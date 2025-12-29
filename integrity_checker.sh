#!/usr/bin/env bash
#
# 
# Purpose:
#   Perform lightweight integrity checks on common NGS file types before they
#   enter downstream pipelines. Exactly ONE line per file is printed:
#     [OK] / [WARNING] / [ERROR] <TYPE> <PATH> - <message>
#
# Checks implemented
#   • FASTQ : extract the first 40 000 lines and run fastQValidator
#   • BAM   : inspect the first 500 header lines and verify the BAM EOF marker
#   • CRAM  : inspect the first 500 header lines and REQUIRE reference
#             MD5 (M5) tags (missing M5 ⇒ ERROR)
#   • VCF   : parse the header plus the first 10 000 variant records with
#             bcftools head (fatal parse ⇒ ERROR)
#           : verify that VCF/BCF files are sorted according to HTSlib rules
#
# Exit status
#   • Any ERROR terminates execution immediately (set -e).
#   • WARNING is non-fatal (e.g., missing helper tool); processing continues.
#   • OK indicates the file passed the implemented checks.
# ---------------------------------------------------------------------------

set +e  # aggregates errors; we do not want to abort afer the first failure
set -o pipefail

declare -a _OKS _FAILS _ERRS
_cur_type=""
_cur_file=""

begin_file() {
    _cur_type="$1"
    _cur_file="$2"
    _OKS=()
    _FAILS=()
    _ERRS=()
}

ok()   { _OKS+=("[OK] $1 $2 - $3"); }
fail() { _FAILS+=("[FAIL] $1 $2 - $3"); }
err()  { _ERRS+=("[ERROR] $1 $2 - $3"); }

end_file() {
  if ((${#_ERRS[@]})); then
    printf '%s\n' "${_ERRS[@]}"
    return 1
  fi

  if ((${#_FAILS[@]})); then
    printf '%s\n' "${_FAILS[@]}"
    return 0
  fi

  printf '[OK] %s %s - file OK\n' "$_cur_type" "$_cur_file"
  return 0
}

FASTQ_LINES=40000      # FASTQ lines inspected
BAM_EOF_BYTES=32768    # bytes read from end of BAM for EOF validation
VCF_RECORDS=10000      # VCF/BCF records parsed
HUM_REF_GENOME=("hg16" "hg17" "hg18" "GRCh37" "GRCh38" "T2T")  # human reference genome names


##############################################################################
# FASTQ
##############################################################################
check_fastq() {
    local f="$1" type="FASTQ"
    begin_file "$type" "$f"

    local tmp
    tmp=$(mktemp)

    if [[ "$f" == *.gz ]]; then
        if ! gunzip -c "$f" >/dev/null | head -n "$FASTQ_LINES" > "$tmp"; then
            err "$type" "$f" "failed to read/decompress FASTQ"
            rm -f "$tmp"
            end_file; return $?
        fi
    else
        if ! head -n "$FASTQ_LINES" "$f" > "$tmp" 2>/dev/null; then
        err "$type" "$f" "failed to read FASTQ"
        rm -f "$tmp"
        end_file; return $?
    fi
  fi

    #validator check 
    if ! command -v fastQValidator >/dev/null 2>&1; then
        rm -f "$tmp"
        fail "$type" "$f" "fastQValidator not found; FASTQ check skipped"
        end_file; return $?
    fi

    if fastQValidator --file "$tmp" --maxErrors 1 --disableSeqIDCheck >/dev/null 2>&1; then
        ok "$type" "$f" "validator passed on first ${FASTQ_LINES} lines"
    else       
        err "$type" "$f" "fastQValidator reported format problems"
    fi

    rm -f "$tmp"
    end_file; return $?
}

##############################################################################
# BAM
##############################################################################
check_bam() {
    local f="$1" type="BAM"
    begin_file "$type" "$f"

    local hdr_tmp tail_tmp sorted species

    # EOF marker
    tail_tmp=$(mktemp)
    tail -c "$BAM_EOF_BYTES" "$f" > "$tail_tmp" 2>/dev/null || true

    if xxd -p "$tail_tmp" | grep -iq '42430200'; then
        ok "$type" "$f" "BAM EOF magic present"
    else
        err "$type" "$f" "BAM EOF magic (42430200) not found"
    fi
    rm -f "$tail_tmp"


    # ----- samtools check -----
    
    if ! command -v samtools >/dev/null 2>&1; then
        fail "$type" "$f" "samtools not found; header/sortedness checks skipped"
        end_file
        return $?
    fi

    # Checks if header is readable 
    hdr_tmp=$(mktemp)

    if ! samtools view -H "$f" > "$hdr_tmp" 2>/dev/null; then
        err "$type" "$f" "BAM header missing or unreadable"
        rm -f "$hdr_tmp"
        end_file
        return $?
    fi
    
    if ! grep -q '^@SQ' "$hdr_tmp"; then
        err "$type" "$f" "BAM header missing @SQ lines"
        rm -f "$hdr_tmp"
        end_file
        return $?
    fi
    
    # Sortedness by coordinate
    sorted=$(grep -m1 '^@HD' "$hdr_tmp" | grep -oE "SO:[^[:space:]]*" | cut -d: -f2)
    if [[ "$sorted" == "coordinate" ]] ; then
        ok "$type" "$f" "BAM file sorted by coordinate"
    else
        err "$type" "$f" "BAM not sorted by coordinate"
    fi

    rm -f "$hdr_tmp"

    # Human reference genome check
    if ! command -v refgenDetector_main.py >/dev/null 2>&1; then
        fail "$type" "$f" "refgenDetector not found"
        end_file; return $?
    fi

    species=$(refgenDetector_main.py -f "$f" -t BAM/CRAM 2>/dev/null \
        | awk -F'Species detected:[[:space:]]*' '/Species detected:/ {print $2}' \
        | xargs)

    if [[ "$species" == "Homo sapiens" ]]; then
        ok "$type" "$f" "BAM file is human"
    else
        err "$type" "$f" "BAM file is not human"
    fi

    end_file; return $?

}


##############################################################################
# CRAM
##############################################################################
check_cram() {
    local f="$1" type="CRAM"
    begin_file "$type" "$f"

    local hdr_tmp sorted species

    #samtools check 
    if ! command -v samtools >/dev/null 2>&1; then
        fail "$type" "$f" "samtools not found; CRAM check skipped"
        end_file; return $?
    fi

    # Read header once 
    hdr_tmp=$(mktemp)
    if ! samtools view -H "$f" > "$hdr_tmp" 2>/dev/null; then
        err "$type" "$f" "CRAM header missing or unreadable"
        rm -f "$hdr_tmp"
        end_file; return $?
    fi

    # Require @SQ (reference dictionary)
    if grep -q '^@SQ' "$hdr_tmp"; then
        ok "$type" "$f" "header readable; @SQ present"
    else
        err  "$type" "$f" "header missing @SQ"
        rm -f "$hdr_tmp"
        end_file; return $?
    fi

    # Sortedness by coordinate (header-based; @HD may be missing)
    sorted=$(grep -m1 '^@HD' "$hdr_tmp" | grep -oE 'SO:[^[:space:]]+' | cut -d: -f2)

    if [[ "$sorted" == "coordinate" ]]; then
        ok "$type" "$f" "sortedness: coordinate"
    else
        err "$type" "$f" "CRAM not sorted by coordinate (SO:${sorted:-missing})"
    fi

    # Reference MD5 tags (now REQUIRED)
    if grep -q 'M5:' "$hdr_tmp"; then
        ok "$type" "$f" "M5 reference MD5 tags present"
    else
        err "$type" "$f" "missing required M5 reference MD5 tags"
    fi

    rm -f "$hdr_tmp"

    # Human reference genome check
    if ! command -v refgenDetector_main.py >/dev/null 2>&1; then
        fail "$type" "$f" "refgenDetector not found"
        end_file; return $?
    fi

    species=$(
        refgenDetector_main.py -f "$f" -t BAM/CRAM 2>/dev/null \
        | awk -F'Species detected:[[:space:]]*' '/Species detected:/ {print $2}' \
        | xargs
    )

    if [[ -z "$species" ]]; then
        err "$type" "$f" "refgenDetector produced no species result"
    elif [[ "$species" == "Homo sapiens" ]]; then
        ok "$type" "$f" "species: Homo sapiens"
    else
        err "$type" "$f" "species is not human ($species)"
    fi

    end_file; return $?
}

##############################################################################
# VCF / BCF
##############################################################################
check_vcf() {
    local f="$1" type="VCF/BCF"

    if ! command -v bcftools >/dev/null 2>&1; then
        fail "$type" "$f" "bcftools not found; VCF/BCF check skipped"
        return
    fi

    if bcftools head -n "$VCF_RECORDS" "$f" >/dev/null 2>&1; then
        ok "$type" "$f" "bcftools parsed header + first ${VCF_RECORDS} records"
    else
        err "$type" "$f" "bcftools failed parsing within first ${VCF_RECORDS} records"
    fi
}

check_variant_sorted() {
    local f="$1" type="VCF/BCF"

    if [[ -z "$f" ]]; then
        fail "$type" "$f" "no file given; sortedness check skipped"
        return
    fi

    # If BCF-like, we need bcftools to convert to VCF text
    if [[ "$f" == *.bcf || "$f" == *.bcf.gz || "$f" == *.bcf.bz2 ]]; then
        if ! command -v bcftools >/dev/null 2>&1; then
            fail "$type" "$f" "bcftools not found; sortedness check skipped for BCF"
            return
        fi
    fi

    # Run the sortedness check; capture first offending locus if unsorted.
    # awk prints "chr:pos" and exits 1 on the first out-of-order record,
    # or prints nothing and exits 0 if everything is sorted.
    local unsorted_at=""
    if ! unsorted_at=$(
        if [[ "$f" == *.bcf || "$f" == *.bcf.gz || "$f" == *.bcf.bz2 ]]; then
            # Any BCF (compressed or not) → VCF stream via bcftools
            bcftools view -Ov "$f" 2>/dev/null
        elif [[ "$f" == *.vcf.bz2 ]]; then
            bzcat "$f" 2>/dev/null
        elif [[ "$f" == *.vcf.gz ]]; then
            zcat "$f" 2>/dev/null
        else
            cat "$f"
        fi | awk -F'\t' '
            BEGIN {
                last_tid = -1
                last_pos = -1
                nctg = 0
            }

            # Collect contigs from header
            /^##contig=/ {
                if (match($0, /ID=([^,>]+)/, a)) {
                    cid = a[1]
                    if (!(cid in tid)) {
                        tid[cid] = nctg
                        nctg++
                    }
                }
                next
            }

            # Skip other header lines
            /^#/ { next }

            # Variant lines
            {
                chrom = $1
                pos   = $2 + 0

                # If chrom not in header, append at end
                if (!(chrom in tid)) {
                    tid[chrom] = nctg
                    nctg++
                }
                t = tid[chrom]

                # Same logic as HTSlib: (tid, pos) must not go backwards
                if (t < last_tid || (t == last_tid && pos < last_pos)) {
                    # Print first offending locus to stdout and exit non-zero
                    printf("%s:%d\n", chrom, pos)
                    exit 1
                }

                last_tid = t
                last_pos = pos
            }
        '
    ); then
        # Non-zero exit from awk / pipeline
        if [[ -n "$unsorted_at" ]]; then
            err "$type" "$f" "variants not sorted by contig/POS (first offending locus: $unsorted_at)"
        else
            err "$type" "$f" "sortedness scan failed (parse or I/O error)"
        fi
    else
        ok "$type" "$f" "variants sorted by contig/POS"
    fi
}





##############################################################################
# Main
##############################################################################
if [[ $# -ne 1 ]]; then
    echo "[ERROR] FILE - usage: $0 <file>"
    exit 1
fi

file="$1"

if [[ ! -f "$file" ]]; then
    echo "[ERROR] FILE $file - not found"
    exit 1
fi



case "$file" in
  *.fastq|*.fastq.gz|*.fq|*.fq.gz)
    check_fastq "$file"
    exit $?
    ;;
  *.bam|*.bam.gz)
    check_bam "$file"
    exit $?
    ;;
  *.cram|*.cram.gz)
    check_cram "$file"
    exit $?
    ;;
  *.vcf|*.vcf.gz|*.bcf|*.bcf.gz|*.vcf.bz2|*.bcf.bz2)
    # keep your current VCF calls for now
    check_vcf "$file"
    rc1=$?
    check_variant_sorted "$file"
    rc2=$?
    # if either fails, fail
    (( rc1 != 0 || rc2 != 0 )) && exit 1
    exit 0
    ;;
  *)
    echo "[WARNING] FILE $file - unsupported extension; skipping"
    exit 0
    ;;
esac

