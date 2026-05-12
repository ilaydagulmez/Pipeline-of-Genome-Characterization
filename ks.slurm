#!/bin/bash

set -euo pipefail

# -----------------------------
# tools
# -----------------------------
MAFFT=/path/to/mafft
PAL2NAL=/path/to/pal2nal.pl
KAKS=path/to/KaKs_Calculator

# -----------------------------
# input files
# -----------------------------
PAIRS=/path/to/pairs
CDS=/path/to/cds
PEP=/path/to/peptide

# sonuç dosyası /arf üstünde kalsın
OUTFILE="$SLURM_SUBMIT_DIR/genelevel.tsv"

# local temp alanı (kota dostu)
WORKDIR="${SLURM_TMPDIR:-/tmp}/${USER}_ks_${SLURM_JOB_ID}"
mkdir -p "$WORKDIR"/tmp "$WORKDIR"/aln "$WORKDIR"/codon "$WORKDIR"/kaks

echo "[INFO] JOB ID   : $SLURM_JOB_ID"
echo "[INFO] HOST     : $(hostname)"
echo "[INFO] PWD      : $(pwd)"
echo "[INFO] SUBMIT   : $SLURM_SUBMIT_DIR"
echo "[INFO] WORKDIR  : $WORKDIR"
echo "[INFO] PAIRS    : $PAIRS"
echo "[INFO] CDS      : $CDS"
echo "[INFO] PEP      : $PEP"

ls -lh "$PAIRS" "$CDS" "$PEP"

# -----------------------------
# functions
# -----------------------------
get_fa () {
  local id="$1"
  local fa="$2"

  awk -v id="$id" '
    BEGIN {keep=0}
    /^>/{
      header=substr($0,2)
      split(header,a,/[\t ]+/)
      keep=(a[1]==id)
      if (keep) print ">" a[1]
      next
    }
    keep {print}
  ' "$fa"
}

fasta_to_axt () {
  local fa="$1"
  local axt="$2"

  awk '
    /^>/ {n++; next}
    {
      gsub(/[[:space:]]+/, "", $0)
      $0=toupper($0)
      gsub(/[^ACGT-]/, "", $0)

      if (n==1) s1=s1 $0
      else if (n==2) s2=s2 $0
    }
    END {
      if (length(s1)>0 && length(s2)>0 && length(s1)==length(s2)) {
        print "1 2"
        print s1
        print s2
      }
    }
  ' "$fa" > "$axt"
}

cleanup_pair_files () {
  rm -f \
    "$WORKDIR/tmp/p1.pep.fa" \
    "$WORKDIR/tmp/p2.pep.fa" \
    "$WORKDIR/tmp/p1.cds.fa" \
    "$WORKDIR/tmp/p2.cds.fa" \
    "$WORKDIR/tmp/pair.pep.fa" \
    "$WORKDIR/tmp/pair.cds.fa" \
    "$WORKDIR/aln/pair.pep.aln" \
    "$WORKDIR/codon/pair.codon.fa" \
    "$WORKDIR/tmp/pair.axt" \
    "$WORKDIR/kaks/pair.kaks"
}

# -----------------------------
# initialize output
# -----------------------------
echo -e "gene1\tgene2\tKs" > "$OUTFILE"

total=0
missing_pep=0
missing_cds=0
bad_axt=0
kaks_fail=0
written=0

# -----------------------------
# main loop
# -----------------------------
while IFS=$'\t' read -r g1 g2; do
  [[ -n "${g1:-}" && -n "${g2:-}" ]] || continue

  total=$((total+1))

  # Windows CRLF varsa temizle
  g1="${g1//$'\r'/}"
  g2="${g2//$'\r'/}"

  if (( total <= 5 )); then
    echo "[DEBUG] pair $total: $g1  $g2"
  fi

  cleanup_pair_files

  # protein
  get_fa "$g1" "$PEP" > "$WORKDIR/tmp/p1.pep.fa"
  if [[ ! -s "$WORKDIR/tmp/p1.pep.fa" ]]; then
    missing_pep=$((missing_pep+1))
    continue
  fi

  get_fa "$g2" "$PEP" > "$WORKDIR/tmp/p2.pep.fa"
  if [[ ! -s "$WORKDIR/tmp/p2.pep.fa" ]]; then
    missing_pep=$((missing_pep+1))
    continue
  fi

  # cds
  get_fa "$g1" "$CDS" > "$WORKDIR/tmp/p1.cds.fa"
  if [[ ! -s "$WORKDIR/tmp/p1.cds.fa" ]]; then
    missing_cds=$((missing_cds+1))
    continue
  fi

  get_fa "$g2" "$CDS" > "$WORKDIR/tmp/p2.cds.fa"
  if [[ ! -s "$WORKDIR/tmp/p2.cds.fa" ]]; then
    missing_cds=$((missing_cds+1))
    continue
  fi

  cat "$WORKDIR/tmp/p1.pep.fa" "$WORKDIR/tmp/p2.pep.fa" > "$WORKDIR/tmp/pair.pep.fa"
  cat "$WORKDIR/tmp/p1.cds.fa" "$WORKDIR/tmp/p2.cds.fa" > "$WORKDIR/tmp/pair.cds.fa"

  # mafft
  "$MAFFT" --leavegappyregion --retree 1 --maxiterate 0 "$WORKDIR/tmp/pair.pep.fa" > "$WORKDIR/aln/pair.pep.aln" || continue

  # pal2nal
  perl "$PAL2NAL" "$WORKDIR/aln/pair.pep.aln" "$WORKDIR/tmp/pair.cds.fa" -output fasta > "$WORKDIR/codon/pair.codon.fa" || continue

  # axt
  fasta_to_axt "$WORKDIR/codon/pair.codon.fa" "$WORKDIR/tmp/pair.axt"

  if [[ ! -s "$WORKDIR/tmp/pair.axt" ]]; then
    bad_axt=$((bad_axt+1))
    continue
  fi

  # kaks
  "$KAKS" -i "$WORKDIR/tmp/pair.axt" -o "$WORKDIR/kaks/pair.kaks" -m YN -c 1 || {
    kaks_fail=$((kaks_fail+1))
    continue
  }

  if [[ ! -s "$WORKDIR/kaks/pair.kaks" ]]; then
    kaks_fail=$((kaks_fail+1))
    continue
  fi

  ks=$(awk '
    NR==1 {
      for(i=1;i<=NF;i++) if($i=="Ks") c=i
    }
    NR==2 && c {
      print $c
    }
  ' "$WORKDIR/kaks/pair.kaks")

  if [[ -n "${ks:-}" && "$ks" != "NA" ]]; then
    echo -e "$g1\t$g2\t$ks" >> "$OUTFILE"
    written=$((written+1))
  fi

  if (( total % 100 == 0 )); then
    echo "[INFO] processed=$total written=$written missing_pep=$missing_pep missing_cds=$missing_cds bad_axt=$bad_axt kaks_fail=$kaks_fail"
  fi

done < "$PAIRS"

cleanup_pair_files

echo "[INFO] DONE"
echo "[INFO] total_pairs   = $total"
echo "[INFO] written_Ks    = $written"
echo "[INFO] missing_pep   = $missing_pep"
echo "[INFO] missing_cds   = $missing_cds"
echo "[INFO] bad_axt       = $bad_axt"
echo "[INFO] kaks_fail     = $kaks_fail"
echo "[INFO] outfile       = $OUTFILE"
