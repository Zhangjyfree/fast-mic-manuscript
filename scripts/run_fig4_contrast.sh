#!/usr/bin/env bash
# Controlled carbon-quality contrast (Figure 4G, Table S10).
#
# Screens both probiotic panels against the UHGG community on three media that
# share the L1–L4 background: L5 (pectin), L5-glc (pectin replaced by free
# glucose) and L6 (resistant starch: free glucose + maltose + maltodextrin).
# L5 and L6 reproduce the corresponding files in results/{akk,lac}_vs_uhgg/.
#
# Run from the repository root after extracting test/UHGG/final_gapseq_xml
# (README). The fast-mic binary and the media (including the L5-glc control)
# come from the engine checkout, $FASTMIC_ENGINE (default ../fast-mic).
#
# Usage:  bash scripts/run_fig4_contrast.sh [threads]      (default: 0 = all cores)
# Output: results/fig4/{akk,lac}_{L5_pectin,L5glc,L6_resistant_starch}.tsv

set -euo pipefail

ENGINE="${FASTMIC_ENGINE:-../fast-mic}"
BIN="$ENGINE/target/release/fast-mic"
THREADS="${1:-0}"
PARTNERS="test/UHGG/final_gapseq_xml"

[[ -x "$BIN" ]] || { echo "fast-mic binary not found at $BIN (set FASTMIC_ENGINE)" >&2; exit 1; }
[[ -d "$PARTNERS" ]] || { echo "$PARTNERS not found: extract the UHGG archive first (README)" >&2; exit 1; }

mkdir -p results/fig4
for sys in akk lac; do
    case "$sys" in
        akk) models="test/akk/akk_gapseq_xml" ;;
        lac) models="test/lac/lac_genomes_faa_gapseq_wdm_xml" ;;
    esac
    for level in L5_pectin L5glc L6_resistant_starch; do
        echo "[$(date +%H:%M:%S)] $sys × UHGG on $level"
        "$BIN" --group1 "$models" --group2 "$PARTNERS" \
               --medium-file "$ENGINE/media/gradient_${level}_gapseq.csv" \
               --threads "$THREADS" -o "results/fig4/${sys}_${level}.tsv" > /dev/null
    done
done
echo "done -> results/fig4/"
