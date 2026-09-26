#!/usr/bin/env bash
# Run fast-mic across the ten-level prebiotic gradient (L0 base -> L9 MOS) and
# collect per-level pairwise interaction results. The raw per-pair tables in
# results/{akk,lac}_vs_uhgg/ of this repository were produced with this script.
#
# Run it from the repository root. The fast-mic binary and the media come from
# the engine checkout, $FASTMIC_ENGINE (default ../fast-mic): the ten levels are
# those listed in its media/gradient_media_list.txt. test/UHGG/final_gapseq_xml
# must be extracted first (README).
#
# Usage:
#   # Akkermansia: models gap-filled on the Akkermansia minimal medium
#   bash scripts/run_gradient.sh \
#        --group1 test/akk/akk_gapseq_xml \
#        --group2 test/UHGG/final_gapseq_xml \
#        --threads 12 --full-tsv --out results/akk_vs_uhgg
#
#   # Lactobacillus-group: models gap-filled on Western diet + mucin
#   bash scripts/run_gradient.sh \
#        --group1 test/lac/lac_genomes_faa_gapseq_wdm_xml \
#        --group2 test/UHGG/final_gapseq_xml \
#        --threads 12 --full-tsv \
#        --out results/lac_vs_uhgg          # <- explicit output directory
#
# Options:
#   --models DIR        all-vs-all within a single group
#   --group1 DIR        first group of a cross-group run
#   --group2 DIR        second group of a cross-group run
#   --threads N         number of threads (0 = all cores)
#   --full-tsv          also write *.full.tsv (per-metabolite cross-feeding)
#   --out DIR           output directory (default: results/<group1>_vs_<group2>, see below)
#
# Output:
#   Without --out, the directory is derived from the group names:
#     group1=akk, group2=UHGG  →  results/akk_vs_uhgg/
#     group1=lac, group2=UHGG  →  results/lac_vs_uhgg/
#     --models akk             →  results/akk/
#   Directory contents:
#     L0_base.tsv ... L9_mos.tsv
#     L*.full.tsv   (with --full-tsv)
#     gradient_run.log

set -euo pipefail

ENGINE="${FASTMIC_ENGINE:-../fast-mic}"   # engine checkout, relative to the repository root
BIN="$ENGINE/target/release/fast-mic"
MEDIA_LIST="$ENGINE/media/gradient_media_list.txt" # L0–L9, paths relative to the engine
OUT_DIR=""                       # empty = derive automatically (see below)

# ── CLI options ─────────────────────────────────────────────────
MODELS=""
GROUP1=""
GROUP2=""
THREADS=0
WANT_FULL_TSV=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --models)      MODELS="$2"; shift 2 ;;
        --group1)      GROUP1="$2"; shift 2 ;;
        --group2)      GROUP2="$2"; shift 2 ;;
        --threads)     THREADS="$2"; shift 2 ;;
        --out)         OUT_DIR="$2"; shift 2 ;;
        --full-tsv)    WANT_FULL_TSV=1; shift ;;
        -h|--help)
            sed -n '2,40p' "$0"; exit 0 ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

# ── derive a short label from a model directory ──────
# If the leaf directory name is generic (contains gapseq / xml / bacteria / final), use its parent.
#   test/akk/akk_gapseq_xml                 -> akk
#   test/UHGG/final_gapseq_xml              -> UHGG
#   test/lac/lac_genomes_faa_gapseq_wdm_xml -> lac
group_label() {
    local d leaf parent
    d="${1%/}"                         # strip trailing slash
    leaf="$(basename "$d")"
    parent="$(basename "$(dirname "$d")")"
    if [[ "$leaf" =~ (gapseq|xml|bacteria|final|sbml) ]]; then
        echo "$parent"
    else
        echo "$leaf"
    fi
}

# ── auto-derive OUT_DIR (only when --out was not given) ──────
if [[ -z "$OUT_DIR" ]]; then
    lc() { tr '[:upper:]' '[:lower:]' <<<"$1"; }
    if [[ -n "$MODELS" ]]; then
        OUT_DIR="results/$(lc "$(group_label "$MODELS")")"
    elif [[ -n "$GROUP1" && -n "$GROUP2" ]]; then
        OUT_DIR="results/$(lc "$(group_label "$GROUP1")")_vs_$(lc "$(group_label "$GROUP2")")"
    else
        OUT_DIR="results/gradient_result"   # fallback (validation below will error out)
    fi
fi

# Validate
if [[ -z "$MODELS" && ( -z "$GROUP1" || -z "$GROUP2" ) ]]; then
    echo "Error: provide either --models DIR or both --group1 DIR --group2 DIR" >&2
    exit 1
fi
if [[ ! -x "$BIN" ]]; then
    echo "Error: fast-mic binary not found at $BIN — build the engine (cargo build --release) or set FASTMIC_ENGINE." >&2
    exit 1
fi
if [[ ! -f "$MEDIA_LIST" ]]; then
    echo "Error: $MEDIA_LIST not found — run from the repository root, or set FASTMIC_ENGINE." >&2
    exit 1
fi

mkdir -p "$OUT_DIR"
LOG="$OUT_DIR/gradient_run.log"
: > "$LOG"
echo "===========================================================" | tee -a "$LOG"
echo " Gradient run started $(date -u +%FT%TZ)" | tee -a "$LOG"
echo " OUTPUT DIR ▶  $OUT_DIR" | tee -a "$LOG"
echo "===========================================================" | tee -a "$LOG"

# ── Species-set summary ─────────────────────────────────────────
if [[ -n "$MODELS" ]]; then
    N=$(find "$MODELS" -maxdepth 1 \( -name '*.xml' -o -name '*.sbml' \) | wc -l | tr -d ' ')
    echo "Mode: all-vs-all on $MODELS  ($N models, $((N*(N-1)/2)) pairs per medium)" | tee -a "$LOG"
else
    N1=$(find "$GROUP1" -maxdepth 1 \( -name '*.xml' -o -name '*.sbml' \) | wc -l | tr -d ' ')
    N2=$(find "$GROUP2" -maxdepth 1 \( -name '*.xml' -o -name '*.sbml' \) | wc -l | tr -d ' ')
    echo "Mode: cross-group $GROUP1 x $GROUP2  ($N1 x $N2 = $((N1*N2)) pairs per medium)" | tee -a "$LOG"
fi
echo "Threads: $THREADS" | tee -a "$LOG"
echo "Output:  $OUT_DIR" | tee -a "$LOG"
echo "" | tee -a "$LOG"

# ── Main loop over the gradient media (L0–L9) ────────────────────────────
while read -r rel; do
    [[ -z "$rel" ]] && continue
    med="$ENGINE/$rel"
    # extract short label
    label=$(basename "$med" .csv | sed 's/^gradient_//;s/_gapseq$//')
    out="$OUT_DIR/${label}.tsv"
    full="$OUT_DIR/${label}.full.tsv"

    n_cpd=$(($(wc -l < "$med") - 1))
    echo "[$(date +%H:%M:%S)] Running $label ($n_cpd compounds) → $out" | tee -a "$LOG"

    args=( --medium-file "$med" --threads "$THREADS" -o "$out" )
    if [[ $WANT_FULL_TSV -eq 1 ]]; then
        args+=( --full-tsv "$full" )
    fi
    if [[ -n "$MODELS" ]]; then
        # let fast-mic glob every SBML in the directory
        # shellcheck disable=SC2206
        models_arr=( "$MODELS"/*.xml "$MODELS"/*.sbml )
        # drop unmatched globs
        models_clean=()
        for m in "${models_arr[@]}"; do [[ -f "$m" ]] && models_clean+=( "$m" ); done
        args+=( --summary "${models_clean[@]}" )
    else
        args+=( --group1 "$GROUP1" --group2 "$GROUP2" --summary )
    fi

    # time + log
    /usr/bin/time -p "$BIN" "${args[@]}" >>"$LOG" 2>&1
    echo "    done." | tee -a "$LOG"
done < "$MEDIA_LIST"

echo "" | tee -a "$LOG"
echo "=== Gradient run finished $(date -u +%FT%TZ) ===" | tee -a "$LOG"
echo ""
echo "Next step (in this repository):"
echo "  python3 scripts/extract_fig5_crossfeed.py   # cross-feeding prevalence"
echo "  python3 scripts/extract_fig6_enrich.py      # strain bootstrap + partner phyla"
echo "  python3 scripts/extract_interception.py     # higher-order interception bound"
