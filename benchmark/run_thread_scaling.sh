#!/usr/bin/env bash
# ================================================================
#  run_benchmark.sh — fast-mic benchmark suite
#
#  Part 1: Correctness  — fast-mic vs COBRApy growth rates (Pearson r, R², MAE)
#  Part 2: Speed / thread scaling — wall time, memory, load/FBA breakdown
#
#  Usage:
#    ./run_benchmark.sh
#    MODEL_DIR=./models MEDIUM=WesternDiet N_LIST=100,1000 ./run_benchmark.sh
#
#  Environment variables:
#    MODEL_DIR        directory with .xml models        (default: test/benchmark_uhgg/benchmark_9)
#    MEDIA_DB         media database TSV (BiGG-style)   (default: media/media_db.tsv)
#    MEDIUM           medium name in MEDIA_DB           (default: WesternDiet)
#    MEDIUM_FILE      gapseq SEED-format CSV (compounds,name,maxFlux) — when set,
#                     overrides MEDIA_DB+MEDIUM; required for ModelSEED/gapseq models.
#                     e.g. MEDIUM_FILE=media/western_diet_mucin_gapseq.csv
#    OUTDIR           output directory                  (default: benchmark_results)
#    REPEATS          timed runs per (thread, n) config (default: 3)
#    N_LIST           comma-separated model counts      (default: MAX_MODELS)
#                     e.g. N_LIST=100,1000 runs two curves in the scaling figure
#    MAX_MODELS       fallback single model count, 0=all (default: 0)
#    CORRECTNESS_N    models for correctness check      (default: 50)
#    SKIP_COBRA       skip COBRApy speed baseline       (default: 0)
#    SKIP_CORRECTNESS skip correctness check            (default: 0)
# ================================================================

set -euo pipefail
trap 'echo "ERROR: script exited at line $LINENO (exit code $?)" >&2' ERR

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Script lives in <project>/benchmark/ — step up to the Cargo workspace root
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

MODEL_DIR="${MODEL_DIR:-${PROJECT_ROOT}/test/benchmark_uhgg/benchmark_9}"
MEDIA_DB="${MEDIA_DB:-${PROJECT_ROOT}/media/media_db.tsv}"
MEDIUM="${MEDIUM:-WesternDiet}"
MEDIUM_FILE="${MEDIUM_FILE:-}"            # gapseq-style CSV; overrides MEDIA_DB+MEDIUM
MEDIA_LIST="${MEDIA_LIST:-}"             # file listing several medium CSV paths (one per
                                        # line); each model is loaded once and evaluated
                                        # under EVERY medium (e.g. L0-L9 gradient).
                                        # Overrides MEDIUM_FILE / MEDIA_DB+MEDIUM.
OUTDIR="${OUTDIR:-${SCRIPT_DIR}/benchmark_results}"
REPEATS="${REPEATS:-3}"
N_LIST="${N_LIST:-}"
MAX_MODELS="${MAX_MODELS:-0}"
CORRECTNESS_N="${CORRECTNESS_N:-50}"
COBRA_WORKERS="${COBRA_WORKERS:-4}"   # parallel COBRApy workers for correctness check
SKIP_COBRA="${SKIP_COBRA:-0}"
SKIP_CORRECTNESS="${SKIP_CORRECTNESS:-0}"

RUST_BIN="${PROJECT_ROOT}/target/release/bench-single-fba"
COBRA_SCRIPT="${SCRIPT_DIR}/benchmark_cobra.py"
ACCURACY_SCRIPT="${PROJECT_ROOT}/scripts/accuracy_check.py"
SCALING_DIR="${OUTDIR}/thread_scaling"
CORRECTNESS_DIR="${OUTDIR}/correctness"

# ── Build ──
echo "Building bench-single-fba ..."
(cd "$PROJECT_ROOT" && cargo build --release --bin bench-single-fba 2>&1 | tail -3)
echo ""

mkdir -p "$SCALING_DIR" "$CORRECTNESS_DIR"

# ── Build invocation argument arrays for fast-mic / COBRApy ──
# When MEDIUM_FILE is set, the gapseq CSV defines the medium directly
# (no <media_db> <medium_name> positionals). Otherwise fall back to the
# legacy TSV+name pair.
# Tiered uptake limits are fixed defaults inside the binaries (carbon 10 /
# aa 1 / nuc 0.5 / cofactor 0.1 mmol/gDW/h). Per-compound maxFlux from a
# MEDIUM_FILE CSV still overrides them.
if [ -n "$MEDIA_LIST" ]; then
    if [ ! -f "$MEDIA_LIST" ]; then
        echo "ERROR: MEDIA_LIST='$MEDIA_LIST' not found" >&2; exit 1
    fi
    # Each model is loaded once and evaluated under every medium in the list
    # (one output row per model×medium). Correctness then compares fast-mic vs
    # COBRApy per medium; the rest of the benchmark is unchanged.
    RUST_MEDIUM_ARGS=( --media-list "$MEDIA_LIST" )
    COBRA_MEDIUM_ARGS=( --media-list "$MEDIA_LIST" )
    N_MEDIA=$(grep -cve '^[[:space:]]*$' -e '^[[:space:]]*#' "$MEDIA_LIST" 2>/dev/null || wc -l < "$MEDIA_LIST")
    MEDIUM_DESC="MEDIA_LIST=$MEDIA_LIST (${N_MEDIA} media, per-model loop)"
elif [ -n "$MEDIUM_FILE" ]; then
    if [ ! -f "$MEDIUM_FILE" ]; then
        echo "ERROR: MEDIUM_FILE='$MEDIUM_FILE' not found" >&2; exit 1
    fi
    RUST_MEDIUM_ARGS=( --medium-file "$MEDIUM_FILE" )
    COBRA_MEDIUM_ARGS=( --medium-file "$MEDIUM_FILE" )
    MEDIUM_DESC="MEDIUM_FILE=$MEDIUM_FILE (SEED-format, gapseq-compatible)"
else
    RUST_MEDIUM_ARGS=( "$MEDIA_DB" "$MEDIUM" )
    COBRA_MEDIUM_ARGS=( "$MEDIA_DB" "$MEDIUM" )
    MEDIUM_DESC="MEDIA_DB=$MEDIA_DB  MEDIUM=$MEDIUM"
fi

# ── Collect full model list ──
FULL_LIST="$SCALING_DIR/full_model_list.txt"
find "$MODEL_DIR" -maxdepth 1 -name '*.xml' -type f | sort > "$FULL_LIST"
N_TOTAL=$(wc -l < "$FULL_LIST" | tr -d ' ')

if [ "$N_TOTAL" -eq 0 ]; then
    echo "ERROR: No .xml files found in $MODEL_DIR" >&2
    exit 1
fi

# ── Determine thread counts ──
MAX_THREADS=$(nproc 2>/dev/null || sysctl -n hw.logicalcpu 2>/dev/null || echo 8)
THREAD_COUNTS="1"
t=2
while [ "$t" -le "$MAX_THREADS" ]; do
    THREAD_COUNTS="$THREAD_COUNTS $t"
    t=$((t * 2))
done
last=$(echo "$THREAD_COUNTS" | awk '{print $NF}')
if [ "$last" -ne "$MAX_THREADS" ]; then
    THREAD_COUNTS="$THREAD_COUNTS $MAX_THREADS"
fi

# ── Determine sizes to benchmark ──
if [ -n "$N_LIST" ]; then
    SIZES=$(echo "$N_LIST" | tr ',' ' ')
elif [ "$MAX_MODELS" -gt 0 ] && [ "$MAX_MODELS" -lt "$N_TOTAL" ]; then
    SIZES="$MAX_MODELS"
else
    SIZES="$N_TOTAL"
fi

echo "╔══════════════════════════════════════════════════════╗"
echo "║  fast-mic Benchmark Suite                          ║"
echo "╚══════════════════════════════════════════════════════╝"
echo ""
echo "  Model dir:     $MODEL_DIR  ($N_TOTAL models available)"
echo "  Medium:        $MEDIUM_DESC"
echo "  Sizes (N):     $SIZES"
echo "  Correctness N: $CORRECTNESS_N"
echo "  Thread counts: $THREAD_COUNTS"
echo "  Repeats:       $REPEATS"
echo "  Output:        $OUTDIR"
echo ""

# ── Results TSV — append-safe; write header only if file is new or empty ──
RESULTS_TSV="$SCALING_DIR/thread_scaling_results.tsv"
# Always start fresh: overwrite the header each run so results never accumulate
# stale corpus sizes (e.g. a previous default benchmark_9 run leaving an n=9 curve).
printf "tool\tthreads\trun\twall_s\tn_models\tmodels_per_s\tpeak_mem_mb\tavg_load_s\tavg_fba_s\n" \
    > "$RESULTS_TSV"

# ── Helpers ──
extract_wall_s() {
    grep "Total wall time" "$1" 2>/dev/null \
        | awk '{for(i=1;i<=NF;i++) if($i=="time:") print $(i+1)}' | head -1 || echo "0"
}
extract_mem_mb() {
    grep "Peak memory" "$1" 2>/dev/null \
        | grep -o '[0-9.]*[[:space:]]*MB' | grep -o '[0-9.]*' | head -1 || echo "0"
}
elapsed_hms() { printf "%dm%02ds" $(($1/60)) $(($1%60)); }

# avg of column $1 (1-based) in a header-having TSV
avg_col() {
    local col=$1 tsv=$2
    awk -F'\t' -v c="$col" 'NR>1 && $c!="" && $c+0>0 {s+=$c; n++} END {
        if(n>0) printf "%.6f", s/n; else print "0"
    }' "$tsv"
}

# ================================================================
# PART 1: CORRECTNESS
# ================================================================
if [ "$SKIP_CORRECTNESS" -eq 0 ]; then
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "  Part 1: Correctness — fast-mic vs COBRApy"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

    CORR_LIST="$CORRECTNESS_DIR/model_list.txt"
    # Sample CORRECTNESS_N models deterministically
    N_CORR_MAX=$(( CORRECTNESS_N < N_TOTAL ? CORRECTNESS_N : N_TOTAL ))
    awk -v n="$N_CORR_MAX" -v seed=99 \
        'BEGIN{srand(seed)} {lines[NR]=$0} END{
            for(i=1;i<=NR;i++){j=int(rand()*(i))+1; t=lines[i]; lines[i]=lines[j]; lines[j]=t}
            for(i=1;i<=n;i++) print lines[i]
        }' "$FULL_LIST" | sort > "$CORR_LIST"
    N_CORR=$(wc -l < "$CORR_LIST" | tr -d ' ')
    echo "  Using $N_CORR models"
    echo ""

    echo "  Running fast-mic (1 thread) ..."
    "$RUST_BIN" "${RUST_MEDIUM_ARGS[@]}" \
        --model-list "$CORR_LIST" --threads 1 \
        > "$CORRECTNESS_DIR/fastmic.tsv" \
        2> "$CORRECTNESS_DIR/fastmic.log"
    echo "  Done."

    echo "  Running COBRApy ($COBRA_WORKERS workers, ${MODEL_TIMEOUT_S:-120}s timeout/model) ..."
    python3 "$COBRA_SCRIPT" "${COBRA_MEDIUM_ARGS[@]}" \
        --model-list "$CORR_LIST" \
        --workers "$COBRA_WORKERS" \
        > "$CORRECTNESS_DIR/cobra.tsv" \
        2> "$CORRECTNESS_DIR/cobra.log"
    echo "  Done."

    echo ""
    python3 "$ACCURACY_SCRIPT" \
        "$CORRECTNESS_DIR/fastmic.tsv" \
        "$CORRECTNESS_DIR/cobra.tsv" \
        "$CORRECTNESS_DIR/stats.tsv" \
        "$CORRECTNESS_DIR/scatter.tsv"
    echo ""
fi

# ================================================================
# PART 2: SPEED / THREAD-SCALING
# ================================================================
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Part 2: Speed — thread scaling"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

for SIZE in $SIZES; do
    echo "┌── N = $SIZE models ─────────────────────────────────────"

    # Build model list for this size
    MODEL_LIST="$SCALING_DIR/model_list_n${SIZE}.txt"
    if [ "$SIZE" -lt "$N_TOTAL" ]; then
        awk -v n="$SIZE" -v seed=42 \
            'BEGIN{srand(seed)} {lines[NR]=$0} END{
                for(i=1;i<=NR;i++){j=int(rand()*(i))+1; t=lines[i]; lines[i]=lines[j]; lines[j]=t}
                for(i=1;i<=n;i++) print lines[i]
            }' "$FULL_LIST" | sort > "$MODEL_LIST"
        echo "  Sampled $SIZE from $N_TOTAL models"
    else
        cp "$FULL_LIST" "$MODEL_LIST"
        echo "  Using all $N_TOTAL models"
    fi

    # ── Probe for time estimate ──
    PROBE_LIST="$SCALING_DIR/probe_n${SIZE}.txt"
    PROBE_N=$(( SIZE < 10 ? SIZE : 10 ))
    head -"$PROBE_N" "$MODEL_LIST" > "$PROBE_LIST"
    PROBE_START=$(date +%s)
    "$RUST_BIN" "${RUST_MEDIUM_ARGS[@]}" --model-list "$PROBE_LIST" --threads 1 \
        > /dev/null 2> "$SCALING_DIR/probe_n${SIZE}.log"
    PROBE_END=$(date +%s)
    PROBE_S=$((PROBE_END - PROBE_START))
    EST_1T=$(awk "BEGIN {printf \"%d\", ($SIZE / $PROBE_N.0) * $PROBE_S}")
    printf "  Probe: %ds for %d models → est. %s at 1 thread\n" \
        "$PROBE_S" "$PROBE_N" "$(elapsed_hms $EST_1T)"
    echo ""

    # ── fast-mic thread scaling ──
    echo "  fast-mic:"
    for T in $THREAD_COUNTS; do
        printf "    threads=%-3s  " "$T"
        for run in $(seq 1 "$REPEATS"); do
            LOG="$SCALING_DIR/fastmic_t${T}_n${SIZE}_r${run}.log"
            # Save per-model TSV for run=1 (used for stacked bar in plot)
            if [ "$run" -eq 1 ]; then
                PERMODEL="$SCALING_DIR/fastmic_t${T}_n${SIZE}_r1_permodel.tsv"
            else
                PERMODEL="/dev/null"
            fi

            "$RUST_BIN" "${RUST_MEDIUM_ARGS[@]}" \
                --model-list "$MODEL_LIST" --threads "$T" \
                > "$PERMODEL" \
                2> "$LOG"

            WALL=$(extract_wall_s "$LOG")
            MEM=$(extract_mem_mb "$LOG")
            MOPS=$(awk "BEGIN {printf \"%.2f\", $SIZE / ($WALL + 1e-9)}")

            # avg load/fba from per-model TSV (run 1 only; others use 0)
            if [ "$run" -eq 1 ] && [ "$PERMODEL" != "/dev/null" ]; then
                AVG_LOAD=$(avg_col 7 "$PERMODEL")
                AVG_FBA=$(avg_col 8 "$PERMODEL")
            else
                AVG_LOAD="0"
                AVG_FBA="0"
            fi

            printf "fast-mic\t$T\t$run\t$WALL\t$SIZE\t$MOPS\t$MEM\t$AVG_LOAD\t$AVG_FBA\n" \
                >> "$RESULTS_TSV"
            printf "r%d:%.1fs  " "$run" "$WALL"
        done
        echo ""
    done

    # ── COBRApy thread scaling (parallel via multiprocessing, HiGHS solver) ──
    if [ "$SKIP_COBRA" -eq 0 ]; then
        echo ""
        echo "  COBRApy (HiGHS, parallel — same thread counts as fast-mic):"
        for T in $THREAD_COUNTS; do
            printf "    threads=%-3s  " "$T"
            for run in $(seq 1 "$REPEATS"); do
                LOG="$SCALING_DIR/cobra_t${T}_n${SIZE}_r${run}.log"
                if [ "$run" -eq 1 ]; then
                    PERMODEL="$SCALING_DIR/cobra_t${T}_n${SIZE}_r1_permodel.tsv"
                else
                    PERMODEL="/dev/null"
                fi
                python3 "$COBRA_SCRIPT" "${COBRA_MEDIUM_ARGS[@]}" \
                    --model-list "$MODEL_LIST" \
                    --threads "$T" \
                    > "$PERMODEL" \
                    2> "$LOG"
                WALL=$(extract_wall_s "$LOG")
                MEM=$(extract_mem_mb "$LOG")
                MOPS=$(awk "BEGIN {printf \"%.3f\", $SIZE / ($WALL + 1e-9)}")
                if [ "$run" -eq 1 ] && [ "$PERMODEL" != "/dev/null" ]; then
                    AVG_LOAD=$(avg_col 7 "$PERMODEL")
                    AVG_FBA=$(avg_col 8 "$PERMODEL")
                else
                    AVG_LOAD="0"; AVG_FBA="0"
                fi
                printf "COBRApy\t$T\t$run\t$WALL\t$SIZE\t$MOPS\t$MEM\t$AVG_LOAD\t$AVG_FBA\n" \
                    >> "$RESULTS_TSV"
                printf "r%d:%.1fs  " "$run" "$WALL"
            done
            echo ""
        done
    fi

    echo "└────────────────────────────────────────────────────"
    echo ""
done

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Benchmark complete"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Results: $RESULTS_TSV"
if [ "$SKIP_CORRECTNESS" -eq 0 ]; then
    echo "  Correctness: $CORRECTNESS_DIR/stats.tsv"
fi
echo ""
echo "  Generate figure:"
echo "    Rscript ${PROJECT_ROOT}/scripts/plot_fig1_benchmark.R \\"
echo "      $RESULTS_TSV \\"
echo "      $CORRECTNESS_DIR/scatter.tsv \\"
echo "      $CORRECTNESS_DIR/stats.tsv"
