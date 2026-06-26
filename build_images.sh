#!/usr/bin/env bash
# Build every per-tool benchmark image, one by one. Run from the repo root.
# polygenic builds from the repo root (it vendors models/pgx); the rest from docker/.
set -uo pipefail
LOG=build_images.log; : > "$LOG"

declare -A TAG=(
  [polygenic]=pgxbench/polygenic:local
  [pharmcat]=pgxbench/pharmcat:2.15.5  [pypgx]=pgxbench/pypgx:0.25.0
  [aldy]=pgxbench/aldy:4.6             [stellarpgx]=pgxbench/stellarpgx:1.2.7
  [cyrius]=pgxbench/cyrius:1.1.1       [ursapgx]=pgxbench/ursapgx:1.0.0
  [pgxpop]=pgxbench/pgxpop:1.0         [panno]=pgxbench/panno:0.2.0
  [t1k]=pgxbench/t1k:1.0.5             [stargazer]=pgxbench/stargazer:2.0.2
  [pbstarphase]=pgxbench/pbstarphase:2.1.0  [pangu]=pgxbench/pangu:0.2.8
  [decypher]=pgxbench/decypher:0.1.0   [chinook]=pgxbench/wf-pgx:latest
  [specimmune]=pgxbench/specimmune:1.0.0
)
ORDER=(polygenic pharmcat cyrius aldy pypgx pgxpop panno stargazer t1k stellarpgx ursapgx pbstarphase pangu specimmune decypher chinook)

for t in "${ORDER[@]}"; do
  df="docker/$t.Dockerfile"
  [ -f "$df" ] || { echo "MISSING-DOCKERFILE $t" | tee -a "$LOG"; continue; }
  # polygenic vendors models/pgx -> needs repo-root context; the rest use docker/.
  ctx="docker"; [ "$t" = polygenic ] && ctx="."
  echo "=== building $t (${TAG[$t]}) ===" | tee -a "$LOG"
  if docker build -f "$df" -t "${TAG[$t]}" "$ctx" >>"$LOG" 2>&1; then
    echo "OK   $t" | tee -a "$LOG"
  else
    echo "FAIL $t (see $LOG)" | tee -a "$LOG"
  fi
done
grep -E "^(OK|FAIL|MISSING) " "$LOG"
