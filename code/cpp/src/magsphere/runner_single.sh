#!/usr/bin/env bash
set -euo pipefail

exe="./magsphere.out"
base_input="input.base"
outdir="runs"
mkdir -p "$outdir"

# Sweep settings
start=0.0
end=1.0
n=100

mapfile -t mags < <(python3 - "$start" "$end" "$n" <<'PY'
import sys
start = float(sys.argv[1])
end = float(sys.argv[2])
n = int(sys.argv[3])

if n == 1:
    print(start)
else:
    step = (end - start) / (n - 1)
    for i in range(n):
        print(start + i * step)
PY
)

i=0
for mx in "${mags[@]}"; do
  tmp_input="$(mktemp)"

  awk -v mx="$mx" '
    NR==9 {$1 = mx}
    {print}
  ' "$base_input" > "$tmp_input"

  "$exe" < "$tmp_input"

  printf -v idx "%03d" "$i"   # 000, 001, ..., 099
  mv -f data.csv "$outdir/data_mz_${idx}.csv"

  echo "$idx $mx" >> "$outdir/index_map.txt"

  rm -f "$tmp_input"
  ((++i))
done