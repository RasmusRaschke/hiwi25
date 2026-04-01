#!/usr/bin/env bash
set -euo pipefail

exe="./magsphere.out"
base_input="input.base"
outdir="runs"
mkdir -p "$outdir"

# Generate values from start to end with a fixed step.
vals() {
  awk -v start="$1" -v end="$2" -v step="$3" '
    BEGIN {
      for (x = start; x <= end + step/1000; x += step) {
        printf "%.10g\n", x
      }
    }'
}

tagify() {
  printf '%s' "$1" | sed 's/-/m/g; s/\./p/g'
}

run_case() {
  local mx="$1"
  local my="$2"
  local tmp_input
  tmp_input="$(mktemp)"

  # Replace line 7 and line 8
  awk -v mx="$mx" -v my="$my" '
    NR==7 {$1 = mx}
    NR==8 {$1 = my}
    {print}
  ' "$base_input" > "$tmp_input"

  "$exe" < "$tmp_input"

  mv -f data.csv "$outdir/data_mx_$(tagify "$mx")_my_$(tagify "$my").csv"
  rm -f "$tmp_input"
}

# Example:
# x-moment from 0 to 1 in steps of 0.01
# y-moment from -0.5 to 0.5 in steps of 0.1
mapfile -t mx_vals < <(vals -1 1 0.1)
mapfile -t my_vals < <(vals -1 1 0.1)

for mx in "${mx_vals[@]}"; do
  for my in "${my_vals[@]}"; do
    run_case "$mx" "$my"
  done
done
