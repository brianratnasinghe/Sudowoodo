#!/bin/bash
set -euo pipefail

for dir in pp_eps_*/; do
    echo "Running in $dir"
    (cd "$dir" && bash run.sh)
done
