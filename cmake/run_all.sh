#!/bin/bash
set -e
for e in 0*_* 1*_* 2*_* 3*_*; do
  echo "=== $e ==="
  ./"$e"
  echo ""
done
echo "=== All examples completed ==="
