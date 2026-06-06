#!/bin/bash
failures=0
for e in 0*_* 1*_* 2*_* 3*_*; do
  echo "=== $e ==="
  ./"$e" || failures=$((failures + 1))
  echo ""
done
echo "=== All examples completed ==="
exit $failures
