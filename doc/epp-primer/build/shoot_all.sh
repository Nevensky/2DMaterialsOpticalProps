#!/bin/bash
# Shoot figures section by section from isolated assemblies. Usage: shoot_all.sh <outdir> [fragment files...]
S=$(cd "$(dirname "$0")/.." && pwd)
OUT=$1; shift
mkdir -p $OUT
FILES=${@:-$(ls $S/sections/*.html)}
one() {
  f=$1; b=$(basename $f .html); d=$S/test/vis/$b; mkdir -p $d/frag; cp $f $d/frag/
  python3 $S/build/assemble.py $d/page.html $d/frag >/dev/null
  node $S/build/shoot.js $d/page.html $2/$b >/dev/null 2>&1 && echo "ok $b" || echo "FAIL $b"
}
export -f one; export S
printf "%s\n" $FILES | xargs -P 3 -I{} bash -c "one {} $OUT"
