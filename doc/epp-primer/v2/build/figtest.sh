#!/bin/bash
# Test one figure file in isolation. Usage: build/figtest.sh sections/B05/fig-xxx.html
# The section is read from the folder name (A02 -> A.2, X01 -> X.1, B10 -> B.10), so figure numbers come out right.
# Writes test/<name>/page.html, runs check.js (errors, TeX, overflow, blank canvases, dead links) and shoots the figure.
V=$(cd "$(dirname "$0")/.." && pwd)
F=$(cd "$(dirname "$1")" && pwd)/$(basename "$1")
N=$(basename "$1" .html)
SEC=$(basename "$(dirname "$F")")
if [[ $SEC =~ ^([A-Z])0*([0-9]+)$ ]]; then NUM="${BASH_REMATCH[1]}.${BASH_REMATCH[2]}"; else SEC=A01; NUM=A.1; fi
D=$V/test/$SEC-$N; rm -rf "$D"; mkdir -p "$D/frag"
cat > "$D/frag/$SEC-test.html" <<HTML
<section class="chapter" id="figtest" data-num="$NUM" data-title="Figure test">
<header class="ch-head"><h2>Figure test: $SEC/$N</h2></header>
<!-- include: $F -->
</section>
HTML
python3 "$V/build/assemble.py" "$D/page.html" "$D/frag" >/dev/null || exit 1
node "$V/build/check.js" "$D/page.html" "$D/check" --eager > "$D/check.json" 2>&1; head -40 "$D/check.json"
node "$V/build/shoot.js" "$D/page.html" "$D/shots" >/dev/null 2>&1 && ls "$D/shots" | head
echo "page: $D/page.html  shots: $D/shots"
