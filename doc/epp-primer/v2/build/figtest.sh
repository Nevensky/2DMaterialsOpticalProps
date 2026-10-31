#!/bin/bash
# Test one figure file in isolation. Usage: build/figtest.sh sections/A01/fig-xxx.html
# Writes test/<name>/page.html, runs check.js (errors, TeX, overflow, blank canvases) and shoots the figure.
V=$(cd "$(dirname "$0")/.." && pwd)
F=$(cd "$(dirname "$1")" && pwd)/$(basename "$1")
N=$(basename "$1" .html)
D=$V/test/$N; rm -rf "$D"; mkdir -p "$D/frag"
cat > "$D/frag/A01-test.html" <<HTML
<section class="chapter" id="a-lattice" data-num="A.1" data-title="Figure test">
<header class="ch-head"><h2>Figure test: $N</h2></header>
<!-- include: $F -->
</section>
HTML
python3 "$V/build/assemble.py" "$D/page.html" "$D/frag" >/dev/null || exit 1
node "$V/build/check.js" "$D/page.html" "$D/check" --eager > "$D/check.json" 2>&1; cat "$D/check.json" | head -40
node "$V/build/shoot.js" "$D/page.html" "$D/shots" >/dev/null 2>&1 && ls "$D/shots" | head
