#!/bin/bash
# Whole-page check of primer v2: assemble, lint every fragment, run the page checker, dump the text.
# Usage: build/fullcheck.sh [--shots]   (writes out/epp-primer-v2.html, test/full/, out/primer-text.md)
V=$(cd "$(dirname "$0")/.." && pwd)
python3 "$V/build/assemble.py" 2>&1 | grep -iv '^wrote' | head -20
echo "== lint"; python3 "$V/build/lint.py" "$V/sections" --inventory "$V/test/inventory.md" | tail -40
echo "== check.js"; mkdir -p "$V/test/full"
node "$V/build/check.js" "$V/out/epp-primer-v2.html" "$V/test/full/check" --eager > "$V/test/full/check.json" 2>&1
python3 - "$V/test/full/check.json" <<'PY'
import json, sys
r = json.load(open(sys.argv[1]))
for e in r['errors']: print('ERROR', e[:1500])
print('math errors:', r['mathErrors'], '| overflow:', {k: (v['px'], v['offenders'][:5]) for k, v in r['overflow'].items()}, '| blank canvases:', r['blankCanvases'][:10])
print('sections:', len(r.get('sections', [])), '| height:', r.get('height'))
PY
python3 "$V/build/textdump.py" >/dev/null 2>&1 && echo "text dump: $V/out/primer-text.md ($(wc -c < "$V/out/primer-text.md") bytes)"
if [ "$1" = "--shots" ]; then node "$V/build/shoot.js" "$V/out/epp-primer-v2.html" "$V/test/full/shots" >/dev/null 2>&1; echo "shots: $(ls "$V/test/full/shots" | wc -l)"; fi
