#!/bin/bash
# Test one v2 section in context. Usage: build/sectest.sh A02 [--shots]
# Assembles C00 + A01 (finished) + this section + empty stubs for every other section, then runs
# lint (this section's files) and check.js. Writes test/sec-<ID>/{page.html,check.json,shots/}.
# Dead links into other sections' glossary, equation or figure ids are listed separately as EXTERNAL:
# they resolve only in the full page, so each must be an id from NOTATION.md section 12 or a brief's exports.
V=$(cd "$(dirname "$0")/.." && pwd)
ID=$1; shift
F=$(ls "$V/sections/$ID"-*.html 2>/dev/null | head -1)
[ -z "$F" ] && { echo "no section file sections/$ID-*.html"; exit 1; }
D=$V/test/sec-$ID; rm -rf "$D"; mkdir -p "$D/frag"
for s in C00-conventions A01-lattice; do echo "<!-- include: $V/sections/$s.html -->" > "$D/frag/$s.html"; done
echo "<!-- include: $F -->" > "$D/frag/$(basename "$F")"
while read -r sid num key; do
  [ "$key" = "$ID" ] && continue; [ "$key" = "A01" ] && continue
  echo "<section class=\"chapter\" id=\"$sid\" data-num=\"$num\" data-title=\"(not in this test)\"><header class=\"ch-head\"><h2>$num (stub)</h2></header></section>" > "$D/frag/$key-stub.html"
done <<'LIST'
a-coupling A.2 A02
a-dfpt A.3 A03
a-polar A.4 A04
a-hedin-baym A.5 A05
a-phonon-se A.6 A06
a-electron-se A.7 A07
a-temperature A.8 A08
a-applications A.9 A09
x-light-cone X.1 X01
x-dictionary X.2 X02
b-foundations B.1 B01
b-lagrangian B.2 B02
b-dielectric B.3 B03
b-fields B.4 B04
b-photon-se B.5 B05
b-propagator B.6 B06
b-phonon-se B.7 B07
b-frohlich B.8 B08
b-electrostatic B.9 B09
b-2d B.10 B10
LIST
python3 "$V/build/assemble.py" "$D/page.html" "$D/frag" 2>&1 | grep -i warn
echo "== lint ($ID)"; L=$(mktemp -d); echo "<!-- include: $F -->" > "$L/$(basename "$F")"
python3 "$V/build/lint.py" "$L" | grep -v 'hero.html' ; rm -rf "$L"
echo "== check.js"
node "$V/build/check.js" "$D/page.html" "$D/check" --eager > "$D/check.json" 2>&1
PFX=$(echo "$ID" | sed -E 's/^([A-Z])0*([0-9]+)$/\1\2/' | tr 'A-Z' 'a-z')
python3 - "$D/check.json" "$PFX" <<'PY'
import json, sys, re
r = json.load(open(sys.argv[1])); pfx = sys.argv[2]
own = re.compile(rf'^(viz-{pfx}-|gl-{pfx}-|fn-{pfx}-|{pfx}-)')
ext = set()
for e in r['errors']:
    m = re.match(r'\[(\w+)\] dead in-page links: (.*)', e)
    if m:
        ids = [x.strip() for x in m.group(2).split(',')]
        bad = [i for i in ids if own.match(i) or not re.match(r'^(viz-|gl-|[abx]\d+-eq-|fn-)', i)]
        ext.update(i for i in ids if i not in bad)
        if bad: print('ERROR', f'[{m.group(1)}] dead links inside this section or to unknown ids:', ', '.join(bad))
    else: print('ERROR', e)
if ext: print('EXTERNAL (must be registry/export ids):', ', '.join(sorted(ext)))
for k, v in (r.get('inlineClip') or {}).items():
    print(f'CLIPPED INLINE MATH [{k}] (wider than the text column; shorten it or move it to a display equation): {len(v)}')
    for x in v[:12]: print('   ', x)
print('math errors:', r['mathErrors'], '| overflow px:', {k: v['px'] for k, v in r['overflow'].items()}, '| blank canvases:', r['blankCanvases'])
PY
if [ "$1" = "--shots" ]; then node "$V/build/shoot.js" "$D/page.html" "$D/shots" >/dev/null 2>&1; ls "$D/shots" | grep -v '^index' | grep -i "$PFX\|$(echo $ID | tr A-Z a-z)" | head -40; fi
echo "page: $D/page.html"
