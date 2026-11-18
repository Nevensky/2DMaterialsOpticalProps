#!/usr/bin/env python3
"""Build the research roadmap of the open EPP topics.

Primer v2's v2/build/head.html (retitled "EPP Research Roadmap"; fonts, tokens, CSS and the MathJax config with
the v2 macros and lazy loading kept as they are) + v2/build/lib.js inlined exactly as v2/build/assemble.py does
+ roadmap.html, wrapped in <div class="roadmap-page"> (no contents rail). The page uses the notation of primer v2
(v2/NOTATION.md); rules the v2 head lacks or that the page overrides live in roadmap.html's own <style> block.

Typesetting: typeset() adds no-break spaces and word joiners to the prose at build time; roadmap.html stays plain.

Usage: python3 open/build_roadmap.py [out.html]   (default: open/out/epp-roadmap.html)
"""
import os, re, sys

H = os.path.dirname(os.path.abspath(__file__))   # open/
S = os.path.dirname(H)                           # doc/epp-primer/
B = os.path.join(S, 'v2', 'build')               # primer v2: head.html (v2 macros, tokens, CSS) and lib.js
out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(H, 'out', 'epp-roadmap.html')

def read(p):
    with open(p, encoding='utf-8') as f:
        return f.read()

head = read(os.path.join(B, 'head.html'))
head, n = re.subn(r'<title>.*?</title>', '<title>EPP Research Roadmap</title>', head, count=1, flags=re.S)
if n != 1:
    sys.exit('build_roadmap: no <title> found in v2/build/head.html')
# The v2 head always typesets the primer's glossary panel and popovers (lazyAlwaysTypeset). This page has neither, and
# MathJax 3.2.2's lazy extension throws on an empty container list, so no formula would be typeset: drop the option.
head, n = re.subn(r",\s*lazyAlwaysTypeset:\s*\[[^\]]*\]", '', head, count=1)
if n != 1:
    sys.exit('build_roadmap: lazyAlwaysTypeset not found in v2/build/head.html; check the MathJax options')
lib = read(os.path.join(B, 'lib.js'))
body = read(os.path.join(H, 'roadmap.html'))

# Line-break hygiene, applied to prose only (style/script blocks, tags and TeX are stashed first):
# a no-break space between a number and its unit (a nowrap span when the number is inline math), a word joiner after the en dash of
# numeric and equation ranges and of LO–TO / L–T, nowrap spans around math ranges, "$x$-word", "word-$x$" and a formula with the
# punctuation that closes it, and a no-break space in "topic 7" and "step 2".
UNIT = r'(?:eV|meV|μeV|Å[⁻¹²³]*|K|cm[⁻¹²³]+(?:/Vs)?|nm|μm|C/m²|amu|u|weeks?|days?|months?)(?!\w)'
WJ = '\u2060'
def typeset(html):
    keep = []
    def stash(m):
        t = m.group(0)
        keep.append(t)
        return ('\x02%d\x03' if t.startswith('$') else '\x00%d\x01') % (len(keep) - 1)
    s = re.sub(r'<(style|script)\b.*?</\1>|<!--.*?-->|<[^>]+>|\$\$.*?\$\$|\$[^$]+\$', stash, html, flags=re.S)
    # inline math is an inline-block, so a no-break space beside it does not hold: wrap math + unit instead
    s = re.sub(r'((?:\x02\d+\x03–)?\x02\d+\x03) (' + UNIT + ')', '<span class="nw">\\1\u00a0\\2</span>', s)
    s = re.sub(r'(\d) (?=' + UNIT + ')', '\\1\u00a0', s)
    s = s.replace('eV Å', 'eV\u00a0Å')
    s = re.sub(r'(\d)–(?=\d|[EG]\.)', '\\1–' + WJ, s)
    s = re.sub(r'(?<!>)\x02\d+\x03–\x02\d+\x03', lambda m: '<span class="nw">' + m.group(0) + '</span>', s)
    s = s.replace('LO–TO', 'LO–' + WJ + 'TO')
    s = re.sub(r'\bL–T\b', 'L–' + WJ + 'T', s)
    s = re.sub(r'“?\x02\d+\x03-[A-Za-z]+|[A-Za-z]+-\x02\d+\x03', lambda m: '<span class="nw">' + m.group(0) + '</span>', s)
    # a formula keeps the bracket or punctuation that closes it, and "topic 7", "steps 1" stay on one line
    s = re.sub(r'\(?\x02\d+\x03(?:\)[.,;:]?|[.,;:])', lambda m: '<span class="nw">' + m.group(0) + '</span>', s)
    s = re.sub(r'\b(topics?|steps?) (?=\d)', '\\1\u00a0', s)
    return re.sub(r'[\x00\x02](\d+)[\x01\x03]', lambda m: keep[int(m.group(1))], s)

doc = (head + '\n<script>\n' + lib + '\n</script>\n'
       + '<div class="roadmap-page">\n<main id="main">\n' + typeset(body) + '\n</main>\n</div>\n')
os.makedirs(os.path.dirname(out), exist_ok=True)
with open(out, 'w', encoding='utf-8') as f:
    f.write(doc)

# prose word count as rendered: drop style/script blocks, display math and tags (SVG text included);
# inline math counts as one word
def words(html):
    txt = re.sub(r'<(style|script)\b.*?</\1>', ' ', html, flags=re.S)
    txt = re.sub(r'\$\$.*?\$\$', ' ', txt, flags=re.S)
    txt = re.sub(r'\$[^$]+\$', ' x ', txt)
    txt = re.sub(r'<[^>]+>', ' ', txt)
    return sum(1 for w in txt.split() if re.search(r'[A-Za-z0-9Ͱ-Ͽ]', w))

total = words(body)
print(f'wrote {out}: {len(doc)/1024:.0f} KB, {total} words of prose (tables, captions and reference lists included)')
for m in re.finditer(r'<(section|header|footer)\b[^>]*\bid="([^"]+)"[^>]*>(.*?)</\1>', body, re.S):
    print(f'  {m.group(2):<18} {words(m.group(3)):>5} words')
