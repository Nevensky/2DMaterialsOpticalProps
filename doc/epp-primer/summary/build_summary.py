#!/usr/bin/env python3
"""Build the executive summary of the primer.

head.html (retitled "Electron–Phonon–Photon Brief"; fonts, tokens, CSS and the MathJax config with lazy
loading kept as they are) + lib.js inlined exactly as build/assemble.py does + summary.html, wrapped in
<div class="summary-page"> (no contents rail).

Usage: python3 summary/build_summary.py [out.html]   (default: summary/out/epp-brief.html)
"""
import os, re, sys

H = os.path.dirname(os.path.abspath(__file__))   # summary/
S = os.path.dirname(H)                           # doc/epp-primer/
B = os.path.join(S, 'build')
out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(H, 'out', 'epp-brief.html')

def read(p):
    with open(p, encoding='utf-8') as f:
        return f.read()

head = read(os.path.join(B, 'head.html'))
head, n = re.subn(r'<title>.*?</title>', '<title>Electron–Phonon–Photon Brief</title>', head, count=1, flags=re.S)
if n != 1:
    sys.exit('build_summary: no <title> found in build/head.html')
lib = read(os.path.join(B, 'lib.js'))
body = read(os.path.join(H, 'summary.html'))

doc = (head + '\n<script>\n' + lib + '\n</script>\n'
       + '<div class="summary-page">\n<main id="main">\n' + body + '\n</main>\n</div>\n')
os.makedirs(os.path.dirname(out), exist_ok=True)
with open(out, 'w', encoding='utf-8') as f:
    f.write(doc)

# prose word count as rendered: drop style/script blocks, display math and tags; keep "Read more" lines and add
# the vertex caption that the script injects on load (the finite-c one); inline math counts as one word
m = re.search(r"full: '(.*?)',\s*\n", body, re.S)
txt = re.sub(r'<(style|script)\b.*?</\1>', ' ', body, flags=re.S)
if m:
    txt += ' ' + re.sub(r"'\s*\+\s*link\('[^']*',\s*'([^']*)'\)\s*\+\s*'", r'\1', m.group(1))
txt = re.sub(r'\$\$.*?\$\$', ' ', txt, flags=re.S)
txt = re.sub(r'\$[^$]+\$', ' x ', txt)
txt = re.sub(r'<[^>]+>', ' ', txt)
words = sum(1 for w in txt.split() if re.search(r'[A-Za-z0-9\u0370-\u03FF]', w))
print(f'wrote {out}: {len(doc)/1024:.0f} KB, {words} words of prose (captions, injected caption and "Read more" lines included)')
