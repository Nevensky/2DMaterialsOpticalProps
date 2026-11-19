#!/usr/bin/env python3
"""Assemble the primer: head + lib + toc + hero + section fragments (sorted by filename) + footer."""
import re, sys, glob, os, html
B = os.path.dirname(os.path.abspath(__file__))
S = os.path.dirname(B)
out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(S, 'out', 'epp-primer.html')
secdir = sys.argv[2] if len(sys.argv) > 2 else os.path.join(S, 'sections')
head = open(os.path.join(B, 'head.html')).read()
lib = open(os.path.join(B, 'lib.js')).read()
hero = open(os.path.join(B, 'hero.html')).read() if os.path.exists(os.path.join(B, 'hero.html')) else ''
foot = open(os.path.join(B, 'foot.html')).read() if os.path.exists(os.path.join(B, 'foot.html')) else ''
ORDER = {'A': 0, 'X': 1, 'B': 2, 'Z': 3}  # Part A, Bridge, Part B, closing
files = sorted(glob.glob(os.path.join(secdir, '*.html')), key=lambda f: (ORDER.get(os.path.basename(f)[0], 9), os.path.basename(f)))
body = []
toc = []  # (kind, id, num, title)
for f in files:
    s = open(f).read()
    body.append(f'\n<!-- ===== {os.path.basename(f)} ===== -->\n' + s)
    for m in re.finditer(r'<(div|section)\b([^>]*)>', s):
        attrs = dict(re.findall(r'([\w-]+)="([^"]*)"', m.group(2)))
        cls = attrs.get('class', '')
        if 'part-head' in cls.split():
            toc.append(('part', attrs.get('id', ''), '', attrs.get('data-part-title', '')))
        elif 'chapter' in cls.split():
            toc.append(('ch', attrs.get('id', ''), attrs.get('data-num', ''), attrs.get('data-title', '')))
nav = ['<nav class="toc" aria-label="Contents"><div class="toc-inner"><details open id="toc-details"><summary>Contents</summary>']
open_ol = False
for kind, id_, num, title in toc:
    if kind == 'part':
        if open_ol: nav.append('</ol></div>')
        nav.append(f'<div class="toc-part"><span><a href="#{id_}" style="display:inline;padding:0;margin:0;color:inherit">{html.escape(title)}</a></span><ol>')
        open_ol = True
    else:
        if not open_ol:
            nav.append('<div class="toc-part"><ol>'); open_ol = True
        nav.append(f'<li><a href="#{id_}"><span class="n">{html.escape(num)}</span><span>{html.escape(title)}</span></a></li>')
if open_ol: nav.append('</ol></div>')
nav.append('</details></div></nav>')
doc = (head + '\n<script>\n' + lib + '\n</script>\n' + '<div class="page">\n' + '\n'.join(nav) + '\n<main id="main">\n' + hero + ''.join(body) + '\n' + foot + '\n</main>\n</div>\n')
os.makedirs(os.path.dirname(out), exist_ok=True)
open(out, 'w').write(doc)
print(f'wrote {out}: {len(doc)/1024:.0f} KB, {len(files)} fragments, {sum(1 for t in toc if t[0]=="ch")} chapters')
