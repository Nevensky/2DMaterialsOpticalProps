#!/usr/bin/env python3
"""Assemble primer v2: head + lib + contents rail + hero + section fragments + foot.

Usage: assemble.py [out.html] [sections_dir]
Fragments are the *.html files directly inside sections_dir, ordered C (conventions), A, X, B, Z.
A fragment may pull in other files with <!-- include: relative/or/absolute/path.html -->.
"""
import re, sys, glob, os, html
B = os.path.dirname(os.path.abspath(__file__))
S = os.path.dirname(B)
out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(S, 'out', 'epp-primer-v2.html')
secdir = sys.argv[2] if len(sys.argv) > 2 else os.path.join(S, 'sections')
rd = lambda p: open(p, encoding='utf-8').read()
head, lib, lib2 = rd(os.path.join(B, 'head.html')), rd(os.path.join(B, 'lib.js')), ''
hero = rd(os.path.join(B, 'hero.html')) if os.path.exists(os.path.join(B, 'hero.html')) else ''
foot = rd(os.path.join(B, 'foot.html')) if os.path.exists(os.path.join(B, 'foot.html')) else ''
INC = re.compile(r'<!--\s*include:\s*([^\s>]+)\s*-->')
def expand(text, base, depth=0):
    if depth > 5: raise RuntimeError('include depth')
    def sub(m):
        p = m.group(1) if os.path.isabs(m.group(1)) else os.path.join(base, m.group(1))
        if not os.path.exists(p):
            print(f'warning: include {os.path.relpath(p, S)} not found; placeholder inserted', file=sys.stderr)
            return f'<figure class="viz"><p class="viz-title">Figure pending: {html.escape(os.path.basename(p))}</p></figure>'
        return f'<!-- from {os.path.relpath(p, S)} -->\n' + expand(rd(p), os.path.dirname(p), depth + 1)
    return INC.sub(sub, text)
ORDER = {'C': -1, 'A': 0, 'X': 1, 'B': 2, 'Z': 3}
files = sorted(glob.glob(os.path.join(secdir, '*.html')), key=lambda f: (ORDER.get(os.path.basename(f)[0], 9), os.path.basename(f)))
body, toc = [], []
for f in files:
    s = expand(rd(f), os.path.dirname(f))
    body.append(f'\n<!-- ===== {os.path.basename(f)} ===== -->\n' + s)
    for m in re.finditer(r'<(div|section)\b([^>]*)>', s):
        attrs = dict(re.findall(r'([\w-]+)="([^"]*)"', m.group(2)))
        cls = attrs.get('class', '').split()
        if 'part-head' in cls: toc.append(('part', attrs.get('id', ''), '', attrs.get('data-part-title', '')))
        elif 'chapter' in cls: toc.append(('ch', attrs.get('id', ''), attrs.get('data-num', ''), attrs.get('data-title', '')))
nav = ['<nav class="toc" aria-label="Contents"><div class="toc-inner"><details open id="toc-details"><summary>Contents</summary>']
open_ol = False
for kind, id_, num, title in toc:
    if kind == 'part':
        if open_ol: nav.append('</ol></div>')
        nav.append(f'<div class="toc-part"><span><a href="#{id_}" style="display:inline;padding:0;margin:0;color:inherit">{html.escape(title)}</a></span><ol>')
        open_ol = True
    else:
        if not open_ol: nav.append('<div class="toc-part"><ol>'); open_ol = True
        nav.append(f'<li><a href="#{id_}"><span class="n">{html.escape(num)}</span><span>{html.escape(title)}</span></a></li>')
if open_ol: nav.append('</ol></div>')
nav.append('</details></div></nav>')
doc = (head + '\n<script>\n' + lib + '\n</script>\n<div class="page">\n' + '\n'.join(nav) + '\n<main id="main">\n'
       + hero + ''.join(body) + '\n' + foot + '\n</main>\n</div>\n')
os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
open(out, 'w', encoding='utf-8').write(doc)
print(f'wrote {out}: {len(doc)/1024:.0f} KB, {len(files)} fragments, {sum(1 for t in toc if t[0]=="ch")} chapters')
