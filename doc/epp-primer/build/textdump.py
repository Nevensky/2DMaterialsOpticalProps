#!/usr/bin/env python3
"""Text-only dump of the assembled primer (TeX kept, scripts dropped) for reviewers."""
import re, sys, html, glob, os
S = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ORDER = {'A': 0, 'X': 1, 'B': 2, 'Z': 3}
files = sorted(glob.glob(S + '/sections/*.html'), key=lambda f: (ORDER.get(os.path.basename(f)[0], 9), os.path.basename(f)))
out = []
hero = open(S + '/build/hero.html').read()
for name, s in [('hero.html (build/hero.html)', hero)] + [(os.path.basename(f), open(f).read()) for f in files]:
    s = re.sub(r'<script.*?</script>', '', s, flags=re.S)
    s = re.sub(r'<style.*?</style>', '', s, flags=re.S)
    s = re.sub(r'<h([1-4])[^>]*>', lambda m: '\n' + '#' * (int(m.group(1)) + 1) + ' ', s)
    s = re.sub(r'</h[1-4]>', '\n', s)
    s = re.sub(r'<(p|li|figcaption|summary|tr|div)[^>]*>', '\n', s)
    s = re.sub(r'<td[^>]*>|<th[^>]*>', ' | ', s)
    s = re.sub(r'<[^>]+>', '', s)
    s = html.unescape(s)
    s = re.sub(r'\n\s*\n+', '\n\n', s)
    out.append(f'\n\n======== FILE: {name} ========\n' + s.strip())
open(S + '/out/primer-text.md', 'w').write(''.join(out))
print(len(''.join(out).split()), 'words')
