#!/usr/bin/env python3
"""Text-only dump of the assembled v2 page (TeX kept, scripts dropped) for reviewers: writes out/primer-text.md."""
import re, sys, html, glob, os
V = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(V, 'build'))
INC = re.compile(r'<!--\s*include:\s*([^\s>]+)\s*-->')
def expand(text, base):
    def sub(m):
        p = m.group(1) if os.path.isabs(m.group(1)) else os.path.join(base, m.group(1))
        return expand(open(p).read(), os.path.dirname(p)) if os.path.exists(p) else '[figure pending: ' + os.path.basename(p) + ']'
    return INC.sub(sub, text)
ORDER = {'C': -1, 'A': 0, 'X': 1, 'B': 2, 'Z': 3}
files = sorted(glob.glob(V + '/sections/*.html'), key=lambda f: (ORDER.get(os.path.basename(f)[0], 9), os.path.basename(f)))
out = []
for name, s in [('hero', open(V + '/build/hero.html').read())] + [(os.path.basename(f), expand(open(f).read(), os.path.dirname(f))) for f in files]:
    s = re.sub(r'<script.*?</script>', '', s, flags=re.S)
    s = re.sub(r'<style.*?</style>', '', s, flags=re.S)
    s = re.sub(r'<h([1-4])[^>]*>', lambda m: '\n' + '#' * (int(m.group(1)) + 1) + ' ', s)
    s = re.sub(r'</h[1-4]>', '\n', s)
    s = re.sub(r'<sup class="fn"><a href="#([^"]+)">\d+</a></sup>', r'[^\1]', s)
    s = re.sub(r'<li id="(fn-[^"]+)">', r'\n[^\1]: ', s)
    s = re.sub(r'<summary>', '\n[box] ', s)
    s = re.sub(r'<span class="why">', '\n  WHY: ', s)
    s = re.sub(r'<dt[^>]*>', '\n  - ', s); s = re.sub(r'<dd>', ' : ', s)
    s = re.sub(r'<(p|li|figcaption|tr|div|ol|ul|dl|details)[^>]*>', '\n', s)
    s = re.sub(r'<td[^>]*>|<th[^>]*>', ' | ', s)
    s = re.sub(r'<[^>]+>', '', s)
    s = html.unescape(s)
    s = re.sub(r'[ \t]+\n', '\n', s); s = re.sub(r'\n\s*\n+', '\n\n', s)
    out.append(f'\n\n======== {name} ========\n' + s.strip())
os.makedirs(V + '/out', exist_ok=True)
open(V + '/out/primer-text.md', 'w').write(''.join(out))
print(len(''.join(out).split()), 'words ->', V + '/out/primer-text.md')
