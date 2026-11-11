#!/usr/bin/env python3
"""Notation linter for primer v2 (see NOTATION.md, section 11).

Usage: lint.py [sections_dir] [--inventory out.md]

Checks the main text of every fragment, that is everything outside the footnotes (div.notes),
the glossary box, source callouts, <span class="gl-src"> and <script>:
  - source-only symbols (Giustino/EPP notation) that belong in footnotes only;
  - italic labels and italic exponentials ({\\rm ...}, e^{i...}, ^{ion} ...);
  - the old forms of the electron linewidth (NOTATION §12.3): a bare \\Gamma with an electron subscript
    (\\Gamma_{n\\bk}, \\Gamma_k), the review's \\Gamma^{\\rm ZP} (with or without a subscript) or v1's \\Gamma(T), instead of
    \\Gamma^{\\tau}_{n\\bk}, and, in A.7 and A.8, y as the imaginary part of the pole (\\ii y, |y|, y = Im Sigma,
    y <= 0); the zone centre \\Gamma stays legal;
  - footnote markers whose target is missing, and footnotes never referenced;
  - glossary entries with a missing data-first anchor.
With --inventory it also writes the symbol inventory: every distinct symbol in reading order with
the sentence where it first appears, for the naive-reader check.
"""
import re, sys, os, glob, html

V = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
args, _skip = [], False
for _a in sys.argv[1:]:
    if _skip: _skip = False; continue
    if _a == '--inventory': _skip = True; continue
    if not _a.startswith('--'): args.append(_a)
secdir = args[0] if args else os.path.join(V, 'sections')
inv_out = sys.argv[sys.argv.index('--inventory') + 1] if '--inventory' in sys.argv else None

INC = re.compile(r'<!--\s*include:\s*([^\s>]+)\s*-->')
def expand(text, base):
    def sub(m):
        p = m.group(1) if os.path.isabs(m.group(1)) else os.path.join(base, m.group(1))
        return expand(open(p).read(), os.path.dirname(p)) if os.path.exists(p) else ''
    return INC.sub(sub, text)

def strip_blocks(s):
    """Remove the parts where source notation is allowed, keeping everything else in place."""
    s = re.sub(r'<script\b.*?</script>', ' ', s, flags=re.S)
    s = re.sub(r'<span class="gl-src">.*?</span>', ' ', s, flags=re.S)
    s = re.sub(r'<div class="notes">.*?</ol>\s*</div>', ' ', s, flags=re.S)
    s = re.sub(r'<details class="glossary-box">.*?</details>', ' ', s, flags=re.S)
    s = re.sub(r'<div class="callout source">.*?</div>', ' ', s, flags=re.S)
    s = re.sub(r'<!--\s*srcmap\s*-->.*?<!--\s*/srcmap\s*-->', ' ', s, flags=re.S)  # dictionary tables (X.2)
    return s

MATH = re.compile(r'\$\$(.+?)\$\$|\$(.+?)\$', re.S)

FORBIDDEN = [
    (r'\\Delta\s*\\tau', 'Giustino displacement Delta tau: use u_{l kappa alpha}'),
    (r'(?<![A-Za-z\\])C_\{\\kappa', 'Giustino force constants C: use Phi'),
    (r'D\^\{\\(rm|mathrm)\s*\{?dm', 'Giustino dynamical matrix D^dm: use tilde Phi'),
    (r'(?<![A-Za-z\\])N_p\b', 'N_p: use N_ell'),
    (r'(?<![A-Za-z\\])M_0\b', 'reference mass M_0 (Giustino): drop it'),
    (r'(?<![A-Za-z\\])l_\{\\bq', "Giustino's zero-point length l_q nu: write (hbar/2 M omega)^(1/2)"),
    (r'\\varepsilon_\{n', 'band energy epsilon_nk: use E_nk'),
    (r'\\Omega\^2', 'EPP dynamical matrix Omega^2: use Phi(q) or tilde Phi'),
    (r'\\Omega_\{\\rm uc\}', 'EPP cell volume Omega_uc: use \\OV'),
    (r'(?<![A-Za-z\\])Q_\{i', 'EPP displacement Q_{ia mu}: use u'),
    (r'(?<![A-Za-z\\])m_\{?\\rm\s*a|(?<![A-Za-z\\])m_a\b', 'EPP mass m_a: use M_kappa'),
    (r'_\{[^}]*\\mu[^}]*\}', 'mu as an index: Cartesian indices are alpha, beta (mu is the chemical potential)'),
    (r'\\epsilon', 'dielectric functions are \\varepsilon (the sources write \\epsilon)'),
    (r'(?<![A-Za-z\\])D_\{?11\b', 'EPP sector D_11: use D_{\\rho\\rho}'),
    (r'\\Pi_\{?(jj|\\rho\s*j|j\s*\\rho)', 'EPP correlator labels: use Pi_{alpha beta}, Pi_{rho beta}, Pi_{alpha rho}'),
    (r'\\Omega\^\{?\\(rm|mathrm)\s*\{?e|\\Omega\^\{?e\b', 'EPP phonon self-energy Omega^e: use \\Pi^{\\ph}'),
    (r'W_\{\\(rm|mathrm)\s*\{?(e|ph)', "Giustino's W_e, W_ph: use W^{\\rme}, W^{\\ph}"),
    (r'(?<![A-Za-z\\])S_\{?\\(rm|mathrm)\s*\{?uc', 'EPP cell area S_uc: use \\OS'),
]
ITALIC = [
    (r'[\^_]\{\\rm\s', 'label with {\\rm ...}: use \\mathrm{...} (or a v2 macro) so the whole word is upright'),
    (r'(?<![\\a-zA-Z])e\^\{\s*-?\s*(\\ii|i)\b', 'italic e in an exponential: use \\ee^{\\ii ...}'),
    (r'\\ee\^\{\s*-?\s*i\b', 'italic i in an exponential: use \\ii'),
    (r'[\^_]\{(ion|ph|e|para|dia|SC|KS|xc|dd|sr|lr|FM|DW|qp|NA|irr|ZPR|mob|abs|2D|en|H|el|ep|L|T|LT|TL|LO|TO|F|ph,A|ph,NA)\}', 'italic label: use \\mathrm{...} or the v2 macro'),
]

SCRIPT_RESERVED = [
    (r'D\^\{?dm|Δτ|N_p\b|M_0\b|Ω\^2|Ω_\{?uc', 'source-only symbol in a figure string'),
    (r'σ(?![_]?\{?xx|_\{?αβ)', 'sigma is reserved for the conductivity (use u_rms for an rms displacement)'),
    (r'Γ(?:_|\(|\^\{?ZP)', 'electron linewidth written as a bare Γ: write Γ^{τ} (Γ alone is the zone centre)'),
]

# The electron linewidth is \Gamma^{\tau}_{n\bk} (NOTATION §12.3); \Gamma alone is the zone centre.
GAMMA = re.compile(r'\\Gamma(?![A-Za-z])')
def _script(tex, i):
    """The ^ or _ script that starts at tex[i] (spaces skipped): (kind, body, end), or None."""
    j = i
    while j < len(tex) and tex[j].isspace(): j += 1
    if j >= len(tex) or tex[j] not in '^_': return None
    kind = tex[j]; j += 1
    while j < len(tex) and tex[j].isspace(): j += 1
    if j >= len(tex): return None
    if tex[j] == '{':
        depth = 0
        for k in range(j, len(tex)):
            depth += {'{': 1, '}': -1}.get(tex[k], 0)
            if depth == 0: return kind, tex[j + 1:k], k + 1
        return kind, tex[j + 1:], len(tex)
    m = re.match(r'\\[A-Za-z]+', tex[j:])
    return (kind, m.group(0), j + len(m.group(0))) if m else (kind, tex[j], j + 1)
def _electron_sub(sub):
    """A subscript that names an electron state: \\bk or \\mathbf k, or a plain k (\\Gamma_k, \\Gamma_{nk}), with or
    without a band index. Upright labels (\\mathrm{...}) and other commands are not read as a k."""
    if re.search(r'\\bk(?![A-Za-z])|\\mathbf\s*\{?\s*k(?![A-Za-z])', sub): return True
    s = re.sub(r'\\(?:mathrm|text|operatorname|mathsf)\s*\{[^{}]*\}|\{\\rm\s[^{}]*\}', ' ', sub)
    return 'k' in re.sub(r'\\[A-Za-z]+', ' ', s)
def bare_linewidths(tex):
    """Each \\Gamma that is an electron linewidth without the superscript \\tau: a subscript that names an electron state
    (\\Gamma_{n\\bk}, \\Gamma_{nn\\bk}, \\Gamma_{\\bk}, \\Gamma_k), the review's zero-point label with or without a subscript
    (\\Gamma^{\\mathrm{ZP}}, \\Gamma^{\\mathrm{ZP}}_{n\\bk}), or a temperature argument (\\Gamma(T), the form of primer v1).
    The zone centre (\\Gamma alone, or with any other script) and \\Gamma^{\\tau}_{n\\bk}, \\Gamma_{n\\bk}^{\\tau},
    \\Gamma^{\\tau}(T) pass."""
    out = []
    for m in GAMMA.finditer(tex):
        sup = sub = None; end = m.end()
        for _ in range(2):
            r = _script(tex, end)
            if not r: break
            kind, body, end = r
            if kind == '^': sup = body
            else: sub = body
        if sup is not None and re.search(r'\\tau(?![A-Za-z])|τ', sup): continue
        temp = re.match(r'(?:\s|\\[,;:!])*\(\s*T(?![A-Za-z])', tex[end:])
        if not ((sup is not None and 'ZP' in sup) or temp or (sub is not None and _electron_sub(sub))): continue
        out.append(tex[m.start():end + (temp.end() if temp else 0)])
    return out
# y as the imaginary part of the pole z = E^qp + i y (A.7, A.8 only; y stays a coordinate or a dummy variable elsewhere)
POLE_Y = [
    (r'(?:\\ii|\\mathrm\s*\{?\s*i\s*\}?)\s*(?:\\[,;:!]\s*)*y(?![A-Za-z])', 'i y'),
    (r'\|\s*y\s*\|', '|y|'),
    (r'(?<![A-Za-z\\])y\s*=\s*\\mathrm\s*\{?\s*Im', 'y = Im'),
    (r'(?<![A-Za-z\\])y\s*\\(?:le|leq)(?![A-Za-z])', 'y <= 0'),
]

def sentence_at(text, pos):
    a = max(text.rfind('. ', 0, pos), text.rfind('\n', 0, pos), 0)
    b = text.find('. ', pos); b = len(text) if b < 0 else b + 1
    return re.sub(r'\s+', ' ', text[a:b]).strip()[:240]

TOKEN = re.compile(r'\\(?:b[a-zA-Z]+|tPhi|hat\s*\\?[a-zA-Z]+|tilde\s*\\?[a-zA-Z]+|boldsymbol\\[a-zA-Z]+|[A-Za-z]+)(?:_\{[^{}]*\}|_\\[A-Za-z]+|_[A-Za-z0-9])?|(?<![\\A-Za-z])[A-Za-z](?:_\{[^{}]*\}|_\\[A-Za-z]+|_[A-Za-z0-9])?')
SKIP_CMDS = {'frac', 'sum', 'qquad', 'quad', 'text', 'left', 'right', 'big', 'Big', 'begin', 'end', 'cdot', 'mathrm', 'dots', 'cdots',
             'partial', 'sqrt', 'ee', 'ii', 'dd', 'le', 'ge', 'lt', 'gt', 'neq', 'simeq', 'equiv', 'pm', 'to', 'infty', 'Longleftrightarrow',
             'underbrace', 'displaystyle', 'tag', 'pmatrix', 'gathered', 'aligned', 'substack', 'smallmatrix', 'bmatrix', 'langle', 'rangle',
             'cos', 'sin', 'coth', 'det', 'Re', 'Im', 'exp', 'ln', 'log', 'delta', 'dagger', 'ell', 'kappa', 'alpha', 'beta', 'nu', 'mu', 'gamma',
             'over', 'mathbf', 'boldsymbol', 'operatorname', 'sinh', 'cosh', 'tanh', 'tan', 'sim', 'leftrightarrow', 'rightarrow', 'leftarrow', 'iint', 'iiint', 'oint', 'nabla', 'max', 'min', 'lim', 'sup', 'inf', 'mathcal', 'mathsf', 'bar', 'overline', 'widehat', 'widetilde', 'quad', 'lbrace', 'rbrace', 'Bigg', 'bigl', 'bigr', 'Bigl', 'Bigr', 'cot', 'csc', 'sec', 'arctan', 'arcsin', 'arccos', 'mapsto', 'Rightarrow', 'Leftarrow', 'perp', 'parallel', 'cap', 'cup', 'forall', 'exists', 'colon', 'phantom', 'hphantom', 'vphantom', 'qquad', 'textstyle', 'scriptstyle', 'limits', 'nolimits', 'stackrel', 'overset', 'underset', 'times', 'approx', 'hbar', 'pi', 'kB', 'hat', 'tilde', 'circ', 'prime', 'in', 'notin',
             'tfrac', 'mathsf', 'ddot', 'dot', 'll', 'gg', 'mp', 'propto', 'arg', 'Delta', 'int', 'lesssim', 'gtrsim', 'quad', 'big', 'Big', 'bigg', 'vert', 'mid'}

def main():
    files = sorted(glob.glob(os.path.join(secdir, '*.html'))) + [os.path.join(V, 'build', 'hero.html')]
    problems, inventory, seen = [], [], set()
    for f in files:
        raw = expand(open(f).read(), os.path.dirname(f))
        name = os.path.basename(f)
        # footnote markers and targets
        refs = re.findall(r'<sup class="fn"><a href="#([^"]+)"', raw)
        targets = re.findall(r'<li id="(fn-[^"]+)"', raw)
        for r in refs:
            if r not in targets: problems.append((name, 'footnote', f'marker points to missing note #{r}'))
        for t in targets:
            if t not in refs: problems.append((name, 'footnote', f'note #{t} is never referenced'))
        order = [t for t in targets if t in refs]
        firsts = []
        for r in refs:
            if r not in firsts: firsts.append(r)
        if order != firsts: problems.append((name, 'footnote', 'notes are not listed in order of first reference: ' + ', '.join(firsts)))
        ids = set(re.findall(r'\bid="([^"]+)"', raw))
        for gid, first in re.findall(r'<dt id="(gl-[^"]+)" data-first="([^"]*)"', raw):
            if first and first not in ids: problems.append((name, 'glossary', f'{gid}: data-first="{first}" has no element with that id'))
        for sc in re.findall(r'<script\b[^>]*>(.*?)</script>', raw, flags=re.S):
            for lit in re.findall(r"'(?:[^'\\\n]|\\.)*'|\"(?:[^\"\\\n]|\\.)*\"|`[^`]*`", sc):
                for pat, msg in SCRIPT_RESERVED:
                    if re.search(pat, lit): problems.append((name, 'script string', f'{msg}  [{lit[:60]}]'))
        main_text = strip_blocks(raw)
        plain = html.unescape(re.sub(r'<[^>]+>', ' ', main_text))
        for m in MATH.finditer(main_text):
            tex = m.group(1) or m.group(2)
            line = main_text.count('\n', 0, m.start()) + 1
            for pat, msg in FORBIDDEN:
                if re.search(pat, tex): problems.append((name, f'line {line}', f'source notation in main text: {msg}  [{tex.strip()[:70]}]'))
            for pat, msg in ITALIC:
                if re.search(pat, tex): problems.append((name, f'line {line}', f'{msg}  [{tex.strip()[:70]}]'))
            for g in bare_linewidths(tex):
                problems.append((name, f'line {line}', f'electron linewidth written as a bare \\Gamma: use \\Gamma^{{\\tau}}_{{n\\bk}}, \\Gamma^{{\\tau}}_{{n\\bk}}(T), \\Gamma^{{\\tau}}_{{n\\bk}}(T=0) (\\Gamma alone is the zone centre)  [{g[:70]}]'))
            if name.startswith(('A07', 'A08')):
                for pat, form in POLE_Y:
                    if re.search(pat, tex): problems.append((name, f'line {line}', f'y as the imaginary part of the pole ({form}): write z = E^qp - i\\Gamma^{{\\tau}}_{{n\\bk}} with \\Gamma^{{\\tau}}_{{n\\bk}} >= 0  [{tex.strip()[:70]}]'))
        if inv_out:
            for m in MATH.finditer(main_text):
                tex = m.group(1) or m.group(2)
                tex = re.sub(r'\\(begin|end|text|mathrm|tag|operatorname)\{[^{}]*\}', ' ', tex)
                for t in TOKEN.finditer(tex):
                    tok = re.sub(r'\s+', '', t.group(0))
                    tok = re.sub(r'^\\hat\s*', '', tok); base = re.sub(r'_.*', '', tok).lstrip('\\')
                    if base in SKIP_CMDS or tok in seen: continue
                    seen.add(tok)
                    ctx_plain = html.unescape(re.sub(r'<[^>]+>', ' ', main_text[max(0, m.start() - 300):m.end() + 200]))
                    inventory.append((name, tok, sentence_at(ctx_plain, min(len(ctx_plain) - 1, 300))))
    for p in problems: print(' | '.join(p))
    print(f'{len(problems)} problem(s) in {len(files)} fragment(s)')
    if inv_out:
        gl_tex = ' '.join(t for f in files for t in re.findall(r'<dt id="gl-[^"]*"[^>]*>(.*?)</dt>', expand(open(f).read(), os.path.dirname(f)), flags=re.S))
        gl_tex = re.sub(r'\s+', '', gl_tex)
        ALLOW = {'x', 'y', 'z', 't', 'f', 'i', 'A', 'B', 'C', 'X', 'P', 'L', 'T_x', 'O'}
        uncovered = []
        for name, tok, ctx in inventory:
            base = re.sub(r'_.*', '', tok)
            if tok in gl_tex or base in gl_tex or base.lstrip('\\') in ALLOW: continue
            uncovered.append((name, tok, ctx))
        with open(inv_out, 'w') as fh:
            fh.write('# Symbol inventory (first appearance in reading order)\n\n')
            for name, tok, ctx in inventory: fh.write(f'- `{tok}` ({name}): {ctx}\n')
            fh.write('\n# Warnings: symbols in the main text with no glossary entry\n\n')
            for name, tok, ctx in uncovered: fh.write(f'- `{tok}` ({name}): {ctx[:160]}\n')
        print(f'{len(uncovered)} symbol(s) without a glossary entry (warnings, see {inv_out})')
        print(f'wrote {inv_out}: {len(inventory)} symbols')
    return 1 if problems else 0

sys.exit(main())
