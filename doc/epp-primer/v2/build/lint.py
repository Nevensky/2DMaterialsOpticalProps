#!/usr/bin/env python3
"""Notation linter for primer v2 (see NOTATION.md, section 11).

Usage: lint.py [sections_dir] [--inventory out.md]

Checks the main text of every fragment, that is everything outside the footnotes (div.notes),
the glossary box, source callouts, <span class="gl-src"> and <script>:
  - source-only symbols (Giustino/EPP notation) that belong in footnotes only;
  - italic labels and italic exponentials ({\\rm ...}, e^{i...}, ^{ion} ...);
  - footnote markers whose target is missing, and footnotes never referenced;
  - glossary entries with a missing data-first anchor.
With --inventory it also writes the symbol inventory: every distinct symbol in reading order with
the sentence where it first appears, for the naive-reader check.
"""
import re, sys, os, glob, html

V = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
args = [a for a in sys.argv[1:] if not a.startswith('--')]
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
]
ITALIC = [
    (r'[\^_]\{\\rm\s', 'label with {\\rm ...}: use \\mathrm{...} (or a v2 macro) so the whole word is upright'),
    (r'(?<![\\a-zA-Z])e\^\{\s*-?\s*(\\ii|i)\b', 'italic e in an exponential: use \\ee^{\\ii ...}'),
    (r'\\ee\^\{\s*-?\s*i\b', 'italic i in an exponential: use \\ii'),
    (r'\^\{(ion|ph|e|para|dia|SC|KS|xc|dd|sr)\}', 'italic label superscript: use \\mathrm{...} or the v2 macro'),
]

SCRIPT_RESERVED = [
    (r'D\^\{?dm|Δτ|N_p\b|M_0\b|Ω\^2|Ω_\{?uc', 'source-only symbol in a figure string'),
    (r'σ(?![_]?\{?xx|_\{?αβ)', 'sigma is reserved for the conductivity (use u_rms for an rms displacement)'),
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
             'over', 'mathbf', 'boldsymbol', 'operatorname', 'times', 'approx', 'hbar', 'pi', 'kB', 'hat', 'tilde', 'circ', 'prime', 'in', 'notin',
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
