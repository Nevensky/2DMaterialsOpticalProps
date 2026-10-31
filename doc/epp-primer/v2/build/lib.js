/* EPP helper library: theme-aware canvas plots, controls, animation loop,
   Feynman-diagram SVG builder, small numerics. Exposed as window.EPP. */
(function () {
  'use strict';
  const EPP = {};
  const root = document.documentElement;

  // ---------- theme ----------
  const cssCache = {};
  EPP.css = (name) => getComputedStyle(root).getPropertyValue(name).trim();
  const NAMED = { electron: '--c-electron', phonon: '--c-phonon', photon: '--c-photon', c4: '--c-4', c5: '--c-5',
    'electron-t': '--t-electron', 'phonon-t': '--t-phonon', 'photon-t': '--t-photon', 'c4-t': '--t-4', 'c5-t': '--t-5',
    ink: '--fg', fg: '--fg', fg2: '--fg-2', muted: '--muted', grid: '--grid', rule: '--rule', accent: '--accent',
    surface: '--surface', surface2: '--surface-2', bg: '--bg' };
  EPP.color = (c) => {
    if (!c) return EPP.css('--fg');
    if (NAMED[c]) return EPP.css(NAMED[c]);
    if (c.startsWith('--')) return EPP.css(c);
    return c;
  };
  // text color for an annotation drawn next to a mark of color c: the hue's text variant, else secondary ink
  EPP.textColor = (c) => (c && NAMED[c + '-t']) ? EPP.css(NAMED[c + '-t']) : (c && (c === 'muted' || c === 'ink' || c === 'fg' || c === 'fg2') ? EPP.color(c) : (c && !NAMED[c] ? EPP.color(c) : EPP.css('--fg-2')));
  EPP.alpha = (c, a) => {
    const col = EPP.color(c);
    if (col.startsWith('#')) {
      let h = col.slice(1);
      if (h.length === 3) h = h.split('').map(x => x + x).join('');
      const n = parseInt(h, 16);
      return `rgba(${(n >> 16) & 255},${(n >> 8) & 255},${n & 255},${a})`;
    }
    const m = col.match(/rgba?\(([^)]+)\)/);
    if (m) { const p = m[1].split(',').map(s => s.trim()); return `rgba(${p[0]},${p[1]},${p[2]},${a})`; }
    return col;
  };
  // parse any CSS color to [r,g,b] via a 1px canvas (cached per string)
  const rgbCache = new Map();
  const pc = document.createElement('canvas'); pc.width = pc.height = 1; const pctx = pc.getContext('2d', { willReadFrequently: true });
  EPP.rgb = (c) => { const col = EPP.color(c); if (rgbCache.has(col)) return rgbCache.get(col); pctx.clearRect(0, 0, 1, 1); pctx.fillStyle = '#000'; pctx.fillStyle = col; pctx.fillRect(0, 0, 1, 1); const d = pctx.getImageData(0, 0, 1, 1).data; const out = [d[0], d[1], d[2]]; rgbCache.set(col, out); return out; };
  const mix = (a, b, t) => [a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t, a[2] + (b[2] - a[2]) * t].map(Math.round);
  // ramp(color, color2?, gamma) -> t in [0,1] -> [r,g,b]. Sequential: surface -> color -> deeper toward --fg.
  // Diverging (color2 given): color2 at 0, surface-gray at 0.5, color at 1.
  EPP.ramp = (color, color2, gamma = 1) => {
    const s = EPP.rgb('--surface'), c = EPP.rgb(color), f = EPP.rgb('--fg');
    if (color2) { const c2 = EPP.rgb(color2), mid = mix(s, EPP.rgb('--muted'), 0.18); return (t) => { t = EPP.clamp(t, 0, 1); return t < 0.5 ? mix(c2, mid, Math.pow(t / 0.5, 1 / gamma)) : mix(mid, c, Math.pow((t - 0.5) / 0.5, gamma)); }; }
    return (t) => { t = Math.pow(EPP.clamp(t, 0, 1), gamma); return t < 0.82 ? mix(s, c, t / 0.82) : mix(c, f, (t - 0.82) / 0.18 * 0.38); };
  };
  const themeCbs = new Set();
  EPP.onTheme = (cb) => { themeCbs.add(cb); return () => themeCbs.delete(cb); };
  const fireTheme = () => { rgbCache.clear(); themeCbs.forEach(cb => { try { cb(); } catch (e) { console.error(e); } }); };
  try { matchMedia('(prefers-color-scheme: dark)').addEventListener('change', fireTheme); } catch (e) {}
  new MutationObserver(fireTheme).observe(root, { attributes: true, attributeFilter: ['data-theme', 'class', 'style'] });
  try { document.fonts.ready.then(fireTheme); } catch (e) {}
  EPP.reducedMotion = (() => { try { return matchMedia('(prefers-reduced-motion: reduce)').matches; } catch (e) { return false; } })();

  // ---------- numerics ----------
  EPP.linspace = (a, b, n) => { const out = new Array(n); for (let i = 0; i < n; i++) out[i] = a + (b - a) * i / (n - 1); return out; };
  EPP.clamp = (x, a, b) => Math.min(b, Math.max(a, x));
  EPP.const = {
    hbar_eVs: 6.582119569e-16, kB_eV: 8.617333262e-5, c_ms: 2.99792458e8, e_C: 1.602176634e-19,
    me_kg: 9.1093837015e-31, amu_kg: 1.66053906660e-27, eps0: 8.8541878128e-12, a0_A: 0.529177210903,
    Ry_eV: 13.605693122994, Ha_eV: 27.211386245988,
    meV_per_cm1: 0.1239841984, // 1 cm^-1 = 0.12398 meV
    THz_per_meV: 0.2417989242 // 1 meV = 0.2418 THz (frequency nu, not omega)
  };
  // Bose-Einstein occupation, E and kT in the same units (or E in eV, T in K via boseT)
  EPP.bose = (E, kT) => (kT <= 0 ? 0 : 1 / Math.expm1(E / kT));
  EPP.boseT = (E_eV, T_K) => EPP.bose(E_eV, EPP.const.kB_eV * T_K);
  EPP.fermi = (E, mu, kT) => (kT <= 0 ? (E < mu ? 1 : (E === mu ? 0.5 : 0)) : 1 / (1 + Math.exp((E - mu) / kT)));
  EPP.fmt = (v, d = 3) => {
    if (!isFinite(v)) return v > 0 ? '∞' : (v < 0 ? '−∞' : '—');
    const a = Math.abs(v);
    let s;
    if (a !== 0 && (a < 1e-3 || a >= 1e5)) s = v.toExponential(Math.max(0, d - 1)).replace('e', '×10^').replace('+', '');
    else s = Number(v.toPrecision(d)).toString();
    return s.replace('-', '−');
  };
  EPP.fix = (v, d = 2) => (isFinite(v) ? v.toFixed(d).replace('-', '−') : '—');

  // minimal complex numbers [re, im]
  const C = EPP.C = {
    c: (re, im = 0) => [re, im],
    add: (a, b) => [a[0] + b[0], a[1] + b[1]],
    sub: (a, b) => [a[0] - b[0], a[1] - b[1]],
    mul: (a, b) => [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]],
    scale: (a, s) => [a[0] * s, a[1] * s],
    div: (a, b) => { const d = b[0] * b[0] + b[1] * b[1]; return [(a[0] * b[0] + a[1] * b[1]) / d, (a[1] * b[0] - a[0] * b[1]) / d]; },
    abs: (a) => Math.hypot(a[0], a[1]),
    arg: (a) => Math.atan2(a[1], a[0]),
    conj: (a) => [a[0], -a[1]],
    inv: (a) => C.div([1, 0], a),
    // principal sqrt with branch cut on negative real axis
    sqrt: (a) => { const r = Math.hypot(a[0], a[1]); const re = Math.sqrt(Math.max(0, (r + a[0]) / 2)); let im = Math.sqrt(Math.max(0, (r - a[0]) / 2)); if (a[1] < 0) im = -im; return [re, im]; },
    exp: (a) => { const e = Math.exp(a[0]); return [e * Math.cos(a[1]), e * Math.sin(a[1])]; },
    log: (a) => [Math.log(Math.hypot(a[0], a[1])), Math.atan2(a[1], a[0])]
  };

  // eigen-decomposition of a real symmetric 2x2 [[a,b],[b,d]] -> {vals:[l1<=l2], vecs:[[x,y],[x,y]]}
  EPP.eigSym2 = (a, b, d) => {
    const tr = a + d, det = a * d - b * b, disc = Math.sqrt(Math.max(0, tr * tr / 4 - det));
    const l1 = tr / 2 - disc, l2 = tr / 2 + disc;
    const vec = (l) => { let x, y; if (Math.abs(b) > 1e-14) { x = b; y = l - a; } else { x = (Math.abs(l - a) < Math.abs(l - d)) ? 1 : 0; y = 1 - x; } const n = Math.hypot(x, y) || 1; return [x / n, y / n]; };
    return { vals: [l1, l2], vecs: [vec(l1), vec(l2)] };
  };
  // Hermitian 2x2 [[a, b],[b*, d]] with a,d real, b complex [re,im] -> eigenvalues and complex eigenvectors
  EPP.eigHerm2 = (a, d, b) => {
    const tr = a + d, bb = b[0] * b[0] + b[1] * b[1];
    const disc = Math.sqrt(Math.max(0, (a - d) * (a - d) / 4 + bb));
    const vals = [tr / 2 - disc, tr / 2 + disc];
    const vecs = vals.map(l => {
      // (a-l) x + b y = 0  -> x = b, y = l - a (complex), fallback for b=0
      let x, y;
      if (bb > 1e-28) { x = [b[0], b[1]]; y = [l - a, 0]; } else { const first = Math.abs(l - a) < Math.abs(l - d); x = [first ? 1 : 0, 0]; y = [first ? 0 : 1, 0]; }
      const n = Math.sqrt(x[0] * x[0] + x[1] * x[1] + y[0] * y[0] + y[1] * y[1]) || 1;
      return [[x[0] / n, x[1] / n], [y[0] / n, y[1] / n]];
    });
    return { vals, vecs };
  };
  // Jacobi eigenvalue algorithm for a real symmetric n x n matrix (array of arrays). Returns ascending.
  EPP.jacobi = (A0, maxSweeps = 60) => {
    const n = A0.length, A = A0.map(r => r.slice()), V = A0.map((r, i) => r.map((_, j) => (i === j ? 1 : 0)));
    for (let s = 0; s < maxSweeps; s++) {
      let off = 0;
      for (let p = 0; p < n; p++) for (let q = p + 1; q < n; q++) off += A[p][q] * A[p][q];
      if (off < 1e-22) break;
      for (let p = 0; p < n; p++) for (let q = p + 1; q < n; q++) {
        if (Math.abs(A[p][q]) < 1e-300) continue;
        const th = (A[q][q] - A[p][p]) / (2 * A[p][q]);
        const t = Math.sign(th || 1) / (Math.abs(th) + Math.sqrt(th * th + 1));
        const c = 1 / Math.sqrt(t * t + 1), sn = t * c;
        for (let k = 0; k < n; k++) { const akp = A[k][p], akq = A[k][q]; A[k][p] = c * akp - sn * akq; A[k][q] = sn * akp + c * akq; }
        for (let k = 0; k < n; k++) { const apk = A[p][k], aqk = A[q][k]; A[p][k] = c * apk - sn * aqk; A[q][k] = sn * apk + c * aqk; }
        for (let k = 0; k < n; k++) { const vkp = V[k][p], vkq = V[k][q]; V[k][p] = c * vkp - sn * vkq; V[k][q] = sn * vkp + c * vkq; }
      }
    }
    const idx = [...Array(n).keys()].sort((i, j) => A[i][i] - A[j][j]);
    return { vals: idx.map(i => A[i][i]), vecs: idx.map(i => V.map(r => r[i])) };
  };
  // ---------- v2 numerics ----------
  // Hermitian n x n eigenproblem by cyclic complex Jacobi rotations. re, im: arrays of rows (im antisymmetric).
  // Returns ascending eigenvalues and eigenvectors vecs[k] = { re: [...], im: [...] }.
  EPP.eigHerm = (re0, im0, maxSweeps = 60) => {
    const n = re0.length;
    const Ar = re0.map(r => r.slice()), Ai = im0.map(r => r.slice());
    const Vr = Ar.map((_, i) => Ar.map((__, j) => (i === j ? 1 : 0))), Vi = Ar.map(() => new Array(n).fill(0));
    for (let sw = 0; sw < maxSweeps; sw++) {
      let off = 0, dg = 0;
      for (let p = 0; p < n; p++) { dg += Ar[p][p] * Ar[p][p]; for (let q = p + 1; q < n; q++) off += Ar[p][q] * Ar[p][q] + Ai[p][q] * Ai[p][q]; }
      if (off <= 1e-28 * (dg + 1e-300)) break;
      for (let p = 0; p < n - 1; p++) for (let q = p + 1; q < n; q++) {
        const xr = Ar[p][q], xi = Ai[p][q], r = Math.hypot(xr, xi);
        if (r < 1e-300) continue;
        // phase step: D = diag(1 at p, e^{-i phi} at q); A <- D^dagger A D makes A_pq = r real
        const cr = xr / r, ci = -xi / r;
        for (let k = 0; k < n; k++) { const a = Ar[k][q], b = Ai[k][q]; Ar[k][q] = a * cr - b * ci; Ai[k][q] = a * ci + b * cr; }
        for (let k = 0; k < n; k++) { const a = Ar[q][k], b = Ai[q][k]; Ar[q][k] = a * cr + b * ci; Ai[q][k] = -a * ci + b * cr; }
        for (let k = 0; k < n; k++) { const a = Vr[k][q], b = Vi[k][q]; Vr[k][q] = a * cr - b * ci; Vi[k][q] = a * ci + b * cr; }
        Ar[p][q] = Ar[q][p] = r; Ai[p][q] = Ai[q][p] = 0; Ai[q][q] = 0;
        // real Jacobi rotation in the (p, q) plane that zeroes A_pq
        const th = (Ar[q][q] - Ar[p][p]) / (2 * r);
        const t = Math.sign(th || 1) / (Math.abs(th) + Math.sqrt(th * th + 1));
        const c = 1 / Math.sqrt(t * t + 1), sn = t * c;
        for (let k = 0; k < n; k++) {
          const ar = Ar[k][p], ai = Ai[k][p], br = Ar[k][q], bi = Ai[k][q];
          Ar[k][p] = c * ar - sn * br; Ai[k][p] = c * ai - sn * bi; Ar[k][q] = sn * ar + c * br; Ai[k][q] = sn * ai + c * bi;
        }
        for (let k = 0; k < n; k++) {
          const ar = Ar[p][k], ai = Ai[p][k], br = Ar[q][k], bi = Ai[q][k];
          Ar[p][k] = c * ar - sn * br; Ai[p][k] = c * ai - sn * bi; Ar[q][k] = sn * ar + c * br; Ai[q][k] = sn * ai + c * bi;
        }
        for (let k = 0; k < n; k++) {
          const ar = Vr[k][p], ai = Vi[k][p], br = Vr[k][q], bi = Vi[k][q];
          Vr[k][p] = c * ar - sn * br; Vi[k][p] = c * ai - sn * bi; Vr[k][q] = sn * ar + c * br; Vi[k][q] = sn * ai + c * bi;
        }
        Ar[p][q] = Ar[q][p] = 0; Ai[p][q] = Ai[q][p] = 0; Ai[p][p] = 0; Ai[q][q] = 0;
      }
    }
    const idx = [...Array(n).keys()].sort((i, j) => Ar[i][i] - Ar[j][j]);
    return { vals: idx.map(i => Ar[i][i]), vecs: idx.map(i => ({ re: Vr.map(r => r[i]), im: Vi.map(r => r[i]) })) };
  };
  // global phase gauge for display: the largest component becomes real and positive
  EPP.fixPhase = (v) => {
    let k = 0, m = -1;
    for (let i = 0; i < v.re.length; i++) { const a = v.re[i] * v.re[i] + v.im[i] * v.im[i]; if (a > m + 1e-12) { m = a; k = i; } }
    const a = Math.sqrt(m) || 1, cr = v.re[k] / a, ci = -v.im[k] / a;
    return { re: v.re.map((x, i) => x * cr - v.im[i] * ci), im: v.re.map((x, i) => x * ci + v.im[i] * cr) };
  };
  // simple bisection root finder on [a,b]
  EPP.bisect = (f, a, b, tol = 1e-10, maxIt = 200) => {
    let fa = f(a), fb = f(b);
    if (!isFinite(fa) || !isFinite(fb) || fa * fb > 0) return NaN;
    for (let i = 0; i < maxIt; i++) { const m = 0.5 * (a + b), fm = f(m); if (fa * fm <= 0) { b = m; fb = fm; } else { a = m; fa = fm; } if (Math.abs(b - a) < tol) break; }
    return 0.5 * (a + b);
  };
  // trapezoid integration of samples y on grid x
  EPP.trapz = (x, y) => { let s = 0; for (let i = 1; i < x.length; i++) s += 0.5 * (y[i] + y[i - 1]) * (x[i] - x[i - 1]); return s; };

  // ---------- rich text: a_b, a_{bc}, a^b, a^{bc} ----------
  function parseRich(s) {
    const out = []; let i = 0, buf = '';
    const flush = () => { if (buf) { out.push({ t: buf, m: 0 }); buf = ''; } };
    while (i < s.length) {
      const ch = s[i];
      if ((ch === '_' || ch === '^') && i + 1 < s.length) {
        flush();
        let j = i + 1, txt;
        if (s[j] === '{') { const k = s.indexOf('}', j); txt = s.slice(j + 1, k < 0 ? s.length : k); i = (k < 0 ? s.length : k + 1); }
        else { txt = s[j]; i = j + 1; }
        out.push({ t: txt, m: ch === '_' ? -1 : 1 });
      } else { buf += ch; i++; }
    }
    flush();
    return out;
  }
  EPP.parseRich = parseRich;
  // draw rich text on canvas; returns width
  EPP.drawRich = (ctx, s, x, y, opt = {}) => {
    const size = opt.size || 12, fam = opt.font || EPP.css('--font-ui'), style = opt.italic ? 'italic ' : '';
    // opt.bold: bold base text (vector symbols); sub- and superscripts stay regular
    const fontOf = (p) => `${style}${opt.bold && !p.m ? '700 ' : ''}${p.m ? Math.round(size * 0.72) : size}px ${fam}`;
    const parts = parseRich(String(s));
    const widths = parts.map(p => { ctx.font = fontOf(p); return ctx.measureText(p.t).width; });
    const total = widths.reduce((a, b) => a + b, 0);
    let x0 = x;
    const align = opt.align || 'left';
    if (align === 'center') x0 = x - total / 2; else if (align === 'right') x0 = x - total;
    const save = ctx.textAlign; ctx.textAlign = 'left';
    parts.forEach((p, k) => {
      ctx.font = fontOf(p);
      const dy = p.m === -1 ? size * 0.28 : (p.m === 1 ? -size * 0.38 : 0);
      ctx.fillText(p.t, x0, y + dy);
      x0 += widths[k];
    });
    ctx.textAlign = save;
    return total;
  };
  EPP.measureRich = (ctx, s, size = 12, fam, opt = {}) => {
    fam = fam || EPP.css('--font-ui');
    const style = opt.italic ? 'italic ' : '';
    return parseRich(String(s)).reduce((w, p) => { ctx.font = `${style}${opt.bold && !p.m ? '700 ' : ''}${p.m ? Math.round(size * 0.72) : size}px ${fam}`; return w + ctx.measureText(p.t).width; }, 0);
  };
  // rich text -> HTML (for legends/readouts)
  EPP.richHTML = (s) => parseRich(String(s)).map(p => {
    const t = p.t.replace(/&/g, '&amp;').replace(/</g, '&lt;');
    return p.m === -1 ? `<sub>${t}</sub>` : p.m === 1 ? `<sup>${t}</sup>` : t;
  }).join('');

  // ---------- nice ticks ----------
  function niceTicks(a, b, target = 6) {
    if (!(isFinite(a) && isFinite(b)) || a === b) return [a];
    const lo = Math.min(a, b), hi = Math.max(a, b);
    const raw = (hi - lo) / target, mag = Math.pow(10, Math.floor(Math.log10(raw))), r = raw / mag;
    const step = (r < 1.5 ? 1 : r < 3 ? 2 : r < 7 ? 5 : 10) * mag;
    const out = []; for (let v = Math.ceil(lo / step - 1e-9) * step; v <= hi + step * 1e-9; v += step) out.push(Math.abs(v) < step * 1e-9 ? 0 : v);
    return out;
  }
  EPP.niceTicks = niceTicks;
  const tickLabel = (v) => { const a = Math.abs(v); if (a === 0) return '0'; if (a >= 1e4 || a < 1e-3) return v.toExponential(0).replace('-', '−'); return String(Number(v.toPrecision(4))).replace('-', '−'); };

  // ---------- responsive canvas ----------
  EPP.canvas = (container, opt = {}) => {
    const el = typeof container === 'string' ? document.querySelector(container) : container;
    el.classList.add('plot');
    const cv = document.createElement('canvas');
    el.appendChild(cv);
    const ctx = cv.getContext('2d');
    const self = { el, canvas: cv, ctx, w: 0, h: 0, dpr: 1, _cbs: [] };
    const aspect = opt.aspect || 0.6, minH = opt.minHeight || 180, maxH = opt.maxHeight || 520;
    const resize = () => {
      const w = Math.max(1, el.clientWidth);
      const h = opt.height ? opt.height : Math.round(EPP.clamp(w * aspect, minH, maxH));
      const dpr = Math.min(3, window.devicePixelRatio || 1);
      if (w === self.w && h === self.h && dpr === self.dpr) return false;
      self.w = w; self.h = h; self.dpr = dpr;
      cv.width = Math.round(w * dpr); cv.height = Math.round(h * dpr);
      cv.style.height = h + 'px';
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      return true;
    };
    self.onResize = (cb) => { self._cbs.push(cb); };
    self.clear = () => { ctx.setTransform(self.dpr, 0, 0, self.dpr, 0, 0); ctx.clearRect(0, 0, self.w, self.h); };
    resize();
    let raf = 0;
    const ro = new ResizeObserver(() => { cancelAnimationFrame(raf); raf = requestAnimationFrame(() => { if (resize()) self._cbs.forEach(cb => cb()); }); });
    ro.observe(el);
    EPP.onTheme(() => self._cbs.forEach(cb => cb()));
    return self;
  };

  // ---------- Plot ----------
  // new EPP.Plot(el, {xlabel, ylabel, xlim, ylim, aspect, height, legend, hover, xfmt, yfmt, margin})
  // plot.draw(layers) ; layers: {type:'line'|'fill'|'points'|'hline'|'vline'|'band'|'hband'|'text'|'arrow'|'segments', ...}
  class Plot {
    constructor(container, opt = {}) {
      this.opt = Object.assign({ aspect: 0.58, legend: true, hover: true, grid: true }, opt);
      const el = typeof container === 'string' ? document.querySelector(container) : container;
      this.wrap = el;
      this.holder = document.createElement('div');
      el.appendChild(this.holder);
      this.cv = EPP.canvas(this.holder, { aspect: this.opt.aspect, height: this.opt.height, minHeight: this.opt.minHeight || 200, maxHeight: this.opt.maxHeight || 460 });
      this.ctx = this.cv.ctx;
      this.layers = [];
      this.xlim = opt.xlim || [0, 1]; this.ylim = opt.ylim || [0, 1];
      if (this.opt.legend) { this.legendEl = document.createElement('div'); this.legendEl.className = 'legend'; el.appendChild(this.legendEl); }
      this.tip = document.createElement('div'); this.tip.className = 'tip'; this.tip.hidden = true; this.holder.appendChild(this.tip);
      this.cv.onResize(() => this.render());
      this._hoverX = null;
      const cvs = this.cv.canvas;
      const toData = (e) => { const r = cvs.getBoundingClientRect(); const px = e.clientX - r.left, py = e.clientY - r.top; return { px, py, x: this.invX(px), y: this.invY(py) }; };
      this._dragging = false;
      cvs.addEventListener('pointerdown', (e) => { const d = toData(e); if (this.onPointer) { this._dragging = true; try { cvs.setPointerCapture(e.pointerId); } catch (_) {} this.onPointer(d.x, d.y, 'down', d); } });
      cvs.addEventListener('pointermove', (e) => {
        const d = toData(e);
        if (this._dragging && this.onPointer) this.onPointer(d.x, d.y, 'move', d);
        if (this.opt.hover && e.pointerType !== 'touch') { const inside = d.px >= this.m.l && d.px <= this.cv.w - this.m.r && d.py >= this.m.t && d.py <= this.cv.h - this.m.b; this._hoverX = inside || (d.px >= this.m.l && d.px <= this.cv.w - this.m.r) ? d.x : null; this._hoverY = inside ? d.y : null; this.render(); }
      });
      const end = (e) => { if (this._dragging && this.onPointer) { const d = toData(e); this.onPointer(d.x, d.y, 'up', d); } this._dragging = false; };
      cvs.addEventListener('pointerup', end); cvs.addEventListener('pointercancel', end);
      cvs.addEventListener('pointerleave', () => { this._hoverX = null; this.tip.hidden = true; this.render(); });
      if (this.opt.onPointer) this.onPointer = this.opt.onPointer;
      if (this.opt.ariaLabel) { cvs.setAttribute('role', 'img'); cvs.setAttribute('aria-label', this.opt.ariaLabel); }
    }
    get m() {
      const mm = this.opt.margin || {};
      return { l: mm.l ?? (this.opt.ylabel ? 58 : 44), r: mm.r ?? 14, t: mm.t ?? 12, b: mm.b ?? (this.opt.xlabel ? 44 : 28) };
    }
    X(x) { const m = this.m, [a, b] = this.xlim; if (this.opt.logx) return m.l + (Math.log10(x) - Math.log10(a)) / (Math.log10(b) - Math.log10(a)) * (this.cv.w - m.l - m.r); return m.l + (x - a) / (b - a) * (this.cv.w - m.l - m.r); }
    Y(y) { const m = this.m, [a, b] = this.ylim; if (this.opt.logy) return this.cv.h - m.b - (Math.log10(y) - Math.log10(a)) / (Math.log10(b) - Math.log10(a)) * (this.cv.h - m.t - m.b); return this.cv.h - m.b - (y - a) / (b - a) * (this.cv.h - m.t - m.b); }
    invX(px) { const m = this.m, [a, b] = this.xlim, f = (px - m.l) / (this.cv.w - m.l - m.r); return this.opt.logx ? Math.pow(10, Math.log10(a) + f * (Math.log10(b) - Math.log10(a))) : a + f * (b - a); }
    invY(py) { const m = this.m, [a, b] = this.ylim, f = (this.cv.h - m.b - py) / (this.cv.h - m.t - m.b); return this.opt.logy ? Math.pow(10, Math.log10(a) + f * (Math.log10(b) - Math.log10(a))) : a + f * (b - a); }
    setLimits({ xlim, ylim, xlabel, ylabel } = {}) { if (xlim) this.xlim = xlim; if (ylim) this.ylim = ylim; if (xlabel !== undefined) this.opt.xlabel = xlabel; if (ylabel !== undefined) this.opt.ylabel = ylabel; }
    draw(layers, lim) { if (lim) this.setLimits(lim); this.layers = layers.filter(Boolean); this.render(); this.renderLegend(); }
    renderLegend() {
      if (!this.legendEl) return;
      const items = this.layers.filter(L => L.label && !L.nolegend && L.type !== 'heatmap');
      this.legendEl.innerHTML = items.map(L => {
        const col = EPP.color(L.color);
        const cls = L.type === 'fill' || L.type === 'band' || L.type === 'hband' ? 'box' : (L.dash ? (L.dash[0] <= 2 ? 'dot' : 'dash') : '');
        return `<span class="item" style="color:${col}"><span class="sw ${cls}"></span><span style="color:var(--fg-2)">${EPP.richHTML(L.label)}</span></span>`;
      }).join('');
      this.legendEl.hidden = items.length === 0;
    }
    render() {
      const ctx = this.ctx, w = this.cv.w, h = this.cv.h, m = this.m, o = this.opt;
      if (!w || !h) return;
      this.cv.clear();
      const fUI = EPP.css('--font-ui'), cGrid = EPP.css('--grid'), cMuted = EPP.css('--muted'), cFg = EPP.css('--fg'), cRule = EPP.css('--rule');
      const [xa, xb] = this.xlim, [ya, yb] = this.ylim;
      // ticks + grid
      const xt = o.xticks || (o.logx ? this._logTicks(xa, xb) : niceTicks(xa, xb, Math.max(3, Math.round((w - m.l - m.r) / 90))));
      const yt = o.yticks || (o.logy ? this._logTicks(ya, yb) : niceTicks(ya, yb, Math.max(3, Math.round((h - m.t - m.b) / 55))));
      ctx.lineWidth = 1; ctx.font = `11px ${fUI}`; ctx.fillStyle = cMuted;
      ctx.textAlign = 'center'; ctx.textBaseline = 'top';
      xt.forEach(t => { const v = typeof t === 'object' ? t.v : t, lab = typeof t === 'object' ? t.label : (o.xfmt ? o.xfmt(v) : tickLabel(v)); const X = this.X(v); if (X < m.l - 0.5 || X > w - m.r + 0.5) return; if (o.grid) { ctx.strokeStyle = cGrid; ctx.beginPath(); ctx.moveTo(Math.round(X) + 0.5, m.t); ctx.lineTo(Math.round(X) + 0.5, h - m.b); ctx.stroke(); } EPP.drawRich(ctx, lab, X, h - m.b + 6, { size: 11, align: 'center', font: fUI }); });
      ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
      yt.forEach(t => { const v = typeof t === 'object' ? t.v : t, lab = typeof t === 'object' ? t.label : (o.yfmt ? o.yfmt(v) : tickLabel(v)); const Y = this.Y(v); if (Y < m.t - 0.5 || Y > h - m.b + 0.5) return; if (o.grid) { ctx.strokeStyle = cGrid; ctx.beginPath(); ctx.moveTo(m.l, Math.round(Y) + 0.5); ctx.lineTo(w - m.r, Math.round(Y) + 0.5); ctx.stroke(); } EPP.drawRich(ctx, lab, m.l - 6, Y, { size: 11, align: 'right', font: fUI }); });
      // zero lines
      ctx.strokeStyle = cRule; ctx.lineWidth = 1;
      if (ya < 0 && yb > 0 && !o.logy) { const Y = Math.round(this.Y(0)) + 0.5; ctx.beginPath(); ctx.moveTo(m.l, Y); ctx.lineTo(w - m.r, Y); ctx.stroke(); }
      if (xa < 0 && xb > 0 && !o.logx) { const X = Math.round(this.X(0)) + 0.5; ctx.beginPath(); ctx.moveTo(X, m.t); ctx.lineTo(X, h - m.b); ctx.stroke(); }
      // frame
      ctx.strokeStyle = cRule; ctx.strokeRect(m.l + 0.5, m.t + 0.5, w - m.l - m.r - 1, h - m.t - m.b - 1);
      // axis labels
      ctx.fillStyle = cFg; ctx.textBaseline = 'alphabetic';
      if (o.xlabel) EPP.drawRich(ctx, o.xlabel, m.l + (w - m.l - m.r) / 2, h - 8, { size: 12.5, align: 'center', font: fUI });
      if (o.ylabel) { ctx.save(); ctx.translate(14, m.t + (h - m.t - m.b) / 2); ctx.rotate(-Math.PI / 2); EPP.drawRich(ctx, o.ylabel, 0, 0, { size: 12.5, align: 'center', font: fUI }); ctx.restore(); }
      // layers (clipped)
      ctx.save(); ctx.beginPath(); ctx.rect(m.l, m.t, w - m.l - m.r, h - m.t - m.b); ctx.clip();
      for (const L of this.layers) this._layer(L, false);
      ctx.restore();
      for (const L of this.layers) if (L.type === 'text' || L.type === 'label') this._layer(L, true);
      // colorbars for heatmaps (drawn in the right margin; give margin.r >= 54)
      for (const L of this.layers) if (L.type === 'heatmap' && L.colorbar !== false && m.r >= 40) {
        const ramp = EPP.ramp(L.color || 'electron', L.color2, L.gamma || 1), x = w - m.r + 10, bw = 9, top = m.t, bot = h - m.b, n = 64;
        for (let k = 0; k < n; k++) { const c = ramp(k / (n - 1)); ctx.fillStyle = `rgb(${c[0]},${c[1]},${c[2]})`; const y1 = bot - (k + 1) * (bot - top) / n; ctx.fillRect(x, y1, bw, (bot - top) / n + 0.6); }
        ctx.strokeStyle = cRule; ctx.strokeRect(x + 0.5, top + 0.5, bw, bot - top - 1);
        ctx.fillStyle = cMuted; ctx.textAlign = 'left'; ctx.textBaseline = 'middle';
        const f = L.cfmt || ((v) => EPP.fmt(v, 2));
        EPP.drawRich(ctx, f(L.vmax), x + bw + 3, top + 5, { size: 10 }); EPP.drawRich(ctx, f(L.vmin), x + bw + 3, bot - 5, { size: 10 });
        if (L.colorbar && typeof L.colorbar === 'string') { ctx.save(); ctx.translate(w - 6, (top + bot) / 2); ctx.rotate(-Math.PI / 2); ctx.textAlign = 'center'; EPP.drawRich(ctx, L.colorbar, 0, 0, { size: 10.5 }); ctx.restore(); }
      }
      // hover
      this._drawHover();
    }
    _logTicks(a, b) { const out = []; for (let e = Math.floor(Math.log10(a)); e <= Math.ceil(Math.log10(b)); e++) { const v = Math.pow(10, e); if (v >= a * 0.999 && v <= b * 1.001) out.push({ v, label: e >= -2 && e <= 3 ? String(v) : `10^{${e}}` }); } return out; }
    _layer(L, unclipped) {
      const ctx = this.ctx, m = this.m, w = this.cv.w, h = this.cv.h;
      const col = EPP.color(L.color || 'ink');
      ctx.setLineDash(L.dash || []);
      ctx.lineWidth = L.width || 2; ctx.strokeStyle = col; ctx.fillStyle = col; ctx.globalAlpha = L.opacity ?? 1;
      const t = L.type;
      if (t === 'line') {
        const x = L.x, y = L.y, ylo = this.ylim[0], yhi = this.ylim[1], span = Math.abs(yhi - ylo);
        const jump = L.breakJump ?? 3 * span;
        ctx.beginPath(); let pen = false, prev = null;
        for (let i = 0; i < x.length; i++) {
          const yi = y[i];
          if (!isFinite(yi) || !isFinite(x[i]) || (prev !== null && Math.abs(yi - prev) > jump)) { pen = false; prev = isFinite(yi) ? yi : null; if (!isFinite(yi)) continue; }
          const X = this.X(x[i]); let Y = this.Y(yi); Y = EPP.clamp(Y, -1e4, 1e4);
          if (!pen) { ctx.moveTo(X, Y); pen = true; } else ctx.lineTo(X, Y);
          prev = yi;
        }
        ctx.lineJoin = 'round'; ctx.lineCap = 'round'; ctx.stroke();
      } else if (t === 'heatmap') {
        // L.z: row-major, j over y (y0 -> y1), i over x (x0 -> x1), length nx*ny
        const sig = [EPP.css('--surface'), col, EPP.css('--fg'), L.vmin, L.vmax, L.color2 || ''].join('|');
        if (!L._img || L._sig !== sig || L._z !== L.z) {
          const oc = document.createElement('canvas'); oc.width = L.nx; oc.height = L.ny;
          const g = oc.getContext('2d'), im = g.createImageData(L.nx, L.ny);
          const ramp = EPP.ramp(L.color || 'electron', L.color2, L.gamma || 1);
          const lo = L.vmin, hi = L.vmax, div = !!L.color2;
          for (let j = 0; j < L.ny; j++) for (let i = 0; i < L.nx; i++) {
            const v = L.z[j * L.nx + i]; let tt;
            if (!isFinite(v)) tt = div ? 0.5 : 0; else tt = EPP.clamp((v - lo) / (hi - lo), 0, 1);
            const c = ramp(tt), k = ((L.ny - 1 - j) * L.nx + i) * 4;
            im.data[k] = c[0]; im.data[k + 1] = c[1]; im.data[k + 2] = c[2]; im.data[k + 3] = 255;
          }
          g.putImageData(im, 0, 0); L._img = oc; L._sig = sig; L._z = L.z;
        }
        ctx.imageSmoothingEnabled = L.smooth !== false;
        const X0 = this.X(L.x0), X1 = this.X(L.x1), Y0 = this.Y(L.y0), Y1 = this.Y(L.y1);
        ctx.drawImage(L._img, Math.min(X0, X1), Math.min(Y0, Y1), Math.abs(X1 - X0), Math.abs(Y1 - Y0));
      } else if (t === 'segments') { // L.segs: [[x0,y0,x1,y1],...]
        ctx.beginPath(); for (const s of L.segs) { ctx.moveTo(this.X(s[0]), this.Y(s[1])); ctx.lineTo(this.X(s[2]), this.Y(s[3])); } ctx.stroke();
      } else if (t === 'fill') {
        ctx.globalAlpha = L.alpha ?? 0.18; ctx.beginPath();
        const x = L.x, y1 = L.y1, y2 = L.y2 || L.x.map(() => (L.base ?? 0));
        let started = false;
        for (let i = 0; i < x.length; i++) { if (!isFinite(y1[i])) continue; const X = this.X(x[i]), Y = EPP.clamp(this.Y(y1[i]), -1e4, 1e4); if (!started) { ctx.moveTo(X, Y); started = true; } else ctx.lineTo(X, Y); }
        for (let i = x.length - 1; i >= 0; i--) { if (!isFinite(y2[i])) continue; ctx.lineTo(this.X(x[i]), EPP.clamp(this.Y(y2[i]), -1e4, 1e4)); }
        ctx.closePath(); ctx.fill();
      } else if (t === 'points') {
        const r = L.r || 3.5;
        for (let i = 0; i < L.x.length; i++) { if (!isFinite(L.y[i])) continue; ctx.beginPath(); ctx.arc(this.X(L.x[i]), this.Y(L.y[i]), r, 0, 2 * Math.PI); if (L.hollow) { ctx.lineWidth = 1.6; ctx.stroke(); } else ctx.fill(); }
      } else if (t === 'hline' || t === 'vline') {
        ctx.lineWidth = L.width || 1.2; if (!L.dash) ctx.setLineDash([5, 4]);
        ctx.beginPath();
        if (t === 'hline') { const Y = this.Y(L.at); ctx.moveTo(m.l, Y); ctx.lineTo(w - m.r, Y); } else { const X = this.X(L.at); ctx.moveTo(X, m.t); ctx.lineTo(X, h - m.b); }
        ctx.stroke(); ctx.setLineDash([]);
        if (L.text) { ctx.fillStyle = EPP.textColor(L.color); ctx.font = `11px ${EPP.css('--font-ui')}`; if (t === 'hline') { ctx.textBaseline = 'bottom'; EPP.drawRich(ctx, L.text, (L.textAt === 'left' ? m.l + 6 : w - m.r - 6), this.Y(L.at) - 3, { size: 11, align: L.textAt === 'left' ? 'left' : 'right' }); } else { ctx.textBaseline = 'top'; EPP.drawRich(ctx, L.text, this.X(L.at) + 4, m.t + 4 + (L.textDy || 0), { size: 11 }); } }
      } else if (t === 'band' || t === 'hband') {
        ctx.globalAlpha = L.alpha ?? 0.12;
        if (t === 'band') { const X0 = this.X(L.x0), X1 = this.X(L.x1); ctx.fillRect(Math.min(X0, X1), m.t, Math.abs(X1 - X0), h - m.t - m.b); }
        else { const Y0 = this.Y(L.y0), Y1 = this.Y(L.y1); ctx.fillRect(m.l, Math.min(Y0, Y1), w - m.l - m.r, Math.abs(Y1 - Y0)); }
        ctx.globalAlpha = 1;
        if (L.text) { ctx.fillStyle = EPP.textColor(L.color); ctx.textBaseline = 'top'; if (t === 'band') EPP.drawRich(ctx, L.text, (this.X(L.x0) + this.X(L.x1)) / 2, m.t + 4, { size: 11, align: 'center' }); else EPP.drawRich(ctx, L.text, w - m.r - 6, Math.min(this.Y(L.y0), this.Y(L.y1)) + 3, { size: 11, align: 'right' }); }
      } else if ((t === 'text' || t === 'label') && unclipped) {
        ctx.globalAlpha = L.opacity ?? 1; ctx.textBaseline = L.baseline || 'middle'; ctx.fillStyle = EPP.textColor(L.color || 'fg2');
        const X = this.X(L.x) + (L.dx || 0), Y = this.Y(L.y) + (L.dy || 0);
        EPP.drawRich(ctx, L.text, X, Y, { size: L.size || 12, align: L.align || 'left', italic: L.italic, bold: L.bold });
      } else if (t === 'arrow') {
        const X0 = this.X(L.x0), Y0 = this.Y(L.y0), X1 = this.X(L.x1), Y1 = this.Y(L.y1);
        ctx.lineWidth = L.width || 1.6; ctx.beginPath(); ctx.moveTo(X0, Y0); ctx.lineTo(X1, Y1); ctx.stroke();
        const a = Math.atan2(Y1 - Y0, X1 - X0), s = L.head || 8; ctx.setLineDash([]);
        ctx.beginPath(); ctx.moveTo(X1, Y1); ctx.lineTo(X1 - s * Math.cos(a - 0.4), Y1 - s * Math.sin(a - 0.4)); ctx.lineTo(X1 - s * Math.cos(a + 0.4), Y1 - s * Math.sin(a + 0.4)); ctx.closePath(); ctx.fill();
      }
      ctx.globalAlpha = 1; ctx.setLineDash([]);
    }
    _drawHover() {
      const x = this._hoverX;
      if (x === null || x === undefined) { this.tip.hidden = true; return; }
      const lines = this.layers.filter(L => L.type === 'line' && L.label && L.hover !== false);
      const heats = this.layers.filter(L => L.type === 'heatmap' && L.hover !== false);
      if (!lines.length && !(heats.length && this._hoverY != null)) { this.tip.hidden = true; return; }
      const ctx = this.ctx, m = this.m, X = this.X(x);
      ctx.strokeStyle = EPP.css('--muted'); ctx.lineWidth = 1; ctx.setLineDash([2, 3]);
      ctx.beginPath(); ctx.moveTo(X, m.t); ctx.lineTo(X, this.cv.h - m.b); ctx.stroke(); ctx.setLineDash([]);
      const xf = this.opt.hoverX || ((v) => EPP.fmt(v, 3));
      const yf = this.opt.hoverY || ((v) => EPP.fmt(v, 4));
      const rows = [`<span style="color:var(--muted)">${EPP.richHTML(this.opt.xname || 'x')} = ${xf(x)}</span>`];
      for (const L of lines) {
        // nearest sample by x (assumes monotonic x)
        const xs = L.x; let lo = 0, hi = xs.length - 1; const asc = xs[hi] >= xs[0];
        while (hi - lo > 1) { const mid = (lo + hi) >> 1; if ((xs[mid] < x) === asc) lo = mid; else hi = mid; }
        const i = Math.abs(xs[lo] - x) < Math.abs(xs[hi] - x) ? lo : hi;
        const y = L.y[i];
        if (!isFinite(y) || Math.abs(xs[i] - x) > Math.abs(this.xlim[1] - this.xlim[0]) * 0.03) continue;
        const Y = this.Y(y);
        if (Y >= m.t && Y <= this.cv.h - m.b) { ctx.fillStyle = EPP.color(L.color); ctx.beginPath(); ctx.arc(this.X(xs[i]), Y, 3.5, 0, 2 * Math.PI); ctx.fill(); }
        rows.push(`<span style="color:${EPP.color(L.color)}">■</span> ${EPP.richHTML(L.label)}: ${yf(y)}`);
      }
      if (this._hoverY != null) for (const L of heats) {
        const i = Math.round((x - L.x0) / (L.x1 - L.x0) * (L.nx - 1)), j = Math.round((this._hoverY - L.y0) / (L.y1 - L.y0) * (L.ny - 1));
        if (i >= 0 && i < L.nx && j >= 0 && j < L.ny) { const yf2 = this.opt.hoverX2 || ((v) => EPP.fmt(v, 3)); rows.push(`<span style="color:var(--muted)">${EPP.richHTML(this.opt.yname || 'y')} = ${yf2(this._hoverY)}</span>`); rows.push(`${EPP.richHTML(L.label || 'value')}: ${(L.cfmt || ((v) => EPP.fmt(v, 3)))(L.z[j * L.nx + i])}`); }
      }
      if (rows.length <= 1) { this.tip.hidden = true; return; }
      this.tip.innerHTML = rows.join('<br>');
      this.tip.hidden = false;
      const tw = this.tip.offsetWidth, left = X + 12 + tw > this.cv.w ? X - 12 - tw : X + 12;
      this.tip.style.left = Math.max(0, left) + 'px'; this.tip.style.top = (m.t + 4) + 'px';
    }
  }
  EPP.Plot = Plot;
  EPP.sample = (f, a, b, n = 400) => { const x = EPP.linspace(a, b, n); return { x, y: x.map(f) }; };

  // ---------- controls ----------
  // EPP.controls(el, specs, onChange) ; spec: {key, type:'range'|'select'|'toggle'|'seg'|'button', label, min, max, step, value, fmt, options:[[value,label],...], onClick}
  EPP.controls = (container, specs, onChange) => {
    const el = typeof container === 'string' ? document.querySelector(container) : container;
    el.classList.add('controls');
    const fig = el.closest('figure, section, div[id]');
    const base = (fig && fig.id) ? fig.id : ('c' + Math.random().toString(36).slice(2, 7));
    const state = {}, setters = {};
    let scheduled = false;
    const fire = () => { if (scheduled) return; scheduled = true; requestAnimationFrame(() => { scheduled = false; onChange && onChange(state); }); };
    specs.forEach(s => {
      const id = `${base}-${s.key}`;
      const box = document.createElement('div'); box.className = 'ctrl' + (s.type === 'toggle' ? ' toggle' : '');
      const longSeg = s.type === 'seg' && s.options.reduce((n, o) => n + String(o[1]).length + 3, 0) > 30;
      if (s.wide || longSeg) box.style.gridColumn = '1 / -1';
      if (s.type === 'range') {
        const fmt = s.fmt || ((v) => String(v));
        box.innerHTML = `<label for="${id}"><span>${EPP.richHTML(s.label)}</span><output for="${id}"></output></label><input type="range" id="${id}" min="${s.min}" max="${s.max}" step="${s.step ?? 'any'}" value="${s.value}">`;
        const inp = box.querySelector('input'), out = box.querySelector('output');
        const map = s.log ? { to: (v) => Math.pow(10, v), from: (v) => Math.log10(v) } : null;
        if (map) { inp.min = Math.log10(s.min); inp.max = Math.log10(s.max); inp.step = s.step ?? 'any'; inp.value = Math.log10(s.value); }
        const read = () => { const v = parseFloat(inp.value); state[s.key] = map ? map.to(v) : v; out.textContent = fmt(state[s.key]); };
        inp.addEventListener('input', () => { read(); fire(); });
        setters[s.key] = (v) => { inp.value = map ? Math.log10(v) : v; read(); };
        read();
      } else if (s.type === 'select') {
        box.innerHTML = `<label for="${id}"><span>${EPP.richHTML(s.label)}</span></label><select id="${id}">${s.options.map(o => `<option value="${o[0]}">${o[1]}</option>`).join('')}</select>`;
        const sel = box.querySelector('select'); sel.value = s.value ?? s.options[0][0];
        const read = () => { state[s.key] = sel.value; };
        sel.addEventListener('change', () => { read(); fire(); });
        setters[s.key] = (v) => { sel.value = v; read(); };
        read();
      } else if (s.type === 'toggle') {
        box.innerHTML = `<label for="${id}"><input type="checkbox" id="${id}" ${s.value ? 'checked' : ''}><span>${EPP.richHTML(s.label)}</span></label>`;
        const cb = box.querySelector('input');
        const read = () => { state[s.key] = cb.checked; };
        cb.addEventListener('change', () => { read(); fire(); });
        setters[s.key] = (v) => { cb.checked = !!v; read(); };
        read();
      } else if (s.type === 'seg') {
        box.innerHTML = `<span class="lab"><span>${EPP.richHTML(s.label)}</span></span><div class="seg" role="group" aria-label="${s.label.replace(/[_^{}]/g, '')}">${s.options.map(o => `<button type="button" id="${id}-${o[0]}" data-v="${o[0]}" aria-pressed="false">${EPP.richHTML(o[1])}</button>`).join('')}</div>`;
        const btns = [...box.querySelectorAll('button')];
        const setv = (v) => { state[s.key] = v; btns.forEach(b => b.setAttribute('aria-pressed', String(b.dataset.v === String(v)))); };
        btns.forEach(b => b.addEventListener('click', () => { setv(b.dataset.v); fire(); }));
        setters[s.key] = setv; setv(s.value ?? s.options[0][0]);
      } else if (s.type === 'button') {
        box.innerHTML = `<span class="lab"><span>&nbsp;</span></span><button type="button" class="btn" id="${id}">${EPP.richHTML(s.label)}</button>`;
        box.querySelector('button').addEventListener('click', () => { s.onClick && s.onClick(state); fire(); });
      }
      el.appendChild(box);
    });
    const api = { state, set: (k, v, silent) => { setters[k] && setters[k](v); if (!silent) fire(); }, refresh: () => onChange && onChange(state), el };
    // initial draw after fonts/layout
    requestAnimationFrame(() => onChange && onChange(state));
    return api;
  };
  // readout: EPP.readout(el).set({ 'ω_LO': '45 meV', ... })
  EPP.readout = (container) => {
    const el = typeof container === 'string' ? document.querySelector(container) : container;
    el.classList.add('readout'); el.setAttribute('aria-live', 'polite');
    return { set: (obj) => { el.innerHTML = Object.entries(obj).map(([k, v]) => `<span>${EPP.richHTML(k)}<b>${v}</b></span>`).join(''); } };
  };

  // ---------- animation ----------
  // EPP.anim(figureEl, frame(t, dt), {buttonHost}) -> {play, pause, toggle, running}. Pauses offscreen; starts paused under reduced motion.
  EPP.anim = (host, frame, opt = {}) => {
    const el = typeof host === 'string' ? document.querySelector(host) : host;
    let running = false, visible = true, raf = 0, last = 0, t = 0;
    const btn = document.createElement('button'); btn.type = 'button'; btn.className = 'btn';
    const label = () => { btn.textContent = running ? 'Pause' : 'Play'; btn.setAttribute('aria-pressed', String(running)); };
    const loop = (ts) => { if (!running || !visible) { raf = 0; return; } const dt = last ? Math.min(0.05, (ts - last) / 1000) : 0; last = ts; t += dt; frame(t, dt); raf = requestAnimationFrame(loop); };
    const kick = () => { if (running && visible && !raf) { last = 0; raf = requestAnimationFrame(loop); } };
    const api = {
      play() { running = true; label(); kick(); }, pause() { running = false; label(); },
      toggle() { running ? api.pause() : api.play(); }, get running() { return running; }, get t() { return t; },
      redraw() { frame(t, 0); }, button: btn
    };
    btn.addEventListener('click', api.toggle);
    if (opt.buttonHost) { const bh = typeof opt.buttonHost === 'string' ? el.querySelector(opt.buttonHost) : opt.buttonHost; bh.appendChild(btn); }
    try { new IntersectionObserver((ents) => { visible = ents[0].isIntersecting; kick(); }).observe(el); } catch (e) {}
    label();
    requestAnimationFrame(() => { frame(0, 0); if (!EPP.reducedMotion && opt.autoplay !== false) api.play(); });
    return api;
  };

  // ---------- Feynman diagrams (SVG) ----------
  // EPP.fd.draw(el, {w, h, items:[...], label}) ; item types:
  //   {t:'fermion'|'fermion2'|'phonon'|'phonon2'|'photon'|'photon2'|'coulomb'|'screened', a:[x,y], b:[x,y], bend:0, arrow:true, label, lpos:'+'|'-', ldist}
  //   {t:'vertex', at:[x,y], r}   {t:'blob', at:[x,y], r, label}   {t:'text', at:[x,y], text, color, size, anchor}
  const NS = 'http://www.w3.org/2000/svg';
  function curvePoints(a, b, bend, n) {
    const pts = []; const dx = b[0] - a[0], dy = b[1] - a[1], L = Math.hypot(dx, dy);
    if (!bend) { for (let i = 0; i <= n; i++) { const s = i / n; pts.push([a[0] + dx * s, a[1] + dy * s, dx / L, dy / L]); } return pts; }
    // circular arc through a, b with sagitta = bend * L (positive bulges to the screen-left of a->b: up for a left-to-right line)
    const sag = bend * L, R = (L * L / 4 + sag * sag) / (2 * Math.abs(sag));
    const mx = (a[0] + b[0]) / 2, my = (a[1] + b[1]) / 2, nx = dy / L, ny = -dx / L; // screen-left normal (y down)
    const sgn = Math.sign(sag), d = R - Math.abs(sag);
    const cx = mx - sgn * nx * d, cy = my - sgn * ny * d;
    let t0 = Math.atan2(a[1] - cy, a[0] - cx), t1 = Math.atan2(b[1] - cy, b[0] - cx);
    // choose the arc passing through the sagitta point
    const px = mx + sgn * nx * Math.abs(sag), py = my + sgn * ny * Math.abs(sag), tp = Math.atan2(py - cy, px - cx);
    const norm = (x) => { while (x < 0) x += 2 * Math.PI; while (x >= 2 * Math.PI) x -= 2 * Math.PI; return x; };
    let span = norm(t1 - t0); const tpRel = norm(tp - t0);
    if (tpRel > span) span = span - 2 * Math.PI;
    for (let i = 0; i <= n; i++) { const th = t0 + span * i / n; const x = cx + R * Math.cos(th), y = cy + R * Math.sin(th); const tx = -Math.sin(th) * Math.sign(span), ty = Math.cos(th) * Math.sign(span); pts.push([x, y, tx, ty]); }
    return pts;
  }
  const pathFrom = (pts) => 'M' + pts.map(p => p[0].toFixed(2) + ' ' + p[1].toFixed(2)).join(' L');
  function offsetPath(pts, f) { return pts.map((p, i) => { const s = i / (pts.length - 1); const o = f(s); return [p[0] - p[3] * o, p[1] + p[2] * o, p[2], p[3]]; }); }
  function arcLen(pts) { let L = 0; for (let i = 1; i < pts.length; i++) L += Math.hypot(pts[i][0] - pts[i - 1][0], pts[i][1] - pts[i - 1][1]); return L; }
  function svgEl(tag, attrs, parent) { const e = document.createElementNS(NS, tag); for (const k in attrs) e.setAttribute(k, attrs[k]); if (parent) parent.appendChild(e); return e; }
  function richSVG(textEl, s) {
    parseRich(String(s)).forEach(p => {
      const ts = svgEl('tspan', {}, textEl); ts.textContent = p.t;
      if (p.m) { ts.setAttribute('font-size', '72%'); ts.setAttribute('baseline-shift', p.m === -1 ? 'sub' : 'super'); ts.setAttribute('font-style', 'normal'); }
    });
  }
  EPP.fd = {
    draw(container, spec) {
      const el = typeof container === 'string' ? document.querySelector(container) : container;
      const w = spec.w || 200, h = spec.h || 100;
      const svg = svgEl('svg', { viewBox: `${spec.x0 || 0} ${spec.y0 || 0} ${w} ${h}`, width: w, height: h, role: 'img', 'aria-label': spec.label || 'Feynman diagram' });
      const paint = () => {
        svg.innerHTML = '';
        const col = { fermion: EPP.color('electron'), phonon: EPP.color('phonon'), photon: EPP.color('photon'), ink: EPP.color('ink'), muted: EPP.color('muted') };
        const tcol = { fermion: EPP.color('electron-t'), phonon: EPP.color('phonon-t'), photon: EPP.color('photon-t') };
        for (const it of spec.items) {
          const t = it.t;
          if (['fermion', 'fermion2', 'phonon', 'phonon2', 'photon', 'photon2', 'coulomb', 'screened', 'plain'].includes(t)) {
            const kind = t.replace('2', '');
            const c = it.color ? EPP.color(it.color) : (kind === 'fermion' ? col.fermion : kind === 'phonon' ? col.phonon : (kind === 'photon' || kind === 'coulomb' || kind === 'screened') ? col.photon : col.ink);
            const base = curvePoints(it.a, it.b, it.bend || 0, 120);
            const L = arcLen(base);
            // dressed lines: double for fermions/Coulomb, heavier stroke for wavy/zigzag lines
            const dressed = t.endsWith('2') || kind === 'screened';
            const dbl = dressed && (kind === 'fermion' || kind === 'screened' || kind === 'plain');
            const strokes = dbl ? [-2.2, 2.2] : [0];
            for (const off of strokes) {
              let pts = off ? offsetPath(base, () => off) : base;
              let d, attrs = { fill: 'none', stroke: c, 'stroke-width': it.width || (dressed && !dbl ? 3.2 : 1.8), 'stroke-linecap': 'round', 'stroke-linejoin': 'round' };
              if (kind === 'photon') { const nw = Math.max(2, Math.round(L / (it.wl || 9))); pts = offsetPath(pts, s => (it.amp || 3.6) * Math.sin(2 * Math.PI * nw * s)); d = pathFrom(pts); }
              else if (kind === 'phonon') { const nz = Math.max(2, Math.round(L / (it.wl || 7))); const zig = curvePoints(it.a, it.b, it.bend || 0, nz * 2); const zz = (off ? offsetPath(zig, () => off) : zig).map((p, i) => { const o = (i === 0 || i === zig.length - 1) ? 0 : (i % 2 ? 1 : -1) * (it.amp || 3.6); return [p[0] - p[3] * o, p[1] + p[2] * o]; }); d = 'M' + zz.map(p => p[0].toFixed(2) + ' ' + p[1].toFixed(2)).join(' L'); }
              else { d = pathFrom(pts); if (kind === 'coulomb' || kind === 'screened') attrs['stroke-dasharray'] = '5 4'; }
              attrs.d = d; svgEl('path', attrs, svg);
            }
            if (kind === 'fermion' && it.arrow !== false) {
              const mid = base[Math.floor(base.length / 2)], ang = Math.atan2(mid[3], mid[2]) * 180 / Math.PI;
              const g = svgEl('g', { transform: `translate(${mid[0]},${mid[1]}) rotate(${ang})` }, svg);
              svgEl('path', { d: 'M5 0 L-4 -4.5 L-2 0 L-4 4.5 Z', fill: c }, g);
            }
            if (it.label) {
              // label on the screen-left of the direction a->b by default (above a left-to-right line); lpos:'right' (or '+') puts it on the other side
              const mid = base[Math.floor(base.length / 2)], dist = it.ldist ?? (dressed ? 14 : 12), sgn = (it.lpos === '+' || it.lpos === 'right') ? -1 : 1;
              const lx = mid[0] + mid[3] * dist * sgn, ly = mid[1] - mid[2] * dist * sgn;
              const lc = it.color ? EPP.color(it.color) : (kind === 'fermion' ? tcol.fermion : kind === 'phonon' ? tcol.phonon : (kind === 'photon' || kind === 'coulomb' || kind === 'screened') ? tcol.photon : col.ink);
              const te = svgEl('text', { x: lx, y: ly + 4, 'text-anchor': 'middle', class: 'fd-lab', fill: lc }, svg); richSVG(te, it.label);
            }
          } else if (t === 'vertex') {
            svgEl('circle', { cx: it.at[0], cy: it.at[1], r: it.r || 3, fill: it.color ? EPP.color(it.color) : col.ink }, svg);
          } else if (t === 'blob') {
            const c = it.color ? EPP.color(it.color) : col.ink;
            svgEl('circle', { cx: it.at[0], cy: it.at[1], r: it.r || 14, fill: EPP.color('surface') }, svg); // opaque base hides lines underneath
            svgEl('circle', { cx: it.at[0], cy: it.at[1], r: it.r || 14, fill: EPP.alpha(it.color || 'ink', 0.12), stroke: c, 'stroke-width': 1.6 }, svg);
            if (it.label) { const te = svgEl('text', { x: it.at[0], y: it.at[1] + 5, 'text-anchor': 'middle', class: 'fd-lab', fill: c, 'font-style': 'normal' }, svg); richSVG(te, it.label); }
          } else if (t === 'text') {
            const te = svgEl('text', { x: it.at[0], y: it.at[1], 'text-anchor': it.anchor || 'middle', class: 'fd-lab', fill: it.color ? EPP.color(it.color) : col.ink, 'font-size': it.size || 14 }, svg);
            if (it.upright) te.setAttribute('font-style', 'normal');
            richSVG(te, it.text);
          }
        }
      };
      paint();
      EPP.onTheme(paint);
      el.appendChild(svg);
      return svg;
    }
  };

  // ---------- small DOM helpers ----------
  EPP.$ = (sel, rootEl = document) => rootEl.querySelector(sel);
  EPP.$$ = (sel, rootEl = document) => [...rootEl.querySelectorAll(sel)];
  window.EPP = EPP;
})();
