// Usage: node check.js <page.html> <outdir> [--shots] [--sections id1,id2]
// Serves vendored MathJax for the cdnjs URL, wraps the fragment in a doc skeleton like the artifact host,
// reports console errors, MathJax errors, horizontal overflow, empty canvases; optional screenshots.
const { chromium } = require('playwright');
const fs = require('fs'), path = require('path');
const S = path.resolve(__dirname, '..');
const [, , pagePath, outDir = path.join(S, 'shots'), ...rest] = process.argv;
const shots = rest.includes('--shots');
const secArg = rest.indexOf('--sections') >= 0 ? rest[rest.indexOf('--sections') + 1].split(',') : null;
const MJ = path.resolve(__dirname, '..', '..', 'vendor', 'mathjax', 'es5');
fs.mkdirSync(outDir, { recursive: true });
let frag = fs.readFileSync(pagePath, 'utf8');
if (rest.includes('--eager')) frag = frag.replace("loader: { load: ['ui/lazy'] },", '');
const doc = `<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover"><style>:root{color-scheme:light;padding-top:env(safe-area-inset-top);padding-bottom:env(safe-area-inset-bottom)}body{margin:0;font:14px system-ui}img{max-width:100%}[hidden]{display:none!important}</style></head><body>${frag}</body></html>`;
(async () => {
  const browser = await chromium.launch();
  const report = { errors: [], warnings: [], mathErrors: 0, overflow: {}, blankCanvases: [], timings: {} };
  for (const cfg of [{ name: 'desktop', w: 1280, h: 900, scheme: 'light' }, { name: 'phone', w: 390, h: 844, scheme: 'dark' }]) {
    const ctx = await browser.newContext({ viewport: { width: cfg.w, height: cfg.h }, colorScheme: cfg.scheme, deviceScaleFactor: 1 });
    const page = await ctx.newPage();
    page.on('console', m => { if (m.type() === 'error') report.errors.push(`[${cfg.name}] ${m.text()}`); else if (m.type() === 'warning') report.warnings.push(`[${cfg.name}] ${m.text()}`); });
    page.on('pageerror', e => report.errors.push(`[${cfg.name}] pageerror: ${e.message}`));
    await page.route('https://cdnjs.cloudflare.com/ajax/libs/mathjax/3.2.2/es5/**', r => {
      const rel = r.request().url().split('/es5/')[1].split('?')[0];
      const f = path.join(MJ, rel);
      if (fs.existsSync(f)) r.fulfill({ path: f, contentType: 'application/javascript' }); else r.fulfill({ status: 404, body: '' });
    });
    await page.route('https://fonts.googleapis.com/**', r => r.fulfill({ status: 200, contentType: 'text/css', body: '' }));
    await page.route('https://fonts.gstatic.com/**', r => r.fulfill({ status: 404, body: '' }));
    const t0 = Date.now();
    await page.setContent(doc, { waitUntil: 'load' });
    try { await page.waitForFunction(() => window.MathJax && MathJax.startup && MathJax.startup.promise, null, { timeout: 20000 }); await page.evaluate(() => MathJax.startup.promise); } catch (e) { report.errors.push(`[${cfg.name}] MathJax did not start: ${e.message}`); }
    report.timings[cfg.name] = Date.now() - t0;
    await page.waitForTimeout(800);
    const r = await page.evaluate(() => {
      const de = document.documentElement;
      const over = de.scrollWidth - de.clientWidth;
      const wide = [];
      if (over > 0) {
        document.querySelectorAll('main *').forEach(el => { const b = el.getBoundingClientRect(); if (b.right > de.clientWidth + 1 && b.width > 0) { let p = el.parentElement, scroll = false; while (p) { const cs = getComputedStyle(p); if (/(auto|scroll|hidden)/.test(cs.overflowX)) { scroll = true; break; } p = p.parentElement; } if (!scroll) wide.push((el.closest('section')?.id || '?') + ' ' + el.tagName + '.' + el.className + ' right=' + Math.round(b.right)); } });
      }
      const mErr = [...document.querySelectorAll('mjx-merror, [data-mjx-error]')].map(e => (e.closest('section')?.id || '?') + ': ' + (e.getAttribute('data-mjx-error') || e.textContent).slice(0, 120));
      const blank = [];
      document.querySelectorAll('canvas').forEach(c => { try { const g = c.getContext('2d'); const d = g.getImageData(0, 0, c.width, c.height).data; let n = 0; for (let i = 3; i < d.length; i += 4 * 97) if (d[i] > 0) n++; if (n === 0) blank.push(c.closest('figure')?.id || c.closest('section')?.id || '?'); } catch (e) {} });
      const sections = [...document.querySelectorAll('section.chapter')].map(s => s.id);
      return { over, wide: wide.slice(0, 15), mErr, blank, sections, h: de.scrollHeight };
    });
    report.overflow[cfg.name] = { px: r.over, offenders: r.wide };
    report.mathErrors += r.mErr.length; if (r.mErr.length) report.errors.push(...r.mErr.map(x => `[${cfg.name}] TeX: ${x}`));
    if (r.blank.length) report.blankCanvases.push(`[${cfg.name}] ${r.blank.join(', ')}`);
    report.sections = r.sections; report.height = report.height || {}; report.height[cfg.name] = r.h;
    if (shots) {
      const ids = secArg || r.sections;
      for (const id of ids) {
        const el = await page.$(`[id="${id}"]`);
        if (!el) continue;
        await el.scrollIntoViewIfNeeded(); await page.waitForTimeout(900);
        await el.screenshot({ path: path.join(outDir, `${cfg.name}-${id}.png`) });
      }
      await page.evaluate(() => window.scrollTo(0, 0)); await page.waitForTimeout(200);
      await page.screenshot({ path: path.join(outDir, `${cfg.name}-top.png`) });
    }
    await ctx.close();
  }
  await browser.close();
  console.log(JSON.stringify(report, null, 1));
})();
