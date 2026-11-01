// Measure display equations that overflow their box at phone width (390 px), with all details boxes opened.
// Usage: node eqwidth.js <page.html> <out.json> [--min 40]   (writes a list sorted by section and overflow)
const { chromium } = require('playwright');
const fs = require('fs'), path = require('path');
const [, , pagePath, outPath, ...rest] = process.argv;
const min = rest.indexOf('--min') >= 0 ? +rest[rest.indexOf('--min') + 1] : 1;
const MJ = path.resolve(__dirname, '..', '..', 'vendor', 'mathjax', 'es5');
let frag = fs.readFileSync(pagePath, 'utf8').replace("loader: { load: ['ui/lazy'] },", '');
const doc = `<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><style>:root{color-scheme:light}body{margin:0}[hidden]{display:none!important}</style></head><body>${frag}</body></html>`;
(async () => {
  const browser = await chromium.launch();
  const ctx = await browser.newContext({ viewport: { width: 390, height: 844 }, deviceScaleFactor: 1 });
  const page = await ctx.newPage();
  await page.route('https://cdnjs.cloudflare.com/ajax/libs/mathjax/3.2.2/es5/**', r => { const f = path.join(MJ, r.request().url().split('/es5/')[1].split('?')[0]); fs.existsSync(f) ? r.fulfill({ path: f, contentType: 'application/javascript' }) : r.fulfill({ status: 404, body: '' }); });
  await page.route('https://fonts.googleapis.com/**', r => r.fulfill({ status: 200, contentType: 'text/css', body: '' }));
  await page.setContent(doc, { waitUntil: 'load', timeout: 240000 });
  await page.waitForFunction(() => window.MathJax && MathJax.startup && MathJax.startup.promise, null, { timeout: 180000 });
  await page.evaluate(() => MathJax.startup.promise);
  await page.evaluate(() => document.querySelectorAll('details').forEach(d => d.open = true));
  await page.waitForTimeout(1500);
  const res = await page.evaluate((min) => {
    const tex = new Map(); for (const m of MathJax.startup.document.math) if (m.typesetRoot) tex.set(m.typesetRoot, m.math);
    const out = [];
    document.querySelectorAll('.eq, mjx-container[display="true"]').forEach(el => {
      const box = el.classList.contains('eq') ? el : el.parentElement;
      if (!box || (el.tagName.toLowerCase() === 'mjx-container' && box.closest('.eq'))) return;
      const over = box.scrollWidth - box.clientWidth;
      if (over < min) return;
      const c = el.tagName.toLowerCase() === 'mjx-container' ? el : el.querySelector('mjx-container');
      const t = c ? (tex.get(c) || '') : '';
      const tag = (t.match(/\\tag\{([^}]*)\}/) || [])[1] || '';
      out.push({ section: el.closest('section')?.id || '?', id: el.id || box.id || '', tag, over, tex: t.replace(/\s+/g, ' ').slice(0, 160) });
    });
    return out;
  }, min);
  res.sort((a, b) => a.section.localeCompare(b.section) || b.over - a.over);
  fs.writeFileSync(outPath, JSON.stringify(res, null, 1));
  const by = {}; res.forEach(r => by[r.section] = (by[r.section] || 0) + 1);
  console.log(`${res.length} display equations overflow by >= ${min} px at 390 px:`, JSON.stringify(by));
  await browser.close();
})();
