// Screenshot every figure-like element of the assembled page (eager MathJax).
// Usage: node shoot.js <page.html> <outdir> [--only A0,B1 (filename prefixes of sections)]
const { chromium } = require('playwright');
const fs = require('fs'), path = require('path');
const S = path.resolve(__dirname, '..');
const [, , pagePath, outDir, ...rest] = process.argv;
const only = rest.indexOf('--only') >= 0 ? rest[rest.indexOf('--only') + 1].split(',') : null;
const MJ = path.resolve(__dirname, '..', '..', 'vendor', 'mathjax', 'es5');
fs.mkdirSync(outDir, { recursive: true });
let frag = fs.readFileSync(pagePath, 'utf8').replace("loader: { load: ['ui/lazy'] },", '');
const doc = `<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover"><style>:root{color-scheme:light}body{margin:0}[hidden]{display:none!important}</style></head><body>${frag}</body></html>`;
(async () => {
  const browser = await chromium.launch();
  const index = [];
  for (const cfg of [{ name: 'desk', w: 1280, h: 900, scheme: 'light' }, { name: 'phone', w: 390, h: 844, scheme: 'dark' }]) {
    const ctx = await browser.newContext({ viewport: { width: cfg.w, height: cfg.h }, colorScheme: cfg.scheme });
    const page = await ctx.newPage();
    await page.route('https://cdnjs.cloudflare.com/ajax/libs/mathjax/3.2.2/es5/**', r => { const f = path.join(MJ, r.request().url().split('/es5/')[1].split('?')[0]); fs.existsSync(f) ? r.fulfill({ path: f, contentType: 'application/javascript' }) : r.fulfill({ status: 404, body: '' }); });
    await page.route('https://fonts.googleapis.com/**', r => r.fulfill({ status: 200, contentType: 'text/css', body: '' }));
    await page.setContent(doc, { waitUntil: 'load' });
    await page.waitForFunction(() => window.MathJax && MathJax.startup && MathJax.startup.promise, null, { timeout: 30000 });
    await page.evaluate(() => MathJax.startup.promise);
    await page.waitForTimeout(800);
    const items = await page.evaluate(() => {
      const out = []; let k = 0;
      document.querySelectorAll('.hero, figure.viz, div.fd, div.tbl, .part-head').forEach(el => {
        const sec = el.closest('section.chapter');
        const id = el.id || ('el' + (k++));
        el.setAttribute('data-shot', id);
        out.push({ shot: id, section: sec ? sec.id : (el.classList.contains('hero') ? 'hero' : 'part-head'), kind: el.className.split(' ')[0] });
      });
      return out;
    });
    for (const it of items) {
      if (only && !only.some(p => it.section.startsWith(p))) continue;
      const el = await page.$(`[data-shot="${it.shot}"]`);
      if (!el || !(await el.isVisible())) { index.push({ ...it, cfg: cfg.name, skipped: 'hidden' }); continue; }
      try { await el.scrollIntoViewIfNeeded({ timeout: 3000 }); } catch (e) { index.push({ ...it, cfg: cfg.name, skipped: 'noscroll' }); continue; }
      await page.waitForTimeout(350);
      const file = `${cfg.name}--${it.section}--${it.shot}.png`;
      try { await el.screenshot({ path: path.join(outDir, file) }); index.push({ ...it, cfg: cfg.name, file }); } catch (e) { index.push({ ...it, cfg: cfg.name, error: e.message }); }
    }
    await ctx.close();
  }
  await browser.close();
  fs.writeFileSync(path.join(outDir, 'index.json'), JSON.stringify(index, null, 1));
  console.log(index.length + ' screenshots');
})();
