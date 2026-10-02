"""One-off: add the 'pixel variable' colouring (pixattr/ files from build_pixel_attrs.py)
to fire_pcb_map.html in place. Every replacement must match exactly once."""
import sys

path = "fire_pcb_map.html"
t = open(path, "r", encoding="utf-8").read()
if "PIX_VAR_IX" in t:
    sys.exit("already patched")


def sub(old, new):
    global t
    n = t.count(old)
    assert n == 1, f"{n} matches for: {old[:70]!r}"
    t = t.replace(old, new)


# ---- 1. globals: loader + ramps (before BurnLayer) ---------------------------------
sub("const BURN_MODE_WORDS = {", r'''// ---- pixel variables (pixattr/, built by build_pixel_attrs.py from pixels.parquet) --
// Every burn-grid cell can carry V uint8 attributes (255 = none), in the SAME order as the
// grid's di/dj/dt, so a variable is just another colour source for the same canvas.
// Files are per year and fetched only when a variable is actually picked: *_pcb (the
// pyroCb-linked whole-fire grids + the +/-24 h windows) and *_all (background fires).
let PIX_VAR_IX = -1;                       // index into PIX_META.vars, -1 = off
let PIX_META = null, _pixMetaP = null;
const _pixFiles = {};                      // "2021_pcb" -> Promise<{idx, raw}>
const PIX_NB = 32, PIX_NODATA_COL = "#9a9a9a";
const PIX_SEQ = ["#440154","#414487","#2a788e","#22a884","#7ad151","#fde725"];   // viridis
const PIX_DIV = ["#2166ac","#67a9cf","#f7f7f7","#ef8a62","#b2182b"];
const PIX_CAT = ["#1b5e20","#9ccc3c","#f2b705","#6aa84f","#c27c3a",
                 "#d9c07a","#bdbdbd","#3a7ca5","#d95f9e","#4b4b4b"];
const _pixLut = {};
function pixMeta(){
  if (!_pixMetaP) _pixMetaP = fetch("pixattr/meta.json")
    .then(r => { if (!r.ok) throw new Error("pixattr/meta.json " + r.status); return r.json(); })
    .then(m => (PIX_META = m)).catch(e => { _pixMetaP = null; throw e; });
  return _pixMetaP;
}
async function _gz(url){
  const r = await fetch(url); if (!r.ok) throw new Error(url + " " + r.status);
  return new Uint8Array(await new Response(r.body.pipeThrough(new DecompressionStream("gzip"))).arrayBuffer());
}
function pixFile(y, kind){
  const k = y + "_" + kind;
  if (!_pixFiles[k]) _pixFiles[k] = Promise.all([
    fetch(`pixattr/${k}.json`).then(r => { if (!r.ok) throw new Error(k + ".json " + r.status); return r.json(); }),
    _gz(`pixattr/${k}.bin.gz`)]).then(([idx, raw]) => ({idx, raw}))
    .catch(e => { delete _pixFiles[k]; throw e; });
  return _pixFiles[k];
}
// -> Uint8Array [V*n] for this grid payload (tagged with _y/_sec/_key/_all at load), or null
async function pixBlock(payload){
  if (!payload || payload._y == null) return null;
  const [meta, f] = await Promise.all([pixMeta(), pixFile(payload._y, payload._all ? "all" : "pcb")]);
  const e = (f.idx[payload._sec] || {})[payload._key];
  if (!e || e[1] !== payload.di.length) return null;
  const V = meta.vars.length;
  return f.raw.subarray(e[0]*V, (e[0] + e[1])*V);
}
// 256-entry value -> colour table for variable vi (continuous vars quantised to PIX_NB bands)
function pixLut(vi){
  if (_pixLut[vi]) return _pixLut[vi];
  const m = PIX_META.vars[vi], lut = new Array(256).fill(PIX_NODATA_COL);
  if (m.kind === "cat") for (let v = 0; v <= m.hi; v++) lut[v] = PIX_CAT[v % PIX_CAT.length];
  else {
    const cols = _bandArray(m.kind === "div" ? PIX_DIV : PIX_SEQ, PIX_NB);
    for (let v = 0; v < 255; v++) lut[v] = cols[Math.min(PIX_NB-1, (v*PIX_NB/255)|0)];
  }
  return (_pixLut[vi] = lut);
}
const pixUnscale = (m, f) => { const v = m.lo + f*(m.hi - m.lo); return m.kind === "log" ? Math.pow(10, v) : v; };
const pixFmt = v => Math.abs(v) >= 100 ? v.toFixed(0) : Math.abs(v) >= 10 ? v.toFixed(1) : String(+v.toFixed(2));
function pixKeyHTML(vi){
  const m = PIX_META.vars[vi];
  if (m.kind === "cat")
    return `Colour = ${esc(m.label)}: ` + m.labels.map((l, i) =>
      `<span style="color:${PIX_CAT[i % PIX_CAT.length]};text-shadow:0 0 1px #999">&#9632;</span> ${esc(l)}`).join(" &middot; ");
  const stops = (m.kind === "div" ? PIX_DIV : PIX_SEQ).join(",");
  const u = m.unit ? " " + esc(m.unit) : "";
  return `Colour = ${esc(m.label)}${m.kind === "log" ? " (log scale)" : ""}: `
    + `<span style="display:inline-block;width:90px;height:8px;vertical-align:middle;border:1px solid #bbb;`
    + `background:linear-gradient(to right,${stops})"></span> `
    + `${pixFmt(pixUnscale(m, 0))} &rarr; ${pixFmt(pixUnscale(m, 1))}${u}`
    + ` <span style="color:${PIX_NODATA_COL}">&#9632;</span> no value`;
}
// Attach (or clear) the attribute block of the grid the canvas is showing, then repaint.
async function applyPixAttr(payload, g){
  if (!burn) return;
  if (PIX_VAR_IX < 0){ if (burn._attr){ burn.setAttr(null, 0); burn.repaint(); } return; }
  try {
    const blk = await pixBlock(payload);
    if (burn._g !== g) return;               // the canvas moved on while we were fetching
    if (!blk){ PIX_VAR_IX = -1; burn.setAttr(null, 0); burn.repaint();
      showNotice("info", "No pixel attributes for this grid."); return; }
    burn.setAttr(blk, g.n); burn.repaint();
  } catch (e) {
    console.warn(e); PIX_VAR_IX = -1;
    if (burn._g === g){ burn.setAttr(null, 0); burn.repaint(); }
    showNotice("warn", "Could not load pixel attributes (" + esc(e.message) + ").");
  }
}
const BURN_MODE_WORDS = {''')

# ---- 2. BurnLayer ------------------------------------------------------------------
sub("this._absT0=0; this._sp=null;         // grid t=0 in ms UTC; decoded spread series",
    "this._absT0=0; this._sp=null;         // grid t=0 in ms UTC; decoded spread series\n"
    "    this._attr=null; this._attrN=0;       // pixel-variable block for the current grid")
sub("this._absT0 = absT0||0; this._sp = decodeSpread(spread);\n    this._range = [0,0]; this._drawn = null; this._blank();",
    "this._absT0 = absT0||0; this._sp = decodeSpread(spread); this._attr = null; this._attrN = 0;\n"
    "    this._range = [0,0]; this._drawn = null; this._blank();")
sub("  clear(){ this._g = null; this._sp = null;",
    "  setAttr(blk, n){ this._attr = blk; this._attrN = n; },\n  clear(){ this._g = null; this._sp = null; this._attr = null;")
sub("    const mode = BURN_COLOR_MODE, sp = this._sp;\n    const bands = (mode",
    """    if (PIX_VAR_IX >= 0){
      // pixel-variable mode: colour = this cell's value of the chosen variable. Until the
      // block has arrived draw nothing (applyPixAttr repaints) rather than flash burn-age.
      if (this._attr && this._attrN === g.n) this._paintAttr(i0, i1, S, PX, PY, w, x0, x1, y0, y1);
      return;
    }
    const mode = BURN_COLOR_MODE, sp = this._sp;
    const bands = (mode""")
sub("  gridBounds(){\n    const g = this._g;\n    if (!g || !g.n) return null;\n    let a=g.ix[0]",
    """  _paintAttr(i0, i1, S, PX, PY, w, x0, x1, y0, y1){
    const g = this._g, ctx = this._ctx, ix = g.ix, iy = g.iy;
    const a = this._attr, off = PIX_VAR_IX * g.n, lut = pixLut(PIX_VAR_IX);
    ctx.globalAlpha = 0.85;
    let cur = null;
    for (let k = i0; k < i1; k++){
      const x = ix[k]*S + PX, y = PY - (iy[k]+1)*S;
      if (x < x0 || y < y0 || x > x1 || y > y1) continue;
      const c = lut[a[off + k]];
      if (c !== cur){ cur = c; ctx.fillStyle = c; }
      ctx.fillRect(x, y, w, w);
    }
    ctx.globalAlpha = 1;
  },
  gridBounds(){
    const g = this._g;
    if (!g || !g.n) return null;
    let a=g.ix[0]""")

# ---- 3. showGrid attaches the block when the canvas switches grid -------------------
sub("  if (burn._g !== g) burn.setGrid(g, span, rampT0, o.absT0Ms||0, spr);\n  else burn.setContext",
    "  if (burn._g !== g){ burn.setGrid(g, span, rampT0, o.absT0Ms||0, spr);"
    " if (PIX_VAR_IX >= 0) applyPixAttr(payload, g); }\n  else burn.setContext")

# ---- 4. tag every grid payload with where its attributes live ----------------------
sub("  Object.assign(byFireSim, sim);         // mergid  -> whole-fire series",
    """  Object.entries(pix || {}).forEach(([k, g]) => { g._y = +y; g._sec = "pix"; g._key = k; });
  Object.entries(sim || {}).forEach(([k, v]) => { if (v && v.g){ v.g._y = +y; v.g._sec = "sim"; v.g._key = k; } });
  Object.assign(byFireSim, sim);         // mergid  -> whole-fire series""")
sub("  Object.assign(byFireSim, sim || {});",
    """  Object.entries(sim || {}).forEach(([k, v]) => { if (v && v.g){ v.g._y = +y; v.g._sec = "sim"; v.g._key = k; v.g._all = true; } });
  Object.assign(byFireSim, sim || {});""")

# ---- 5. panel UI ------------------------------------------------------------------
sub("""+`<button data-m="frp95"${bySpread[p.mergid]?"":" disabled"}>FRP95</button></span></label>`""",
    """+`<button data-m="frp95"${bySpread[p.mergid]?"":" disabled"}>FRP95</button></span></label>`
      +`<label>pixel variable <select id="pvSel" title="colour each 375 m cell by an attribute from pixels.parquet">`
        +`<option value="">${PIX_META ? "none - use colour above" : "loading..."}</option></select></label>`""")
sub("""    $("#tierkey").innerHTML = (BURN_COLOR_MODE === "age" ? AGE_KEY : SPD_KEY(BURN_COLOR_MODE))""",
    """    $("#tierkey").innerHTML = (PIX_VAR_IX >= 0 && PIX_META ? pixKeyHTML(PIX_VAR_IX)
      : BURN_COLOR_MODE === "age" ? AGE_KEY : SPD_KEY(BURN_COLOR_MODE))""")
sub("""      BURN_COLOR_MODE = b.dataset.m;
      cmBtns.forEach(x => x.classList.toggle("on", x === b));
      setColourKey();
      burn.setContext(spanOf(), t0Of(), absT0Of(), bySpread[p.mergid]); burn.repaint();
    };
  });""", """      BURN_COLOR_MODE = b.dataset.m;
      cmBtns.forEach(x => x.classList.toggle("on", x === b));
      PIX_VAR_IX = -1; const pv0 = $("#pvSel"); if (pv0) pv0.value = "";
      setColourKey();
      burn.setContext(spanOf(), t0Of(), absT0Of(), bySpread[p.mergid]); burn.repaint();
    };
  });
  // pixel variable: a dropdown over every attribute shipped in pixattr/. It overrides the
  // colour buttons above while set; picking a button (or "none") hands colour back.
  (async () => {
    const sel = $("#pvSel"); if (!sel) return;
    try { await pixMeta(); } catch (e) { sel.innerHTML = `<option value="">unavailable</option>`; sel.disabled = true; return; }
    const grp = {};
    PIX_META.vars.forEach((m, i) => (grp[m.group] = grp[m.group] || []).push(
      `<option value="${i}">${esc(m.label)}${m.unit ? " (" + esc(m.unit) + ")" : ""}</option>`));
    sel.innerHTML = `<option value="">none - use colour above</option>`
      + Object.entries(grp).map(([g, o]) => `<optgroup label="${esc(g)}">${o.join("")}</optgroup>`).join("");
    sel.value = PIX_VAR_IX >= 0 ? String(PIX_VAR_IX) : "";
    setColourKey();
    sel.onchange = () => {
      PIX_VAR_IX = sel.value === "" ? -1 : +sel.value;
      if (PIX_VAR_IX >= 0) cmBtns.forEach(x => x.classList.remove("on"));
      else cmBtns.forEach(x => x.classList.toggle("on", x.dataset.m === BURN_COLOR_MODE));
      setColourKey();
      const gp = gridOf();
      applyPixAttr(gp, decodeGrid(gp));
    };
  })();""")

open(path, "w", encoding="utf-8").write(t)
print("patched", len(t))
