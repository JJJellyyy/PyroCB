"""One-off, runs after patch_map_pixattr.py: colour ALL polygons by a pixel variable
(fire-level aggregates, pixattr/fires_{year}.json) and add a toggle that colours the
non-pyroCb ("All fire activity") fires instead of leaving them grey."""
import sys

path = "fire_pcb_map.html"
t = open(path, "r", encoding="utf-8").read()
if "PAINT_ALL" in t:
    sys.exit("already patched")


def sub(old, new):
    global t
    n = t.count(old)
    assert n == 1, f"{n} matches for: {old[:70]!r}"
    t = t.replace(old, new)


# ---- 1. fire-level loader + hue ----------------------------------------------------
sub("// Attach (or clear) the attribute block of the grid the canvas is showing, then repaint.",
    r'''// ---- whole-polygon colouring by variable (MODE === "var") --------------------------
// pixattr/fires_{year}.json: {mergid: [V codes]}, one value per fire per variable, quantised
// on the FIRE-level scale (meta flo/fhi) - separate from the per-cell scale above.
let PAINT_ALL = false;                     // also colour fires that never made a pyroCb
let VAR_MODE_IX = -1;                      // variable colouring polygons when MODE === "var"
const fireAggData = {}, _fireAggP = {};
function ensureFireAgg(y){
  if (!_fireAggP[y]) _fireAggP[y] = pixMeta()
    .then(() => fetch(`pixattr/fires_${y}.json`))
    .then(r => { if (!r.ok) throw new Error(`fires_${y}.json ${r.status}`); return r.json(); })
    .then(d => (fireAggData[y] = d))
    .catch(e => { delete _fireAggP[y]; throw e; });
  return _fireAggP[y];
}
function varHue(y, p){
  const d = fireAggData[y];
  if (!d){                                  // year not loaded yet: grey now, recolour on arrival
    ensureFireAgg(y).then(() => forceFilter()).catch(e => console.warn(e));
    return PIX_NODATA_COL;
  }
  const r = d[String(p.mergid)];
  return pixLut(VAR_MODE_IX, null)[r ? r[VAR_MODE_IX] : 255];
}
function varLegendHTML(){
  const m = PIX_META.vars[VAR_MODE_IX];
  const nod = `<span style="color:${PIX_NODATA_COL}">&#9632;</span> no value`;
  if (m.kind === "cat")
    return `<div class="sub" style="margin:6px 0 3px">Colour = ${esc(m.label)} (fire's dominant class)</div>`
      + m.labels.map((l, i) => `<div class="row"><i class="sq" style="background:${PIX_CAT[i % PIX_CAT.length]}"></i>${esc(l)}</div>`).join("")
      + `<div class="sub">${nod}</div>`;
  const fm = Object.assign({}, m, {lo:m.flo, hi:m.fhi});
  const stops = (m.kind === "div" ? PIX_DIV : PIX_SEQ).join(",");
  const u = m.unit ? " " + esc(m.unit) : "";
  return `<div class="sub" style="margin:6px 0 3px">Colour = ${esc(m.label)}${m.kind === "log" ? " (log)" : ""},`
    + ` mean over each fire's pixels</div>`
    + `<div id="bar" style="background:linear-gradient(90deg,${stops})"></div>`
    + `<div id="barlab"><span>${pixFmt(pixUnscale(fm, 0))}</span><span>${pixFmt(pixUnscale(fm, 0.5))}</span>`
    + `<span>${pixFmt(pixUnscale(fm, 1))}${u}</span></div><div class="sub">${nod}</div>`;
}
// Attach (or clear) the attribute block of the grid the canvas is showing, then repaint.''')

# ---- 2. hueOf / bndStyle ----------------------------------------------------------
sub("""const hueOf = (y,f) => f.properties.pcb_linked === false ? ALLFIRE_HUE
  : MODE==="year" ? YCOL[y] : doyColor(f.properties.doy);""",
    """const hueOf = (y,f) => (f.properties.pcb_linked === false && !PAINT_ALL) ? ALLFIRE_HUE
  : MODE==="year" ? YCOL[y] : MODE==="var" ? varHue(y, f.properties) : doyColor(f.properties.doy);""")
sub("""  const hue = MODE==="year" ? YCOL[y] : doyColor(p.doy);
  if (OUTLINE.bnd)""",
    """  const hue = MODE==="year" ? YCOL[y] : MODE==="var" ? varHue(y, p) : doyColor(p.doy);
  if (OUTLINE.bnd)""")

# ---- 3. all-fire layers join recolour(); their ignitions follow hueOf ---------------
sub("""      color:"#5d5d5d", weight:1, fillColor:ALLFIRE_HUE, fillOpacity:0.95}),""",
    """      color:"#5d5d5d", weight:1, fillColor:hueOf(y, ft), fillOpacity:0.95}),""")
sub("""      else if (l.getElement && l.getElement()) l.getElement().firstChild.style.background = hueOf(y,l.feature);
    }));
  }
  drawLegend();""",
    """      else if (l.getElement && l.getElement()) l.getElement().firstChild.style.background = hueOf(y,l.feature);
    }));
  }
  for (const y of Object.keys(allFires)){           // the "All fire activity" layers
    const a = allFires[y];
    a.per.setStyle(ft => perStyle(y, ft));
    a.ign.eachLayer(l => { if (l.setStyle) l.setStyle({fillColor: hueOf(y, l.feature)}); });
  }
  drawLegend();""")

# ---- 4. legend ---------------------------------------------------------------------
sub("""  } else {
    const stops = Array.from({length:24},(_,i)=>ramp(i/23)).join(",");""",
    """  } else if (MODE==="var" && PIX_META && VAR_MODE_IX >= 0){
    legendBody.innerHTML = common + varLegendHTML();
  } else {
    const stops = Array.from({length:24},(_,i)=>ramp(i/23)).join(",");""")

# ---- 5. main-panel controls --------------------------------------------------------
sub("""<div class="seg"><button id="mYear" class="on">Year</button><button id="mDoy">Day of year</button></div>""",
    """<div class="seg"><button id="mYear" class="on">Year</button><button id="mDoy">Day of year</button><button id="mVar">Variable</button></div>
    <div id="varBox" hidden style="margin-top:4px"><select id="mVarSel" style="width:100%"></select></div>
    <label style="display:block;margin:4px 0 0;font-size:12px" title="Fires without a pyroCb are normally grey (All fire activity layer)"><input type="checkbox" id="paintAll"> colour non-pyroCb fires too</label>""")
sub("""  const my=d.querySelector("#mYear"), mdy=d.querySelector("#mDoy");
  my.onclick  = ()=>{MODE="year"; my.classList.add("on"); mdy.classList.remove("on"); recolour();};
  mdy.onclick = ()=>{MODE="doy";  mdy.classList.add("on"); my.classList.remove("on"); recolour();};""",
    """  const my=d.querySelector("#mYear"), mdy=d.querySelector("#mDoy"), mvr=d.querySelector("#mVar");
  const varBox=d.querySelector("#varBox"), varSel=d.querySelector("#mVarSel"), paintAll=d.querySelector("#paintAll");
  const setMode = m => { MODE = m;
    [[my,"year"],[mdy,"doy"],[mvr,"var"]].forEach(([b,k]) => b.classList.toggle("on", k===m));
    varBox.hidden = m !== "var"; forceFilter(); };
  my.onclick  = () => setMode("year");
  mdy.onclick = () => setMode("doy");
  // Variable mode: polygons coloured by the fire-level mean of a pixels.parquet variable.
  // Load every loaded year's aggregates BEFORE switching, so nothing flashes grey.
  const loadAggs = () => Promise.all(Object.keys(cache).map(ensureFireAgg));
  let varTouched = false;
  mvr.onclick = async () => {
    try { await pixMeta(); await loadAggs(); }
    catch (e){ showNotice("warn", "Pixel variables are unavailable: " + esc(e.message)
      + ". The pixattr/ folder must be deployed next to the map and the map opened over http(s), not file://."); return; }
    if (!varSel.options.length){
      const grp = {};
      PIX_META.vars.forEach((m, i) => (grp[m.group] = grp[m.group] || []).push(
        `<option value="${i}">${esc(m.label)}${m.unit ? " (" + esc(m.unit) + ")" : ""}</option>`));
      varSel.innerHTML = Object.entries(grp).map(([g, o]) => `<optgroup label="${esc(g)}">${o.join("")}</optgroup>`).join("");
    }
    if (VAR_MODE_IX < 0){ VAR_MODE_IX = 0; varSel.value = "0"; }
    if (!varTouched){ varTouched = true; PAINT_ALL = true; paintAll.checked = true; }
    setMode("var");
  };
  varSel.onchange = async () => { VAR_MODE_IX = +varSel.value; await loadAggs().catch(console.warn); forceFilter(); };
  paintAll.onchange = e => { PAINT_ALL = e.target.checked; forceFilter(); };""")

# ---- 6. clearer failure text in the panel dropdown ----------------------------------
sub("""catch (e) { sel.innerHTML = `<option value="">unavailable</option>`; sel.disabled = true; return; }""",
    """catch (e) { sel.innerHTML = `<option value="">unavailable - pixattr/ not found</option>`; sel.title = e.message; sel.disabled = true; return; }""")

open(path, "w", encoding="utf-8").write(t)
print("patched", len(t))
