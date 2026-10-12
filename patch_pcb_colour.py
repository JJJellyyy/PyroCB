"""One-off patch (2026-10-05): colour pyroCb markers by injection-height variables.

Adds "Colour pyroCbs by" under "Colour fires by": Same as fires | Injection height | Tropopause height |
Injection - tropopause | Injection above fire ground | Share of troposphere above fire reached. Values come from
the pcb payload (max_inject_alt, trop_alt, trop_relative_alt, fire_elev_m, inj_above_fire, inj_frac_trop - see
build_active_pcb.py add_heights). Legend key, popup rows and info-panel rows included.
    python patch_pcb_colour.py   then   python map_tool.py join && node map_check.js
"""
from pathlib import Path

TAIL = Path(__file__).resolve().parent / "map_src" / "tail.html"
t = TAIL.read_text(encoding="utf-8")
if "PCB_CMODE" in t:
    raise SystemExit("already patched")


def sub(old, new, count=1):
    global t
    assert t.count(old) == count, (t.count(old), old[:90])
    t = t.replace(old, new)


# ---- state, scales and the hue function (right after hueOf) ----------------------------------------
sub('''const hueOf = (y,f) => (f.properties.pcb_linked === false && !PAINT_ALL) ? ALLFIRE_HUE
  : MODE==="year" ? YCOL[y] : MODE==="var" ? varHue(y, f.properties) : doyColor(f.properties.doy);''',
'''const hueOf = (y,f) => (f.properties.pcb_linked === false && !PAINT_ALL) ? ALLFIRE_HUE
  : MODE==="year" ? YCOL[y] : MODE==="var" ? varHue(y, f.properties) : doyColor(f.properties.doy);

// ---- pyroCb markers coloured by injection height (independent of the fire colouring) -------------
// Fixed scales (about the 2nd-98th percentile of the 2013-2023 inventory) so colours compare across years.
// Heights: viridis. Signed / ratio quantities: RdBu diverging about their neutral value (tropopause = 0 / 1).
const PCB_CVARS = {
  inj:     {key:"max_inject_alt",    label:"Injection height",               unit:"m ASL", lo:8000,  hi:15000, stops:RAMPS.index,
            def:"Maximum smoke injection altitude from the inventory (brightness-temperature based), metres above sea level."},
  trop:    {key:"trop_alt",          label:"Tropopause height",              unit:"m ASL", lo:8500,  hi:13700, stops:RAMPS.index,
            def:"MERRA-2 tropopause altitude at the event (inventory), metres above sea level."},
  rel:     {key:"trop_relative_alt", label:"Injection - tropopause",         unit:"m",     lo:-3500, hi:3500,  stops:RAMPS.delta, mid:"0 = at the tropopause",
            def:"Injection height minus tropopause height: below zero = stayed in the troposphere, above zero = reached the stratosphere."},
  above:   {key:"inj_above_fire",    label:"Injection above fire ground",    unit:"m",     lo:7500,  hi:14000, stops:RAMPS.index,
            def:"Injection height minus the fire's ground elevation (median Copernicus DEM elevation of the event's active pyroCb pixels; DEM at the plotted point when there are none): how high the plume rose above the fire."},
  frac:    {key:"inj_frac_trop",     label:"Share of troposphere reached",   unit:"",      lo:0.7,   hi:1.3,   stops:RAMPS.delta, mid:"1 = tropopause",
            def:"(injection - fire ground) / (tropopause - fire ground): the fraction of the tropospheric column above the fire that the plume crossed. 1 = reached the tropopause, > 1 = overshot into the stratosphere. Both heights relative to the fire's own altitude."},
};
let PCB_CMODE = "none";                         // "none" = same colour as the fire | one of PCB_CVARS
const PCB_NOVAL = "#c8c8c8";
function _rampAt(stops, f){
  f = Math.max(0, Math.min(1, f));
  const x = f*(stops.length-1), i = Math.min(stops.length-2, x|0), g = x-i;
  const a = hex2rgb(stops[i]), b = hex2rgb(stops[i+1]);
  return rgb2hex(a.map((v,q)=>v+(b[q]-v)*g));
}
function pcbHue(y, ft){
  if (PCB_CMODE === "none") return hueOf(y, ft);
  const v = PCB_CVARS[PCB_CMODE], x = ft.properties[v.key];
  return (x == null || !isFinite(x)) ? PCB_NOVAL : _rampAt(v.stops, (x - v.lo) / (v.hi - v.lo));
}
const fmtM = v => v == null ? "–" : Math.round(v).toLocaleString("en").replace(/,/g, " ") + " m";''')

# ---- markers use pcbHue ------------------------------------------------------------------------------
sub('''               html:`<div class="pcblbl" style="background:${hueOf(y,ft)};border:1.5px solid #111;''',
    '''               html:`<div class="pcblbl" style="background:${pcbHue(y,ft)};border:1.5px solid #111;''')
sub('''                   weight:1.5,fillColor:hueOf(y,ft),fillOpacity:1});},''',
    '''                   weight:1.5,fillColor:pcbHue(y,ft),fillOpacity:1});},''')
sub('''    ["pcb","pcbC"].forEach(k=>c[k].eachLayer(l=>{
      if (l.setStyle) l.setStyle({fillColor:hueOf(y,l.feature)});
      else if (l.getElement && l.getElement()) l.getElement().firstChild.style.background = hueOf(y,l.feature);''',
    '''    ["pcb","pcbC"].forEach(k=>c[k].eachLayer(l=>{
      if (l.setStyle) l.setStyle({fillColor:pcbHue(y,l.feature)});
      else if (l.getElement && l.getElement()) l.getElement().firstChild.style.background = pcbHue(y,l.feature);''')
# divIcon markers built while their element is not yet in the DOM keep the build-time colour: recolour on add
sub('''      l.on("click", () => DETAIL_MODE==="panel" ? openInfoPanel(p)
                                               : renderActiveWindow(p, 12, 12));}''',
    '''      l.on("click", () => DETAIL_MODE==="panel" ? openInfoPanel(p)
                                               : renderActiveWindow(p, 12, 12));
      l.on("add", () => { const el = l.getElement && l.getElement();
        if (el && el.firstChild) el.firstChild.style.background = pcbHue(y, ft); });}''')

# ---- popup rows --------------------------------------------------------------------------------------
sub('''          ["Max injection alt (m)",p.max_inject_alt],["Rel. tropopause (m)",p.trop_relative_alt],''',
    '''          ["Max injection alt (m)",p.max_inject_alt],["Tropopause (m)",p.trop_alt],
          ["Rel. tropopause (m)",p.trop_relative_alt],["Fire ground elevation (m)",p.fire_elev_m],
          ["Injection above fire (m)",p.inj_above_fire],["Share of troposphere reached",p.inj_frac_trop],''')

# ---- info panel row ----------------------------------------------------------------------------------
sub('''    ["New coordinate (matched)", p.np_lat!=null''',
    '''    ["Injection height", p.max_inject_alt!=null
      ? `${fmtM(p.max_inject_alt)} ASL &middot; tropopause ${fmtM(p.trop_alt)} &middot; `
        + `${p.trop_relative_alt>=0?"+":""}${fmtM(p.trop_relative_alt)} vs tropopause`
        + (p.inj_above_fire!=null ? `<br>${fmtM(p.inj_above_fire)} above the fire (ground ${fmtM(p.fire_elev_m)}) &middot; `
          + `${Math.round(p.inj_frac_trop*100)}% of the troposphere above the fire` : "")
      : "not in the inventory"],
    ["New coordinate (matched)", p.np_lat!=null''')

# ---- control under "Colour fires by" -----------------------------------------------------------------
sub('''    <label class="row" title="Fires without a pyroCb are grey by default (the All fire activity layer)"><input type="checkbox" id="paintAll"> Include fires without a pyroCb</label>''',
    '''    <label class="row" title="Fires without a pyroCb are grey by default (the All fire activity layer)"><input type="checkbox" id="paintAll"> Include fires without a pyroCb</label>
    <span class="lbl">Colour pyroCbs by</span>
    <div style="display:flex;gap:6px;align-items:center"><select class="fld" id="pcbColSel">
      <option value="none">Same as fires</option>
      <optgroup label="Injection height">
        <option value="inj">1 &middot; Injection height (m ASL)</option>
        <option value="trop">2 &middot; Tropopause height (m ASL)</option>
        <option value="rel">3 &middot; Injection &minus; tropopause (m)</option>
        <option value="above">4 &middot; Injection above fire ground (m)</option>
        <option value="frac">5 &middot; Share of troposphere above the fire reached</option>
      </optgroup></select><button class="ibtn" id="pcbColInfo" title="What is this and how was it calculated?">i</button></div>
    <div class="defbox" id="pcbColDef" hidden></div>''')
sub('''  const setMode = m => { MODE = m;''',
    '''  const pcbColSel = d.querySelector("#pcbColSel"), pcbColDef = d.querySelector("#pcbColDef");
  const syncPcbDef = () => { const v = PCB_CVARS[PCB_CMODE];
    pcbColDef.innerHTML = v ? `<b>${v.label}</b> &ndash; ${v.def}` : "pyroCb markers take the colour of their fire (year / day of year / variable above)."; };
  pcbColSel.onchange = () => { PCB_CMODE = pcbColSel.value; syncPcbDef(); recolour(); };
  d.querySelector("#pcbColInfo").onclick = () => { pcbColDef.hidden = !pcbColDef.hidden; syncPcbDef(); };
  const setMode = m => { MODE = m;''')

# ---- legend ------------------------------------------------------------------------------------------
sub('''  legendBody.innerHTML = features + colour + lc + apcKey + note + notes;''',
    '''  const pv = PCB_CVARS[PCB_CMODE];
  const pcbKey = pv ? `<span class="lbl">pyroCb colour</span>` + rampKey(pv.label + (pv.unit ? ` (${pv.unit})` : ""),
      pv.stops, (pv.unit ? Math.round(pv.lo).toLocaleString("en") : pv.lo) + " or less",
      (pv.unit ? Math.round(pv.hi).toLocaleString("en") : pv.hi) + " or more", pv.mid || null)
    + `<div class="keylab" style="justify-content:flex-start;gap:6px;margin-top:4px">`
    + `<span class="swatch" style="background:${PCB_NOVAL}"></span>no value</div>` : "";
  legendBody.innerHTML = features + colour + pcbKey + lc + apcKey + note + notes;''')

TAIL.write_text(t, encoding="utf-8", newline="\n")
print("patched")
