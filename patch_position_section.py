"""One-off patch of map_src/tail.html + head.html (2026-10-03): the "pyroCb position" section now offers only
the reported point and the NEW coordinate (pcbcube's matching, np_* properties); v2/v3/v4 are gone. The active
pyroCb pixels tickbox lives under it, each pyroCb gets its own colour + badge, the pyroCb info panel can be
minimised.   python patch_position_section.py   then   python map_tool.py join && node map_check.js
"""
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
TAIL, HEAD = HERE / "map_src" / "tail.html", HERE / "map_src" / "head.html"
t = TAIL.read_text(encoding="utf-8")
h = HEAD.read_text(encoding="utf-8")
if "np_lat" in t:
    raise SystemExit("already patched")


def sub(old, new, count=1):
    global t
    assert t.count(old) == count, (t.count(old), old[:90])
    t = t.replace(old, new)


def cut(start, end, new=""):
    """Replace t[start_marker : end_marker) (end marker kept)."""
    global t
    a = t.index(start)
    b = t.index(end, a)
    t = t[:a] + new + t[b:]


# ---------------------------------------------------------------- state
cut('let PCB_POS  = "fire4";', 'const regionOK',
    '''let PCB_POS  = "new";                  // "new" = pcbcube's matched coordinate (np_lat/np_lon, falls back to the report) | "report" = atlas lat/lon
let SHOW_ERR = false;                  // draw the +/-0.1 deg report uncertainty
let SHOW_MOVE = false;                 // arrow: as reported -> new coordinate
''')

# ---------------------------------------------------------------- info panel rows
cut('    ["Matched fire v4 (4 overpasses, time-weighted)"', '  ];\n  // one segmented control',
    '''    ["New coordinate (matched)", p.np_lat!=null
      ? `${p.np_lat.toFixed(3)}, ${p.np_lon.toFixed(3)} &middot; ${p.np_move_km} km from the report &middot; `
        + `${p.np_n} px within ${p.np_r} km${p.np_wk ? " &middot; <b>weak match</b>" : ""}`
        + (p.np_src === "pixels_grown" ? " &middot; position from the credited fire's connected pixels" : "")
      : "no detections in range &ndash; stays at the reported point"],
''')

# ---------------------------------------------------------------- position section HTML
cut('    <details class="sec extOnly" hidden><summary>pyroCb position</summary><div>', '    <details class="sec extOnly" hidden><summary>pyroCb click</summary>',
    '''    <details class="sec" open><summary>pyroCb position</summary><div>
      <div class="seg sm"><button id="posNew" class="on">New coordinate</button><button id="posRep">As reported</button></div>
      <div class="hint" style="margin:6px 0 4px">New coordinate = the dataset's matched position: VIIRS pixels of all fires within 15 km
      (widened to 25 km if empty) of the reported point in the 2 overpasses before and 2 after onset, weighted by
      exp(&minus;|&Delta;t|/12&nbsp;h); thin matches use the credited fire's connected pixels instead. Unmatched events stay at the report.</div>
      <label class="row"><input type="checkbox" id="showMove"> Arrow: as reported &rarr; new coordinate</label>
      <label class="row"><input type="checkbox" id="showErr"> Report uncertainty box</label>
      <div class="hint">93 % of atlas coordinates are rounded to 0.1&deg; (~11 km N-S); the box is drawn at that grid
      (tighter for the few finer ones).</div>
      <label class="row" title="The dataset's label = 1 pixels: every VIIRS pixel within 10 km of the new coordinate, in the overpass just before onset and the one just after"><input type="checkbox" data-kind="apc"> Active pyroCb pixels <span class="opt">10 km, before + after onset</span></label>
      <div class="hint">Each pyroCb has its own colour: ring, badge and pixels match. Select a pyroCb to highlight its own pixels.</div>
    </div></details>
''')
# the tickbox no longer sits under Layers
cut('      <label class="row" title="The dataset\'s label = 1 pixels: every VIIRS pixel within 20 km', '      <label class="row"><input type="checkbox" data-kind="lc">', '')

# ---------------------------------------------------------------- handlers
cut('  // ---- pyroCb position (advanced) ---', '  const dp=d.querySelector("#dPanel")',
    '''  // ---- pyroCb position --------------------------------------------------------------------
  const pn = d.querySelector("#posNew"), pr = d.querySelector("#posRep");
  const setPos = (mode, on) => {
    PCB_POS = mode; [pn, pr].forEach(b => b.classList.toggle("on", b === on));
    repositionPcb(); drawLegend();
  };
  pn.onclick = () => setPos("new", pn);
  pr.onclick = () => setPos("report", pr);
  d.querySelector("#showErr").onchange = e => { SHOW_ERR = e.target.checked; repositionPcb(); };
  d.querySelector("#showMove").onchange = e => { SHOW_MOVE = e.target.checked; repositionPcb(); drawLegend(); };
''')

# ---------------------------------------------------------------- legend notes
cut('  const pos = (c, txt) => `<div class="hint" style="color:${c};margin-top:6px">${txt}</div>`;', '  const apcKey = ',
    '''  const notes = PCB_POS === "new"
    ? `<div class="hint" style="color:#b45309;margin-top:6px">pyroCbs plotted on the NEW coordinate (dataset matching).`
      + `${SHOW_MOVE ? " Line + triangle = as reported &rarr; new." : ""} Unmatched events stay at the reported point.</div>` : "";
''')
cut('  const apcKey = showKind.apc', '  legendBody.innerHTML',
    '''  const apcKey = showKind.apc ? `<span class="lbl">Active pyroCb pixels</span>
    <div class="row"><i class="sq" style="background:#8a8a8a"></i>Solid = overpass before onset</div>
    <div class="row"><i class="sq" style="background:#d4d4d4;border:1px solid #8a8a8a"></i>Pale = overpass after onset</div>
    <div class="row"><i style="border:1.5px dashed #555;background:none"></i>10 km around the new coordinate</div>
    <div class="hint">One colour per pyroCb: its ring, its <b>#n</b> badge (same n as the numbered marker) and its pixels.
    Dark outline = a pixel of another fire than the credited one. All are label&nbsp;=&nbsp;1 in the dataset.</div>` : "";
''')

# ---------------------------------------------------------------- position function + the nb2/nb3/nb4 machinery
cut('function pcbDisplayLatLng(p, coords){', '// fsim.pcb_ll is baked',
    '''function pcbDisplayLatLng(p, coords){
  if (PCB_POS === "new" && p.np_lat != null && p.np_lon != null) return [p.np_lat, p.np_lon];
  return [coords[1], coords[0]];
}

''')
cut('const v2Key = (base)', '// Where a pyroCb sits on the map right now', '')
cut('  rebuildNb2();\n  rebuildNb3();\n  rebuildNb4();\n}', '// Visibility filter: country toggles',
    '''  rebuildMove();
}

// ---- arrow: as reported -> new coordinate -----------------------------------------------------
function bearingDeg(lat1, lon1, lat2, lon2){
  const rad = d => d * Math.PI / 180;
  const dLon = rad(lon2 - lon1), la1 = rad(lat1), la2 = rad(lat2);
  const y = Math.sin(dLon) * Math.cos(la2);
  const x = Math.cos(la1) * Math.sin(la2) - Math.sin(la1) * Math.cos(la2) * Math.cos(dLon);
  return (Math.atan2(y, x) * 180 / Math.PI + 360) % 360;
}
let moveLayer = null;
function rebuildMove(){
  if (moveLayer){ map.removeLayer(moveLayer); moveLayer = null; }
  if (PCB_POS !== "new" || !SHOW_MOVE) return;
  moveLayer = L.layerGroup().addTo(map);
  for (const y of Object.keys(cache)){
    if (!activeYears.has(+y) || !cache[y].pcb) continue;
    cache[y].pcb.eachLayer(l => {
      const p = l.feature.properties, g = l.feature.geometry.coordinates;
      if (!pcbShown(l.feature) || p.np_lat == null) return;
      const from = [g[1], g[0]], to = [p.np_lat, p.np_lon];
      if (from[0] === to[0] && from[1] === to[1]) return;
      L.polyline([from, to], {pane:"nb2", renderer:R.nb2, color:"#b45309", weight:2, opacity:0.85,
        interactive:false}).addTo(moveLayer);
      const brg = bearingDeg(from[0], from[1], to[0], to[1]);
      L.marker(to, {pane:"nb2", interactive:false, icon:L.divIcon({className:"", iconSize:[14,14], iconAnchor:[7,7],
        html:`<div style="width:0;height:0;border-left:6px solid transparent;border-right:6px solid transparent;`
          +`border-bottom:11px solid #b45309;transform:rotate(${brg}deg);transform-origin:50% 70%"></div>`})}).addTo(moveLayer);
    });
  }
}

''')
# the old block between the old repositionPcb tail and filterKey is now fully replaced (cut ends at the filterKey comment)
sub('''  if (errLayer || SHOW_NB2_CIRCLE || SHOW_NB2_FIRE || SHOW_NB3_CIRCLE || SHOW_NB3_FIRE
      || SHOW_NB4_CIRCLE || SHOW_NB4_FIRE) repositionPcb();''', '''  if (errLayer || SHOW_MOVE) repositionPcb();''')

# ---------------------------------------------------------------- active pyroCb pixels: one colour per pyroCb
cut('const APC_COL = ["#f59e0b", "#dc2626"];', 'const R_APC = () => R.apc;',
    '''// One colour per pyroCb (Okabe-Ito + 3 extras, assigned in time order within the year so neighbours differ).
const APC_HUES = ["#E69F00", "#56B4E9", "#009E73", "#D55E00", "#CC79A7", "#0072B2", "#7b3294", "#a6761d", "#1b9e77", "#e7298a"];
const apcByPid = {};                                   // pcb_id -> {hue, layers:[...], fillA:[...]}
function apcLayer(y, apc){
  const grp = L.layerGroup();
  let built = false;
  grp.on("add", () => {
    if (built) return; built = true;
    const pids = Object.keys(apc || {}).sort();      // pcb_id starts with the date: time order
    pids.forEach((pid, k) => {
      const e = apc[pid], hue = APC_HUES[k % APC_HUES.length], R = R_APC();
      const [lo0, la0] = e.c, own = {hue, layers: [], solid: [], r: e.r};
      apcByPid[pid] = own;
      const nSide = [0, 0], nFire = {};
      for (const q of e.p){
        const lo = lo0 + q[0]*1e-5, la = la0 + q[1]*1e-5, hx = q[2]*1e-5, hy = q[3]*1e-5;
        nSide[q[4]]++; nFire[q[5]] = (nFire[q[5]]||0) + 1;
        const credited = e.fires[q[5]] === e.mg, before = q[4] === 0;
        const r = L.rectangle([[la-hy, lo-hx], [la+hy, lo+hx]], {pane:"apc", renderer:R, interactive:false,
          stroke:true, color: credited ? hue : "#111", weight: credited ? 0.5 : 0.9,
          fillColor:hue, fillOpacity: before ? 0.9 : 0.4});
        r._apcBase = {fo: before ? 0.9 : 0.4, so: 1};
        r.addTo(grp); own.layers.push(r);
      }
      // 10 km ring in the pyroCb's colour, with a #n badge matching the numbered pyroCb marker
      const ring = L.circle([la0, lo0], {pane:"apc", renderer:R, radius:e.r*1000, interactive:false,
        color:hue, weight:2.2, dashArray:"7 5", fill:false});
      ring.addTo(grp); own.ring = ring;
      const ft = pcbFeatureById(pid), pr = ft ? ft.properties : {};
      const badge = L.marker([la0 + e.r/111.32, lo0], {pane:"apc", interactive:true, icon:L.divIcon({className:"",
        iconSize:[1,1], iconAnchor:[0,10],
        html:`<div class="apcBadge" style="background:${hue}">#${pr.pcb_index != null ? pr.pcb_index : "?"}`
          + `<span>${esc(String(pid).slice(9, 13))}-${esc(String(pid).slice(14, 16))}${esc(String(pid).slice(16, 18))}Z</span></div>`})});
      const fires = e.fires.map((m, i) => `${esc(m)}${m === e.mg ? " (credited)" : ""}: ${nFire[i]||0}`).join("<br>");
      badge.bindPopup(`<b style="color:${hue}">&#9632;</b> <b>Active pyroCb &middot; ${esc(pid)}</b><table>${rows([
          ["New coordinate", `${la0.toFixed(4)}, ${lo0.toFixed(4)}`],
          ["Position from", e.src === "pixels_grown" ? "credited fire's connected pixels (thin match)" : "time-weighted pixel centroid"],
          ["Weak match", e.wk ? "yes" : "no"],
          ["Overpass(es) at/before onset", (e.s[0] || "-") + ` &middot; ${nSide[0]} px`],
          ["Overpass(es) after onset", (e.s[1] || "-") + ` &middot; ${nSide[1]} px`],
          ["Pixels per fire", fires]])}</table>
          <div class="hint">All pixels within ${e.r} km of the new coordinate in these overpasses are the dataset's
          label&nbsp;=&nbsp;1 pixels, whichever fire they belong to.</div>`);
      badge.on("click", () => apcHighlight(pid));
      badge.addTo(grp);
      L.circleMarker([la0, lo0], {pane:"apc", renderer:R, radius:4, color:hue, weight:2.5, fillColor:"#fff",
        fillOpacity:1, interactive:false}).addTo(grp);
    });
  });
  return grp;
}
// Highlight one pyroCb's pixels (others fade); null clears. Called when a pyroCb is opened or its badge clicked.
let apcFocus = null;
function apcHighlight(pid){
  apcFocus = pid;
  for (const [k, o] of Object.entries(apcByPid)){
    const dim = pid != null && k !== pid;
    o.layers.forEach(r => r.setStyle({fillOpacity: dim ? 0.06 : r._apcBase.fo, opacity: dim ? 0.12 : 1,
                                      weight: !dim && pid != null ? 1.1 : (r.options.color === "#111" ? 0.9 : 0.5)}));
    o.ring && o.ring.setStyle({opacity: dim ? 0.2 : 1, weight: !dim && pid != null ? 3.4 : 2.2});
  }
}
''')

# open pyroCb -> highlight; close -> clear
sub('function closeInfo(){ stopAnim(); stopGoes(); infoEl.style.display = "none"; clearHighlight(); clearPops(); syncFlow();',
    'function closeInfo(){ stopAnim(); stopGoes(); infoEl.style.display = "none"; clearHighlight(); clearPops(); syncFlow(); if (apcFocus) apcHighlight(null);')
sub('  infoEl.style.display = "block";\n  infoEl.querySelector(".x").onclick = closeInfo;',
    '''  infoEl.style.display = "block";
  infoEl.classList.remove("min");
  infoEl.querySelector(".x").onclick = closeInfo;
  const mn = infoEl.querySelector(".mn");
  mn.onclick = () => { const on = infoEl.classList.toggle("min"); mn.innerHTML = on ? "&#9650;" : "&#9660;";
                       mn.title = on ? "Expand" : "Minimise"; };
  if (!FIRE_ONLY && p.pcb_id && showKind.apc) apcHighlight(p.pcb_id);''')
# minimise button in the header
sub('''+`</div><div class="s">${sub}</div><span class="x" title="Close">&times;</span></div>`''',
    '''+`</div><div class="s">${sub}</div><span class="mn" title="Minimise">&#9660;</span><span class="x" title="Close">&times;</span></div>`''')
# keep highlight in step with the tickbox
sub('''    if (i.dataset.kind==="lc" || i.dataset.kind==="apc") drawLegend();});''',
    '''    if (i.dataset.kind==="lc" || i.dataset.kind==="apc") drawLegend();
    if (i.dataset.kind==="apc" && !i.checked) apcFocus = null;});''')
TAIL.write_text(t, encoding="utf-8", newline="\n")

# ---------------------------------------------------------------- CSS
css = '''
  /* minimise the pyroCb info panel: header only */
  #info .hd .mn{position:absolute;right:36px;top:8px;width:24px;height:24px;display:flex;align-items:center;
        justify-content:center;border-radius:6px;cursor:pointer;font-size:11px;line-height:1;opacity:.8}
  #info .hd .mn:hover{background:rgba(255,255,255,.14);opacity:1}
  #info.min{overflow:hidden;max-height:none}
  #info.min .bd{display:none}
  #info.min .hd{border-radius:var(--radius)}
  #info .hd{padding-right:68px}
  /* active-pyroCb badge: same colour as the ring and pixels */
  .apcBadge{display:inline-block;transform:translate(-50%,0);color:#fff;font:700 11px/1 var(--font);padding:3px 6px;
        border-radius:9px;border:1.5px solid #fff;box-shadow:0 1px 4px rgba(0,0,0,.45);white-space:nowrap;cursor:pointer;
        text-shadow:0 0 2px rgba(0,0,0,.5)}
  .apcBadge span{font-weight:500;opacity:.9;margin-left:5px}
'''
assert "</style>" in h
i = h.index("</style>")
h = h[:i] + css + h[i:]
HEAD.write_text(h, encoding="utf-8", newline="\n")
print("patched")
