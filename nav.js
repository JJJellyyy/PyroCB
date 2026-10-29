/* Shared page navigation for the pyroCb site: one small pill with the sub-pages, current page highlighted.
   Include with <script src="nav.js" defer data-top="10"></script> (data-top = distance from the top in px).
   To add a page, add it to PAGES below - every page picks it up. Wide screens: a row of links; narrow screens:
   one button that opens the list. The year / fire in the URL hash is not carried over (each page keeps its own). */
(function () {
  const PAGES = [
    ["index.html", "PyroCb map", "Fires and pyroCbs 2013-2023 with growth bands, per-fire animation and variables"],
    ["compare.html", "Atlas comparison", "Published vs new atlas versions, S-NPP + MODIS v2, MODIS burned area, differences"],
    ["fire3d.html", "Fires in 3D", "The 10 largest fires as terrain blocks with their burn history"],
  ];
  const me = document.currentScript;
  const top = (me && me.dataset.top) || "10";
  let here = location.pathname.split("/").pop() || "index.html";
  if (!PAGES.some(p => p[0] === here)) here = "index.html";

  const css = `
  .pnav{position:fixed;z-index:10000;top:${top}px;left:50%;transform:translateX(-50%);display:flex;gap:2px;padding:3px;
        background:var(--pnav-bg);border:1px solid var(--pnav-line);border-radius:999px;box-shadow:0 1px 2px rgba(16,24,40,.06),0 4px 14px rgba(16,24,40,.14);
        font:500 12.5px/1 Inter,ui-sans-serif,system-ui,-apple-system,"Segoe UI",Roboto,sans-serif}
  .pnav a{color:var(--pnav-fg);text-decoration:none;padding:7px 12px;border-radius:999px;white-space:nowrap}
  .pnav a:hover{background:var(--pnav-hover)}
  .pnav a[aria-current="page"]{background:var(--pnav-fg);color:var(--pnav-bg);font-weight:600}
  .pnav button{display:none;border:none;background:none;color:var(--pnav-fg);font:inherit;padding:7px 12px;cursor:pointer}
  .pnav{--pnav-bg:#fff;--pnav-fg:#1c2430;--pnav-line:#e4e7ec;--pnav-hover:#f2f4f7}
  @media (prefers-color-scheme: dark){:root:not([data-theme="light"]) .pnav{--pnav-bg:#161b22;--pnav-fg:#e6e9ee;--pnav-line:#2b3340;--pnav-hover:#232a35}}
  :root[data-theme="dark"] .pnav{--pnav-bg:#161b22;--pnav-fg:#e6e9ee;--pnav-line:#2b3340;--pnav-hover:#232a35}
  @media (max-width: 760px){
    .pnav{flex-direction:column;border-radius:14px;left:auto;right:16px;transform:none;gap:0}
    .pnav button{display:block}
    .pnav a{display:none;border-radius:10px}
    .pnav.open a{display:block}
  }`;
  const st = document.createElement("style");
  st.textContent = css;
  document.head.appendChild(st);

  const nav = document.createElement("nav");
  nav.className = "pnav";
  nav.setAttribute("aria-label", "Pages");
  const cur = PAGES.find(p => p[0] === here);
  const btn = document.createElement("button");
  btn.type = "button";
  btn.textContent = (cur ? cur[1] : "Pages") + " ▾";
  btn.setAttribute("aria-expanded", "false");
  btn.onclick = () => { const o = nav.classList.toggle("open"); btn.setAttribute("aria-expanded", String(o)); };
  nav.appendChild(btn);
  for (const [href, label, title] of PAGES) {
    const a = document.createElement("a");
    a.href = href; a.textContent = label; a.title = title;
    if (href === here) a.setAttribute("aria-current", "page");
    nav.appendChild(a);
  }
  document.addEventListener("click", e => { if (!nav.contains(e.target)) nav.classList.remove("open"); });
  (document.body ? Promise.resolve() : new Promise(r => addEventListener("DOMContentLoaded", r))).then(() => document.body.appendChild(nav));
})();
