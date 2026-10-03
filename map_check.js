// syntax-check the map's inline script: node map_check.js   (DATA is stubbed out)
const fs = require("fs");
const tail = fs.readFileSync("map_src/tail.html", "utf8");
const js = "const DATA={}" + tail.slice(0, tail.lastIndexOf("</script>"));
fs.writeFileSync(require("os").tmpdir() + "/map_check_script.js", js);
try { new (require("vm").Script)(js, {filename: "map.js"}); console.log("syntax OK"); }
catch (e) { console.log(e.stack.split("\n").slice(0, 6).join("\n")); process.exit(1); }
