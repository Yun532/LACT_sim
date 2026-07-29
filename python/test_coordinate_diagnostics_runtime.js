const fs = require("fs");
const vm = require("vm");

const html = fs.readFileSync("docs/assets/lact-coordinate-system-3d.html", "utf8");
if (!html.includes('<details id="infoPanel" class="panel">') ||
    html.includes('<details id="infoPanel" class="panel" open>')) {
  throw new Error("left information panel must exist and be collapsed by default");
}
if (!html.includes('<details id="controlPanel" class="panel">') ||
    html.includes('<details id="controlPanel" class="panel" open>')) {
  throw new Error("view and parameter controls must exist and be collapsed by default");
}
for (const marker of [
  '<details class="textfold"><summary>页面说明</summary>',
  '<details class="textfold"><summary>图像说明</summary>',
  '<details class="textfold"><summary>图表说明</summary>',
  '固定地图：北上东右',
  '北 N',
  '南 S',
  '东 E',
  '西 W',
  '屏幕方向（随视角）',
  'title="左拖旋转；Shift 或右拖平移；滚轮以光标为中心缩放">拖动 / 缩放',
  'cctx.fillStyle="#0b1925"',
  'const PARALLEL_COLORS=["#ff6d7a","#64c7ff","#ffd45e","#70e39c"]',
  'caseDensityColor(plotColor',
  'id="cameraCoords"',
  'id="corsikaPanelMode"',
  'function corePlotColor(t)',
  'q.fillStyle="#000000"',
  'id="deformMetric"',
  'function pylastReaderCoordinates(u,v)',
  'function pylastPlotCoordinates(pix_x,pix_y)',
  'function drawPylastCameraPanel(',
  'function drawUnifiedPixelCamera(',
  'function drawPylastCorePlot(',
  'function corsikaPrimaryViews()',
  'function corsikaPreviewViews()',
  'id="returnShower"',
  'camera_geometry_lact_uv',
  'LACT +u = camera x_m = local +x',
  'LACT +v = camera y_m = local +y',
  '相机面：ROOT x_m≡LACT u；ROOT y_m≡LACT v',
  '+u (=x_m)',
  '+v (=y_m)',
  'camera_x_m=hit.u_m',
  'camera_y_m=hit.v_m',
  '不是 CORSIKA 阵列 NWU 的 x/y',
  '矩阵是 <code>[[0,-1],[-1,0]]</code>',
  '正方位基底朝 West',
  'pyLAST +pix_x = -v',
  'pyLAST +pix_y = -u',
  'pyLAST 画布：向右 = -u；向上 = -v',
  'NWU：x=N，y=W，E=−y',
  '显示：+x=N，+y=E',
  'function corsikaAngleDisplay()',
  'function beamDirectionFromAngles(thetaDeg,phiDeg)',
]) {
  if (!html.includes(marker)) throw new Error(`generated page is missing UI marker: ${marker}`);
}
const scriptStart = html.lastIndexOf("<script>") + "<script>".length;
const scriptEnd = html.indexOf("</script>", scriptStart);
if (scriptStart < "<script>".length || scriptEnd < 0) throw new Error("runtime script missing");
const source = html.slice(scriptStart, scriptEnd);

function canvasContext() {
  return new Proxy({}, {
    get(target, key) {
      if (key === "measureText") return text => ({width: String(text).length * 7});
      if (!(key in target)) target[key] = () => {};
      return target[key];
    },
    set(target, key, value) { target[key] = value; return true; },
  });
}

function element(id, dataset = {}) {
  const context = canvasContext();
  const initialValues = {cameraScale: "spot", cameraCoords: "lact", corsikaPanelMode: "camera", deformMetric: "normal", globalAz: "0", globalEl: "70", photonTel: "19"};
  return {
    id, dataset, style: {}, classList: {add() {}, remove() {}, toggle() {}},
    value: initialValues[id] || "0", hidden: false,
    width: 900, height: 700, innerHTML: "", textContent: "", max: 0,
    listeners: {}, listenerLists: {}, addEventListener(type, handler) {
      (this.listenerLists[type] ||= []).push(handler);
      this.listeners[type] = event => this.listenerLists[type].forEach(fn => fn(event));
    }, setPointerCapture() {},
    getBoundingClientRect() { return {left: 0, top: 0, width: id === "scene" ? 900 : 400, height: id === "scene" ? 700 : 336}; },
    getContext() { return context; },
  };
}

function run(search) {
  const ids = ["scene", "camera", "info", "cameraPanel", "cameraPanelTitle", "cameraScale", "cameraCoords", "corsikaPanelMode", "cameraHead",
    "cameraCaption", "deformPanel", "deformChart", "deformTitle", "deformStat", "colorMetric", "colorMin", "colorMax",
    "globalToolbar", "parallelToolbar", "elevationToolbar", "corsikaToolbar", "globalAz",
    "globalAzValue", "globalEl", "globalElValue", "angleValue", "bunch", "photonTel",
    "bunchValue", "elevation", "deformMetric", "arrowScale", "corsikaCase", "groundView", "showerView",
    "zoomOut", "zoomIn", "subtitle", "legend", "returnShower"];
  const elements = new Map(ids.map(id => [id, element(id)]));
  const pageButtons = ["global", "parallel", "elevation", "corsika"].map(p => element("page-" + p, {page: p}));
  const viewButtons = ["horizon", "iso", "top", "north", "east"].map(view => element("view-" + view, {view}));
  const document = {
    documentElement: {},
    querySelector(selector) {
      if (selector.startsWith("#")) return elements.get(selector.slice(1));
      if (selector === '[data-page="corsika"]') return pageButtons.find(button => button.dataset.page === "corsika");
      return null;
    },
    querySelectorAll(selector) {
      if (selector === "[data-page]") return pageButtons;
      if (selector === "[data-view]") return viewButtons;
      if (selector === ".reset") return [];
      return [];
    },
  };
  const sandbox = {
    console, document, location: {search}, URLSearchParams,
    devicePixelRatio: 1, innerWidth: 1280, innerHeight: 800,
    getComputedStyle() { return {getPropertyValue() { return "#ffffff"; }}; },
    addEventListener() {}, Math, Map, Infinity,
  };
  vm.runInNewContext(source, sandbox, {filename: "lact-coordinate-system-3d.html"});
  elements.sandbox = sandbox;
  return elements;
}

const checks = [
  ["?page=global", "当前程序角度"],
  ["?page=parallel", "完整输出绘图"],
  ["?page=elevation", "当前角度"],
  ["?page=corsika&case=event1909", "zem 与发射点"],
];
for (const [query, expected] of checks) {
  const elements = run(query);
  if (!elements.get("info").innerHTML.includes(expected)) {
    throw new Error(`${query}: info panel is missing ${expected}`);
  }
  console.log("runtime OK", query);
}

const globalElements = run("?page=global");
const globalAxisNames = globalElements.sandbox.globalScene().lines.map(line => line.name);
for (const expected of ["pyLAST +pix_x = -v", "pyLAST +pix_y = -u"]) {
  if (!globalAxisNames.includes(expected)) {
    throw new Error(`global definition scene is missing ${expected}`);
  }
}
const globalInfo = globalElements.get("info").innerHTML;
for (const expected of ["pix_x=-v", "横轴=pix_y=-u", "root/LactEventSource.cpp:245-265"]) {
  if (!globalInfo.includes(expected)) {
    throw new Error(`global pyLAST definition is missing ${expected}`);
  }
}
console.log("runtime OK global page contains pyLAST field and canvas definitions");

const elevationWhiteboardElements = run("?page=elevation");
const elevationCase = elevationWhiteboardElements.sandbox.currentParallelCase();
const elevationHistogram = elevationWhiteboardElements.sandbox.spotHistogram(elevationCase, 140);
if (!elevationWhiteboardElements.get("cameraHead").textContent.includes("完整 output u/v") ||
    elevationHistogram.n <= 0 ||
    elevationHistogram.n !== elevationCase.full_output_uv_m.length ||
    elevationHistogram.n !== elevationCase.validation.full_output_uv_rows ||
    elevationHistogram.n !== elevationCase.summary.hit_output_plane) {
  throw new Error("elevation page does not reuse the complete program-style whiteboard output");
}
console.log("runtime OK elevation uses complete program-style whiteboard output");

const normalScale = elevationWhiteboardElements.sandbox.deformationColorScale();
const normalScene = elevationWhiteboardElements.sandbox.buildScene();
if (elevationWhiteboardElements.get("deformMetric").value !== "normal" ||
    normalScale.min !== 0 || !Number.isFinite(normalScale.max) || normalScale.max <= 0 ||
    elevationWhiteboardElements.get("colorMin").textContent !== "0.000 mdeg" ||
    !elevationWhiteboardElements.get("colorMax").textContent.includes(normalScale.max.toFixed(3)) ||
    !elevationWhiteboardElements.get("colorMetric").textContent.includes("Δnormal") ||
    !normalScene.lines.some(line => line.name.includes("法向偏转"))) {
  throw new Error("elevation page does not default to a shared normal-angle scale");
}
elevationWhiteboardElements.get("deformMetric").listeners.change({target: {value: "position"}});
const positionScale = elevationWhiteboardElements.sandbox.deformationColorScale();
const positionScene = elevationWhiteboardElements.sandbox.buildScene();
if (positionScale.min !== 0 || Math.abs(positionScale.max - 2.2563430904) > 1e-9 ||
    elevationWhiteboardElements.get("colorMin").textContent !== "0.000 mm" ||
    !elevationWhiteboardElements.get("colorMetric").textContent.includes("Δcenter") ||
    !positionScene.lines.some(line => line.name.includes("中心位移向量"))) {
  throw new Error("elevation page cannot switch to the shared center-displacement scale");
}
console.log("runtime OK selectable deformation metrics", {normalScale, positionScale});

const corsikaFullCameraElements = run("?page=corsika&case=event1909");
if (corsikaFullCameraElements.get("cameraScale").value !== "full" ||
    !corsikaFullCameraElements.get("cameraHead").textContent.includes("ROOT PE 前四") ||
    vm.runInNewContext("globalScale", corsikaFullCameraElements.sandbox) !== "shower") {
  throw new Error("CORSIKA page must default to the full shower and four ROOT camera previews");
}
const previewIds = vm.runInNewContext("corsikaPreviewViews().map(function(v){return v.telescope_id})", corsikaFullCameraElements.sandbox);
const previewColors = vm.runInNewContext("corsikaPreviewViews().map(function(v){return corsikaRayColor(v.telescope_id)})", corsikaFullCameraElements.sandbox);
if (previewIds.join(",") !== "19,16,21,20" ||
    previewColors.join(",") !== "#ff6d7a,#64c7ff,#ffd45e,#70e39c" ||
    !corsikaFullCameraElements.get("cameraCaption").textContent.includes("triggered 标志均为 false")) {
  throw new Error(`CORSIKA preview cameras/colors do not reproduce the four original highest-PE views: ${previewIds}; ${previewColors}`);
}
const groundPointBeforeZoom = corsikaFullCameraElements.sandbox.project([0, 0, 0]);
corsikaFullCameraElements.sandbox.zoomAt(120, groundPointBeforeZoom[0], groundPointBeforeZoom[1]);
const groundPointAfterZoom = corsikaFullCameraElements.sandbox.project([0, 0, 0]);
if (Math.hypot(groundPointAfterZoom[0] - groundPointBeforeZoom[0],
               groundPointAfterZoom[1] - groundPointBeforeZoom[1]) > 1e-7) {
  throw new Error("cursor-centered zoom failed to keep the selected ground point fixed");
}
console.log("runtime OK CORSIKA defaults to full shower, four ROOT camera previews and cursor-centered zoom");

const corsikaCoreElements = run("?page=corsika&case=event1909");
corsikaCoreElements.get("corsikaPanelMode").value = "cores";
corsikaCoreElements.get("corsikaPanelMode").listeners.change({target: {value: "cores"}});
const coreLayout = vm.runInNewContext("corsikaCorePlotLayout", corsikaCoreElements.sandbox);
const tel19Core = vm.runInNewContext("D.event1909.array.find(function(t){return t.telescope_id===19})", corsikaCoreElements.sandbox);
const tel19Plot = coreLayout.points.find(point => point.telescope_id === 19);
if (!coreLayout || coreLayout.points.length !== 32 || !tel19Plot ||
    corsikaCoreElements.get("cameraCoords").style.display !== "none" ||
    corsikaCoreElements.get("cameraScale").style.display !== "none" ||
    !corsikaCoreElements.get("cameraHead").textContent.includes("pyLAST 阵列 / true core") ||
    !corsikaCoreElements.get("cameraCaption").textContent.includes("East=−West") ||
    !corsikaCoreElements.get("cameraCaption").textContent.includes("高对比蓝—青—绿—黄—红") ||
    !corsikaCoreElements.get("cameraCaption").textContent.includes("include_non_triggered=True") ||
    !corsikaCoreElements.get("info").innerHTML.includes("plot_x=East=−array_y_west")) {
  throw new Error("pyLAST core plot does not expose all 32 ROOT telescope positions, high-contrast styling and the audited NWU-to-East/North rule");
}
corsikaCoreElements.get("camera").listeners.click({
  clientX: tel19Plot.x / coreLayout.devicePixelRatio,
  clientY: tel19Plot.y / coreLayout.devicePixelRatio,
});
const coreClickedTel = vm.runInNewContext("state.photonTel", corsikaCoreElements.sandbox);
const coreClickedCenter = vm.runInNewContext("center", corsikaCoreElements.sandbox);
if (coreClickedTel !== 19 || vm.runInNewContext("globalScale", corsikaCoreElements.sandbox) !== "telescope" ||
    Math.hypot(coreClickedCenter[0] - tel19Core.position_nwu_m[0], coreClickedCenter[1] - tel19Core.position_nwu_m[1]) > 1e-6) {
  throw new Error("clicking the pyLAST core plot did not focus the matching ROOT telescope model");
}
console.log("runtime OK high-contrast pyLAST core plot maps NWU to East/North and selects the matching 3D telescope");

const rootCameraCheck = vm.runInNewContext(`(() => {
  const root = D.event1909.camera_geometry_lact_uv;
  const configured = new Map(D.camera_geometry.map(p => [p.id, p]));
  let maxError = 0;
  for (const p of root) {
    const q = configured.get(p.id);
    maxError = Math.max(maxError, Math.abs(p.u-q.u), Math.abs(p.v-q.v), Math.abs(p.size-q.size));
  }
  return {count:root.length,maxError};
})()`, corsikaFullCameraElements.sandbox);
if (rootCameraCheck.count !== 1616 || rootCameraCheck.maxError > 1e-9 ||
    !corsikaFullCameraElements.get("info").innerHTML.includes("camera_pixels") ||
    !corsikaFullCameraElements.get("info").innerHTML.includes("x_m≡u") ||
    !corsikaFullCameraElements.get("info").innerHTML.includes("camera_x_m=hit.u_m") ||
    !corsikaFullCameraElements.get("info").innerHTML.includes("不是 CORSIKA 阵列 NWU 的 x/y") ||
    !corsikaFullCameraElements.get("info").innerHTML.includes("v=-0.571319 m") ||
    !corsikaFullCameraElements.get("info").innerHTML.includes("pix_y=0.571243 m")) {
  throw new Error(`CORSIKA LACT u/v does not come directly from the ROOT camera geometry: ${JSON.stringify(rootCameraCheck)}`);
}
console.log("runtime OK CORSIKA LACT u/v uses all 1616 ROOT camera_pixels coordinates");

const corsikaArrayElements = run("?page=corsika&case=event1909");
const corsikaArrayScene = corsikaArrayElements.sandbox.buildScene();
const allTelescopeFacets = corsikaArrayScene.polys.filter(polygon => /^tel\d+ 镜片 /.test(polygon.name));
const telescopeOptions = corsikaArrayElements.get("photonTel").innerHTML.match(/<option /g) || [];
if (allTelescopeFacets.length !== 32 * 54 || telescopeOptions.length !== 32 ||
    !corsikaArrayElements.get("info").innerHTML.includes("source.coordinate_frame=corsika_nwu_relative") ||
    !corsikaArrayElements.get("info").innerHTML.includes("32 台望远镜都可选择")) {
  throw new Error(`CORSIKA full-array geometry/config audit is incomplete: facets=${allTelescopeFacets.length}; options=${telescopeOptions.length}`);
}
corsikaArrayElements.get("photonTel").listeners.change({target: {value: "0"}});
const focusedScene = corsikaArrayElements.sandbox.buildScene();
const focusedFacets = focusedScene.polys.filter(polygon => polygon.name.startsWith("tel0 镜片 "));
const otherFocusedFacets = focusedScene.polys.filter(polygon => /^tel(?!0 )\d+ 镜片 /.test(polygon.name));
const focusedPhotonPaths = focusedScene.lines.filter(line => line.name.includes("raw zem ray"));
const focusedAxisNames = focusedScene.lines.map(line => line.name);
const focusedPreviewIds = vm.runInNewContext("corsikaPreviewViews().map(function(v){return v.telescope_id})", corsikaArrayElements.sandbox);
const focusedPreviewSignals = vm.runInNewContext("corsikaPreviewViews().map(function(v){return [v.telescope_id,v.cherenkov_pe_sum,v.camera_signal.length]})", corsikaArrayElements.sandbox);
const focusedCenter = vm.runInNewContext("center", corsikaArrayElements.sandbox);
const tel0 = vm.runInNewContext("D.event1909.array.find(function(t){return t.telescope_id===0})", corsikaArrayElements.sandbox);
if (focusedFacets.length !== 54 || otherFocusedFacets.length !== 31 * 54 || focusedPhotonPaths.length !== 32 * 120 ||
    focusedPreviewIds.join(",") !== "19,16,21,0" ||
    focusedPreviewSignals[3][1] !== 40 || focusedPreviewSignals[3][2] !== 7 ||
    Math.hypot(focusedCenter[0] - tel0.position_nwu_m[0], focusedCenter[1] - tel0.position_nwu_m[1]) > 1e-6 ||
    !focusedAxisNames.includes("LACT +u = camera x_m = local +x") ||
    !focusedAxisNames.includes("LACT +v = camera y_m = local +y") ||
    !focusedScene.labels.some(label => label.t.includes("ROOT x_m≡LACT u")) ||
    !focusedScene.labels.some(label => label.t === "+u (=x_m)") ||
    !focusedScene.labels.some(label => label.t === "+v (=y_m)") ||
    !corsikaArrayElements.get("cameraHead").textContent.includes("PE 前三 + 所选 tel0") ||
    !corsikaArrayElements.get("cameraCaption").textContent.includes("右下角立即替换为该台真实 image_cherenkov_pe") ||
    !corsikaArrayElements.get("info").innerHTML.includes("当前所选")) {
  throw new Error(`CORSIKA telescope selection filtered the complete scene or missed the tel0 rotation center: center=${focusedCenter}`);
}
console.log("runtime OK telescope focus preserves all 32 structures and 3840 real sampled photon paths");

const corsikaHierarchyElements = run("?page=corsika&case=event1909");
corsikaHierarchyElements.get("camera").listeners.click({clientX: 300, clientY: 250});
const clickedScale = vm.runInNewContext("globalScale", corsikaHierarchyElements.sandbox);
const clickedTel = vm.runInNewContext("state.photonTel", corsikaHierarchyElements.sandbox);
const clickedScene = corsikaHierarchyElements.sandbox.buildScene();
const clickedCenter = vm.runInNewContext("center", corsikaHierarchyElements.sandbox);
const tel20 = vm.runInNewContext("D.event1909.array.find(function(t){return t.telescope_id===20})", corsikaHierarchyElements.sandbox);
if (clickedScale !== "telescope" || clickedTel !== 20 ||
    clickedScene.polys.filter(polygon => /^tel\d+ 镜片 /.test(polygon.name)).length !== 32 * 54 ||
    clickedScene.lines.filter(line => line.name.includes("raw zem ray")).length !== 32 * 120 ||
    Math.hypot(clickedCenter[0] - tel20.position_nwu_m[0], clickedCenter[1] - tel20.position_nwu_m[1]) > 1e-6 ||
    corsikaHierarchyElements.get("returnShower").style.display !== "inline-block") {
  throw new Error(`CORSIKA preview click did not enter the selected telescope view: scale=${clickedScale}, tel=${clickedTel}`);
}
corsikaHierarchyElements.get("returnShower").listeners.click({});
const returnedScene = corsikaHierarchyElements.sandbox.buildScene();
if (vm.runInNewContext("globalScale", corsikaHierarchyElements.sandbox) !== "shower" ||
    returnedScene.polys.filter(polygon => /^tel\d+ 镜片 /.test(polygon.name)).length !== 32 * 54 ||
    vm.runInNewContext("viewAzDeg", corsikaHierarchyElements.sandbox) !== 0 ||
    vm.runInNewContext("viewElDeg", corsikaHierarchyElements.sandbox) !== -22 ||
    corsikaHierarchyElements.get("returnShower").style.display !== "none" ||
    !corsikaHierarchyElements.get("cameraHead").textContent.includes("ROOT PE 前四")) {
  throw new Error("CORSIKA return control did not restore the initial full-shower view");
}
console.log("runtime OK CORSIKA camera click only changes fit/rotation center and explicit return restores full-shower scale");

const parallelElements = run("?page=parallel");
const fullDensityColors = ["#ff6d7a", "#64c7ff", "#ffd45e", "#70e39c"]
  .map(color => parallelElements.sandbox.caseDensityColor(color, 1));
if (fullDensityColors.join("|") !== "rgb(255,109,122)|rgb(100,199,255)|rgb(255,212,94)|rgb(112,227,156)") {
  throw new Error(`parallel density colors do not match ray colors: ${fullDensityColors}`);
}
const parallelCase = parallelElements.sandbox.currentParallelCase();
const rawCentroid = parallelCase.camera_summary.output_uv_centroid_m;
const repeatedHistogram = parallelElements.sandbox.spotHistogram(parallelCase, 140);
if (Math.abs(repeatedHistogram.centroid[0] - rawCentroid[0]) > 1e-12 ||
    Math.abs(repeatedHistogram.centroid[1] - rawCentroid[1]) > 1e-12 ||
    repeatedHistogram.n !== parallelCase.full_output_uv_m.length ||
    parallelElements.get("cameraCoords").style.display !== "inline-block") {
  throw new Error("parallel whiteboard must preserve raw u/v and expose the audited pyLAST redraw selector");
}
parallelElements.get("camera").listeners.click({clientX: 300, clientY: 250});
if (!parallelElements.get("info").innerHTML.includes("offset=4°")) {
  throw new Error("parallel camera click did not select the bottom-right sky case");
}
if (!parallelElements.get("cameraCaption").textContent.includes("真实夹角都是 4°") ||
    !parallelElements.get("cameraCaption").textContent.includes("不是简单的方位角 ±4°")) {
  throw new Error("parallel whiteboard plots do not preserve absolute u/v ticks while auto-framing");
}
console.log("runtime OK parallel whiteboard selection and auto-framed absolute axes");

const offsetCases = vm.runInNewContext("D.parallel.four_direction_cases", parallelElements.sandbox);
for (const item of offsetCases) {
  const direction = item.beam.direction_local;
  const localOffset = Math.acos(Math.max(-1, Math.min(1, -direction[2]))) * 180 / Math.PI;
  const pointingAz = item.pointing.az_deg * Math.PI / 180;
  const pointingEl = item.pointing.el_deg * Math.PI / 180;
  const sourceAz = item.source_sky.az_deg * Math.PI / 180;
  const sourceEl = item.source_sky.el_deg * Math.PI / 180;
  const skyDot = Math.sin(pointingEl) * Math.sin(sourceEl) +
    Math.cos(pointingEl) * Math.cos(sourceEl) * Math.cos(sourceAz - pointingAz);
  const skyOffset = Math.acos(Math.max(-1, Math.min(1, skyDot))) * 180 / Math.PI;
  if (Math.abs(localOffset - 4) > 1e-6 || Math.abs(skyOffset - 4) > 1e-6 ||
      Math.abs(item.beam.theta_deg - 4) > 1e-9) {
    throw new Error(`${item.id}: offset is not exactly 4 degrees: local=${localOffset}; sky=${skyOffset}`);
  }
}
const westOffset = offsetCases.find(item => item.id === "parallel_sky_west");
const eastOffset = offsetCases.find(item => item.id === "parallel_sky_east");
if (Math.abs(westOffset.beam.direction_local[1]) > 1e-12 ||
    Math.abs(eastOffset.beam.direction_local[1]) > 1e-12 ||
    westOffset.beam.direction_local[0] <= 0 || eastOffset.beam.direction_local[0] >= 0 ||
    Math.abs(Math.abs(eastOffset.source_sky.az_deg) - 11.5550088042) > 1e-8) {
  throw new Error("East/West 4-degree offsets are not pure local ±u offsets with zenith-angle correction");
}
console.log("runtime OK all four parallel sources are exact 4-degree optical-axis offsets");

const completePointElements = run("?page=parallel");
const upPointClouds = completePointElements.sandbox.buildScene().pointClouds;
const upPointCounts = upPointClouds.map(cloud => cloud.points.length);
const upSummary = completePointElements.sandbox.selectedParallelCase().summary;
const upExpectedCounts = [
  upSummary.hit_mirror,
  upSummary.blocked_incoming,
  upSummary.blocked_reflected,
];
if (upPointCounts.join(",") !== upExpectedCounts.join(",")) {
  throw new Error(`parallel up full mirror/obstruction points are incomplete: actual=${upPointCounts}; expected=${upExpectedCounts}`);
}
completePointElements.get("camera").listeners.click({clientX: 300, clientY: 250});
const rightPointCounts = completePointElements.sandbox.buildScene().pointClouds
  .map(cloud => cloud.points.length);
const rightSummary = completePointElements.sandbox.selectedParallelCase().summary;
const rightExpectedCounts = [
  rightSummary.hit_mirror,
  rightSummary.blocked_incoming,
  rightSummary.blocked_reflected,
];
if (rightPointCounts.join(",") !== rightExpectedCounts.join(",")) {
  throw new Error(`parallel right full mirror/obstruction points are incomplete: actual=${rightPointCounts}; expected=${rightExpectedCounts}`);
}
console.log("runtime OK complete parallel reflection and obstruction point clouds");

const globalStructureElements = run("?page=global");
const parallelStructureElements = run("?page=parallel");
const globalFacets = globalStructureElements.sandbox.buildScene().polys
  .filter(polygon => polygon.name.startsWith("理想镜片 "));
const parallelFacets = parallelStructureElements.sandbox.buildScene().polys.slice(0, 54);
if (globalFacets.length !== 54 || parallelFacets.length !== 54) {
  throw new Error(`ideal mirror count mismatch: global=${globalFacets.length}, parallel=${parallelFacets.length}`);
}
let maxStructureVertexError = 0;
for (let facet = 0; facet < 54; ++facet) {
  for (let vertex = 0; vertex < globalFacets[facet].points.length; ++vertex) {
    const a = globalFacets[facet].points[vertex];
    const b = parallelFacets[facet].points[vertex];
    maxStructureVertexError = Math.max(
      maxStructureVertexError,
      Math.hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2])
    );
  }
}
if (maxStructureVertexError > 1e-9) {
  throw new Error(`global and parallel ideal mirrors disagree: ${maxStructureVertexError}`);
}
console.log("runtime OK global/parallel ideal mirror orientation", maxStructureVertexError);

const corsikaStructureElements = run("?page=corsika&case=event1909");
const corsikaStructureScene = corsikaStructureElements.sandbox.buildScene();
const corsikaFacets = corsikaStructureScene.polys
  .filter(polygon => polygon.name.startsWith("tel19 镜片 "));
const corsikaOrigin = corsikaStructureScene.marks.find(item => item.name === "tel19").point;
const parallelMirrorVertex = parallelStructureElements.sandbox.genericPoint(
  [0, 0, -16], parallelStructureElements.sandbox.currentParallelCase()
);
let maxParallelCorsikaStructureError = 0;
if (corsikaFacets.length !== 54) {
  throw new Error(`CORSIKA ideal mirror count mismatch: ${corsikaFacets.length}`);
}
for (let facet = 0; facet < 54; ++facet) {
  for (let vertex = 0; vertex < parallelFacets[facet].points.length; ++vertex) {
    const a = parallelFacets[facet].points[vertex].map((value, axis) =>
      value - parallelMirrorVertex[axis]
    );
    const bNwu = corsikaFacets[facet].points[vertex].map((value, axis) =>
      value - corsikaOrigin[axis]
    );
    // Compare physical vectors, not raw component arrays: generic is North-East-Up,
    // while the CORSIKA page preserves North-West-Up, so East = -West.
    const b = [bNwu[0], -bNwu[1], bNwu[2]];
    maxParallelCorsikaStructureError = Math.max(
      maxParallelCorsikaStructureError,
      Math.hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2])
    );
  }
}
if (maxParallelCorsikaStructureError > 1e-9) {
  throw new Error(`parallel/CORSIKA telescope-local geometry disagrees: ${maxParallelCorsikaStructureError}`);
}
console.log("runtime OK parallel/CORSIKA telescope-local mirror orientation", maxParallelCorsikaStructureError);

const elevationStructureElements = run("?page=elevation");
const parallelAxes = parallelStructureElements.sandbox.buildScene().lines
  .filter(item => item.name.startsWith("通用 local +"));
const elevationAxes = elevationStructureElements.sandbox.buildScene().lines
  .filter(item => item.name.startsWith("通用 local +"));
let maxParallelElevationAxisError = 0;
for (let axis = 0; axis < 3; ++axis) {
  for (let endpoint = 0; endpoint < 2; ++endpoint) {
    const a = parallelAxes[axis].points[endpoint];
    const b = elevationAxes[axis].points[endpoint];
    maxParallelElevationAxisError = Math.max(
      maxParallelElevationAxisError,
      Math.hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2])
    );
  }
}
if (parallelAxes.length !== 3 || elevationAxes.length !== 3 || maxParallelElevationAxisError > 1e-9) {
  throw new Error(`parallel/elevation trace frame mismatch: ${maxParallelElevationAxisError}`);
}
console.log("runtime OK parallel/elevation trace frame", maxParallelElevationAxisError);

const horizonElements = run("?page=global");
const horizonYs = [[-10, -10, 0], [-10, 10, 0], [10, -10, 0], [10, 10, 0]]
  .map(point => horizonElements.sandbox.project(point)[1]);
if (Math.max(...horizonYs) - Math.min(...horizonYs) > 1e-9) {
  throw new Error("global default ground plane is not projected to one horizon line");
}
console.log("runtime OK horizontal ground projection");

const elevationHorizonElements = run("?page=elevation");
const elevationHorizonYs = [[-10, -10, -17], [-10, 10, -17], [10, -10, -17], [10, 10, -17]]
  .map(point => elevationHorizonElements.sandbox.project(point)[1]);
if (Math.max(...elevationHorizonYs) - Math.min(...elevationHorizonYs) > 1e-9) {
  throw new Error("elevation default ground plane is not projected to one horizon line");
}
console.log("runtime OK elevation horizontal ground projection");

for (const query of ["?page=global", "?page=parallel", "?page=elevation", "?page=corsika&case=event1909"]) {
  const elements = run(query);
  const compass = elements.sandbox.buildScene().compass;
  if (!compass || !compass.frame) {
    throw new Error(`${query}: visible ground compass metadata is missing`);
  }
  console.log("runtime OK ground compass", query, compass.frame);

  const sceneData = elements.sandbox.buildScene();
  const layout = elements.sandbox.sceneLayout();
  const projectedYs = sceneData.points.map(point => elements.sandbox.project(point)[1]);
  if (Math.min(...projectedYs) < -1e-7 || Math.max(...projectedYs) > layout.usableHeight + 1e-7) {
    throw new Error(`${query}: fitted 3D content exceeds the full canvas`);
  }
}
console.log("runtime OK full-height 3D projection after compact compass restoration");

const fixedGlobal = run("?page=global").sandbox.fixedSkyReference();
const fixedParallel = run("?page=parallel").sandbox.fixedSkyReference();
const fixedElevation = run("?page=elevation").sandbox.fixedSkyReference();
const fixedCorsika = run("?page=corsika&case=event1909").sandbox.fixedSkyReference();
if (fixedGlobal.az_deg !== 0 || fixedGlobal.alt_deg !== 70 ||
    fixedParallel.az_deg !== 0 || fixedParallel.alt_deg !== 74 ||
    fixedElevation.az_deg !== 0 || fixedElevation.alt_deg !== 70 ||
    Math.abs(fixedCorsika.az_deg - 357.6867155597) > 1e-9 ||
    Math.abs(fixedCorsika.alt_deg - 70.6202124943) > 1e-9) {
  throw new Error("the unified fixed angle reference does not use each page's real sky direction");
}
console.log("runtime OK unified fixed sky-angle reference across pages");

for (const query of ["?page=global", "?page=parallel", "?page=elevation", "?page=corsika&case=event1909"]) {
  const elements = run(query);
  const forward = elements.sandbox.displayBasis().forward;
  if (!(forward[0] > 0)) {
    throw new Error(`${query}: default view must put North into the screen and South out of the screen`);
  }
}
console.log("runtime OK all four default views place South screen-out and North screen-in");

function topDirection(elements, direction) {
  elements.sandbox.setView("top");
  const a = elements.sandbox.rot([0, 0, 0]);
  const b = elements.sandbox.rot(direction);
  return [b[0] - a[0], b[1] - a[1]];
}
const genericTop = run("?page=parallel");
const genericNorth = topDirection(genericTop, [1, 0, 0]);
const genericEast = topDirection(genericTop, [0, 1, 0]);
if (!(genericNorth[1] > 0 && Math.abs(genericNorth[0]) < 1e-12 &&
      genericEast[0] > 0 && Math.abs(genericEast[1]) < 1e-12)) {
  throw new Error(`generic map top is not north-up/east-right: N=${genericNorth}, E=${genericEast}`);
}
const corsikaTop = run("?page=corsika&case=event1909");
const corsikaNorth = topDirection(corsikaTop, [1, 0, 0]);
const corsikaEast = topDirection(corsikaTop, [0, -1, 0]);
if (!(corsikaNorth[1] > 0 && Math.abs(corsikaNorth[0]) < 1e-12 &&
      corsikaEast[0] > 0 && Math.abs(corsikaEast[1]) < 1e-12)) {
  throw new Error(`CORSIKA map top is not north-up/east-right: N=${corsikaNorth}, E=${corsikaEast}`);
}
console.log("runtime OK map top north-up/east-right for generic and CORSIKA NWU");

const pylastCorsikaElements = run("?page=corsika&case=event1909");
const corsikaGeometryBefore = pylastCorsikaElements.sandbox.buildScene().polys
  .flatMap(poly => poly.points.flat());
pylastCorsikaElements.get("cameraCoords").value = "pylast";
pylastCorsikaElements.sandbox.cameraDraw();
const corsikaGeometryAfter = pylastCorsikaElements.sandbox.buildScene().polys
  .flatMap(poly => poly.points.flat());
let maxPylastGeometryMutation = 0;
for (let i = 0; i < corsikaGeometryBefore.length; i++) {
  maxPylastGeometryMutation = Math.max(
    maxPylastGeometryMutation,
    Math.abs(corsikaGeometryBefore[i] - corsikaGeometryAfter[i])
  );
}
if (maxPylastGeometryMutation > 0 ||
    !pylastCorsikaElements.get("cameraHead").textContent.includes("pyLAST 规则重画") ||
    !pylastCorsikaElements.get("cameraCaption").textContent.includes("同一 ROOT event 1909") ||
    pylastCorsikaElements.sandbox.pylastReaderCoordinates(0.2, -0.3).pix_x !== 0.3 ||
    pylastCorsikaElements.sandbox.pylastReaderCoordinates(0.2, -0.3).pix_y !== -0.2 ||
    pylastCorsikaElements.sandbox.pylastPointFromLact(0.2, -0.3).join(",") !== "-0.2,0.3") {
  throw new Error(`pyLAST browser redraw mutated 3D geometry or broke the audited mapping: ${maxPylastGeometryMutation}`);
}
console.log("runtime OK pyLAST browser redraw preserves CORSIKA 3D geometry and exact reader/renderer mapping");

const pylastParallelElements = run("?page=parallel");
pylastParallelElements.get("cameraCoords").value = "pylast";
pylastParallelElements.sandbox.cameraDraw();
if (!pylastParallelElements.get("cameraHead").textContent.includes("pyLAST 坐标") ||
    !pylastParallelElements.get("cameraCaption").textContent.includes("精确 offset=4°") ||
    !pylastParallelElements.get("cameraCaption").textContent.includes("绘图风格、暗色背景、强度和缩放与 LACT 像素图统一")) {
  throw new Error("parallel-light page does not redraw its real 4-degree output with unified pyLAST styling");
}
console.log("runtime OK parallel-light pyLAST redraw uses the real 4-degree cases and unified styling");

function maxBasisError(a, b) {
  return Math.max(...["x", "y", "z"].flatMap(key =>
    a[key].map((value, i) => Math.abs(value - b[key][i]))
  ));
}
const genericDisplayElements = run("?page=parallel");
const genericDisplay = genericDisplayElements.sandbox.displayBasis();
const genericProgram = genericDisplayElements.sandbox.genericBasis({az_deg: 0, el_deg: -28});
const genericDisplayError = maxBasisError(
  {x: genericDisplay.right, y: genericDisplay.up, z: genericDisplay.forward}, genericProgram
);
const corsikaDisplayElements = run("?page=corsika&case=event1909");
const corsikaDisplay = corsikaDisplayElements.sandbox.displayBasis();
const corsikaProgram = corsikaDisplayElements.sandbox.corsikaBasis({az_deg: 0, el_deg: -22});
const corsikaDisplayError = maxBasisError(
  {x: corsikaDisplay.right, y: corsikaDisplay.up, z: corsikaDisplay.forward},
  corsikaProgram
);
if (genericDisplayError > 1e-12 || corsikaDisplayError > 1e-12) {
  throw new Error(`display/program azimuth basis mismatch: generic=${genericDisplayError}, CORSIKA=${corsikaDisplayError}`);
}
const genericAz90 = genericDisplayElements.sandbox.genericBasis({az_deg: 90, el_deg: 0}).z;
const corsikaAz90 = corsikaDisplayElements.sandbox.corsikaBasis({az_deg: 90, el_deg: 0}).z;
if (Math.abs(genericAz90[0]) > 1e-12 || Math.abs(genericAz90[1] - 1) > 1e-12 ||
    Math.abs(corsikaAz90[0]) > 1e-12 || Math.abs(corsikaAz90[1] + 1) > 1e-12) {
  throw new Error(`az=90 must point East: generic=${genericAz90}, CORSIKA=${corsikaAz90}`);
}
console.log("runtime OK display basis follows program North-to-East azimuth", genericDisplayError, corsikaDisplayError);

const corsikaNorthFrame = corsikaDisplayElements.sandbox.corsikaBasis({az_deg: 0, el_deg: 70});
const sin70 = Math.sin(70 * Math.PI / 180);
const cos70 = Math.cos(70 * Math.PI / 180);
const expectedCorsikaNorthFrame = {
  x: [0, -1, 0],
  y: [-sin70, 0, cos70],
  z: [cos70, 0, sin70],
};
const corsikaNorthFrameError = maxBasisError(corsikaNorthFrame, expectedCorsikaNorthFrame);
if (corsikaNorthFrameError > 1e-12) {
  throw new Error(`CORSIKA HTML basis differs from latest C++ formula: ${corsikaNorthFrameError}`);
}
console.log("runtime OK CORSIKA HTML basis matches latest C++ numeric frame", corsikaNorthFrameError);

// The same physical source one degree East of a North-pointing telescope must
// enter both source adapters with the same telescope-local u/v signs.
const toLocal = (direction, basis) => [basis.x, basis.y, basis.z].map(axis =>
  direction.reduce((sum, value, index) => sum + value * axis[index], 0)
);
const el70 = 70 * Math.PI / 180;
const az1 = Math.PI / 180;
const genericEastSource = [Math.cos(el70) * Math.cos(az1), Math.cos(el70) * Math.sin(az1), Math.sin(el70)];
const corsikaEastSource = [genericEastSource[0], -genericEastSource[1], genericEastSource[2]];
const genericNorthFrame = genericDisplayElements.sandbox.genericBasis({az_deg: 0, el_deg: 70});
const genericEastLocal = toLocal(genericEastSource.map(value => -value), genericNorthFrame);
const corsikaEastLocal = toLocal(corsikaEastSource.map(value => -value), corsikaNorthFrame);
const horizontalEntryError = Math.max(...genericEastLocal.map((value, index) => Math.abs(value - corsikaEastLocal[index])));
if (horizontalEntryError > 1e-12 || genericEastLocal[0] >= 0) {
  throw new Error(`generic/CORSIKA East-offset local u mismatch: generic=${genericEastLocal}, CORSIKA=${corsikaEastLocal}`);
}
console.log("runtime OK generic and CORSIKA East-offset sources share local u/v signs", horizontalEntryError);

const showerAngles = corsikaDisplayElements.sandbox.corsikaAngleDisplay();
const expectedShowerAzimuth = ((showerAngles.array_rotation_deg - showerAngles.phi_deg + 180) % 360 + 360) % 360;
if (Math.abs(showerAngles.altitude_deg - (90 - showerAngles.theta_deg)) > 1e-12 ||
    Math.abs(showerAngles.derived_azimuth_deg - expectedShowerAzimuth) > 1e-12 ||
    Math.abs(showerAngles.azimuth_north_to_east_deg - showerAngles.derived_azimuth_deg) > 1e-5) {
  throw new Error(`CORSIKA theta/phi direction conversion mismatch: ${JSON.stringify(showerAngles)}`);
}
const beamPhi0 = corsikaDisplayElements.sandbox.beamDirectionFromAngles(1, 0);
const beamPhi90 = corsikaDisplayElements.sandbox.beamDirectionFromAngles(1, 90);
if (beamPhi0[0] <= 0 || Math.abs(beamPhi0[1]) > 1e-12 || beamPhi0[2] >= 0 ||
    Math.abs(beamPhi90[0]) > 1e-12 || beamPhi90[1] <= 0 || beamPhi90[2] >= 0) {
  throw new Error("synthetic beam theta/phi does not follow local -z, +x toward +y convention");
}
console.log("runtime OK CORSIKA shower and synthetic-beam theta/phi directions", showerAngles);

const scaleElements = run("?page=corsika&case=event1909");
const origin = scaleElements.sandbox.rot([0, 0, 0]);
const axisLengths = [[1, 0, 0], [0, 1, 0], [0, 0, 1]].map(point => {
  const projected = scaleElements.sandbox.rot(point);
  return Math.hypot(
    projected[0] - origin[0], projected[1] - origin[1], projected[2] - origin[2]
  );
});
if (Math.max(...axisLengths) - Math.min(...axisLengths) > 1e-12) {
  throw new Error("CORSIKA 3D projection does not use one isotropic x/y/z scale");
}
console.log("runtime OK isotropic CORSIKA x/y/z scale");
