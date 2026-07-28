const fs = require("fs");
const vm = require("vm");

const html = fs.readFileSync("docs/assets/lact-coordinate-system-3d.html", "utf8");
if (!html.includes('<details id="infoPanel" class="panel">') ||
    html.includes('<details id="infoPanel" class="panel" open>')) {
  throw new Error("left information panel must exist and be collapsed by default");
}
for (const marker of [
  '<details class="textfold"><summary>页面说明</summary>',
  '<details class="textfold"><summary>图像说明</summary>',
  '<details class="textfold"><summary>图表说明</summary>',
  '固定地图：北上东右',
  '当前屏幕投影（随视角旋转）',
  'cctx.fillStyle="#0b1925"',
  'const PARALLEL_COLORS=["#ff6d7a","#64c7ff","#ffd45e","#70e39c"]',
  'caseDensityColor(plotColor',
  'id="cameraCoords"',
  'cameraCoordinateMode()==="pylast"?[-p[1],-p[0]]',
  'pyLAST +pix_x = -u',
  'pyLAST +pix_y = -v',
  'pyLAST 画布：向右 = -v；向上 = -u',
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
  const initialValues = {cameraScale: "spot", cameraCoords: "lact", globalAz: "0", globalEl: "70", photonTel: "19"};
  return {
    id, dataset, style: {}, classList: {add() {}, remove() {}, toggle() {}},
    value: initialValues[id] || "0", hidden: false,
    width: 900, height: 700, innerHTML: "", textContent: "", max: 0,
    listeners: {}, addEventListener(type, handler) { this.listeners[type] = handler; }, setPointerCapture() {},
    getBoundingClientRect() { return {left: 0, top: 0, width: id === "scene" ? 900 : 400, height: id === "scene" ? 700 : 336}; },
    getContext() { return context; },
  };
}

function run(search) {
  const ids = ["scene", "camera", "info", "cameraPanel", "cameraPanelTitle", "cameraScale", "cameraCoords", "cameraHead",
    "cameraCaption", "deformPanel", "deformChart", "deformStat", "colorMin", "colorMax",
    "globalToolbar", "parallelToolbar", "elevationToolbar", "corsikaToolbar", "globalAz",
    "globalAzValue", "globalEl", "globalElValue", "angleValue", "bunch", "photonTel",
    "bunchValue", "elevation", "arrowScale", "corsikaCase", "groundView", "showerView",
    "zoomOut", "zoomIn", "subtitle", "legend"];
  const elements = new Map(ids.map(id => [id, element(id)]));
  const pageButtons = ["global", "parallel", "elevation", "corsika"].map(p => element("page-" + p, {page: p}));
  const viewButtons = ["horizon", "iso", "top", "north", "east"].map(view => element("view-" + view, {view}));
  const document = {
    documentElement: {},
    querySelector(selector) {
      if (selector.startsWith("#")) return elements.get(selector.slice(1));
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
  ["?page=corsika&case=event1909", "zem 原值"],
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
for (const expected of ["pyLAST +pix_x = -u", "pyLAST +pix_y = -v"]) {
  if (!globalAxisNames.includes(expected)) {
    throw new Error(`global definition scene is missing ${expected}`);
  }
}
const globalInfo = globalElements.get("info").innerHTML;
for (const expected of ["pix_x=-u", "横轴=pix_y=-v", "LactEventSource.cpp:245-246", "visualize.py:82-83"]) {
  if (!globalInfo.includes(expected)) {
    throw new Error(`global pyLAST definition is missing ${expected}`);
  }
}
console.log("runtime OK global page contains pyLAST field and canvas definitions");

const elevationWhiteboardElements = run("?page=elevation");
const elevationCase = elevationWhiteboardElements.sandbox.currentParallelCase();
const elevationHistogram = elevationWhiteboardElements.sandbox.spotHistogram(elevationCase, 140);
if (!elevationWhiteboardElements.get("cameraHead").textContent.includes("程序同风格完整光斑") ||
    elevationHistogram.n <= 0 ||
    elevationHistogram.n !== elevationCase.full_output_uv_m.length ||
    elevationHistogram.n !== elevationCase.validation.full_output_uv_rows ||
    elevationHistogram.n !== elevationCase.summary.hit_output_plane) {
  throw new Error("elevation page does not reuse the complete program-style whiteboard output");
}
console.log("runtime OK elevation uses complete program-style whiteboard output");

const corsikaFullCameraElements = run("?page=corsika&case=event1909");
if (corsikaFullCameraElements.get("cameraScale").value !== "full" ||
    !corsikaFullCameraElements.get("cameraHead").textContent.includes("完整相机")) {
  throw new Error("CORSIKA page must default to the complete physical camera view");
}
const groundPointBeforeZoom = corsikaFullCameraElements.sandbox.project([0, 0, 0]);
corsikaFullCameraElements.sandbox.zoomAt(120, groundPointBeforeZoom[0], groundPointBeforeZoom[1]);
const groundPointAfterZoom = corsikaFullCameraElements.sandbox.project([0, 0, 0]);
if (Math.hypot(groundPointAfterZoom[0] - groundPointBeforeZoom[0],
               groundPointAfterZoom[1] - groundPointBeforeZoom[1]) > 1e-7) {
  throw new Error("cursor-centered zoom failed to keep the selected ground point fixed");
}
console.log("runtime OK CORSIKA defaults to full cameras and can zoom into the ground anchor");

const parallelElements = run("?page=parallel");
const fullDensityColors = ["#ff6d7a", "#64c7ff", "#ffd45e", "#70e39c"]
  .map(color => parallelElements.sandbox.caseDensityColor(color, 1));
if (fullDensityColors.join("|") !== "rgb(255,109,122)|rgb(100,199,255)|rgb(255,212,94)|rgb(112,227,156)") {
  throw new Error(`parallel density colors do not match ray colors: ${fullDensityColors}`);
}
const parallelCase = parallelElements.sandbox.currentParallelCase();
const rawCentroid = parallelCase.camera_summary.output_uv_centroid_m;
parallelElements.get("cameraCoords").value = "pylast";
const pylastPoint = parallelElements.sandbox.cameraDisplayPoint([0.12, -0.34]);
const pylastHistogram = parallelElements.sandbox.spotHistogram(parallelCase, 140);
if (Math.abs(pylastPoint[0] - 0.34) > 1e-12 || Math.abs(pylastPoint[1] + 0.12) > 1e-12 ||
    Math.abs(pylastHistogram.centroid[0] + rawCentroid[1]) > 1e-12 ||
    Math.abs(pylastHistogram.centroid[1] + rawCentroid[0]) > 1e-12 ||
    pylastHistogram.n !== parallelCase.full_output_uv_m.length) {
  throw new Error("pyLAST camera view does not implement horizontal=-v, vertical=-u on the complete point set");
}
parallelElements.get("cameraCoords").value = "lact";
parallelElements.get("camera").listeners.click({clientX: 300, clientY: 250});
if (!parallelElements.get("info").innerHTML.includes("az=+1°")) {
  throw new Error("parallel camera click did not select the bottom-right sky case");
}
if (!parallelElements.get("cameraCaption").textContent.includes("中心不是 (0,0)") ||
    !parallelElements.get("cameraCaption").textContent.includes("绝对 u/v 原值")) {
  throw new Error("parallel whiteboard plots do not preserve absolute u/v ticks while auto-framing");
}
console.log("runtime OK parallel whiteboard selection and auto-framed absolute axes");

const completePointElements = run("?page=parallel");
const upPointClouds = completePointElements.sandbox.buildScene().pointClouds;
const upPointCounts = upPointClouds.map(cloud => cloud.points.length);
if (upPointCounts.join(",") !== "14803,2973,127") {
  throw new Error(`parallel up full mirror/obstruction points are incomplete: ${upPointCounts}`);
}
completePointElements.get("camera").listeners.click({clientX: 300, clientY: 250});
const rightPointCounts = completePointElements.sandbox.buildScene().pointClouds
  .map(cloud => cloud.points.length);
if (rightPointCounts.join(",") !== "14795,2904,121") {
  throw new Error(`parallel right full mirror/obstruction points are incomplete: ${rightPointCounts}`);
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
}

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
    !pylastCorsikaElements.get("cameraHead").textContent.includes("pyLAST")) {
  throw new Error(`pyLAST camera display mutated 3D geometry or failed to label itself: ${maxPylastGeometryMutation}`);
}
console.log("runtime OK pyLAST camera view preserves CORSIKA 3D geometry");

function maxBasisError(a, b) {
  return Math.max(...["x", "y", "z"].flatMap(key =>
    a[key].map((value, i) => Math.abs(value - b[key][i]))
  ));
}
const genericDisplayElements = run("?page=parallel");
const genericDisplay = genericDisplayElements.sandbox.displayBasis();
const genericProgram = genericDisplayElements.sandbox.genericBasis({az_deg: 225, el_deg: -28});
const genericDisplayError = maxBasisError(
  {x: genericDisplay.right, y: genericDisplay.up, z: genericDisplay.forward}, genericProgram
);
const corsikaDisplayElements = run("?page=corsika&case=event1909");
const corsikaDisplay = corsikaDisplayElements.sandbox.displayBasis();
const corsikaProgram = corsikaDisplayElements.sandbox.corsikaBasis({az_deg: 220, el_deg: -22});
const corsikaDisplayError = maxBasisError(
  {x: corsikaDisplay.up, y: corsikaDisplay.right, z: corsikaDisplay.forward}, corsikaProgram
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
