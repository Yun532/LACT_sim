# LACT 坐标系说明

建议先打开[三维坐标网页](assets/lact-coordinate-system-3d.html)，再对照[逐代码审计笔记](coordinate_code_audit_zh.md)。

最短的数据链是：

```text
CORSIKA NWU (+N,+W,+Up)
  -> mirror/obstruction local (+x=West at az=0, +y=sky-up, +z=boresight)
  -> LACT focal plane (+u=-mirror x=East, +v=mirror y)
  -> ROOT camera x_m/y_m = u/v
  -> pyLAST pix_x=-v, pix_y=-u
```

望远镜方位角从北开始，顺时针向东为正。网页四页都使用同一 NWU 地图：俯视时上北、下南、左西、右东；三维旋转只改变观察相机，不改变这些物理方向。

平行光、不同仰角和 CORSIKA 相机图均直接读取本次程序输出。LACT 模式画原始 `u/v`；pyLAST 模式只执行 reader 的 `(-v,-u)` 和 pyLAST 自己的画布轴顺序，不重新模拟、不移动光斑中心。

镜片/遮挡 `x/y` 与相机 `u/v` 不能混叫成一套轴：前者是相机看镜面，后者是镜片看相机，因此 `u=-mirror x`。这也是输出面必须满足 `u×v=normal` 的原因。
