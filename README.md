<div align="center">
  <img src="https://github.com/user-attachments/assets/15c3e997-0e62-44f2-90a8-bdc786a2412d" alt="Phased Array Dashboard Preview" width="850">
  <p><i>Interactive Phased Array Antenna Dashboard</i></p>
</div>
----
# 📡 Phased Array Antenna Pattern Demo
An interactive 2D phased array antenna pattern simulation tool based on Python. Utilizing `matplotlib`'s animation engine, this program visually demonstrates the far-field radiation characteristics of phased array antennas under beam scanning, parameter tuning, and various amplitude taper (window) functions.
## ✨ Main Features
- **Dynamic Beam Scanning**: Automatically simulates the continuous scanning process of the antenna's main beam in space (Theta and Phi directions).
- **Multi-dimensional Visualization**:
  - 📊 **E-plane / H-plane Cartesian Plots**: Real-time display of beam shape, mainlobe, and sidelobe distributions.
  - 🎯 **Polar View**: Displays the E-plane beam pattern in polar coordinates.
  - 🗺️ **2D Spatial Radar Heatmap**: Uses contour plotting (`contourf`) to render the true 2D spatial radiation energy distribution, simulating a radar scanning screen.
  - 🌐 **3D Radiation Pattern**: Clicking "Show 3D" pops up an interactive, rotatable 3D radiation surface.
- **Real-time Interactive Parameter Tuning**:
  - Supports real-time modification via input boxes at the bottom of the interface: **Number of elements (M/N)**, **Element spacing (dy/dz)**, and **Operating frequency (f0)**.
  - The wavelength (\(\lambda\)) and radiation pattern are redrawn immediately upon parameter modification.
- **13 Taper (Amplitude Weighting) Window Functions**:
  - Built-in support for 13 window functions including Uniform, Hamming, Hann, Blackman, Kaiser, and Chebyshev.
  - Used to analyze the impact of different weighting methods on the **Half-Power Beamwidth (HPBW)** and **Sidelobe Level (SLL)**.
- **Smart Grating Lobe Warning**:
  - Monitors in real-time whether the element spacing is greater than half a wavelength (\(\lambda/2\)). If the condition is not met, a red warning pops up on the interface indicating the presence of grating lobes.
- **Real-time Performance Metrics Calculation**:
  - Dynamically calculates and plots the curves of HPBW and the first SLL as they vary with the scan angle.

## 🛠️ Requirements
------

# 📡 Phased Array Antenna Pattern Demo (相控阵天线方向图动态仿真)

这是一个基于 Python 的交互式二维相控阵天线方向图仿真工具。该程序利用 `matplotlib` 的动画引擎，直观地展示了相控阵天线在波束扫描、参数调节以及不同加权窗函数下的远场辐射特性。

## ✨ 主要功能

- **动态波束扫描 (Beam Scanning)**：自动模拟天线主波束在空间（Theta 和 Phi 方向）的连续扫描过程。
- **多维度可视化展示**：
  - 📊 **E面 / H面直角坐标图**：实时显示波束形状、主瓣与旁瓣分布。
  - 🎯 **极坐标图 (Polar View)**：以极坐标形式展示 E 面波束方向图。
  - 🗺️ **2D 空间雷达热力图**：利用等高线图 (`contourf`) 绘制真实的二维空间辐射能量分布，模拟雷达扫描屏效果。
  - 🌐 **3D 立体方向图**：点击 "Show 3D" 可弹出可拖拽旋转的 3D 辐射曲面。
- **实时参数交互调节**：
  - 支持在界面下方输入框实时修改：**阵元数量 (M/N)**、**阵元间距 (dy/dz)**、**工作频率 (f0)**。
  - 修改参数后，波长（λ）及方向图会立即重绘。
- **13 种窗函数（幅度加权）切换**：
  - 内置 Uniform, Hamming, Hann, Blackman, Kaiser, Chebyshev 等 13 种窗函数。
  - 用于分析不同加权方式对**主瓣宽度 (HPBW)** 和 **旁瓣电平 (SLL)** 的影响。
- **智能栅瓣警告 (Grating Lobe Warning)**：
  - 实时监测阵元间距是否大于半波长 (\(\lambda/2\))，若条件不满足，界面会弹出红色警告，提示出现栅瓣。
- **实时性能指标计算**：
  - 动态计算并绘制半功率波束宽度 (HPBW) 和第一旁瓣电平 (SLL) 随扫描角度变化的曲线。

## 🛠️ 环境依赖

运行此代码需要 Python 3.7+ 环境。请确保安装以下依赖库：

```bash
pip install numpy matplotlib scipy


Running this code requires a Python 3.7+ environment. Please ensure the following dependencies are installed:

```bash
pip install numpy matplotlib scipy
