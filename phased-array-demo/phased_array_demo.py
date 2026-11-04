#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phased Array Antenna Pattern Demo v25 (Fixed Matplotlib API & Radar Heatmap)
- Fixed AttributeError: 'QuadContourSet' object has no attribute 'collections'.
- Real 2D Spatial Radiation Mapping (Contourf) with 8-figure scanning.
"""
import sys, os
import matplotlib

try:
    import tkinter
    matplotlib.use('TkAgg')
    HAS_GUI = True
except ImportError:
    matplotlib.use('Agg')
    HAS_GUI = False

import numpy as np  # 【已修复】原先为 import numpy as n
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Button, TextBox
from matplotlib.font_manager import FontProperties, fontManager

# ---------------------------------------------------------------
# Font setup
# ---------------------------------------------------------------
def detect_font():
    cands = ['Microsoft YaHei', 'SimHei', 'SimSun', 'Microsoft JhengHei',
             'PingFang SC', 'Hiragino Sans GB', 'Heiti SC', 'STHeiti',
             'Noto Sans CJK SC', 'WenQuanYi Zen Hei', 'WenQuanYi Micro Hei']
    avail = {f.name for f in fontManager.ttflist}
    # 【已优化】将循环变量 n 改为 font_name 避免混淆
    for font_name in cands:
        if font_name in avail:
            try:
                return FontProperties(family=font_name)
            except:
                pass
    for p in [r'C:\Windows\Fonts\msyh.ttc', r'C:\Windows\Fonts\simhei.ttf',
              r'C:\Windows\Fonts\simsun.ttc',
              '/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc',
              '/usr/share/fonts/truetype/wqy/wqy-zenhei.ttc',
              '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf']:
        if os.path.exists(p):
            try:
                return FontProperties(fname=p)
            except:
                pass
    return FontProperties()


FP = detect_font()
print(f"[Info] Using font: {FP.get_name()}")

plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['font.size'] = 9

# ---------------------------------------------------------------
# Physical constants & state
# ---------------------------------------------------------------
c0 = 299_792_458.0
DEFAULTS = dict(M=16, N=32, dy=14.275831333333334e-3, dz=16.915e-3, f0=10.5e9, taper_idx=6)

TAPER_LABELS = [
    'Uniform', 'Hamming', 'Hann', 'Blackman',
    'Bartlett', 'Blackman-Harris', 'Kaiser(b=1.65)', 'Chebyshev(-25dB)',
    'Taylor(-25dB,appx)', 'Gaussian', 'Rayleigh(appx)', 'Binomial', 'Nuttall'
]

class State:
    def __init__(self):
        self.reset_all(hard=True)
        self.running = True
        self.hist_a = []
        self.hist_e = []
        self.hist_h = []
        self.hist_s = []

    def reset_all(self, hard=False):
        for k, v in DEFAULTS.items():
            setattr(self, k, v)
        self.need_rebuild = True
        self._rebuild()

    def _make_window(self, n, idx):
        if idx == 0:
            w = np.ones(n)
        elif idx == 1:
            w = np.hamming(n)
        elif idx == 2:
            w = np.hanning(n)
        elif idx == 3:
            w = np.blackman(n)
        elif idx == 4:
            w = np.bartlett(n)
        elif idx == 5:
            a = [0.35875, 0.48829, 0.14128, 0.01168]
            k = np.arange(n)
            w = a[0] - a[1] * np.cos(2 * np.pi * k / (n - 1)) + a[2] * np.cos(4 * np.pi * k / (n - 1)) - a[3] * np.cos(
                6 * np.pi * k / (n - 1))
        elif idx == 6:
            w = np.kaiser(n, 1.65)
        elif idx == 7:
            try:
                from scipy.signal import windows as _sw
                w = _sw.chebwin(n, at=25)
            except Exception:
                try:
                    from scipy.signal import chebwin as _cw
                    w = _cw(n, at=25)
                except Exception:
                    w = np.kaiser(n, 3.0)
        elif idx == 8:
            w = np.kaiser(n, 2.5)
        elif idx == 9:
            k = np.arange(n) - (n - 1) / 2
            s = 0.4 * (n - 1)
            w = np.exp(-0.5 * (k / s) ** 2)
        elif idx == 10:
            w = np.kaiser(n, 0.5)
        elif idx == 11:
            from math import comb
            w = np.array([comb(n - 1, k) for k in range(n)], dtype=float)
            w = np.maximum(w, 1e-30)
        elif idx == 12:
            a = [0.355768, 0.487396, 0.144232, 0.012604]
            k = np.arange(n)
            w = a[0] - a[1] * np.cos(2 * np.pi * k / (n - 1)) + a[2] * np.cos(4 * np.pi * k / (n - 1)) - a[3] * np.cos(
                6 * np.pi * k / (n - 1))
        else:
            w = np.ones(n)
        w = np.abs(w).astype(float)
        return w / w.max() if w.max() > 0 else w

    def _rebuild(self):
        self.lam = c0 / self.f0
        self.k0 = 2 * np.pi / self.lam
        mi = np.arange(self.M) - (self.M - 1) / 2
        ni = np.arange(self.N) - (self.N - 1) / 2
        self.y = mi * self.dy
        self.z = ni * self.dz
        self.wy = self._make_window(self.M, self.taper_idx)
        self.wz = self._make_window(self.N, self.taper_idx)
        self.need_rebuild = False

    def set_taper(self, i):
        self.taper_idx = i
        self.need_rebuild = True

    def set_M(self, v):
        try:
            v = max(2, min(128, int(round(float(v)))))
        except:
            v = self.M
        if v != self.M: self.M = v; self.need_rebuild = True
        return v

    def set_N(self, v):
        try:
            v = max(2, min(128, int(round(float(v)))))
        except:
            v = self.N
        if v != self.N: self.N = v; self.need_rebuild = True
        return v

    def set_dy(self, v):
        try:
            v = max(1.0, min(100.0, float(v))) / 1000.0
        except:
            v = self.dy
        if abs(v - self.dy) > 1e-10: self.dy = v; self.need_rebuild = True
        return v * 1000

    def set_dz(self, v):
        try:
            v = max(1.0, min(100.0, float(v))) / 1000.0
        except:
            v = self.dz
        if abs(v - self.dz) > 1e-10: self.dz = v; self.need_rebuild = True
        return v * 1000

    def set_f0(self, v):
        try:
            v = max(0.5, min(100.0, float(v))) * 1e9
        except:
            v = self.f0
        if abs(v - self.f0) > 1e3: self.f0 = v; self.need_rebuild = True
        return v / 1e9

    def apply(self):
        if self.need_rebuild:
            self._rebuild()
            self.hist_a.clear()
            self.hist_e.clear()
            self.hist_h.clear()
            self.hist_s.clear()


st = State()


# ---------------------------------------------------------------
# Pattern computation
# ---------------------------------------------------------------
def array_factor(th, ph, y, z, wy, wz, th0, ph0):
    ky = st.k0 * np.sin(th) * np.sin(ph)
    kz = st.k0 * np.sin(th) * np.cos(ph)
    ky0 = st.k0 * np.sin(th0) * np.sin(ph0)
    kz0 = st.k0 * np.sin(th0) * np.cos(ph0)
    dy_ = ky[..., None] - ky0
    dz_ = kz[..., None] - kz0
    Ay_c = (wy * np.cos(y * dy_)).sum(-1)
    Ay_s = (wy * np.sin(y * dy_)).sum(-1)
    Az_c = (wz * np.cos(z * dz_)).sum(-1)
    Az_s = (wz * np.sin(z * dz_)).sum(-1)
    return np.hypot(Ay_c * Az_c - Ay_s * Az_s, Ay_c * Az_s + Ay_s * Az_c)


def cut_curve(ph_deg, th0, ph0, n=361):
    th_d = np.linspace(-90, 90, n)
    th = np.deg2rad(th_d)
    ph = np.deg2rad(ph_deg)
    AF = array_factor(th, ph, st.y, st.z, st.wy, st.wz, th0, ph0)
    pat = np.abs(AF) * np.abs(np.cos(th))
    pat /= pat.max() + 1e-15
    return th_d, 20 * np.log10(pat + 1e-15)


def hpbw(th_d, pat_db):
    p = np.argmax(pat_db)
    pk = pat_db[p]
    idx = np.where(pat_db >= pk - 3.0)[0]
    return float(th_d[idx[-1]] - th_d[idx[0]]) if len(idx) >= 2 else np.nan


def peak_pos(th_d, pat_db):
    return float(th_d[np.argmax(pat_db)])


# ---------------------------------------------------------------
# 2D Radar Grid
# ---------------------------------------------------------------
r_line = np.linspace(0, 90, 101)  # Theta in degrees
ph_line = np.linspace(-90, 90, 181)  # Phi in degrees
R_MESH, PH_MESH = np.meshgrid(r_line, ph_line)
TH_RAD = np.deg2rad(R_MESH)
PH_RAD = np.deg2rad(PH_MESH)

# Project to Cartesian for plotting (X-Y plane)
X_MESH = np.sin(TH_RAD) * np.sin(PH_RAD)
Y_MESH = np.sin(TH_RAD) * np.cos(PH_RAD)

# ---------------------------------------------------------------
# Figure & Layout
# ---------------------------------------------------------------
FW, FH = 18.0, 11.5
fig = plt.figure(figsize=(FW, FH), dpi=100, facecolor='#f5f5f5')

fig.suptitle('Phased Array Antenna Pattern Demo | Beam Scanning Animation',
             fontsize=18, fontweight='bold', y=0.97)

# ---- Row 1: E / H / 2D Beam ----
PLOT_Y, PLOT_H = 0.66, 0.26
ax_e = fig.add_axes([0.04, PLOT_Y, 0.28, PLOT_H])
ax_h = fig.add_axes([0.37, PLOT_Y, 0.28, PLOT_H])
ax_beam = fig.add_axes([0.70, 0.33, 0.25, 0.50])

# ---- Row 2: Performance / Info / Polar ----
ax_perf = fig.add_axes([0.04, 0.32, 0.28, 0.25])
ax_info = fig.add_axes([0.36, 0.30, 0.16, 0.27])
ax_info.axis('off')
ax_pol = fig.add_axes([0.53, 0.32, 0.11, 0.23], projection='polar')

# ---- Taper buttons ----
TAP_Y_TOP = 0.25
TAP_HBTN = 0.022
TAP_DY = 0.024
TAP_W = 0.14
TAP_DX = 0.15
TAP_X0 = 0.04
taper_title = fig.text(TAP_X0 + (4 * TAP_DX - TAP_DX + TAP_W) / 2, TAP_Y_TOP + 0.005,
                       'Tapering Window Selection', ha='center', va='bottom',
                       fontsize=11, fontweight='bold')
taper_axes = []
taper_btns = []
for i, lbl in enumerate(TAPER_LABELS):
    col = i % 4
    row = i // 4
    x = TAP_X0 + col * TAP_DX
    y = TAP_Y_TOP - row * TAP_DY - TAP_HBTN
    ax_b = fig.add_axes([x, y, TAP_W, TAP_HBTN])
    b = Button(ax_b, lbl, color='#eeeeee', hovercolor='#dddddd')
    taper_axes.append(ax_b)
    taper_btns.append(b)
# Hint cell
ax_hint = fig.add_axes([TAP_X0 + 3 * TAP_DX, TAP_Y_TOP - 3 * TAP_DY - TAP_HBTN, TAP_W, TAP_HBTN])
ax_hint.text(0.5, 0.5, 'Keys: 1~9 / a~d', ha='center', va='center', fontsize=7, color='gray')
ax_hint.axis('off')

# ---- Bottom controls ----
TTL_Y = 0.13
INP_Y, INP_H = 0.08, 0.03
BTN_Y, BTN_H = 0.01, 0.03
LAM_Y = 0.06

ttl1 = fig.text(0.16, TTL_Y, 'Array Elements (M / N)', ha='center', fontsize=11, fontweight='bold')
ttl2 = fig.text(0.46, TTL_Y, 'Element Spacing dy/dz (mm)', ha='center', fontsize=11, fontweight='bold')
ttl3 = fig.text(0.81, TTL_Y, 'Frequency f0 / Wavelength', ha='center', fontsize=11, fontweight='bold')

# Input boxes
BW = 0.06
BH = INP_H
ax_m = fig.add_axes([0.06, INP_Y, BW, BH])
ax_n = fig.add_axes([0.18, INP_Y, BW, BH])
ax_dy = fig.add_axes([0.36, INP_Y, BW, BH])
ax_dz = fig.add_axes([0.48, INP_Y, BW, BH])
ax_f0 = fig.add_axes([0.76, INP_Y, 0.10, BH])

tb_m = TextBox(ax_m, 'M (Y):', initial=str(st.M), textalignment='center')
tb_n = TextBox(ax_n, 'N (Z):', initial=str(st.N), textalignment='center')
tb_dy = TextBox(ax_dy, 'dy(mm):', initial=f'{st.dy * 1000:.2f}', textalignment='center')
tb_dz = TextBox(ax_dz, 'dz(mm):', initial=f'{st.dz * 1000:.2f}', textalignment='center')
tb_f0 = TextBox(ax_f0, 'f0(GHz):', initial=f'{st.f0 / 1e9:.2f}', textalignment='center')

lam_text = fig.text(0.81, LAM_Y, f'lambda = {st.lam * 1000:.2f} mm',
                    ha='center', va='center', fontsize=11, color='#0066cc', fontweight='bold')

# Action buttons
BBW = 0.07
ax_p = fig.add_axes([0.47, BTN_Y, BBW, BTN_H])
ax_r = fig.add_axes([0.55, BTN_Y, BBW, BTN_H])
ax_x = fig.add_axes([0.63, BTN_Y, BBW, BTN_H])
ax_3 = fig.add_axes([0.71, BTN_Y, BBW, BTN_H])

btn_p = Button(ax_p, 'Pause', color='#ffe08a', hovercolor='#ffd060')
btn_r = Button(ax_r, 'Run', color='#a5d8a0', hovercolor='#7fc479')
btn_x = Button(ax_x, 'Reset', color='#d9d9d9', hovercolor='#bbbbbb')
btn_3 = Button(ax_3, 'Show 3D', color='#a4c2f4', hovercolor='#76a5ea')
status_text = fig.text(0.46, BTN_Y + BTN_H / 2, 'Status: Running', ha='right', va='center',
                       fontsize=11, color='darkgreen', fontweight='bold')


# ---------------------------------------------------------------
# Plot styling
# ---------------------------------------------------------------
def style_cart(ax, title, color, lbl, leg_loc):
    l, = ax.plot([], [], color=color, lw=2, label=lbl)
    vl = ax.axvline(0, color='g', ls='--', lw=1.4, alpha=0.8, label='Mainlobe')
    ax.set_xlim(-90, 90)
    ax.set_ylim(-40, 2)
    ax.grid(True, ls='--', alpha=0.6)
    ax.set_xlabel('theta (deg)', fontsize=10)
    ax.set_ylabel('Normalized Gain (dB)', fontsize=10)
    ax.legend(loc=leg_loc, fontsize=8)
    ax.set_title(title, fontsize=12, fontweight='bold', pad=6)
    return l, vl


line_e, vline_e = style_cart(ax_e, 'E-Plane Pattern (phi=0)', 'blue', 'E-Plane', leg_loc='upper right')
line_h, vline_h = style_cart(ax_h, 'H-Plane Pattern (phi=90)', 'red', 'H-Plane', leg_loc='upper right')

# Polar
line_pol, = ax_pol.plot([], [], 'm-', lw=1.4)
vline_pol, = ax_pol.plot([0, 0], [0, 1], 'g--', lw=1.4)
ax_pol.set_theta_zero_location('N')
ax_pol.set_thetamin(-90)
ax_pol.set_thetamax(90)
ax_pol.set_rlim(-40, 2)
ax_pol.set_rlabel_position(135)
ax_pol.set_rticks([-40, -30, -20, -10, 0])
ax_pol.set_xticks(np.deg2rad([-90, -60, -30, 0, 30, 60, 90]))
ax_pol.set_xticklabels(['-90', '-60', '-30', '0', '30', '60', '90'], fontsize=8)
ax_pol.grid(True, ls='--', alpha=0.6)
ax_pol.set_title('E-Plane Polar View', fontsize=11, fontweight='bold', pad=12)

# Performance plot
ax_perf.set_xlim(-90, 90)
ax_perf.set_ylim(-40, 22)
ax_perf.grid(True, ls='--', alpha=0.6)
ax_perf.axhline(-3, color='gray', ls=':', lw=1)
pe_e, = ax_perf.plot([], [], 'b-', lw=2, label='HPBW_E (deg)')
pe_h, = ax_perf.plot([], [], 'r-', lw=2, label='HPBW_H (deg)')
pe_s, = ax_perf.plot([], [], 'm-.', lw=2, label='SLL (dB)')
ax_perf.set_xlabel('Beam Point theta0 (deg)', fontsize=10)
ax_perf.set_ylabel('Value', fontsize=10)
ax_perf.legend(loc='upper left', fontsize=8, ncol=3)
ax_perf.set_title('Performance vs Scan Angle', fontsize=12, fontweight='bold', pad=6)

# 2D Solid Beam View (利用真实能量分布)
ax_beam.set_title('2D Spatial Beam Scanning', fontsize=12, fontweight='bold', pad=6)
ax_beam.set_facecolor('#0a0a1a')
ax_beam.set_xlim(-1.15, 1.15)
ax_beam.set_ylim(-0.1, 1.15)
ax_beam.set_xticks([-1.0, -0.5, 0.0, 0.5, 1.0])
ax_beam.set_yticks([0.0, 0.5, 1.0])
ax_beam.tick_params(colors='white', labelsize=10)
ax_beam.set_xlabel('Cross-Plane (X)', color='white', fontsize=11)
ax_beam.set_ylabel('Forward (Z)', color='white', fontsize=11)
ax_beam.grid(True, color='#333355', ls='--', alpha=0.5)
ax_beam.spines['bottom'].set_color('white')
ax_beam.spines['left'].set_color('white')
ax_beam.spines['top'].set_visible(False)
ax_beam.spines['right'].set_visible(False)

# 绘制雷达辅助半圆
for r in [0.3, 0.6, 0.9]:
    ax_beam.add_patch(plt.Circle((0, 0), r, color='#444466', fill=False, lw=1, ls='--'))
ax_beam.add_patch(plt.Circle((0, 0), 1.0, color='#555577', fill=False, lw=1.5))

# 添加角度刻度文字
for angle in [-60, -30, 0, 30, 60]:
    x_pos = 1.08 * np.sin(np.deg2rad(angle))
    y_pos = 1.08 * np.cos(np.deg2rad(angle))
    ax_beam.text(x_pos, y_pos, f'{angle}°', color='white', fontsize=10, fontweight='bold', ha='center', va='center')

# 天线原点标识
ax_beam.scatter([0], [0], c='white', s=80, marker='s', zorder=10, label='Array')

# 初始化热力图 contourf
radar_cf = None

# 动态扫描线
line_sweep, = ax_beam.plot([], [], 'w-', lw=1.5, alpha=0.8, zorder=7)

# 警告文字
warn_text = fig.text(0.81, 0.26, '', ha='center', va='top', color='red', fontsize=9, fontweight='bold')

# Info text
info_text = ax_info.text(0.02, 0.94, '', transform=ax_info.transAxes, fontsize=7.0,
                         va='top', family='monospace',
                         bbox=dict(boxstyle='round,pad=0.4', fc='white', alpha=0.95))
ax_info.set_title('Real-Time Parameters', fontsize=11, fontweight='bold', pad=10)


# ---------------------------------------------------------------
# UI callbacks
# ---------------------------------------------------------------
def highlight(i):
    for k, b in enumerate(taper_btns):
        b.color = '#b6e3b6' if k == i else '#eeeeee'
        b.ax.set_facecolor(b.color)


highlight(st.taper_idx)


def cb_pause(_):
    st.running = False
    status_text.set_text('Status: Paused')
    status_text.set_color('darkred')


def cb_run(_):
    st.apply()
    tb_m.set_val(str(st.M))
    tb_n.set_val(str(st.N))
    tb_dy.set_val(f'{st.dy * 1000:.2f}')
    tb_dz.set_val(f'{st.dz * 1000:.2f}')
    tb_f0.set_val(f'{st.f0 / 1e9:.2f}')
    lam_text.set_text(f'lambda = {st.lam * 1000:.2f} mm')
    st.running = True
    status_text.set_text('Status: Running')
    status_text.set_color('darkgreen')


def cb_reset(_):
    st.reset_all()
    highlight(st.taper_idx)
    tb_m.set_val(str(st.M))
    tb_n.set_val(str(st.N))
    tb_dy.set_val(f'{st.dy * 1000:.2f}')
    tb_dz.set_val(f'{st.dz * 1000:.2f}')
    tb_f0.set_val(f'{st.f0 / 1e9:.2f}')
    lam_text.set_text(f'lambda = {st.lam * 1000:.2f} mm')
    st.running = True
    status_text.set_text('Status: Running')
    status_text.set_color('darkgreen')


def cb_m(v):  v = st.set_M(v);  tb_m.set_val(str(v))


def cb_n(v):  v = st.set_N(v);  tb_n.set_val(str(v))


def cb_dy(v): v = st.set_dy(v); tb_dy.set_val(f'{v:.2f}')


def cb_dz(v): v = st.set_dz(v); tb_dz.set_val(f'{v:.2f}')


def cb_f0(v):
    v = st.set_f0(v)
    tb_f0.set_val(f'{v:.2f}')
    lam_text.set_text(f'lambda = {c0 / (v * 1e9) * 1000:.2f} mm')


tb_m.on_submit(cb_m)
tb_n.on_submit(cb_n)
tb_dy.on_submit(cb_dy)
tb_dz.on_submit(cb_dz)
tb_f0.on_submit(cb_f0)
btn_p.on_clicked(cb_pause)
btn_r.on_clicked(cb_run)
btn_x.on_clicked(cb_reset)


def taper_cb(i):
    def _f(_): st.set_taper(i); highlight(i)
    return _f


for i, b in enumerate(taper_btns):
    b.on_clicked(taper_cb(i))


def on_key(ev):
    k = (ev.key or '').lower()
    mp = {'1': 0, '2': 1, '3': 2, '4': 3, '5': 4, '6': 5, '7': 6, '8': 7, '9': 8, '0': 9, 'a': 10, 'b': 11, 'c': 12}
    if k in mp: st.set_taper(mp[k]); highlight(mp[k])


fig.canvas.mpl_connect('key_press_event', on_key)


def show_3d(_):
    f3 = plt.figure(figsize=(10, 8), dpi=100)
    ax3 = f3.add_subplot(111, projection='3d')
    th = np.linspace(0, np.pi / 2, 61)
    ph = np.linspace(-np.pi, np.pi, 121)
    TH, PH = np.meshgrid(th, ph)
    AF = array_factor(TH, PH, st.y, st.z, st.wy, st.wz, np.deg2rad(30), 0.0)
    pat = np.abs(AF) * np.abs(np.cos(TH))
    pat /= pat.max() + 1e-15
    R = pat
    X = R * np.sin(TH) * np.cos(PH)
    Y = R * np.sin(TH) * np.sin(PH)
    Z = R * np.cos(TH)
    G = 20 * np.log10(pat + 1e-15)
    ax3.plot_surface(X, Y, Z, facecolors=plt.cm.jet((G + 40) / 40), rstride=2, cstride=2, shade=False)
    ax3.set_xlabel('X')
    ax3.set_ylabel('Y')
    ax3.set_zlabel('Z')
    ax3.set_title(f'3D Pattern (theta0=+30deg, {TAPER_LABELS[st.taper_idx]})', fontweight='bold')
    ax3.set_xlim(-1, 1)
    ax3.set_ylim(-1, 1)
    ax3.set_zlim(0, 1)
    plt.show()


btn_3.on_clicked(show_3d)


# ---------------------------------------------------------------
# Animation update
# ---------------------------------------------------------------
def update(frame):
    global radar_cf
    if st.need_rebuild and st.running: st.apply()
    if not st.running:
        return (line_e, line_h, line_pol, vline_e, vline_h, vline_pol,
                info_text, pe_e, pe_h, pe_s, warn_text, lam_text, line_sweep)

    # 扫描范围：-40 到 40 (theta), -30 到 30 (phi)
    scan_th = 40 * np.sin(np.deg2rad(frame * 10 - 90))
    scan_ph = 30 * np.sin(np.deg2rad(frame * 7.5))

    th0 = np.deg2rad(scan_th)
    ph0 = np.deg2rad(scan_ph)

    the, pate = cut_curve(0.0, th0, ph0)
    bwe = hpbw(the, pate)
    pke = peak_pos(the, pate)
    line_e.set_data(the, pate)
    vline_e.set_xdata([pke, pke])

    thh, path = cut_curve(90.0, th0, ph0)
    bwh = hpbw(thh, path)
    pkh = peak_pos(thh, path)
    line_h.set_data(thh, path)
    vline_h.set_xdata([pkh, pkh])

    # Polar
    th_p = np.deg2rad(the)
    r_p = pate
    line_pol.set_data(th_p, r_p)
    vline_pol.set_data([np.deg2rad(pke), np.deg2rad(pke)], [-40, 0])

    # --- 2D Radar Heatmap Update (Contourf) ---
    AF_2d = array_factor(TH_RAD, PH_RAD, st.y, st.z, st.wy, st.wz, th0, ph0)
    pat_2d = np.abs(AF_2d) * np.cos(TH_RAD)
    pat_2d_db = 20 * np.log10(pat_2d / (pat_2d.max() + 1e-15) + 1e-15)

    # 修复：直接调用 remove() 方法，兼容新版 Matplotlib
    if radar_cf is not None:
        radar_cf.remove()

    radar_cf = ax_beam.contourf(X_MESH, Y_MESH, pat_2d_db, levels=np.linspace(-40, 0, 20), cmap='turbo', zorder=2)

    # 更新指向线（投射到2D平面的x和y坐标）
    line_sweep.set_data([0, np.sin(th0) * np.sin(ph0)], [0, np.sin(th0) * np.cos(ph0)])

    # SLL 计算
    if not np.isnan(bwe):
        m = (np.abs(the - pke) > bwe * 2)
        sll = float(np.max(pate[m])) if np.any(m) else -25.0
    else:
        sll = -25.0

    # history
    if len(st.hist_a) == 0 or abs(scan_th - st.hist_a[-1]) > 0.2:
        st.hist_a.append(scan_th)
        st.hist_e.append(bwe if not np.isnan(bwe) else np.nan)
        st.hist_h.append(bwh if not np.isnan(bwh) else np.nan)
        st.hist_s.append(sll)
        if len(st.hist_a) > 600:
            for L in (st.hist_a, st.hist_e, st.hist_h, st.hist_s): L.pop(0)
    pe_e.set_data(st.hist_a, st.hist_e)
    pe_h.set_data(st.hist_a, st.hist_h)
    pe_s.set_data(st.hist_a, st.hist_s)

    # Grating lobe warning
    hl = st.lam / 2
    parts = []
    if st.dy > hl: parts.append(f'dy={st.dy * 1000:.2f}mm > lambda/2={hl * 1000:.2f}mm')
    if st.dz > hl: parts.append(f'dz={st.dz * 1000:.2f}mm > lambda/2={hl * 1000:.2f}mm')

    if parts:
        warn_str = "WARNING: " + ", ".join(parts) + " -> grating lobes!"
    else:
        warn_str = ""
    warn_text.set_text(warn_str)

    # Info panel
    info = (
        "+---- Phased Array Params ----+\n"
        f"| Beam Point theta0 = {scan_th:+6.2f} deg |\n"
        f"| Beam Point phi    = {scan_ph:+6.2f} deg |\n"
        f"| Frequency  f0     = {st.f0 / 1e9:5.2f} GHz  |\n"
        f"| Wavelength lam    = {st.lam * 1000:5.2f} mm   |\n"
        f"| Array Size MxN    = {st.M} x {st.N}        |\n"
        f"| Spacing dy,dz     = {st.dy * 1000:.2f}, {st.dz * 1000:.2f} mm |\n"
        f"| Taper Window      = {TAPER_LABELS[st.taper_idx][:18]:<18s}|\n"
        "+---- Performance Metrics ---+\n"
        f"| E-plane HPBW      = {bwe:5.2f} deg   |\n"
        f"| H-plane HPBW      = {bwh:5.2f} deg   |\n"
        f"| E-peak / H-peak   = {pke:5.1f} / {pkh:5.1f}|\n"
        f"| 1st Sidelobe SLL  = {sll:5.1f} dB    |\n"
        "+-----------------------------+\n"
        "Green dash = mainlobe peak\n"
        "Keys 1~9 / a~d switch taper"
    )
    info_text.set_text(info)
    lam_text.set_text(f'lambda = {st.lam * 1000:.2f} mm')

    return (line_e, line_h, line_pol, vline_e, vline_h, vline_pol,
            info_text, pe_e, pe_h, pe_s, warn_text, lam_text, line_sweep)


# ---------------------------------------------------------------
# Launch
# ---------------------------------------------------------------
print("Phased array animation started. Close window to exit.")
ani = FuncAnimation(fig, update, frames=np.arange(0, 360, 1), interval=10,
                    blit=False, repeat=True)

if not HAS_GUI:
    update(45)
    fig.canvas.draw()
    fig.savefig('v25_real_radar_fixed.png', dpi=95, bbox_inches='tight')
    print("Saved: v25_real_radar_fixed.png")
    try:
        print("Saving GIF (may take a moment)...")
        ani.save('phased_array_v25.gif', writer='pillow', fps=30, dpi=75)
        print("Saved: phased_array_v25.gif")
    except Exception as e:
        print(f"GIF save skipped: {e}")
else:
    plt.show()