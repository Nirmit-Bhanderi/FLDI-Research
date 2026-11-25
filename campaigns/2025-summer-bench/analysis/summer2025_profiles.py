"""Sensitivity profiles for the summer 2025 FLDI bench campaign (hot-air gun).

Two modes:

1. From the raw recordings (one t<RNo>.mat per run, written by
   acquisition/DAQ_FLDI_Setup.m), recompute the per-position statistics:

       python analysis/summer2025_profiles.py --raw "/path/to/FLDI Expts"

   The raw folder must contain the day folders "3 June 2025", "4 June 2025"
   and "12 June 2025". The raw data is not stored in this repository.

2. Without --raw, re-plot from the processed CSV files already in
   data/processed/:

       python analysis/summer2025_profiles.py

Outputs: data/processed/traverse_stats.csv and
figures/summer2025_overview.png.

Requires numpy, scipy and matplotlib.
"""
import argparse
import csv
import glob
import os
import re

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROC = os.path.join(ROOT, "data", "processed")
FIG = os.path.join(ROOT, "figures", "summer2025_overview.png")
FS = 10_000.0  # Hz, sampling rate used on 3, 4 and 12 June

# Fringe extremes from each day's Wollaston test (runs 1000-1002), in volts.
# Near quadrature the photodiode responds (Vmax - Vmin)/2 volts per radian, so
# phase_rms = 2 * V_rms / (Vmax - Vmin).
# 12 June: the log sheet lists Vmin as 0.124 V, but the recordings show
# 0.0124 V (a dropped zero); the corrected value is used here.
FRINGE = {
    "3 June 2025": (0.2362, 0.3309),
    "4 June 2025": (0.2340, 0.3316),
    "12 June 2025": (0.0124, 0.3759),
}
COLOURS = {"3 June 2025": "#9a6200", "4 June 2025": "#1d6a86", "12 June 2025": "#c42b23"}


def position_mm(label):
    """'30 mm (diode)' -> +30, '12 mm (laser)' -> -12, anything else -> None."""
    m = re.match(r"([\d.]+) mm \((laser|diode)\)", label)
    if not m:
        return None
    return float(m.group(1)) * (-1 if m.group(2) == "laser" else 1)


def stats_from_raw(raw):
    from scipy.io import loadmat
    from scipy.signal import welch

    rows = []
    for day in FRINGE:
        per_z = {}
        for path in glob.glob(os.path.join(raw, day, "t*.mat")):
            m = loadmat(path)
            run = int(m["RNo"].ravel()[0])
            z = position_mm(str(m["loc"][0]) if m["loc"].size else "")
            # Runs >= 1000 are Wollaston tests / bias checks. t1004.mat on
            # 4 June carries a stale '0 mm (laser)' label, so filter by run number.
            if z is None or run >= 1000:
                continue
            v = m["data"].ravel()
            v = v - v.mean()
            f, p = welch(v, fs=FS, nperseg=2048)
            df = f[1] - f[0]
            per_z.setdefault(z, []).append(
                (v.var(), p[(f >= 5) & (f < 500)].sum() * df, p[f >= 1000].sum() * df)
            )
        vmin, vmax = FRINGE[day]
        for z in sorted(per_z):
            a = np.array(per_z[z])
            rms = np.sqrt(a[:, 0])
            rows.append({
                "day": day,
                "position_mm": z,
                "n_runs": len(a),
                "rms_V": rms.mean(),
                "rms_V_std": rms.std(),
                "phase_rms_rad": 2 * rms.mean() / (vmax - vmin),
                "phase_rms_below_500Hz_rad": 2 * np.sqrt(a[:, 1].mean()) / (vmax - vmin),
                "phase_rms_above_1kHz_rad": 2 * np.sqrt(a[:, 2].mean()) / (vmax - vmin),
            })
    out = os.path.join(PROC, "traverse_stats.csv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()})
    print("wrote", out)


def read_csv(name):
    with open(os.path.join(PROC, name), newline="") as fh:
        return list(csv.DictReader(fh))


def fwhm(z, y):
    """Full width at half maximum of y(z) by linear interpolation (nan if not bracketed)."""
    i = int(np.argmax(y))
    h = y[i] / 2
    if y[: i + 1].min() >= h or y[i:].min() >= h:
        return float("nan")
    left = np.interp(h, y[: i + 1], z[: i + 1])
    right = np.interp(h, y[i:][::-1], z[i:][::-1])
    return right - left


def plot():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    stats = read_csv("traverse_stats.csv")
    vert = read_csv("2025-06-30_vertical_traverse.csv")
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
    for day, col in COLOURS.items():
        d = [r for r in stats if r["day"] == day]
        z = np.array([float(r["position_mm"]) for r in d])
        vmin, vmax = FRINGE[day]
        ax[0].plot(z, [1e3 * float(r["rms_V"]) for r in d], "o-", ms=3.5, color=col,
                   label=f"{day} (fringe half-amplitude {1e3 * (vmax - vmin) / 2:.0f} mV)")
        hi = np.array([float(r["phase_rms_above_1kHz_rad"]) for r in d])
        ax[1].plot(z, hi, "o-", ms=3.5, color=col, label=day)
        ax[1].plot(z, [float(r["phase_rms_below_500Hz_rad"]) for r in d], "--", lw=1.2, color=col)
        if day == "12 June 2025":
            var = np.array([float(r["rms_V"]) ** 2 for r in d])
            print(f"12 June: variance FWHM = {fwhm(z, var):.1f} mm, "
                  f"above-1 kHz variance FWHM = {fwhm(z, hi ** 2):.1f} mm")
    ax[0].set_title("Raw RMS voltage vs traverse position")
    ax[0].set_ylabel("RMS (mV)")
    ax[0].legend(fontsize=8, frameon=False)
    ax[1].set_title("Phase RMS (rad): solid = above 1 kHz, dashed = below 500 Hz")
    ax[1].set_ylabel("rad")
    for a in ax[:2]:
        a.set_xlabel("position along beam (mm)   - laser side  |  + diode side")
        a.grid(alpha=0.3)

    # 30 June: average the repeats at each vertical position.
    pos = {}
    for r in vert:
        pos.setdefault(float(r["vertical_position_mm"]), []).append(float(r["rms_V"]))
    y = sorted(pos)
    ax[2].plot(y, [1e3 * np.mean(pos[k]) for k in y], "o-", ms=3.5, color="#c42b23")
    ax[2].set_title("30 June: vertical traverse across the hot jet\n"
                    "(10 mm from nozzle, 200 °C; RMS values from the log sheet)")
    ax[2].set_xlabel("vertical position (mm)")
    ax[2].set_ylabel("RMS (mV)")
    ax[2].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIG, dpi=150)
    print("wrote", FIG)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--raw", help='folder containing "3 June 2025", "4 June 2025", "12 June 2025"')
    args = ap.parse_args()
    if args.raw:
        stats_from_raw(args.raw)
    plot()
