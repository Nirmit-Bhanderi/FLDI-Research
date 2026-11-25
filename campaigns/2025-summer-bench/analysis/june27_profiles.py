"""27 June 2025: FLDI sensitivity with the hot-air gun, sampled at 250 kHz.

Four traverses were recorded that day (folder names as on the lab drive):
  Horizontal - T1, T2, T3 : gun moved along the optical axis
  Vertical - T1           : gun moved vertically across the beam near the focus
Gun nozzle 10 mm from the beam, 200 C. Positions come from the log sheet
(27th june.xlsx); the 'loc' labels saved in the .mat files are not used,
because in Horizontal - T1 they kept saying 'laser' after the traverse
crossed to the photodiode side.

Two modes:

  python analysis/june27_profiles.py --raw "/path/to/27 June 2025"
      recompute data/processed/2025-06-27_*.csv from the raw recordings
      (needs the .mat files and 27th june.xlsx; about 2 GB, not in the repo)

  python analysis/june27_profiles.py
      re-plot figures/june27_overview.png from the processed CSVs

Requires numpy, scipy, matplotlib and (for --raw) openpyxl.
"""
import argparse
import csv
import glob
import os

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROC = os.path.join(ROOT, "data", "processed")
FIG = os.path.join(ROOT, "figures", "june27_overview.png")

# log sheet tab -> data folder
SHEETS = {
    "Horizontal run 1": "Horizontal - T1",
    "Horizontal run 2": "Horizontal - T2",
    "Horizontal run 3": "Horizontal - T3",
    "Vertical run 1": "Vertical - T1",
}
# Fringe extremes (V) from the Wollaston sweeps recorded that day (t1000-t1001).
# Horizontal - T2 has no sweep files; it used the same alignment as T1.
FRINGE = {
    "Horizontal - T1": (0.0053, 0.0637),
    "Horizontal - T2": (0.0053, 0.0637),
    "Horizontal - T3": (0.0052, 0.0648),
    "Vertical - T1": (0.0052, 0.0648),
}
BANDS = [(5, 1e3), (1e3, 1e4), (1e4, 5e4), (5e4, 1.25e5)]  # Hz
BAND_NAMES = ["below 1 kHz", "1-10 kHz", "10-50 kHz", "50-125 kHz"]
MAP_FOLDER = "Horizontal - T2"  # widest traverse; used for the frequency map
MAP_EDGES = np.geomspace(100, 125e3, 61)


def positions_from_log(xlsx):
    """{(folder, run): position_mm} from the 'R no | location (mm) | V' table of each tab."""
    import openpyxl

    wb = openpyxl.load_workbook(xlsx, data_only=True)
    pos = {}
    for sheet, folder in SHEETS.items():
        loc, started = None, False
        for row in wb[sheet].iter_rows(values_only=True):
            c = row[11:14]
            if c[0] == "R no" and c[1] == "location (mm)":
                started = True
                continue
            if started and isinstance(c[0], (int, float)):
                if c[1] is not None:
                    loc = float(c[1])
                pos[(folder, int(c[0]))] = loc
    return pos


def from_raw(raw):
    from scipy.io import loadmat
    from scipy.signal import welch

    pos = positions_from_log(os.path.join(raw, "27th june.xlsx"))
    rows, map_rows = [], []
    for folder, (vmin, vmax) in FRINGE.items():
        span = vmax - vmin  # volts per 2 rad near quadrature
        mid = (vmin + vmax) / 2
        spectra, means, noise = {}, {}, []
        for path in glob.glob(os.path.join(raw, folder, "t*.mat")):
            m = loadmat(path)
            run = int(m["RNo"].ravel()[0])
            v = m["data"].ravel()
            f, p = welch(v - v.mean(), fs=float(m["fs"].ravel()[0]), nperseg=16384)
            if run <= 5:  # gun off
                noise.append(p)
            elif run < 1000 and pos.get((folder, run)) is not None:
                z = pos[(folder, run)]
                spectra.setdefault(z, []).append(p)
                means.setdefault(z, []).append(v.mean())
        df = f[1] - f[0]
        p_off = np.mean(noise, axis=0)
        for z in sorted(spectra):
            p = np.clip(np.mean(spectra[z], axis=0) - p_off, 0, None)  # remove instrument noise
            band_var = [p[(f >= a) & (f < b)].sum() * df for a, b in BANDS]
            mean_v = float(np.mean(means[z]))
            rows.append({
                "traverse": folder,
                "position_mm": z,
                "n_runs": len(spectra[z]),
                "phase_rms_rad": 2 * np.sqrt(sum(band_var)) / span,
                **{f"phase_rms_{n.replace(' ', '_')}_rad": 2 * np.sqrt(b) / span
                   for n, b in zip(BAND_NAMES, band_var)},
                "mean_V": mean_v,
                # mean bias phase relative to mid-fringe; shifts when the gun changes the mean density gradient
                "mean_phase_shift_rad": float(np.arcsin(np.clip((mean_v - mid) / (span / 2), -1, 1))),
            })
            if folder == MAP_FOLDER:
                for lo, hi in zip(MAP_EDGES[:-1], MAP_EDGES[1:]):
                    sel = (f >= lo) & (f < hi)
                    if sel.any():
                        map_rows.append({"position_mm": z, "f_Hz": np.sqrt(lo * hi),
                                         "psd_rad2_per_Hz": p[sel].mean() * (2 / span) ** 2})
        rows.append({"traverse": folder, "position_mm": "gun off",
                     "n_runs": len(noise),
                     "phase_rms_rad": 2 * np.sqrt(p_off[f >= 5].sum() * df) / span})
    write_csv("2025-06-27_traverse_stats.csv", rows)
    write_csv("2025-06-27_T2_frequency_map.csv", map_rows)


def write_csv(name, rows):
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    with open(os.path.join(PROC, name), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()})
    print("wrote", os.path.join(PROC, name))


def read_csv(name):
    with open(os.path.join(PROC, name), newline="") as fh:
        return list(csv.DictReader(fh))


def fwhm(z, y):
    i = int(np.argmax(y))
    h = y[i] / 2
    if y[: i + 1].min() >= h or y[i:].min() >= h:
        return float("nan")
    return np.interp(h, y[i:][::-1], z[i:][::-1]) - np.interp(h, y[: i + 1], z[: i + 1])


def plot():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    stats = [r for r in read_csv("2025-06-27_traverse_stats.csv") if r["position_mm"] != "gun off"]
    noise = {r["traverse"]: float(r["phase_rms_rad"])
             for r in read_csv("2025-06-27_traverse_stats.csv") if r["position_mm"] == "gun off"}

    def series(folder, key):
        d = [r for r in stats if r["traverse"] == folder]
        return (np.array([float(r["position_mm"]) for r in d]),
                np.array([float(r[key]) for r in d]))

    fig, ax = plt.subplots(2, 2, figsize=(13, 9))

    a = ax[0, 0]
    for folder, col in zip(["Horizontal - T1", "Horizontal - T2", "Horizontal - T3"],
                           ["#c42b23", "#1d6a86", "#9a6200"]):
        z, phi = series(folder, "phase_rms_rad")
        a.plot(z, phi, "o-", ms=3.5, color=col,
               label=f"{folder} (variance FWHM {fwhm(z, phi ** 2):.0f} mm)")
    a.axhline(np.mean(list(noise.values())), color="gray", ls=":", lw=1, label="gun off (instrument noise)")
    a.set_title("Along the beam: phase RMS vs gun position\n(each traverse has its own scale zero)")
    a.set_xlabel("position along beam (mm)   - laser side  |  + photodiode side")
    a.set_ylabel("phase RMS (rad)")
    a.legend(fontsize=8, frameon=False)

    a = ax[0, 1]
    # 50-125 kHz is left out: after removing the gun-off noise it is at the noise level.
    for name, col in zip(BAND_NAMES[:3], ["#9a6200", "#c42b23", "#1d6a86"]):
        z, phi = series(MAP_FOLDER, f"phase_rms_{name.replace(' ', '_')}_rad")
        v = phi ** 2
        a.plot(z, v / v.max(), "o-", ms=3, color=col, label=f"{name}: FWHM {fwhm(z, v):.0f} mm")
    a.set_title(f"{MAP_FOLDER}: variance in each frequency band,\nnormalised to its peak")
    a.set_xlabel("position along beam (mm)")
    a.set_ylabel("normalised variance")
    a.legend(fontsize=8, frameon=False)

    a = ax[1, 0]
    m = read_csv("2025-06-27_T2_frequency_map.csv")
    zs = np.array(sorted({float(r["position_mm"]) for r in m}))
    fs = np.array(sorted({float(r["f_Hz"]) for r in m}))
    grid = np.full((len(fs), len(zs)), np.nan)
    for r in m:
        grid[np.searchsorted(fs, float(r["f_Hz"])), np.searchsorted(zs, float(r["position_mm"]))] = float(r["psd_rad2_per_Hz"])
    grid = grid / np.nanmax(grid, axis=1, keepdims=True)  # normalise each frequency to its own peak
    pc = a.pcolormesh(zs, fs / 1e3, np.clip(grid, 1e-2, 1), norm=LogNorm(1e-2, 1), cmap="magma", shading="nearest")
    a.set_yscale("log")
    a.set_ylim(0.1, 40)  # above ~40 kHz only electronic lines and noise remain
    fig.colorbar(pc, ax=a, label="PSD / peak PSD at that frequency")
    a.set_title(f"{MAP_FOLDER}: where along the beam each frequency is sensed")
    a.set_xlabel("position along beam (mm)")
    a.set_ylabel("frequency (kHz)")

    a = ax[1, 1]
    z, phi = series("Vertical - T1", "phase_rms_rad")
    _, shift = series("Vertical - T1", "mean_phase_shift_rad")
    a.plot(z, phi, "o-", ms=3.5, color="#c42b23", label="phase RMS (fluctuations)")
    a.set_xlabel("vertical gun position (mm)")
    a.set_ylabel("phase RMS (rad)", color="#c42b23")
    b = a.twinx()
    b.plot(z, shift - shift[0], "s--", ms=3, color="#1d6a86", label="mean phase shift")
    b.set_ylabel("mean phase shift from gun-edge value (rad)", color="#1d6a86")
    a.set_title("Across the beam (vertical traverse near the focus):\ntwo shear layers, mean shifts of opposite sign")

    for a in ax.ravel():
        a.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIG, dpi=150)
    print("wrote", FIG)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--raw", help='the "27 June 2025" folder (with 27th june.xlsx and the four traverse folders)')
    args = ap.parse_args()
    if args.raw:
        from_raw(args.raw)
    plot()
