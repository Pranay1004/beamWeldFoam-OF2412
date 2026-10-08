#!/usr/bin/env python3
"""Post-process a beamWeldFoam ArcCase run.

Reads the function-object time series written by system/controlDict
(temperatureMax, meltVolume, metalVolume, temperatureProbes) and the
field snapshots of the written time directories, then produces

  * <out>_results.png    four-panel figure (time series, probe
                         histories, temperature map, liquid-fraction map)
  * <out>_timeseries.csv combined function-object time series
  * a printed summary    peak temperature, melt volume, penetration
                         depth, pool length, metal-volume drift

Usage:
    python3 postProcess.py [caseDir] [--out PREFIX]

Requires numpy and matplotlib.  The cell centres are taken from the
Cx/Cy/Cz fields written by the cellCentres function object; if those
are absent the script falls back to rebuilding them from
constant/polyMesh.
"""

import argparse
import re
import sys
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# --------------------------------------------------------------------------
# OpenFOAM file helpers
# --------------------------------------------------------------------------

def read_scalar_field(path):
    """Read an OpenFOAM ASCII volScalarField internalField as a 1-D array."""
    text = Path(path).read_text(errors="replace")
    m = re.search(
        r"internalField\s+nonuniform\s+List<scalar>\s*\(?\s*(\d+)\s*\)?\s*\((.*?)\)\s*;",
        text, re.S)
    if m:
        return np.array(m.group(2).split(), dtype=float)
    m = re.search(
        r"internalField\s+uniform\s+([-+0-9.eE]+)", text)
    if m:
        return np.array([float(m.group(1))])
    raise ValueError(f"cannot parse internalField in {path}")


def read_fo_dat(path):
    """Read a volFieldValue .dat file -> (times, values[ntime, ncol])."""
    rows = []
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        rows.append([float(x) for x in line.split()])
    if not rows:
        return np.zeros(0), np.zeros((0, 0))
    arr = np.array(rows)
    return arr[:, 0], arr[:, 1:]


def read_probes(path):
    """Read a probes .dat file -> (locations[5,3], times, values[ntime, 5])."""
    locs, rows = [], []
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if line.startswith("# Probe"):
            m = re.search(r"\(\s*([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s*\)", line)
            if m:
                locs.append([float(g) for g in m.groups()])
        elif line and not line.startswith("#"):
            rows.append([float(x) for x in line.split()])
    if not rows:
        return np.zeros((0, 3)), np.zeros(0), np.zeros((0, 0))
    arr = np.array(rows)
    return np.array(locs), arr[:, 0], arr[:, 1:]


def time_dirs(case):
    """Numeric time directories of a case, sorted ascending (0 included)."""
    out = []
    for p in case.iterdir():
        if p.is_dir() and re.fullmatch(r"\d+(\.\d+)?([eE][-+]?\d+)?", p.name):
            try:
                out.append((float(p.name), p))
            except ValueError:
                pass
    return sorted(out)


# --------------------------------------------------------------------------
# cell centres and volumes
# --------------------------------------------------------------------------

def read_polymesh(case):
    """points, faces, owner, neighbour arrays from constant/polyMesh."""
    pm = case / "constant" / "polyMesh"

    def body_of(text):
        text = re.sub(r"FoamFile\s*\{.*?\}", "", text, flags=re.S)
        m = re.search(r"\d+\s*\((.*)\)", text, re.S)
        return m.group(1)

    pts_body = body_of((pm / "points").read_text(errors="replace"))
    points = np.array(
        [p.split() for p in re.findall(r"\(([^)]*)\)", pts_body)], dtype=float)

    owner = np.array(body_of((pm / "owner").read_text(errors="replace")).split(),
                     dtype=np.int64)
    neighbour = np.array(
        body_of((pm / "neighbour").read_text(errors="replace")).split(),
        dtype=np.int64)

    faces_body = body_of((pm / "faces").read_text(errors="replace"))
    faces = [(int(n), [int(i) for i in ids.split()])
             for n, ids in re.findall(r"(\d+)\s*\(([^)]*)\)", faces_body)]
    return points, faces, owner, neighbour


def cell_data(case, times):
    """Cell centres (xc, yc, zc) and volumes.

    Centres come from the Cx/Cy/Cz fields written by the cellCentres
    function object when present, otherwise from constant/polyMesh
    (vertex mean, exact for hexahedra).  Volumes always come from
    constant/polyMesh via the divergence theorem, which is also exact
    for directions with a single cell.
    """
    xc = yc = zc = None
    for _, td in reversed(times):
        if all((td / f"C{c}").exists() for c in "xyz"):
            xc = read_scalar_field(td / "Cx")
            yc = read_scalar_field(td / "Cy")
            zc = read_scalar_field(td / "Cz")
            break

    points, faces, owner, neighbour = read_polymesh(case)
    ncells = int(max(owner.max(), neighbour.max())) + 1

    if xc is None:  # fallback: rebuild centres from polyMesh
        cell_pts = [[] for _ in range(ncells)]
        for fi, (_, pts) in enumerate(faces):
            cell_pts[owner[fi]].extend(pts)
            if fi < len(neighbour):
                cell_pts[neighbour[fi]].extend(pts)
        centres = np.array([points[np.unique(cell_pts[c])].mean(axis=0)
                            for c in range(ncells)])
        xc, yc, zc = centres[:, 0], centres[:, 1], centres[:, 2]

    vols = np.zeros(ncells)
    for fi, (_, ids) in enumerate(faces):
        P = points[ids]
        fc = P.mean(axis=0)
        sf = np.zeros(3)
        for i in range(1, len(P) - 1):
            sf += np.cross(P[i] - P[0], P[i + 1] - P[0])
        contrib = np.dot(fc, sf) / 6.0   # (1/3) * (Cf . (Sf/2))
        vols[owner[fi]] += contrib
        if fi < len(neighbour):          # internal face: outward for neighbour
            vols[neighbour[fi]] -= contrib
    return xc, yc, zc, vols


def grid_index(coords, tol=1e-9):
    """Map unstructured coordinates of a rectilinear grid to (i, n, pos)."""
    pos = np.unique(np.round(coords, 12))
    idx = np.abs(coords[:, None] - pos[None, :]).argmin(axis=1)
    return idx, len(pos), pos


# --------------------------------------------------------------------------
# case parameters
# --------------------------------------------------------------------------

def case_parameters(case):
    """Extract process/material numbers needed for the report."""
    p = {"v": None, "lg": None, "tsol": None, "tliq": None}
    fvs = (case / "system" / "fvSolution").read_text(errors="replace")
    m = re.search(r"HS_velocity\s+([-+0-9.eE]+)", fvs)
    if m:
        p["v"] = float(m.group(1))
    m = re.search(r"HS_lg\s+([-+0-9.eE]+)", fvs)
    if m:
        p["lg"] = float(m.group(1))
    m = re.search(r"HS_Q\s+([-+0-9.eE]+)", fvs)
    if m:
        p["Q"] = float(m.group(1))
    tp = (case / "constant" / "transportProperties").read_text(errors="replace")
    m = re.search(r"Tsolidus\s+([-+0-9.eE]+)", tp)
    if m:
        p["tsol"] = float(m.group(1))
    m = re.search(r"Tliquidus\s+([-+0-9.eE]+)", tp)
    if m:
        p["tliq"] = float(m.group(1))
    return p


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("case", nargs="?", default=str(Path(__file__).resolve().parent),
                    help="case directory (default: directory of this script)")
    ap.add_argument("--out", default="postProcess",
                    help="output prefix for .png/.csv (default: postProcess)")
    ap.add_argument("--time", type=float, default=None,
                    help="snapshot time for the field maps (default: latest)")
    args = ap.parse_args()

    case = Path(args.case).resolve()
    pp = case / "postProcessing"
    if not pp.is_dir():
        sys.exit(f"no postProcessing directory in {case} - run the case first")

    params = case_parameters(case)
    tsol, tliq = params["tsol"] or 1673.0, params["tliq"] or 1723.0
    v, lg = params["v"] or 0.0, params["lg"] or 0.0

    # ---- function-object time series (merged over restarts) ---------------
    # On resume (startFrom latestTime) each function object writes into a new
    # postProcessing/<name>/<startTime>/ directory, so all of them are read
    # and concatenated.
    def merged_fo(name, fname):
        rows = []
        for d in sorted((pp / name).glob("*")):
            # nb: on restarts OpenFOAM may append _<startTime> to the file
            # name (volFieldValue_0.04.dat), so read every .dat file
            for f in sorted(d.glob("*.dat")):
                t, y = read_fo_dat(f)
                if y.size:
                    rows.append(np.column_stack([t, y]))
        if not rows:
            return None, None
        arr = np.vstack(rows)
        arr = arr[np.argsort(arr[:, 0])]
        keep = np.ones(len(arr), bool)   # drop duplicate stamps, keep last
        keep[:-1] = arr[1:, 0] != arr[:-1, 0]
        arr = arr[keep]
        return arr[:, 0], arr[:, 1:]

    series, labels = [], []
    for name, col in (("temperatureMax", "T max [K]"),
                      ("meltVolume", "melt volume [m3]"),
                      ("metalVolume", "metal volume [m3]")):
        t, y = merged_fo(name, "volFieldValue.dat")
        if t is not None and y.size:
            series.append((name, t, y[:, 0]))
            labels.append(col)

    def merged_probes(field):
        locs, all_t, all_v = None, [], []
        for d in sorted((pp / "temperatureProbes").glob("*")):
            f = d / field
            if f.is_file():
                locs_, t, v = read_probes(f)
                if locs is None and len(locs_):
                    locs = locs_
                if t.size:
                    all_t.append(t)
                    all_v.append(v)
        if not all_t:
            return None, np.zeros(0), np.zeros((0, 0))
        t = np.concatenate(all_t)
        v = np.vstack(all_v)
        order = np.argsort(t)
        return locs, t[order], v[order]

    probe_times, probe_vals, probe_locs = {}, {}, None
    probe_locs, probe_times["T"], probe_vals["T"] = merged_probes("Temperature")
    _, probe_times["eps"], probe_vals["eps"] = merged_probes("epsilon1")

    times = time_dirs(case)
    times = [td for td in times
             if (td[1] / "Temperature").exists()
             and (td[1] / "alpha.phase1").exists()]
    if not times:
        sys.exit(f"no usable (field) time directories found in {case}")
    if args.time is not None:
        earlier = [td for td in times if td[0] <= args.time + 1e-12]
        if not earlier:
            sys.exit(f"no snapshot at or before t = {args.time:g}")
        t_end, td_end = earlier[-1]
    else:
        t_end, td_end = times[-1]

    # ---- final-time fields ----------------------------------------------
    need = ["Temperature", "alpha.phase1"]
    missing = [n for n in need if not (td_end / n).exists()]
    if missing:
        sys.exit(f"missing field(s) {missing} in {td_end}")

    T = read_scalar_field(td_end / "Temperature")
    alpha = read_scalar_field(td_end / "alpha.phase1")
    eps_path = td_end / "epsilon1"
    eps = read_scalar_field(eps_path) if eps_path.exists() else None
    xc, yc, zc, vols = cell_data(case, times)

    metal = alpha > 0.5
    if eps is not None:
        frac = np.clip(eps, 0, 1)          # solver liquid fraction
    else:                                   # fall back to temperature
        frac = np.clip((T - tsol) / max(tliq - tsol, 1e-9), 0, 1)
    frac = np.where(metal, frac, np.nan)    # gas cells -> NaN (epsilon1 == 1 there)
    liquid = np.nan_to_num(frac) > 0.5
    liquid_vol = float(np.nansum(frac * vols))

    # rectilinear grid (single cell in x for this case)
    iy, ny, yp = grid_index(yc)
    iz, nz, zp = grid_index(zc)
    Tmap = np.full((ny, nz), np.nan)
    Fmap = np.full((ny, nz), np.nan)
    Mmap = np.zeros((ny, nz))
    Tmap[iy, iz] = T
    Fmap[iy, iz] = frac
    Mmap[iy, iz] = liquid

    # plate geometry: arc-face datum = lowest metal layer that is part of
    # the contiguous body (ignores isolated spatter above the pool)
    y_metal = yc[metal]
    if y_metal.size:
        edges = np.arange(0.0, y_metal.max() + 0.0002, 0.0002)
        counts, _ = np.histogram(y_metal, bins=edges)
        body = np.where(counts > 0.02 * metal.sum())[0]
        plate_top = edges[body[0]] if len(body) else y_metal.min()
        plate_bot = y_metal.max()
    else:
        plate_top, plate_bot = (0.0, 0.0)

    melt_cells = liquid & (yc >= plate_top - 0.0001)  # pool, not spatter
    if melt_cells.any():
        depth = float(yc[melt_cells].max() - plate_top)  # from arc face downwards
        zlo, zhi = zc[melt_cells].min(), zc[melt_cells].max()
        pool_len = float(zhi - zlo)
        ymax = float(T.max())
    else:
        depth = pool_len = 0.0
        ymax = float(T.max())

    # ---- summary numbers -------------------------------------------------
    tmax_t, tmax_v = (series[0][1], series[0][2]) if series else (np.zeros(0), np.zeros(0))
    melt_t, melt_v = (series[1][1], series[1][2]) if len(series) > 1 else (np.zeros(0), np.zeros(0))
    metal_t, metal_v = (series[2][1], series[2][2]) if len(series) > 2 else (np.zeros(0), np.zeros(0))

    print("=" * 66)
    print(f" beamWeldFoam ArcCase post-processing   ({case})")
    print("=" * 66)
    print(f" last written time          : {t_end:g} s")
    if params.get("Q"):
        print(f" arc power HS_Q             : {params['Q']:.2f} W")
    print(f" peak temperature           : {ymax:.1f} K "
          f"({ymax - 273.15:.0f} degC)")
    print(f" melt volume (epsilon1)     : {liquid_vol * 1e9:.4f} mm3")
    if melt_v.size:
        print(f"   FO meltVolume last value : {melt_v[-1] * 1e9:.4f} mm3")
    print(f" penetration depth (from arc face y = {plate_top*1e3:.3f} mm) : "
          f"{depth*1e3:.3f} mm")
    print(f" pool length along z        : {pool_len*1e3:.3f} mm")
    if metal_v.size:
        drift = 100 * (metal_v.max() - metal_v.min()) / metal_v.mean()
        print(f" metal volume               : {metal_v[-1] * 1e9:.4f} mm3 "
              f"(drift {drift:.2e} %)")
    if (probe_locs is not None and "T" in probe_vals
            and probe_vals["T"] is not None and probe_vals["T"].size):
        print(" probe temperatures at last sample [K]:")
        for i, loc in enumerate(probe_locs):
            print(f"   P{i} (y={loc[1]*1e3:6.3f} mm, z={loc[2]*1e3:7.3f} mm) : "
                  f"{probe_vals['T'][-1, i]:8.2f}")
    print("=" * 66)

    # ---- CSV -------------------------------------------------------------
    csv = Path(args.out).with_suffix(".timeseries.csv")
    with open(csv, "w") as fh:
        fh.write("time," + ",".join(labels) + "\n")
        if series:
            grid = np.unique(np.concatenate([s[1] for s in series]))
            cols = []
            for _, tt, vv in series:
                cols.append(np.interp(grid, tt, vv))
            for k, t in enumerate(grid):
                fh.write(f"{t:g}," + ",".join(f"{c[k]:.10e}" for c in cols) + "\n")
    print(f"wrote {csv}")

    # ---- figure ----------------------------------------------------------
    fig, ax = plt.subplots(2, 2, figsize=(14, 9.5))

    a = ax[0, 0]
    if tmax_t.size:
        a.plot(tmax_t * 1e3, tmax_v, "r-", label="max T")
    a.set_xlabel("time [ms]")
    a.set_ylabel("peak temperature [K]", color="r")
    a.tick_params(axis="y", colors="r")
    a.grid(alpha=0.3)
    a2 = a.twinx()
    if melt_t.size:
        a2.plot(melt_t * 1e3, melt_v * 1e9, "b-", label="melt volume")
    a2.set_ylabel("melt volume [mm$^3$]", color="b")
    a2.tick_params(axis="y", colors="b")
    a.set_title("(a) peak temperature and melt volume")

    a = ax[0, 1]
    if (probe_locs is not None and "T" in probe_vals
            and probe_vals["T"] is not None and probe_vals["T"].size):
        names = ["arc face", "mid", "back face", "gas above", "arc start"]
        for i in range(probe_vals["T"].shape[1]):
            lbl = names[i] if i < len(names) else f"P{i}"
            a.plot(probe_times["T"] * 1e3, probe_vals["T"][:, i],
                   label=f"P{i} {lbl} (y={probe_locs[i,1]*1e3:.2f} mm)")
        a.legend(fontsize=8)
    a.set_xlabel("time [ms]")
    a.set_ylabel("temperature [K]")
    a.set_title("(b) history probes")
    a.grid(alpha=0.3)

    a = ax[1, 0]
    zm, ym = np.meshgrid(zp * 1e3, yp * 1e3)
    im = a.pcolormesh(zm, ym, Tmap, shading="nearest", cmap="inferno")
    fig.colorbar(im, ax=a, label="T [K]")
    if liquid.any():
        a.contour(zm, ym, Mmap, levels=[0.5], colors="cyan", linewidths=1.0)
    if v and lg is not None:
        a.axvline((v * t_end + lg) * 1e3, color="lime", ls="--", lw=1,
                  label=f"torch z={v*t_end+lg:.4f} m")
        a.legend(loc="upper right", fontsize=8)
    a.set_xlabel("z [mm]")
    a.set_ylabel("y [mm]  (arc side up)")
    a.invert_yaxis()          # y+ is gravity direction -> depth downwards
    a.set_title(f"(c) temperature at t = {t_end:g} s")

    a = ax[1, 1]
    im = a.pcolormesh(zm, ym, Fmap, shading="nearest", cmap="viridis",
                      vmin=0, vmax=1)
    fig.colorbar(im, ax=a, label="liquid fraction $\\epsilon_1$")
    a.contour(zm, ym, Mmap, levels=[0.5], colors="red", linewidths=1.0)
    a.set_xlabel("z [mm]")
    a.set_ylabel("y [mm]")
    a.invert_yaxis()
    a.set_title("(d) liquid fraction in the metal "
                f"(depth {depth*1e3:.2f} mm, length {pool_len*1e3:.2f} mm)")

    fig.suptitle("beamWeldFoam ArcCase - GMAW SS308L "
                 f"(HS_Q = {params.get('Q', float('nan')):.0f} W, "
                 f"v = {v*1e3*60:.0f} mm/min)", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    out = Path(args.out).with_suffix(".png")
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
