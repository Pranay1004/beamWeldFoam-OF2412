# beamWeldFoam (OpenFOAM v2412)

Open-source volume-of-fluid (VOF) solver for high-energy-density advanced
manufacturing processes — laser welding, electron-beam welding, arc welding and
additive manufacturing.

This repository is a port of
[tomflint22/beamWeldFoam](https://github.com/tomflint22/beamWeldFoam)
(upstream developed against OpenFOAM 6) to **OpenFOAM v2412**
(the openfoam.com release, `v2412`).

---

## Overview

In this implementation the metallic substrate and the shielding gas phase are
treated as incompressible. The solver fully captures the fusion/melting state
transition of the metallic substrate. For the vapourisation of the substrate the
explicit volumetric dilation due to the vapourisation state transition is
neglected; instead a phenomenological recoil-pressure term captures the
contribution to the momentum and energy fields due to vapourisation events.

beamWeldFoam also captures:

* surface tension effects, and their temperature dependence (Marangoni),
* latent heat effects due to melting/fusion (and vapourisation),
* buoyancy due to thermal expansion of the phases (Boussinesq approximation),
* momentum damping due to solidification,
* a representative heat source for an incident laser/electron beam, which can
  also be configured to represent arc-welding processes.

The solver approach is based on the adiabatic two-phase `interFoam` code
developed by [OpenCFD Ltd.](https://www.openfoam.com/).

---

## Requirements

| Component | Notes |
|---|---|
| OpenFOAM **v2412** | openfoam.com release, e.g. the `openfoam2412` apt package, installed under `/usr/lib/openfoam/openfoam2412` |
| C++17 compiler | `g++` 9 or newer recommended |
| `make`, `git` | standard build tools |
| OpenMPI | only needed for parallel runs (`decomposePar` / `mpirun`) |

The OpenFOAM environment must be sourced in every shell you use:

```bash
source /usr/lib/openfoam/openfoam2412/etc/bashrc
```

(Adjust the path if OpenFOAM v2412 is installed elsewhere; add the line to
`~/.bashrc` to make it permanent.)

---

## Building

```bash
git clone https://github.com/Pranay1004/beamWeldFoam-OF2412.git
cd beamWeldFoam-OF2412
./Allwmake
```

`Allwmake` will:

1. stop with a clear message if OpenFOAM has not been sourced,
2. link-test the compiler found first on `PATH` against the OpenFOAM libraries,
   and automatically fall back to a system compiler (e.g. `g++-9`) if the active
   toolchain cannot link them — this happens with conda toolchains, which ship
   an old glibc sysroot that cannot resolve OpenFOAM's symbols,
3. compile the solver with `wmake` into `$FOAM_USER_APPBIN`.

Clean the build with:

```bash
./Allwclean
```

Check the result:

```bash
beamWeldFoam -help
```

### Optional: install for all users / any directory

The build lands in your user directory
(`$FOAM_USER_APPBIN`, which is on `PATH`). To also make it available from the
main OpenFOAM installation for every user:

```bash
sudo cp $FOAM_USER_APPBIN/beamWeldFoam $FOAM_APPBIN/
```

`$FOAM_APPBIN` is `/usr/lib/openfoam/openfoam2412/platforms/linux64GccDPInt32Opt/bin`,
where all standard OpenFOAM solvers live, so afterwards `beamWeldFoam` can be
called from any case directory without any extra setup.

---

## Running a case

Every tutorial keeps its initial fields in `initial/` rather than `0/`.
From inside a case directory:

```bash
cp -r initial 0        # create the initial time directory
blockMesh              # generate the mesh
setFields              # initialise the phase field
beamWeldFoam           # run in serial
```

Parallel (6 ranks in this example):

```bash
cp -r initial 0
blockMesh
setFields
decomposePar
mpirun -np 6 beamWeldFoam -parallel
reconstructPar
```

Redirect any of these commands if you want to keep a log, e.g.
`blockMesh > log.blockMesh 2>&1`.

Remove previous results before re-running (`rm -rf 0* [1-9]* processor* log`),
or use the provided `Allclean` scripts.

> **Note — MPI and conda environments.** If `mpirun` segfaults immediately at
> startup (before any OpenFOAM output), a conda environment is usually
> injecting its own OpenMPI libraries through `LD_LIBRARY_PATH` and clashing
> with the system OpenMPI. Launch with a cleaned library path instead:
>
> ```bash
> LD_LIBRARY_PATH=$(echo "$LD_LIBRARY_PATH" | tr ':' '\n' \
>     | grep -viE "miniforge|conda|envs" | paste -sd:) \
>   mpirun -np 6 beamWeldFoam -parallel
> ```

`PowderBed2D` and `PowderBed3D` ship ready-made run scripts:

```bash
cd tutorials/PowderBed2D
./Allrun              # serial
./Allrun parallel     # parallel (decomposePar / runParallel / reconstructPar)
```

---

## Tutorial cases

| Case | Description | Mesh |
|---|---|---|
| `Test` | Minimal smoke-test case — start here to verify the build | 1 × 60 × 60 (2D) |
| `GalluimCase` | Gallium melting in a cavity: melting, buoyancy-driven flow and solidification — classic validation case [1] | 1 × 140 × 100 |
| `GalluimCase3D` | 3D version of the gallium melting case | 120 × 280 × 200 |
| `SenDavies` | Marangoni (thermocapillary) flow in a partially filled cavity with a flat interface and imposed temperature gradient; analytical steady-state solution available [2] | 1 × 240 × 480 |
| `ArcCase` | Arc welding of SS308L built around a specific GMAW process (see table below): Marangoni-driven surface flow, solidification and weld-pool penetration over an 8 mm pass | 1 × 60 × 120 |
| `SingleParticleMelting2D` | Melting of a single particle | 1 × 60 × 60 |
| `PowderBed2D` | Laser powder-bed (additive manufacturing) case with `LaserProperties` heat-source dictionary | 1 × 240 × 480 |
| `PowderBed3D` | 3D powder-bed case | 60 × 60 × 60 |
| `EB_3D` | 3D high-energy-beam welding case | 60 × 60 × 240 |
| `Tancase` | Welding benchmark case | 1 × 50 × 150 |

The gallium melting and Sen & Davies cases reproduce experimental/analytical
data from the literature and serve as validation of the implementation; the
welding cases exercise the heat source, Marangoni flow, solidification damping
and phase-change treatment (a Ti6Al4V laser butt-weld validation is reported in
[3]).

### ArcCase process parameters

`ArcCase` is built around a single GMAW pass on **SS308L**:

| Parameter | Value | Set in |
|---|---|---|
| Welding current | 151 A | folded into `HS_Q` |
| Arc voltage | 19.7 V | folded into `HS_Q` |
| Arc efficiency | 0.8 | folded into `HS_Q` |
| Arc power into the workpiece | 2379.76 W (= 151 × 19.7 × 0.8) | `system/fvSolution` → `MELTING → HS_Q` |
| Travel speed | 500 mm/min (8.3333e-3 m/s) | `MELTING → HS_velocity` |
| Torch start / direction | z = −4 mm, travelling to +z (8 mm domain, 0.96 s to cross) | `MELTING → HS_lg` |
| Arc spot radius | 2 mm | `MELTING → HS_a` |
| Power profile | full power for the whole pass | `MELTING → Q_ramp_rate 0`, `tshift 0` |
| Pass duration | 1.0 s | `system/controlDict → endTime` |
| Wire feed speed | 6 m/min — **not modelled** | — |
| Material | SS308L, see below | `constant/transportProperties` |

Material (`phase1`): ρ = 7900 kg/m³, c_p = 800/500 J/kg·K (liquid/solid),
κ = 30/20 W/m·K, T_sol = 1673 K, T_liq = 1723 K (1400/1450 °C),
L_f = 270 kJ/kg, ν = 7.5e-7 m²/s, β = 1.6e-5 /K; interface/vapour:
σ = 1.85 N/m, dσ/dT = −0.4e-3 N/m·K, T_vap = 3135 K, M = 0.0558 kg/mol,
L_vap = 6.1e6 J/kg.

Notes:

* **Wire feed (6 m/min).** beamWeldFoam's source term adds heat only — there
  is no filler-wire mass source, so deposited filler metal is not added to the
  pool, and the energy that would melt the wire is not removed either. For a
  1.0 mm wire at 6 m/min that sink is ≈ 0.62 g/s × ≈ 1.41 MJ/kg ≈ **875 W
  (~37 % of the arc power)**, so the modelled net heat input — and hence
  penetration — is high relative to the real process unless it is considered
  part of the 0.8 arc efficiency.
* **Cost.** With `maxCo`/`maxAlphaCo` = 0.5/0.25 and `maxDeltaT` = 5e-6 the step
  pins at Δt = 5e-6 s, so the full 1.0 s pass is ≈ 200k steps (≈ 7 h on
  20 ranks, 2 outer correctors).
* **Stability.** The case as first configured blew up at ≈ 14–18 ms
  (spurious velocities → Δt collapse → continuity errors ~1e16), twice, just
  as the pool approached T_vap. Three changes, all documented at the entries
  themselves, stabilised it through 40 ms and counting: gas density
  ρ₂ = 0.1 → 1.0 kg/m³ (the 79000:1 metal/gas ratio destabilised the
  pressure–velocity coupling), gas viscosity ν₂ = 1.48e-5 → 1e-3 m²/s (the
  hot gas plume, not the pool, was limiting Δt — pool velocities stay
  sub-m/s), upwind convection for `U`, `nOuterCorrectors` 1 → 2, and
  `damperSwitch true` so interface forces/sinks partition by density
  (without it the explicit evaporative sink drives low-heat-capacity gas
  cells to negative T and the run dies in a floating-point exception).
  A second crash mode (SIGFPE in `pow(2πMmRT, 0.5)` when an overshoot cell
  hits T ≤ 0) is guarded by a minimal solver patch — `TEqn.H`/`UEqn.H`
  evaluate the evaporation/recoil exponentials with
  `Tpos = max(Temperature, SMALL)`, identical to before for all T > 0
  (see the `Tpos` comments in the sources).
* **Mass loss.** The pool sheds metal as spatter through the top patch from
  ≈ 0.05 s on (top-patch flux ≈ volume drift rate, so outflow — not
  numerical clipping): Marangoni shear tears flotsam off the surface, the
  single-cell beam locks onto the lowest dense droplets instead of the pool
  (starving it, T_max sags) until they exit. Fixes so far:
  `HS_deposition_cutoff` 0.5 → 0.99 (conduction mode — rays ignore smeared
  cells), `cAlpha` back to 1.0 (0.5 smeared the interface and leaked mass
  through MULES clipping), `epsilonRelaxation` 0.95 → 0.5, `scAlpha`
  0.5 → 0, `recoilScale` 0.2 (calms recoil digging; Qv cooling untouched),
  plus a `-1` guard for empty rays (`DivergingOscillatingGaussian.H`).
  Status at 0.24 s: near-full-penetration pool (depth ≈ 2 mm, length
  ≈ 6.5 mm), T_max plateau ≈ 3100–3300 K, sub-m/s pool flow, metal drift
  ≈ 5% and accelerating — the bleed persists through every
  numerics-only fix, pointing at the driving physics itself (Marangoni
  strength / power density), see below.
* **Open point — Marangoni strength.** The `dsigmadT = -0.4e-3` N/m/K used
  here is the pure-iron end of the range (my choice; upstream used -0.1e-4).
  For SS308L the effective value depends on S/O content and can be much
  smaller (even positive). The pool violence — and hence the spatter rate —
  scales directly with it, so the full 1.0 s pass needs either a
  literature-backed value for this wire/plate combination or acceptance of
  a spattery pool with documented mass loss.

### ArcCase post-processing

`system/controlDict` defines function objects that write during the run:

| Function object | Output | Purpose |
|---|---|---|
| `temperatureMax` | `postProcessing/temperatureMax/0/volFieldValue.dat` | peak temperature history |
| `meltVolume` | `postProcessing/meltVolume/0/volFieldValue.dat` | melted metal volume, ε₁ weighted by α₁ (raw `volIntegrate(epsilon1)` would include the gas, where ε₁ ≡ 1) |
| `metalVolume` | `postProcessing/metalVolume/0/volFieldValue.dat` | total metal volume — must stay constant (no mass sources) |
| `temperatureProbes` | `postProcessing/temperatureProbes/0/{Temperature,epsilon1}` | five history probes: arc face, mid-thickness, back face, gas above the plate, arc start point |
| `cellCentres` | `Cx`, `Cy`, `Cz` in the last time directory | cell coordinates for the script below |

A single Python script turns these into a figure and a summary:

```bash
python3 tutorials/ArcCase/postProcess.py      # writes postProcess.png and postProcess.timeseries.csv
python3 tutorials/ArcCase/postProcess.py myCase --out results
```

It reads the function-object time series plus the field snapshots
(`Temperature`, `alpha.phase1`, `epsilon1`) of every written time directory
(cell centres are rebuilt from `constant/polyMesh` if `Cx`/`Cy`/`Cz` are not
available), prints peak temperature, melt volume, penetration depth and pool
length, and saves a four-panel figure: peak-T/melt-volume time series, probe
histories, the temperature map, and the liquid-fraction map with the melt-pool
outline. Requires `numpy` and `matplotlib`.

For parallel runs the snapshots live in `processor*/` until reconstructed:

```bash
reconstructPar                # all times
reconstructPar -time 0.01     # a single snapshot
```

---

## Algorithm

Initially the solver loads the mesh, reads in fields and boundary conditions,
reads certain mesh information into arrays (for the heat source application),
and selects the turbulence model (if specified). The main solver loop is then
initiated. First the time step is dynamically modified to ensure numerical
stability. Next, the two-phase fluid mixture properties and turbulence
quantities are updated. The discretized phase-fraction equation is then solved
for a user-defined number of sub-time steps (typically 3) using the
multidimensional universal limiter with explicit solution solver
[MULES](https://openfoam.org/release/2-3-0/multiphase/), which performs
conservative solution of hyperbolic convective transport equations with defined
bounds (0 and 1 for α₁). Once the updated phase field is obtained, the program
enters the pressure–velocity loop, in which *p* and *U* are corrected in an
alternating fashion; *T* is also solved here so that the buoyancy predictions
are consistent with the *U* and *p* fields. The sequential correction of
pressure and velocity is the pressure implicit with splitting of operators
(PISO) algorithm. In OpenFOAM, PISO may be repeated for multiple iterations at
each time step; combined with SIMPLE-type outer iterations this is the merged
PIMPLE algorithm, run for a user-specified number of outer correctors.

The main solver loop iterates until program termination:

* **beamWeldFoam simulation algorithm summary**
  * Initialize simulation data and mesh
  * **WHILE** t < t_end **DO**
    1. Update Δt for stability
    2. Phase-equation sub-cycle
    3. Update interface location for heat-source application
    4. Update fluid properties
    5. PISO loop
       1. Form the *U* equation
       2. Energy transport loop
          1. Solve the *T* equation
          2. Update the fluid-fraction field
          3. Re-evaluate source terms due to latent heat
       3. PISO
          1. Obtain and correct face fluxes
          2. Solve the *p*-Poisson equation
          3. Correct *U*
    6. Write fields

---

## Repository layout

```
beamWeldFoam-OF2412/
├── Allwmake                     build the solver (with compiler fallback)
├── Allwclean                    remove build products
├── applications/
│   └── solvers/beamWeldFoam/    solver sources (beamWeldFoam.C + *.H, VoF/)
│       └── Make/                wmake files/options
└── tutorials/                   example cases (initial/ + system/ + constant/)
```

---

## Changes vs. upstream (OpenFOAM 6 → v2412)

* `dimensionedScalar` lookups, `CorrectPhi.H` inclusion and the VoF
  time-step headers updated for the v2412 API (lowercase `small`/`great` are
  only available under `COMPAT_OPENFOAM_ORG`, so `SMALL`/`GREAT` are used).
* `Allwmake` now link-tests the active compiler and falls back to a system
  compiler when a toolchain (e.g. conda) cannot link the OpenFOAM libraries.
* Mismatched braces in the `relaxationFactors` dictionary of four tutorials
  corrected (they previously aborted `setFields` on v2412).
* Stale solver name in `tutorials/ArcCase/Allrun` corrected.
* `ArcCase` re-parameterised around a 151 A / 19.7 V / η 0.8 / 500 mm/min
  GMAW pass on SS308L (previously an aluminium case with a stationary 780 W
  source and a ramp-down).

---

## License

OpenFOAM, and by extension the beamWeldFoam application, is licensed free and
open source only under the
[GNU General Public Licence version 3](https://www.gnu.org/licenses/gpl-3.0.en.html).
One reason for OpenFOAM's popularity is that its users are granted the freedom
to modify and redistribute the software and have a right of continued free use,
within the terms of the GPL.

## Acknowledgements

The work was generously supported by the Engineering and Physical Sciences
Research Council (EPSRC) under the "Cobalt-free Hard-facing for Reactor
Systems" grant EP/T016728/1, and Science Foundation Ireland (SFI), co-funded
under European Regional Development Fund and by I-Form industry partners,
grant 16/RC/3872.

## Citing this work

If you use beamWeldFoam in your work, please cite:

> Thomas F. Flint, Gowthaman Parivendhan, Alojz Ivankovic, Michael C. Smith,
> Philip Cardiff,
> *beamWeldFoam: Numerical simulation of high energy density fusion and
> vapourisation-inducing processes*,
> SoftwareX, Volume 18, 2022, 101065, ISSN 2352-7110,
> https://doi.org/10.1016/j.softx.2022.101065

## References

1. Kay Wittig and Petr A. Nikrityuk, *IOP Conf. Ser.: Mater. Sci. Eng.* **27**
   012054 (2012).
2. Sen, A., & Davis, S. (1982). Steady thermocapillary flows in two-dimensional
   slots. *Journal of Fluid Mechanics*, 121, 163–186.
   doi:10.1017/S0022112082001840
3. S. L. Campanelli, G. Casalino, M. Mortello, A. Angelastro, A. D. Ludovico,
   *Microstructural Characteristics and Mechanical Properties of Ti6Al4V Alloy
   Fiber Laser Welds*.
