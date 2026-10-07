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
| `ArcCase` | Arc welding: surface heat flux on a substrate between two gas regions, with power ramp-down after 0.25 s and extinguishing at 0.35 s; shows Marangoni-driven surface flow and weld-pool penetration | 1 × 60 × 120 |
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
