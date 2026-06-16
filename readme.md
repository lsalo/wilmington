# Wilmington flow–geomechanics

Multiphase-flow and geomechanics modeling of the Wilmington Oil Field (Los Angeles Basin, California): reservoir flow in MRST, geomechanics in Abaqus, and Julia post-processing for subsidence, stress state, and fault stability.

[![DOI](https://zenodo.org/badge/719353093.svg)](https://doi.org/10.5281/zenodo.20723535)

## Citation

> Saló-Salgado, L., Silva, J. A., Plesch, A., Wolfe, F. D., Shaw, J. H. & Juanes, R. (2026), "Stress state, subsidence, and faulting in the Wilmington Oil Field, California: a multiphase flow–geomechanics modeling assessment (1936–2020)", in review.
>
> Preprint: arXiv:2605.16711 (https://doi.org/10.48550/arXiv.2605.16711).
> This record will be updated with the journal citation and DOI upon acceptance.

## Requirements
- **Don't have Julia?** Download the installer for your OS from
  https://julialang.org/downloads/ and run it — this installs Julia and the
  `juliaup` version manager. Then get the pinned version with `juliaup add 1.11`.

- **Julia 1.11.7** — pinned in `Manifest.toml`. Install with `juliaup add 1.11`.
- **External engines (full pipeline only):** MRST (open-source, MATLAB) and
  Abaqus (commercial). Neither is redistributed here; the model setup is
  specified in the manuscript Methods. The Julia post-processing and figure
  scripts run against the deposited dataset *without* these engines.

## Setup

1. Open this folder in VS Code; the Julia extension auto-activates the
   environment (confirm the status bar shows `wilmington`, not `v1.x`).
2. Install the locked dependencies — in the REPL: `] instantiate`.
3. Point at the data: download the dataset deposit, then create
   `config_local.jl` next to `config.jl`:

```julia
const LOCAL_DATA_DIR = "/path/to/wilmington/data"
```

   It is gitignored; without it, `config.jl` falls back to `../data`.
4. Run a script, e.g. `include("scripts/abaqus_results_fault_tractions.jl")`.
   Every script starts with `include(joinpath(@__DIR__, "..", "config.jl"))`,
   which activates the environment and defines all paths (`DATA_DIR`,
   `GRIDS_ABAQUS`, `OUT_DIR`, …) and the module loader `tool(:name)`.

## Layout

```
config.jl          path resolution + module loading (entry point for all scripts)
config_local.jl    machine-specific data path (gitignored)
Project.toml       dependencies;  Manifest.toml — locked versions, Julia 1.11.7
jl-tools/modules/  vendored analysis modules (provenance in jl-tools/SNAPSHOT.md)
scripts/           mrst_*                                   reservoir flow setup & results
                   abaqus_create_*                          build Abaqus inputs (incl. pressure from MRST)
                   abaqus_results_*                         post-process & plot
                   abaqus_python/                           run inside Abaqus, not Julia (see below)
                   mrst/TransmissibilityWithMultipliers     state fcn to modify face transmissibilities in MRST.  
output/ figures/   generated results — gitignored, created on run
```

Data (grids, raw and processed well rates, MRST/Abaqus simulation outputs) is
archived separately under CC-BY-4.0. Dataset DOI: [10.5281/zenodo.20612866](https://doi.org/10.5281/zenodo.20612866).
The SCEC Community Stress Model is third-party and cited at source, not included.
See the dataset README for the full layout.

## Abaqus Python scripts (`scripts/abaqus_python/`)

These run inside Abaqus (CAE Python terminal), not Julia. They extracted the
displacement/stress fields in the dataset's `simulation_output/abaqus/*/extracted/`
from the solver `.odb` files. They are **not required to reproduce the figures** —
the extracted fields are deposited directly — and are included only to document how
those fields were produced. To run one inside Abaqus:

```
execfile('filename.py')
```

## License

Code released under the MIT License (see `LICENSE`), **except**
`scripts/mrst/TransmissibilityWithMultipliers.m`, which extends the MATLAB Reservoir
Simulation Toolbox (MRST) and is therefore distributed under MRST's GNU GPL v3 license
(see the footer in that file). Data: CC-BY-4.0.