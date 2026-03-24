# openFoamCreateCase

A Python-based preprocessor for automating parametric OpenFOAM simulations. Given a set of inlet velocities and mesh configurations, it generates complete, ready-to-run OpenFOAM cases — handling mesh creation, fluid property lookup, turbulence parameter calculation, and boundary condition setup — then launches each simulation automatically.

---

## What It Does

Setting up many OpenFOAM cases by hand is tedious and error-prone. This tool automates the entire preprocessing pipeline:

1. **Reads parametric inputs** from two plain-text files — `inletData` (inlet velocities) and `meshData` (mesh cell counts) — and generates one case for every combination.
2. **Generates meshes** by copying a template (`headMesh`), substituting cell count values into `blockMeshDict`, and running the mesh script.
3. **Creates OpenFOAM cases** by copying a template (`headCase`), then automatically filling in:
   - Inlet velocity `U`
   - Turbulence parameters: turbulent kinetic energy `k`, dissipation rate `ε`, and specific dissipation rate `ω`
   - Fluid physical properties: density `ρ` and kinematic viscosity `ν` for air, fresh water, or sea water — interpolated from lookup tables at the configured temperature
   - Solver settings: `endTime`, `writeInterval`, `rhoInf`, number of parallel processors
4. **Runs each case** by invoking its `run.sh` bash script via `subprocess`.
5. **Skips completed cases** — if a case directory already contains `postProcessing/` output, it is left untouched. Incomplete cases (directory exists but no output) are cleaned and re-run.

---

## Repository Structure

```
openFoamCreateCase/
├── createCase_v1_3.py       # Main script — loops over all mesh/inlet combinations
├── openFoam_func_lib.py     # Core library: mesh creation, case setup, post-processing
├── initialConditions.py     # User-editable simulation parameters
├── inletData                # Inlet velocity sweep (one value per line)
├── meshData                 # Mesh configurations (cell count per line)
├── headCase/                # Template OpenFOAM case (must be provided by user)
├── headMesh/                # Template mesh case (must be provided by user)
└── meanProp/                # Fluid property lookup tables
    ├── air
    ├── fresh_water
    └── sea_water
```

---

## Configuration

All user-facing settings live in `initialConditions.py`:

| Parameter | Description | Default |
|---|---|---|
| `T` | Fluid temperature (°C) for property lookup | `20` |
| `nutnu` | Eddy viscosity ratio `νt/ν` for turbulence inlet BC | `10` |
| `Tu` | Turbulence intensity (%) | `0.5` |
| `endTime` | Simulation end time (s) | `20` |
| `writeInterval` | Output write interval (s) | `1` |
| `n_of_proc` | Number of MPI processes for `decomposeParDict` | `30` |

Inlet velocities are listed in `inletData` (one per line, after a header row) and mesh configurations in `meshData`.

---

## How the Physics Is Computed

Fluid properties (ρ and ν) are **interpolated from tabulated data** in `meanProp/` at the temperature `T` set in `initialConditions.py`. Three fluids are supported: air, fresh water, and sea water. Currently the tool uses fresh water and air properties; switching fluids requires a one-line change in `openFoam_func_lib.py`.

Turbulence inlet boundary conditions are derived from standard RANS relations:

```
u' = (Tu/100) · U
k  = 3/2 · u'²
ε  = Cμ · k² / (νt/ν · ν)
ω  = ε / (Cμ · k)         where Cμ = 0.09
```

These computed values are written directly into the `0/include/initialConditions` file of each case.

---

## Output

Each run produces a set of numbered case directories:

```
case1_U=0.1/
case1_U=0.2/
mesh1_1.0_c=120/
```

Post-processing utilities are included in `openFoam_func_lib.py`:

- **`prop_forces(path)`** — reads `postProcessing/forces/` output, strips parentheses from OpenFOAM vector notation, and returns time-averaged drag force components and moments.
- **`yPlus(path)`** — reads `postProcessing/yPlus/yPlus.dat` and returns min, max, and average y⁺ values from the last time step.

---

## Getting Started

### 1. Clone the repository

```bash
git clone https://github.com/<your-username>/openFoamCreateCase.git
cd openFoamCreateCase
```

### 2. Install Python dependencies

```bash
pip install numpy
```

### 3. Prepare your templates

The tool requires two OpenFOAM template directories that you provide:

- **`headCase/`** — a complete OpenFOAM case skeleton with placeholder tokens in the relevant files:
  - `0/include/initialConditions` → `U*`, `k*`, `omega*`
  - `constant/transportProperties` → `rho_air*`, `nu_air*`, `rho_water*`, `nu_water*`
  - `system/controlDict` → `endTime*`, `writeInterval*`, `rhoInf*`
  - `system/decomposeParDict` → `n_proc*`
  - `run.sh` → `loc`, `loc1`, `n_of_proc*`
- **`headMesh/`** — a blockMesh template with `cells*` in `system/blockMeshDict`

### 4. Set your simulation parameters

Edit `initialConditions.py` to set the temperature, turbulence settings, end time, and processor count.

Edit `inletData` to list the inlet velocities you want to sweep (one per line, keep the header row):

```
U, m/s
0.5 1
1.0 1
2.0 1
```

Edit `meshData` to list the mesh configurations:

```
Mesh number cells
1 120
1 240
```

### 5. Run the preprocessor

```bash
python createCase_v1_3.py
```

The script will print progress to the terminal, including start time, case names, and elapsed time per case. Each completed case will have its own directory (e.g. `case1_U=0.5/`) containing full OpenFOAM output.

### Re-running

- Cases with existing `postProcessing/` output are **skipped automatically**.
- Cases where the directory exists but output is missing are **deleted and re-run**.
- To force a full re-run, simply delete the case directories manually before running the script.

---

## Requirements

- Python 3.x
- NumPy
- OpenFOAM (accessible in the shell environment)
- A configured `headCase/` template with placeholder strings (e.g. `U*`, `k*`, `omega*`, `rho_water*`, `nu_water*`, `endTime*`, etc.)
- A configured `headMesh/` template with a `cells*` placeholder in `blockMeshDict`

---

## License

GNU General Public License v3.0 — see [LICENSE](LICENSE) for details.
