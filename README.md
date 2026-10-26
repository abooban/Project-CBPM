# Project CBPM · 2-D helium gas

An educational hard-disk simulation of helium atoms, elastic collisions, thermal speeds, and an idealized Maxwell-demon gate. The project began as a **first-year 2024 student submission** for *Introduction to Computer-Based Physical Modeling* in the International Physics Studies Program (IPSP) at Leipzig University. The code and documentation on this branch are a later portfolio improvement; they are not part of the original coursework.

![Short animation of the 2-D gas and central gate](assets/simulation.gif)

## What it demonstrates

- Equal-mass, elastic hard-disk impacts with momentum and kinetic-energy conservation.
- Radius-aware reflecting box walls and a central aperture. A demon rule acts only when a particle crosses the divider: fast particles go left-to-right and slow particles right-to-left.
- A reproducible 2-D thermal initialization and the **2-D Maxwell–Boltzmann (Rayleigh) speed density**, `f(v) = m v/(k_B T) exp[-m v²/(2 k_B T)]`.
- Kinetic temperature from equipartition in two dimensions (`E/N = k_B T` at equilibrium), with center-of-mass motion removed for temperature estimates.
- Adaptive ballistic substeps, a neighboring-cell broad phase, overlap correction, and collision impulses only for approaching pairs.

The demon is an educational idealization. Its measurement, memory, and reset costs are not modeled, so the temperature sorting does not violate the second law.

## Results

The seed-fixed plots below compare a finite random sample with the theoretical 2-D speed density and track energy and regional temperatures in a gate run.

| Speed distribution | Total kinetic energy | Regional temperatures |
|---|---|---|
| ![2-D speed distribution](assets/speed_distribution.png) | ![Kinetic energy over time](assets/kinetic_energy.png) | ![Regional temperature over time](assets/regional_temperature.png) |

## Run on Ubuntu or WSL

```bash
sudo apt update
sudo apt install python3 python3-venv git
git clone https://github.com/abooban/Project-CBPM.git
cd Project-CBPM
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
pytest
jupyter lab Project.ipynb
```

The notebook is intentionally lightweight and has no embedded output. It imports the tested simulation module. To regenerate the GIF and plots, run `python scripts/generate_artifacts.py`.

## What changed since the original notebook

The original notebook remains in Git history as submitted in 2024. Its saved animation and plot outputs are cleared in this branch to reduce checkout size; the current notebook is a concise entry point to the improved model.

The original implementation missed collisions across grid-cell boundaries, could respond repeatedly to separating overlaps, ignored particle radius at outer boundaries, and applied gate rules based on a particle's side rather than a detected crossing. Its stated speed law and temperature factor were 3-D formulas in a 2-D model. It also contained global-state and initialization inconsistencies. These have been replaced with a small importable simulation module and regression tests. The tests check elastic collision conservation, wall reflection, gate crossing, overlap handling, and bounded long-run stepping; they support those specific properties, not a claim of research-grade accuracy.

## Files

```text
Project.ipynb                 # short, reproducible notebook entry point
helium_demon.py               # simulation and 2-D statistical mechanics
tests/test_helium_demon.py    # lightweight regression tests
scripts/generate_artifacts.py # deterministic figures and animation
assets/                       # README visualizations
```
