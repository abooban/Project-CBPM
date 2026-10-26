"""Run paired, seeded demon/no-demon trials and plot temperature separation."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from helium_demon import Gas

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "assets"


def run_trial(seed: int, demon: bool, steps: int, dt: float,
              sample_every: int) -> list[dict[str, float | int | str]]:
    initial = Gas.random(count=55, temperature=300, seed=seed, demon=demon)
    initial_energy = initial.kinetic_energy
    rows = []
    for step in range(1, steps + 1):
        initial.step(dt)
        if step % sample_every and step != steps:
            continue
        tl, tr = initial.regional_temperatures()
        row: dict[str, float | int | str] = {
            "seed": seed,
            "demon": demon,
            "time_ns": step * dt * 1e9,
            "step": step,
            "temperature_left_K": tl,
            "temperature_right_K": tr,
            "delta_temperature_K": tr - tl,
            "count_left": int(np.sum(initial.positions[:, 0] < initial.width / 2)),
            "count_right": int(np.sum(initial.positions[:, 0] >= initial.width / 2)),
            "kinetic_energy_J": initial.kinetic_energy,
            "relative_energy_drift": (initial.kinetic_energy - initial_energy) / initial_energy,
        }
        row.update(initial.gate_counts)
        rows.append(row)
    return rows


def plot_results(rows: list[dict[str, float | int | str]]) -> None:
    conditions = [(True, "Demon"), (False, "No demon")]
    fig, axes = plt.subplots(3, 1, figsize=(9, 10), sharex=True, layout="constrained")
    colors = {"Demon": "#e6550d", "No demon": "#3182bd"}
    def is_demon(row):
        return row["demon"] is True or str(row["demon"]).lower() in ("true", "1")
    for demon, label in conditions:
        selected = [r for r in rows if is_demon(r) is demon]
        seeds = sorted({int(r["seed"]) for r in selected})
        series = []
        for seed in seeds:
            trial = sorted((r for r in selected if int(r["seed"]) == seed),
                           key=lambda r: float(r["time_ns"]))
            times = np.array([float(r["time_ns"]) for r in trial])
            if not series:
                grid = times
            series.append(np.array([[float(r["temperature_left_K"]),
                                     float(r["temperature_right_K"]),
                                     float(r["delta_temperature_K"]),
                                     float(r["count_left"]),
                                     float(r["count_right"])] for r in trial]))
        values = np.stack(series)
        mean = np.mean(values, axis=0)
        sem = np.std(values, axis=0, ddof=1) / np.sqrt(len(series)) if len(series) > 1 else np.zeros_like(mean)
        axes[0].plot(grid, mean[:, 0], color=colors[label], label=f"{label}: left")
        axes[0].fill_between(grid, mean[:, 0] - sem[:, 0], mean[:, 0] + sem[:, 0],
                             color=colors[label], alpha=.12)
        axes[0].plot(grid, mean[:, 1], color=colors[label], ls="--", label=f"{label}: right")
        axes[0].fill_between(grid, mean[:, 1] - sem[:, 1], mean[:, 1] + sem[:, 1],
                             color=colors[label], alpha=.12)
        axes[1].plot(grid, mean[:, 2], color=colors[label], label=label)
        axes[1].fill_between(grid, mean[:, 2] - sem[:, 2], mean[:, 2] + sem[:, 2],
                             color=colors[label], alpha=.18)
        axes[2].plot(grid, mean[:, 3], color=colors[label], label=f"{label}: left")
        axes[2].fill_between(grid, mean[:, 3] - sem[:, 3], mean[:, 3] + sem[:, 3],
                             color=colors[label], alpha=.12)
        axes[2].plot(grid, mean[:, 4], color=colors[label], ls="--", label=f"{label}: right")
        axes[2].fill_between(grid, mean[:, 4] - sem[:, 4], mean[:, 4] + sem[:, 4],
                             color=colors[label], alpha=.12)
    axes[0].axhline(300, color="#777", lw=1, ls=":")
    axes[0].set_ylabel("Regional temperature (K)")
    axes[0].legend(ncol=2, frameon=False)
    axes[0].set_title("Paired seeded trials · mean ± standard error")
    axes[1].axhline(0, color="#777", lw=1, ls=":")
    axes[1].set_ylabel("ΔT = right − left (K)")
    axes[1].legend(frameon=False)
    axes[2].set_ylabel("Particle count")
    axes[2].set_xlabel("Physical simulation time (ns)")
    axes[2].legend(frameon=False)
    OUT.mkdir(exist_ok=True)
    fig.savefig(OUT / "demon_experiment.png", dpi=160)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--steps", type=int, default=4000)
    parser.add_argument("--seeds", type=int, nargs="+", default=[7, 19, 42, 3, 11, 23, 81, 101])
    parser.add_argument("--dt", type=float, default=2e-13)
    parser.add_argument("--sample-every", type=int, default=40)
    parser.add_argument("--append", action="store_true", help="append trials to existing CSV and replot")
    args = parser.parse_args()
    if args.steps < 1 or args.sample_every < 1:
        parser.error("steps and sample-every must be positive")

    rows = []
    csv_path = OUT / "demon_experiment.csv"
    if args.append and csv_path.exists():
        with csv_path.open(newline="", encoding="utf-8") as stream:
            rows.extend(csv.DictReader(stream))
    for seed in args.seeds:
        for demon in (True, False):
            print(f"Running seed={seed}, demon={demon}", flush=True)
            rows.extend(run_trial(seed, demon, args.steps, args.dt, args.sample_every))
    OUT.mkdir(exist_ok=True)
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    plot_results(rows)

    def is_demon(row):
        return row["demon"] is True or str(row["demon"]).lower() in ("true", "1")
    for demon, label in ((True, "demon"), (False, "no demon")):
        finals = [r for r in rows if is_demon(r) is demon and int(r["step"]) == args.steps]
        dts = [float(r["delta_temperature_K"]) for r in finals]
        drifts = [float(r["relative_energy_drift"]) for r in finals]
        print(f"{label}: final ΔT mean={np.mean(dts):.2f} K, SD={np.std(dts, ddof=1):.2f} K; "
              f"max |relative energy drift|={max(map(abs, drifts)):.3g}")
        if demon:
            totals = {key: sum(int(r[key]) for r in finals) for key in finals[0]
                      if key.endswith(("_accepted", "_rejected"))}
            print("demon final cumulative crossing counts:", totals)
    final_by_seed = {(int(r["seed"]), is_demon(r)): r for r in rows
                     if int(r["step"]) == args.steps}
    paired_effects = [float(final_by_seed[(seed, True)]["delta_temperature_K"])
                      - float(final_by_seed[(seed, False)]["delta_temperature_K"])
                      for seed in sorted({key[0] for key in final_by_seed})
                      if (seed, True) in final_by_seed and (seed, False) in final_by_seed]
    print(f"paired demon-minus-control ΔT: mean={np.mean(paired_effects):.2f} K, "
          f"SD={np.std(paired_effects, ddof=1):.2f} K, n={len(paired_effects)}")
    print(f"Saved {csv_path} and {OUT / 'demon_experiment.png'}")


if __name__ == "__main__":
    main()
