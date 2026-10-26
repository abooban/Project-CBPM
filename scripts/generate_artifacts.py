"""Generate deterministic portfolio figures and an actual-trajectory animation."""
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.colors import Normalize
from matplotlib.patches import Rectangle
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from helium_demon import Gas, K_B, maxwell_speed_pdf_2d

OUT = Path(__file__).resolve().parents[1] / "assets"
OUT.mkdir(exist_ok=True)
plt.rcParams.update({"figure.dpi": 130, "font.size": 10, "axes.spines.top": False, "axes.spines.right": False})

# Distribution plot: sampled 2-D Maxwell-Boltzmann velocities.
gas = Gas.random(count=1200, temperature=300, seed=2024, demon=False)
speeds = np.linalg.norm(gas.velocities, axis=1)
fig, ax = plt.subplots(figsize=(7.2, 4.2))
ax.hist(speeds, bins=35, density=True, color="#6baed6", alpha=.72, edgecolor="white", label="Sampled particle speeds")
v = np.linspace(0, speeds.max()*1.05, 400)
ax.plot(v, maxwell_speed_pdf_2d(v, 300, gas.mass), color="#e6550d", lw=2.5, label="2-D Maxwell–Boltzmann density")
ax.set(xlabel="Speed (m s$^{-1}$)", ylabel="Probability density (s m$^{-1}$)", title="Thermal speed distribution · 2-D helium model")
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(OUT/"speed_distribution.png", bbox_inches="tight")
plt.close(fig)

# Demon run and energy histories, with seeded initial conditions.
sim = Gas.random(count=55, temperature=300, seed=7)
dt = 2e-13
steps = 4000
times, energies, left_t, right_t = [], [], [], []
for k in range(steps):
    sim.step(dt)
    if k % 3 == 0:
        times.append((k+1)*dt*1e9)
        energies.append(sim.kinetic_energy)
        tl, tr = sim.regional_temperatures()
        left_t.append(tl)
        right_t.append(tr)
fig, ax = plt.subplots(figsize=(7.2, 4.2))
ax.plot(times, energies, color="#3182bd", lw=1.5)
ax.set(xlabel="Time (ns)", ylabel="Total kinetic energy (J)", title="Total kinetic energy · elastic hard-disk model")
ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
fig.tight_layout()
fig.savefig(OUT/"kinetic_energy.png", bbox_inches="tight")
plt.close(fig)

fig, ax = plt.subplots(figsize=(7.2, 4.2))
ax.plot(times, left_t, color="#3182bd", label="Left")
ax.plot(times, right_t, color="#e6550d", label="Right")
ax.set(xlabel="Time (ns)", ylabel="Kinetic temperature (K)", title="Regional temperatures · idealized demon gate")
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(OUT/"regional_temperature.png", bbox_inches="tight")
plt.close(fig)

# The 16-second loop advances 20 physical simulation steps per displayed frame.
# FuncAnimation records the actual Gas state after stepping; it does not tween
# positions, so particles never receive fabricated paths through a boundary.
anim_gas = Gas.random(count=65, temperature=300, seed=19)
fps = 12.5
frames = 200
steps_per_frame = 20
speed_limit = 3*np.sqrt(2*K_B*300/anim_gas.mass)
speed_norm = Normalize(vmin=0, vmax=speed_limit, clip=True)
cmap = plt.get_cmap("plasma")
fig, ax = plt.subplots(figsize=(8.2, 7.2), layout="constrained")
mid = anim_gas.width/2
lo, hi = (anim_gas.height-anim_gas.gate_opening)/2, (anim_gas.height+anim_gas.gate_opening)/2
ax.add_patch(Rectangle((0, 0), mid, anim_gas.height, color="#4c9ac1", alpha=.055, zorder=0))
ax.add_patch(Rectangle((mid, 0), mid, anim_gas.height, color="#e58b53", alpha=.055, zorder=0))
ax.add_patch(Rectangle((0, 0), anim_gas.width, anim_gas.height, fill=False,
                       edgecolor="#243746", linewidth=2.0, zorder=2))
ax.plot([mid, mid], [0, lo], color="#243746", lw=5, solid_capstyle="butt", zorder=3)
ax.plot([mid, mid], [hi, anim_gas.height], color="#243746", lw=5, solid_capstyle="butt", zorder=3)
ax.add_patch(Rectangle((mid-.10e-9, lo), .20e-9, anim_gas.gate_opening,
                       facecolor="#f2c14e", edgecolor="#b07816", linewidth=1.25,
                       alpha=.78, zorder=2.5))
ax.annotate("GATE OPENING", xy=(mid, (lo+hi)/2), xytext=(mid+1.0e-9, (lo+hi)/2+1.0e-9),
            ha="left", va="bottom", fontsize=7.5, weight="bold", color="#76551b",
            arrowprops={"arrowstyle":"-", "color":"#b07816", "linewidth":1.1},
            bbox={"boxstyle":"round,pad=.22", "facecolor":"#fff8df", "edgecolor":"#e3c36d", "alpha":.96},
            zorder=4)
ax.set(xlim=(0, anim_gas.width), ylim=(0, anim_gas.height), xticks=[], yticks=[])
ax.set_aspect("equal")
ax.set_title("Fast particles pass left → right · slow particles pass right → left",
             fontsize=10, color="#405363", pad=10)
speeds = np.linalg.norm(anim_gas.velocities, axis=1)
sc = ax.scatter(anim_gas.positions[:, 0], anim_gas.positions[:, 1], s=42,
                c=speeds, cmap=cmap, norm=speed_norm, edgecolors="white",
                linewidths=.65, alpha=.96, zorder=5)
sm = plt.cm.ScalarMappable(norm=speed_norm, cmap=cmap)
sm.set_array([])
cbar = fig.colorbar(sm, ax=ax, pad=.035, shrink=.82, aspect=28)
cbar.set_label("Particle speed (m/s)", rotation=90, labelpad=10)
cbar.outline.set_linewidth(.5)
left_label = ax.text(.035, .965, "LEFT  ·  -- K  ·  n=--", transform=ax.transAxes,
                     va="top", ha="left", color="#173b50", fontsize=10, weight="bold",
                     bbox={"boxstyle":"round,pad=.35", "facecolor":"#eaf5fa", "edgecolor":"#9cc5d6", "alpha":.96},
                     zorder=6)
right_label = ax.text(.965, .965, "RIGHT  ·  -- K  ·  n=--", transform=ax.transAxes,
                      va="top", ha="right", color="#713d22", fontsize=10, weight="bold",
                      bbox={"boxstyle":"round,pad=.35", "facecolor":"#fff0e5", "edgecolor":"#e5b28e", "alpha":.96},
                      zorder=6)
path_lengths = np.zeros(anim_gas.count)
initial_energy = anim_gas.kinetic_energy

def validate_simulation_state():
    p = anim_gas.positions
    if not np.all(np.isfinite(p)) or not np.all(np.isfinite(anim_gas.velocities)):
        raise RuntimeError("non-finite particle state while generating animation")
    r = anim_gas.radius
    if (np.any(p[:, 0] < r) or np.any(p[:, 0] > anim_gas.width-r)
            or np.any(p[:, 1] < r) or np.any(p[:, 1] > anim_gas.height-r)):
        raise RuntimeError("particle center crossed an outer wall")
    delta = p[:, None, :] - p[None, :, :]
    distance2 = np.sum(delta*delta, axis=2)
    distance2[np.diag_indices_from(distance2)] = np.inf
    if np.min(distance2) < (2*r)**2*(1-1e-8):
        raise RuntimeError("particle overlap detected while generating animation")

def update(frame):
    for _ in range(steps_per_frame):
        previous = anim_gas.positions.copy()
        anim_gas.step(dt)
        path_lengths[:] += np.linalg.norm(anim_gas.positions-previous, axis=1)
        validate_simulation_state()
    sc.set_offsets(anim_gas.positions)
    sc.set_array(np.linalg.norm(anim_gas.velocities, axis=1))
    tl, tr = anim_gas.regional_temperatures()
    left_count = int(np.sum(anim_gas.positions[:, 0] < anim_gas.width/2))
    right_count = anim_gas.count-left_count
    left_label.set_text(f"LEFT  ·  {tl:.0f} K  ·  n={left_count}" if np.isfinite(tl) else f"LEFT  ·  -- K  ·  n={left_count}")
    right_label.set_text(f"RIGHT  ·  {tr:.0f} K  ·  n={right_count}" if np.isfinite(tr) else f"RIGHT  ·  -- K  ·  n={right_count}")
    return sc, left_label, right_label

animation = FuncAnimation(fig, update, frames=frames, interval=1000/fps, blit=True)
animation.save(OUT/"simulation.gif", writer=PillowWriter(fps=fps), dpi=88)
final_relative_energy_drift = abs(anim_gas.kinetic_energy-initial_energy)/initial_energy
if final_relative_energy_drift > 1e-8:
    raise RuntimeError(f"kinetic energy drift exceeded tolerance: {final_relative_energy_drift:.3g}")
if np.any(path_lengths < 1e-12):
    raise RuntimeError("one or more particles remained effectively stationary in the animation")
print(f"Animation validated: {frames} frames at {fps} FPS, {steps_per_frame} simulation steps/frame, "
      f"{final_relative_energy_drift:.2e} relative energy drift, minimum path {path_lengths.min():.3g} m")
plt.close(fig)
