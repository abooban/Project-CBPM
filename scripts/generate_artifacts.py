"""Generate deterministic portfolio figures and a short simulation animation."""
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
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
steps = 900
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

# Smooth animation at fewer frames; each frame advances ten stable time steps.
anim_gas = Gas.random(count=65, temperature=300, seed=19)
fig, ax = plt.subplots(figsize=(5.2, 5.2))
mid = anim_gas.width/2
lo, hi = (anim_gas.height-anim_gas.gate_opening)/2, (anim_gas.height+anim_gas.gate_opening)/2
ax.plot([mid, mid], [0, lo], color="#253746", lw=4, solid_capstyle="butt")
ax.plot([mid, mid], [hi, anim_gas.height], color="#253746", lw=4, solid_capstyle="butt")
ax.set(xlim=(0, anim_gas.width), ylim=(0, anim_gas.height), xticks=[], yticks=[], title="2-D helium gas · idealized Maxwell demon")
ax.set_aspect("equal")
sc = ax.scatter(anim_gas.positions[:, 0], anim_gas.positions[:, 1], s=17, c="#2b8cbe", edgecolors="white", linewidths=.25)
label = ax.text(.04, .96, "", transform=ax.transAxes, va="top", color="#34495e")
def update(frame):
    for _ in range(10):
        anim_gas.step(dt)
    sc.set_offsets(anim_gas.positions)
    tl, tr = anim_gas.regional_temperatures()
    label.set_text(f"Left: {tl:.0f} K   Right: {tr:.0f} K")
    return sc, label
animation = FuncAnimation(fig, update, frames=100, interval=45, blit=True)
animation.save(OUT/"simulation.gif", writer=PillowWriter(fps=20), dpi=95)
plt.close(fig)
