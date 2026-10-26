"""Small 2-D hard-disk gas model for the Project-CBPM teaching project.

Lengths are metres, times seconds, velocities m/s, masses kg and energies J.
The demon is an explicitly idealized, information-using gate; its work/entropy
cost is not modeled, so it is not a physical counterexample to the second law.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

K_B = 1.380649e-23  # J/K, exact SI value
HELIUM_MASS = 4.002602e-3 / 6.02214076e23  # kg per atom
HELIUM_DIAMETER = 2.6e-10  # m; illustrative effective hard-sphere diameter


def maxwell_speed_pdf_2d(speed: np.ndarray | float, temperature: float,
                         mass: float = HELIUM_MASS) -> np.ndarray:
    """Normalized 2-D Maxwell (Rayleigh) speed density, in (m/s)^-1."""
    v = np.asarray(speed, dtype=float)
    if temperature <= 0 or mass <= 0:
        raise ValueError("temperature and mass must be positive")
    a = mass / (K_B * temperature)
    return a * np.maximum(v, 0.0) * np.exp(-0.5 * a * np.maximum(v, 0.0) ** 2)


@dataclass
class Gas:
    """Equal-mass hard disks in a rectangular box split by a demon gate."""
    positions: np.ndarray
    velocities: np.ndarray
    width: float = 20e-9
    height: float = 20e-9
    radius: float = HELIUM_DIAMETER / 2
    gate_opening: float = 3 * HELIUM_DIAMETER
    mass: float = HELIUM_MASS
    demon: bool = True
    speed_threshold: float | None = None

    def __post_init__(self) -> None:
        self.positions = np.asarray(self.positions, dtype=float).copy()
        self.velocities = np.asarray(self.velocities, dtype=float).copy()
        if self.positions.ndim != 2 or self.positions.shape[1] != 2:
            raise ValueError("positions must have shape (N, 2)")
        if self.velocities.shape != self.positions.shape:
            raise ValueError("velocities must have the same shape as positions")
        if self.width <= 2*self.radius or self.height <= 2*self.radius:
            raise ValueError("box dimensions must exceed the particle diameter")
        if not 0 < self.gate_opening <= self.height - 2*self.radius:
            raise ValueError("gate opening must fit within the box")
        if self.mass <= 0:
            raise ValueError("mass must be positive")
        if np.any(self.positions[:, 0] < self.radius) or np.any(self.positions[:, 0] > self.width-self.radius) or np.any(self.positions[:, 1] < self.radius) or np.any(self.positions[:, 1] > self.height-self.radius):
            raise ValueError("initial particle centers must be inside the radius-adjusted box")
        mid = self.width / 2
        gate_lo = (self.height-self.gate_opening)/2 + self.radius
        gate_hi = (self.height+self.gate_opening)/2 - self.radius
        in_solid_divider = (np.abs(self.positions[:, 0]-mid) < self.radius) & (
            (self.positions[:, 1] < gate_lo) | (self.positions[:, 1] > gate_hi))
        if np.any(in_solid_divider):
            raise ValueError("initial particle centers must not overlap the solid divider")
        if self.speed_threshold is None:
            self.speed_threshold = float(np.sqrt(2*K_B*300/self.mass))
        self.gate_counts = {
            f"{direction}_{speed}_{outcome}": 0
            for direction in ("left_to_right", "right_to_left")
            for speed in ("fast", "slow")
            for outcome in ("accepted", "rejected")
        }

    @classmethod
    def random(cls, count: int = 100, temperature: float = 300.0,
               seed: int = 2024, **kwargs) -> "Gas":
        """Create non-overlapping positions and 2-D thermal Gaussian velocities."""
        if count < 1 or temperature <= 0:
            raise ValueError("count and temperature must be positive")
        rng = np.random.default_rng(seed)
        width = kwargs.get("width", 20e-9)
        height = kwargs.get("height", 20e-9)
        radius = kwargs.get("radius", HELIUM_DIAMETER/2)
        gate_opening = kwargs.get("gate_opening", 3*HELIUM_DIAMETER)
        mid = width/2
        gate_lo = (height-gate_opening)/2 + radius
        gate_hi = (height+gate_opening)/2 - radius
        pos: list[np.ndarray] = []
        min_d2 = (2*radius)**2
        for _ in range(count):
            for attempt in range(100_000):
                candidate = np.array([rng.uniform(radius, width-radius), rng.uniform(radius, height-radius)])
                if abs(candidate[0]-mid) < radius and not gate_lo <= candidate[1] <= gate_hi:
                    continue
                if all(np.sum((candidate-p)**2) >= min_d2 for p in pos):
                    pos.append(candidate)
                    break
            else:
                raise ValueError("could not pack particles; reduce count or particle radius")
        mass = kwargs.get("mass", HELIUM_MASS)
        velocities = rng.normal(0.0, math.sqrt(K_B*temperature/mass), size=(count, 2))
        return cls(np.asarray(pos), velocities, **kwargs)

    @property
    def count(self) -> int:
        return len(self.positions)

    @property
    def kinetic_energy(self) -> float:
        return float(0.5*self.mass*np.sum(self.velocities**2))

    @property
    def temperature(self) -> float:
        """2-D kinetic temperature, excluding center-of-mass drift."""
        if self.count == 0:
            return float("nan")
        dof = max(1, 2*self.count - 2)
        v = self.velocities - np.mean(self.velocities, axis=0)
        return float(self.mass*np.sum(v*v)/(dof*K_B))

    def regional_temperatures(self) -> tuple[float, float]:
        """Return left/right 2-D kinetic temperatures, NaN for empty regions."""
        out = []
        left = self.positions[:, 0] < self.width/2
        for mask in (left, ~left):
            n = int(np.sum(mask))
            if n < 2:
                out.append(float("nan"))
                continue
            v = self.velocities[mask]
            v = v - np.mean(v, axis=0)
            out.append(float(self.mass*np.sum(v*v)/((2*n-2)*K_B)))
        return out[0], out[1]

    def _gate_crossings(self, old_x: np.ndarray, old_y: np.ndarray | None = None) -> None:
        """Resolve actual divider crossings; blocked particles reflect at x=mid."""
        mid = self.width/2
        new_x = self.positions[:, 0]
        if old_y is None:
            old_y = self.positions[:, 1]
        crossed = ((old_x < mid) & (new_x >= mid)) | ((old_x > mid) & (new_x <= mid))
        lo = (self.height-self.gate_opening)/2 + self.radius
        hi = (self.height+self.gate_opening)/2 - self.radius
        for i in np.flatnonzero(crossed):
            left_to_right = old_x[i] < mid
            fast = np.linalg.norm(self.velocities[i]) >= self.speed_threshold
            direction = "left_to_right" if left_to_right else "right_to_left"
            speed_class = "fast" if fast else "slow"
            allowed_speed = (not self.demon or (left_to_right == fast))
            vx = self.velocities[i, 0]
            if vx == 0:
                passage_y = (self.positions[i, 1], self.positions[i, 1])
            else:
                # A disk occupies a 2r-wide strip while passing the divider.
                # Check both edges so it cannot clear the midplane and then
                # clip a gate corner before its center clears the divider.
                y_entry = old_y[i] + ((mid-self.radius)-old_x[i]) * self.velocities[i, 1]/vx
                y_exit = old_y[i] + ((mid+self.radius)-old_x[i]) * self.velocities[i, 1]/vx
                passage_y = (y_entry, y_exit)
            allowed_y = min(passage_y) >= lo and max(passage_y) <= hi
            accepted = allowed_y and allowed_speed
            outcome = "accepted" if accepted else "rejected"
            self.gate_counts[f"{direction}_{speed_class}_{outcome}"] += 1
            if not accepted:
                self.positions[i, 0] = mid - self.radius if old_x[i] < mid else mid + self.radius
                self.velocities[i, 0] *= -1

    def _walls(self, old_x: np.ndarray, old_y: np.ndarray | None = None) -> None:
        r = self.radius
        if old_y is None:
            old_y = self.positions[:, 1].copy()
        # Reflecting outer box, clamped to physically valid center coordinates.
        for axis, extent in ((0, self.width), (1, self.height)):
            low, high = r, extent-r
            below, above = self.positions[:, axis] < low, self.positions[:, axis] > high
            self.positions[below, axis] = 2*low-self.positions[below, axis]
            self.velocities[below, axis] = np.abs(self.velocities[below, axis])
            self.positions[above, axis] = 2*high-self.positions[above, axis]
            self.velocities[above, axis] = -np.abs(self.velocities[above, axis])
            self.positions[:, axis] = np.clip(self.positions[:, axis], low, high)
        mid = self.width/2
        lo = (self.height-self.gate_opening)/2+r
        hi = (self.height+self.gate_opening)/2-r
        # Reflect at the finite-radius faces of the solid divider. A disk
        # entering the strip through the aperture also reflects off its edges.
        vx = self.velocities[:, 0]
        moving_right = vx > 0
        moving_left = vx < 0
        left_face_crossed = moving_right & (old_x <= mid-r) & (self.positions[:, 0] > mid-r)
        right_face_crossed = moving_left & (old_x >= mid+r) & (self.positions[:, 0] < mid+r)
        for crossed_face, face_x in ((left_face_crossed, mid-r), (right_face_crossed, mid+r)):
            indices = np.flatnonzero(crossed_face)
            for i in indices:
                y_at_face = old_y[i] + (face_x-old_x[i])*self.velocities[i, 1]/vx[i]
                if not lo <= y_at_face <= hi:
                    self.positions[i, 0] = face_x
                    self.velocities[i, 0] *= -1
        in_strip = np.abs(self.positions[:, 0]-mid) < r
        below_gate = in_strip & (self.positions[:, 1] < lo)
        above_gate = in_strip & (self.positions[:, 1] > hi)
        self.positions[below_gate, 1] = 2*lo-self.positions[below_gate, 1]
        self.velocities[below_gate, 1] = np.abs(self.velocities[below_gate, 1])
        self.positions[above_gate, 1] = 2*hi-self.positions[above_gate, 1]
        self.velocities[above_gate, 1] = -np.abs(self.velocities[above_gate, 1])
        self._gate_crossings(old_x, old_y)

    def _candidate_pairs(self) -> set[tuple[int, int]]:
        """Uniform-grid broad phase; adjacent cells are included exactly once."""
        cell_size = 2*self.radius
        cells: dict[tuple[int, int], list[int]] = {}
        for i, (x, y) in enumerate(self.positions):
            key = (int(x//cell_size), int(y//cell_size))
            cells.setdefault(key, []).append(i)
        pairs: set[tuple[int, int]] = set()
        for (cx, cy), members in cells.items():
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    for i in members:
                        for j in cells.get((cx+dx, cy+dy), ()):
                            if i < j:
                                pairs.add((i, j))
        return pairs

    def _collide(self) -> None:
        # A few positional-correction sweeps resolve small multi-particle overlaps.
        # Impulses are applied only for approaching pairs, preventing repeat bounces.
        diameter = 2*self.radius
        lo = (self.height-self.gate_opening)/2+self.radius
        hi = (self.height+self.gate_opening)/2-self.radius
        for _ in range(8):
            changed = False
            mid = self.width/2
            left_side = self.positions[:, 0] < mid
            for i, j in sorted(self._candidate_pairs()):
                if (self.positions[i, 0]-mid)*(self.positions[j, 0]-mid) < 0:
                    lo = (self.height-self.gate_opening)/2+self.radius
                    hi = (self.height+self.gate_opening)/2-self.radius
                    if not (lo <= self.positions[i, 1] <= hi and lo <= self.positions[j, 1] <= hi):
                        continue  # the solid portion of the divider separates them
                delta = self.positions[j]-self.positions[i]
                d2 = float(delta@delta)
                if d2 >= diameter*diameter:
                    continue
                dist = math.sqrt(d2)
                if dist < 1e-15*self.radius:
                    # Deterministic fallback for exactly coincident centers.
                    normal = np.array([1.0, 0.0]) if (i+j) % 2 == 0 else np.array([0.0, 1.0])
                else:
                    normal = delta/dist
                overlap = diameter-dist
                at_divider_i = ((self.positions[i, 1] < lo or self.positions[i, 1] > hi)
                                and abs(abs(self.positions[i, 0]-mid)-self.radius) <= 1e-12*self.radius)
                at_divider_j = ((self.positions[j, 1] < lo or self.positions[j, 1] > hi)
                                and abs(abs(self.positions[j, 0]-mid)-self.radius) <= 1e-12*self.radius)
                correction = overlap + 1e-12*self.radius
                if at_divider_i and not at_divider_j:
                    self.positions[j] += normal*correction
                elif at_divider_j and not at_divider_i:
                    self.positions[i] -= normal*correction
                else:
                    self.positions[i] -= normal*(correction/2)
                    self.positions[j] += normal*(correction/2)
                relative_normal = float((self.velocities[j]-self.velocities[i])@normal)
                if relative_normal < 0:
                    impulse = -relative_normal
                    self.velocities[i] -= impulse*normal
                    self.velocities[j] += impulse*normal
                changed = True
            # Positional overlap corrections must not push disks into the solid
            # divider. Ballistic crossings are handled separately by the gate.
            outside_aperture = (self.positions[:, 1] < lo) | (self.positions[:, 1] > hi)
            left_crossed = left_side & (self.positions[:, 0] >= mid)
            right_crossed = ~left_side & (self.positions[:, 0] <= mid)
            self.positions[left_crossed, 0] = mid-self.radius
            self.positions[right_crossed, 0] = mid+self.radius
            left_blocked = left_side & outside_aperture & (self.positions[:, 0] > mid-self.radius)
            right_blocked = ~left_side & outside_aperture & (self.positions[:, 0] < mid+self.radius)
            self.positions[left_blocked, 0] = mid-self.radius
            self.positions[right_blocked, 0] = mid+self.radius
            if not changed:
                break
        # Positional corrections can approach a box edge; keep centers inside.
        np.clip(self.positions[:, 0], self.radius, self.width-self.radius, out=self.positions[:, 0])
        np.clip(self.positions[:, 1], self.radius, self.height-self.radius, out=self.positions[:, 1])

    def step(self, dt: float) -> None:
        """Advance by dt with adaptive ballistic substeps and hard-disk impacts."""
        if not np.isfinite(dt) or dt <= 0:
            raise ValueError("dt must be finite and positive")
        vmax = float(np.max(np.linalg.norm(self.velocities, axis=1), initial=0.0))
        # Limit travel to 1/5 radius so discrete overlap detection cannot skip a hit.
        n_sub = max(1, int(math.ceil(vmax*dt/(0.2*self.radius))))
        if n_sub > 200_000:
            raise ValueError("dt is too large for stable stepping; reduce dt")
        h = dt/n_sub
        for _ in range(n_sub):
            old_positions = self.positions.copy()
            self.positions += self.velocities*h
            self._walls(old_positions[:, 0], old_positions[:, 1])
            self._collide()


def initial_gas(count: int = 100, temperature: float = 300.0,
                seed: int = 2024, **kwargs) -> Gas:
    return Gas.random(count=count, temperature=temperature, seed=seed, **kwargs)
