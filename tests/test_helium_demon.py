import numpy as np
import pytest

from helium_demon import Gas, HELIUM_MASS, K_B, maxwell_speed_pdf_2d


def gas_at(pos, vel, **kwargs):
    return Gas(np.array(pos, float), np.array(vel, float), **kwargs)


def test_equal_mass_elastic_collision_conserves_energy_and_momentum():
    gas = gas_at([[4e-9, 5e-9], [4e-9+2.59e-10, 5e-9]],
                 [[300., 40.], [-100., 10.]], demon=False)
    p0 = gas.mass*gas.velocities.sum(axis=0)
    e0 = gas.kinetic_energy
    gas._collide()
    np.testing.assert_allclose(gas.mass*gas.velocities.sum(axis=0), p0, rtol=1e-14, atol=1e-40)
    assert gas.kinetic_energy == pytest.approx(e0, rel=1e-14)
    assert gas.velocities[0, 0] < gas.velocities[1, 0]


def test_neighboring_grid_cells_are_checked_and_overlap_removed():
    # Radius is 0.13 nm; the overlapping pair straddles a grid-cell boundary.
    r = 0.13e-9
    gas = gas_at([[2*r-0.01*r, 3*r], [2*r+1.9*r, 3*r]],
                 [[0., 0.], [0., 0.]], radius=r, width=10*r, height=10*r, demon=False)
    gas._collide()
    distance = np.linalg.norm(gas.positions[1]-gas.positions[0])
    assert distance >= 2*r*(1-1e-10)


def test_random_initialization_respects_solid_divider_and_rejects_invalid_centers():
    gas = Gas.random(count=200, seed=9)
    mid = gas.width/2
    lo = (gas.height-gas.gate_opening)/2+gas.radius
    hi = (gas.height+gas.gate_opening)/2-gas.radius
    invalid = ((np.abs(gas.positions[:, 0]-mid) < gas.radius)
               & ((gas.positions[:, 1] < lo) | (gas.positions[:, 1] > hi)))
    assert not np.any(invalid)
    with pytest.raises(ValueError, match="solid divider"):
        gas_at([[mid, gas.height/4]], [[0., 0.]])


def test_collision_position_correction_does_not_push_particle_through_divider():
    r = 0.13e-9
    mid = 5e-9
    gas = gas_at([[mid-r-0.26e-9, 5e-9], [mid-r-0.01e-9, 5e-9]],
                 [[0., 0.], [0., 0.]], radius=r, width=10e-9, height=10e-9)
    gas._collide()
    assert np.all(gas.positions[:, 0] <= mid-r)
    assert np.linalg.norm(gas.positions[1]-gas.positions[0]) >= 2*r*(1-1e-10)


def test_solid_divider_face_reflects_without_losing_energy():
    r = 0.13e-9
    mid = 5e-9
    gas = gas_at([[mid-r-0.01e-9, 2.5e-9]], [[1000., 20.]],
                 radius=r, width=10e-9, height=10e-9)
    energy = gas.kinetic_energy
    gas.step(2e-13)
    assert gas.positions[0, 0] <= mid-r
    assert gas.velocities[0, 0] < 0
    assert gas.kinetic_energy == pytest.approx(energy, rel=1e-14)


def test_wall_reflection_uses_particle_radius_and_preserves_speed():
    r = 0.13e-9
    gas = gas_at([[r+0.01e-9, 4e-9]], [[-100., 50.]], radius=r,
                 width=10e-9, height=10e-9, demon=False)
    energy = gas.kinetic_energy
    gas.step(2e-12)
    assert gas.positions[0, 0] >= r
    assert gas.velocities[0, 0] > 0
    assert gas.kinetic_energy == pytest.approx(energy, rel=1e-14)


def test_gate_only_changes_velocity_on_crossing_and_filters_direction():
    r = 0.13e-9
    common = dict(width=10e-9, height=10e-9, radius=r, gate_opening=2e-9,
                  speed_threshold=150., demon=True)
    mid = 5e-9
    y = 5e-9
    # Slow particle from the left is reflected; its speed is compared at crossing.
    slow = gas_at([[mid-0.01e-9, y]], [[100., 0.]], **common)
    slow.step(1e-13)
    assert slow.positions[0, 0] < mid
    assert slow.velocities[0, 0] < 0
    # Fast particle from the left passes through the aperture.
    fast = gas_at([[mid-0.01e-9, y]], [[300., 0.]], **common)
    fast.step(1e-13)
    assert fast.positions[0, 0] > mid
    assert fast.velocities[0, 0] > 0
    # The same fast particle from the right is rejected.
    from_right = gas_at([[mid+0.01e-9, y]], [[-300., 0.]], **common)
    from_right.step(1e-13)
    assert from_right.positions[0, 0] > mid
    assert from_right.velocities[0, 0] > 0
    # The circle center must clear the aperture with a full radius of margin.
    edge = gas_at([[mid-0.14e-9, 6.0e-9]], [[300., 0.]], **common)
    edge.step(5e-13)
    assert edge.positions[0, 0] < mid
    assert edge.velocities[0, 0] < 0


def test_gate_counts_selection_by_direction_and_speed():
    mid = 5e-9
    gas = gas_at(
        [[mid + 0.02e-9, 4e-9], [mid + 0.02e-9, 6e-9],
         [mid - 0.02e-9, 4e-9], [mid - 0.02e-9, 6e-9]],
        [[-300., 0.], [-100., 0.], [300., 0.], [100., 0.]],
        width=10e-9, height=10e-9, radius=0.13e-9,
        gate_opening=4e-9, speed_threshold=150.,
    )
    gas.positions[:, 0] = [mid - 0.02e-9, mid - 0.02e-9,
                           mid + 0.02e-9, mid + 0.02e-9]
    gas._gate_crossings(np.array([mid + 0.02e-9, mid + 0.02e-9,
                                  mid - 0.02e-9, mid - 0.02e-9]))
    assert gas.gate_counts["left_to_right_fast_accepted"] == 1
    assert gas.gate_counts["left_to_right_slow_rejected"] == 1
    assert gas.gate_counts["right_to_left_slow_accepted"] == 1
    assert gas.gate_counts["right_to_left_fast_rejected"] == 1
    assert sum(gas.gate_counts.values()) == 4


def test_substeps_keep_fast_particles_in_box_without_tunneling_or_energy_drift():
    gas = Gas.random(count=30, temperature=300., seed=11)
    e0 = gas.kinetic_energy
    for _ in range(200):
        gas.step(2e-13)
        assert np.all(gas.positions[:, 0] >= gas.radius)
        assert np.all(gas.positions[:, 0] <= gas.width-gas.radius)
        assert np.all(gas.positions[:, 1] >= gas.radius)
        assert np.all(gas.positions[:, 1] <= gas.height-gas.radius)
        d = gas.positions[:, None, :] - gas.positions[None, :, :]
        dist2 = np.sum(d*d, axis=2)
        dist2[np.diag_indices_from(dist2)] = np.inf
        assert np.min(dist2) >= (2*gas.radius)**2*(1-1e-8)
    assert gas.kinetic_energy == pytest.approx(e0, rel=1e-10)


def test_two_dimensional_speed_density_and_temperature_scale():
    t, m = 300., HELIUM_MASS
    v = np.linspace(0, 20_000, 100_000)
    density = maxwell_speed_pdf_2d(v, t, m)
    integral = np.sum((density[:-1]+density[1:])*np.diff(v)/2)
    assert integral == pytest.approx(1., rel=2e-5)
    gas = gas_at([[3e-9, 3e-9], [7e-9, 7e-9]], [[100., 0.], [-100., 0.]], demon=False)
    # 2D kinetic energy convention: E/N = k_B T for an equilibrium ideal gas.
    assert gas.temperature == pytest.approx(HELIUM_MASS*100**2/K_B)
