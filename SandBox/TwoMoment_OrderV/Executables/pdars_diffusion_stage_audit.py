#!/usr/bin/env python3

"""Verify the fixed-coefficient PDARS diffusion-limit stage algebra.

The scalar model retains the distinction present in Thornado between the
centered density flux carried by the evolved first moment and the direct
density-jump stabilization evaluated from the current explicit stage.
"""

from __future__ import annotations

import numpy as np


DT = 0.2
EQUILIBRIUM_FRACTION = 0.4
PRESSURE_GRADIENT = 0.7
DENSITY_DIVERGENCE = 1.2
MOMENTUM_OPACITY = 0.9
FACE_STABILIZATION = -0.3
INITIAL_INVARIANT = 1.3


def equilibrium_density(invariant: float) -> float:
    return EQUILIBRIUM_FRACTION * invariant


def fick_flux(invariant: float) -> float:
    return (
        -PRESSURE_GRADIENT
        * equilibrium_density(invariant)
        / MOMENTUM_OPACITY
    )


def centered_operator(invariant: float) -> float:
    return DENSITY_DIVERGENCE * fick_flux(invariant)


def penalty_operator(invariant: float) -> float:
    return FACE_STABILIZATION * equilibrium_density(invariant)


def transport(state: np.ndarray, epsilon: float) -> np.ndarray:
    density, first_moment, _matter = state
    return np.array(
        [
            DENSITY_DIVERGENCE * first_moment / epsilon
            + FACE_STABILIZATION * density,
            -PRESSURE_GRADIENT * density / epsilon,
            0.0,
        ]
    )


def implicit_collision(input_state: np.ndarray, stage_dt: float, epsilon: float) -> np.ndarray:
    density_input, first_moment_input, matter_input = input_state
    invariant = density_input + matter_input
    stiffness = stage_dt / epsilon**2
    density = (
        density_input + stiffness * equilibrium_density(invariant)
    ) / (1.0 + stiffness)
    first_moment = first_moment_input / (
        1.0 + stiffness * MOMENTUM_OPACITY
    )
    return np.array([density, first_moment, invariant - density])


def pdars_step(epsilon: float) -> tuple[np.ndarray, np.ndarray]:
    density = equilibrium_density(INITIAL_INVARIANT)
    state_1 = np.array(
        [
            density,
            epsilon * fick_flux(INITIAL_INVARIANT),
            INITIAL_INVARIANT - density,
        ]
    )
    transport_1 = transport(state_1, epsilon)

    input_2 = state_1 + DT * transport_1
    state_2 = implicit_collision(input_2, DT, epsilon)
    collision_2 = (state_2 - input_2) / DT
    transport_2 = transport(state_2, epsilon)

    input_3 = state_1 + 0.5 * DT * (
        transport_1 + transport_2 + collision_2
    )
    state_3 = implicit_collision(input_3, 0.5 * DT, epsilon)
    return state_2, state_3


def limiting_stages() -> tuple[float, float, float]:
    p_2 = INITIAL_INVARIANT + DT * (
        centered_operator(INITIAL_INVARIANT)
        + penalty_operator(INITIAL_INVARIANT)
    )
    p_3 = INITIAL_INVARIANT + DT * centered_operator(INITIAL_INVARIANT) \
        + 0.5 * DT * (
            penalty_operator(INITIAL_INVARIANT) + penalty_operator(p_2)
        )
    return p_2, p_3, fick_flux(p_2)


def main() -> None:
    p_2_limit, p_3_limit, flux_3_limit = limiting_stages()
    print("epsilon       |P2-P2lim|       |P3-P3lim|       |I3/eps-F(P2)|")
    for epsilon in (1.0e-1, 3.0e-2, 1.0e-2, 3.0e-3, 1.0e-3, 3.0e-4):
        state_2, state_3 = pdars_step(epsilon)
        p_2 = state_2[0] + state_2[2]
        p_3 = state_3[0] + state_3[2]
        print(
            f"{epsilon:9.1e}  "
            f"{abs(p_2 - p_2_limit):15.7e}  "
            f"{abs(p_3 - p_3_limit):15.7e}  "
            f"{abs(state_3[1] / epsilon - flux_3_limit):15.7e}"
        )

    expected_mismatch = abs(fick_flux(p_2_limit) - fick_flux(p_3_limit))
    print()
    print(f"limiting P2                         = {p_2_limit:.15e}")
    print(f"limiting accepted P                = {p_3_limit:.15e}")
    print(f"|F(P2)-F(Paccepted)|               = {expected_mismatch:.15e}")


if __name__ == "__main__":
    main()
