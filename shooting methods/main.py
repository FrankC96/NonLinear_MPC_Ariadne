import numpy as np

from controller import Controller
from simulator import Simulator
from functools import partial


def robot_model(x: np.array, u: np.array) -> np.array:
    """
    A mathematical model representing the plant, by a set of 1st order nonlinear equations.

    The nonlinear system is normalized by
        -> [L] distance between each wheel
        -> [V] maximum linear velocity of the robot
    """
    L, V = 0.5, 5.0

    T = L / V

    u0_norm = u[0] / V
    u1_norm = u[1] * T

    dxdt = np.cos(x[2]) * u0_norm
    dydt = np.sin(x[2]) * u0_norm
    dthetadt = u1_norm

    dxdt_norm = dxdt / L
    dydt_norm = dydt / L

    return np.array([dxdt_norm, dydt_norm, dthetadt])


def robot_plant(x: np.array, u: np.array, mag: float) -> np.array:
    """
    Adding random uniform noise to the dynamics to approximate the mismatch
    between the dynamics and the actual plant.
    """
    return robot_model(x, u) + np.random.uniform(-mag, mag, (3,))


if __name__ == "__main__":
    x_init = np.array([5.5, 5.5, 0.0])
    x_ref = np.array([0.0, 0.0, 0.0])

    T_MAX = 5
    N_STATES = len(x_init)
    N_INPUTS = 2

    PLANT_NOISE = 0.2
    # FIXME: put something reasonable
    STATE_BOUNDS = np.array([-1, 1])
    INPUT_BOUNDS = np.array([-1, 1])

    Q = np.eye(N_STATES)
    R = np.eye(N_INPUTS)

    # using the dynamics function to perform a prediction step
    mpc_coll = Controller(
        constr_method="COLL",
        model=robot_model,
        n_states=N_STATES,
        n_inputs=N_INPUTS,
        n_pred=20,
        t_max=T_MAX,
        Q=Q,
        R=R,
        state_bounds=STATE_BOUNDS,
        input_bounds=INPUT_BOUNDS,
        minimize_method="SLSQP",
        term_constr=False,
    )
    mpc_dms = Controller(
        constr_method="DMS",
        model=robot_model,
        n_states=N_STATES,
        n_inputs=N_INPUTS,
        n_pred=20,
        t_max=T_MAX,
        Q=Q,
        R=R,
        state_bounds=STATE_BOUNDS,
        input_bounds=INPUT_BOUNDS,
        minimize_method="SLSQP",
        term_constr=False,
    )

    # using the plant function to perform a simulation step
    sim_coll = Simulator(
        controller=mpc_coll,
        x_0=x_init,
        x_r=x_ref,
        sim_steps=20,
        plant=partial(robot_plant, mag=PLANT_NOISE),
        input_bounds=None,
    )

    sim_dms = Simulator(
        controller=mpc_dms,
        x_0=x_init,
        x_r=x_ref,
        sim_steps=20,
        plant=partial(robot_plant, mag=PLANT_NOISE),
        input_bounds=None,
    )

    x_coll, u_coll, t_coll = sim_coll.get_orbit()
    print(f"[{mpc_coll.constr_method}] Final x, u {x_coll[-1, :]} | {u_coll[-1, :]}")

    x_dms, u_dms, t_dms = sim_dms.get_orbit()
    print(f"[{mpc_coll.constr_method}] Final x, u {x_dms[-1, :]} | {u_dms[-1, :]}")
