import time
import numpy as np
from typing import AnyStr, Callable

try:
    from tqdm import trange
except ImportError as e:
    trange = range
    print(f"[WARNING] No progress bar will be displayed, install tqdm")

from plotter import Plotter


class Simulator(Plotter):
    def __init__(
        self,
        controller,
        x_0: np.array,
        x_r: np.array,
        sim_steps: int,
        plant: Callable[[np.array, np.array, float], np.array],
        input_bounds: np.array = None,
    ):
        assert len(x_0) == len(x_r)

        self.controller = controller
        self.x0 = x_0
        self.x_ref = x_r
        self.sim_steps = sim_steps
        self.plant = plant
        self.input_bounds = input_bounds
        if isinstance(self.input_bounds, np.ndarray):
            self.u_max = np.max(self.input_bounds)
            self.u_min = np.min(self.input_bounds)

        self.u0 = controller.optimize(self.x0, self.x_ref)["u"][0]

        self.input_plotter = Plotter(
            n_states=controller.n_inputs,
            title=f"{controller.constr_method} input",
            save_path=f"{controller.constr_method}_input.mp4",
        )
        self.state_plotter = Plotter(
            n_states=controller.n_states,
            title=f"{controller.constr_method} state",
            save_path=f"{controller.constr_method}_state.mp4",
        )

    def get_orbit(self):
        def saturate(u: float):
            # FIXME: provide input bounds for each input variable and not a generic one
            # meaning the user can provide min, max bounds for input variable u[0], u[1]
            # like np.array([[-2, 2], [-5, 5]])
            for input in u:
                if input > self.u_max:
                    u[int(np.where(u == input)[0])] = self.u_max
                elif input < self.u_min:
                    u[int(np.where(u == input)[0])] = self.u_min
            return u

        # Start recording frames at each timestep for plotting
        self.input_plotter.start_recording(fps=15)
        self.state_plotter.start_recording(fps=15)

        # Initialize states [x] and inputs [u] matrices with np.nan values
        x, u = np.full([self.sim_steps + 1, self.controller.n_states], np.nan), np.full(
            [self.sim_steps + 1, self.controller.n_inputs], np.nan
        )

        # Populate the initial state, input in x, u at index [0]
        x[0] = self.x0
        if isinstance(self.input_bounds, np.ndarray):
            u[0] = saturate(self.controller.optimize(x[0], self.x_ref)["u"][0])
        else:
            self.controller.optimize(x[0], self.x_ref)["u"][0]

        running_time = 0.0
        for k in trange(self.sim_steps):
            start = time.perf_counter()

            self.input_plotter.update(k, u[k, :])
            self.state_plotter.update(k, x[k, :])

            x[k + 1] = self.plant(x[k], u[k])
            if isinstance(self.input_bounds, np.ndarray):
                u[k + 1] = saturate(self.controller.optimize(x[k], self.x_ref)["u"][0])
            else:
                u[k + 1] = self.controller.optimize(x[k], self.x_ref)["u"][0]

            end = time.perf_counter()
            running_time += end - start

        # Free memory of writer object in Plotter class
        self.input_plotter.stop_recording()
        self.state_plotter.stop_recording()

        return x, u, running_time
