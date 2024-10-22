### notes on trap. collocation

1. `sim_f` is only for visualization purposes, just simulating the system with the computed vector `u`
2. my `f` in case of `dt=0` returns the righ-hand-side and if `dt` is provided, returns an Euler step in for simulating the system.
    - Is this bad practice in terms each function should serve a single purpose?
3. `v[i+1] - v[i] - (h/2) * (u[i+1] + u[i])` -> `(m/s) - (m/s) - (s) * (m/s^2)` i see it now, I got confused with `m=1` mass and meters...
