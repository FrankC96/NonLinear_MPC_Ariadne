import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

def sim_f(x0: np.array, u: np.array, dt: float, f=None) -> list[list, list]:
    # f callable for an arb system.
    x0_sim, x1_sim = [x0[0]], [x0[1]]
    for i in range(len(u)-1):
        x0 = f(x0, u[i], dt)
        x0_sim.append(x0[0])
        x1_sim.append(x0[1])

    return [x0_sim, x1_sim]

def f(state: np.array, u: float, dt: float=0.0) -> np.array:
    """
    We will model the block as a unit point mass that
    slides without friction in one dimension. The state of the block is its position x and velocity ν, and the
    control is the force u applied to the block.
    """
    x, v = state

    dxdt = v
    dvdt = u

    if dt == 0.0:
        return np.array([dxdt, dvdt])
    else:
        # 2 differrent functions for different functionality
        return state + np.array([dxdt, dvdt])*dt

def constraints(X):
    t = np.linspace(t0, tf, N)

    # decision variables
    x = X[0:len(t)]
    v = X[len(t):2*len(t)]
    u = X[2*len(t):3*len(t)]

    assert len(x) == len(v) == len(u)

    # collocation consts
    ceq = []
    for i in range(N-1):
        h = t[i+1] - t[i]

        # Dimensional analysis
        # m - m - s * (m/s)
        # (m/s) - (m/s) - s * (kg*m/s^2), dividing s with kg to ensure dimensional consistency.
        ceq.append(x[i+1] - x[i] - (h) * ((v[i+1] + v[i]))/2)
        ceq.append(v[i+1] - v[i] - (h) * ((u[i+1] + u[i]))/2)

    ceq.append(x[0] - x0_init[0])
    ceq.append(v[0] - x0_init[1])
    ceq.append(x[-1] - 1)  # meyer term for reaching x=1m at t=1s
    ceq.append(v[-1] - 0)
    # ceq.append(u[-1] - 0) # meyer term for forcing u[end]=0, terrible unfeasible switching
    return ceq

def obj(X: np.array) -> float:
    x = X[0:len(t)]
    v = X[len(t):2*len(t)]
    u = X[2*len(t):3*len(t)]
    return sum([(t[i+1] - t[i])/2 * (u[i]**2 + u[i+1]**2) for i in range(len(u)-1)])

t0, tf = 0, 1
N = 100
dt = (tf - t0) / N

t = np.linspace(t0, tf, N)
# different initization for the problem, 0's also seem to bring good results
x_init = t
v_init = np.ones([len(t)])
u_init = np.zeros([len(t)])

init_guess = np.vstack([x_init, v_init, u_init]).flatten()

global x0_init  # i know, but for some reason scipy doesn't provide extra arguments for the constraints function, only for the objective.
x0_init = np.array([0, 0])

result = minimize(obj, init_guess, constraints={"type": "eq", "fun": constraints})
result_flag = result.success
result_debug = result.message

x_opt, v_opt, u_opt = result.x[0:len(t)], result.x[len(t):2*len(t)], result.x[2*len(t):3*len(t)]

if result_flag:
    print(f'Found solution: {result.success}')
else:
    print(f'Minimization was False, solver exited with\n->{result_debug}\n->{sum(constraints(result.x))}')


x0_sim, x1_sim = sim_f(x0_init, u_opt, dt)

plt.figure(1)
plt.plot(t, x0_sim, '-b', label='x_simulated')
plt.plot(t, x1_sim, '-r', label='v_simulated')
plt.plot(t, x_opt, '--b', label='x')
plt.plot(t, v_opt, '--r', label='v')
plt.legend()
plt.grid("on")
plt.title(rf"State evolution over time {t0, tf} with {N} datapoints")
plt.xlabel("k")
plt.ylabel("x, v")

plt.figure(2)
plt.plot(t, u_opt, '-k', label='input')
plt.legend()
plt.grid("on")
plt.title(rf"Input over discrete time {t0, tf} with {N} datapoints")
plt.xlabel("k")
plt.ylabel("u")

plt.show()

