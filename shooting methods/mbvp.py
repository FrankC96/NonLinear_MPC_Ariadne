import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Define the system of ODEs for y1' = y2, y2' = -y1
def system(t, y):
    return [y[1], -y[0]]

# Boundary value problem conditions
x_start, x_end = 0, np.pi / 2
y1_start, y1_end = 0, 1  # Boundary conditions y(0) = 0, y(pi/2) = 1

# Single Shooting: Solve for a guessed initial condition of y2(0)
def single_shooting(guess_y2_start):
    # Initial condition for IVP
    y0 = [y1_start, guess_y2_start]
    # Solve the ODE using solve_ivp from x_start to x_end
    sol = solve_ivp(system, [x_start, x_end], y0, dense_output=True)
    return sol

# Guess for y2(0) - Single Shooting Approach
guess_y2_single = 1.0
sol_single = single_shooting(guess_y2_single)

# Multiple Shooting: Divide the interval into two segments [0, pi/4] and [pi/4, pi/2]
x_mid = np.pi / 4

# First segment [0, pi/4] - Initial guess for y2(0)
guess_y2_multiple = 1.0
y0_segment1 = [y1_start, guess_y2_multiple]
sol_segment1 = solve_ivp(system, [x_start, x_mid], y0_segment1, dense_output=True)

# Use the result of segment 1 as the starting point for segment 2
y0_segment2 = [sol_segment1.y[0, -1], sol_segment1.y[1, -1]]
sol_segment2 = solve_ivp(system, [x_mid, x_end], y0_segment2, dense_output=True)

# Plot the solutions from single and multiple shooting
x_plot = np.linspace(x_start, x_end, 100)
y_single = sol_single.sol(x_plot)[0]

x_plot_segment1 = np.linspace(x_start, x_mid, 50)
x_plot_segment2 = np.linspace(x_mid, x_end, 50)
y_segment1 = sol_segment1.sol(x_plot_segment1)[0]
y_segment2 = sol_segment2.sol(x_plot_segment2)[0]

# Plot
plt.figure(figsize=(10, 6))
plt.plot(x_plot, y_single, label='Single Shooting', color='blue')
plt.plot(x_plot_segment1, y_segment1, label='Multiple Shooting - Segment 1', color='green')
plt.plot(x_plot_segment2, y_segment2, label='Multiple Shooting - Segment 2', color='red')
plt.scatter([x_start, x_end], [y1_start, y1_end], color='black', zorder=5, label="Boundary conditions")
plt.title("Comparison of Single Shooting and Multiple Shooting Methods")
plt.xlabel("x")
plt.ylabel("y(x)")
plt.legend()
plt.grid(True)
plt.show()
