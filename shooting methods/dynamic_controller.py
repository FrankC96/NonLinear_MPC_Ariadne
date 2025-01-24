import numpy as np
from sympy import symbols, cos, sin, Matrix

# Define symbols for states and inputs
x1, x2, x3, u1, u2 = symbols("x1 x2 x3 u1 u2")  # States and input

states = [x1, x2, x3]
inputs = [u1, u2]

f1 = cos(x3) * u1
f2 = sin(x3) * u1
f3 = u2

rhs = [f1, f2, f3]

# Define the nonlinear system
f = Matrix(rhs)

state = Matrix(states)  # State vector
input = Matrix(inputs)  # Input vector

# Compute Jacobian matrices
A = f.jacobian(state)  # Partial derivatives of f w.r.t. state
B = f.jacobian(input)  # Partial derivatives of f w.r.t. input

# Display symbolic Jacobians
print("Jacobian A (w.r.t state):")
print(A)
print("\nJacobian B (w.r.t input):")
print(B)

# Define operating point
operating_point = {x1: 0, x2: 0, x3: 0, u1: 0, u2: 0}

# Evaluate Jacobians at the operating point
A_eval = A.subs(operating_point)
B_eval = B.subs(operating_point)

print("\nEvaluated Jacobian A:")
print(np.array(A_eval.evalf()))
print("\nEvaluated Jacobian B:")
print(np.array(B_eval.evalf()))
