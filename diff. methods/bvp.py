import time
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


def f(t, x):
    return [x[1], 6*x[0]**2-t]

class BVP:
    def __init__(self, f, steps, a0, a1, y0=1, t=[0, 1]):
        self.steps = steps
        self.f = f
        self.t_span = t
        self.a = [a0, a1]
        self.y0 = y0
        self.hist_a = [self.a]
    
    def next_a(self, y, ref):
        y = [y[i]-ref for i in range(len(y))]
        
        tmp = self.a[1] - ((self.a[1] - self.a[0]) / (y[1] - y[0])) * y[1]
        self.a = [self.a[1], tmp]

        return self.a
    
    def ivp_solve(self, tmp_a):
        return solve_ivp(self.f, t_span=self.t_span, y0=[self.y0, tmp_a], method="RK45")["y"][0][-1]

    def conv_check(self, y, ref):
        return np.abs(y - ref)

    def solve(self, ref=5):
        for i in range(self.steps):
            ivp_a = self.ivp_solve(self.a[0])
            ivp_b = self.ivp_solve(self.a[1])

            a_next = self.next_a(y=[ivp_a, ivp_b], ref=ref)
            self.hist_a.append(a_next)

            ivp_F = self.ivp_solve(self.a[1])
            print(f'Step {i} with error {self.conv_check(ivp_F, ref=ref)}')
            if self.conv_check(ivp_F, ref=ref) < 1e-1:
                self.a = [float(self.a[i]) for i in range(len(self.a))]
                print(f"Converged at a->{self.a} under {i} steps")
                break

if __name__ == "__main__":
    # main loop

    prob = BVP(f, steps=20, a0=1.2, a1=1.5)

    prob.solve(ref=5)
