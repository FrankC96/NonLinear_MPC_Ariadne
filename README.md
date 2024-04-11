# NonLinear_MPC_Ariadne

Design a nonlinear controller using Ariadne for a tanker ship steering model, derived by identification methods by Norrbin.

--compiler GCC 10.5

g++ driver.cpp -I/usr/local/include/ariadne -L/usr/local/lib -std=c++20 -lariadne -w -fcompare-debug-second

**point -I and -L to the headers and binaries respectively for Ariadne**
