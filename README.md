# NonLinear_MPC_Ariadne

Design a nonlinear controller using Ariadne for a tanker ship steering model, derived by identification methods by Norrbin.

--compiler GCC 10.5
g++-10 -I/usr/local/include/ariadne -fdiagnostics-color=always -g driver.cpp -o driver -std=c++20 -lariadne -w -fcompare-debug-second

**point -I to the headers for Ariadne**