#include <iostream>
#include <vector>
#include <string>

#include "ariadne/ariadne.hpp"
#include "ariadne/function/function_patch.hpp"
#include "ariadne/symbolic/expression_patch.hpp"

using namespace std;
using namespace Ariadne;

// PC: Changed matrix variables to a list of vector variables.
List<Variable<Vector<Real>>> list_of_vector_variables(String name, int m, int n)
{
    List<Variable<Vector<Real>>> vars;
    for (int i = 0; i != m; ++i)
    {
        vars.append(Variable<Vector<Real>>(name + to_string(i), n));
        // for (int j=0; j!=n; ++j) {
        //     vars.push_back(RealVariable(name + to_string(i) + to_string(j)));
        // }
    };

    return vars;
};

// PC: It's useful as documentation to put the return type in the function signature rather than use auto.
ValidatedVectorMultivariateFunctionPatch integrate_dynamics(EffectiveVectorMultivariateFunction f, ExactBoxType x_dom, StepSizeType h)
{
    auto integrator = TaylorPicardIntegrator(1E-8);
    /* FC: I cannot understand the reason we need to recompute the integrand (phi),
            if we want to evaluate the integrand in different time ranges for each
            shooting node in the domain [0, N], then the TaylorPicardIntegrator
            isn't suitable for our problem because because it integrates in the
            range [0, h] every time, we need [0, h] -> [h, 2*h].

            N: prediction horizon
            h: time step for the integrator
    */
    auto x = integrator.flow_step(f, x_dom, h);

    return x;
};

ValidatedVectorMultivariateFunctionPatch project_fun(ValidatedVectorMultivariateFunctionPatch f, Range p)
{
    return Array<ValidatedScalarMultivariateFunctionPatch>(p.size(), [&f, &p](SizeType i)
        { return f[p[i]]; });
}

int main()
{

    // PC: Prefer const values to macros for constants in C++
    const int pred_hor = 2;
    const int nx = 2;
    const int nu = 1;

    RealVariable x1("x1"), x2("x2"), u("u");
    RealVariable u0("u0"), u1("u1");

    // PC: Function for test system and Lagrangian
    RealExpression x1dot = 10 * x1 + x2;
    RealExpression x2dot = x1 + u;
    RealExpression L = pow(x1, 2) + pow(x2, 2) + pow(u, 2);
    EffectiveVectorMultivariateFunction rhs = Function({ x1, x2, u }, { x1dot, x2dot, L });

    // integrated dynamics, this is to replace F in ln. 57
    // StepSizeType h = 1 / two;
    // ValidatedVectorMultivariateFunctionPatch phi = integrate_dynamics(rhs, x_dom, h);

    // ValidatedVectorMultivariateFunctionPatch phi_dynamics = project_fun(phi, Range(0, 2));
    // ValidatedVectorMultivariateFunctionPatch phi_lagr = project_fun(phi, Range(2, 3));

    // // PC: I think it's easier to make x and u variables separately
    // RealVectorVariable x_var = RealVectorVariable("x", nx);
    // RealVectorVariable u_var = RealVectorVariable("u", nu);
    // // PC: I think it's easier to make x and u variables separately
    List<RealVectorVariable> x_vars = list_of_vector_variables("x", pred_hor, nx);
    List<RealVectorVariable> u_vars = list_of_vector_variables("u", pred_hor, nu);

    // SizeType x_size = 2;
    // SizeType t_size = 1;
    // SizeType u_size = 1;

    // auto x_range = Range(0, x_size);
    // auto t_index = x_size + u_size;
    // auto u_range = Range(x_size + 1u, x_size + 1 + u_size);

    // // We have a function [y,c]=phi(x,t,u).
    // // Extract y=psi(x,h,u) for constant h
    // BoxDomainType xut_dom = phi_dynamics.domain();
    // // Dyadic h = Dyadic(xut_dom[t_index].upper_bound());
    // BoxDomainType x_dom = project(xut_dom, x_range);
    // IntervalDomainType t_dom = xut_dom[x_size];
    // BoxDomainType u_dom = project(xut_dom, u_range);

    // // FC: empty constructors for the NLP
    ScalarMultivariateFunction<ValidatedTag> f;
    // VectorMultivariateFunction<ValidatedTag> g;

    // ExpressionPatch<ValidatedTag, RealVector> ep1 = phi_lagr({x_vars[0] | x_dom, u_vars[0] | u_dom, h});
    // ExpressionPatch<ValidatedTag, RealVector> ep2 = phi_lagr({x_vars[1] | x_dom, u_vars[1] | u_dom, h});

    cout << x_vars[0] << endl;
    StepSizeType h = 1 / two;
    ExactBoxType xu_dom({ {0, 1}, {0, 1}, {0, 1} });
    for (int k = 1; k <= pred_hor; k++)
    {

        /*
        DIRECT MULTIPLE SHOOTING

        1.DISCRETIZE CONTROLS ON AN EVENLY(OPT) SPREAD GRID
        2.SOLVE ODE FOR EACH INTERVAL [t_i, t_{i+1}] NUMERICALLY, STARTING WITH AN INITIAL ARBITRARY VALUE
        3.COMPUTE THE INTEGRAL LOSS
        */
        // cout << x_dom << endl;
        ValidatedVectorMultivariateFunctionPatch phi = integrate_dynamics(rhs, xu_dom, h);

        // TO PROVLIMA EINAI STIN project_fun, girnaei vector akoma kai stin scalar, ftiaxto to gamidi na paei na gamithei
        ValidatedVectorMultivariateFunctionPatch phi_dynamics = project_fun(phi, Range(0, 2));
        ValidatedVectorMultivariateFunctionPatch phi_lagr = project_fun(phi, Range(2, 3));

        BoxDomainType xut_dom = phi_dynamics.domain();
        auto x_dom = project(xut_dom, Range(0, 2));  // 

        IntervalDomainType t_dom = xut_dom[2];
        BoxDomainType u_dom = project(xut_dom, Range(3, 4));

        FunctionPatch<ValidatedTag, RealVector(RealVector)> psi =
            FunctionPatch<ValidatedTag, RealVector(RealVector)>({ x_vars[0] | x_dom, u_vars[0] | u_dom }, phi_dynamics({ x_vars[0], u_vars[0], h }));
        // x_dom = cast_exact_box(psi.range());
        // x_dom[3] = {0, 1};
        xu_dom = product(x_dom, u_dom);
        // cout << xu_dom << endl;
        //  ~~~~WANTED FUNC
        // g.set(k - 1,
        //   FunctionPatch<ValidatedTag, RealVector(RealVector)>({x_vars[0] | x_dom, u_vars[0] | u_dom},
        //                                                       phi_dynamics({x_vars[0], u_vars[0], h})));

        // TODO : J += phi_lagr({x_vars[k - 1] | x_dom, u_vars[k - 1] | u_dom, k * h});
        /* operator `+` that is friend to ExpressionPatch class and adds 2 patches together as expressions.*/
        // f = f + phi_lagr({ x_vars[0] | x_dom, u_vars[0] | u_dom, h });

        // TODO : assign all computed psi`s to a structure.
        /*
        //  ~~~~WANTED FUNC
        // g.set(k - 1,
        //       FunctionPatch<ValidatedTag, RealVector(RealVector)>({x_vars[k - 1] | x_dom, u_vars[k - 1] | u_dom},
        //   phi_dynamics({x_vars[k - 1], u_vars[k - 1], h})));

        // cout << FunctionPatch<ValidatedTag, RealVector(RealVector)>({x_vars[0] | x_dom, u_vars[0] | u_dom}, phi_dynamics({x_vars[0], u_vars[0], k * h})) << endl;
        */
    }

    return 0;
};

/*
=========QUESTIONS=========
how can I:
    1. add 2 ExpressionPatch objects?
    2. evaluate the integrand in different sub ranges in [0, N]


The function f : x |-> x^2 and the function g : y |-> y^2

are the same function!

We have a function phi : R^n * R |-> R^n

The function psi (x0) := phi(x0,h) and psi (x1):=phi(x1,h)

are the same function.

You define g0 by
g0(x0) := phi(x0,h)

and g1 by
g1(x1):=phi(x1,2*h)

So here, the variable name x0 or x1 is irrelevant

The definition of phi says:
phi(x,t) gives the solution of the differential equation dx/dt=f(x), starting at x, and flowing for time t.

g0 Is the the function taking x to the state at time h

g1 is the function taking x to the state after time 2*h

The function from decision node k to node k+1 (separated by time h) is always psi, where psi(x,u) := phi(x,t,u)

However, the possible states you can start at at time t0=0, may be different from the states you can start at at time t1=h

If the set of possible starting points at time t_k  = k*h is D_k, then for the kth step, we need to compute psi over the domain D_k.

This means you probably need to compute phi over a different domain.

flow_step(f, D0, h) will compute the same function phi as flow_step(f,D1,h)

but will compute an approximation $\hat{\phi} _0$ over domain $x_init \in D_0$(and $t \in[0, h] $)

    wheras flow_step(f, D1, h)
computes an approximation $\bat{\phi} _1$ over domain $x_init \in D_1$ and $t\in[0, h] $.

    Maybe here is a source of confusion :

    We can define a one -
    step function
        psi(x, u) : = phi(x, h, u); equivalently psi(x0,u0) := phi(x0,h,u0)

Then the two-step function has two inputs, and is given by
psi_{(2)}(x0,u0,u1) = psi(x1,u1) where x1=psi(x0,u0).

This is then given by psi_{(2)}(x0,u0,u1)=psi(psi(x0,u0),u1) = phi(phi(x0,h,u0),u1)

Compute the domain of psi for the second step by
  cast_exact_box(range(psi))

cast_exact_box(psi.range())

Your g is a function of [for two steps] x0,x1,x2,u0,u1

Constraints are psi(x0,u0)=x1 and psi(x1,u1)=x2

Rewrite constraints as psi(x0,u0)-x1=0 or x1-psi(x0,u0)=0

So "pieces" of g will be of this form, and the corresponding C will be 0

[Or just write list of constraints psi(x0,u0)-x1=0; but this will require some extensions to the code.]

Extension: allow vector functions and function patches.

*/