#include <iostream>
#include <vector>
#include <string>

#include "ariadne/ariadne.hpp"
#include "ariadne/function/function_patch.hpp"
#include "ariadne/symbolic/expression_patch.hpp"

using namespace std;
using namespace Ariadne;
using namespace Ariadne;

ValidatedVectorMultivariateFunctionPatch integrate_dynamics(EffectiveVectorMultivariateFunction f, ExactBoxType x_dom)
{
    auto integrator = TaylorPicardIntegrator(1E-8);

    StepSizeType h = 1 / two;
    auto x = integrator.flow_step(f, x_dom, h);

    return x;
};

ValidatedVectorMultivariateFunctionPatch project_fun(ValidatedVectorMultivariateFunctionPatch f, Range p)
{
    return Array<ValidatedScalarMultivariateFunctionPatch>(p.size(), [&f, &p](SizeType i)
        { return f[p[i]]; });
};

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


int main()
{
    List<RealVectorVariable> x_vars = list_of_vector_variables("x", 2, 2);
    List<RealVectorVariable> u_vars = list_of_vector_variables("u", 2, 1);

    RealVariable x1("x1"), x2("x2"), u("u");
    RealVariable u0("u0"), u1("u1");

    // PC: Function for test system and Lagrangian
    RealExpression x1dot = 10 * x1 + x2;
    RealExpression x2dot = x1 + u;
    RealExpression L = pow(x1, 2) + pow(x2, 2) + pow(u, 2);
    EffectiveVectorMultivariateFunction rhs = Function({ x1, x2, u }, { x1dot, x2dot, L });


    RealScalar time_horizon = 10;
    SizeType num_shooting_nodes = 5;


    Vector<RealScalar> shooting_nodes(num_shooting_nodes + 1);

    // define the shooting nodes
    for (SizeType i = 0;i <= num_shooting_nodes;i++) {
        shooting_nodes[i] = i * time_horizon / num_shooting_nodes;
    }

    // define the time step
    RealScalar dt = time_horizon / num_shooting_nodes;
    StepSizeType h = 1 / two;

    // define function patches for each shooting node
    Vector<FunctionPatch<ValidatedTag, RealVector(RealVector)>> fp(num_shooting_nodes);
    ExactBoxType domain = { {0, 1}, {0, 1}, {0, 1} };
    for (SizeType i = 0;i < num_shooting_nodes;i++) {
        ValidatedVectorMultivariateFunctionPatch phi = integrate_dynamics(rhs, domain);
        ValidatedVectorMultivariateFunctionPatch phi_dynamics = project_fun(phi, Range(0, 2));
        ValidatedVectorMultivariateFunctionPatch phi_lagr = project_fun(phi, Range(2, 3));

        domain = phi_dynamics.domain();

        auto x_dom = project(domain, Range(0, 2));
        BoxDomainType u_dom = project(domain, Range(3, 4));
        fp[i] = FunctionPatch<ValidatedTag, RealVector(RealVector)>({ x_vars[0] | x_dom, u_vars[0] | u_dom }, phi_dynamics({ x_vars[0], u_vars[0], h }));

        domain = product(x_dom, u_dom);
    }

    return 0;
};