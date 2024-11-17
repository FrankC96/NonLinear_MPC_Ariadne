#include <iostream>
#include "ariadne/ariadne.hpp"
#include "symbolic/expression_patch.hpp"
#include "solvers/nonlinear_programming.hpp"

using namespace std;
using namespace Ariadne;

void pprint(string str, auto var)
{
    std::cout << str << " => \n" << var << std::endl;
    for (int i = 0;i < 50;++i) { std::cout << "-"; };
    std::cout << std::endl;
};
ValidatedVectorMultivariateFunctionPatch integrate_dynamics(EffectiveVectorMultivariateFunction f, ExactBoxType x_dom, StepSizeType h)
{
    auto integrator = TaylorPicardIntegrator(1E-2);

    auto x = integrator.flow_step(f, x_dom, h);

    return x;
};
ValidatedVectorMultivariateFunctionPatch project_fun(ValidatedVectorMultivariateFunctionPatch f, Range p)
{
    return Array<ValidatedScalarMultivariateFunctionPatch>(p.size(), [&f, &p](SizeType i) { return f[p[i]]; });
};
List<Variable<Vector<Real>>> list_of_vector_variables(String name, int m, int n)
{
    List<Variable<Vector<Real>>> vars;
    for (int i = 0; i != m; ++i)
    {
        vars.append(Variable<Vector<Real>>(name + to_string(i), n));
    };

    return vars;
}
List<RealVariable> flatten_vars(List<Variable<Vector<Real>>> x, List<Variable<Vector<Real>>> u)
{
    List<RealVariable> temp_vars = {};
    for (int i = 0;i < x.size();++i)
    {
        for (int j = 0;j < x[i].size(); ++j)
        {
            temp_vars.push_back(x[i][j]);
        }
    }
    for (int i = 0;i < u.size();++i)
    {
        for (int j = 0;j < u[i].size(); ++j)
        {
            temp_vars.push_back(u[i][j]);
        }
    }

    return temp_vars;
}


int main()
{
    RealVectorVariable x("x", 2);
    RealVariable t("t");
    RealVariable u("u"), u0("u0"), u1("u1");
    RealVariable c("c");
    StepSizeType h = 0.0625_x;

    BoxDomainType stateDomain = { { 0.0_x,0.125_x }, { 0.0_x,0.125_x } };
    IntervalDomainType inputDomain = { -0.125_x,+0.125_x };
    IntervalDomainType costDomain = { -0.5_x, 0.5_x };
    IntervalDomainType tDomain = { 0, h };
    BoxDomainType varsDomain = product(stateDomain, inputDomain, costDomain);
    BoxDomainType fullDomain = product(varsDomain, tDomain);

    RealExpression x1_dot = x[1];
    RealExpression x2_dot = u;
    RealExpression L = pow(x[0], 2) + pow(x[1], 2) + pow(u, 2);
    RealExpression u_dot = 0;
    EffectiveVectorMultivariateFunction rhs = Function({ x[0], x[1], u, c }, { x1_dot, x2_dot, u_dot, L });

    ValidatedVectorMultivariateFunctionPatch phi0 = integrate_dynamics(rhs, varsDomain, h);
    ValidatedVectorMultivariateFunctionPatch phi_dynamics0 = project_fun(phi0, Range(0, 2));
    ValidatedScalarMultivariateFunctionPatch phi_lagr0 = project_fun(phi0, Range(3, 4))[0];

    auto xuct_domain = phi_dynamics0.domain();
    stateDomain = project(xuct_domain, Range(0, 2));
    costDomain = xuct_domain[3];
    tDomain = {0, xuct_domain[4].upper_bound()};

    BoxDomainType xut_domain = product(stateDomain, inputDomain, tDomain);
    // 1. restrict the integrated dynamics to the xut domain
    ValidatedVectorMultivariateFunctionPatch psi0 = ValidatedVectorMultivariateRestrictedFunction(phi_dynamics0, xuct_domain);
    // 2. create a symbolic expression patch psi, with symbolic variables x, u, t bounded in xut domain
    ValidatedVectorExpressionPatch phiExp0 = psi0({ x[0], x[1], u, c, t });
    // 3. create a symbolic expression patch for x domain with x bounded in state domain provided by the integrator
    ValidatedVectorExpressionPatch s0 = ValidatedVectorExpressionPatch({ x | stateDomain, c | costDomain }, { x[0], x[1] });
    // 4. create the constraint phi(x, u, c, t) - x
    ValidatedVectorExpressionPatch constExp0 = phiExp0 - s0;
    // // // 5. convert the constraint to a FunctionPatch object
    ValidatedVectorMultivariateFunctionPatch constFP0 = ValidatedVectorMultivariateFunctionPatch({ x | stateDomain, u | inputDomain, c | costDomain, t | tDomain }, constExp0);

    stateDomain = cast_exact_box(constFP0.range());
    varsDomain = product(stateDomain, inputDomain, costDomain);

    ValidatedVectorMultivariateFunctionPatch phi1 = integrate_dynamics(rhs, varsDomain, 2 * h);
    ValidatedVectorMultivariateFunctionPatch phi_dynamics1 = project_fun(phi1, Range(0, 3));
    ValidatedScalarMultivariateFunctionPatch phi_lagr1 = project_fun(phi1, Range(3, 4))[0];
    tDomain = { h, 2 * h };
    xut_domain = product(varsDomain, tDomain);
    // restrict l(x, u, t) to (h, 2h) instead of (0, 2h) the integrator gives.
    phi_lagr1 = restriction(phi_lagr1, xut_domain);

    VectorVariableBoxDomainType state2Domain(x, stateDomain);
    VariableIntervalDomainType input02Domain(u, inputDomain);
    VariableIntervalDomainType input12Domain(u1, inputDomain);
    VariableIntervalDomainType cost2Domain(c, costDomain);
    VariableIntervalDomainType t2Domain(t, tDomain);

    ValidatedVectorMultivariateFunctionPatch psi1 = ValidatedVectorMultivariateFunctionPatch({ state2Domain,input02Domain, input12Domain,cost2Domain,t2Domain }, phi_dynamics1({ phi_dynamics0({ x,u,c,t }), input12Domain,cost2Domain,t2Domain }));
    ValidatedVectorExpressionPatch phiExp1 = psi1({ x[0], x[1], u, u1, c, t });
    ValidatedVectorExpressionPatch s1 = ValidatedVectorExpressionPatch({ state2Domain, cost2Domain }, { x[0], x[1], c });
    ValidatedVectorExpressionPatch constExp1 = phiExp1 - s1;
    ValidatedVectorMultivariateFunctionPatch constFP1 = ValidatedVectorMultivariateFunctionPatch({ state2Domain,input02Domain,input12Domain,cost2Domain,t2Domain }, constExp1);

    stateDomain = cast_exact_box(constFP1.range());

    ValidatedScalarExpressionPatch phi_lagrExp0 = phi_lagr0({ x[0], x[1], u, c, t });
    ValidatedScalarExpressionPatch phi_lagrExp1 = phi_lagr1({ x[0], x[1], u1, c, t });
    ValidatedScalarExpressionPatch runningCost = phi_lagrExp0 + phi_lagrExp1;

    ValidatedScalarMultivariateFunction runningCostFP = ValidatedScalarMultivariateFunctionPatch({ state2Domain,input02Domain,input12Domain,cost2Domain,t2Domain }, runningCost);
    ValidatedVectorMultivariateFunction g = combine(constFP0, constFP1);

    pprint("running cost arg size", runningCostFP.argument_size());
    pprint("constraints g arg size", g.argument_size());
    // // Solve the NLP
    NonlinearInfeasibleInteriorPointOptimiser solver;
    // ExactBoxType D = ExactBoxType{ {0.0_x,inf}, {0.0_x,inf}, {0.0_x,inf}, {0.0_x,inf}, {0.0_x,inf}, {0.0_x,inf} };
    // ExactBoxType C = ExactBoxType{ {0.0_x,inf}, {0.0_x,inf}, {0.0_x,inf}, {0.0_x,inf} };
    // auto uOptimal = solver.minimise(runningCostFP, D, g, C);

    return 0;
};