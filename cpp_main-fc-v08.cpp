#include <iostream>
#include <vector>
#include <string>

#include "ariadne/ariadne.hpp"
#include "ariadne/function/function.hpp"
#include "ariadne/function/function_patch.hpp"
#include "ariadne/symbolic/expression_patch.hpp"
#include "ariadne/solvers/nonlinear_programming.hpp"

using namespace std;
using namespace Ariadne;

class MPC {

public:
    EffectiveVectorMultivariateFunction system;
    MPC(List<RealVariable> x, List<RealVariable> u, Vector<RealExpression> odes)
    {
        // unpack all state and input variables to a list
        List<RealVariable> tmp_vars = {};
        for (RealVariable c : x)
        {
            tmp_vars.push_back(c);
        }
        for (RealVariable c : u)
        {
            tmp_vars.push_back(c);
        }

        // issue 1 - FIXED: from documentation, odes object should have been a vector, not a list
        // EffectiveVectorMultivariateFunction system = Function(tmp_vars, odes);
    };

    ValidatedVectorMultivariateFunctionPatch integrate_dynamics(EffectiveVectorMultivariateFunction f, ExactBoxType x_dom, StepSizeType h)
    {
        auto integrator = TaylorPicardIntegrator(1E-2);

        // StepSizeType h = 1 / (two / 2);
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

};


int main()
{

    RealVariable x1("x1"), x2("x2"), u("u");

    // PC: Function for test system and Lagrangian
    RealExpression x1dot = x2;
    RealExpression x2dot = u;
    RealExpression L = pow(x1, 2);

    List<RealVariable> X = { x1, x2 };
    List<RealVariable> U = { u };
    Vector<RealExpression> system = { x1dot, x2dot, L };

    MPC mpc = MPC(X, U, system);

    EffectiveVectorMultivariateFunction rhs = Function({ x1, x2, u }, { x1dot, x2dot, L });

    List<RealVectorVariable> x_vars = mpc.list_of_vector_variables("x", 2, 2);
    List<RealVectorVariable> u_vars = mpc.list_of_vector_variables("u", 2, 1);

    RealScalar time_horizon = 10;
    SizeType num_shooting_nodes = 5;


    Vector<RealScalar> shooting_nodes(num_shooting_nodes + 1);

    // FC: define the shooting nodes
    for (SizeType i = 0;i <= num_shooting_nodes;i++) {
        shooting_nodes[i] = i * time_horizon / num_shooting_nodes;
    }

    // FC: define the time step
    RealScalar dt = time_horizon / num_shooting_nodes;
    StepSizeType h = 1 / (2 * two);

    List<RealVariable> flat_vars = mpc.flatten_vars(x_vars, u_vars);
    ExactBoxType cost_func_domain = { {0, 1}, {0, 1}, {0, 1}, {0, 1} };

    FunctionFactory<ValidatedTag> factory();
    TaylorFunctionFactory function_factory(ThresholdSweeper<FloatDP>(dp, 1e-8));
    ScalarMultivariateFunction<ValidatedTag> f = function_factory.create_zero(cost_func_domain);

    List<FunctionPatch<ValidatedTag, RealVector(RealVector)>> fp(num_shooting_nodes);
    ExactBoxType domain = { {0, 1}, {0, 1}, {0, 1} };
    for (SizeType k = 0;k < num_shooting_nodes;k++) {
        ValidatedVectorMultivariateFunctionPatch phi = mpc.integrate_dynamics(rhs, domain, h);
        ValidatedVectorMultivariateFunctionPatch phi_dynamics = mpc.project_fun(phi, Range(0, 2));
        ValidatedScalarMultivariateFunctionPatch phi_lagr = mpc.project_fun(phi, Range(2, 3))[0];

        domain = phi_dynamics.domain();

        auto x_dom = project(domain, Range(0, 2));
        BoxDomainType u_dom = project(domain, Range(2, 3));

        fp[k] = FunctionPatch<ValidatedTag, RealVector(RealVector)>({ x_vars[0] | x_dom, u_vars[0] | u_dom }, phi_dynamics({ x_vars[0], u_vars[0], h }));

        f = f + phi_lagr;


        x_dom = cast_exact_box(fp[k].range());
        domain = product(x_dom, u_dom);
    }

    VectorMultivariateFunction<ValidatedTag> g = join(fp[0], fp[1]);

    cout << 'f ' << f.argument_size() << 'g ' << g.argument_size() << endl;

    // Solve the NLP
    NonlinearInfeasibleInteriorPointOptimiser solver;
    ExactBoxType D = ExactBoxType{ {-1.0_x,2.0_x}, {-3.0_x,5.0_x}, {-3.0_x,5.0_x} };
    ExactBoxType C = ExactBoxType{ {0.0_x,inf}, {0.0_x,inf}, {0.0_x,inf} , {0.0_x,inf} };
    FloatDPBoundsVector x_optimal = solver.minimise(f, D, g, C);

    return 0;
};

