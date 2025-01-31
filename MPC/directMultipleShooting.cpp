#include <iostream>
#include "ariadne.hpp"
#include "utils.hpp"
#include "ariadne/solvers/runge_kutta_integrator.hpp"

using namespace Ariadne;
using VectorAndScalarPatch = Pair<ValidatedVectorMultivariateFunctionPatch, ValidatedScalarMultivariateFunctionPatch>;

const int num_shooting_nodes = 2;

RealVectorVariable x("x", 2);
RealVariable u("u");
RealVariable t("t");
RealVariable c("c");
StepSizeType h = 0.0625_x;

BoxDomainType xdom = { { -2.0_x, +1.0_x }, { -2.0_x, +1.0_x } };
IntervalDomainType udom = { -0.125_x, +0.125_x };
IntervalDomainType cdom = { -0.0_x, +0.0_x };
IntervalDomainType tdom = { 0, h};
BoxDomainType x_sim_dom = xdom;
IntervalDomainType c_sim_dom = cdom;

BoxDomainType init_doms = product(xdom, udom, cdom);

Vector<RealExpression> x_dot = { u - 2*x[0] - 1*x[1], x[0]};
RealExpression u_dot = 0;
RealExpression c_dot = pow(x[0],2) + pow(x[1], 2) + pow(u, 2);

EffectiveVectorMultivariateFunction f_xuc = Function({ x[0], x[1], u, c }, { x_dot[0], x_dot[1], u_dot, c_dot });

List<RealVectorVariable> X = {};
List<RealVariable> U;

List<BoxDomainType> x_doms = {};
List<IntervalDomainType> u_doms = {};
List<IntervalDomainType> c_doms = {};

List<VariableIntervalOrBoxDomainType> xu_restr_domains;

List<ValidatedVectorMultivariateFunctionPatch> phi_patches;
List<ValidatedScalarMultivariateFunctionPatch> gamma_patches;

List<ValidatedScalarMultivariateFunctionPatch> constraints;

FloatDPApproximationVector optimization_step(BoxDomainType xdom, IntervalDomainType udom, IntervalDomainType cdom)
{   
    X = {};
    U = {};

    x_doms = {};
    u_doms = {};
    c_doms = {};

    xu_restr_domains = {};

    phi_patches = {};
    gamma_patches = {};
    constraints = {};

    ValidatedVectorMultivariateFunctionPatch step;
    ValidatedVectorMultivariateFunctionPatch state;
    ValidatedVectorMultivariateFunctionPatch g;
    ValidatedScalarMultivariateFunctionPatch objective;

    for(int i=0;i<=num_shooting_nodes;++i)
    {
        RealVectorVariable x = RealVectorVariable("x" + to_string(i), 2);
        RealVariable u = RealVariable("u" + to_string(i));
        RealVariable t = RealVariable("t");

        X.append(x);
        U.append(u);
    }
    x_doms.append(xdom);
    u_doms.append(udom);
    c_doms.append(cdom);
    for(int i=0;i<num_shooting_nodes;++i)
    {
        // std::cout << "Creating " << i << "-th shooting node" << "\n"; 

        VectorVariableBoxDomainType xr = x | x_doms[i];
        VariableIntervalDomainType ur = u | u_doms[i];
        VariableIntervalDomainType tr = t | tdom;

        ValidatedVectorMultivariateFunctionPatch phigamma_ = integrateDynamics(f_xuc, product(x_doms[i], udom, c_doms[i]), h);
        ValidatedVectorMultivariateFunctionPatch phigamma = make_function_patch({ xr,ur,tr }, phigamma_, { x[0],x[1],u,0,t });
        ValidatedVectorMultivariateFunctionPatch phi = project_function(phigamma, Range(0, x.size()));
        ValidatedScalarMultivariateFunctionPatch gamma = phigamma[x.size()+1];
        ValidatedVectorMultivariateFunctionPatch psi = make_function_patch({ xr,ur }, phi, { x[0],x[1],u,h });

        xdom = cast_exact_box(psi.range());
        cdom = cast_exact_interval(gamma.range());

        x_doms.append(xdom);
        u_doms.append(udom);
        c_doms.append(cdom);

        phi_patches.append(phi);
        gamma_patches.append(gamma);
    }

    assert(X.size() == U.size());
    assert(x_doms.size() == u_doms.size());
    assert(X.size() == x_doms.size());
    for(int i=0;i<X.size();++i) { xu_restr_domains.append(X[i] | x_doms[i]); }
    for(int i=0;i<U.size();++i) { xu_restr_domains.append(U[i] | u_doms[i]); }

    step = make_function_patch(xu_restr_domains, phi_patches[0], { X[0][0],X[0][1],U[0],h });
    state = make_function_patch(xu_restr_domains, { X[1][0],X[1][1] });
    ValidatedVectorMultivariateFunctionPatch constraint = step - state;
    for(int j=0;j<constraint.result_size();++j) { constraints.append(constraint[j] );}

    for(int i=0;i<num_shooting_nodes-1;++i)
    {
        step = make_function_patch(xu_restr_domains, phi_patches[i], { X[i][0],X[i][1],U[i],h });
        state = make_function_patch(xu_restr_domains, { X[i+1][0],X[i+1][1] });
        ValidatedVectorMultivariateFunctionPatch constraint = step - state;

        for(int j=0;j<constraint.result_size();++j) { constraints.append(constraint[j] );}
    }
    // std::cout << "Appended " << constraints.size() << " constraints." << "\n"; 

    g = ValidatedVectorMultivariateFunctionPatch(constraints);

    objective = make_function_patch(xu_restr_domains, gamma_patches[0], { X[0][0],X[0][1],U[0],h });
    // std::cout << "Initialized cost" << "\n"; 
    for(int j=1;j<num_shooting_nodes;++j)
    {
        objective += make_function_patch(xu_restr_domains, gamma_patches[j], { X[j][0],X[j][1],U[j],h });
        // std::cout << "Appended " << j << "-th objective." << "\n"; 
    }

    NonlinearInfeasibleInteriorPointOptimiser nlio;
        
    int num_decision_state_vars = (num_shooting_nodes+1) * 2;
    int num_decision_input_vars = (num_shooting_nodes+1) * 1;

    int num_decision_vars = num_decision_state_vars + num_decision_input_vars;
    int num_constraints = (num_shooting_nodes) * 2;

    List<IntervalDomainType> D_temp;
    List<IntervalDomainType> C_temp;
    for(int i=0;i<num_constraints;++i) {C_temp.append({-0.125_x, +0.125_x});}

    assert(objective.result_size() == 1);
    assert(objective.domain() == g.domain());
    auto D = Box<ExactIntervalType>(objective.domain());
    auto C = Box<ExactIntervalType>(Vector<IntervalDomainType>(C_temp));

    return nlio.minimise(objective, D, g, C);
}

int main()
{   
    for(int i=0;i<20;++i)
    {   
        std::cout << "Simulating step ---------------\t[" << i << "]" << "\n";
        
        PRINT(x_sim_dom);
        RealVectorVariable x = RealVectorVariable("x_sim" + to_string(i), 2);
        RealVariable u = RealVariable("u_sim" + to_string(i));

        FloatDPApproximationVector vars_opt = optimization_step(x_sim_dom, udom, c_sim_dom);    
        udom = IntervalDomainType({-abs(vars_opt[2*(1+num_shooting_nodes)]), +abs(vars_opt[2*(1+num_shooting_nodes)])});  
        PRINT(vars_opt);
        PRINT(udom);
        List<FloatDPApproximation> x_opt; for (int i=0;i<2*(1+num_shooting_nodes);++i) { x_opt.append(vars_opt[i]); }            
        List<FloatDPApproximation> u_opt; for (int i=2*(1+num_shooting_nodes);i<vars_opt.size();++i) { u_opt.append(vars_opt[i]); }            

        VectorVariableBoxDomainType xr = x | x_sim_dom;
        VariableIntervalDomainType ur = u | udom;
        VariableIntervalDomainType tr = t | tdom;

        ValidatedVectorMultivariateFunctionPatch PHI_GAMMA_ = integrateDynamics(f_xuc, product(x_sim_dom, udom, cdom), h);
        ValidatedVectorMultivariateFunctionPatch PHI_GAMMA = make_function_patch({ xr,ur,tr }, PHI_GAMMA_, { x[0],x[1],u,0,t });
        ValidatedVectorMultivariateFunctionPatch PHI = project_function(PHI_GAMMA, Range(0, x.size()));
        ValidatedScalarMultivariateFunctionPatch GAMMA = PHI_GAMMA[x.size()+1];
        ValidatedVectorMultivariateFunctionPatch PSI = make_function_patch({ xr,ur }, PHI, { x[0],x[1],u,h });

        x_sim_dom = cast_exact_box(PSI.range());
        c_sim_dom = cast_exact_interval(GAMMA.range());
    }
    return 0;
};