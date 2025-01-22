#include <iostream>
#include <ariadne.hpp>
#include <symbolic/expression.hpp>
#include <solvers/nonlinear_programming.hpp>
#include <function/function_patch.hpp>

#include "utils.hpp"

using namespace std;
using namespace Ariadne;


int main()
{
    const int num_shooting_nodes = 2;
    // FIXME: system class
    //  === VARIABLES ===
    RealVectorVariable x("x", 2);
    RealVariable u("u");
    RealVariable t("t");
    RealVariable c("c");
    List<RealVariable> system_vars = {u, t, c};

    // === GOV. EQUATIONS ===
    Vector<RealExpression> x_dot = { 2+x[1],u };
    RealExpression c_dot = pow(x[0], 2) + pow(x[1], 2) + pow(u, 2);
    RealExpression u_dot = 0;
    // === RHS ===
    EffectiveVectorMultivariateFunction f_xuc = Function({ x[0], x[1], u, c }, { x_dot[0], x_dot[1], u_dot, c_dot });

    // === STEP SIZE ===
    StepSizeType h = 0.0625_x;
    // === INIT. BOX CONDS. (x0, u0, t0, c0) ===
    BoxDomainType xdom = { { -1.0_x, +1.0_x }, { -1.0_x, +1.0_x } };
    IntervalDomainType udom = { -0.125_x, +0.125_x };
    IntervalDomainType tdom = { 0, h };
    IntervalDomainType cdom = { -0.0_x, 0.0_x };
    Tuple<BoxDomainType, IntervalDomainType, IntervalDomainType, IntervalDomainType> restr_doms = make_tuple(xdom, udom, cdom, tdom);
    Tuple<
        BoxDomainType,
        IntervalDomainType,
        ValidatedVectorMultivariateFunctionPatch,
        ValidatedScalarMultivariateFunctionPatch
        > step_result;
    
    ValidatedVectorMultivariateFunctionPatch phi;
    ValidatedScalarMultivariateFunctionPatch gamma;
    List<RealVectorVariable> state_variables = {x};
    List<RealVariable> input_variables = {u};
    List<RealVariable> time_variables = {t};
    List<ValidatedVectorMultivariateFunctionPatch> phi_patches;
    List<ValidatedScalarMultivariateFunctionPatch> gamma_patches;
    List<BoxDomainType> x_doms;
    List<IntervalDomainType> u_doms;
    List<VariableIntervalOrBoxDomainType> xu_doms_restr;

    for(int iter=1; iter<=num_shooting_nodes; ++iter)
    {
        x_doms.append(xdom);
        u_doms.append(udom);

        x = RealVectorVariable("x" + to_string(iter), x.size());
        u = RealVariable("u" + to_string(iter));
        system_vars = {u, t, c};

        step_result = makeStep(iter, f_xuc, x, system_vars, restr_doms, h);
        
        phi = std::get<2>(step_result);
        gamma = std::get<3>(step_result);

        xdom = std::get<0>(step_result);
        tdom = std::get<1>(step_result);

        restr_doms = make_tuple(xdom, udom, cdom, tdom);

        phi_patches.append(phi);
        gamma_patches.append(gamma);
        x_doms.append(xdom);
        u_doms.append(udom);

        state_variables.append(x);
        input_variables.append(u);
        time_variables.append(t);

        for(int i=0;i<x.size();++i) { xu_doms_restr.append(state_variables[i] | x_doms[i]); }
        // FIXME: generalize 
        for(int i=0;i<num_shooting_nodes;++i) { xu_doms_restr.append(input_variables[i] | u_doms[i]); }
        
        // if(iter!=num_shooting_nodes) { xu_doms_restr.clear(); }
        
        ValidatedVectorMultivariateFunctionPatch step0 = make_function_patch(xu_doms_restr, phi_patches[iter-1], { state_variables[iter-1][0],state_variables[iter-1][1],input_variables[iter-1],h });
        ValidatedVectorMultivariateFunctionPatch step1 = make_function_patch(xu_doms_restr, phi_patches[iter], { state_variables[iter][0],state_variables[iter][1],input_variables[iter],h });
        ValidatedVectorMultivariateFunctionPatch state1 = make_function_patch(xu_doms_restr, { state_variables[iter][0],state_variables[iter][1] });

        ValidatedScalarMultivariateFunctionPatch objective0 = make_function_patch(xu_doms_restr, gamma_patches[iter-1], { state_variables[iter-1][0],state_variables[iter-1][1],input_variables[iter-1],h });
        ValidatedScalarMultivariateFunctionPatch objective1 = make_function_patch(xu_doms_restr, gamma_patches[iter], { state_variables[iter][0],state_variables[iter][1],input_variables[iter],h });
        ValidatedScalarMultivariateFunctionPatch objective = objective0 + objective1;

        ValidatedVectorMultivariateFunctionPatch constraint1 = step0 - state1;
        PRINT(constraint1);
    };

    // for(int iter=1;iter<2;++iter)
    // {
    //     ValidatedVectorMultivariateFunctionPatch step0 = make_function_patch(xu_doms_restr, phi_patches[iter-1], { state_variables[iter-1][0],state_variables[iter-1][1],input_variables[iter-1],h });
    //     ValidatedVectorMultivariateFunctionPatch step1 = make_function_patch(xu_doms_restr, phi_patches[iter], { state_variables[iter][0],state_variables[iter][1],input_variables[iter],h });
    //     ValidatedVectorMultivariateFunctionPatch state1 = make_function_patch(xu_doms_restr, { state_variables[iter][0],state_variables[iter][1] });

    //     ValidatedScalarMultivariateFunctionPatch objective0 = make_function_patch(xu_doms_restr, gamma_patches[iter-1], { state_variables[iter-1][0],state_variables[iter-1][1],input_variables[iter-1],h });
    //     ValidatedScalarMultivariateFunctionPatch objective1 = make_function_patch(xu_doms_restr, gamma_patches[iter], { state_variables[iter][0],state_variables[iter][1],input_variables[iter],h });
    //     ValidatedScalarMultivariateFunctionPatch objective = objective0 + objective1;

    //     ValidatedVectorMultivariateFunctionPatch constraint1 = step0 - state1;
    //     PRINT(constraint1);
    // }
    return 0;
};