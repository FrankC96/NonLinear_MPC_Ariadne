#include <iostream>
#include "ariadne.hpp"
#include "utils.hpp"

class problemFormulation
{
public:
    int num_states;
    int num_inputs;

    int num_shooting_nodes;
    StepSizeType time_horizon;

    EffectiveVectorMultivariateFunction dynamics_eqs;

    BoxDomainType x0dom;
    IntervalDomainType u0dom;
    IntervalDomainType c0dom;
    IntervalDomainType t0dom;

    List<BoxDomainType> xdoms;


    problemFormulation(int num_states, int num_inputs, int num_shooting_nodes, StepSizeType time_horizon) : 
                                    num_states(num_states), num_inputs(num_inputs), num_shooting_nodes(num_shooting_nodes), time_horizon(time_horizon) {};

    void setDynamics(EffectiveVectorMultivariateFunction f_xuc, ExactBoxType initial_domains)
    {   
        // RHS
        dynamics_eqs = f_xuc;

        x0dom = project(initial_domains, Range(0, num_states));
        u0dom = initial_domains[2];
        c0dom = initial_domains[3];
        t0dom = initial_domains[4];

        xdoms.append(x0dom);
    };
};

class problemOptimization
{
private:
    // make problemOptimization inherit from problemFormulation
    problemFormulation prob;

public:
    List<RealVectorVariable> X;
    List<RealVariable> U;
    List<ExactBoxType> X_DOMS;
    List<ExactBoxType> U_DOMS;
    List<BoxOrIntervalDomain> XU_DOMAIN;
    Tuple<ValidatedVectorMultivariateFunctionPatch, ValidatedScalarMultivariateFunctionPatch> state_result;
    List<ValidatedVectorMultivariateFunctionPatch> phi_patches;
    List<ValidatedScalarMultivariateFunctionPatch> gamma_patches;

    problemOptimization(problemFormulation prob) : prob(prob) {};

    Tuple<ValidatedVectorMultivariateFunctionPatch, ValidatedScalarMultivariateFunctionPatch> calculateNextState(EffectiveVectorMultivariateFunction f, int current_shooting_node)
    {

        auto x = RealVectorVariable("x" + to_string(current_shooting_node), prob.num_states);
        auto u = RealVariable("u" + to_string(current_shooting_node));
        auto t = RealVariable("t");

        VectorVariableBoxDomainType xr = x | prob.xdoms[current_shooting_node];
        VariableIntervalDomainType ur = u | prob.u0dom;
        VariableIntervalDomainType tr = t | prob.t0dom;

        ValidatedVectorMultivariateFunctionPatch phigamma_ = integrateDynamics(prob.dynamics_eqs, product(prob.xdoms[current_shooting_node], prob.u0dom, prob.c0dom), prob.time_horizon);
        ValidatedVectorMultivariateFunctionPatch phigamma = make_function_patch({ xr,ur,tr }, phigamma_, { x[0],x[1],u,0,t });
        ValidatedVectorMultivariateFunctionPatch phi = project_function(phigamma, Range(0, x.size()));
        ValidatedScalarMultivariateFunctionPatch gamma = phigamma[x.size()+1];
        ValidatedVectorMultivariateFunctionPatch psi = make_function_patch({ xr,ur }, phi, { x[0],x[1],u,prob.time_horizon });

        prob.xdoms.append(cast_exact_box(psi.range()));
        return make_tuple(phi, gamma);
    }

    List<ExactBoxType> getDomains() {return prob.xdoms;}

    ValidatedVectorMultivariateFunctionPatch createConstraints(
        List<ValidatedVectorMultivariateFunctionPatch> phi_patches,
        List<RealVectorVariable> full_x,
        List<RealVariable> full_u,
        List<BoxOrIntervalDomain> full_xu_doms)
    {
        List<VariableIntervalOrBoxDomainType> restricted_domains;
        ValidatedVectorMultivariateFunctionPatch step;
        ValidatedVectorMultivariateFunctionPatch state;
        ValidatedVectorMultivariateFunctionPatch g;
        List<ValidatedScalarMultivariateFunctionPatch> constraints;

        
        for(int i=1;i<full_x.size();++i) {restricted_domains.append(full_x[i] | prob.xdoms[i]);} 
        for(int i=1;i<full_u.size();++i) {restricted_domains.append(full_u[i] | prob.u0dom);} 

        std::cout << "\n" << "Constraint vars" << "\n";
        for(int i=1;i<full_x.size();++i) {std::cout<< full_x[i] << "\t";}
        std::cout << "\n";
        for(int i=1;i<full_x.size();++i) {std::cout<< full_u[i] << "\t";}
        std::cout << "\n";
        std::cout << "with xi size " << full_x[1].size() << "\n";
        std::cout << "and ui size " << full_x[1].size() << "\n";
        for(int i=1;i<prob.num_shooting_nodes;++i)
        {   
            step = make_function_patch(restricted_domains, phi_patches[i-1], { full_x[i][0],full_x[i][1],full_u[i],prob.time_horizon });
            state = make_function_patch(restricted_domains, { full_x[i+1][0],full_x[i+1][1] });
            ValidatedVectorMultivariateFunctionPatch constraint = step - state;

            // append each state constraint to a list
            for(int j=0;j<constraint.result_size();++j) { constraints.append(constraint[j] );}
        }
        
        // convert constraint list to ValidatedVectorMultivariateFunctionPatch
        return ValidatedVectorMultivariateFunctionPatch(constraints);
    }

    ValidatedScalarMultivariateFunctionPatch createObjective(
        List<ValidatedScalarMultivariateFunctionPatch> gamma_patches,
        List<RealVectorVariable> full_x,
        List<RealVariable> full_u,
        List<BoxOrIntervalDomain> full_xu_doms)
    {
        List<VariableIntervalOrBoxDomainType> restricted_domains;
        ValidatedScalarMultivariateFunctionPatch objective;

        for(int i=1;i<full_x.size();++i) {restricted_domains.append(full_x[i] | prob.xdoms[i]);} 
        for(int i=1;i<full_u.size();++i) {restricted_domains.append(full_u[i] | prob.u0dom);} 

        // DEBUG:
        std::cout << "\n" << "Objective vars" << "\n";
        for(int i=1;i<full_x.size();++i) {std::cout<< full_x[i] << "\t";}
        std::cout << "\n";
        for(int i=1;i<full_x.size();++i) {std::cout<< full_u[i] << "\t";}
        std::cout << "\n";
        std::cout << "with xi size " << full_x[1].size() << "\n";
        std::cout << "and ui size " << full_x[1].size() << "\n \n";

        // initial cost
        objective = make_function_patch(restricted_domains, gamma_patches[1], { full_x[1][0],full_x[1][1],full_u[1],prob.time_horizon });
        for(int i=2;i<=prob.num_shooting_nodes;++i)
        {
            objective -= make_function_patch(restricted_domains, gamma_patches[i], { full_x[i][0],full_x[i][1],full_u[i],prob.time_horizon });
        }
        return objective;
    }
    FloatDPBoundsVector optimize()
    {
        ValidatedVectorMultivariateFunctionPatch g;
        ValidatedScalarMultivariateFunctionPatch f;

        for(int i=0;i<=prob.num_shooting_nodes;++i)
        {   
            state_result = calculateNextState(prob.dynamics_eqs, i);

            phi_patches.append(std::get<0>(state_result));
            gamma_patches.append(std::get<1>(state_result));

            X.append(RealVectorVariable("x" + to_string(i), prob.num_states));
            U.append(RealVariable("u" + to_string(i)));
        }
        PRINT(X);
        X_DOMS = getDomains();

        // FIXME: 1 loop
        for(int i=0;i<=prob.num_shooting_nodes;++i) {XU_DOMAIN.append(X_DOMS[i]);}
        for(int i=0;i<=prob.num_shooting_nodes;++i) {XU_DOMAIN.append(prob.u0dom);}
        
        g = createConstraints(phi_patches, X, U, XU_DOMAIN);
        f = createObjective(gamma_patches, X, U, XU_DOMAIN);

        PRINT(g.argument_size());
        PRINT(f.argument_size());
        NonlinearInfeasibleInteriorPointOptimiser nlio;

        int num_decision_vars = prob.num_shooting_nodes * (prob.num_states + prob.num_inputs);
        int num_constraints = (prob.num_shooting_nodes-1) * (prob.num_states);

        List<IntervalDomainType> D_temp;
        List<IntervalDomainType> C_temp;
        // putting  {-0.5_x, +0.5_x} for all state vars and {-0.125_x, +0.125_x} for input vars
        for(int i=0;i<prob.num_shooting_nodes*prob.num_states;i++) {D_temp.append({-0.25_x, +0.25_x});}
        for(int i=prob.num_shooting_nodes*prob.num_states;i<num_decision_vars;i++) {D_temp.append({-0.125_x, +0.125_x});}
        for(int i=0;i<num_constraints;i++) {C_temp.append({-0.0_x, +0.0_x});}

        auto D = Box<ExactIntervalType>(Vector<IntervalDomainType>(D_temp));
        auto C = Box<ExactIntervalType>(Vector<IntervalDomainType>(C_temp));
 
        FloatDPBoundsVector u_optimal = nlio.minimise(f, D, g, C);

        PRINT(u_optimal);
        return u_optimal;
    }
};
int main()
{
    RealVectorVariable x("x", 2);
    RealVariable u("u");
    RealVariable t("t");
    RealVariable c("c");

    BoxDomainType xdom = { { -1.0_x, +1.0_x }, { -1.0_x, +1.0_x } };
    IntervalDomainType udom = { -1.125_x, +1.125_x };
    IntervalDomainType cdom = { -1.0_x, 1.0_x };
    IntervalDomainType tdom = { 0, 0.0625_x};
    auto i_doms = product(xdom, udom, cdom, tdom);

    Vector<RealExpression> x_dot = { x[1],u };
    RealExpression u_dot = 0;
    RealExpression c_dot = pow(x[0], 2) + pow(x[1], 2) + pow(u, 2);
    
    // === RHS ===
    EffectiveVectorMultivariateFunction f_xuc = Function({ x[0], x[1], u, c }, { x_dot[0], x_dot[1], u_dot, c_dot });

    const int NUM_STATE = x.size();
    const int NUM_INPUT = 1;
    const int NUM_SHOOTING_NODES = 2;
    const StepSizeType TIME_HOR = 0.0625_x;

    PRINT(NUM_SHOOTING_NODES);
    problemFormulation problem_setup = problemFormulation(NUM_STATE, NUM_INPUT, NUM_SHOOTING_NODES, TIME_HOR);
    problem_setup.setDynamics(f_xuc, i_doms);

    problemOptimization problem_opt = problemOptimization(problem_setup);
    problem_opt.optimize();

    return 0;
};