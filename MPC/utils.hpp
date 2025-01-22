using namespace Ariadne;

void pprint(std::string str, auto var)
{
    std::cout << str << " => \n" << var << "\n";
    for (int i = 0;i < 50;++i) { std::cout << "-"; };
    std::cout << "\n";
};
#define PRINT(expr) { pprint(#expr,(expr)); }


ValidatedVectorMultivariateFunctionPatch integrate_dynamics(EffectiveVectorMultivariateFunction f, ExactBoxType x_dom, StepSizeType h)
{
    auto integrator = TaylorPicardIntegrator(1E-2);

    ValidatedVectorMultivariateFunctionPatch x_d = integrator.flow_step(f, x_dom, h);
    return x_d;
};
ValidatedVectorMultivariateFunctionPatch project_function(ValidatedVectorMultivariateFunctionPatch f, Range p)
{
    return Array<ValidatedScalarMultivariateFunctionPatch>(p.size(), [&f, &p](SizeType i) { return f[p[i]]; });
};

typedef Variant<RealVectorVariable, RealVariable> RealOrVectorVariable;

class VariableIntervalDomainType : public VariableInterval<FloatDP> {
    using VariableInterval<FloatDP>::VariableInterval;
};
VariableIntervalDomainType operator| (RealVariable v, IntervalDomainType d) {
    return VariableIntervalDomainType(v, d);
}

class VectorVariableBoxDomainType {
    RealVectorVariable _var;
    BoxDomainType _dom;
public:
    RealVectorVariable variable() const { return _var; }
    BoxDomainType set() const { return _dom; }
    VectorVariableBoxDomainType(RealVectorVariable v, BoxDomainType d) : _var(v), _dom(d) {};
    VariableIntervalDomainType operator[](SizeType i) const { return _var[i] | _dom[i]; }
    friend OutputStream& operator<<(OutputStream& os, VectorVariableBoxDomainType const& vd) {
        return os << vd._var << '|' << vd._dom;
    }
};
VectorVariableBoxDomainType operator| (RealVectorVariable vv, BoxDomainType d) {
    return VectorVariableBoxDomainType(vv, d);
}

typedef Variant<VariableIntervalDomainType, VectorVariableBoxDomainType> VariableIntervalOrBoxDomainType;

ValidatedScalarMultivariateFunctionPatch make_function_patch(List<VariableIntervalOrBoxDomainType> sv_arg_doms, RealExpression e) {
    List<VariableIntervalDomainType> arg_doms;
    List<RealVariable> args; for (auto arg_dom : arg_doms) { args.append(arg_dom.variable()); }
    RealSpace spc(args);
    List<IntervalDomainType> doms; for (auto arg_dom : arg_doms) { doms.append(arg_dom.interval()); }
    BoxDomainType dom = BoxDomainType(Array<IntervalDomainType>(doms.begin(), doms.end()));
    EffectiveScalarMultivariateFunction fe = make_function(spc, e);
    ValidatedScalarMultivariateTaylorFunctionModelDP fem(dom, fe, ThresholdSweeperDP(dp, 1e-8));
    return fem;
}

template<class... TS> OutputStream& operator<<(OutputStream& os, Variant<TS...> const& var) {
    std::visit([&os](auto const& t) {os << t;}, var); return os;
}
template<class T> OutputStream& operator<<(OutputStream& os, InitializerList<T> const& lst) {
    return os << List<T>(lst);
}

List<VariableIntervalDomainType> make_scalar_variable_domains(List<VariableIntervalOrBoxDomainType> const& sv_arg_doms) {
    List<VariableIntervalDomainType> arg_doms;
    for (auto sv_arg_dom : sv_arg_doms) {
        if (std::holds_alternative<VariableIntervalDomainType>(sv_arg_dom)) {
            arg_doms.append(std::get<VariableIntervalDomainType>(sv_arg_dom));
        }
        else {
            auto bx_arg_dom = std::get<VectorVariableBoxDomainType>(sv_arg_dom);
            for (SizeType i = 0; i != bx_arg_dom.variable().size(); ++i) {
                arg_doms.append(bx_arg_dom[i]);
            }
        }
    }
    return arg_doms;
}

ValidatedVectorMultivariateFunctionPatch make_function_patch(List<VariableIntervalOrBoxDomainType> sv_arg_doms, List<RealExpression> es) {
    List<VariableIntervalDomainType> arg_doms = make_scalar_variable_domains(sv_arg_doms);
    List<RealVariable> args; for (auto arg_dom : arg_doms) { args.append(arg_dom.variable()); }
    RealSpace spc(args);
    List<IntervalDomainType> doms; for (auto arg_dom : arg_doms) { doms.append(arg_dom.interval()); }
    BoxDomainType dom = BoxDomainType(Array<IntervalDomainType>(doms.begin(), doms.end()));
    Vector<RealExpression> e(es);
    EffectiveVectorMultivariateFunction fe = make_function(spc, e);
    ValidatedVectorMultivariateTaylorFunctionModelDP fem(dom, fe, ThresholdSweeperDP(dp, 1e-8));
    return fem;
}

ValidatedScalarMultivariateFunctionPatch make_function_patch(List<VariableIntervalOrBoxDomainType> sv_arg_doms, ValidatedScalarMultivariateFunctionPatch f, List<RealExpression> es) {
    return compose(cast_unrestricted(f), make_function_patch(sv_arg_doms, es));
}

ValidatedVectorMultivariateFunctionPatch make_function_patch(List<VariableIntervalOrBoxDomainType> sv_arg_doms, ValidatedVectorMultivariateFunctionPatch f, List<RealExpression> es) {
    return compose(cast_unrestricted(f), make_function_patch(sv_arg_doms, es));
}

// INTEGRATE DYNAMICS
// LOSE THE C TERM
// EXTRACT X EQS
// EXTRACT GAMMA EQ 
// RESTRICT [X EQS] [QUESTION] isn't phi0 = psi0? 

// RETURN XDOM, TDOM, PSI, GAMMA

Tuple<
        BoxDomainType,
        IntervalDomainType, 
        ValidatedVectorMultivariateFunctionPatch,
        ValidatedScalarMultivariateFunctionPatch
        > makeStep
(
    int                                 iter,
    EffectiveVectorMultivariateFunction f_xuc,
    RealVectorVariable                  x,
    List<RealVariable>                  vars,  // FIXME: unify with RealOrVectorVariable
    Tuple<BoxDomainType, IntervalDomainType, IntervalDomainType, IntervalDomainType> restr_doms,
    StepSizeType                        h
)
{   
    RealVariable("u" + to_string(iter));

    BoxDomainType       x_dom_next;
    IntervalDomainType  t_dom_next;

    RealVariable u = vars[0];
    RealVariable c = vars[1];
    RealVariable t = vars[2];

    BoxDomainType x_dom = std::get<0>(restr_doms);
    IntervalDomainType u_dom = std::get<1>(restr_doms);
    IntervalDomainType c_dom = std::get<2>(restr_doms);
    IntervalDomainType t_dom = std::get<3>(restr_doms);

    VectorVariableBoxDomainType xr = x | x_dom;
    VariableIntervalDomainType ur = u | u_dom, tr = t | t_dom;

    BoxDomainType xuc_dom = product(x_dom, u_dom, c_dom);

    ValidatedVectorMultivariateFunctionPatch phigamma_ = integrate_dynamics(f_xuc, xuc_dom, h);
    ValidatedVectorMultivariateFunctionPatch phigamma = make_function_patch({ xr,ur,tr }, phigamma_, { x[0],x[1],u,0,t });
    ValidatedVectorMultivariateFunctionPatch phi = project_function(phigamma, Range(0, x.size()));
    ValidatedScalarMultivariateFunctionPatch gamma = phigamma[x.size()+1]; // [PROBLEM +1?] 
    ValidatedVectorMultivariateFunctionPatch psi = make_function_patch({ xr,ur }, phi, { x[0],x[1],u,h });

    x_dom_next = cast_exact_box(psi.range());
    t_dom_next = {0, h};

    
    return make_tuple(x_dom_next, t_dom_next, phi, gamma);
};