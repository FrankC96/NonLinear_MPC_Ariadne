
#include "config.hpp"
#include "algebra/algebra.hpp"
#include "function/function_patch.hpp"
#include "function/domain.hpp"
#include "symbolic/expression.hpp"
#include "symbolic/expression_patch.hpp"
#include "symbolic/expression_set.hpp"
#include "symbolic/space.hpp"

using namespace std;
using namespace Ariadne;

int main()
{
    RealVectorVariable x("x", 2);
    RealVariable t("t");
    RealVariable u("u"), u0("u0"), u1("u1");

    Dyadic h(0.125_x);

    BoxDomainType xdom = { {0.0_x,0.5_x},{0.25_x,0.75_x} };
    BoxDomainType xrdom = { {0.125_x,0.375_x},{0.375_x,0.625_x} };
    IntervalDomainType tdom = { 0,h };
    IntervalDomainType udom = { -1.0_x,+1.0_x };

    BoxDomainType xtudom = product(xdom, tdom, udom);

    VectorVariableBoxDomainType xd(x, xdom);
    xd = x | xdom;
    VariableIntervalDomainType xd0(x[0], xdom[0]);
    VariableIntervalDomainType xd1(x[1], xdom[1]);
    VariableIntervalDomainType td(t, tdom);
    td = t | tdom;
    VariableIntervalDomainType ud(u, udom);
    VariableIntervalDomainType ud0(u0, udom);
    VariableIntervalDomainType ud1(u1, udom);
    // ARIADNE_TEST_PRINT(xd);
    // ARIADNE_TEST_PRINT(td);

    auto xtud = ValidatedVectorExpressionPatch({ x | xdom,t | tdom,u | udom }, { x[0],x[1],t,u });
    // ARIADNE_TEST_PRINT(xtud);
    auto xhud = ValidatedVectorExpressionPatch({ xd[0],xd[1],ud }, { x[0],x[1],h,u });
    // ARIADNE_TEST_PRINT(xhud);

    EffectiveVectorMultivariateFunction phif({ x[0],x[1],t,u }, { x[0] + u * t, (x[1] - x[0] + u) * exp(-t) + u * t + x[0] - u });
    ValidatedVectorMultivariateFunctionPatch phifp = ValidatedVectorMultivariateRestrictedFunction(phif, xtudom);
    // ARIADNE_TEST_PRINT(phifp);
    ValidatedVectorMultivariateFunctionPatch phi1fp = restriction(phifp, product(xrdom, tdom, udom));

    ValidatedVectorExpressionPatch phiep = phifp({ x[0],x[1],t,u0 });
    // ARIADNE_TEST_PRINT(phiep);
    phiep = phifp({ x,t,u0 });
    // ARIADNE_TEST_PRINT(phiep);
    ValidatedVectorExpressionPatch phi1ep = phi1fp({ x[0],x[1],t,u });
    // ARIADNE_TEST_PRINT(phi1ep);
    ValidatedVectorExpressionPatch phieptu({ phi1ep,t | tdom,u | udom });
    // ARIADNE_TEST_PRINT(phieptu);
    ValidatedVectorExpressionPatch phi2ep = phifp(phieptu);
    // ARIADNE_TEST_PRINT(phi2ep);

    phi2ep = phifp({ phi1fp({x,t,u0}),t | tdom,u | udom });
    // ARIADNE_TEST_PRINT(phi2ep);
    auto phi2fp = ValidatedVectorMultivariateFunctionPatch({ xd,td,ud0,ud1 }, phifp({ phi1fp({x,t,u0}),td,ud1 }));
    phi2fp = ValidatedVectorMultivariateFunctionPatch({ x | xdom,t | tdom,u0 | udom,u1 | udom }, phifp({ phi1fp({x,t,u0}),td,ud1 }));
    // ARIADNE_TEST_PRINT(phi2fp);

    ValidatedVectorExpressionPatch psiep = phifp({ x,h,u });

    // ARIADNE_TEST_PRINT(psiep);
    psiep = phifp({ x[0],x[1],h,u });
    // ARIADNE_TEST_PRINT(psiep);

    ValidatedVectorExpressionPatch psi1ep = phi1fp({ x,h,u0 });
    ValidatedVectorExpressionPatch psi2ep = phifp({ phi1fp({x,h,u0}),h,ud1 });
    // ARIADNE_TEST_PRINT(psi1ep);
    // ARIADNE_TEST_PRINT(psi2ep);

    auto psifp = ValidatedVectorMultivariateFunctionPatch({ xd,ud }, psiep);
    psifp = ValidatedVectorMultivariateFunctionPatch({ x | xdom,u | udom }, psiep);
    // ARIADNE_TEST_PRINT(psifp);

    auto psi1fp = restriction(psifp, product(xrdom, udom));
    auto psi2fp = ValidatedVectorMultivariateFunctionPatch({ x | xdom,u0 | udom,u1 | udom }, psifp({ psi1fp({x,u0}),ud1 }));
    //psi2fp = ValidatedVectorMultivariateFunctionPatch({x|xdom,u0|udom,u1|udom},psifp({psifp({xd,ud0}),ud1}));
    // ARIADNE_TEST_PRINT(psi2fp);

    cout << "psi1fp " << psi1fp.domain() << endl;
    cout << "psi2fp " << psi2fp.domain() << endl;
    ValidatedVectorMultivariateFunction psi12ep = join(psi1fp, psi2fp);
    // ARIADNE_TEST_PRINT(psi12ep);

    return 0;
};