// g++-10 *.cpp -I/usr/local/include/ariadne -L/usr/local/lib -std=c++20 -lariadne -w -fcompare-debug-second
#include <iostream>
#include <ariadne.hpp>
#include <ariadne/solvers/nonlinear_programming.hpp>

using namespace Ariadne;
using ValidatedScalarMultivariateFunction = ScalarMultivariateFunction<ValidatedTag>;
using ValidatedVectorMultivariateFunction = VectorMultivariateFunction<ValidatedTag>;

int main() {
    
    ValidatedScalarMultivariateFunction x = ValidatedScalarMultivariateFunction::coordinate(EuclideanDomain(2),0);
    ValidatedScalarMultivariateFunction y = ValidatedScalarMultivariateFunction::coordinate(EuclideanDomain(2),1);

    ValidatedVectorMultivariateFunction g = ValidatedVectorMultivariateFunction{x, y};
    ValidatedScalarMultivariateFunction f = ValidatedScalarMultivariateFunction{1+x*y};

    ExactBoxType D = ExactBoxType({{-4,5},{-4,5}});
    ExactBoxType C = ExactBoxType({{-4,5},{-4,5}});

    PenaltyFunctionOptimiser solver;
    auto sol = solver.minimise(f, D, g, C);
    return 0;
}