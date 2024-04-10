#include "ModelPredictiveController.h"

using namespace Ariadne;

ModelPredictiveController::ModelPredictiveController() {
    // DEFAULT CONSTRUCTOR
}

ModelPredictiveController::ModelPredictiveController(
    DottedRealAssignments eqs,
    unsigned int np,
    unsigned int nc,
    Vector<FloatDP> x0,
    Float<DP> u0,
    Vector<FloatDP> ref) {
        dynamics = eqs;
        n = 6; m = 1; r = 2;

        k = 0;

        unsigned int maxSimulationSamples = ref.size() - np;

        states.resize(n, maxSimulationSamples);
        // ATT: no need to set to 0, it already is? this is the constructor, so it better be.

        inputs.resize(m, maxSimulationSamples-1);  // ATT: discarding 1 element, to add the first controlledInput?
        // same

        outputs.resize(r, maxSimulationSamples-1);  // ATT: discarding 1 element, to add the last computedOutput?
    }

Float<DP> ModelPredictiveController::solveMPC() {
    // solve the NLP problem and return the optimal input
    Float<DP> a = FloatDP(0.0_x, double_precision);
    return a;

}

Vector<FloatDP> ModelPredictiveController::propagateDynamics(Vector<FloatDP> x, Float<DP> u) {
    Vector<FloatDP> b = Vector<FloatDP>::zero(6, double_precision);
    return b;
}