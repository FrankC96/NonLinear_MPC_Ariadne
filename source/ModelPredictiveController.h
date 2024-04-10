#ifndef MODELPREDICTIVECONTROLLER_H
#define MODELPREDICTIVECONTROLLER_H

#include <ariadne.hpp>
#include <ariadne/numeric/floatdp.hpp>
#include <ariadne/algebra/matrix.hpp>
#include <ariadne/algebra/vector.hpp>
#include <ariadne/solvers/runge_kutta_integrator.hpp>

using namespace Ariadne;

using FloatDPVector = Vector<FloatDP>;

class ModelPredictiveController {
    public:
        ModelPredictiveController();

        // the following constructor will initialize all variables required for the MPC object.

        ModelPredictiveController(
            DottedRealAssignments dynamics,
            unsigned int predictionHorizon,
            unsigned int controlHorizon,
            Vector<FloatDP> initialStateVector,
            Float<DP> initialInputVector,
            Vector<FloatDP> ref);

        static Float<DP> solveMPC();
        Vector<FloatDP> propagateDynamics(Vector<FloatDP> x, Float<DP> u);

        // WHY? this function is specified as ModelPredictiveController function, though computeControlInputs is not, also const at the end?
        // void saveData(string desiredControlTrajectoryTotalFile, string inputsFile, string statesFile, string outputsFile, string OFile, string MFile) const;
    
    private:
        unsigned int k;  // index constant for the difference equations
        unsigned int m;  // m, n ,r = input_dimension, state_dimension, output_dimension'
        unsigned int n;
        unsigned int r;

        // TODO: dynamics object to be passed to RK4
        DottedRealAssignments dynamics;
        Matrix<FloatDP> Q = Matrix<FloatDP>::identity(6, double_precision);
        FloatDPVector R;
        FloatDPVector x0;   // initial state vector
        FloatDPVector ref;  // reference vector
        
        unsigned int np, nc;  // prediction_horizon, control horizon

        // ALREADY DEFINED MATRICES ACCUMULATING RESULTS FOR 100 TIMESTEPS
        Matrix<FloatDP> states = Matrix<FloatDP>::zero(n, 100, double_precision);
        Matrix<FloatDP> inputs = Matrix<FloatDP>::zero(m, 100, double_precision);
        Matrix<FloatDP> outputs = Matrix<FloatDP>::zero(r, 100, double_precision);

        static FloatDP newControl;  // optimal computed control
};
#endif