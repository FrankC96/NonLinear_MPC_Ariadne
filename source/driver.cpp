/* 
* Nonlinear Model Predictive Controller
* author: Chalaris Fragkiskos
* date: April 2024 
* Follows the pattern :
* https://github.com/AleksandarHaber/Model-Predictive-Control-for-Linear-Systems-in-Cpp-by-Using-Eigen-Library/
*/
#include <iostream>
#include <ariadne.hpp>
#include <ariadne/algebra/matrix.hpp>
#include <ariadne/algebra/vector.hpp>
#include <ariadne/solvers/runge_kutta_integrator.hpp>

using namespace Ariadne;
using state = DifferentialInclusion::EnclosureType::BoundingBoxType::EuclideanSetType;

class ModelPredictiveController {
    public:
        ModelPredictiveController(
            DottedRealAssignments dynamics,
            unsigned int predictionHorizon,
            unsigned int controlHorizon,
            Matrix<FloatDP> Q,
            Vector<FloatDP> R,
            Vector<FloatDP> initialStateVector,
            Vector<FloatDP> initialInputVector,
            Vector<FloatDP> ref) : dynamics(dynamics),
                predictionHorizon(predictionHorizon), 
                controlHorizon(controlHorizon), 
                Q(Q),
                R(R),
                x0(x0),
                u0(u0),
                ref(ref){}

        Vector<FloatDP> solveMPC() {
            Vector<FloatDP> u({0.0_x}, double_precision);
            
            return u;
            }

        // FloatDPApproximationVector propagateDynamics_RK4(Vector<FloatDP> x, Vector<FloatDP> u) {
        //     THIS WON'T RUN BECAUSE RK4 WON'T ACCEPT A DottedRealAssignments object

        //     RealVariable delta("delta");
        //     RealVariablesBox inputs={u.get(0)<=delta<=u.get(0)}; // ATT: only for 1 input systems

        //     DifferentialInclusion ivf(this->dynamics, inputs);

        //     RungeKutta4Integrator integrator(1e-8);

        //     Point<FloatDP> state_pt = x;  // initial state to continue from 


        //     FloatDPApproximationVector nextState = integrator.step(this->dynamics, state_pt, FloatDP(1e-8_x, double_precision));

        //     return nextState;
        // }
        state propagateDynamics_Picard(Vector<FloatDP> x, Vector<FloatDP> u) {
            
            RealVariable delta("delta");
            RealVariablesBox inputs={u.get(0)<=delta<=u.get(0)}; // ATT: only for 1 input systems

            ThresholdSweeperDP sweeper(DoublePrecision(), 1e-8);
            TaylorPicardIntegrator integrator(
                step_maximum_error = 1e-3,
                sweeper,
                lipschitz_tolerance = 0.5_x,
                minimum_temporal_order = 4,
                maximum_temporal_order = 12);
                
                RealVariable uu("uu"), vv("vv"), rr("rr"), ppsi("ppsi"), xx("xx"), yy("yy");
                RealVariablesBox initial = {
                    {5<=uu<=10},
                    {-5<=vv<=10},
                    {-5<=rr<=2},
                    {-1<=ppsi<=1},
                    {500<=xx<=1000},
                    {500<=yy<=1000}
                    };


            DifferentialInclusion ivf(this->dynamics, inputs);
            SizeType period_of_parameter_reduction = 12;
            ExactDouble ratio_of_parameters_to_keep = 6.0_x;
            LohnerReconditioner reconditioner(initial.variables().size(),inputs.variables().size(),period_of_parameter_reduction,ratio_of_parameters_to_keep);

            double step = 1.0/32;
            auto evolver = DifferentialInclusionEvolver(ivf, integrator, reconditioner);
            evolver.configuration().set_maximum_step_size(step);

            Real evolution_time = 5;
            auto orbit = evolver.orbit(initial,evolution_time);

            auto n = ivf.dimension();
            state graphics_box(n);
            for(auto set:orbit.reach()) {
                graphics_box = hull(graphics_box,set.euclidean_set().bounding_box());
            }

            return graphics_box;
        };

        auto displayValues(auto x) {return x;}
    
    private:

        DottedRealAssignments dynamics = {};
        unsigned int predictionHorizon = 10;
        unsigned int controlHorizon = 5;

        Matrix<FloatDP> Q = Matrix<FloatDP>::identity(6, double_precision);
        Vector<FloatDP> R = Vector<FloatDP>::zero(6, double_precision);
        Vector<FloatDP> x0 = Vector<FloatDP>::zero(6, double_precision);
        Vector<FloatDP> u0 = Vector<FloatDP>::zero(1, double_precision); 
        Vector<FloatDP> ref = Vector<FloatDP>::zero(6, double_precision);
        
        // [IGNORE THIS] DEFINE HISTORY MATRICES ACCUMULATING RESULTS FOR 100 TIMESTEPS
        Matrix<FloatDP> states = Matrix<FloatDP>::zero(6, 100, double_precision);
        Matrix<FloatDP> inputs = Matrix<FloatDP>::zero(1, 100, double_precision);
        Matrix<FloatDP> outputs = Matrix<FloatDP>::zero(6, 100, double_precision); // assuming perfect state knowledge
};

int main() {

    // ** MPC params
    unsigned int Np = 20;  // Prediction horizon
    unsigned int Nc = 18;  // Controll horizon

    RealConstant sampling("sampling", 0.05_dec);
    // end of MPC params

    // ** Model params
    // pass this to a discrete integrator (RK4)
    RealConstant L("L", 350), m("m", 356000), n("n", 87.6_dec), T("T", 3000), g("g", 9.81_dec);

    RealConstant X1("X1", 1.070_dec), X2("X2", -0.0344_dec), X3("X3", 0), X4("X4", 0.0267_dec),
                     X5("X5", 1.492_dec), X6("X6", -300), X7("X7", -0.42_dec), X8("X8", 0.80_dec);
    RealConstant Y1("Y1", 1.820_dec), Y2("Y2", 0), Y3("Y3", -0.7_dec), Y4("Y4", -0.820_dec),
                     Y5("Y5", -0.910_dec), Y6("Y6", 0.234_dec), Y7("Y7", -2.150_dec), Y8("Y8", 0),
                     Y9("Y9", 0.5_dec), Y10("Y10", -0.5_dec), Y11("Y11", 0.176_dec), Y12("Y12", 1.2_dec),
                     Y13("Y13", 0.040_dec);
    RealConstant R1("Z1", 0), R2("Z2", 0.103_dec), R3("Z3", -0.350_dec), R4("Z4", -0.163_dec),
                     R5("Z5", -0.45_dec), R6("Z6", -0.177_dec), R7("Z7", 0.064_dec), R8("Z8", -0.129_dec),
                     R9("Z9", -0.180_dec), R10("Z10", 0.180_dec), R11("Z11", -0.083_dec), R12("Z12", -0.564_dec),
                     R13("Z13", -0.020_dec);
    // end of Model params

    // ** System vars
    RealVariable u("u"), v("v"), r("r"), psi("psi"), x("x"), y("y"), delta("delta");

    DottedRealAssignments equationsOfMotion = {
        dot(u)=
            div((div(1, L) * X2 * pow(u, 2)), X1) + div(div(1, g*pow(L, 2))*X3*pow(u, 4), X1) + 
			div(L*X4*pow(r, 2), X1) + div(X5*v*r, X1) + div(div(1, g*pow(L, 2))*X6*u*abs(v)*pow(v, 2), X1)
            + div(X7*(div(1, L)*Y11*pow(u, 2)*delta + Y12*div(T, m))*delta, X1) + div(div(T, m)*X8, X1),
        dot(v) = 
            div(div(1, L)*Y3*u*v, Y1) + div(Y4*u*r, Y1) + div(div(1, sqrt(g*L)*L)*Y5*pow(u, 2)*v, Y1)
            + div(div(1, sqrt(g*L))*Y6*pow(u, 2)*r, Y1) + div(div(1, L)*Y7*v*abs(v) , Y1) + div(L*Y8*r*abs(r), Y1)
            + div(Y9*abs(v)*r, Y1) + div(Y10*v*abs(r), Y1) + div(div(1, L)*Y11*pow(u, 2)*delta, Y1) +
            div(Y12*div(T, m)*delta, Y1) + div(Y13*div(T, m), Y1),
        dot(r) = div(div(1, L)*R3*u*v, L*R2) + div(R4*u*r, L*R2) + div(div(1, sqrt(g*L)*L)*R5*pow(u, 2)*v , L*R2)
            + div(div(1, sqrt(g*L))*R6*pow(u, 2)*r, L*R2) + div(div(1, L)*R7*v*abs(v), L*R2) + div(L*R8*r*abs(r), L*R2)
            + div(R9*abs(v)*r, L*R2) + div(R10*v*abs(r), L*R2) + div(div(1, L)*R11*pow(u, 2)*delta, L*R2) +
            div(R12*div(T, m)*delta, L*R2) + div(R13*div(T, m), L*R2),
        dot(psi) = r,
        dot(x) = u*cos(psi) - v*sin(psi),
        dot(y) = u*sin(psi) + v*cos(psi)
    };
    // end of system vars

    Matrix<FloatDP> Q = Matrix<FloatDP>::identity(6, double_precision);
    Vector<FloatDP> R = Vector<FloatDP>::one(6, double_precision);

    // **Define the desired trajectory
    unsigned int timeSteps = 300;
    // just a stabilization trajectory, drive everything to 0
    Vector<FloatDP> desiredTraj = Vector<FloatDP>::zero(timeSteps, double_precision);
    // end of desird trajectory

    // ** Main MPC loop

    // these will will form the MPC object to make the calculations for timestep k = 0,
    // they need to be initialized with the MPC object.
    Vector<FloatDP> initialInput = Vector<FloatDP>::zero(1, double_precision);  // initial input
    Vector<FloatDP> initialState = Vector<FloatDP>::zero(6, double_precision);  // initial state
    
    Vector<FloatDP> tempInput = Vector<FloatDP>::zero(1, double_precision);  // temp input
    Vector<FloatDP> tempState = Vector<FloatDP>::zero(6, double_precision);  // temp state
    tempState = initialState;

    ModelPredictiveController mpc(equationsOfMotion, Np, Nc, Q, R, initialState, initialInput, desiredTraj);
    for(int i=0;i<1;i++) {
        Vector<FloatDP> tempInput = mpc.solveMPC();  // this should return the next optimal input to be fed again to the system.
        // FloatDPApproximationVector newState = mpc.propagateDynamics_RK4(tempState, tempInput);  // this should accept the previous state, and the next computed input.
        state  newState = mpc.propagateDynamics_Picard(tempState, tempInput);  // this should accept the previous state, and the next computed input.
        // std::cout << B << std::endl;
    
    }
    // end of main MPC loop

    return 0;
}
