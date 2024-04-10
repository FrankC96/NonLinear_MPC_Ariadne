/* 
* Nonlinear Model Predictive Controller
* author: Chalaris Fragkiskos
* date: April 2024 
* Follows the pattern :
* https://github.com/AleksandarHaber/Model-Predictive-Control-for-Linear-Systems-in-Cpp-by-Using-Eigen-Library/
*/

#include <ariadne_main.hpp>

#include <ModelPredictiveController.h>

using namespace Ariadne;

void ariadne_main() {

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
    Vector<FloatDP> V({1, 2}, double_precision); // define a symbolic state vector

    DottedRealAssignments equationsOfMotion = {
        dot(u)=
            div((div(1, L) * X2 * pow(u, 2)), X1) + div(div(1, g*pow(L, 2))*X3*pow(u, 4), X1) + 
			div(L*X4*pow(r, 2), X1) + div(X5*v*r, X1) + div(div(1, g*pow(L, 2))*X6*u*pos(v)*pow(v, 2), X1)
            + div(X7*(div(1, L)*Y11*pow(u, 2)*delta + Y12*div(T, m))*delta, X1) + div(div(T, m)*X8, X1),
        dot(v) = 
            div(div(1, L)*Y3*u*v, Y1) + div(Y4*u*r, Y1) + div(div(1, sqrt(g*L)*L)*Y5*pow(u, 2)*v, Y1)
            + div(div(1, sqrt(g*L))*Y6*pow(u, 2)*r, Y1) + div(div(1, L)*Y7*v*pos(v) , Y1) + div(L*Y8*r*pos(r), Y1)
            + div(Y9*pos(v)*r, Y1) + div(Y10*v*pos(r), Y1) + div(div(1, L)*Y11*pow(u, 2)*delta, Y1) +
            div(Y12*div(T, m)*delta, Y1) + div(Y13*div(T, m), Y1),
        dot(r) = div(div(1, L)*R3*u*v, L*R2) + div(R4*u*r, L*R2) + div(div(1, sqrt(g*L)*L)*R5*pow(u, 2)*v , L*R2)
            + div(div(1, sqrt(g*L))*R6*pow(u, 2)*r, L*R2) + div(div(1, L)*R7*v*pos(v), L*R2) + div(L*R8*r*pos(r), L*R2)
            + div(R9*pos(v)*r, L*R2) + div(R10*v*pos(r), L*R2) + div(div(1, L)*R11*pow(u, 2)*delta, L*R2) +
            div(R12*div(T, m)*delta, L*R2) + div(R13*div(T, m), L*R2),
        dot(psi) = r,
        dot(x) = u*cos(psi) - v*sin(psi),
        dot(y) = u*sin(psi) + v*cos(psi)
    };
    // end of system vars

    // **Define the desired trajectory
    unsigned int timeSteps = 300;
    Vector<FloatDP> desiredTraj = Vector<FloatDP>::zero(timeSteps, double_precision);
    // end of desird trajectory

    // ** Main MPC loop
    Float<DP> temporaryInput = FloatDP(0.0_x, double_precision);  // initial input
    Vector<FloatDP> temporaryState = Vector<FloatDP>::zero(6, double_precision);  // initial state
    ModelPredictiveController mpc(equationsOfMotion, Np, Nc, temporaryState, temporaryInput, desiredTraj);
    
    for(int i=0;i<20;i++) {
        mpc.solveMPC();  // this should return the next optimal input to be fed again to the system.
        mpc.propagateDynamics(temporaryState, temporaryInput);  // this should accept previous state, and the next computed input.

    }

    // mpc.saveData();  // Save state trajectories to a csv file
    // end of main MPC loop
    // return 0;
}