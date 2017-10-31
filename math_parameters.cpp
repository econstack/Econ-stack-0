#include <math.h>
#include <fstream>
#include <string>
#include <conio.h>
#include <iostream>
#include <iomanip>

using namespace std;

#include "math_parameters.h"

ofstream log_fout("log_data.out");	//	need to define log_fout in main.cpp
//  math_grunts.h -------------------------------------------------------------------------------
//const int TWO_SIDE_FLAG = 0;          //  = 1 for 2 sided numerical differentiation
const int UNUSED = -999;
const int OUT_OF_RANGE_FLAG = -1;
int ON = 1;
int OFF = 0;

//	numerical differentiation step size
double stepJacob = 0.001;
//  Nlin_solvers.h parameters --------------------------------------------------------------------
int Nmaxiter_xnext = 10;          //  max attempts to push candidate x into constraint set
int MaxIter_Nlin = 1000;          //  max iterations for non-linear solver
double eps_Nlin = 0.0001;         //  define zero for non-linear solver
double eps_MatInv = 0.0000000001;   //  define zero for matrice inversion
double stepNewton = 0.1;         //  Newton step, s, in Nlin solver: x' = x - s*Jinv(x)*x
double stepNewton_dp = 0.05; // option: max percent change in x
double step_NlinNewton = 0.1;     //  Newton step, s, in Nlin solver: x' = x - s*Jinv(x)*x
double step_Broyden = 1.00;       //  Broyden step, s, x' = x + s*step
double e_math = 2.718281828459;   //  mathematical constant
double pi_math = 3.14159265358;   //  mathematical constant

//	ERROR FLAGS --------------------------------------------------------------------------------
int FAILURE_FLAG = -1;
int NONCONVERGENCE_FLAG = -2;
int X_STUCK_FLAG = -3;
int EXCEEDED_MAXITER_FLAG = -4;
int SINGULAR_MATRIX_FLAG = -5;
int NO_PIVOT_FLAG = -5;
int CONSTRAINT_SET_VIOLATION_FLAG = -6;
int X_SUBOPTIMAL_FLAG = -7;
int NOT_LOCALLY_CONVEX_FLAG = -8;
int NEWTON_POWELL_FLAG = -9;
int OVERFLOW_FLAG = -10;

//	FOUT FLAGS --------------------------------------------------------------------------------
//	rescaling next guess x into constraint set
int fout_rescale_xnext_FLAG = OFF;
//	newton and broyden method
int fout_nlin_final_FLAG = ON;
int fout_nlin_iter_FLAG = OFF;
int fout_nlin_final_NMIN = 10; // does not write out final summary if (#iter < X)
//	fminBFGS and polytope in R(N) and fminR1 in R(1)
int fout_fmin_iter_FLAG = OFF;
int fout_fmin_diag_FLAG = OFF;
int fout_fmin_final_FLAG = ON;
int fout_update_hessian_fmin_FLAG = OFF;
//	linemin newton and golden
int fout_linemin_iter_FLAG = OFF;
int fout_linemin_final_FLAG = OFF;
int FOUT_GET_GRADIENT_DIAG_FLAG = ON;
int FOUT_NLIN_DIAG_FLAG = OFF;
int FOUT_BROYDEN_DIAG_FLAG = OFF;
//	calib_mom2mom
int fout_calib_mom2mom_diag_FLAG = ON;
//	misc functions in math_grunts

//	misc identifiers and flags  ------------------------------------------------------------
double HELL = -1000000;
int DIFF_3PT_FLAG = 1;
int DIFF_5PT_FLAG = 2;
int fminBFGS_FLAG = ON;
int update_x0_FLAG = ON;

int fout_xnext_broyden_diag_FLAG = OFF;
int fout_jxnext_broyden_FLAG = OFF;
int FOUT_GET_ERGODIC_FINAL_FLAG = ON;
int FOUT_GET_ERGODIC_ITER_FLAG = ON;
int Nquad = 20;
const int N_fout_nlin = 50;             //  no. iterations between nlin solver loops for fout

const int Ncalibration_moments = 10;

//	fminBFGS ----------------------------------------------------------------------------
double epsFmin = 0.001;
double deltaFmin = 10*epsFmin;
int NmaxiterFmin = 200;
double FMIN_GRADIENT_STEP_PERCENT = 0.01;
//  fmin.h - polytope -------------------------------------------------------------------
double gamma_polytope = 1.0;          //  init vertices parameter
double alpha_polytope = 1.0;          //  x_reflected scale parameter
double beta_polytope = 2.0;           //  x_expanded scale parameter
double zeta_polytope = 0.5;           //  x_contracted scale parameter
double eta_polytope = 0.5;            //  shrink simplex scale parameter
//  linemin ----------------------------------------------------------------------------
int NmaxiterLinemin = 100;
double epsLinemin = 0.001;
double stepNewtonLinemin = 0.5;
int LINEMIN_GOLDEN_FLAG = 1;
int LINEMIN_NEWTON_FLAG = 2;
//	linemin_newton
double linemin_stepNewton = 0.2;
double linemin_stepJacob = stepJacob;
//  linemin_golden
double alpha_linemin = 2.0;         //  scaling parameter in seeking init triple
double delta_linemin = 1.0;         //  step size in seeking init triple
int NmaxiterLM_triple = 10;			//  max no. iter in seeking init triple


//const int Nfmin = 2;

//  fmin.h - fmin_BFGS ---------------------------------------------------


double eps_polytope = 0.01;           //  stopping condition tolerance


//  diagnostics flags for math.h    (=1 for diagnostics turned on)
int fout_get_d2fdx2_along_s_FLAG = ON;
int fout_get_dfdx_along_s_FLAG = ON;
int fout_get_gradient_FLAG = ON;
int fout_swap_row_FLAG = OFF;

//  diagnostics flags for fmin.h (=1 for diagnostics turned on)
const int fmin_BFGS_DIAG_FLAG = 1;
const int polytope_DIAG_FLAG = 0;

//  diagnostics flags for nlin_solvers.h (=1 for diagnostics turned on)
int newton_method_DIAG_FLAG = 1;
int broyden_method_DIAG_FLAG = 0;
int fout_nlin_labelled_eqn_resids_FLAG = 0;		//	calls fcn to write out eqn resids if ON

// diagnostics flags for newton_calib solvers
int FOUT_CALIB_STEPWISE_FINAL_FLAG = 1;

void toggle_fout_nlin_FLAGS()
{
	FOUT_GET_GRADIENT_DIAG_FLAG = 1 - FOUT_GET_GRADIENT_DIAG_FLAG;
	fout_get_d2fdx2_along_s_FLAG = 1 - fout_get_d2fdx2_along_s_FLAG;
	fout_get_dfdx_along_s_FLAG = 1 - fout_get_dfdx_along_s_FLAG;
	fout_get_gradient_FLAG = 1 - fout_get_gradient_FLAG;
	return;
}