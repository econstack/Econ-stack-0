//  math_grunts.h -------------------------------------------------------------------------------
//	global log file 
extern int temp;
extern std::ofstream log_fout;
//	basic swtiches
extern int ON;
extern int OFF;
//	math constants
extern double e_math; 
extern double pi_math; 
//	numerical differentiation step size
extern double stepJacob;         
//  Nlin_solvers.h --------------------------------------------------------------------
extern int Nmaxiter_xnext; //  max attempts to push candidate x into constraint set
extern int MaxIter_Nlin; //  max iterations for non-linear solver
extern double eps_Nlin; //  define zero for non-linear solver
extern double eps_MatInv; //  define zero for matrice inversion
extern double stepNewton; //  Newton step, s, in Nlin solver: x' = x - s*Jinv(x)*x
extern double step_NlinNewton; //  Newton step, s, in Nlin solver: x' = x - s*Jinv(x)*x
extern double step_Broyden; //  Broyden step, s, x' = x + s*step
extern double stepNewton_dp;
//	fminBFGS ----------------------------------------------------------------------------
extern double epsFmin; //	convergence tolerance
extern double deltaFmin; //	changes between x(k) and x(k+1) tolerance
extern int NmaxiterFmin; //	max iterations
extern double FMIN_GRADIENT_STEP_PERCENT; 
//  fmin.h - polytope -------------------------------------------------------------------
extern double eps_polytope; //  stopping condition tolerance
extern double gamma_polytope; //  init vertices parameter
extern double alpha_polytope; //  x_reflected scale parameter
extern double beta_polytope; //  x_expanded scale parameter
extern double zeta_polytope; //  x_contracted scale parameter
extern double eta_polytope; //  shrink simplex scale parameter
//  linemin ----------------------------------------------------------------------------
extern int NmaxiterLinemin; 
extern double epsLinemin; //	convergence tolerance 
extern double stepNewtonLinemin;
extern int LINEMIN_GOLDEN_FLAG; //	flags for newton vs golden section in linemin call from fmin
extern int LINEMIN_NEWTON_FLAG;
//	linemin_newton
extern double linemin_stepNewton; // step size from x(k) to x(k+1)
extern double linemin_stepJacob; //	step size in numerical diff
//  linemin_golden 
extern double alpha_linemin; //  scaling parameter in seeking init triple
extern double delta_linemin; //  step size in seeking init triple
extern int NmaxiterLM_triple; //  max no. iter in seeking init triple

//	ERROR FLAGS --------------------------------------------------------------------------------
extern int NO_PIVOT_FLAG;
extern int FAILURE_FLAG;
extern int NONCONVERGENCE_FLAG;
extern int X_STUCK_FLAG;
extern int EXCEEDED_MAXITER_FLAG;
extern int SINGULAR_MATRIX_FLAG;
extern int CONSTRAINT_SET_VIOLATION_FLAG;
extern int X_SUBOPTIMAL_FLAG;
extern int NEWTON_POWELL_FLAG;
extern int OVERFLOW_FLAG;

//	FOUT FLAGS --------------------------------------------------------------------------------
//	rescaling next guess x into constraint set
extern int fout_rescale_xnext_FLAG;		
//	newton and broyden method
extern int fout_nlin_final_FLAG;			
extern int fout_nlin_iter_FLAG;
extern int fout_nlin_final_NMIN;
//	broyden diag
extern int fout_xnext_broyden_diag_FLAG;
extern int fout_jxnext_broyden_FLAG;
//	fminBFGS and polytope in R(N) and fminR1 in R(1)
extern int fout_fmin_iter_FLAG;      
extern int fout_fmin_diag_FLAG;
extern int fout_fmin_final_FLAG;
extern int fout_update_hessian_fmin_FLAG;
extern const int fmin_BFGS_DIAG_FLAG; //?????
extern const int polytope_DIAG_FLAG; //??????
//	linemin newton and golden
extern int fout_linemin_iter_FLAG;
extern int fout_linemin_final_FLAG;
extern int FOUT_GET_GRADIENT_DIAG_FLAG;
extern int FOUT_NLIN_DIAG_FLAG;
extern int FOUT_BROYDEN_DIAG_FLAG;
//	misc output
extern int FOUT_GET_ERGODIC_FINAL_FLAG;
extern int FOUT_GET_ERGODIC_ITER_FLAG;
//	gradient and jacobian - numerical diff
extern int fout_get_d2fdx2_along_s_FLAG;
extern int fout_get_dfdx_along_s_FLAG;
extern int fout_get_gradient_FLAG;
//	calib_mom2mom
extern int fout_calib_mom2mom_diag_FLAG;
//	mat inv 
extern int fout_swap_row_FLAG;
//	misc identifiers and flags  ------------------------------------------------------------
extern double HELL;
extern int DIFF_3PT_FLAG;
extern int DIFF_5PT_FLAG;
extern int fminBFGS_FLAG;
extern int update_x0_FLAG;

extern int Nquad;
extern const int N_fout_nlin;             //  no. iterations between nlin solver loops for fout

//extern const int Ncalibration_moments;




//const int Nfmin = 2;

//  fmin.h - fmin_BFGS ---------------------------------------------------




//  diagnostics flags for math.h    (=1 for diagnostics turned on)


//  diagnostics flags for fmin.h (=1 for diagnostics turned on)

//  diagnostics flags for nlin_solvers.h (=1 for diagnostics turned on)
extern int newton_method_DIAG_FLAG;
extern int broyden_method_DIAG_FLAG;
extern int fout_nlin_labelled_eqn_resids_FLAG;		//	calls fcn to write out eqn resids if ON

// diagnostics flags for newton_calib solvers
extern int FOUT_CALIB_STEPWISE_FINAL_FLAG;

void toggle_fout_nlin_FLAGS();
