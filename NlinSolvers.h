//	computation parameters for newton_method 
class ParamNewton {
public:
	int Nf; // dimension (size) of system of equations
	double loss; // L2 norm residual for |f(x*)| at solution x*
	double time_used; // total time used in call to _solver()
	double step_newton; // size of newton step st x(n+1) = x(n) + s*J(n)^(-1) where s = step_newton
	double eps_newton; // stopping tolerance 
	double step_jacob; // zero for numerical differentiation (f(x+h)-f(x))/h where h = step_jacob
	double step_dpercent; // max percentage change (represented as decimal)
	int powell_flag; // ON = use powell hybrid method, check if (fx)^2 decreases from x(n) to x(n+1)
	int newton_flag; //	(=0) for newton_method and (!=0) for broyden method
	int Nmaxiter; // max # of iterations
	int Nfout; // # of iterations between write to outfile 
	int step_dpercent_flag; // ON = max change from x(n) to x(n+1) of step_dpercent
	//	member functions =========================================================================
	ParamNewton(); // constructor - init to default values (non-zero)
	// use default copy and "=" constructors
};

/*  TEMPLATE. Solves a system of equations of dim Nx by Newton's method. Function pointers to the system of  equations (*fct) and to the system of constraints (*fct_constraints) must be passed */
int newton_method(double *x, int Nx, void (*fct)(double *fx, double *x),  
				  int (*fct_constraints)(double *x), ParamNewton& pm_newt);
int newton_method(double *x, int Nx, void (*fct)(double *fx, double *x),  
				  int (*fct_constraints)(double *x),
				  double& loss, int Nmaxiter_tmp, double step_newton, 
				  int Nfout_tmp, double &time_used, 
				  double eps_newton = eps_Nlin, double step_jacob = stepJacob);
/*No constraints version*/
int newton_method(double *x, int Nx, void (*fct)(double *fx, double *x), 
				  double& loss, int Nmaxiter_tmp, double step_newton, 
				  int Nfout_tmp, double &time_used, double eps_newton = eps_Nlin);
int newton_method(double *x, int Nx, void (*fct)(double *fx, double *x), 
				  int (*fct_constraints)(double *x) /*default version*/);
//	newton_method() utilities --------------------------------------------------------
int get_xnext_newton_powell(double *xnext, double *x, double *dx, double *fx, 
							double *fxnext, int Nx, void (*fct)(double *fx, double *x));
/* TEMPLATE. Solves a set of nonlinear eqns in R(Nx) using broyden method. Returns a flag for the number of iterations taken (=1 if failure); returns info on the error term st loss = sup_norm(x) and Nmaxiter_tmp = max number of iterations, Nfout_tmp = no iterations between fout to log file and to console, time_used = total time used*/
int broyden_method(double *x /*solution*/, int Nx,	void (*fct)(double *fx, double *x),		
				   int (*fct_constraints)(double *x), double& loss /*residual at x*/,	
				   int Nmaxiter_tmp, int Nfout_tmp, double& time_used,
				   double step_broyden_tmp = step_Broyden,
				   double eps_broyden_tmp = eps_Nlin);
/* No constrainst version*/
int broyden_method(double *x, int Nx, void (*fct)(double *fx, double *x), 
				   double &loss, int Nmaxiter_tmp, int Nfout_tmp, double& time_used);
int broyden_method(double *x, int Nx, void (*fct)(double *fx, double *x),
				   int (*fct_constraints)(double *x) /*default version*/);
//	solves for the eqm of a system as a fixed point of a nlin set of eqns
int get_eqm(double *x0, double *x, double *fx, void (*fct)(double *fx, double *x), 
			int Nx_tmp, int Nx0_tmp, int (*fct_constraints)(double *x),
			double& loss, double& time_used, int Nmaxiter_tmp, int Nfout, 
			double step_newton, double eps_newton, double step_jacob = stepJacob);
/*	 update Jacobian in broyden_method: Set y(k) = f(x(k+1)) - f(x(k)) and set J(k+1) = J(k) + {(y(k)-J(k)s(k))s(k)'}/{s(k)'s(k)} where J(k+1) = Jxnext, f(x(k+1)) = fxnext, y(k) = y_temp */
void get_jxnext_broyden(double *fx, double *fxnext, double *xnext, double *y_temp, 
						double *dir, double *Jx, double *Jxnext, int Nx_tmp, 
						void (*fct)(double *fx, double *x));
/*	gets next candidate x(k+1) in broyden_method, also (i) gets the inverse jacobian, Jinv(k), (ii) evaluates the eqns at x, f(k), and (iii) computes the search direction, dir(k). */
int get_xnext_broyden(double *Jx, double *Jxinv, double *fx, double *dir, double *x, 
					  double *xnext, void (*fct)(double *fx, double *x), int Nx_tmp, 
					  double step_broyden_tmp = step_Broyden);
void get_xnext_newton_dpercent(double *xnext, double *x, int Nx, double step_newton);
//	------------------------ DIAGNOSTICS -----------------------------------------------------
void fout_xnext_newton_dpercent(ofstream& fout, double *x, double *xnext, double *xnext_old, 
								double *dx_percent, double max_x_dpercent, int jmax_x_dpercent, 
								double max_dpercent, int Nx);
void fout_get_xnext_newton_powell(ofstream& fout, double *x, double *xnext, double *xnext_tmp, 
								  int Nx, double loss_fx, double loss_fxnext, double loss_fxnext_tmp);
void fout_ParamNewton(ofstream& fout, ParamNewton prm);

