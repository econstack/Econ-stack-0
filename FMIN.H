class ParamFmin {
public:
	double step_jacob;		//	step size for numerical Jacobian
	double eps_fmin;		//	tolerance for fmin
	double step_newton_LM;	//	scalar for step size in line_min call
	double eps_LM;			//	tolerance for line_min call
	double percentstep_LM;	//	max change in percentage terms per step in line_min call
	int Nmax_iter_LM;		//	max # iterations in line_min call
	int Nmax_iter_fmin;		//	max # iterations for fmin 
	int linemin_flag;		//	switch for type of line_min (golden section vs newton)
	int flag_fout_final;
	int Nfout_fmin;			//	# iterations between write to log for fmin
	int Nfout_LM;			//	# iterations between write to log in line_min call
	double loss;
	//	member function ===============================================================
	ParamFmin(); // constructor and initialization
	//	TODO write to outfile 

	// use default copy and "=" constructors
};
//	used by linemin_newton and contract_linemin_newton_step
extern int NOT_LOCALLY_CONVEX_FLAG;
extern int Nmaxiter_linemin_xnext;
extern int fout_contract_linemin_newton_step_FLAG;
/*	Revised April 24, 2006. Added option to limit step size in each iteration. Set to zero to turn off. */
int linemin_newton(double *x, double *dir, double *xmin, int Nx_tmp,
				   double& fx, double& df, double& d2f,
				   double (*fct)(double *x), int (*fct_constraints)(double *x),
				   ParamFmin pm_fmin);
/*	Version: Aug 7, 2005. Solves constrained minimization problem along given direction. Finds the minimum *x of a function *fct along a direction (one-dimensional subsp) *dir using Newton's method, where updating is: x(k+1) = x(k) - dir*(f'(x(k))/f"(x(k))). The derivatives are approximated by numerical differentiation. Returns #iterations if converged, o/w returns non-convergence flag. Also returns the 1st and 2nd derivative at xmin along dir. Nfout_tmp = #iterations between write-to-file diagnostics. Nmaxiter_tmp = max # iteration. eps_tmp = convergence tolerance. */
int linemin_newton(double *x, double *dir, double *xmin, int Nx_tmp, double& fx, 
				   double& dfdx, double& d2fdx, double (*fct)(double *x),
				   int (*fct_constraints)(double *x), int Nfout_tmp, int Nmaxiter_tmp,
				   double step_newton_tmp = linemin_stepNewton, 
				   double step_jacob = stepJacob, double eps_tmp = epsLinemin,
				   int flag_fout_final = OFF);
double linemin_newton(double *x, double *dir, double *xmin, int Nx,
					  double (*fct)(double *x), double step_jacob_tmp = linemin_stepJacob);
/*	Finds the local min near x0 in R(Nx) for a function *fct along a line in R(Nx) using bracketing. The function's value at the local min will be returned. Given function *fct, candidate x0, and search dir, the min is stored in xmin. */
double linemin_golden(double *x0, double *dir, double *xmin, int Nx,
					  double (*fct)(double *x));
//	constrained minimization problem in R1
int fmin_R1(double x, double& xmin, double df, double d2f, double loss_fx,	
			double (*fct)(double x), int (*fct_constraints)(double x), ParamFmin prm);
//	old fct call in fmin_v3.h and older
int fmin_R1(double x, double& xmin, double step_newton_tmp, double df, double d2f,
			double loss_fx,	int Nfout, double (*fct)(double x), 
			int (*fct_constraints)(double x), double step_df, 
			double epsFmin_tmp = epsFmin, int Nmax_iter_tmp = NmaxiterFmin);
/* version: Sept 14, 2005. Solves constrained minimization problem. Finds the min of a function (*fct) from R(N) to R using the BFGS method which uses gradient info. The gradient is computed numerically. The subrountine for finding a minimum along a direction is set by the linemin_flag. Latest change is the bundling of computation parameters into struct ParamFmin. Note that current version makes obsolete older versions of the function call as default values are now built in. */
int fmin_BFGS(double *x, double *xmin, const int Nx, double *Hx, 
			  double (*fct)(double *x), int (*fct_constraints)(double *x),
			  ParamFmin& prm);
//	old fct call from fmin_v3.h
int fmin_BFGS(double *x, double *xmin, const int Nx, double step_BFGS, 
			  double *Hx, int Nfout, double (*fct)(double *x), 
			  int (*fct_constraints)(double *x),
			  double step_jacob = stepJacob, double eps_tmp = epsFmin,
			  double step_newton_LM = linemin_stepNewton,
			  double eps_LM = epsLinemin, int Nmax_iter_LM = NmaxiterLinemin,
			  int Nmax_iter_tmp = NmaxiterFmin, int linemin_flag = LINEMIN_NEWTON_FLAG,
			  int flag_fout_final = OFF);
double fmin_BFGS(double *x0 /*initial_candidate*/, int Nx, 
				 double *xmin /*solution*/, double (*fct)(double *x), 
				 double step_jacob_tmp = stepJacob,
				 double epsFmin_tmp = epsFmin,
				 int Nmax_iter_tmp = NmaxiterFmin,
				 int linemin_flag = LINEMIN_NEWTON_FLAG);
double fmin_BFGS(double *x0 /*initial_candidate*/, int Nx, 
				 double *xmin /*solution*/, double (*fct)(double *x), 
				 int linemin_flag = LINEMIN_NEWTON_FLAG);
//	fmin with DFP updating for approx hessian
int fmin_DFP(double *x, double *xmin, const int Nx, double *Hx, 
			double (*fct)(double *x), int (*fct_constraints)(double *x), ParamFmin prm);
/*  Implements Powell's conjugate direction method. First initialize all direction u(i) = e(i) for i \in {1,Nx} where e(i) are the basis vectors in R(Nx). Then while the stopping condition is not satisfied: (step-i) Save starting point x(0), (step-ii) for i = 1,..,Nx, move x(i-1) to the minimum along direction u(i) and call this point x(i), (step-iii) for i = 1,..,(Nx-1), set u(i) to u(i+1) (step-iv) set u(Nx) to x(Nx)-x(0, and (step-v) move x(Nx) to the minimum along direction u(Nx) and call this point p(0) */
int fmin_Powell(double *x, int Nx, 
				void (*fct)(double *fx, double *x), 
				double (*fct_constraints)(double *x));
/*	created: feb 24, 2005. This version implements a more efficient version for costly function evaluations where the function keeps track of the changes in vertices and does not re-evaluate the function at a vertice known from previous iteration. */
void update_gradient_fminBFGS(double *y, double *z, double *x, double *xnext,
							  double *Jx, double *Jxnext, double& fx, 
							  double (*fct)(double *x), int Nx);
/* fmin subfunctions. Compute (i) H(k)z(k)z(k)H(k) and (ii) z(k)H(k)z(k) where H1(k) = (i) and H1a(k) = H(k)z(k), H1b(k) = z(k)H(k) and H2 = z(k)H(k)z(k). They are calculated as: H1b = transpose(H1a), H2 = H1b(k)z(k). Note that H2(k) is a scalar. Compute (iii) y(k)y(k) and (iv) y(k)z(k) where H3(k) = (iii) and H4(k) = (iv) note that H4(k) is a scalar */
void update_hessian_fminBFGS(double *Hx, double *Hxnext, double *x,
							 double *y, double *z, int Nx);
void update_hessian_fminBFGS_NR(double *Hx, double *Hxnext, double *x, double *y,
								 double *z, int Nx);
void update_hessian_fminDFP(double *Hx, double *Hxnext, double *x, double *y, double *z, int Nx);
//	linemin_newton subfunctions 
int contract_linemin_newton_step(double& fx, double *x, double *xnext, double *dir, int Nx_tmp, 
								 double (*fct)(double *x), double step_jacob, double& df, 
								 double& d2f, double& fnext, double& dfnext, double& d2fnext,
								 double step_newton_tmp, int DIFF_3PT_FLAG);
//	--------------------------- DIAGNOSTICS --------------------------------------------------
void fout_update_hessian_fmin(ofstream& fout, double *Hpos, double *Hneg, 
							  double *Hx, double *Hxnext, int Nx);
