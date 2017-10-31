class ParamLinemin {
public:
	int Nx;						//	dimension of system dim(x) 
	double step_newton_linemin;	//	scalar for step size in linemin() call
	double eps_linemin;			//	tolerance for linemin() call
	double percentstep_linemin;	//	max change in percentage terms per step in line_min call
	int Nmax_iter_linemin;		//	max # iterations in line_min call
	int Nfout_linemin;			//	# iterations between write to log in line_min call
	double step_jacob;			//	step size for numerical differentation 
	double flag_locally_convex;	//	=ON for breaking if cannot find locally convex candidate
	//	golden section parameters
	double golden_alpha;		//  scaling parameter in seeking init triple
	double golden_delta;		//  step size in seeking init triple
	int golden_Nmaxiter_triple; //  max no. iter in seeking init triple
	//	process diagnostics
	double loss_x;
	double fx;
	double dfx;
	double d2fx;
	double time_used;
	int process_flag;
	//	member functions =================================================================
	ParamLinemin();	// constructor
	//	write computation parameters to outfile
	void _fout(ofstream& fout);
	//	solvers: searches for min along given direction *dir 
	int _solver_newton(double *x, double *xmin, double *dir, double& fx, double& df, double& d2f,
					   double (*fct)(double *x), int (*fct_constraints)(double *x), ofstream *fout = 0);
	int _solver_golden(double *x, double *xmin, double *dir, double &fx, 
					   double (*fct)(double *x), int (*fct_constraints)(double *x), ofstream *fout = 0);
	//	solver utils 
	int _solver_golden_get_init_triple(double *x0, double *dir, double *xa, double *xb, double *xc, 
		  							   double &fa, double &fb, double &fc, double (*fct)(double *), 
									   int (*fct_constraints)(double *), std::ofstream *fout = 0);
	//	write out diagnostics and summary data
	void _fout_iter(ofstream &fout, double *x, double *xnext, const int& i_lmin, const double& fx, 
					const double& fxnext);
	void _fout_exit(ofstream &fout, double *x0, double *xstar, const double& fx0, 
					const double& fxstar);
	void _fout_error(ofstream &fout, double *x, double *xnext, const int& i_nlin, 
					 const double& fx, string label);
	//	checks that x(k+1) is within percentage of x(k), rescales x(k+1) if needed
	int _get_xnext_dpercent(double *xnext, double *x, ofstream *fout = 0);
	int _get_xnext_convex(double *xnext, double *x, double *dir, double& fnext, double &dfnext, 
						   double& d2fnext, double (*fct)(double *x), ofstream *fout = 0);
	//	use default copy and assign operators
	
};
class ParamFxmin {
public:
	ParamLinemin pm_lmin;		//	computation parameters for linemin
	double step_jacob;			//	step size for numerical Jacobian 
	double eps_fmin;			//	tolerance for fmin
	int Nmax_iter_fmin;			//	max # iterations for fmin 
	int Nfout_fmin;				//	# iterations between write to log for fmin
	int linemin_newton_flag;	//	switch for type of line_min (golden section vs newton)
	int BFGS_flag;				//	switch for type of solver and hessian update BFGS = default, DFP = alt
	//	process diagnostics
	double loss_x;
	double fx;
	int process_flag;
	double time_used;
	//	member function ====================================================================
	ParamFxmin(); // constructor and initialization
	//	write parameters to outfile
	void _fout(ofstream& fout);
	//	solvers : solves for min of function f given starting point x ----------------------
	int _solver_BFGS(double *x, double *xmin, double (*fct)(double *x), 
					 int (*fct_constraints)(double *x), ofstream *fout = 0);
	int _solver_DFP(double *x, double *xmin, double (*fct)(double *x), 
					int (*fct_constraints)(double *x), ofstream *fout = 0);
	int _solver_R1(const double& x0, double& xmin, double& fx, double& df, double& d2f,
				   double (*fct)(double x), int (*fct_constraints)(double x), ofstream *fout = 0);
	//	solver utils ----------------------------------------------------------------------
	void _solver_update_gradient(double *y, double *z, double *x, double *xnext, double *Jx, double *Jxnext,
								 double& fxnext, double (*fct)(double *x), ofstream *fout = 0);
	void _solver_update_hessian_BFGS(double *Hx, double *Hxnext, double *x, double *y, double *z, 
									 ofstream *fout = 0);
	void _solver_update_hessian_DFP(double *Hx, double *Hxnext, double *x, double *y, double *z, 
									ofstream *fout = 0);
	void _solver_update_hessian_fout(ofstream& fout, double *Hx, double *Hxnext);
	void _solver_R1_fout_iter(ofstream& fout, const double& x, const double& xnext, 
							  const double& fx, const double& fxnext, const int& j_fmin);
	void _solver_R1_fout_exit(ofstream& fout, const double& x0, const double& xmin, 
							  const double& fx0, const double& fxmin);
	//	write out diagnostics and summary data ---------------------------------------------
	void _fout_iter(ofstream &fout, double *x, double *xnext, const int& i_lmin, const double& fx, 
					const double& fxnext);
	void _fout_exit(ofstream &fout, double *x0, double *xstar, const double& fx0, 
					const double& fxstar);
	void _fout_error(ofstream &fout, double *x, double *xnext, const int& i_nlin, 
					 const double& fx, string label);
	//	misc utils --------------------------------------------------------------------


};
