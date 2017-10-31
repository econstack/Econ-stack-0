
class ParamNlin {
public:
	int Nf; // dimension (size) of system of equations
	double loss_fx; // L2 norm residual for |f(x*)| at solution x*
	double time_used; // total time used in call to _solver()
	double step_nlin; // size of broyden or newton step st x(n+1) = x(n) + s*J(n)^(-1) where s = step_newton 
	double eps_nlin; // stopping tolerance 
	double step_jacob; // zero for numerical differentiation (f(x+h)-f(x))/h where h = step_jacob
	double step_dpercent; // max percentage change (represented as decimal)
	int powell_flag; // ON = use powell hybrid method, check if (fx)^2 decreases from x(n) to x(n+1)
	int Nmaxiter; // max # of iterations
	int Nfout; // # of iterations between write to outfile 
	int flag; // main process flag, returns # iters if successful or (negative) error code
	int step_dpercent_flag; // ON = max change from x(n) to x(n+1) of step_dpercent
	//	member functions =========================================================================
	ParamNlin(); // constructor - init to default values (non-zero)
	//	write parameters to outfile
	void _fout(ofstream& fout);
	/*	------------------------------------------------------------------------------------------
		Solves for zero of a system of equations given an initial guess as double *x using 
		Newton method or Broyden method
	---------------------------------------------------------------------------------------------- */
	void _newton_method(double *x, void (*fct)(double *fx, double *x), int (*fct_constraints)(double *x),
						ofstream *fout = 0);
	void _broyden_method(double *x, void (*fct)(double *fx, double *x), int (*fct_constraints)(double *x),
						 ofstream *fout = 0);
	/*	-----------------------------------------------------------------------------------------
		Broyden_method utilites:
		(i) get_jxnext() : update Jacobian in broyden_method: 
			Set y(k) = f(x(k+1)) - f(x(k)) and set J(k+1) = J(k) + {(y(k)-J(k)s(k))s(k)'}/{s(k)'s(k)} 
			where J(k+1) = Jxnext, f(x(k+1)) = fxnext, y(k) = y_temp 
		(ii) get_xnext() : gets next candidate x(k+1) in broyden_method, also (i) gets the 
			inverse jacobian, Jinv(k), (ii) evaluates the eqns at x, f(k), and (iii) computes 
			the search direction, dir(k)
	--------------------------------------------------------------------------------------------- */
	void _broyden_get_jxnext(double *fx, double *fxnext, double *xnext, double *y_temp, 
							 double *dir, double *Jx, double *Jxnext, 
							void (*fct)(double *fx, double *x), ofstream *fout = 0);
	int _broyden_get_xnext(double *Jx, double *Jxinv, double *fx, double *dir, double *x, 
						   double *xnext, void (*fct)(double *fx, double *x), ofstream *fout = 0);
	/*	-------------------------------------------------------------------------------------------
		Iteration utilies:
		(i) get_xnext_dpercent(): checks if x(k+1) is within step_dpercent% of x(k), i.e
			|x(k+1)-x(k)| < step_percent 
		(ii) get_xnext_powell(): implements powell hybrid method where x(k+1) must be st the SSR 
		decreases, ie 
			|f(k+1)| < |f(k)|
		In both cases (i) and (ii), contracts the candidate x(k+1) back toward x(k) if possible 
	--------------------------------------------------------------------------------------------- */
	//	TODO write member fucntion
	int _get_xnext_dpercent(double *xnext, double *x, ofstream *fout = 0);
	void _get_xnext_dpercent_fout(ofstream& fout, double *x, double *xnext, double *dx, double *dxpercent,
								  const double& max_dxpercent, const int& j_max_dxpercent);
	//	TODO write member function
	int _get_xnext_powell(double *xnext, double *x, double *fx, double *fxnext, 
						  void (*fct)(double *fx, double *x), ofstream *fout = 0);
	//	Iteration diagnostics - write to outfile --------------------------------------------------
	void _fout_iter(ofstream& fout, const int& i_nlin, const double& loss_x, const double& loss_fx,
					const double& time_used, double *x, double *xnext, double *fx, double *fxnext);
	//	Exit solver diagnostics - write to outfile
	void _fout_exit(ofstream& fout, const int& i_nlin, const double& loss_fx, 
					double *x0, double *xstar, double *fx0, double *fxstar);
	//	Generic diagnostic for error - write to outfile
	void _fout_error(ofstream& fout, const int& i_nlin, const double& loss_fx,
					 double *x, double *xnext, double *fx, const string& label_tmp);
	// use default copy and "=" constructors
};