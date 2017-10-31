/*	------------- GENERIC SUMMARY AND ERROR DIAGNOSTICS ------------------------------				Inverse matrix error diagnostics */
void fout_inv_mat_error(ofstream& , int , int , double *, double *);
//	generic constraint set error diagnostics
void fout_xnext_error(ofstream& , double *, double *, int , string );
void fout_xnext_error(ofstream& , double *, double *, int );
//	generic exceeded max # iteration diagnostics
void fout_max_iter_error(ofstream& , int , double *, double *, int , string );
//	generic iterative solvers (nlin or fmin) function call summary if successful
void fout_iter_final_diag(ofstream& fout, double *x0_tmp, double *xstar_tmp, int Nx_tmp,
						  int jiter, int Nmaxiter_tmp, double time_used_tmp, double eps_tmp,
						  string label_tmp);
void fout_iter_final_diag(ofstream& , double *, double *, int , int , int , double , string );
//	generic iterative solvers (nlin or fmin) function exit if unsuccessful
void fout_iter_exit(ofstream& , double *, double *, int , int , int , string );
//	generic nlin iteration statistics
void fout_nlin_iter(ofstream& , int , int , double , double , double , double *,
					double *, double *, double *, string );
//	generic nlin_solver function call stats summary
void fout_nlin_final(ofstream& fout, int Nx, int process_flag, double loss_x,
					 double loss_fx, double *x0, double *x, double *fx, double time_used);
/*	------------- LINEMIN NEWTON DIAGNOSTICS AND ERRORS	---------------------------     linemin_newton iteration diagnostics (note final exit summary is not fct specific). */
void fout_linemin_newton_iter(ofstream& fout, int jiter, double timeused,
							  double dfdx, double d2fdx, double dfnext,
							  double d2fnext, double *dir, double *xnext,
							  double *x, double fx, double fxnext, int Nx, double loss_x);
void fout_linemin_newton_iter(ofstream& fout, int jiter, double timeused,
							  double dfdx, double d2fdx, double *dir,
							  double *xnext, double *x, double fx, int Nx,
							  double loss_x);
/*	-------------- LINEMIN GOLDEN DIAGNOSTICS AND ERRORS ----------------------------               Cannot find initial triple for golden section. */
void fout_linemin_golden_findtriple_error(ofstream& fout, double *x0, double *dir, int Nx);
/*	iteration and function call summary stats. */
void fout_linemin_golden_iter(ofstream& fout, int jiter, double timeused,
							  double fa, double *xa, double fc, double *xc, int N);
void fout_linemin_golden_final(ofstream& fout, int process_flag, double timeused,
							   double fa, double *xa, double fc, double *xc, int N);
/*	------------- GET DERIVATIVE AND GRADIENT DIAGNOSTICS ------------------------- */
void fout_get_dfdx_diag(ofstream& fout, double *dir_tmp, double fx, double fxabove,
						double fxbelow,	double *x, double *xabove, double *xbelow,
						double dfdx, int Nx, double step_jacob);
void fout_get_d2fdx2_diag(ofstream& fout, double *dir_tmp, double dfdx_above, double dfdx_below,
						  double *xabove, double *xbelow, int Nx, double d2fdx2,
						  double step_jacob, double dfdx);
void fout_get_d2fdx2_diag(ofstream& fout, double *sdir_tmp, double ffabove, double fabove,
						  double fx, double fbelow, double ffbelow, int Nx, double h,
						  double *xxabove, double *xabove, double *x,
						  double *xbelow, double *xxbelow, double d2fdx2, double dfdx);
void fout_get_gradient_diag(ofstream& fout, double *x, double *Jx, int Nx,
							double fx, double step_jacob);
void fout_get_gradient_diag(ofstream& fout, double *x, double *Jx, int Nx, double fx,
							double *fx_above, double *fx_below, double step_jacob);
/*	---------------- BROYDEN METHOD DIAGNOSTICS AND ERRORS -------------------------------- get_xnext_broyden general diagnostics */
void fout_xnext_broyden_diag(ofstream& fout, double step_broyden_tmp, double *Jx, double *Jxinv,
							 double *dir, double *x, double *xnext, double *fx, int Nx_tmp);
//	get_jxnext_broyden general diagnostics
void fout_jxnext_broyden_diag(ofstream& fout, double *Jx, double *Jxnext,
							  double *xnext, double *fx, double *fxnext,
							  double *y_temp, int Nx_tmp);
/*	------------------ FMIN IN R(1) DIAG AND ERRORS -------------------------------------    iteration stats for fmin_R1 */
void fout_fminR1_iter(ofstream& fout, int jfmin, int Nmax_iter_tmp, double x,
					  double xnext, double loss_x, double loss_fx, double loss_df,
					  double df, double d2f, double time_used);
//	fct call summary
void fout_fminR1_final(ofstream& fout, int process_flag, int Nmax_iter_tmp,
					   double loss_x, double loss_fx, double loss_df,
					   double time_used, double x0_tmp, double xmin);
/*	------------------- FMIN_BFGS DIAG AND ERRORS ----------------------------------------- iteration statistics for fminBFGS */
void fout_fminBFGS_iter(ofstream& fout, int jfmin, int Nmaxiter_tmp, double timeused,
						double fx, double fxnext, double *x, double *xnext, int N,
						double norm_Jxnext, double norm_dx, double eps_fmin);
void fout_fminBFGS_iter(ofstream& fout, int jfmin, int Nmaxiter_tmp, double timeused,
						double fx, double fxnext, double *x, double *xnext, int N,
						double norm_Jxnext, double norm_dx);
void fout_fminBFGS_iter(ofstream& fout, int jfmin, int Nmaxiter_tmp, double timeused,
						double fx, double fxnext, double *x, double *xnext, int N);
void fout_fminBFGS_final(ofstream& fout, double *x0, double fx0, double *xmin, double fxmin,
						 int Nx, int process_flag, int Nmax_iter_fmin, double time_used,
						 double eps_fmin);
//	general debugging diagnostics for fminBFGS
void fout_fminBFGS_diag(ofstream& fout, double *x, double *xnext, double *y, double *z,
						double *Jxnext, double *Hxnext, int N, int jfmin, int Nmaxiter_tmp,
						double loss_x, double loss_Jx);
//	suboptimal point error (non-zero jacobian at solution)
void fout_fminBFGS_badx_error(ofstream& fout, double *x, double *xnext, double *Jxnext, int Nx,
							  int jfmin, int Nmaxiter_tmp, double norm_x, double norm_Jx);
//	writes out jacobian and hessian (read in from file) on resumption of fminBFGS
void fout_fminBFGS_resume(ofstream& fout, double *Jx, double *Hx, double *x,
						  double *x_linemin, int Nx);
//	----------------------- POLYTOPE MINIMIZATION DIAG ------------------------------------------ iteration statistics
void fout_polytope_iter(ofstream& fout, double timeused, double *vertices,
						double *fvertices, int Nx);
//	fct call summary
void fout_polytope_final(ofstream& fout, double timeused, double *xmax,
						 double fxmax, int iter_ctr, int Nx);
/*	------------------ MISC -----------------------------------------------	*/
//	error in interpolation of 2 arrays st m(theta) = theta*m1 + (1-theta)*m2
void fout_mat_interpolation_error(ofstream& fout, double *m1, double *m2,
								  double theta_m1, int N);
//	constraint set call <rescale_xnext> summary
void fout_xnext_diag(ofstream& fout, double *x, double *xnext, double *xnext0, int Nx);
void fout_vf_final(ofstream& fout, int Niters, double loss, double time_used);
/*	Writes out the diagnostics for get_eqm() */
void fout_get_eqm_diag(ofstream& fout, double *x0, double *x, int Nx,
					   void (*fct)(double *fx, double *x));
/*	--------------------- ECON NLIN ----------------------------------------                        diagnostics for eqm2eqm recursion call */
void fout_eqm2eqm_recursion_diag(ofstream& fout, double *param_BM_tmp, double *x0_BM_tmp,
								 double *param_XP_tmp, int Nparam_tmp, int Nx_tmp);
