extern int fout_calib_mom2mom_iter_FLAG;
extern int mom2mom_step_geometric_FLAG;
extern int fout_get_optimal_step_mom2mom_FLAG;


//	calibration function given benchmark (moments, parameters, eqm) and target moments - this function uses fmin() as the inner loop
int calib_mom2mom(double *momBM, double *paramBM, double *xeqmBM,
				  double *momXP, double *paramXP, double *xeqmXP,
				  double (*fct_calib)(double *p), int (*constraint_set)(double *p),
				  void (*mom_extract)(double *m), double *gxeqm,
				  int Nmom, int Nxeqm, int Nparam, ParamFmin pm_fmin, double percent_step);
//	obsolete
int calib_mom2mom(double *momBM, double *paramBM, double *xeqmBM,
				  double *momXP, double *paramXP, double *xeqmXP,
				  double (*fct_calib)(double *p), int (*constraint_set)(double *p),
				  void (*mom_extract)(double *m), double *gxeqm,
				  int Nmom, int Nxeqm, int Nparam, int Nsteps,
				  int Nmaxrecr, int& jrecr, ParamFmin pm_fmin, double percent_Nstep);
//	calibration using newton_method() for inner loop
int calib_mom2mom(double *momBM, double *paramBM, double *xeqmBM,
				  double *momXP, double *paramXP, double *xeqmXP,
				  void (*fct_calib)(double *fp, double *p), int (*constraint_set)(double *p),
				  void (*mom_extract)(double *m), double *gxeqm,
				  int Nmom, int Nxeqm, int Nparam, ParamNewton pm_newt, double percent_step);
//	writes diag for calib_mom2mom recursion
void fout_calib_mom2mom_recr_diag(ofstream& fout, double *momBM_tmp, double *paramBM_tmp,
								  double *xeqmBM_tmp, double *momXP_tmp,
								  int Nmom, int Nxeqm, int Nparam,
								  int Nsteps, int Nmaxrecr, int jrecr);
//	writes fct call summary upon exit
void fout_mom2mom_final(ofstream& fout, double *momBM_tmp, double *paramCur,
						double *xeqmCur, double *momCur, double *paramXP_tmp,
						double *xeqmXP_tmp, int Nmom, int Nxeqm, int Nparam,
						int Nsteps, double percent_step, double time_used);
//	fct call iter diag
void fout_mom2mom_iter(ofstream& fout, double *momBM, double *paramBM,
					   double *xeqmBM, double *momCur, double *paramCur,
					   double *xeqmCur, int Nmom, int Nxeqm, int Nparam,
					   int i_step, int Nsteps, double time_used);
void get_optimal_step(double *x_BM, double *x_XP, int Nx, double percent_Nstep,
					  double *dx, int& Nsteps);
void fout_get_optimal_step(ofstream& fout, double max_change, double min_change,
						   double Napprox, int Nsteps, double percent_Nstep,
						   double *x_BM, double *x_XP, double *change, double *dmom,
						   int Nmom);
//	given a benchmark set of parameters and its eqm, computes the eqm for another set of parameters by proceeding to new parameters in constant percentage changes (assumes all parameters are positive)
int eqm2eqm(double *param_bm, double *eqm_bm, double *param_target, double *eqm_target,
			double percent_step, int Nparam, int Neqm, void (*param_extract)(double *p),
			void (*fct)(double *fx, double *x), int (*fct_constraints)(double *x),
			ParamNewton pm_newt, ofstream& fout);
int eqm2eqm(double *param_bm, double *eqm_bm, double *param_target, double *eqm_target,
			double percent_step, int Nparam, int Neqm, void (*param_extract)(double *p),
			void (*fct)(double *fx, double *x), int (*fct_constraints)(double *x),
			ParamNewton pm_newt);
//	summary/diagnostics for eqm2eqm()
void fout_eqm2eqm_iter(ofstream& fout, double *param_bm, double *param_cur, double *param_target,
					   int Nparam, double *eqm_bm, double *eqm_cur, int Neqm, double timeused,
					   int i, int Nsteps);

void fout_eqm2eqm_iter_error(ofstream& fout, double *param_bm, double *param_cur,
							 double *param_target, int Nparam,
							 double *eqm_bm, double *eqm_cur, int Neqm, int Nsteps);
//	fct summary for eqm2eqm_solver
void fout_eqm2eqm_final(ofstream& fout, double *param_BM, double *xeqm_BM, double *param_XP,
						double *xeqm_XP, int Nparam_tmp, int Nx_tmp, int Nsteps, double timeused);
