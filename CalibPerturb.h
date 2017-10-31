
/*
	Calibrates a model using perturbation methods where a baseline eqm and its associated parameters and moments are known. Proceeds in geometrically equal steps to target moments to be calibrated. 
*/
class CalibPerturb {
public:
	//	computation parameters
	ParamFmin pm_fmin;		//	computation parameters for inner call to fmin()
	int Nmom;				//	# moments in calibration
	int Nparam;				//	# parameters to be calibrated
	double percent_step;	//	step size as fraction from baseline moments to target moments 
	int Nsteps;				//	# steps needed of size percent_step between baseline and target moments
	int solver_flag;		//	returns # iters if successful; else returns negative error code
	//	member functions ====================================================
	//	constructor
	CalibPerturb();
	//	write to outfile
	void fout(ofstream& fout);
	//	compute # steps to be used in solver
	void get_optimal_step(double *xmom_bm, double *xmom_xp, double *dxmom, ofstream *fout = 0);
	void get_optimal_step_fout(std::ofstream &fout, const double& max_change, const double& min_change, 
		double Napprox, double *xmom_bm, double *xmom_xp, double *factor_change, double *dxmom);
	//	main solver 
	int solver(double *xpm_bm, double *xmom_bm, double *xpm_xp, double *xmom_xp, 
			   double (*fct_calib)(double *p), int (*constraint_set)(double *p),
			   void (*mom_extract)(double *m), ofstream *fout = 0);
	//	misc diagnostics 
	void fout_mom2mom_iter(ofstream& fout, double *xmom_bm, double *xpm_bm, double *xmom_cur, 
		double *xpm_cur, int i_cp, double time_used);
	void fout_mom2mom_final(ofstream& fout, double *xmom_bm, double *xpm_bm, double *xmom_xp, 
		double *xpm_xp, double time_used);
};