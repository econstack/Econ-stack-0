#include <math.h>
#include <fstream>
#include <string>
#include <conio.h>
#include <iostream>
#include <iomanip>
#include <time.h>

using namespace std;

//	math library -------------------------------------------------------
#include "math_parameters.h"
#include "io_grunts.h"
#include "math_grunts.h"
#include "NlinSolvers.h"
#include "fmin.h"
#include "calibration.h"

//	globals ----------------------------------------------------------------
int fout_mom2mom_iter_FLAG = ON;
int mom2mom_step_geometric_FLAG = ON;	// obsolete
int fout_get_optimal_step_mom2mom_FLAG = ON;
int fout_eqm2eqm_iter_FLAG = ON;
int fout_eqm2eqm_final_FLAG = ON;

//	step-wise recursive calibration call. Given a benchmark triple (parameters, moments, eqm) denoted by BM, this function solves for the parameters and eqm associated with a target set of moments (denoted XP). This is done by continuation method approach where it proceeds from the benchmark moments to the target moments in equally spaced steps. At each step, it can call itself if the step is so big that it cannot solve for an eqm or cannot solve for the associated parameters. There MUST be a fct ptr to a fct that extracts moments into globals.  Also, the ptr gxeqm MUST be to a global that stores the eqm (index) so that it can be updated.
int calib_mom2mom(double *xmomBM, double *xpmBM, double *xeqmBM,
				  double *xmomXP, double *xpmXP, double *xeqmXP,
				  double (*fct_calib)(double *p), int (*constraint_set)(double *p),
				  void (*mom_extract)(double *m), double *gxeqm,
				  int Nmom, int Nxeqm, int Nparam, ParamFmin pm_fmin, double percent_step)
{
	double *xdmom; xdmom = new double [Nmom];		//	factor change for moments for each step
	double *xmomCur; xmomCur = new double [Nmom];	//	current (target) moments
	double *xpmCur; xpmCur = new double [Nparam];	//	current parameters 
	double *xeqmCur; xeqmCur = new double [Nxeqm];	//	current eqm 
	double *Hp;	Hp = new double [Nparam*Nparam];	// Hessian used in fmin
	double *xpm_star; xpm_star = new double [Nparam]; // solution for current target moments

	//	time use vars
	double time_used, starttime, stoptime, stoptimef, starttimef = clock();

	//	compute step size for moments
	int Nsteps;
	mat_add(xdmom, xmomXP, xmomBM, -1, Nmom);
	get_optimal_step(xmomBM, xmomXP, Nmom, percent_step, xdmom, Nsteps);

	//	set current moments and its associated parameters and eqm
	copy_array(xmomCur, xmomBM, Nmom);
	copy_array(xpmCur, xpmBM, Nparam);
	copy_array(xeqmCur, xeqmBM, Nxeqm);

	//	make copy of existing global eqm and copy benchmark eqm into global
	double *gxeqm_old; gxeqm_old = new double [Nxeqm];
	copy_array(gxeqm_old, gxeqm, Nxeqm);
	copy_array(gxeqm, xeqmBM, Nxeqm);

	//	misc vars in main loop
	int fmin_flag = 0; // flag for whether calibrated successfully
	double step_newton_LM_old;
	double loss_fp = 0;

	//	loop over steps
	for (int i=1; i<=Nsteps; ++i) {

		//	start iter clock
		starttime = clock();
		//	set current moments
		mat_mult(xmomCur, xdmom, xmomCur, Nmom);
		//	extract current moments into globals
		mom_extract(xmomCur);

		//	calibrate to current moments using previous calibrated params as initial guess and also use its associated eqm as the initial guess for current eqm - speed up newton step for this part
		fmin_flag = fmin_BFGS(xpmCur, xpm_star, Nparam, Hp, fct_calib, constraint_set, pm_fmin);
		//	test for failure of fmin 
		if (fmin_flag == NOT_LOCALLY_CONVEX_FLAG) {
			log_fout << endl << "Exiting calib_mom2mom(). Current candidate not locally convex!";
			break;
		}

		//	store info on results of calibration
		loss_fp = pm_fmin.loss;
		//	update parameters and eqm for initial guess
		copy_array(xeqmCur, gxeqm, Nxeqm);
		copy_array(xpmCur, xpm_star, Nparam);

		//	restart again from solution to ensure that it is a local min
		log_fout << endl << "calib_mom2mom(): Restarting fmin from solution:";
		//	accelerate - increase step_newton in linemin
		step_newton_LM_old = pm_fmin.step_newton_LM;
		pm_fmin.step_newton_LM *= 2.0;
		fout_array(log_fout, xpm_star, 1, Nparam, 6);
		fmin_flag = fmin_BFGS(xpm_star, xpm_star, Nparam, Hp, fct_calib,
			constraint_set, pm_fmin);

		//	if restarted solution is better, then save it; o/w keep old params and info for next calibration step
		double loss_tmp = pm_fmin.loss;
		if ( ((fmin_flag > 0) || (fmin_flag == -3)) && (pm_fmin.loss < loss_fp) ) {
			//	update parameters and eqm for initial guess
			copy_array(xeqmCur, gxeqm, Nxeqm);
			copy_array(xpmCur, xpm_star, Nparam);
		}
		//	restore old newton_step in linemin
		pm_fmin.step_newton_LM = step_newton_LM_old;

		//	write out eqm and moments if needed
		if (fout_mom2mom_iter_FLAG == ON) {
			stoptime = clock();
			time_used = (stoptime - starttime)/1000;
			fout_mom2mom_iter(log_fout, xmomBM, xpmBM, xeqmBM, xmomCur, xpmCur, xeqmCur,
				Nmom, Nxeqm, Nparam, i, Nsteps, time_used);
		}

	}
	//	At final step, start again from solution but decrease the tolerance for movement in candidate and smaller steps along linemin_newton.
	pm_fmin.step_newton_LM *= 0.5;
	pm_fmin.eps_LM *= 0.05;
	fmin_flag = fmin_BFGS(xpm_star, xpm_star, Nparam, Hp, fct_calib,
		constraint_set, pm_fmin);

	//	return solution for parameters and its associated eqm
	copy_array(xpmXP, xpm_star, Nparam);
	copy_array(xeqmXP, gxeqm, Nxeqm);

	//	write out summary
	stoptimef = clock();
	time_used = (stoptimef - starttimef)/1000;
	fout_mom2mom_final(log_fout, xmomBM, xpmBM, xeqmBM, xmomXP, xpm_star, xeqmXP,
		Nmom, Nxeqm, Nparam, Nsteps, percent_step, time_used);
	log_fout << endl << "Exiting: calib_mom2mom()" << endl;

	//	restore previous global eqm
	copy_array(gxeqm, gxeqm_old, Nxeqm);
	//	free mem
	delete[] xdmom; delete[] xmomCur; delete[] xpmCur; delete[] xeqmCur;
	delete[] Hp; delete[] xpm_star; delete[] gxeqm_old;
	return fmin_flag;
}
//	newton_method() version: The older versions used fmin(). This one implements step-wise calibration using modified newton_method() to zero the system of equations characterizing the calibrated model. The modified newton_method() checks if the total loss given by the sum of squares of the residuals from feval(calib). If it cannot locate a point
//	old function call
int calib_mom2mom(double *momBM, double *paramBM, double *xeqmBM,
				  double *momXP, double *paramXP, double *xeqmXP,
				  double (*fct_calib)(double *p), int (*constraint_set)(double *p),
				  void (*mom_extract)(double *m), double *gxeqm,
				  int Nmom, int Nxeqm, int Nparam, int Nsteps,
				  int Nmaxrecr, int& jrecr, ParamFmin pm_fmin, double percent_step)
{
	log_fout << endl << "OBSOLETE FCT CALL : calib_mom2mom() *************************" << endl;

	int flag = calib_mom2mom(momBM, paramBM, xeqmBM, momXP, paramXP, xeqmXP,
		fct_calib, constraint_set, mom_extract, gxeqm, Nmom, Nxeqm, Nparam, pm_fmin,
		percent_step);
	return flag;
}
//	This version is similar except uses newton_method() instead of fmin() at each calibration step
int calib_mom2mom(double *momBM, double *paramBM, double *xeqmBM,
				  double *momXP, double *paramXP, double *xeqmXP,
				  void (*fct_calib)(double *fp, double *p), int (*constraint_set)(double *p),
				  void (*mom_extract)(double *m), double *gxeqm,
				  int Nmom, int Nxeqm, int Nparam, ParamNewton pmn, double percent_step)
{
	//	time use vars
	double time_used, starttime, stoptime, stoptimef, starttimef = clock();
	//	compute step size for moments
	double *dmom; dmom = new double [Nmom];
	mat_add(dmom, momXP, momBM, -1, Nmom);
	int Nsteps;
	get_optimal_step(momBM, momXP, Nmom, percent_step, dmom, Nsteps);

	//	current moments and its associated parameters and eqm
	double *xmomCur; xmomCur = new double [Nmom];
	copy_array(xmomCur, momBM, Nmom);
	double *xparamCur; xparamCur = new double [Nparam];
	copy_array(xparamCur, paramBM, Nparam);
	double *xeqmCur; xeqmCur = new double [Nxeqm];
	copy_array(xeqmCur, xeqmBM, Nxeqm);
	//	model parameters solution to current moments
	double *xparamstar; xparamstar = new double [Nparam];
	copy_array(xparamstar, xparamCur, Nparam);

	//	make copy of existing global eqm and copy benchmark eqm into global
	double *gxeqm_old; gxeqm_old = new double [Nxeqm];
	copy_array(gxeqm_old, gxeqm, Nxeqm);
	copy_array(gxeqm, xeqmBM, Nxeqm);

	//	misc vars in main loop
	int process_flag  = 0; // flag for whether calibrated successfully
	double *paramstar;	paramstar = new double [Nparam];
//	double step_newton_old;
//	double loss_fp;

	//	loop over steps
	for (int i=1; i<=Nsteps; ++i) {

		//	start iter clock
		starttime = clock();
		//	set current moments
		mat_mult(xmomCur, dmom, xmomCur, Nmom);
		//	extract current moments into globals
		mom_extract(xmomCur);

		//	calibrate to current moments using previous calibrated params as initial guess and also use its associated eqm as the initial guess for current eqm
		log_fout << endl << "calib_mom2mom(): main newton loop, " << i << "/"
				<< Nsteps << endl;
        log_fout << setw(20) << "xparamstar";
        fout_array(log_fout, xparamstar, 1, Nparam, 4, 0);

		process_flag = newton_method(xparamstar, Nparam, fct_calib, constraint_set, pmn);
        log_fout << endl << "exiting newton" << endl;

		//	if unsucessful, break
		if (process_flag < 0) {
			log_fout << endl << "calib_mom2mom(newton): Error, step " << i << "/"
				<< Nsteps << endl;
			break;
		}
		//	if successful, update model parameters for init guess and associated eqm
		else {
			//	update model parameters and eqm for initial guess
//			log_fout << endl << "calib_mom2mom(newton): Success, step " << i << "/" << Nsteps << endl;
			copy_array(xeqmCur, gxeqm, Nxeqm);
			copy_array(xparamCur, xparamstar, Nparam);
		}

		//	write out eqm and moments if needed
		if (fout_mom2mom_iter_FLAG == ON) {
			stoptime = clock();
			time_used = (stoptime - starttime)/1000;
			fout_mom2mom_iter(log_fout, momBM, paramBM, xeqmBM, xmomCur, xparamCur, xeqmCur,
				Nmom, Nxeqm, Nparam, i, Nsteps, time_used);
		}

	}

	//	update solution if successful
	if (process_flag > 0) {
		copy_array(momXP, xmomCur, Nmom);
		copy_array(xeqmXP, xeqmCur, Nxeqm);
		copy_array(paramXP, xparamCur, Nparam);
	}
	//	write out summary
	stoptimef = clock();
	time_used = (stoptimef - starttimef)/1000;
	fout_mom2mom_final(log_fout, momBM, paramBM, xeqmBM, momXP, paramXP, xeqmXP,
		Nmom, Nxeqm, Nparam, Nsteps, percent_step, time_used);

	//	restore previous global eqm
	copy_array(gxeqm, gxeqm_old, Nxeqm);

	//	free mem
	delete[] dmom; delete[] xmomCur; delete[] xparamCur; delete[] xeqmCur;
	delete[] paramstar; delete[] gxeqm_old; delete[] xparamstar;

	return process_flag;
}
//	recursion diagnostics (upon entry) for calib_mom2mom()
void fout_calib_mom2mom_recr_diag(ofstream& fout, double *momBM_tmp, double *paramBM_tmp,
								  double *xeqmBM_tmp, double *momXP_tmp,
								  int Nmom, int Nxeqm, int Nparam,
								  int Nsteps, int Nmaxrecr, int jrecr)
{
	fout << endl << "calib_mom2mom_recursion_diag: " << endl;
	fout << "N(params)=" << Nparam << ", N(moments)=" << Nmom << ", N(eqm)=" << Nxeqm
		<< ", recursion depth=" << jrecr << " of Max level " << Nmaxrecr << endl;

	fout << setw(12) << "Benchmarks";
	fout << setw(12) << "moments";
	fout_array(fout, momBM_tmp, 1, Nmom, 6, 0);
	fout << setw(24) << "params";
	fout_array(fout, paramBM_tmp, 1, Nparam, 6, 0);
	fout << setw(24) << "eqm";
	fout_array(fout, xeqmBM_tmp, 1, Nxeqm, 6, 0);

	fout << setw(12) << "Target";
	fout << setw(12) << "moments";
	fout_array(fout, momXP_tmp, 1, Nmom, 6, 0);

	return;
}
//	iteration diagnostics/summary for calib_mom2mom()
void fout_mom2mom_iter(ofstream& fout, double *momBM, double *paramBM,
					   double *xeqmBM, double *momCur, double *paramCur,
					   double *xeqmCur, int Nmom, int Nxeqm, int Nparam,
					   int i_step, int Nsteps, double time_used)
{
	int Nprec = 6;
	fout << endl << "calib_mom2mom_iter: step " << i_step << "/" << Nsteps
		<< ", time_used=" << time_used << endl;

	fout << setw(12) << "Benchmarks" << setw(12) << "moments";
	fout_array(fout, momBM, 1, Nmom, Nprec, 0);

	fout << setw(24) << "parameters";
	fout_array(fout, paramBM, 1, Nparam, Nprec, 0);

	fout << setw(24) << "eqm";
	fout_array(fout, xeqmBM, 1, Nxeqm, Nprec, 0);

	fout << setw(12) << "Current" << setw(12) << "moments";
	fout_array(fout, momCur, 1, Nmom, Nprec, 0);

	fout << setw(24) << "parameters";
	fout_array(fout, paramCur, 1, Nparam, Nprec, 0);

	fout << setw(24) << "eqm";
	fout_array(fout, xeqmCur, 1, Nxeqm, Nprec, 0);

	return;
}
void fout_mom2mom_final(ofstream& fout, double *momBM, double *paramBM,
						double *xeqmBM, double *momXP, double *paramXP,
						double *xeqmXP, int Nmom, int Nxeqm, int Nparam,
						int Nsteps, double percent_step, double time_used)
{
	int Nprec = 6;
	//	header
	log_fout << endl << "calib_mom2mom_final: time_used=" << time_used
		<< ", #steps=" << Nsteps << ", step_size=" << percent_step << endl;

	//	benchmark data
	fout << setw(12) << "Benchmarks" << setw(12) << "moments";
	fout_array(fout, momBM, 1, Nmom, Nprec, 0);

	fout << setw(24) << "parameters";
	fout_array(fout, paramBM, 1, Nparam, Nprec, 0);

	fout << setw(24) << "eqm";
	fout_array(fout, xeqmBM, 1, Nxeqm, Nprec, 0);

	//	target data
	fout << setw(12) << "Target" << setw(12) << "moments";
	fout_array(fout, momXP, 1, Nmom, Nprec, 0);

	fout << setw(24) << "parameters";
	fout_array(fout, paramXP, 1, Nparam, Nprec, 0);

	fout << setw(24) << "eqm";
	fout_array(fout, xeqmXP, 1, Nxeqm, Nprec, 0);
	return;
}

void get_optimal_step(double *x_BM, double *x_XP, int Nx, double percent_Nstep,
					  double *dx, int& Nsteps)
{
	double *change; change = new double [Nx];

	//	find factor of change between benchmark and target
	for (int i=0; i<Nx; ++i)
		*(change+i) = *(x_XP+i)/ *(x_BM+i);
	//	find largest factor change
	double max_change = get_max(change, Nx);	//	incr in moments
	double min_change = get_min(change, Nx);	//	decr in moments
	//	case biggest change is increasing
	double Napprox;
	if (max_change > 1/min_change)
		Napprox = log(max_change)/log(1+percent_Nstep);
	else //	case biggest change is decreasing
		Napprox = log(min_change)/log(1-percent_Nstep);
	//	find the number of steps - round up to next integer
	Nsteps = ( int ) Napprox + 1;
	//	given number of steps - find factor change per step
	for (int i=0; i<Nx; ++i)
		*(dx+i) = pow(*(change+i), 1.0/Nsteps);
	//	write out diag
	if (fout_get_optimal_step_mom2mom_FLAG == ON) {
		fout_get_optimal_step(log_fout, max_change, min_change, Napprox, Nsteps,
			percent_Nstep, x_BM, x_XP, change, dx, Nx);
	}
	//	free mem
	delete[] change;
	return;
}
void fout_get_optimal_step(ofstream& fout, double max_change, double min_change,
						   double Napprox, int Nsteps, double percent_Nstep,
						   double *x_BM, double *x_XP, double *change, double *dmom,
						   int Nmom)
{
	//	write header
	fout << endl << "get_optimal_step_diag : max change(+) = " << max_change
		<< ", max_change(-)=" << min_change << ", Napprox=" << Napprox << ", Nsteps="
		<< Nsteps << ", percent_change=" << percent_Nstep << endl;

	fout << setw(12) << "x(BM)=";
	fout_array(fout, x_BM, 1, Nmom, 6, 0);
	fout << setw(12) << "x(XP)=";
	fout_array(fout, x_XP, 1, Nmom, 6, 0);
	log_fout << setw(12) << "factor";
	fout_array(log_fout, change, 1, Nmom, 6, 0);
	log_fout << setw(12) << "step";
	fout_array(log_fout, dmom, 1, Nmom, 6, 0);
	return;
}
//	given a benchmark set of parameters and its eqm, computes the eqm for another set of parameters by proceeding to new parameters in constant percentage changes (assumes all parameters are positive)
//  added parameters: (i) output file option
int eqm2eqm(double *param_bm, double *eqm_bm, double *param_target, double *eqm_target,
			double percent_step, int Nparam, int Neqm, void (*param_extract)(double *p),
			void (*fct)(double *fx, double *x), int (*fct_constraints)(double *x),
			ParamNewton pm_newt, ofstream& fout)
{
	//	compute number of steps given percent change per step
	double *dparam; dparam = new double [Nparam];
	int Nsteps;
	get_optimal_step(param_bm, param_target, Nparam, percent_step, dparam, Nsteps);

	//	current moments and its associated parameters and eqm
	double *param_cur; param_cur = new double [Nparam];
	copy_array(param_cur, param_bm, Nparam);
	double *eqm_cur; eqm_cur = new double [Neqm];
	copy_array(eqm_cur, eqm_bm, Neqm);

	//	control vars
	double timeused, starttime, stoptime, stoptimef, starttimef = clock();
	int process_flag = 0;

	//	loop over steps
	for (int i=1; i<=Nsteps; ++i) {

		//	start iter clock
		starttime = clock();

		//	compute current parameters
		mat_mult(param_cur, dparam, param_cur, Nparam);
		param_extract(param_cur);

		//	compute eqm for current parameters, note that at each iteration, solution to previous iteration is used as candidate (init guess) for current iteration
		process_flag = newton_method(eqm_cur, Neqm, fct, fct_constraints, pm_newt);

		//	option: write out diagnostics for current step
		if (fout_eqm2eqm_iter_FLAG == ON) {
			stoptime = clock();
			timeused = (stoptime-starttime)/1000;
			fout_eqm2eqm_iter(fout, param_bm, param_cur, param_target, Nparam,
				eqm_bm, eqm_cur, Neqm, timeused, i, Nsteps);
		}

		//	case: cannot solve for eqm at some intermediate step, write out diag and exit
		if (process_flag < 0) {
			fout_eqm2eqm_iter_error(fout, param_bm, param_cur, param_target,
				Nparam, eqm_bm, eqm_cur, Neqm, Nsteps);
			break;
		}

	}

	copy_array(eqm_target, eqm_cur, Neqm);

	//	option: write out diagnostics
	if ((process_flag > 0) && (fout_eqm2eqm_final_FLAG == ON)) {
		stoptimef = clock();
		timeused = (stoptimef - starttimef)/1000;
		fout_eqm2eqm_final(fout, param_bm, eqm_bm, param_target, eqm_target, Nparam,
			Neqm, Nsteps, timeused);
	}
	//	case: failure at some step - copy current parameters to target parametres and current eqm candidate t
	if (process_flag < 0) {
		copy_array(param_target, param_cur, Nparam);
		copy_array(eqm_target, eqm_cur, Neqm);
	}
	//	free mem
	delete[] dparam; delete[] eqm_cur; delete[] param_cur;
	return process_flag;
}
//  old function call
int eqm2eqm(double *param_bm, double *eqm_bm, double *param_target, double *eqm_target,
			double percent_step, int Nparam, int Neqm, void (*param_extract)(double *p),
			void (*fct)(double *fx, double *x), int (*fct_constraints)(double *x),
			ParamNewton pm_newt)
{
    int flag = eqm2eqm(param_bm, eqm_bm, param_target, eqm_target,
        percent_step, Nparam, Neqm, param_extract, fct, fct_constraints, pm_newt, log_fout);

    return flag;
}
//	summary/diagnostics for eqm2eqm()
void fout_eqm2eqm_iter(ofstream& fout, double *param_bm, double *param_cur, double *param_target,
					   int Nparam, double *eqm_bm, double *eqm_cur, int Neqm, double timeused,
					   int i, int Nsteps)
{
	int Nprec = 4;

	fout << "eqm2eqm_iter_diag: " << i << "/" << Nsteps << ", timeused=" << timeused;

	fout << endl << setw(15) << "param(bm)";
	fout_array(fout, param_bm, 1, Nparam, Nprec, 0);

	fout << setw(15) << "param(current)";
	fout_array(fout, param_cur, 1, Nparam, Nprec, 0);

	fout << setw(15) << "param(target)";
	fout_array(fout, param_target, 1, Nparam, Nprec, 0);

	fout << setw(15) << "eqm(bm)";
	fout_array(fout, eqm_bm, 1, Neqm, Nprec, 0);

	fout << setw(15) << "eqm(current)";
	fout_array(fout, eqm_cur, 1, Neqm, Nprec, 0);

	return;
}

void fout_eqm2eqm_iter_error(ofstream& fout, double *param_bm, double *param_cur,
							 double *param_target, int Nparam,
							 double *eqm_bm, double *eqm_cur, int Neqm, int Nsteps)
{
	int Nprec = 4;

	fout << endl << "eqm2eqm_error: Nsteps=" << Nsteps;

	fout << endl << setw(15) << "param(bm)";
	fout_array(fout, param_bm, 1, Nparam, Nprec, 0);

	fout << setw(15) << "param(target)";
	fout_array(fout, param_target, 1, Nparam, Nprec, 0);

	fout << setw(15) << "param(current)";
	fout_array(fout, param_cur, 1, Nparam, Nprec, 0);

	fout << setw(15) << "eqm(bm)";
	fout_array(fout, eqm_bm, 1, Neqm, Nprec, 0);

	fout << setw(15) << "eqm(current)";
	fout_array(fout, eqm_cur, 1, Neqm, Nprec, 0);

	return;
}

//	fct summary for eqm2eqm_solver
void fout_eqm2eqm_final(ofstream& fout, double *param_BM, double *xeqm_BM, double *param_XP,
						double *xeqm_XP, int Nparam_tmp, int Nx_tmp, int Nsteps, double timeused)
{
	//	write header
	fout << endl << "eqm2eqm_solver summary: Nsteps=" << Nsteps << ", timeused=" <<
		timeused << endl;
	//	write out benchmark
	fout << setw(12) << "params(BM)";
	fout_array(fout, param_BM, 1, Nparam_tmp, 6, 0);
	fout << setw(12) << "xstar(BM)";
	fout_array(fout, xeqm_BM, 1, Nx_tmp, 6, 0);
	//	write out XP
	fout << setw(12) << "params(XP)";
	fout_array(fout, param_XP, 1, Nparam_tmp, 6, 0);
	fout << setw(12) << "xstar(XP)";
	fout_array(fout, xeqm_XP, 1, Nx_tmp, 6, 0);

	return;
}
