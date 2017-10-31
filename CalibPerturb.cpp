#include <math.h>
#include <fstream>
#include <string>
#include <conio.h>
#include <iostream>
#include <iomanip>
#include <time.h>

using namespace std;

#include "math_parameters.h"
#include "io_grunts.h"
#include "io_diag.h"
#include "math_grunts.h"
#include "NlinSolvers.h"
#include "Fmin.h"

#include "CalibPerturb.h"

//	=============================================================================================
//	Constructor 
CalibPerturb::CalibPerturb()
{
	Nmom = 1;
	Nparam = 1;
	percent_step = 0.1;
	Nsteps = 1;
	solver_flag = 0;
}
//	write to outfile
void CalibPerturb::fout(ofstream &fout)
{	
	//	write out parameters specific to CalibPerturb
	fout << endl << setw(20) << "CalibPerturb():" 
		<< setw(12) << "flag" << setw(6) << solver_flag
		<< setw(12) << "#mom" << setw(6) << Nmom
		<< setw(12) << "#param" << setw(6) << Nparam
		<< setw(12) << "%step" << setw(6) << percent_step
		<< setw(12) << "#steps" << setw(6) << Nsteps << endl;

	//	TODO - create fmin.fout() and write out fmin parameters
	
	return;	
}
//	computes the # steps for calibration by perturbation method()
void CalibPerturb::get_optimal_step(double *xmom_bm, double *xmom_xp, double *dxmom, ofstream *fout)
{
	double *factor_change; factor_change = new double [this->Nmom];

	//	find factor of change between benchmark and target
	for (int i=0; i<(this->Nmom); ++i)
		factor_change[i] = xmom_xp[i] / xmom_bm[i];

	//	find largest factor change
	double max_change = get_max( factor_change, this->Nmom );	//	incr in moments
	double min_change = get_min( factor_change, this->Nmom );	//	decr in moments
	
	//	case biggest change is increasing
	double Napprox;
	if ( max_change > (1.0/min_change) )
		Napprox = log(max_change) / log(1.0+this->percent_step);
	else //	case biggest change is decreasing
		Napprox = log(min_change)/log(1.0-this->percent_step);
	//	find the number of steps - round up to next integer
	Nsteps = (int) Napprox + 1;
	//	given number of steps - find factor change per step
	for (int i=0; i<(this->Nmom); ++i)
		dxmom[i] = pow( factor_change[i], 1.0/Nsteps);
	//	option - write to outfile
	if ( fout != 0 ) {
		this->get_optimal_step_fout( *fout, max_change, min_change, Napprox, xmom_bm, xmom_xp, 
			factor_change, dxmom);
	}
	//	free heap
	delete[] factor_change;
	return;
}
void CalibPerturb::get_optimal_step_fout(std::ofstream &fout, const double& max_change, const double& min_change, 
										 double Napprox, double *xmom_bm, double *xmom_xp, 
										 double *factor_change, double *dxmom)
{
	//	set precision
	int Nprec = 4;
	fout << setiosflags(ios::fixed) << setprecision(Nprec);
	//	write header
	fout << endl << "get_optimal_step(): max change(+)= " << max_change
		<< ", max_change(-)=" << min_change << ", Napprox=" << Napprox << ", Nsteps="
		<< this->Nsteps << ", percent_change=" << this->percent_step << endl;
	//	write out data
	fout << setw(12) << "x(BM)=";
	fout_array( fout, xmom_bm, 1, this->Nmom, Nprec, OFF );
	fout << setw(12) << "x(XP)=";
	fout_array( fout, xmom_xp, 1, this->Nmom, Nprec, OFF );
	log_fout << setw(12) << "factor";
	fout_array( log_fout, factor_change, 1, this->Nmom, Nprec, OFF );
	log_fout << setw(12) << "dstep";
	fout_array( log_fout, dxmom, 1, this->Nmom, Nprec, OFF );
	return;
}
//	
int CalibPerturb::solver(double *xpm_bm, double *xmom_bm, double *xpm_xp, double *xmom_xp, 
						 double (*fct_calib)(double *), int (*constraint_set)(double *), 
						 void (*mom_extract)(double *), ofstream *fout)
{	
	double *xdmom; xdmom = new double [Nmom];		//	factor change for moments for each step
	double *xmom_cur; xmom_cur = new double [Nmom];	//	current (target) moments
	double *xpm_cur; xpm_cur = new double [Nparam];	//	current parameters 
	double *xHp; xHp = new double [Nparam*Nparam];	// Hessian used in fmin
	double *xpm_star; xpm_star = new double [Nparam]; // solution for current target moments

	//	time use vars
	double time_used, starttime, stoptime, stoptimef, starttimef = clock();

	//	compute step size for each perturbation step
	this->get_optimal_step( xmom_bm, xmom_xp, xdmom );

	//	set current moments and its associated parameters and eqm
	copy_array(xmom_cur, xmom_bm, this->Nmom);
	copy_array(xpm_cur, xpm_bm, this->Nparam);

	//	misc vars in main loop 
	int fmin_flag = 0; // flag for whether calibrated successfully
	double step_newton_LM_old;
	double loss_fp = 0;

	//	loop over steps
	for (int i_cp=1; i_cp<=(this->Nsteps); ++i_cp) {

		//	start iter clock
		starttime = clock();
		//	set current moments
		mat_mult( xmom_cur, xdmom, xmom_cur, this->Nmom );
		//	extract current moments into globals
		mom_extract( xmom_cur );

		//	calibrate to current moments using previous calibrated params as initial guess and also use its associated eqm as the initial guess for current eqm - speed up newton step for this part
		fmin_flag = fmin_BFGS( xpm_cur, xpm_star, this->Nparam, xHp, fct_calib, constraint_set, pm_fmin);
		//	test for failure of fmin 
		if (fmin_flag == NOT_LOCALLY_CONVEX_FLAG) {
			log_fout << endl << "Exiting calib_mom2mom(). Current candidate not locally convex!";
			break;
		}

		//	store info on results of calibration
		loss_fp = pm_fmin.loss;

		//	update parameters 
		copy_array( xpm_cur, xpm_star, this->Nparam );

		//	restart again from solution to ensure that it is a local min
		log_fout << endl << "CalibPerturb_solver(): Restarting fmin from solution:";
		//	accelerate - increase step_newton in linemin and write summary to outfile
		step_newton_LM_old = pm_fmin.step_newton_LM;
		pm_fmin.step_newton_LM *= 2.0;
		log_fout << setw(12) << "xpm(star)=";
		fout_array( log_fout, xpm_star, 1, this->Nparam, 6, OFF );
		fmin_flag = fmin_BFGS(xpm_star, xpm_star, this->Nparam, xHp, fct_calib, constraint_set, pm_fmin);

		//	if restarted solution is better, then save it; o/w keep old params and info for next calibration step
		double loss_tmp = pm_fmin.loss;
		if ( ((fmin_flag > 0) || (fmin_flag == -3)) && (pm_fmin.loss < loss_fp) ) {
			//	update parameters for initial guess
			copy_array( xpm_cur, xpm_star, this->Nparam );
		}
		//	restore old newton_step in linemin
		pm_fmin.step_newton_LM = step_newton_LM_old;

		//	option - write current model parameters and moments to outfile
		if ( fout != 0 ) {
			stoptime = clock();
			time_used = (stoptime - starttime)/1000;
			this->fout_mom2mom_iter( *fout, xmom_bm, xpm_bm, xmom_cur, xpm_cur, i_cp, time_used );
		}

	}

	//	At final step, start again from solution but decrease the tolerance for movement in candidate and smaller steps along linemin_newton.
	pm_fmin.step_newton_LM *= 0.5;
	pm_fmin.eps_LM *= 0.10;
	fmin_flag = fmin_BFGS( xpm_star, xpm_star, this->Nparam, xHp, fct_calib, constraint_set, pm_fmin);

	//	return solution for parameters and its associated eqm
	copy_array( xpm_xp, xpm_star, this->Nparam );

	//	write out summary
	stoptimef = clock();
	time_used = (stoptimef - starttimef)/1000;
	this->fout_mom2mom_final( log_fout, xmom_bm, xpm_bm, xmom_xp, xpm_xp, time_used );
	log_fout << endl << "Exiting: CalibPerturb::solver()" << endl;

	//	free heap
	delete[] xdmom; delete[] xmom_cur; delete[] xpm_cur; 
	delete[] xHp; delete[] xpm_star;
	//	return with process flag
	return fmin_flag;
}
//	misc diagnostics for solver()
void CalibPerturb::fout_mom2mom_iter(std::ofstream &fout, double *xmom_bm, double *xpm_bm, 
									 double *xmom_cur, double *xpm_cur, int i_cp, double time_used)
{
	//	set precision 
	int Nprec = 6;
    fout << setiosflags(ios::fixed) << setprecision(Nprec);
	//	write header
	fout << endl << "CalibPerturb_solver: step " << i_cp << "/" << Nsteps
		<< ", time_used=" << time_used << endl;
	//	write date for baseline
	fout << setw(18) << "Benchmarks" << setw(12) << "moments";
	fout_array( fout, xmom_bm, 1, this->Nmom, Nprec, OFF );
	fout << setw(30) << "parameters";
	fout_array( fout, xpm_bm, 1, this->Nparam, Nprec, OFF );
	//	write data for calibrated 
	fout << setw(18) << "Current Target" << setw(12) << "moments";
	fout_array( fout, xmom_cur, 1, this->Nmom, Nprec, OFF );
	fout << setw(30) << "parameters";
	fout_array(fout, xpm_cur, 1, this->Nparam, Nprec, OFF );
	return;
}
void CalibPerturb::fout_mom2mom_final(std::ofstream &fout, double *xmom_bm, double *xpm_bm, 
									  double *xmom_xp, double *xpm_xp, double time_used)
{
	//	set precision 
	int Nprec = 6;
    fout << setiosflags(ios::fixed) << setprecision(Nprec);
	//	header
	log_fout << endl << "CalibPerturb_solver(): time_used=" << time_used
		<< ", #steps=" << this->Nsteps << ", step_size=" << this->percent_step << endl;
	//	benchmark data
	fout << setw(18) << "Baseline" << setw(12) << "moments";
	fout_array( fout, xmom_bm, 1, this->Nmom, Nprec, OFF );
	fout << setw(30) << "parameters";
	fout_array( fout, xpm_bm, 1, this->Nparam, Nprec, OFF );
	//	target data
	fout << setw(18) << "Target" << setw(12) << "moments";
	fout_array( fout, xmom_xp, 1, this->Nmom, Nprec, OFF );
	fout << setw(30) << "parameters";
	fout_array( fout, xpm_xp, 1, this->Nparam, Nprec, OFF );
	return;

}
//	=============================================================================================
