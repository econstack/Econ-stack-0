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
#include "math_grunts.h"
#include "io_diag.h"
#include "NlinSolvers.h"
#include "fmin.h"

int Nmaxiter_linemin_xnext = 5;
int fout_contract_linemin_newton_step_FLAG = OFF;

//	initializes ParamFmin to default parameters
ParamFmin::ParamFmin()
{
	step_jacob = stepJacob;
	eps_fmin = 0.001;
	step_newton_LM = 0.5;
	eps_LM = 0.005;
	percentstep_LM = 0.05;
	Nmax_iter_LM = 20;
	Nmax_iter_fmin = 20;
	linemin_flag = LINEMIN_NEWTON_FLAG;
	flag_fout_final = ON;
	Nfout_fmin = 0; // no iteration output
	Nfout_LM = 0; // no iteration output
	loss = -1; // sets it to a number that cannot arise when used in function call
}
void fout_ParamFmin(ofstream& fout, ParamFmin pm)
{
	fout << endl << setw(20) << " -- ParamFmin -- " << endl;
	fout << setw(20) << "step_jacob" << setw(12) << pm.step_jacob << endl;
	fout << setw(20) << "eps_fmin" << setw(12) << pm.eps_fmin << endl;
	fout << setw(20) << "step_newton_LM" << setw(12) << pm.step_newton_LM << endl;
	fout << setw(20) << "eps_LM" << setw(12) << pm.eps_LM << endl;
	fout << setw(20) << "Nmax_iter_LM" << setw(12) << pm.Nmax_iter_LM << endl;
	fout << setw(20) << "Nmax_iter_fmin" << setw(12) << pm.Nmax_iter_fmin << endl;
	fout << setw(20) << "linemin_flag" << setw(12) << pm.linemin_flag << endl;
	fout << setw(20) << "Nfout_fmin" << setw(12) << pm.Nfout_fmin << endl;
	fout << setw(20) << "Nfout_LM" << setw(12) << pm.Nfout_LM << endl;
	return;
}
/*	Revised April 24, 2006. Added option to limit step size in each iteration. Set to zero to turn off. */
int linemin_newton(double *x, double *dir, double *xmin, int Nx_tmp,
				   double& fx, double& df, double& d2f,
				   double (*fct)(double *x), int (*fct_constraints)(double *x),
				   ParamFmin pm_fmin)
{
	//	hotfix : prevent crashes from restarting fmin in calib_mom2mom in case of not locally convex errors
	copy_array(xmin, x, Nx_tmp);
	//	extract computation parameters into locals
	int Nfout_tmp = pm_fmin.Nfout_LM;
	int Nmaxiter_tmp = pm_fmin.Nmax_iter_LM;
	double step_newton_tmp = pm_fmin.step_newton_LM;
	double step_jacob = pm_fmin.step_jacob;
	double eps_tmp = pm_fmin.eps_LM;
	int flag_fout_final = pm_fmin.flag_fout_final;
	
	//	candidate solutions
	double *xnext;		xnext = new double [Nx_tmp];
	double *x0_tmp;		x0_tmp = new double [Nx_tmp];
	//	properties at candidate solution
	double dfnext, d2fnext, fxnext = 1000;
	double fx0; 
	//	make copy of init guess
	copy_array(x0_tmp, x, Nx_tmp);
	//	scale dir into unit vector
	double norm_dir_tmp = norm_L2(dir, Nx_tmp);
	mat_rescale(dir, dir, 1/norm_dir_tmp, Nx_tmp);
	//	local control-flow vars
	int i = 0; // current number of iterations
	int xnext_flag, flag_rescale, process_flag = 0;
	double loss_x = 1;
	//	time use vars
	double stoptime_tmp, starttime_tmp, stoptimef_tmp, timeused = 0;
	double starttimef_tmp = clock();

	//	option : write search data to log 
	if ((fout_linemin_iter_FLAG==ON)||(Nfout_tmp==0)) {
		log_fout << endl << "linemin_newton():";
		log_fout << endl << "  x= "; fout_array(log_fout, x, 1, Nx_tmp, 6, 0);
		log_fout << "dir= "; fout_array(log_fout, dir, 1, Nx_tmp, 6, 0);
	}
	//	compute first derivative and 2nd derivative along dir
	d2f = get_d2fdx2_along_s(fx0, x, dir, Nx_tmp, fct, step_jacob, df, DIFF_3PT_FLAG);
	fx = fx0;

	//	main loop
	while ((loss_x > 0.10*eps_tmp*fabs(1.0+norm_L2(x,Nx_tmp))) && (fabs(df) > eps_tmp)) {

		//	dirty error flagging
		if (d2f < 0) {
			log_fout << endl << "linemin error!! not locally convex!!" << endl;
			process_flag = NOT_LOCALLY_CONVEX_FLAG;
			break;
		}
		//	start iteration clock
		starttime_tmp = clock();
		//	update candidate for xmin
		mat_add(xnext, x, dir, -step_newton_tmp*df/d2f, Nx_tmp);
		//	option: check if x(k+1) is within percent change limit, contract step if needed
		if (pm_fmin.percentstep_LM != 0) {
			get_xnext_newton_dpercent(xnext, x, Nx_tmp, pm_fmin.percentstep_LM);
		}
		//	rescale candidate so that it is in constraint set if needed
		flag_rescale = rescale_xnext(x, xnext, Nx_tmp, fct_constraints);
		if (flag_rescale<0) {
			fout_xnext_error(log_fout, x, xnext, Nx_tmp, "linemin_newton");
			process_flag = CONSTRAINT_SET_VIOLATION_FLAG;
			break;
		}
		//	check for local convexity at next candidate, break if unable to find locally convex xnext
		d2fnext = get_d2fdx2_along_s(fxnext, xnext, dir, Nx_tmp, fct, step_jacob,
			dfnext, DIFF_3PT_FLAG);
		if (d2fnext < 0) {
			//	find locally convex xnext and get f, df, d2f at new candidate, xnext
			xnext_flag = contract_linemin_newton_step(fx, x, xnext, dir, Nx_tmp, fct, step_jacob,
				df, d2f, fxnext, dfnext, d2fnext, step_newton_tmp, DIFF_3PT_FLAG);
			//	break if unable to rescale next guess st it is locally convex
			if (xnext_flag < 0) {
//				fout_linemin_error_nonconvexity(log_fout, x, xnext, Nx_tmp, df, d2f, fxnext, dfnext);
				log_fout << endl << "not locally convex!!!!!" << endl;
				cout << endl << "not locally convex!!!" << endl;
				process_flag = NOT_LOCALLY_CONVEX_FLAG;
				break;
			}
		}

		//	compute change in x in normed space
		loss_x = norm_sup(x, xnext, Nx_tmp);
		//	write out iteration diagnostics
		if (Nfout_tmp!=0) { //	iteration diagnostics
			if ((i%Nfout_tmp)==0)  {
				stoptime_tmp = clock();
				timeused = (stoptime_tmp - starttime_tmp) / 1000;
				log_fout << endl << "linemin: step_newton= " << step_newton_tmp << ", |dx|=" << loss_x;
				fout_linemin_newton_iter(log_fout, i, timeused, df, d2f, dfnext, d2fnext,
					dir, xnext, x, fx, fxnext, Nx_tmp, loss_x);
			}
		}
		if ((fout_linemin_iter_FLAG==ON)&&(Nfout_tmp==0)) {
			stoptime_tmp = clock();
			timeused = (stoptime_tmp - starttime_tmp) / 1000;
			fout_linemin_newton_iter(log_fout, i, timeused, df, d2f, dfnext, d2fnext,
					dir, xnext, x, fx, fxnext, Nx_tmp, loss_x);
		}
		//	check for convergence
		if (loss_x < 0.10*eps_tmp*fabs(1.0 + norm_L2(x, Nx_tmp)))	{
			stoptimef_tmp = clock();
			timeused = (stoptimef_tmp - starttimef_tmp) / 1000;
			copy_array(xmin, xnext, Nx_tmp);
			//	successful convergence
			if (fabs(df)<eps_tmp) { process_flag = i; }
			//	x stuck but non-zero derivative
			else { process_flag = X_STUCK_FLAG; }
			break;
		}

		// increment iteration counter and update guess and function at new guess
		++i;
		copy_array(x, xnext, Nx_tmp);
		df = dfnext;
		d2f = d2fnext;
		fx = fxnext;
		//	Break if exceeded max number of iterations
		if (i >= Nmaxiter_tmp) {
			fout_max_iter_error(log_fout, Nmaxiter_tmp, xnext, x0_tmp, Nx_tmp,
				"line_min_newton");
			copy_array(xmin, xnext, Nx_tmp);
			process_flag = EXCEEDED_MAXITER_FLAG;
			break;
		}

	}

	//	update and return process_flag and return solution xmin
	if ((process_flag >= 0) || (process_flag == X_STUCK_FLAG)) {
		process_flag = i;
		if (fx0 > fxnext) {
			copy_array(xmin, xnext, Nx_tmp);
			fx = fxnext;
			df = dfnext;
			d2f = d2fnext;
		} else {
			copy_array(xmin, x0_tmp, Nx_tmp);
			fx = fx0;
		}
	}
	//	write out function call summary
	if ((flag_fout_final == ON)||(fout_linemin_final_FLAG == ON)) {
		fout_iter_final_diag(log_fout, x0_tmp, xmin, Nx_tmp, process_flag, Nmaxiter_tmp,
			timeused, eps_tmp, "line_min_newton");
		loss_x = norm_L2(x0_tmp, xmin, Nx_tmp);
		log_fout << "|dx|=" << loss_x << ", eps(LM)=" << eps_tmp;
	}
	//	restore original init guess
	copy_array(x, x0_tmp, Nx_tmp);
	delete[] xnext; delete[] x0_tmp;
	return process_flag;
}
/*	Version: Aug 7, 2005. Solves constrained minimization problem along given direction. Finds the minimum *x of a function *fct along a direction (one-dimensional subsp) *dir using Newton's method, where updating is: x(k+1) = x(k) - dir*(f'(x(k))/f"(x(k))). The derivatives are approximated by numerical differentiation. Returns #iterations if converged, o/w returns non-convergence flag. Also returns the 1st and 2nd derivative at xmin along dir. Nfout_tmp = #iterations between write-to-file diagnostics. Nmaxiter_tmp = max # iteration. eps_tmp = convergence tolerance. */

int linemin_newton(double *x, double *dir, double *xmin, int Nx_tmp,
				   double& fx, double& df, double& d2f,
				   double (*fct)(double *x), int (*fct_constraints)(double *x),
				   int Nfout_tmp, int Nmaxiter_tmp, double step_newton_tmp,
				   double step_jacob, double eps_tmp, int flag_fout_final)
{
	ParamFmin pm_fmin;
	pm_fmin.Nfout_LM = Nfout_tmp;
	pm_fmin.Nmax_iter_LM = Nmaxiter_tmp;
	pm_fmin.step_newton_LM = step_newton_tmp;
	pm_fmin.step_jacob = step_jacob;
	pm_fmin.eps_LM = eps_tmp;
	pm_fmin.flag_fout_final = flag_fout_final;
	//pm_fmin.percentstep_LM = 0; // turn on step size constraint

	int flag = linemin_newton(x, dir, xmin, Nx_tmp, fx, df, d2f, fct, fct_constraints, pm_fmin);
	return flag;
}
double linemin_newton(double *x, double *dir, double *xmin, int Nx,
					  double (*fct)(double *x), double step_jacob_tmp)
{
	double fx, dfdx, d2fdx;
	//	set computation parameters for linemin_newton
	int Nfout_tmp = 0;
	int Nmaxiter_tmp = NmaxiterLinemin;
	double eps_tmp = epsLinemin;

	int process_flag = linemin_newton(x, dir, xmin, Nx, fx, dfdx, d2fdx, fct, no_constraints,
		Nfout_tmp, Nmaxiter_tmp, step_jacob_tmp, eps_tmp, OFF);

	return fx;
}
/*	Finds the local min near x0 in R(Nx) for a function *fct along a line in R(Nx) using bracketing. The function's value at the local min will be returned. Given function *fct, candidate x0, and search dir, the min is stored in xmin. */
double linemin_golden(double *x0, double *dir, double *xmin, int Nx,
					  double (*fct)(double *x))
{
	double loss = 1;
	double *xa;		xa = new double [Nx];
	double *xb;		xb = new double [Nx];
	double *xc;		xc = new double [Nx];
	double *xd;		xd = new double [Nx];

	double fa, fb, fc, fd;

	int NOT_FOUND = 1;
	int find_triple_flag = NOT_FOUND;
	int process_flag = 0;
	int jiter = 0;
	int d_LESSTHAN_b = 1;
	int d_GREATERTHAN_b = 0;
	int d_flag;			//	flag for d<b or d>b for candiate next point d

	//	rescale dir into unit vector
	double rescale_dir = norm_L2(dir,Nx);
	mat_rescale(dir,dir,1/rescale_dir,Nx);

	//	find points xa<xc<xb st f(xa), f(xb) > f(xc) ------------------
	copy_array(xb,x0,Nx);
	fb = fct(xb);
	int jfind = 0;

	double starttime, stoptime, stoptimef, starttimef = clock();

	while (find_triple_flag==NOT_FOUND) {

		mat_add(xa,xb,dir,pow(alpha_linemin,jfind)*delta_linemin,Nx);
		mat_add(xc,xb,dir,-pow(alpha_linemin,jfind)*delta_linemin,Nx);
		fa = fct(xa);
		fc = fct(xc);

		//	case f(a), f(c) > f(b)
		if ((fa>fb)&&(fc>fb)) {
			find_triple_flag = 1 - NOT_FOUND;
			break;
		}	//	end if()

		++jfind;

		if (jfind == NmaxiterLM_triple) {
			fout_linemin_golden_findtriple_error(log_fout,x0,dir,Nx);
			process_flag = FAILURE_FLAG;
			break;
		}

	}	//	end while()


	//	main loop ------------------------------------
	while ((process_flag!=FAILURE_FLAG)&&(loss>epsLinemin))	{

		if (fout_fmin_iter_FLAG ==1) { starttime = clock(); }

		//	choose next candidate xd and evaluate f(xd)
		if (norm_L2(xb,xa,Nx)<norm_L2(xc,xb,Nx))	{	//	set d = 0.5*(b+c)
			double temp1 = norm_L2(xb,xa,Nx);
			double temp2 = norm_L2(xc,xb,Nx);
			mat_add(xd,xc,xb,1.0,Nx);
			mat_rescale(xd,xd,0.5,Nx);
			d_flag = d_GREATERTHAN_b;
		}
		else {									//	set d = 0.5*(a+b)
			mat_add(xd,xa,xb,1.0,Nx);
			mat_rescale(xd,xd,0.5,Nx);
			d_flag = d_LESSTHAN_b;
		}

		//	set new triple (xa,xb,xc)
		fd = fct(xd);


		if ((d_flag==d_LESSTHAN_b)&&(fd>=fb)) {		//	case: replace (a,b,c) with (d,b,c)
			copy_array(xa,xd,Nx);
			fa = fd;
		}
		else if ((d_flag==d_LESSTHAN_b)&&(fd<fb)) { //	case: replace (a,b,c) with (a,d,b)
			copy_array(xc,xb,Nx);
			copy_array(xb,xd,Nx);
			fc = fb;
			fb = fd;
		}
		else if ((d_flag==d_GREATERTHAN_b)&&(fd<=fb)) { //case: replace (a,b,c) with (b,d,c)
			copy_array(xa,xb,Nx);
			copy_array(xb,xd,Nx);
			fa = fb;
			fb = fd;
		}
		else {										//	case: replace (a,b,c) with (a,b,d)
			copy_array(xc,xd,Nx);
			fc = fd;
		}


		//	check stopping condition
		loss = norm_L2(xc,xa,Nx);
		loss += pow(fabs(fa-fc),2);

		++jiter;
		if (jiter == NmaxiterLinemin) {
			log_fout << endl << "Linemin - exceed max number of iterations " << NmaxiterLinemin
				<< " Exiting...";
			log_fout << endl << "f(x_lower)=" << fa << ", x_lower=";
			fout_array(log_fout,xa,1,Nx,4,0);
			log_fout << "f(x_upper)=" << fc << ", x_upper=";
			fout_array(log_fout,xc,1,Nx,4,0);
			process_flag = FAILURE_FLAG;
			break;
		}

		if (fout_linemin_iter_FLAG == ON)	{	//	write out iteration diag
			stoptime = clock();
			double timeused = (stoptime - starttime) / 1000;
			fout_linemin_golden_iter(log_fout,jiter,timeused,fa,xa,fc,xc,Nx);
		}

	}	//	end while()

	if (process_flag!=FAILURE_FLAG) process_flag = jiter;

	if (fout_linemin_final_FLAG == ON) {	//	write out function diagnostics
		stoptimef = clock();
		double timeused = (stoptimef - starttimef) / 1000;
		fout_linemin_golden_final(log_fout,process_flag,timeused,fa,xa,fc,xc,Nx);
	}

	copy_array(xmin,xb,Nx);

	delete[] xa; delete[] xb; delete[] xc; delete[] xd;
	return loss;
}

//	constrained minimization problem in R1
int fmin_R1(double x, double& xmin, double df, double d2f, double loss_fx,
			double (*fct)(double x), int (*fct_constraints)(double x), ParamFmin prm)
{
	//	extract computation parameters
	int Nfout = prm.Nfout_fmin;
	int Nmax_iter_tmp = prm.Nmax_iter_LM;
	double step_newton_tmp = prm.step_newton_LM;
	double step_df = prm.step_jacob;
	double eps_fmin_tmp = prm.eps_LM;
	//	store and later restore global step_Jacob
	double step_jacob_old = stepJacob;
	stepJacob = step_df;

	int constraint_flag, process_flag = 0;
	//	init loss vars
	double loss_x, loss_df, loss = 1;
	//	init iteration counter
	int jfmin = 0;

	double fx, xnext;

	//	make copy of initial guess for diagnostics
	double x0_tmp = x;
	//	init clock
	double time_used, starttime, stoptime, stoptimef, starttimef = clock();

	while (loss>eps_fmin_tmp) {

		starttime = clock();

		//	step 1: get next candidate x(k+1)
		df = get_dfdx(fx, x, fct, step_df);
		d2f = get_d2fx2(fx, x, fct, step_df);
		xnext = x - step_newton_tmp*df/d2f;

		//	check if xnext in constraint set, rescale if needed
		constraint_flag = fct_constraints(xnext);
		if (constraint_flag>0)
			constraint_flag = rescale_xnext(x, xnext, fct_constraints);

		//	step 2: check stopping condition
		loss_fx = fct(x);
		loss_x = fabs(xnext - x);
		loss_df = fabs(df);
		if ((loss_x/(1.0+fabs(x))) < eps_fmin_tmp) {
			//	check for non-optimal solution
			if (loss_df > eps_fmin_tmp) {
				log_fout << endl << "Error fmin_R1, non-optimal solution" << endl;
				xmin = xnext;
				process_flag = X_SUBOPTIMAL_FLAG;
				break;
			}
			//	break for successful convergence
			else {
				xmin = xnext;
				process_flag = jfmin;
				break;
			}
		}
		//	increment iteration flag
		++jfmin;
		//	check if exceeded max # iterations
		if (jfmin >= Nmax_iter_tmp) {
			log_fout << endl << "Exceeded max #iter, fmin_R1" << endl;
			xmin = xnext;
			process_flag = EXCEEDED_MAXITER_FLAG;
			break;
		}
		//	write out iteration diagnostics
		if (Nfout!=0)
			if ((jfmin%Nfout)==0) {
				stoptime = clock();
				time_used = (stoptime - starttime)/1000.0;
				fout_fminR1_iter(log_fout, jfmin, Nmax_iter_tmp, x, xnext,
					loss_x, loss_fx, loss_df, df, d2f, time_used);
			}
		if ((Nfout == 0)&&(fout_fmin_iter_FLAG == ON)) {
			stoptime = clock();
			time_used = (stoptime - starttime)/1000.0;
			fout_fminR1_iter(log_fout, jfmin, Nmax_iter_tmp, x, xnext,
				loss_x, loss_fx, loss_df, df, d2f, time_used);
		}
		//	update candidate
		x = xnext;
	}
	//	write out function call summary
	if (fout_fmin_final_FLAG == ON) {
		stoptimef = clock();
		time_used = (stoptimef - starttimef)/1000.0;
		fout_fminR1_final(log_fout, process_flag, Nmax_iter_tmp, loss_x, loss_fx,
			loss_df, time_used, x0_tmp, xmin);
		log_fout << endl << "Exit fmin_B1, process_flag " << process_flag << "/"
			<< Nmax_iter_tmp << ", loss(x)=" << loss_x << ", loss(fx)=" << loss_fx
			<< ", loss(df)=" << loss_df << ", time_used=" << time_used;
		log_fout << endl << "x=" << x0_tmp << ", x*=" << xmin;
	}
	//	restore old step_Jacob
	stepJacob = step_jacob_old;
	return process_flag;
}
//	old fct call in fmin_v3.h and older
int fmin_R1(double x, double& xmin, double step_newton_tmp, double df, double d2f,
			double loss_fx,	int Nfout, double (*fct)(double x),
			int (*fct_constraints)(double x), double step_df,
			double eps_fmin_tmp, int Nmax_iter_tmp)
{
	ParamFmin param_tmp;
	//	put computation parameters into struct
	param_tmp.Nfout_fmin = Nfout;
	param_tmp.Nmax_iter_LM = Nmax_iter_tmp;
	param_tmp.step_newton_LM = step_newton_tmp;
	param_tmp.step_jacob = step_df;
	param_tmp.eps_LM = eps_fmin_tmp;

	int process_flag = fmin_R1(x, xmin, df, d2f, loss_fx, fct, fct_constraints, param_tmp);

	return process_flag;
}
/* version: Sept 14, 2005. Solves constrained minimization problem. Finds the min of a function (*fct) from R(N) to R using the BFGS method which uses gradient info. The gradient is computed numerically. The subrountine for finding a minimum along a direction is set by the linemin_flag. Latest change is the bundling of computation parameters into struct ParamFmin. Note that current version makes obsolete older versions of the function call as default values are now built in. */
int fmin_BFGS(double *x, double *xmin, const int Nx, double *Hx,
			  double (*fct)(double *x), int (*fct_constraints)(double *x),
			  ParamFmin& prm)
{
	//	hotfix - set xmin = x (prevent crashes from restarting fmin in calib_mom2mom)
	copy_array(xmin, x, Nx);
	//	extract computation parameters into locals
	double step_jacob_old = stepJacob;
	stepJacob = prm.step_jacob;
	double eps_fmin = prm.eps_fmin;
	double step_newton_LM = prm.step_newton_LM;
	double eps_LM = prm.eps_LM;
	int Nmax_iter_LM = prm.Nmax_iter_LM;
	int Nmax_iter_fmin = prm.Nmax_iter_fmin;
	int linemin_flag = prm.linemin_flag;
	int flag_fout_final = prm.flag_fout_final;
	int Nfout_fmin = prm.Nfout_fmin;
	int Nfout_LM = prm.Nfout_LM;
	//	rescale tolerance for dimension Nx
	eps_fmin /= pow(1.0*Nx, 0.5);

	//	local flags, counters and loss vars
	int process_flag = 0, process_flag_linemin;
	double norm_Jxnext, loss_x, loss = 1;
	//	init iteration counter
	int jfmin = 0;
	//	declare computation parameters linemin_newton
	double dfdx, d2fdx;
	double fxnext, fx, fx0;
	//	computation arrays, hessian, gradient and aux vars used in main loop
	double *Hxinv; Hxinv = new double [Nx*Nx];	// hessian current guess
	double *Hxnext; Hxnext = new double [Nx*Nx]; //	hessian next guess
	double *Jx;	Jx = new double [Nx]; // gradient current guess
	double *Jxnext; Jxnext = new double [Nx]; // gradient next guess
	double *xnext; xnext = new double [Nx];	 //	next guess
	double *z; z = new double [Nx];	//	aux var - used in updating hessian
	double *y; y = new double [Nx];	//	aux var - used in updating hessian
	double *s; s = new double [Nx];	//	current search direction

	//	make copy of initial guess for diagnostics
	double *x0_tmp;		x0_tmp = new double [Nx];
	copy_array(x0_tmp, x, Nx);
	//	set initial guess for H(k) to identity mat
	make_identity_matrix(Hx, Nx);
	//	compute jacobian at current init candidate
	get_gradient(Jx, x, fx0, Nx, fct, stepJacob);
	//	clocking and timeuse vars
	double time_used, starttime, stoptime, stoptimef, starttimef = clock();
	fx = fx0;

	//	main loop
	while (loss>eps_fmin) {

		//	start iteration clock
		starttime = clock();
		// step 1: Solve for search direction s(k): where H(k)s(k) = - J(x(k)) 

		//	compute inverse Hessian at current candidate
		inv_mat(Hxinv, Hx, Nx);
		//	compute search direction and rescale into unit vector
		mat_mult(Hxinv, Jx, s, Nx, Nx, 1);
		mat_rescale(s, s, -1.0, Nx);

		//	step 2: Find min along search direction
		if (linemin_flag == LINEMIN_GOLDEN_FLAG)
			fxnext = linemin_golden(x, s, xnext, Nx, fct);
		if (linemin_flag == LINEMIN_NEWTON_FLAG) {
			process_flag_linemin = linemin_newton(x, s, xnext, Nx, fxnext, dfdx, d2fdx,
				fct, fct_constraints, Nfout_LM, Nmax_iter_LM, step_newton_LM,
				stepJacob, eps_LM, OFF);
		}
		//	test for non-local convex error
		if (process_flag_linemin == NOT_LOCALLY_CONVEX_FLAG) {
			log_fout << endl << "Exiting fmin. Candidate not locally convex!";
			process_flag = process_flag_linemin;
			break;
		}

		//	check stopping condition
		loss_x = norm_sup(x, xnext, Nx);
		norm_Jxnext = norm_L2(Jx, Nx);	//	changed from Jxnext to Jx - sept26/05
		if (loss_x < 0.1*(eps_fmin + eps_fmin*norm_sup(x, Nx))) {
			//	check for x stuck at suboptimal point
			if (norm_Jxnext > eps_fmin*(1.0+fabs(fxnext))) {
				fout_fminBFGS_badx_error(log_fout, x, xnext, Jx, Nx, jfmin,
					Nmax_iter_fmin, loss_x, norm_Jxnext);
				process_flag = X_STUCK_FLAG;
				break;
			}
			//	break for successful convergence to xmin
			break;
		}
		//	stopping condition : break if |Jx|<eps
		if (norm_Jxnext < eps_fmin*(1.0+fabs(fxnext))) break;
		//	check for exceed max # iterations
		if (jfmin >= Nmax_iter_fmin) {
			fout_max_iter_error(log_fout, Nmax_iter_fmin, xnext, x0_tmp, Nx,
				"fmin_BFGS");
			process_flag = EXCEEDED_MAXITER_FLAG;
			break;
		}

		//	update gradient and hessian (either using BFGS=default or DFP)
		update_gradient_fminBFGS(y, z, x, xnext, Jx, Jxnext, fxnext, fct, Nx);
		if (fminBFGS_FLAG == ON) {	//	BFGS updating rule = default
			update_hessian_fminBFGS(Hx, Hxnext, x, y, z, Nx);
		}
		else {	//	DFP updating rule
			update_hessian_fminDFP(Hx, Hxnext, x, y, z, Nx);
		}

		//	write out general iteration diagnostics
		if (fout_fmin_diag_FLAG == ON) {
			fout_fminBFGS_diag(log_fout, x, xnext, y, z, Jxnext, Hxnext, Nx, jfmin,
				Nmax_iter_fmin, loss_x, norm_Jxnext);
		}

		//	write out iteration diagnostics
		if ((Nfout_fmin != 0) && (fout_fmin_iter_FLAG == OFF))
			if ((jfmin%Nfout_fmin) == 0) {
				stoptime = clock();
				time_used = (stoptime - starttime) / 1000;
				fout_fminBFGS_iter(log_fout, jfmin, Nmax_iter_fmin, time_used, fx,
					fxnext, x, xnext, Nx, norm_Jxnext, loss_x, prm.eps_fmin);
			}
		if (fout_fmin_iter_FLAG == ON) {
			stoptime = clock();
			time_used = (stoptime - starttime) / 1000;
			fout_fminBFGS_iter(log_fout, jfmin, Nmax_iter_fmin, time_used, fx,
				fxnext, x, xnext, Nx, norm_Jxnext, loss_x, prm.eps_fmin);
		}
		//	increment fmin iteration counter
		++jfmin;
		//	update Jacobian, Hessian and candidate solution
		copy_array(Hx, Hxnext, Nx*Nx);
		copy_array(x, xnext, Nx);
		copy_array(Jx, Jxnext, Nx);
		fx = fxnext;

	}

	//	save xmin and write out function call diagnostics
	copy_array(xmin, xnext, Nx);
	prm.loss = fx;
	if (process_flag >= 0) process_flag = jfmin;
	if ((flag_fout_final == ON) || (fout_fmin_final_FLAG == ON)) {
		stoptimef = clock();
		time_used = (stoptimef - starttimef)/ 1000.0;
		fout_fminBFGS_final(log_fout, x0_tmp, fx0, xmin, fxnext, Nx, process_flag,
			Nmax_iter_fmin, time_used, eps_fmin);
	}

	//	free up memory
	delete[] x0_tmp; delete[] Hxinv; delete[] Hxnext; delete[] z;
	delete[] Jx; delete[] Jxnext; delete[] xnext;
	delete[] y; delete[] s;
	//	restore stepJacob
	stepJacob = step_jacob_old;
	return process_flag;
}
//	old fct call from fmin_v3.h
int fmin_BFGS(double *x, double *xmin, const int Nx, double step_BFGS,
			  double *Hx, int Nfout, double (*fct)(double *x),
			  int (*fct_constraints)(double *x),
			  double step_jacob, double eps_tmp, double step_newton_LM,
			  double eps_LM, int Nmax_iter_LM, int Nmax_iter_tmp, int linemin_flag,
			  int flag_fout_final)
{
	//	put computation parameters into struct
	ParamFmin param_tmp;
	param_tmp.step_jacob = step_jacob;
	param_tmp.eps_fmin = eps_tmp;
	param_tmp.step_newton_LM = step_newton_LM;
	param_tmp.eps_LM = eps_LM;
	param_tmp.Nmax_iter_LM = Nmax_iter_LM;
	param_tmp.Nmax_iter_fmin = Nmax_iter_tmp;
	param_tmp.linemin_flag = linemin_flag;
	param_tmp.flag_fout_final = flag_fout_final;
	param_tmp.Nfout_fmin = Nfout;
	param_tmp.Nfout_LM = 0;
	//	call template function
	int process_flag = fmin_BFGS(x, xmin, Nx, Hx, fct, fct_constraints, param_tmp);
	return process_flag;
}
double fmin_BFGS(double *x0 /*initial_candidate*/, int Nx,
				 double *xmin /*solution*/, double (*fct)(double *x),
				 double step_jacob_tmp, double eps_fmin_tmp, int Nmax_iter_tmp,
				 int linemin_flag)
{
	//	make copy of initial candidate
	double *x;		x = new double [Nx];
	copy_array(x, x0, Nx);
	//	declare Hessian
	double *Hx;		Hx = new double [Nx*Nx];	//	Hessian

	//	set computation parameters for fmin_BFGS
	double step_BFGS_tmp = 1;
	int Nfout_tmp = 0;
	int Nmax_iter_LM = NmaxiterLinemin;

	int process_flag = fmin_BFGS(x, xmin, Nx, step_BFGS_tmp, Hx, Nfout_tmp, fct,
		no_constraints, step_jacob_tmp, eps_fmin_tmp, Nmax_iter_LM, Nmax_iter_tmp,
		linemin_flag, OFF);

	double loss = fct(xmin);

	delete[] x; delete[] Hx;
	return loss;
}
double fmin_BFGS(double *x0 /*initial_candidate*/, int Nx,
				 double *xmin /*solution*/, double (*fct)(double *x), int linemin_flag)
{
	//	make copy of initial guess
	double *x;		x = new double [Nx];
	copy_array(x, x0, Nx);
	//	declare Hessian
	double *Hx;		Hx = new double [Nx*Nx];
	//	set computation parameters for fmin_BFGS
	double eps_fmin_tmp = epsLinemin;
	int Nfout_tmp = 0;
	int Nmaxiter_tmp = NmaxiterFmin;
	int Nmax_iter_LM_tmp = NmaxiterLinemin;
	double step_jacob_tmp = stepJacob;
	double step_BFGS_tmp = 1.0;
	//	make fct call to template
	int process_flag = fmin_BFGS(x, xmin, Nx, step_BFGS_tmp, Hx, Nfout_tmp, fct, no_constraints,
		step_jacob_tmp, eps_fmin_tmp, Nmax_iter_LM_tmp, Nmaxiter_tmp, linemin_flag, OFF);

	double loss = fct(xmin);

	delete[] x; delete[] Hx;
	return loss;
}

//	fmin_DFP - is the same as the BFGS method except for the updating of the Hessian so this implementation switches a flag in fminBFGS where the default is to use the BFGS updating rule. The flag is switched back to the default upon exit of this fct call
int fmin_DFP(double *x, double *xmin, const int Nx, double *Hx,
			double (*fct)(double *x), int (*fct_constraints)(double *x), ParamFmin prm)
{
	fminBFGS_FLAG = OFF;
	int flag = fmin_BFGS(x, xmin, Nx, Hx, fct, fct_constraints, prm);
	fminBFGS_FLAG = ON;
	return flag;
}
/*  Implements Powell's conjugate direction method. First initialize all direction u(i) = e(i) for i \in {1,Nx} where e(i) are the basis vectors in R(Nx). Then while the stopping condition is not satisfied: (step-i) Save starting point x(0), (step-ii) for i = 1,..,Nx, move x(i-1) to the minimum along direction u(i) and call this point x(i), (step-iii) for i = 1,..,(Nx-1), set u(i) to u(i+1) (step-iv) set u(Nx) to x(Nx)-x(0, and (step-v) move x(Nx) to the minimum along direction u(Nx) and call this point p(0) */
int fmin_Powell(double *x, int Nx,
				void (*fct)(double *fx, double *x),
				double (*fct_constraints)(double *x))
{
    int process_flag = 0;
	double step_change = 1;
	double loss_current, loss_best = -HELL;
	int iiter = 0;				//	init counter for num of iterations
	int FAILURE_FLAG = -1;		//	set FAILURE_FLAG

    //  initialization and step (i)
    double *u;  u = new double [Nx*Nx];		//	search directions
    make_basis_vectors(u,Nx);       //  init all directions u(i) = e(i) to basis vectors

    double *x0;		x0 = new double [Nx];
    double *xx;		xx = new double [Nx*Nx];
	double *xnext;	xnext = new double [Nx];
	double *xbest;	xbest = new double [Nx];
    copy_array(x0,x,Nx);    //  save starting point in row 0

	while (step_change>epsFmin) {

		//  step (ii) For i=1..Nx, move x(i-1) to the min along dir u(i) and call this point x(i)
		for (int i=0; i<Nx; ++i) {
			//	Need to write the function below and uncomment step
			//line_min(x,u+i*Nx,Nx,xnext,fct,fct_constraints);        //  get x(i)
			copy_array(xx+i*Nx,x,Nx);   //  save each iterate on x(i) i = 1,..,Nx
			copy_array(x,xnext,Nx);     //  set current x to x(i)
		}   //  end for(i)

		//  step (iii) - For i=1..(Nx-1), set u(i) to u(i+1)
		for (int j=0; j<Nx; ++j) copy_array(u+j*Nx,u+(j+1)*Nx,Nx);

		//  step (iv) - Set u(Nx) to x(Nx)-x(0)
		for (int ii=0; ii<Nx; ++ii) *(u+(Nx-1)*Nx+ii) = *(x+ii) - *(x0+ii);

		//	step (v) - Move x(Nx) to the minimum along direction u(Nx) and call this point p(0)
		//loss_current = line_min(x,u+(Nx-1),Nx,xnext);
		copy_array(x,xnext,Nx);

		//	save best guess and check stopping condition
		if (loss_current < loss_best) { //	need to be fixed !!!!!!!!!!!!
			loss_best = loss_current;
			copy_array(xbest,x,Nx);
		}
		step_change = norm_L2(x,xx,Nx);

		//	check for too many iterations
		if (iiter>NmaxiterFmin) {
			log_fout << endl << "FMIN: Exceeded max num of iterations!!! Iter=" << iiter
				<< " Exiting...";
			process_flag = FAILURE_FLAG;
			break;
		}	//	end if()

	}	//	end while(loss)

    delete[] u; delete[] x0; delete[] xx;

	if (process_flag!=FAILURE_FLAG) process_flag = iiter;

    return process_flag;
}

/*	gets next candidate x(k+1) in broyden_method, also (i) gets the inverse jacobian, Jinv(k), (ii) evaluates the eqns at x, f(k), and (iii) computes the search direction, dir(k). */
void update_gradient_fminBFGS(double *y, double *z, double *x, double *xnext,
							   double *Jx, double *Jxnext, double& fx,
							   double (*fct)(double *x), int Nx)
{
	mat_add(z, xnext, x, -1.0, Nx);				//	find z(k) = x(k+1) - x(k)
	get_gradient(Jxnext, xnext, fx, Nx, fct);	//	get J(x(k+1))
	mat_add(y, Jxnext, Jx, -1.0, Nx);			//	find y(k) = J(x(k+1)) - J(x(k))
	return;
}
/* Updating the approx Hessian using BFGS updating formula. Compute (i) H(k)z(k)z(k)H(k) and (ii) z(k)H(k)z(k) where H1(k) = (i) and H1a(k) = H(k)z(k), H1b(k) = z(k)H(k) and H2 = z(k)H(k)z(k). They are calculated as: H1b = transpose(H1a), H2 = H1b(k)z(k). Note that H2(k) is a scalar. Compute (iii) y(k)y(k) and (iv) y(k)z(k) where H3(k) = (iii) and H4(k) = (iv) note that H4(k) is a scalar */
void update_hessian_fminBFGS(double *Hx, double *Hxnext, double *x, double *y, double *z, int Nx)
{
	//	declare local vars
	double *H1;			H1 = new double [Nx*Nx];
	double *H1a;		H1a = new double [Nx];
	double *H1b;		H1b = new double [Nx];
	double *H3;			H3 = new double [Nx*Nx];
	double H2, H4;
	//	compute H1a(k) = H(k)z(k)
	mat_mult(Hx, z, H1a, Nx, Nx, 1);
	//	compute H1b(k) = tranpose(H1a(k)) = z(k)H(k)
	transpose_mat(H1b, H1a, Nx, 1);
	//	compute H1(k) = H1a(k)H1b(k) = H(k)z(k)z(k)H(k)
	mat_mult(H1a, H1b, H1, Nx, 1, Nx);
	//	compute H2(k) = z(k)H(k)z(k)
	H2 = mat_mult(z, H1a, Nx);
	//	compute H1(k) = H1(k)/H2(k) = H(k)z(k)z(k)H(k)/z(k)H(k)z(k)
	mat_rescale(H1, H1, 1/H2, Nx*Nx);
	//	compute H3(k) = y(k)y(k)
	mat_mult(y, y, H3, Nx, 1, Nx);
	//	compute H4(k) = y(k)z(k)
	H4 = mat_mult(y, z, Nx);
	//	compute H3(k) = H3(k)/H4(k) = y(k)y(k)/y(k)z(k)
	mat_rescale(H3, H3, 1/H4, Nx*Nx);
	//	update hessian, H(k+1) = H(k) + H3(k) - H1(k)
	mat_add(Hxnext, H3, H1, -1.0, Nx*Nx);
	mat_add(Hxnext, Hx, Hxnext, 1.0, Nx*Nx);
	//	write out diagnostics
	if (fout_update_hessian_fmin_FLAG == ON)
		fout_update_hessian_fmin(log_fout, H3, H1, Hx, Hxnext, Nx);
	//	free mem
	delete[] H1; delete[] H1a; delete[] H1b; delete[] H3;
	return;
}
//	using formulation from numerical recipes which yields the same updated hessian as the original one (based on Judd)
void update_hessian_fminBFGS_NR(double *Hx, double *Hxnext, double *x, double *y,
								 double *z, int Nx)
{
	int Nprec = 6;
	//	compute y*Ty where T=transpose operator
	double *M1;	M1 = new double [Nx*Nx];
	mat_mult(y, y, M1, Nx, Nx);
	log_fout << "y*y outer product";
	fout_array(log_fout, M1, Nx, Nx, Nprec);
	//	compute Tz*y
	double M1_denom = mat_mult(z, y, Nx);
	//	compute (y*Ty)/(Tx*y)
	mat_rescale(M1, M1, 1/M1_denom, Nx*Nx);
	//	compute z*H
	double *M2; M2 = new double [Nx];
	mat_mult(z, Hx, M2, 1, Nx, Nx);
	//	compute outer product of (z*H) X (z*H) where X = outer product operator
	double *M3; M3 = new double [Nx*Nx];
	mat_mult(M2, M2, M3, Nx, Nx);
	//	compute Ty*(Hy)
	double M3_denom = mat_mult(z, M2, Nx);
	mat_rescale(M3, M3, 1/M3_denom, Nx*Nx);
	//	updating terms
	double *dM;	dM = new double [Nx*Nx];
	mat_add(dM, M1, M3, -1.0, Nx*Nx);
	//	update Hessian
	mat_add(Hxnext, Hx, dM, 1.0, Nx*Nx);
	//	write out diagnostics
	if (fout_update_hessian_fmin_FLAG == ON)
		fout_update_hessian_fmin(log_fout, M1, M3, Hx, Hxnext, Nx);
	//	free mem
	delete[] M1; delete[] M2; delete[] M3; delete[] dM;
	return;
}
//	Updating the approx Hessian using the DFP updating formula.
void update_hessian_fminDFP(double *Hx, double *Hxnext, double *x, double *y, double *z, int Nx)
{
	//	compute H*z (note it is of dim (N,1) --> dim(H*z) = (N,1) = col vector
	double *M1;	M1 = new double [Nx];
	mat_mult(Hx, z, M1, Nx, Nx, 1);
	//	compute M1 = y - H*z and dim(M1) = (N,1)
	mat_add(M1, y, M1, 1, Nx);
	//	compute M2 = (y - H*z) Ty where T = transpose operator and dim(M2) = (N,N)
	double *M2; M2 = new double [Nx*Nx];
	mat_mult(M1, y, M2, Nx, 1, Nx);
	//	compute M3 = y*T(y - H*z) and dim(M3) = (N,N)
	double *M3; M3 = new double [Nx*Nx];
	mat_mult(y, M1, M3, Nx, 1, Nx);
	//	compute M2 + M3
	double *M4; M4 = new double [Nx*Nx];
	mat_add(M4, M2, M3, 1, Nx*Nx);
	//	compute M4 = (M3 + M4)/(Ty*z) and dim(M4) = (N,N)
	double scalar1 = mat_mult(y, z, Nx);
	mat_rescale(M4, M4, 1/scalar1, Nx*Nx);
	//	compute Tz*(y - H*z) which is of dim (1,1)
	double scalar2 = mat_mult(z, M1, Nx);
	//	compute y*Ty which is of dim (N,N)
	double *M5; M5 = new double [Nx*Nx];
	mat_mult(y, y, M5, Nx, 1, Nx);
	//	compute M5 = Tz*(y - H*z)*y*Ty which of dim (N,N)
	mat_rescale(M5, M5, scalar2, Nx);
	//	divide M5 by (Ty*z)*(Ty*z) which is a scalar
	mat_rescale(M5, M5, 1/pow(scalar1, 2.0), Nx);
	//	get update part M0 = M4 - M5
	double *M0; M0 = new double [Nx*Nx];
	mat_add(M0, M4, M5, -1.0, Nx*Nx);
	//	get updated Hessian
	mat_add(Hxnext, Hx, M0, 1.0, Nx*Nx);
	//	write out diagnostics
	if (fout_update_hessian_fmin_FLAG == ON)
		fout_update_hessian_fmin(log_fout, M4, M5, Hx, Hxnext, Nx);
	//	free mem
	delete[] M0; delete[] M1; delete[] M2; delete[] M3; delete[] M4; delete[] M5;
	return;
}
//	implements the updating scheme for DFP as per numerical recipes
void update_hessian_fminDFP_NR(double *Hx, double *Hxnext, double *x, double *y, double *z, int Nx)
{
	//	get outer product of z X z where X = outer product
	double *M1; M1 = new double [Nx*Nx];
	mat_mult(z, z, M1, Nx, Nx);
	//	get dot product of z y
	double M2_tmp;
	M2_tmp = mat_mult(z, y, Nx);
	mat_rescale(M1, M1, 1/M2_tmp, Nx*Nx);
	//	compute H*y
	double *M3;	M3 = new double [Nx];
	mat_mult(Hx, y, M3, Nx, Nx, 1);
	//	get outer product of (H*y) X (H*y)
	double *M4; M4 = new double [Nx*Nx];
	mat_mult(M3, M3, M4, Nx, Nx);
	//	compute y*H*y
	double M5_tmp;
	double *M5; M5 = new double [Nx];
	mat_mult(y, Hx, M5, 1, Nx, Nx);
	M5_tmp = mat_mult(M5, y, Nx);
	mat_rescale(M4, M4, 1/M5_tmp, Nx*Nx);
	//	putting the updating additive terms together
	double *dH; dH = new double [Nx*Nx];
	mat_add(dH, M1, M4, -1.0, Nx*Nx);
	//	update hessian
	mat_add(Hxnext, Hx, dH, 1.0, Nx*Nx);

	return;
}
void fout_update_hessian_fmin(ofstream& fout, double *Hpos, double *Hneg,
							  double *Hx, double *Hxnext, int Nx)
{
	int Nprec = 6;
	//	write header
	fout << endl << "update_hessian_fmin_diag: dim=" << Nx << endl;

	fout << "Hessian: H(k)";
	fout_array(fout, Hx, Nx, Nx, Nprec);

	fout << "subtractive update term";
	fout_array(fout, Hneg, Nx, Nx, Nprec);

	fout << "additive update term";
	fout_array(fout, Hpos, Nx, Nx, Nprec);

	fout << "Hessian: H(k+1)";
	fout_array(fout, Hxnext, Nx, Nx, Nprec);

	return;
}
//	used in linemin_newton. Contracts xnext back toward x until xnext is locally convex. Rescales the size of change (along dir) by one-half with each attempt upto a max of 5 attempts. Exits with a failure flag if unable to find such a candidate with Nmaxiter_linemin_xnext such attempts. Note that for Nmaxiter_linemin_xnext = 5, this means that the decrease in movement along dir is reduced by 97 percent. Note that in linemin_newton, this fct is called when xnext is locally concave (ie d2f < 0)
int contract_linemin_newton_step(double& fx, double *x, double *xnext, double *dir, int Nx_tmp,
								 double (*fct)(double *x), double step_jacob, double& df,
								 double& d2f, double& fnext, double& dfnext, double& d2fnext,
								 double step_newton_tmp, int DIFF_3PT_FLAG)
{
	//	control-flow vars
	int flag = NOT_LOCALLY_CONVEX_FLAG;
	//	make copy of original candidate
	double *xnext_old; xnext_old = new double [Nx_tmp];
	copy_array(xnext_old, xnext, Nx_tmp);
	//	compute total change
	double *dx; dx = new double [Nx_tmp];
	mat_add(dx, xnext, x, -1, Nx_tmp);
	//	local vars
	double rescale_factor = 0.5;

	for (int i=0; i<Nmaxiter_linemin_xnext; ++i) {

		//	compute new candidate xnext
		mat_add(xnext, x, dx, rescale_factor, Nx_tmp);
		//	compute derivative info at new candidate
		d2fnext = get_d2fdx2_along_s(fnext, xnext, dir, Nx_tmp, fct, step_jacob,
			dfnext, DIFF_3PT_FLAG);
		//	case: xnext is locally convex - copy info and break
		if (d2fnext > 0) {
			flag = 0;
			break;
		}
		rescale_factor /= 2.0;
	}
	//	write out diagnostics if needed
	if (fout_contract_linemin_newton_step_FLAG == ON) {
//		fout_contract_linemin_newton_step( );
	}
	delete[] xnext_old; delete[] dx;
	return flag;
}
