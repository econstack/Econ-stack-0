#include <math.h>
#include <fstream>
#include <string>
#include <conio.h>
#include <iostream>
#include <iomanip>
#include <time.h>

using namespace std;
//	own libraries
#include "math_parameters.h"
#include "io_grunts.h"
#include "io_diag.h"
#include "math_grunts.h"
#include "NlinSolvers.h"

//	globals
int fout_xnext_newton_dpercent_FLAG = OFF;
int fout_get_xnext_newton_powell_FLAG = ON;

//	initialization of parameter values for class ParamNewton
ParamNewton::ParamNewton()
{
	Nf = 0; // default; solver will exit for system of dim 0
	eps_newton = 0.0001;
	Nfout = 0; // default = no iteration output to log
	Nmaxiter = 1000; // max number of iterations
	step_jacob = stepJacob; // step size in computing numerical derivatives
	step_newton = 0.5; // step size in newton step from x to x'
	step_dpercent_flag = OFF; // option: constrain max percent change from x to x'
	step_dpercent = 0.1; // max percent change for option
	powell_flag = OFF;
	newton_flag = 0; // use newton_method by default
}
/*  TEMPLATE. Solves a system of equations of dim Nx by Newton's method. Function pointers to the system of  equations (*fct) and to the system of constraints (*fct_constraints) must be passed */
int newton_method(double *x, int Nx, void (*fct)(double *fx, double *x),
				  int (*fct_constraints)(double *x), ParamNewton& pm)
{
	//	store old stepJacob and put current newton parameter as global
	double stepJacob_old = stepJacob;
	stepJacob = pm.step_jacob;

	//	make copy of initial guess for diagnostics
	double *x0;	x0 = new double [Nx];
	copy_array(x0, x, Nx);

	//	time use vars
	double time_used, starttime, stoptime, stoptimef, starttimef = clock();
    //  local math vars
	double *step_x; step_x = new double [Nx];	// step from x(k) to x(k+1)
    double *xnext; xnext = new double [Nx];		// next guess
    double *fx;	fx = new double [Nx];			// error array at x(k)
	double *fxnext;	fxnext = new double [Nx];	// error array at x(k+1)
    double *Jx;	Jx = new double [Nx*Nx];		// Jacobian at x
    double *Jxinv; Jxinv = new double [Nx*Nx];	// Jacobian inverse at x

    //  control-flow vars
    double loss, loss_fx, loss_x = 0;
    int i_Nlin = 0; // Nlin solver loop counter
	int process_flag = 0; // flag signaling (non)convergence to be returned
	int inv_mat_flag; // flag for mat inv step

	//	prelim testing of candidate solution
	fct(fx, x);
	loss_fx = norm_sup(fx, Nx);

    //  main loop ---------------------------------------------------------------------------
	while ( (loss_fx > pm.eps_newton) && (i_Nlin < pm.Nmaxiter) )   {

		/*	newton iteration diagnostics - start iteration clock */
		starttime = clock();

		//  get Jacobian at x
		get_jacob(x, Jx, fx, Nx, fct, stepJacob);
		//  get inverse Jacobian at x, break if unable to invert
        inv_mat_flag = inv_mat(Jxinv, Jx, Nx);
		if (inv_mat_flag == SINGULAR_MATRIX_FLAG)  {
			fout_iter_exit(log_fout, x0, x, Nx, i_Nlin, pm.Nmaxiter,
				"newton_method: singular jacobian");
			process_flag = SINGULAR_MATRIX_FLAG;
			break;
		}

        //  compute next guess for x: x(k+1) = x(k) - Jinv(x)f(x)
		mat_mult(Jxinv, fx, step_x, Nx, Nx, 1);
		mat_rescale(step_x, pm.step_newton, Nx);
		mat_add(xnext, x, step_x, -1.0, Nx);

		//	option: check if x(k+1) is within percent change limit, contract step if needed
		if (pm.step_dpercent_flag == ON) {
			get_xnext_newton_dpercent(xnext, x, Nx, pm.step_dpercent);
		}

        //  check if next candidate, x(k+1) is in constraint set, rescale if needed
        int flag_constraint = fct_constraints(xnext);
        if (flag_constraint>0)  {
            int flag_xnext = rescale_xnext(x, Jxinv, fx, xnext, fct_constraints, Nx);
            if (flag_xnext > 0) {
				fout_xnext_error(log_fout, x, xnext, Nx);
				process_flag = NONCONVERGENCE_FLAG;
                break;
            }
        }
		/*  step 3: check stopping condition. Compute error at current candidate x */
		loss_x = norm_sup(xnext, x, Nx);
        loss_x = loss_x / (1.0 + norm_sup(x, Nx));
        fct(fxnext, xnext);
        loss_fx = norm_sup(fxnext, Nx);
		loss = loss_x + loss_fx;

		//	option: Powell hybrid, check if the sum of residuals decreases from x(k) to x(k+1)
		if (pm.powell_flag == ON) {
			int powell_flag = get_xnext_newton_powell(xnext, x, step_x, fx, fxnext, Nx, fct);
			//	break if unable to find candidate st sum of residuals decrease
			if (powell_flag < 0) {
				process_flag = NEWTON_POWELL_FLAG;
				break;
			}
		}
		//  option: newton iteration diagnostics
		if (pm.Nfout != 0)
			if ((i_Nlin%pm.Nfout)==0) {
				stoptime = clock();
				time_used = (stoptime-starttime)/1000.0;
				fout_nlin_iter(log_fout, Nx, i_Nlin, loss_x, loss_fx, time_used,
					x, xnext, fx, fxnext, "newton_method");
				//log_fout << setw(12) << "Jx";
				//fout_array(log_fout, Jx, Nx, Nx, 6, 1);
			}
		if (fout_nlin_iter_FLAG == ON) {
			stoptime = clock();
			time_used = (stoptime-starttime)/1000.0;
			fout_nlin_iter(log_fout, Nx, i_Nlin, loss_x, loss_fx, time_used,
				x, xnext, fx, fxnext, "newton_method");
		}

		//	increment iteration counter
        ++i_Nlin;
		//	break if exceeded max # iterations
        if (i_Nlin >= pm.Nmaxiter)    {
			fout_max_iter_error(log_fout, pm.Nmaxiter, xnext, x0, Nx, "newton_method");
			process_flag = EXCEEDED_MAXITER_FLAG;
			break;
        }

		//	break if x does not change
		if ((norm_sup(xnext, x, Nx) < (pm.eps_newton/100))&&(loss_fx > pm.eps_newton)) {
			fout_iter_exit(log_fout, x0, x, Nx, i_Nlin, pm.Nmaxiter, "x_stuck");
			log_fout << " x'=";
			fout_array(log_fout, xnext, 1, Nx, 6, 0);
			process_flag = X_STUCK_FLAG;
			break;
		}


        //  update candidate x and fx
        copy_array(x, xnext, Nx);
		copy_array(fx, fxnext, Nx);
	}   //  end while()

	//  compute and write out summary stats for function call
	if (process_flag >= 0) process_flag = i_Nlin;
	loss = loss_fx;
	stoptimef = clock();
	time_used = (stoptimef - starttimef)/1000.0;
	if ((fout_nlin_final_FLAG == ON) && (process_flag >= fout_nlin_final_NMIN)) {
		fout_nlin_final(log_fout, Nx, process_flag, loss_x, loss_fx, x0, x,
			fx, time_used);
		cout << endl << "Exit: newton_method, size=" << Nx << ", iter=" << process_flag
			<< "/" << pm.Nmaxiter << ", loss(x)=" << loss_x << ", loss(fx)=" << loss_fx;
	}
	//	always write out failures to converge
	if ((process_flag < 0)&&(fout_nlin_final_FLAG == OFF)) {
		fout_nlin_final(log_fout, Nx, process_flag, loss_x, loss_fx, x0, x,
			fx, time_used);
		cout << endl << "Exit: newton_method, size=" << Nx << ", iter=" << process_flag
			<< "/" << pm.Nmaxiter << ", loss(x)=" << loss_x << ", loss(fx)=" << loss_fx;
	}
	//	return loss and time_used
	pm.loss = loss;
	pm.time_used = time_used;

	//	free up memory
    delete[] xnext; delete[] fx; delete[] Jx; delete[] Jxinv; delete[] x0;
	delete[] fxnext; delete[] step_x;
    return process_flag;
}
int newton_method_old(double *x, int Nx, void (*fct)(double *fx, double *x),
				  int (*fct_constraints)(double *x), ParamNewton& pm_newt)
{
	//	extract computation parameters
	int Nmaxiter_tmp = pm_newt.Nmaxiter;
	double step_newton = pm_newt.step_newton;
	int Nfout_tmp = pm_newt.Nfout;
	double eps_newton = pm_newt.eps_newton;
	double step_jacob = pm_newt.step_jacob;
	int step_dpercent_flag = pm_newt.step_dpercent_flag;
	double step_newton_dp = pm_newt.step_dpercent;
	//	make copy of initial guess for diagnostics
	double *x0;		x0 = new double [Nx];
	copy_array(x0, x, Nx);
	//	time use vars
	double time_used, starttime, stoptime, stoptimef, starttimef = clock();
    //  local math vars
	double *step_x; step_x = new double [Nx]; // step from x(k) to x(k+1)
    double *xnext; xnext = new double [Nx];	//  next guess
    double *fx;	fx = new double [Nx]; //  error array at x(k)
	double *fxnext;	fxnext = new double [Nx]; //	error array at x(k+1)
    double *Jx;	Jx = new double [Nx*Nx]; //  Jacobian at x
    double *Jxinv; Jxinv = new double [Nx*Nx]; //  Jacobian inverse at x
    //  control-flow vars
    double loss, loss_fx, loss_x = 0;
    int i_Nlin = 0; //  Nlin solver loop counter
	int process_flag = 0; // flag signaling (non)convergence to be returned
	int inv_mat_flag;
	double dp_incr, dp_decr, dp_max;
	//	prelim testing of candidate solution
	fct(fx, x);
	loss_fx = norm_sup(fx, Nx);

    //  main loop ---------------------------------------------------------------------------
    while (((loss_x+loss_fx)>eps_newton) && (i_Nlin<Nmaxiter_tmp))   {

		/*	newton iteration diagnostics - start iteration clock */
		starttime = clock();
		//  get Jacobian at x
        get_jacob(x, Jx, fx, Nx, fct, step_jacob);
		/*  get inverse Jacobian at x, break if unable to invert */
        inv_mat_flag = inv_mat(Jxinv, Jx, Nx);
		if (inv_mat_flag == SINGULAR_MATRIX_FLAG)  {
			fout_iter_exit(log_fout, x0, x, Nx, i_Nlin, Nmaxiter_tmp, "newton_method");
			process_flag = SINGULAR_MATRIX_FLAG;
			break;
		}
        /*  compute next guess for x */
		mat_mult(Jxinv, fx, step_x, Nx, Nx, 1);
		mat_rescale(step_x, step_x, step_newton, Nx);
		mat_add(xnext, x, step_x, -1.0, Nx);
		//	option: check if x(k+1) is within percent change limit, contract step if needed
		if (step_dpercent_flag == ON) {
			get_xnext_newton_dpercent(xnext, x, Nx, step_newton_dp);
		}

        /*  check if guess is in constraint set, rescale if needed */
        int flag_constraint = fct_constraints(xnext);
        if (flag_constraint>0)  {
            int flag_xnext = rescale_xnext(x, Jxinv, fx, xnext, fct_constraints, Nx);
            if (flag_xnext > 0) {
				fout_xnext_error(log_fout, x, xnext, Nx);
				process_flag = NONCONVERGENCE_FLAG;
                break;
            }
        }
		/*  step 3: check stopping condition. Compute error at current candidate x */
		loss_x = norm_sup(xnext, x, Nx);
        loss_x = loss_x / (1.0 + norm_sup(x, Nx));
        fct(fxnext, xnext);
        loss_fx = norm_sup(fxnext, Nx);
		loss = loss_x + loss_fx;
		//  option: newton iteration diagnostics
		if (Nfout_tmp!=0)
			if ((i_Nlin%Nfout_tmp)==0) {
				stoptime = clock();
				time_used = (stoptime-starttime)/1000.0;
				fout_nlin_iter(log_fout, Nx, i_Nlin, loss_x, loss_fx, time_used,
					x, xnext, fx, fxnext, "newton_method");
				//log_fout << setw(12) << "Jx";
				//fout_array(log_fout, Jx, Nx, Nx, 6, 1);
			}
		if (fout_nlin_iter_FLAG == ON) {
			stoptime = clock();
			time_used = (stoptime-starttime)/1000.0;
			fout_nlin_iter(log_fout, Nx, i_Nlin, loss_x, loss_fx, time_used,
				x, xnext, fx, fxnext, "newton_method");
		}

		//	increment iteration counter
        ++i_Nlin;
		//	break if exceeded max # iterations
        if (i_Nlin>=Nmaxiter_tmp)    {
			fout_max_iter_error(log_fout, Nmaxiter_tmp, xnext, x0, Nx, "newton_method");
			process_flag = EXCEEDED_MAXITER_FLAG;
			break;
        }

		//	break if x does not change
		if (step_dpercent_flag != ON) {
			if (norm_sup(xnext, x, Nx)<(eps_newton/1000)) {
				fout_iter_exit(log_fout, x0, x, Nx, i_Nlin, Nmaxiter_tmp, "x_stuck");
				log_fout << " x'=";
				fout_array(log_fout, xnext, 1, Nx, 6, 0);
				process_flag = X_STUCK_FLAG;
				break;
			}
		}
		else {
			//	compute largest percent change
			mat_div(step_x, xnext, x, Nx);
			dp_incr = get_max(step_x, Nx);
			dp_decr = get_min(step_x, Nx);
			dp_max = max(dp_incr, 1/dp_decr);
			if (dp_max < eps_newton) {
				fout_iter_exit(log_fout, x0, x, Nx, i_Nlin, Nmaxiter_tmp, "x_stuck");
				process_flag = X_STUCK_FLAG;
				break;
			}
		}
        //  update candidate x and fx
        copy_array(x, xnext, Nx);
		copy_array(fx, fxnext, Nx);
    }   //  end while()

	//  compute and write out summary stats for function call
	if (process_flag >= 0) process_flag = i_Nlin;
	loss = loss_fx;
	stoptimef = clock();
	time_used = (stoptimef - starttimef)/1000.0;
	if (fout_nlin_final_FLAG == ON)  {
		fout_nlin_final(log_fout, Nx, process_flag, loss_x, loss_fx, x0, x,
			fx, time_used);
		cout << endl << "Exit: newton_method, size=" << Nx << ", iter=" << process_flag
			<< "/" << Nmaxiter_tmp << ", loss(x)=" << loss_x << ", loss(fx)=" << loss_fx;
	}
	//	always write out failures to converge
	if ((process_flag < 0)&&(fout_nlin_final_FLAG == OFF)) {
		fout_nlin_final(log_fout, Nx, process_flag, loss_x, loss_fx, x0, x,
			fx, time_used);
		cout << endl << "Exit: newton_method, size=" << Nx << ", iter=" << process_flag
			<< "/" << Nmaxiter_tmp << ", loss(x)=" << loss_x << ", loss(fx)=" << loss_fx;
	}
	//	return loss and time_used
	pm_newt.loss = loss;
	pm_newt.time_used = time_used;

	//	free up memory
    delete[] xnext; delete[] fx; delete[] Jx; delete[] Jxinv; delete[] x0;
	delete[] fxnext; delete[] step_x;
    return process_flag;
}
//	old template (v3 and older used this as main fct call
int newton_method(double *x, int Nx, void (*fct)(double *fx, double *x),
				  int (*fct_constraints)(double *x), double& loss, int Nmaxiter_tmp,
				  double step_newton, int Nfout_tmp, double &time_used,
				  double eps_newton, double step_jacob)
{
	//	zero for supnorm
	double *zero_tmp;	zero_tmp = new double [Nx];
	const_array(zero_tmp, Nx, 0);
	//	make copy of initial guess for diagnostics
	double *x0;		x0 = new double [Nx];
	copy_array(x0, x, Nx);

	int process_flag = 0;			//	flag signaling (non)convergence to be returned
	int inv_mat_flag;

	double step_x;
	double starttime, stoptime, stoptimef, starttimef = clock();

    //  local vars
    double *xnext;		xnext = new double [Nx];	//  next guess
    double *fx;			fx = new double [Nx];       //  error array at x(k)
	double *fxnext;		fxnext = new double [Nx];	//	error array at x(k+1)
    double *Jx;			Jx = new double [Nx*Nx];    //  Jacobian at x
    double *Jxinv;		Jxinv = new double [Nx*Nx]; //  Jacobian inverse at x

    //  declare and init while-loop controls
    double loss_x = 0;
	fct(fx, x);
	double loss_fx = norm_sup(fx, Nx);		//  init loss
    int i_Nlin = 0;                         //  Nlin solver loop counter

    //  main loop ---------------------------------------------------------------------------
    while (((loss_x+loss_fx)>eps_newton) && (i_Nlin<Nmaxiter_tmp))   {

		/*	newton iteration diagnostics - start iteration clock */
		starttime = clock();

		//  get Jacobian at x
        get_jacob(x, Jx, fx, Nx, fct, step_jacob);

		/*  get inverse Jacobian at x, break if unable to invert */
        inv_mat_flag = inv_mat(Jxinv, Jx, Nx);
		if (inv_mat_flag == SINGULAR_MATRIX_FLAG)  {
			fout_iter_exit(log_fout, x0, x, Nx, i_Nlin, Nmaxiter_tmp, "newton_method");
			process_flag = SINGULAR_MATRIX_FLAG;
			break;
		}

        /*  compute next guess for x */
        for (int i=0; i<Nx; ++i) {
			step_x = mat_mult(Jxinv+i*Nx, fx, Nx);
            step_x = step_newton*step_x;
            *(xnext+i) = *(x+i)-step_x;
        }   //  end for(i)

        /*  check if guess is in constraint set, rescale if needed */
        int flag_constraint = fct_constraints(xnext);
        if (flag_constraint>0)  {
            int flag_xnext = rescale_xnext(x,Jxinv,fx,xnext,fct_constraints,Nx);
            if (flag_xnext<0) {
				fout_xnext_error(log_fout, x, xnext, Nx);
				process_flag = NONCONVERGENCE_FLAG;
                break;
            }
        }
		//	write out iteration diagnostics
		if (fout_nlin_iter_FLAG == ON) {
			log_fout << endl << " x="; fout_array(log_fout, x, 1, Nx, 6, 0);
			log_fout << "x'="; fout_array(log_fout, xnext, 1, Nx, 6, 0);
		}
        /*  step 3: check stopping condition. Compute error at current candidate x */
		loss_x = norm_sup(xnext, x, Nx);
        loss_x = loss_x/(1.0 + norm_sup(x, zero_tmp, Nx));

        fct(fxnext, xnext);
        loss_fx = norm_sup(fxnext, zero_tmp, Nx);
		loss = loss_x + loss_fx;

		//  newton iteration diagnostics
		if (Nfout_tmp!=0)
			if ((i_Nlin%Nfout_tmp)==0) {
				stoptime = clock();
				time_used = (stoptime-starttime)/1000.0;
				fout_nlin_iter(log_fout, Nx, i_Nlin, loss_x, loss_fx, time_used,
					x, xnext, fx, fxnext, "newton_method");
				log_fout << setw(12) << "Jx";
				fout_array(log_fout, Jx, Nx, Nx, 6, 1);
			}
		if (fout_nlin_iter_FLAG == ON) {
			stoptime = clock();
			time_used = (stoptime-starttime)/1000.0;
			fout_nlin_iter(log_fout, Nx, i_Nlin, loss_x, loss_fx, time_used,
				x, xnext, fx, fxnext, "newton_method");
		}

		//	increment iteration counter
        ++i_Nlin;
		//	break if exceeded max # iterations
        if (i_Nlin>=Nmaxiter_tmp)    {
			fout_max_iter_error(log_fout, Nmaxiter_tmp, xnext, x0, Nx, "newton_method");
			process_flag = EXCEEDED_MAXITER_FLAG;
			break;
        }

		//	break if x does not change
		if (norm_sup(xnext, x, Nx)<(eps_newton/100)) {
			fout_iter_exit(log_fout, x0, x, Nx, i_Nlin, Nmaxiter_tmp, "x_stuck");
			process_flag = X_STUCK_FLAG;
			break;
		}
        //  update candidate x and fx
        copy_array(x, xnext, Nx);
		copy_array(fx, fxnext, Nx);

    }   //  end while()

	//  compute and write out summary stats for function call
	if (process_flag>=0) process_flag = i_Nlin;
	loss = loss_fx;
	stoptimef = clock();
	time_used = (stoptimef - starttimef)/1000.0;
	if (fout_nlin_final_FLAG == ON)  {
		fout_nlin_final(log_fout, Nx, process_flag, loss_x, loss_fx, x0, x,
			fx, time_used);
		cout << endl << "Exit: newton_method, size=" << Nx << ", iter=" << process_flag
			<< "/" << Nmaxiter_tmp << ", loss(x)=" << loss_x << ", loss(fx)=" << loss_fx;
	}
	//	always write out failures to converge
	if ((process_flag < 0)&&(fout_nlin_final_FLAG == OFF)) {
		fout_nlin_final(log_fout, Nx, process_flag, loss_x, loss_fx, x0, x,
			fx, time_used);
		cout << endl << "Exit: newton_method, size=" << Nx << ", iter=" << process_flag
			<< "/" << Nmaxiter_tmp << ", loss(x)=" << loss_x << ", loss(fx)=" << loss_fx;
	}
	//	free up memory
    delete[] xnext; delete[] fx; delete[] Jx; delete[] Jxinv; delete[] x0;
	delete[] fxnext; delete[] zero_tmp;
    return process_flag;
}

/*No constraints version*/
int newton_method(double *x, int Nx, void (*fct)(double *fx, double *x),
				  double& loss, int Nmaxiter_tmp, double step_newton,
				  int Nfout_tmp, double &time_used, double eps_newton)
{
	int process_flag = newton_method(x, Nx, fct, no_constraints, loss, Nmaxiter_tmp,
		step_newton, Nfout_tmp, time_used, eps_newton);

	return process_flag;
}
int newton_method(double *x, int Nx, void (*fct)(double *fx, double *x),
				  int (*fct_constraints)(double *x) /*default version*/)
{
	double loss, time_used;
	int Nmax_iter_newton_method = MaxIter_Nlin;
	double step_newton_method = step_NlinNewton;
	int Nfout_newton_method = 1;

	int flag = newton_method(x, Nx, fct, fct_constraints, loss, Nmax_iter_newton_method,
		step_newton_method, Nfout_newton_method, time_used);

	return flag;
}
//	----------------- newton_method() utilities --------------------------------------
int get_xnext_newton_powell(double *xnext, double *x, double *dx, double *fx,
							double *fxnext, int Nx, void (*fct)(double *fx, double *x))
{
	int Nmaxiter_tmp = 10;
	int jiter = 1;
	double rescale_factor = 0.5;

	int process_flag = 0;

	double loss_fx = norm_L2(fx, Nx);
	double loss_fxnext = norm_L2(fxnext, Nx);
	double loss_fxnext_tmp = loss_fxnext;

	double *fxnext_old; fxnext_old = new double [Nx];
	copy_array(fxnext_old, fxnext, Nx);
	double *xnext_old; xnext_old = new double [Nx];
	copy_array(xnext_old, xnext, Nx);

	//	 rescale if residuals increase by factor half with each iter upto Nmaxiter=10
	double *fxnext_tmp; fxnext_tmp = new double [Nx];
	double *xnext_tmp; xnext_tmp = new double [Nx];

	while (loss_fxnext_tmp > loss_fx) {

		mat_add(xnext_tmp, x, dx, -rescale_factor, Nx);
		fct(fxnext_tmp, xnext_tmp);
		loss_fxnext_tmp = norm_L2(fxnext_tmp, Nx);
		rescale_factor *= rescale_factor;
		++jiter;

		//	break if residual decreases
		if (loss_fxnext_tmp < loss_fx) {
            copy_array(xnext, xnext_tmp, Nx);
			copy_array(fxnext, fxnext_tmp, Nx);
			break;
		}
		//	break if exceeded max number of attempts
		if (jiter == Nmaxiter_tmp) {
			process_flag = NEWTON_POWELL_FLAG;
			break;
		}

	}

	//	option: write out diagnostics
	if (fout_get_xnext_newton_powell_FLAG == ON) {
		fout_get_xnext_newton_powell(log_fout, x, xnext_old, xnext_tmp, Nx, loss_fx, loss_fxnext,
			loss_fxnext_tmp);
	}

	delete[] fxnext_tmp; delete[] xnext_tmp; delete[] fxnext_old; delete[] xnext_old;

	return process_flag;
}
/* TEMPLATE. Solves a set of nonlinear eqns in R(Nx) using broyden method. Returns a flag for the number of iterations taken (=1 if failure); returns info on the error term st loss = sup_norm(x) and Nmaxiter_tmp = max number of iterations, Nfout_tmp = no iterations between fout to log file and to console, time_used = total time used*/
int broyden_method(double *x /*solution*/, int Nx,	void (*fct)(double *fx, double *x),
				   int (*fct_constraints)(double *x), double& loss /*residual at x*/,
				   int Nmaxiter_tmp, int Nfout_tmp, double& time_used,
				   double step_broyden_tmp, double eps_broyden_tmp)
{
	double starttime, stoptime, stoptimef, starttimef = clock();
    int flag_xnext_tmp, process_flag = 0;

    //  Initialization: Set A(0) = I  (identity matrix)
    double *Jx;		Jx = new double [Nx*Nx];
	make_identity_matrix(Jx, Nx);
	//	make copy of init guess
	double *x0_tmp;	x0_tmp = new double [Nx];
	copy_array(x0_tmp, x, Nx);

    double loss_x = 1; double loss_fx = 1;
    loss = loss_x + loss_fx;

    //  Local vars
    double *fx;         fx = new double [Nx];
    double *fxnext;		fxnext = new double [Nx];
    double *Jxinv;		Jxinv = new double [Nx*Nx];
    double *Jxnext;		Jxnext = new double [Nx*Nx];
    double *xnext;      xnext = new double [Nx];
    double *s_step;     s_step = new double [Nx];
    double *y_temp;     y_temp = new double [Nx];

	int i_Nlin = 0;

    //  Main loop ******************
    while (loss>eps_broyden_tmp)   {

		//	start iteration clock
		starttime = clock();

		/*  step 1. Compute next iterate: (i) Solve J(k)s(k) = -f(x(k)) for s(k), (ii) set x(k+1) = x(k) + s(k). In (i), solve for s(k) = - J(k)^{-1}f(x(k)) where J(k) = Jx, s(k) = s_step, x(k) = x, x(k+1) = xnext, f(x(k)) = fx */
		flag_xnext_tmp = get_xnext_broyden(Jx, Jxinv, fx, s_step, x, xnext, fct, Nx,
			step_broyden_tmp);
		//	break if unable to invert jacobian
		if (flag_xnext_tmp < 0) {
			fout_iter_exit(log_fout, x0_tmp, x, Nx, i_Nlin, Nmaxiter_tmp,
				"broyden, singular jacobian");
			return flag_xnext_tmp;
		}

		/*  step 2. Check for candidate x(k+1) in constraint set, rescale if needed */
        int flag_constraint = fct_constraints(xnext);
        if (flag_constraint>0)  {
            int flag_xnext = rescale_xnext(x, Jxinv, fx, xnext, fct_constraints, Nx);
            if (flag_xnext < 0) {
				fout_xnext_error(log_fout, x, xnext, Nx);
				process_flag = CONSTRAINT_SET_VIOLATION_FLAG;
                break;
            }   //  end if(flag_xnext)
        }   //  end if(flag_constraint)


		/*  step 3. Update Jacobian guess: Set y(k) = f(x(k+1)) - f(x(k)) and set J(k+1) = J(k) + {(y(k)-J(k)s(k))s(k)'}/{s(k)'s(k)} where J(k+1) = Jxnext, f(x(k+1)) = fxnext, y(k) = y_temp */
		get_jxnext_broyden(fx, fxnext, xnext, y_temp, s_step, Jx, Jxnext, Nx, fct);
		/*  step 4. Check stopping conditions, write out diagnostics if needed, update candidate x */
        loss_x = norm_sup(x, xnext, Nx);
        loss_x = loss_x/(1.0 + norm_sup(x, Nx));
        loss_fx = norm_sup(fxnext, Nx);
        loss = loss_x + loss_fx;
		//	break if successful convergence
		if (loss<eps_broyden_tmp) {
			process_flag = i_Nlin;
			break;
		}
		//	write out general diagnostics - EMPTY
		if (fout_nlin_iter_FLAG == ON) {

		}
		//	write out iteration diagnostics
		if (Nfout_tmp!=0)
			if ((i_Nlin%Nfout_tmp)==0) {
				stoptime = clock();
				time_used = (stoptime-starttime)/1000;
				fout_nlin_iter(log_fout, Nx, i_Nlin, loss_x, loss_fx, time_used,
					x, xnext, fx, fxnext, "broyden_method");
        }
		if (fout_nlin_iter_FLAG==ON) {
			stoptime = clock();
			time_used = (stoptime-starttime)/1000;
			fout_nlin_iter(log_fout, Nx, i_Nlin, loss_x, loss_fx, time_used,
				x, xnext, fx, fxnext, "broyden_method");
		}
		//	increment iteration counter
        ++i_Nlin;
		//	break if exceeded max # iterations
        if (i_Nlin>=Nmaxiter_tmp)    {
			fout_max_iter_error(log_fout, Nmaxiter_tmp, xnext, x0_tmp, Nx, "broyden_method");
			return EXCEEDED_MAXITER_FLAG;
        }

        //  update candidate x and jacobian
        copy_array(x, xnext, Nx);
        copy_array(Jx, Jxnext, Nx*Nx);

    }   //  end while(loss)

    //  write out summary stats on computation ----------------------------------------------
	stoptimef = clock();
	time_used = (stoptimef - starttimef)/1000.0;		//	save cpu time used

	if (fout_nlin_final_FLAG == ON) {
		fout_iter_final_diag(log_fout, x0_tmp, xnext, Nx, process_flag, Nmaxiter_tmp,
			time_used, "broyden_method");
		cout << endl << "Exit: broyden_method, size=" << Nx << ", iter=" << process_flag
			<< "/" << Nmaxiter_tmp << ", loss(x)=" << loss_x << ", loss(fx)=" << loss_fx;
	}

    //  delete local pointers
    delete[] fx; delete[] fxnext; delete[] Jx; delete[] Jxinv;
	delete[] Jxnext; delete[] xnext; delete[] s_step;
	delete[] y_temp; delete[] x0_tmp;

    return process_flag;
}


/* No constrainst version*/
int broyden_method(double *x, int Nx, void (*fct)(double *fx, double *x),
				   double &loss, int Nmaxiter_tmp, int Nfout_tmp, double& time_used)
{
	int process_flag = broyden_method(x, Nx, fct, no_constraints, loss, Nmaxiter_tmp,
		Nfout_tmp, time_used);
	return process_flag;
}
int broyden_method(double *x, int Nx, void (*fct)(double *fx, double *x),
				   int (*fct_constraints)(double *x) /*default version*/)
{
	int Nfout_tmp = 1;
	double loss, time_used;

	int process_flag = broyden_method(x, Nx, fct, fct_constraints, loss, MaxIter_Nlin, Nfout_tmp,
		time_used);
	return process_flag;
}
//	solves for the eqm of a system as a fixed point of a nlin set of eqns
int get_eqm(double *x0, double *x, double *fx, void (*fct)(double *fx, double *x),
			int Nx_tmp, int Nx0_tmp, int (*fct_constraints)(double *x),
			double& loss, double& time_used, int Nmaxiter_tmp, int Nfout,
			double step_newton, double eps_newton, double step_jacob)
{
	int process_flag = 0;
	double starttime, stoptime, time_used_newton;
	starttime = clock();

	for (int ix0=0; ix0<Nx0_tmp; ++ix0) {

		copy_array(x, x0+ix0*Nx_tmp, Nx_tmp);
		process_flag = newton_method(x, Nx_tmp, fct, fct_constraints, loss, Nmaxiter_tmp,
			step_newton, Nfout, time_used_newton, eps_newton, step_jacob);

		if (process_flag == 0) break;
		else fout_get_eqm_diag(log_fout, x0+ix0*Nx_tmp, x, Nx_tmp, fct);

	}

	fct(fx, x);		//	evaluate equations at solution and return residuals
	loss = norm_sup(fx, Nx_tmp);

	stoptime = clock();
	time_used = (stoptime - starttime)/1000.0;

	return process_flag;
}
/*	gets next candidate x(k+1) in broyden_method, also (i) gets the inverse jacobian, Jinv(k), (ii) evaluates the eqns at x, f(k), and (iii) computes the search direction, dir(k). */
int get_xnext_broyden(double *Jx, double *Jxinv, double *fx, double *dir, double *x,
					  double *xnext, void (*fct)(double *fx, double *x), int Nx_tmp,
					  double step_broyden_tmp)
{
	//	make copy of Jacobian
	double *Jxcopy;		Jxcopy = new double [Nx_tmp*Nx_tmp];
    copy_array(Jxcopy, Jx, Nx_tmp*Nx_tmp);
	//	get inverse Jacobian, break if singular
    int mat_inv_flag = inv_mat(Jxinv, Jxcopy, Nx_tmp);
	if (mat_inv_flag < 0)
		return SINGULAR_MATRIX_FLAG;
	//	evaluate equations at current candidate
    fct(fx, x);
	//	compute direction vector
    mat_mult(Jxinv, fx, dir, Nx_tmp, Nx_tmp, 1);
	mat_rescale(dir, dir, -1.0*step_broyden_tmp, Nx_tmp);
	//	solve for next candidate
    mat_add(xnext, x, dir, 1, Nx_tmp);
	//	write out function call diagnostics
	if (fout_xnext_broyden_diag_FLAG == ON)
		fout_xnext_broyden_diag(log_fout, step_broyden_tmp, Jx, Jxinv, dir, x, xnext, fx, Nx_tmp);
	//	free up memory
	delete[] Jxcopy;

	return 0;
}
/*	 update Jacobian in broyden_method: Set y(k) = f(x(k+1)) - f(x(k)) and set J(k+1) = J(k) + {(y(k)-J(k)s(k))s(k)'}/{s(k)'s(k)} where J(k+1) = Jxnext, f(x(k+1)) = fxnext, y(k) = y_temp */
void get_jxnext_broyden(double *fx, double *fxnext, double *xnext, double *y_temp,
						double *dir, double *Jx, double *Jxnext, int Nx_tmp,
						void (*fct)(double *fx, double *x))
{
	//	evaluate eqns at next candidate, get f(x(k+1))
    fct(fxnext, xnext);
	//	compute y(k) = f(k+1) - f(k)
    mat_add(y_temp, fxnext, fx, -1, Nx_tmp);
	//	compute (y(k) + f(k))
	double *temp1;		temp1 = new double [Nx_tmp];
	mat_add(temp1, y_temp, fx, 1, Nx_tmp);
	//	compute (y(k) + f(k))s(k)T  where s(k)T = transpose s(k)
	double *temp2;		temp2 = new double [Nx_tmp*Nx_tmp];
	mat_mult(temp1, dir, temp2, Nx_tmp, 1, Nx_tmp);
	//	compute s(k)T s(k)
	double temp3;
	temp3 = mat_mult(dir, dir, Nx_tmp);
	//	compute J(k+1) = J(k) + (y(k) + f(k))s(k)T /s(k)T s(k)
	mat_add(Jxnext, Jx, temp2, 1/temp3, Nx_tmp*Nx_tmp);
	//	write out diagnostics
	if (fout_jxnext_broyden_FLAG == ON) {
		fout_jxnext_broyden_diag(log_fout, Jx, Jxnext, xnext, fx, fxnext, y_temp, Nx_tmp);
	}
	//	free up memory
	delete[] temp1; delete[] temp2;
	return;
}
//	checks for and rescales the next guess for newton_method to with pchange percent of current. step_newton is taken as the max percentage change from x to xnext - note that this fct is limited to checks for candidate vars which have a positivity restriction.
void get_xnext_newton_dpercent(double *xnext, double *x, int Nx, double max_step_dpercent)
{
	//	make copy of original next candidate
	double *xnext_old; xnext_old = new double [Nx];
	copy_array(xnext_old, xnext, Nx);

	//	get percentage change from *x to *xnext
	double *x_dpercent; x_dpercent = new double [Nx];
	get_percentage_change(x_dpercent, x, xnext, Nx);
	//	rescale into decimal from percentage
	mat_rescale(x_dpercent, 0.01, Nx);

	//	check if any of changes from x(k) to x(k+1) exceeds the max percentage change
	double max_x_dpercent; // value of max percentage change
	int j_max_x_dpercent; // index of max percentage change
	get_max_fabs(x_dpercent, Nx, max_x_dpercent, j_max_x_dpercent);

	//	if max percentage change exceeded, rescale xnext
	if (fabs(max_x_dpercent) > max_step_dpercent) {
		double rescale_factor = fabs(max_step_dpercent)/max_x_dpercent;
		double *dx; dx = new double [Nx];
		mat_add(dx, xnext, x, -1, Nx);
		mat_add(xnext, x, dx, rescale_factor, Nx);
		delete[] dx;
	}
	//	option: write diagnostics to log
	if (fout_xnext_newton_dpercent_FLAG == ON) {
		fout_xnext_newton_dpercent(log_fout, x, xnext, xnext_old, x_dpercent, max_x_dpercent,
			j_max_x_dpercent, max_step_dpercent, Nx);
	}
	delete[] x_dpercent; delete[] xnext_old;
	return;
}
//	--------------------------------- DIAGNOSTICS ---------------------------------------------
void fout_xnext_newton_dpercent(ofstream& fout, double *x, double *xnext, double *xnext_old, double *dx_percent,
								double max_x_dpercent, int jmax_x_dpercent, double max_dpercent,
								int Nx)
{
	int Nprec = 6;

	fout << endl << "--- xnext_newton_dpercent_diag: max_percent_change_allowed=" << max_dpercent;
	fout << endl << "max_actual_percent_change=" << max_x_dpercent << ", index(max_actual)="
		<< jmax_x_dpercent << endl;

	fout << setw(12) << "x(k)=";
	fout_array(log_fout, x, 1, Nx, Nprec, 0);

	fout << setw(12) << "new:x(k+1)=";
	fout_array(log_fout, xnext, 1, Nx, Nprec, 0);

	fout << setw(12) << "old:x(k+1)=";
	fout_array(log_fout, xnext_old, 1, Nx, Nprec, 0);

	fout << setw(12) << "dx=";
	fout_array(log_fout, dx_percent, 1, Nx, Nprec, 0);

	return;
}
void fout_get_xnext_newton_powell(ofstream& fout, double *x, double *xnext, double *xnext_tmp,
								  int Nx, double loss_fx, double loss_fxnext, double loss_fxnext_tmp)
{
	int Nprec = 4;

	fout << endl << "get_xnext_newton_powell:" << endl;

	fout << setw(15) << "|fx(k)|" << setw(12) << loss_fx << setw(12) << "x(k)";
	fout_array(fout, x, 1, Nx, Nprec, 0);

	fout << setw(15) << "old:|fx(k+1)|" << setw(12) << loss_fxnext << setw(12) << "x(k+1)";
	fout_array(fout, xnext, 1, Nx, Nprec, 0);

	fout << setw(15) << "new:|fx(k+1)|" << setw(12) << loss_fxnext_tmp << setw(12) << "x(k+1)";
	fout_array(fout, xnext_tmp, 1, Nx, Nprec, 0);

	return;
}
void fout_ParamNewton(ofstream& fout, ParamNewton pm)
{
	fout << endl << "newton parameters:" << endl;
	fout << setw(16) << "eps_newton" << setw(12) << pm.eps_newton << endl;
	fout << setw(16) << "Maxiter" << setw(12) << pm.Nmaxiter << endl;
	fout << setw(16) << "step_jacob" << setw(12) << pm.step_jacob << endl;
	fout << setw(16) << "step_newton" << setw(12) << pm.step_newton << endl;
	fout << setw(16) << "Nfout" << setw(12) << pm.Nfout << endl;
	return;
}
