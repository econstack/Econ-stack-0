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
#include "MathGrunts20.h"
#include "NlinSolver20.h"
#include "Fmin20.h"

/*	-----------------------------------------------------------------------------------------
	Linemin: Solves for the mininmum of a function f along given direction *dir. Two methods 
	are used. 
		(i) solver_newton : (default method) uses the newton method and involves 
		numerical differentiation and uses gradient information. 
		(ii) solver_golden : Golden section method which does not require gradient info	
--------------------------------------------------------------------------------------------- */
//	constructor 
ParamLinemin::ParamLinemin()
{
	Nx = 0;						//	dimension of system dim(x) 
	step_newton_linemin = 0.5;	//	scalar for step size in linemin() call
	eps_linemin = 0.001;		//	tolerance for linemin() call
	percentstep_linemin = 0;	//	max change in percentage terms per step in line_min call
								//	( OFF if = 0), o/w is ON
	Nmax_iter_linemin = 50;		//	max # iterations in line_min call
	Nfout_linemin = 1;			//	# iterations between write to log in line_min call
								//	no writes to outfile if (=0)
	step_jacob = 0.001;			//	step size for numerical differentation 
	//	process diagnostics
	loss_x = 1.0;
	fx = 1.0;
	dfx = 1.0;
	d2fx = 1.0;
	time_used = 0;
	process_flag = 0;
}
//	write to outfile
void ParamLinemin::_fout(std::ofstream &fout)
{
	fout << endl << setw(15) << "ParamLinemin():"
		<< setw(12) << "step(nwt)"
		<< setw(12) << "eps"
		<< setw(12) << "pct_step"
		<< setw(12) << "step(jacob)"
		<< setw(12) << "Nmaxiter"
		<< setw(12) << "Nfout"
		<< setw(12) << "Nx" << endl;
	fout << setw(27) << this->step_newton_linemin
		<< setw(12) << this->eps_linemin
		<< setw(12) << this->percentstep_linemin
		<< setw(12) << this->step_jacob
		<< setw(12) << this->Nmax_iter_linemin
		<< setw(12) << this->Nfout_linemin
		<< setw(12) << this->Nx << endl;
	return;
}
/*	----------------------------------------------------------------------------------------
	Linemin solvers : searches for minimum along direction *dir. The initial guess, *x, MUST
	be provided along with the search direction *dir. 
	---------------------------------------------------------------------------------------- */
//	solver_newton() : searches for minimum along direction *dir using newton method
int ParamLinemin::_solver_newton(double *x, double *xmin, double *dir, double &fx, double &df, double &d2f, 
								 double (*fct)(double *), int (*fct_constraints)(double *), ofstream *fout)
{
	//	set xmin = x (prevents crashes when exiting without being able to compute next candidate x(k+1))
	copy_array( xmin, x, this->Nx );
	//	local control-flow vars
	int i_lmin = 0; // current number of iterations
	int xnext_convex_flag, rescale_flag, process_flag = 0;
	double loss_x = 1;
	
	//	candidate solutions
	double *xnext; xnext = new double [this->Nx];
	double *x0_tmp; x0_tmp = new double [this->Nx];
	//	make copy of init guess and eval system f at init guess
	copy_array( x0_tmp, x, this->Nx );
	double fx0 = fct( x0_tmp ); 

	//	f properties at candidate solution
	double dfnext, d2fnext, fxnext = 1000;

	//	scale direction vector *dir into unit vector
	double norm_dir_tmp = norm_L2( dir, this->Nx );
	mat_rescale( dir, 1.0/norm_dir_tmp, this->Nx );

	//	time use vars
	double stoptime, starttime, stoptimef, starttimef = clock();

	//	compute first derivative and 2nd derivative along dir
	d2f = get_d2fdx2_along_s( fx0, x, dir, this->Nx, fct, this->step_jacob, 
		df, DIFF_3PT_FLAG);
	fx = fx0;

	//	main loop
	double eps_dx = 0.10 * (this->eps_linemin) * fabs( 1.0 + norm_L2( x, this->Nx ) );
	while ( ( loss_x > eps_dx ) && ( fabs(df) > (this->eps_linemin) ) ) {

		//	start iteration clock
		starttime = clock();

		//	Break if exceeded max number of iterations
		if ( i_lmin >= this->Nmax_iter_linemin ) {
			this->_fout_error( log_fout, x, xnext, i_lmin, fx, 
				"Error! Linemin_newton(): Exceeded max # iterations. Exiting. " );
			copy_array( xmin, x, this->Nx );
			this->process_flag = EXCEEDED_MAXITER_FLAG;
			break;
		}

		/*	----------------------------------------------------------------------------
			step 1. Get next candidate for xmin; check if it satisfies the following:
			(i) is in constraint set, (ii) (optional) is within percentage change limit,
		-------------------------------------------------------------------------------- */
		//	Get candidate x(k+1)
		mat_add( xnext, x, dir, -(this->step_newton_linemin)*df/d2f, this->Nx );
		//	verify candidate x(k+1) is in constraint set, rescale if needed
		rescale_flag = rescale_xnext(x, xnext, this->Nx, fct_constraints);
		if ( rescale_flag < 0 ) {
			this->_fout_error( log_fout, x, xnext, i_lmin, fx, 
				"Error! Linemin(): Cannot find candidate in constraint set. Exiting." );
			this->process_flag = CONSTRAINT_SET_VIOLATION_FLAG;
			break;
		}
		//	option: check if x(k+1) is within percent change limit, contract step if needed
		if ( (this->percentstep_linemin) != 0 ) {
			int xnext_percent_flag = this->_get_xnext_dpercent( xnext, x );
			//	TODO write error log
			if ( xnext_percent_flag < 0 ) {

			}
		}

		/*	----------------------------------------------------------------------------
			step 2. Get 2nd derivative at candidate x(k+1); Option : Check if candidate 
			if locally convex, rescale if needed. 
		-------------------------------------------------------------------------------- */
		//	compute 2nd derivative at x(k+1)
		d2fnext = get_d2fdx2_along_s(fxnext, xnext, dir, this->Nx, fct, this->step_jacob,
			dfnext, DIFF_3PT_FLAG);
		//	option: check for local convexity at next candidate, rescale if needed
		if ( (this->flag_locally_convex == ON) && (d2fnext < 0) ) {
			//	find locally convex xnext and get f, df, d2f at new candidate, xnext
			xnext_convex_flag = this->_get_xnext_convex( x, xnext, dir, fxnext, dfnext, d2fnext, fct );
			//	break if unable to rescale next guess st it is locally convex
			if ( xnext_convex_flag < 0 ) {
				this->_fout_error( log_fout, x, xnext, i_lmin, fx, 
					"Error! Linemin_newton(): Candidate not locally convex. Exiting. " );
				copy_array( xmin, x, this->Nx );
				this->process_flag = -88;
				break;
			}
		}

		//	compute change in x in normed space
		this->loss_x = norm_sup( x, xnext, this->Nx );
		eps_dx = 0.10 * (this->eps_linemin) * fabs( 1.0 + norm_L2( x, this->Nx ) );
		
		//	option : write out iteration diagnostics data to file
		if ( this->Nfout_linemin != 0 ) {
			if ( (i_lmin%this->Nfout_linemin) == 0 ){
				stoptime = clock();
				this->time_used = ( stoptime - starttime ) / 1000.0;
				if ( fout == 0 ) {	//	write to default logdata outfile
					this->_fout_iter( log_fout, x, xnext, i_lmin, fx, fxnext );
				} else { // write to specified logdata outfile
					this->_fout_iter( *fout, x, xnext, i_lmin, fx, fxnext );
				}	
			}
		}
		//	check for x-xtuck and/or convergence 
		if ( ((this->loss_x) < eps_dx ) || ( fabs(df) < (this->eps_linemin)) )	{
			stoptimef = clock();
			this->time_used = ( stoptimef - starttimef ) / 1000.0;
			copy_array( xmin, xnext, this->Nx );
			//	successful convergence
			if ( fabs(df) < (this->eps_linemin) ) 
				this->process_flag = i_lmin;
			else //	x stuck but non-zero derivative
				this->process_flag = X_STUCK_FLAG; 
			break;
		}

		// increment iteration counter and update guess and function at new guess
		++i_lmin;
		copy_array( x, xnext, this->Nx );
		df = dfnext;
		d2f = d2fnext;
		fx = fxnext;

	}

	//	update and return process_flag and return solution xmin ==================================
	if ( (process_flag >= 0) || (process_flag == X_STUCK_FLAG) ) {
		process_flag = i_lmin;
		//	update if improvement upon init guess
		if (fx0 > fxnext) {
			copy_array( xmin, xnext, this->Nx );
			fx = fxnext;
			df = dfnext;
			d2f = d2fnext;
		} else { //	 o/w set solution to init guess
			copy_array( xmin, x0_tmp, this->Nx );
			fx = fx0;
		}
	}
	//	write out function call summary
	stoptimef = clock();
	this->time_used = ( stoptimef - starttimef ) / 1000.0;
	if ( fout == 0 ) {
		this->_fout_exit( log_fout, x0_tmp, xmin, fx0, fx );
	} else {
		this->_fout_exit( *fout, x0_tmp, xmin, fx0, fx );
	}
	//	restore original init guess
	copy_array( x, x0_tmp, this->Nx );
	//	free heap
	delete[] xnext; delete[] x0_tmp;
	return process_flag;
}
//	solver_golden() : solves for min along direction *dir by golden section method
int ParamLinemin::_solver_golden(double *x, double *xmin, double *dir, double &fx, 
								 double (*fct)(double *), int (*fct_constraints)(double *), ofstream *fout)
{
	int Nx = this->Nx;
	//	declare local vars for x(i), f(x(i)), i = {a,b,c,d}
	double *xa;	xa = new double [Nx];
	double *xb;	xb = new double [Nx];
	double *xc;	xc = new double [Nx];
	double *xd;	xd = new double [Nx];
	double fa, fb, fc, fd;
	//	control-flow vars
	int j_lmin = 0;
	double loss = 1.0;

	//	rescale direction vector *dir into unit vector
	double rescale_dir = norm_L2( dir,  Nx);
	mat_rescale( dir, dir, 1.0/rescale_dir, Nx );

	/*	--------------------------------------------------------------------------
		step 1: Find points x(a)<x(c)<x(b) st f(x(a)), f(x(b)) > f(x(c))
	----------------------------------------------------------------------------- */
	//	labels and control-flow vars for step 1
	int NOT_FOUND = 1;
	int d_LESSTHAN_b = 1;
	int d_GREATERTHAN_b = 0;
	int xnext_flag;		//	flag for d<b or d>b for candiate next point d

	double starttime, stoptime, stoptimef, starttimef = clock();

	//	find initial triple x(a), x(c) such that f(x(a)), f(x(c)) > f(x(b))
	this->_solver_golden_get_init_triple( x, dir, xa, xb, xc, fa, fb, fc, fct, fct_constraints );
	//	save initial guess data 
	double *x0; x0 = new double [Nx];
	copy_array( x0, xb, Nx );
	double fx0 = fb;

	//	main loop ------------------------------------
	while ( (process_flag!=FAILURE_FLAG) && (loss>epsLinemin) )	{

		starttime = clock();

		/*	-------------------------------------------------------------------------
			Choose next candidate xd and evaluate f(xd). 
			(i)	If current candidate x(b) is closer to x(a), then set next candidate
					x(d) = 0.5*(x(b)+x(c))
				such that x(d) > x(b)
			(ii) If current candidate x(b) is closer to x(c), then set next candidate
				to midpoint between x(a) and x(b)
					x(d) = 0.5*(x(a)+x(b))
				such that x(d) < x(b) 
		----------------------------------------------------------------------------- */
		if ( norm_L2( xb, xa, Nx ) < norm_L2( xc, xb, Nx ) )	{	
			double temp1 = norm_L2( xb, xa, Nx );
			double temp2 = norm_L2( xc, xb, Nx );
			mat_add( xd, xc, xb, 1.0, Nx );
			mat_rescale( xd, xd, 0.5, Nx );
			xnext_flag = d_GREATERTHAN_b;
		} else {								
			mat_add( xd, xa, xb, 1.0, Nx );
			mat_rescale( xd, xd, 0.5, Nx);
			xnext_flag = d_LESSTHAN_b;
		}
		//	evaluate f at new candidate
		fd = fct( xd );

		/*	------------------------------------------------------------------------
			Set new triple ( xa, xb, xc). 
			(i)		If x(d) < x(b) along direction *dir and f(x(d)) >= f(x(b)), replace 
					x(a,b,c) with x(d,b,c)
			(ii)	If x(d) < x(b) but f(x(d)) < f(x(b)), replace x(a,b,c) with x(a,d,b)
			(iii)	If x(d) > x(b) along direction *dir and f(x(d)) <= f(x(b)), replace 
					x(a,b,c) with x(b,d,c)
			(iv)	Else if x(d) <= x(b) but f(x(d)) > f(x(b)), replace x(a,b,c)
					with x(a,b,d)
		---------------------------------------------------------------------------- */
		if ( (xnext_flag == d_LESSTHAN_b) && (fd >= fb) ) {	//	case (i)
			copy_array( xa, xd, Nx );
			fa = fd;
		} else if ( (xnext_flag == d_LESSTHAN_b) && (fd < fb) ) { // case(ii) 
			copy_array( xc, xb, Nx );
			copy_array( xb, xd, Nx );
			fc = fb;
			fb = fd;
		} else if ( (xnext_flag == d_GREATERTHAN_b) && (fd<=fb) ) { // case(iii)
			copy_array( xa, xb, Nx );
			copy_array( xb, xd, Nx );
			fa = fb;
			fb = fd;
		} else { //	case(iv)
			copy_array( xc, xd, Nx );
			fc = fd;
		}


		//	update control-flow variables
		loss = norm_L2( xc, xa, Nx );
		loss += pow( fabs( fa-fc ), 2 );
		++j_lmin;
		//	check stopping condition, break if exceeded max # iterations
		if ( j_lmin >= (this->Nmax_iter_linemin) ) {
			this->_fout_error( log_fout, xb, xd, j_lmin, fb, 
				"Error! Golden_section_Linemin : exceed max number of iterations. Exiting. " );
			this->process_flag = FAILURE_FLAG;
			break;
		}
		//	write out iteration summary 
		if ( (this->Nfout_linemin) != 0 ) {
			if ( j_lmin%(this->Nfout_linemin) == 0 ) {
				stoptime = clock();
				this->time_used = ( stoptime - starttime ) / 1000;
				if ( fout != 0 ) {
					this->_fout_iter( *fout, xb, xd, j_lmin, fb, fd );
				} else {
					this->_fout_iter( log_fout, xb, xd, j_lmin, fb, fd );
				}
			}
		}

	}	

	//	write out fct call summary ========================================================
	stoptimef = clock();
	this->time_used = ( stoptimef - starttimef ) / 1000;
	if ( (this->process_flag) != FAILURE_FLAG ) 
		this->process_flag = j_lmin;

	if ( fout != 0 ) {
		this->_fout_exit( *fout, x0, xb, fx0, fb );
	} else {

	}
	//	return solution in *xmin
	copy_array( xmin, xb, Nx );

	delete[] xa; delete[] xb; delete[] xc; delete[] xd; delete[] x0;
	return (this->process_flag);
}
/*	------------------------------------------------------------------------------------------------
	Golden section linemin search requires an init triple { x(a), x(b), x(c) } along direction *dir 
	such that:
			f(x(a)), f(x(c)) > f(x(b)) 
	where the initial value of x(b) = x(0) MUST be given
	------------------------------------------------------------------------------------------------ */
int ParamLinemin::_solver_golden_get_init_triple(double *x0, double *dir, double *xa, double *xb, double *xc, 
												 double &fa, double &fb, double &fc, double (*fct)(double *), 
												 int (*fct_constraints)(double *), std::ofstream *fout)
{
	//	control-flow vars
	int j_iter = 0;
	int NOT_FOUND = -1;
	int found_xa_flag = NOT_FOUND;
	int found_xc_flag = NOT_FOUND;
	int process_flag = 0;
	//	set x(b) = x(0) and evalute system f at x(b)
	copy_array( xb, x0, this->Nx );
	fb = fct( xb );

	while ( (found_xa_flag == NOT_FOUND) || (found_xc_flag == NOT_FOUND) ) {

		//	select points x(a) and x(c) from x(0) along direction *dir, evaluate f at candidates
		if ( found_xa_flag == NOT_FOUND ) {
			mat_add( xa, xb, dir, pow( this->golden_alpha, j_iter ) * this->golden_delta, this->Nx );
			fa = fct( xa );
		}
		if ( found_xc_flag == NOT_FOUND ) {
			mat_add( xc, xb, dir, -pow( this->golden_alpha, j_iter ) * this->golden_delta, this->Nx );
			fc = fct( xc );
		}

		//	 check whether f(a), f(c) > f(b), upate flags for updating x(a), x(c)
		if ( fa > fb ) {
			found_xa_flag = 1 - NOT_FOUND;
		}
		if ( fc > fb ) {
			found_xc_flag = 1 - NOT_FOUND;
		}	

		//	increment iteration counter, break if exceed max # iterations
		++j_iter;
		if ( j_iter >= (this->Nmax_iter_linemin) ) {
			process_flag = NOT_FOUND;
			break;
		}

	}	
	//	option : write to outfile
	if ( fout != 0 ) {
		int Nprec = 4;
		( *fout ) << endl << "Golden_section_solver() : process_flag = " << process_flag << endl;
		( *fout ) << setw(12) << "f(x(a))=" << setw(12) << fa;
		( *fout ) << "x(a)=";
		fout_array( *fout, xa, 1, this->Nx, Nprec, OFF );
		( *fout ) << setw(12) << "f(x(b))=" << setw(12) << fb;
		( *fout ) << "x(b)=";
		fout_array( *fout, xb, 1, this->Nx, Nprec, OFF );
		( *fout ) << setw(12) << "f(x(c))=" << setw(12) << fc;
		( *fout ) << "x(c)=";
		fout_array( *fout, xc, 1, this->Nx, Nprec, OFF );
	}
	return process_flag;
}
int ParamLinemin::_get_xnext_dpercent(double *xnext, double *x, std::ofstream *fout)
{
	int flag = get_xnext_dpercent( xnext, x, this->Nx, this->percentstep_linemin, fout );
	return flag;
}
int ParamLinemin::_get_xnext_convex(double *xnext, double *x, double *dir, double &fnext, double &dfnext, 
									double &d2fnext, double (*fct)(double *x), std::ofstream *fout)
{
	//	control-flow vars
	int Nmaxiter_tmp = 10;
	int flag = -88;
	//	local vars
	double rescale_factor = 0.5;

	//	make copy of original candidate
	double *xnext_old; xnext_old = new double [this->Nx];
	copy_array( xnext_old, xnext, this->Nx );

	//	get change vector
	double *dx; dx = new double [this->Nx];
	mat_add( dx, xnext, x, -1, this->Nx );
	mat_rescale( dx, rescale_factor, this->Nx );

	//	search for candidate that is locally convex
	for ( int i=0; i<Nmaxiter_tmp; ++i ) {

		//	compute new candidate xnext
		mat_add( xnext, x, dx, 1.0, this->Nx);
		//	compute derivative info at new candidate
		d2fnext = get_d2fdx2_along_s( fnext, xnext, dir, this->Nx, fct, this->step_jacob,
			dfnext, DIFF_3PT_FLAG );
		//	case: xnext is locally convex - update info and break
		if ( d2fnext > 0 ) {
			//	candidate is automatically updated in call get_d2fdx2_along_s(), 
			flag = 0;
			break;
		}
		//	update new change vector
		mat_rescale( dx, rescale_factor, this->Nx ); 
	}

	//	option : write out diagnostics 
	if ( fout != 0 ) {
		(*fout) << endl << "_get_xnext_convex(): d2fnext = " << d2fnext;
		this->_fout_iter( *fout, x, xnext, -88, this->fx, fnext );
	}
	//	free heap
	delete[] xnext_old; delete[] dx;
	return flag;
}
//	write out diagnostics and summary data
void ParamLinemin::_fout_error(ofstream &fout, double *x, double *xnext, const int& i_lmin, 
							   const double& fx, string label)
{
	int Nprec = 6;
	//	write header and summary stats (to console and to logfile)
	cout << endl << label << ", flag=" << this->process_flag;
	cout << endl << "ParamLinemin::_solver(): iter= " << i_lmin << ", |fx|=" << fx;
	fout << endl << label << ", flag=" << this->process_flag;
	fout << endl << "ParamLinemin::_solver(): iter= " << i_lmin 
		<< ", |fx|=" << fx << ", time_used=" << this->time_used << endl;
	//	write out current candidate 
	fout << setw(16) << "x(k)"; 
	fout_array( fout, x, 1, this->Nx, Nprec, OFF );
	//	write out next candidate
	fout << setw(16) << "x(k+1)"; 
	fout_array( fout, xnext, 1, this->Nx, Nprec, OFF );
	return;

}
void ParamLinemin::_fout_exit(std::ofstream &fout, double *x0, double *xstar, const double &fx0, 
							  const double &fxstar)
{
	int Nprec = 6;
	//	write header 
	fout << endl << "Exit ParamNlin::_solver(): iter = " 
		<< this->process_flag << "/" << this->Nmax_iter_linemin 
		<< ", time_used = " << this->time_used 
		<< ", f(x0)= " << fx0 << ", f(x*)" << fxstar << endl;
	//	write out initial data: guess x0 and f(x0)
	fout << setw(14) << "x0=";
	fout_array( fout, x0, 1, this->Nx, Nprec, OFF );
	//	write out solution (or latest candidate) data 
	fout << setw(14) << "x*=";
	fout_array( fout, xstar, 1, this->Nx, Nprec, OFF );
	return;
}
void ParamLinemin::_fout_iter(std::ofstream &fout, double *x, double *xnext, const int &i_lmin, 
							  const double &fx, const double &fxnext)
{
	int Nprec = 6;
	//	write header and summary stats
	cout << endl << "ParamNlin::_solver(): iter= " << i_lmin << ", loss(x)="
		<< loss_x << ", f(k)=" << fx << ", f(k+1)=" << fxnext;
	fout << endl << "ParamNlin::_solver(): iter= " << i_lmin << ", loss(x)="
		<< loss_x << ", f(k)=" << fx << ", f(k+1)=" << fxnext << ", time_used=" << time_used << endl;
	//	write out current candidate 
	fout << setw(16) << "x(k)"; 
	fout_array( fout, x, 1, this->Nx, Nprec, 0 );
	//	write out next candidate
	fout << setw(16) << "x(k+1)"; 
	fout_array( fout, xnext, 1, this->Nx, Nprec, 0 );
	return;
}
//	---------------------------------------------------------------------------------------------
//	Constructor
ParamFxmin::ParamFxmin()
{
	step_jacob = 0.001;			//	step size for numerical Jacobian 
	eps_fmin = 0.001;			//	tolerance for fmin
	Nmax_iter_fmin = 200;		//	max # iterations for fmin 
	Nfout_fmin = 1;				//	# iterations between write to log for fmin
	linemin_newton_flag = ON; 	//	switch for type of line_min (golden section vs newton)
	BFGS_flag = ON;				//	switch for type of fmin solver, default = BFGS, alt = DFP 
								//	differs only in update of hessian
	//	process diagnostics
	loss_x = 1.0;
	fx = 1.0;
	process_flag = 0;
	time_used = 0;
}
//	---------------------------------------------------------------------------------------------
void ParamFxmin::_fout(std::ofstream &fout)
{
	fout << setw(15) << "ParamFxmin: " << setw(12) << "step(jacob)"
		<< setw(12) << "eps(fmin)" << setw(12) << "Max#iter" 
		<< setw(12) << "Nfout(fmin)" << endl;
	fout << setw(27) << step_jacob << setw(12) << eps_fmin
		<< setw(12) << Nmax_iter_fmin << setw(12) << Nfout_fmin << endl;
	return;
}
int ParamFxmin::_solver_BFGS(double *x, double *xmin, double (*fct)(double *x), 
							 int (*fct_constraints)(double *x), ofstream *fout)
{
	//	init xmin to be equal to init guess 
	copy_array( xmin, x, this->pm_lmin.Nx );

	//	store old global stepJacob and set current as global
	double step_jacob_old = stepJacob;
	stepJacob = this->step_jacob;
	//	rescale tolerance for # of dimensions, Nx
	eps_fmin /= pow( 1.0*(this->pm_lmin.Nx), 0.5);
	//	local flags, counters and loss vars
	int process_flag = 0, process_flag_linemin = 0;
	int hessian_inv_flag = 0;
	double norm_Jxnext = 1.0, loss_x = 1.0, loss = 1.0;
	//	init iteration counter
	int i_fmin = 0;
	//	computation arrays, hessian, gradient and aux vars used in main loop
	int Nx = this->pm_lmin.Nx;
	double *Hx; Hx = new double [Nx*Nx];		// hessian at current guess
	double *Hxinv; Hxinv = new double [Nx*Nx];	// inv of hessian at current guess
	double *Hxnext; Hxnext = new double [Nx*Nx]; //	hessian at next guess
	double *Jx;	Jx = new double [Nx];			// gradient at current guess
	double *Jxnext; Jxnext = new double [Nx];	// gradient at next guess
	double *xnext; xnext = new double [Nx];		//	next guess
	double *z; z = new double [Nx];				//	aux var - used in updating hessian
	double *y; y = new double [Nx];				//	aux var - used in updating hessian
	double *s; s = new double [Nx];				//	current search direction
	const_array( xnext, Nx );	// init candidate to zeroes - by construction, x0 should never be an 
								// array of zeros so that the init check for stopping conditions is not satd 

	//	make copy of initial guess, x0, and jacobian and f at init guess, for diagnostics
	double *x0_tmp;	x0_tmp = new double [Nx];
	copy_array( x0_tmp, x, Nx );
	double fx0;
	//	compute f(x0) and jacobian, J(x0) 
	get_gradient( Jx, x, fx0, Nx, fct, stepJacob );

	//	set initial guess for H(k) to identity mat
	make_identity_matrix( Hx, Nx );

	//	clocking and timeuse vars
	double starttime, stoptime, stoptimef, starttimef = clock();
	//	declare computation parameters linemin_newton
	double dfdx, d2fdx;
	double fxnext, fx = fx0;

	//	main loop ====================================================================
	while ( loss> (this->eps_fmin) ) {

		//	start iteration clock
		starttime = clock();

		/*	----------------------------------------------------------------------- 
			step 0: Check stopping/exit conditions
			(i) check if x is stuck at suboptimal, (ii) for successful convergence
			(iii) whether max # iterations have been exceeded. 
		--------------------------------------------------------------------------- */
		//	check stopping conditions for convergence or x-stuck 
		loss_x = norm_sup( x, xnext, Nx );
		norm_Jxnext = norm_L2( Jx, Nx );	//	changed from Jxnext to Jx - sept26/05
		if (loss_x < 0.1*(eps_fmin + eps_fmin*norm_sup(x, Nx))) {
			//	check for x stuck at suboptimal point
			if ( norm_Jxnext > eps_fmin*(1.0+fabs(fxnext)) ) {
				this->_fout_error( log_fout, x, xnext, i_fmin, fx, 
					"Error! ParamFxmin_solver(): x-stuck. Exiting. " );
				this->process_flag = X_STUCK_FLAG;
				break;
			}
			//	break for successful convergence to xmin
			this->process_flag = i_fmin;
			break;
		}

		//	check for exceed max # iterations
		if ( i_fmin >= this->Nmax_iter_fmin ) {
			this->_fout_error( log_fout, x, xnext, i_fmin, fx, 
				"Error! ParamFxmin_solver(): Exceeded max # iterations. Exiting." );
			this->process_flag = EXCEEDED_MAXITER_FLAG;
			break;
		}
		//	increment fmin iteration counter
		++i_fmin;

		/*	----------------------------------------------------------------------- 
			step 1: 
			(i) Solve for search direction s(k): where H(k)s(k) = - J(x(k)). 
				Involves the inversion of current hessian. 
			(ii) Find min along search direction s(k).
		--------------------------------------------------------------------------- */
		//	compute inverse Hessian at current candidate, break if singular
		hessian_inv_flag = inv_mat( Hxinv, Hx, Nx );
		if ( hessian_inv_flag < 0 ) {
			this->_fout_error( log_fout, x, xnext, i_fmin, fx, 
				"Error! ParamFxmin::_solver_BFGS(): Singular hessian. Exiting." );
			this->process_flag = hessian_inv_flag;
			break;
		}
		//	compute search direction and rescale into unit vector
		mat_mult( Hxinv, Jx, s, Nx, Nx, 1 );
		mat_rescale( s, -1.0, Nx );
		//	Find min along search direction
		if ( linemin_newton_flag == ON ) {	//	default : use newton_method()
			this->pm_lmin.process_flag = this->pm_lmin._solver_newton( x, xnext, s, fx, dfdx, 
				d2fdx, fct, fct_constraints );
		} else { //	use golden section 
			this->pm_lmin.process_flag 
				= this->pm_lmin._solver_golden( x, xnext, s, fx, fct, fct_constraints );
		}

		/*	----------------------------------------------------------------------- 
			step 2: 
				(i) Check that candidate x{k+1} is locally convex, rescale x(k+1)
				if needed. 
				(ii) Update gradient and hessian
		--------------------------------------------------------------------------- */
		//	Test for non-local convex 
		if ( (this->pm_lmin.process_flag) == -88 ) {
			this->_fout_error( log_fout, x, xnext, i_fmin, fx, 
				"Error! ParamFxmin::_solver(): Cannot find locally convex candidate. Exiting." );
			this->process_flag = this->pm_lmin.process_flag;
			break;
		}
		//	update gradient and hessian (either using BFGS=default or DFP)
		this->_solver_update_gradient( y, z, x, xnext, Jx, Jxnext, fxnext, fct );
		if ( (this->BFGS_flag) == ON) {	//	BFGS updating rule = default
			this->_solver_update_hessian_BFGS( Hx, Hxnext, x, y, z );
		} else {	//	DFP updating rule
			this->_solver_update_hessian_DFP( Hx, Hxnext, x, y, z );
		}

		/*	----------------------------------------------------------------------- 
			step 3: write out iteration diagnostics, update current guess and data
				(i) option: write to logdata outfile
				(ii) update current candidate and its associated data
		--------------------------------------------------------------------------- */
		//	option : write out general iteration diagnostics
		if ( this->Nfout_fmin != 0 ) {
			stoptime = clock();
			this->time_used = ( stoptime - starttime ) / 1000.0;
			if ( (i_fmin%this->Nfout_fmin) == 0 ) {
				if ( fout == 0 ) {	//	write to default global logdata outfile
					this->_fout_iter( log_fout, x, xnext, i_fmin, fx, fxnext );
				} else { // write to specified logdata outfile
					this->_fout_iter( *fout, x, xnext, i_fmin, fx, fxnext );
				}
			}
		}
		//	update Jacobian, Hessian and candidate solution
		copy_array( Hx, Hxnext, Nx*Nx );
		copy_array( x, xnext, Nx );
		copy_array( Jx, Jxnext, Nx );
		fx = fxnext;
	}

	//	save xmin and write out function call diagnostics ==========================================
	copy_array(xmin, xnext, Nx);
	this->fx = fx;
	stoptimef = clock();
	this->time_used = ( stoptimef - starttimef ) / 1000.0;
	if ( process_flag >= 0 ) 
		this->process_flag = i_fmin;
	//	write out data to logdata outfile
	if ( fout == 0 ) {
		this->_fout_exit( log_fout, x0_tmp, xmin, fx0, this->fx );
	} else {
		this->_fout_exit( *fout, x0_tmp, xmin, fx0, this->fx );
	}

	//	free up memory
	delete[] x0_tmp; delete[] Hxinv; delete[] Hxnext; delete[] z;
	delete[] Jx; delete[] Jxnext; delete[] xnext;
	delete[] y; delete[] s;
	//	restore stepJacob
	stepJacob = step_jacob_old;
	return process_flag;
}
//	alterative fmin using DFP method which only differs in the updating of the hessian
int ParamFxmin::_solver_DFP(double *x, double *xmin, double (*fct)(double *), 
							int (*fct_constraints)(double *), std::ofstream *fout)
{
	this->BFGS_flag = OFF;
	this->_solver_BFGS( x, xmin, fct, fct_constraints, fout );
	return process_flag;
}
//	solver utils -------------------------------------------------------------------------
void ParamFxmin::_solver_update_gradient(double *y, double *z, double *x, double *xnext, 
										 double *Jx, double *Jxnext, double &fxnext, 
										 double (*fct)(double *), std::ofstream *fout)
{
	mat_add( z, xnext, x, -1.0, this->pm_lmin.Nx );				//	find z(k) = x(k+1) - x(k)
	get_gradient( Jxnext, xnext, fx, this->pm_lmin.Nx, fct );	//	get J(x(k+1))
	mat_add( y, Jxnext, Jx, -1.0, this->pm_lmin.Nx );			//	find y(k) = J(x(k+1)) - J(x(k))
	//	option write out diagnostics 
	if ( fout != 0 ) {
		int Nprec = 4;
		(*fout) << endl << "ParamFxmin::solver_update_gradient. Nx = " << this->pm_lmin.Nx << endl;
		(*fout) << setw(12) << "z(k)=";
		fout_array( *fout, z, 1, this->pm_lmin.Nx, Nprec, OFF );
		(*fout) << setw(12) << "y(k)=";
		fout_array( *fout, y, this->pm_lmin.Nx, this->pm_lmin.Nx, Nprec );
	}
	return;
}
void ParamFxmin::_solver_update_hessian_BFGS(double *Hx, double *Hxnext, double *x, 
											 double *y, double *z, std::ofstream *fout)
{
	int Nx = this->pm_lmin.Nx;
	//	declare local vars
	double *H1;	H1 = new double [Nx*Nx];
	double *H1a; H1a = new double [Nx];
	double *H1b; H1b = new double [Nx];
	double *H3;	H3 = new double [Nx*Nx];
	double H2, H4;
	//	compute H1a(k) = H(k)z(k)
	mat_mult( Hx, z, H1a, Nx, Nx, 1 );
	//	compute H1b(k) = tranpose(H1a(k)) = z(k)H(k)
	transpose_mat( H1b, H1a, Nx, 1 );
	//	compute H1(k) = H1a(k)H1b(k) = H(k)z(k)z(k)H(k)
	mat_mult( H1a, H1b, H1, Nx, 1, Nx );
	//	compute H2(k) = z(k)H(k)z(k)
	H2 = mat_mult( z, H1a, Nx );
	//	compute H1(k) = H1(k)/H2(k) = H(k)z(k)z(k)H(k)/z(k)H(k)z(k)
	mat_rescale( H1, H1, 1.0/H2, Nx*Nx );
	//	compute H3(k) = y(k)y(k)
	mat_mult( y, y, H3, Nx, 1, Nx );
	//	compute H4(k) = y(k)z(k)
	H4 = mat_mult( y, z, Nx );
	//	compute H3(k) = H3(k)/H4(k) = y(k)y(k)/y(k)z(k)
	mat_rescale( H3, H3, 1.0/H4, Nx*Nx );
	//	update hessian, H(k+1) = H(k) + H3(k) - H1(k)
	mat_add( Hxnext, H3, H1, -1.0, Nx*Nx );
	mat_add( Hxnext, Hx, Hxnext, 1.0, Nx*Nx );
	//	write out diagnostics
	if ( fout != 0 ) {
		this->_solver_update_hessian_fout( *fout, Hx, Hxnext );
	}
	//	free mem
	delete[] H1; delete[] H1a; delete[] H1b; delete[] H3;
	return;
}
void ParamFxmin::_solver_update_hessian_DFP(double *Hx, double *Hxnext, double *x, 
											double *y, double *z, std::ofstream *fout)
{
	int Nx = this->pm_lmin.Nx;
	//	get outer product of z X z where X = outer product
	double *M1; M1 = new double [Nx*Nx];
	mat_mult( z, z, M1, Nx, Nx );
	//	get dot product of z y
	double M2_tmp;
	M2_tmp = mat_mult( z, y, Nx );
	mat_rescale( M1, M1, 1.0/M2_tmp, Nx*Nx );
	//	compute H*y
	double *M3;	M3 = new double [Nx];
	mat_mult( Hx, y, M3, Nx, Nx, 1 );
	//	get outer product of (H*y) X (H*y)
	double *M4; M4 = new double [Nx*Nx];
	mat_mult( M3, M3, M4, Nx, Nx );
	//	compute y*H*y
	double M5_tmp;
	double *M5; M5 = new double [Nx];
	mat_mult( y, Hx, M5, 1, Nx, Nx );
	M5_tmp = mat_mult( M5, y, Nx );
	mat_rescale( M4, M4, 1.0/M5_tmp, Nx*Nx );
	//	putting the updating additive terms together
	double *dH; dH = new double [Nx*Nx];
	mat_add( dH, M1, M4, -1.0, Nx*Nx );
	//	update hessian
	mat_add( Hxnext, Hx, dH, 1.0, Nx*Nx );
	//	write out diagnostics
	if ( fout != 0 ) {
		this->_solver_update_hessian_fout( *fout, Hx, Hxnext );
	}
	//	free heap
	delete[] M1; delete[] M3; delete[] M4; delete[] M5; delete[] dH;
	return;
}
//	writes out diagnostics for _solver_update_hessian to outfile
void ParamFxmin::_solver_update_hessian_fout(std::ofstream &fout, double *Hx, double *Hxnext)
{
	int Nprec = 4;
	int Nx = this->pm_lmin.Nx;
	fout << endl << "ParamFxmin_solver(): update_hessianBFGS. Nx = " << Nx << endl;
	fout << "H(k)";
	fout_array( fout, Hx, Nx, Nx, Nprec );
	fout << "H(k+1)";
	fout_array( fout, Hxnext, Nx, Nx, Nprec );
	return;
}
//	function miniminization in 1 dimension ----------------------------------------------
int ParamFxmin::_solver_R1(const double &x0, double &xmin, double& fx, double& df, double& d2f,
						   double (*fct)(double), int (*fct_constraints)(double), std::ofstream *fout)
{
	//	store and later restore global step_Jacob
	double step_jacob_old = stepJacob;
	stepJacob = this->step_jacob;
	double step_df = stepJacob;

	//	control-flow vars
	int xnext_flag, constraint_flag = 0;
	//	init loss vars
	double loss_x = 1, loss_df = 1, loss = 1;
	//	init iteration counter
	int j_fmin = 0;

	//	store initial guess, x0, and eval the system f at x0
	double x = x0;
	double fx0 = fct( x0 );
	fx = fx0;
	double xnext, fxnext;

	//	clock vars
	double starttime, stoptime, stoptimef, starttimef = clock();

	//	main function call ==============================================================
	while ( loss > (this->eps_fmin) ) {

		starttime = clock();

		//	step 1: get next candidate x(k+1)
		df = get_dfdx( fx, x, fct, step_df );
		d2f = get_d2fx2( fx, x, fct, step_df );
		xnext = x - (this->pm_lmin.step_newton_linemin) * df / d2f;

		//	check if xnext in constraint set, rescale if needed
		constraint_flag = fct_constraints( xnext );
		if ( constraint_flag > 0)
			xnext_flag = rescale_xnext( x, xnext, fct_constraints );
		//	break if unable to find candidate in constraint set
		if ( xnext_flag < 0 ) {
			log_fout << "Error! Fxmin_solverR1(): Cannot find candidate in constraint set. Exiting." << endl;
			this->_solver_R1_fout_exit( log_fout, x, xnext, fx, -88 );
			this->process_flag = CONSTRAINT_SET_VIOLATION_FLAG;
			break;
		}


		//	step 2: check stopping condition
		fxnext = fct( xnext );
		loss_x = fabs( xnext - x );
		loss_df = fabs( df );
		if ( (loss_x / (1.0+fabs(x)) ) < (this->eps_fmin) ) {
			//	check for non-optimal solution
			if ( loss_df > (this->eps_fmin) ) {
				log_fout << endl << "Error fmin_R1, non-optimal solution" << endl;
				xmin = xnext;
				process_flag = X_SUBOPTIMAL_FLAG;
				break;
			} else { //	break for successful convergence
				xmin = xnext;
				this->process_flag = j_fmin;
				break;
			}
		}
		//	increment iteration flag
		++j_fmin;
		//	check if exceeded max # iterations
		if ( j_fmin >= (this->Nmax_iter_fmin) ) {
			log_fout << endl << "Exceeded max #iter, fmin_R1" << endl;
			xmin = xnext;
			process_flag = EXCEEDED_MAXITER_FLAG;
			break;
		}
		//	write out iteration diagnostics
		if ( (this->Nfout_fmin) !=0) {
			if ( ( j_fmin%(this->Nfout_fmin)) == 0 ) {
				stoptime = clock();
				this->time_used = ( stoptime - starttime ) / 1000.0;
				if ( fout != 0 ) {
					this->_solver_R1_fout_iter( *fout, x, xnext, fx, fxnext, j_fmin );
				} else {
					this->_solver_R1_fout_iter( log_fout, x, xnext, fx, fxnext, j_fmin );
				}
			}
		}
		//	update candidate
		x = xnext;
		fx = fxnext;
	}
	//	write out function call summary ==============================================
	stoptimef = clock();
	this->time_used = ( stoptimef - starttimef ) / 1000.0;
	if ( fout != 0 ) {	//	write summary data to specified outfile
		this->_solver_R1_fout_exit( *fout, x0, xmin, fx0, fx );
	} else {	//	write summary data to default outfile
		this->_solver_R1_fout_exit( log_fout, x0, xmin, fx0, fx );
	}

	//	restore old step_Jacob
	stepJacob = step_jacob_old;
	return (this->process_flag);
}
void ParamFxmin::_solver_R1_fout_iter(std::ofstream &fout, const double &x, const double &xnext, 
									  const double &fx, const double &fxnext, const int& j_fmin )
{
	fout << endl << "ParamFxmin::_solver_R1(): iter = " << j_fmin 
		<< ", time_used = " << this->time_used << endl;
	fout << setw(12) << "x(k)" << setw(12) << "x(k+1)" << setw(12) << "f(k)" << setw(12) << "f(k+1)" << endl;
	fout << setw(12) << x << setw(12) << xnext << setw(12) << fx << setw(12) << fxnext << endl;
	return;
}
void ParamFxmin::_solver_R1_fout_exit(std::ofstream &fout, const double &x0, const double &xmin, 
									  const double &fx0, const double &fxmin)
{
	fout << endl << "ParamFxmin::_solver_R1(): flag = " << this->process_flag 
		<< ", time_used = " << this->time_used << endl;
	fout << setw(12) << "x0" << setw(12) << "x*" << setw(12) << "f(x0)" << setw(12) << "f(x*)" << endl;
	fout << setw(12) << x0 << setw(12) << xmin << setw(12) << fx0 << setw(12) << fxmin << endl;
	return;
}	
//	======================================================================================

