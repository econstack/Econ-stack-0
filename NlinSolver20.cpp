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
#include "math_grunts.h"
#include "MathGrunts20.h"
#include "NlinSolver20.h"

//	==================================================================================================
//	constructor 
ParamNlin::ParamNlin()
{
	Nf = 0; // default; solver will exit for system of dim 0
	eps_nlin = 0.0001;
	Nfout = 0; //  # iterations between output to log (default = no output)
	Nmaxiter = 1000; // max number of iterations
	step_jacob = stepJacob; // step size in computing numerical derivatives
	step_nlin = 0.5; // step size in broyden|newton step from x to x'
	step_dpercent_flag = OFF; // option: constrain max percent change from x to x'
	step_dpercent = 0.1; // max percent change for option
	powell_flag = OFF;
	flag = 0;
}
//	Write parameters to outfile
void ParamNlin::_fout(ofstream& fout)
{
	fout << endl << setw(15) << "ParamNlin(): ";
	fout << setw(12) << "N(f)" << setw(12) << "eps(nlin)" << setw(12) << "step(nlin)" 
		<< setw(12) << "N(#iter)" << setw(12) << "N(fout)" << setw(12) << "step(jacob)" << endl;
	fout << setw(27) << Nf << setw(12) << eps_nlin << setw(12) << step_nlin 
		<< setw(12) << Nmaxiter << setw(12) << Nfout << setw(12) << step_jacob << endl;
	return;
}
/*	-----------------------------------------------------------------------------------------
	Broyden_method utilites:
	(i) get_jxnext() : update Jacobian in broyden_method: 
		Set y(k) = f(x(k+1)) - f(x(k)) and set J(k+1) = J(k) + {(y(k)-J(k)s(k))s(k)'}/{s(k)'s(k)} 
		where J(k+1) = Jxnext, f(x(k+1)) = fxnext, y(k) = y_temp 
	(ii) get_xnext() : gets next candidate x(k+1) in broyden_method, also (i) gets the 
		inverse jacobian, Jinv(k), (ii) evaluates the eqns at x, f(k), and (iii) computes 
		the search direction, dir(k)
--------------------------------------------------------------------------------------------- */
void ParamNlin::_broyden_method(double *x, void (*fct)(double *fx, double *x), 
								int (*fct_constraints)(double *x), ofstream *fout)
{
	//	init check st dimension of f > 0
	if ( this->Nf == 0 ) {
		cout << endl << "Error! Broyden_method(), dim(f) = 0. Exiting" << endl;
		log_fout << endl << "Error! Broyden_method(), dim(f) = 0. Exiting" << endl;
		this->flag = -1;
		return;
	}

	double starttime, stoptime, stoptimef, starttimef = clock();
	const int Nx = this->Nf;

	//	store old stepJacob and put current parameter as global
	double stepJacob_old = stepJacob;
	stepJacob = this->step_jacob;

    //  Initialization: Set Jacobian to identity matrix
    double *Jx;	Jx = new double [Nx*Nx];
	make_identity_matrix( Jx, Nx );

	//	make copy of init guess
	double *x0;	x0 = new double [Nx];
	copy_array( x0, x, Nx );
	//	get system evaluated at init guess f(x0) 
	double *fx0; fx0 = new double [Nx];
	fct( fx0, x0 );

	//	initialize control-flow vars in main while() loop
	int i_nlin = 0;
    double loss_x = 1; double loss_fx = 1;
    double loss = loss_x + loss_fx;

    //  Local vars
    double *fx; fx = new double [Nx];
    double *fxnext; fxnext = new double [Nx];
    double *Jxinv; Jxinv = new double [Nx*Nx];
    double *Jxnext; Jxnext = new double [Nx*Nx];
    double *xnext; xnext = new double [Nx];
    double *dir; dir = new double [Nx]; 
    double *y_temp; y_temp = new double [Nx];

    //  Main loop ===========================================================================
	while ( loss > (this->eps_nlin) )   {

		//	start iteration clock
		starttime = clock();

		/*  -------------------------------------------------------------------------------
			step 1. Compute next iterate: (i) Solve J(k)s(k) = -f(x(k)) for s(k), 
			(ii) set x(k+1) = x(k) + s(k). In (i), solve for s(k) = - J(k)^{-1}f(x(k)) 
			where J(k) = Jx, s(k) = s_step, x(k) = x, x(k+1) = xnext, f(x(k)) = fx 
		--------------------------------------------------------------------------------- */
		int xnext_flag = this->_broyden_get_xnext( Jx, Jxinv, fx, dir, x, xnext, fct );
		//	exit if unable to invert jacobian 
		if ( xnext_flag < 0) {
			this->_fout_error( log_fout, i_nlin, loss_fx, x, xnext, fx, 
				"Error! broyden_method(): singular Jacobian. Exiting." );
			this->flag = xnext_flag;
			return;
		}

		/*  -----------------------------------------------------------------------------
			step 2. Check for candidate x(k+1) in constraint set, rescale if needed 
		--------------------------------------------------------------------------------- */
        int flag_constraint = fct_constraints( xnext );
        if ( flag_constraint > 0 )  {
			//	TODO rewrite rescale_xnext as member function
            int flag_xnext = rescale_xnext(x, Jxinv, fx, xnext, fct_constraints, Nx);
            if (flag_xnext < 0) {
				this->_fout_error( log_fout, i_nlin, loss_fx, x, xnext, fx, 
					"Error! broyden_method(): Unable to find x(k+1) in constraint set. Exiting." );
				this->flag = CONSTRAINT_SET_VIOLATION_FLAG;
                break;
            }   //  end if(flag_xnext)
        }   //  end if(flag_constraint)


		/*  -----------------------------------------------------------------------------
			step 3. Update Jacobian guess: Set y(k) = f(x(k+1)) - f(x(k)) and,
			set J(k+1) = J(k) + {(y(k)-J(k)s(k))s(k)'}/{s(k)'s(k)} 
			where J(k+1) = Jxnext, f(x(k+1)) = fxnext, y(k) = y_temp 
		--------------------------------------------------------------------------------- */
		this->_broyden_get_jxnext( fx, fxnext, xnext, y_temp, dir, Jx, Jxnext, fct );

		/*  -----------------------------------------------------------------------------
			step 4. Check stopping conditions; break if converged,  write out diagnostics
			if needed, update candidate; 
		--------------------------------------------------------------------------------- */
        loss_x = norm_sup( x, xnext, Nx );
        loss_x = loss_x / ( 1.0 + norm_sup( x, Nx ) );
        loss_fx = norm_sup( fxnext, Nx );
        loss = loss_x + loss_fx;
		//	break if successful convergence
		if ( loss < (this->eps_nlin) ) {
			this->flag = i_nlin;
			this->loss_fx = loss_fx;
			break;
		}
		//	option : write to default (global) logdata file 
		stoptime = clock();
		time_used = (stoptime-starttime) / 1000.0;
		if (( (this->Nfout) != 0 ) && ( fout == 0 ) ) {
			this->_fout_iter( log_fout, i_nlin, loss_x, loss_fx, time_used, x, xnext, 
				fx, fxnext );
        }
		//	option : write to specified outfile
		if (( (this->Nfout) != 0 ) && ( fout != 0 ) ) {
			this->_fout_iter( *fout, i_nlin, loss_x, loss_fx, time_used, x, xnext, 
				fx, fxnext );
		}
		//	increment iteration counter
        ++i_nlin;
		//	break if exceeded max # iterations
		if ( i_nlin >= (this->Nmaxiter) ) {
			this->_fout_error( log_fout, i_nlin, loss_fx, x, xnext, fx, 
				"Error! Broyden_method(): Exceeded max # iterations. Exiting." );
			this->flag = EXCEEDED_MAXITER_FLAG;
			return; 
        }

        //  update candidate x and jacobian
        copy_array( x, xnext, Nx );
        copy_array( Jx, Jxnext, Nx*Nx );

    }   //  end while(loss)

    //  write out summary stats on computation ----------------------------------------------
	stoptimef = clock();
	this->time_used = ( stoptimef - starttimef ) / 1000.0;
	//	write to default logfout outfile if no outfile is specified
	if ( fout == 0 ) {
		this->_fout_exit( log_fout, i_nlin, loss_fx, x0, x, fx0, fx );
	} else {	// write to specified outfile
		this->_fout_exit( *fout, i_nlin, loss_fx, x0, x, fx0, fx );
	}

	//	restore old stepJacob
	stepJacob = stepJacob_old;
    //  delete local pointers
    delete[] fx0; delete[] fx; delete[] fxnext; 
	delete[] Jx; delete[] Jxinv; delete[] Jxnext; 
	delete[] x0; delete[] xnext; delete[] dir; delete[] y_temp; 

    return;

}
void ParamNlin::_newton_method(double *x, void (*fct)(double *fx, double *x), 
							   int (*fct_constraints)(double *x), std::ofstream *fout)
{
	//	init check st dimension of f > 0
	if ( this->Nf == 0 ) {
		cout << endl << "Error! Newton_method(), dim(f) = 0. Exiting" << endl;
		log_fout << endl << "Error! Newton_method(), dim(f) = 0. Exiting" << endl;
		this->flag = -1;
		return;
	}

	//	store old stepJacob and put current newton parameter as global
	double stepJacob_old = stepJacob;
	stepJacob = this->step_jacob;

    //  local vars used in main while() loop
	int Nx = this->Nf;
	double *step_x; step_x = new double [Nx];	// step from x(k) to x(k+1)
    double *xnext; xnext = new double [Nx];		// next guess
    double *fx;	fx = new double [Nx];			// error array at x(k)
	double *fxnext;	fxnext = new double [Nx];	// error array at x(k+1)
    double *Jx;	Jx = new double [Nx*Nx];		// Jacobian at x
    double *Jxinv; Jxinv = new double [Nx*Nx];	// Jacobian inverse at x
	const_array( xnext, Nx );
	
	//	make copy of initial guess and evaluate system f at guess
	double *x0;	x0 = new double [Nx];
	copy_array( x0, x, Nx );
	double *fx0; fx0 = new double [Nx];
	fct( fx0, x0 );

	//	time use vars
	double time_used, starttime, stoptime, stoptimef, starttimef = clock();

    //  init control-flow vars in main while() loop 
    double loss_x = 0, loss_fx = norm_sup( fx0, Nx );
	double loss = loss_x + loss_fx;
    int i_nlin = 0; // Nlin solver loop counter
	int inv_mat_flag; // flag for mat inv step

    //  main loop ===========================================================================
	while ( loss > (this->eps_nlin) )   {

		//	newton iteration diagnostics - start iteration clock 
		starttime = clock();

		/*	---------------------------------------------------------------------------------
			step 1: Get Jacobian at x, J(x(k)) = J(k) and compute its inverse,
			J^(-1)(x(k)) = Jinv(k); Compute the next candidate solution x(k+1) st 
			x(k+1) = x(k) - Jinv(x)f(x)
		------------------------------------------------------------------------------------- */
		//	get jacobian
		get_jacob( x, Jx, fx, Nx, fct, stepJacob );
		//  get inverse Jacobian at x, break if unable to invert
        inv_mat_flag = inv_mat( Jxinv, Jx, Nx );
		if ( inv_mat_flag == SINGULAR_MATRIX_FLAG )  {
			this->_fout_error( log_fout, i_nlin, loss_fx, x, xnext, fx,
				"Error! Newton_method(): singular jacobian. Exiting." );
			this->flag = SINGULAR_MATRIX_FLAG;
			return;
		}
        //  compute next guess, x(k+1)
		mat_mult( Jxinv, fx, step_x, Nx, Nx, 1 );
		mat_rescale( step_x, this->step_nlin, Nx );
		mat_add( xnext, x, step_x, -1.0, Nx );

		/*	---------------------------------------------------------------------------------
			step 2 : Check if x(k+1) is acceptable, contract step if needed. (i) check if x(k+1) 
			is in the constraint set, and (ii) (optional) check if x(k+1) is within percent 
			change limit where dpercent_step = (0, 1) is represented as a fraction 
		------------------------------------------------------------------------------------- */
		//	option : check if x(k+1) is within percent change limits, rescale if needed
		if ( (this->step_dpercent_flag == ON) ) {
			int flag_dpercent_tmp = this->_get_xnext_dpercent( xnext, x );
			if ( flag_dpercent_tmp < 0 ) {
				stoptimef = clock();
				this->time_used = (stoptimef - starttimef) / 1000.0;
				this->_fout_error( log_fout, i_nlin, loss_fx, x, xnext, fx, 
					"Error! Newton_method(): Cannot find candidate within dpercent. Exiting." );
				this->flag = -1;
				break;
			}
		}
        //  check if x(k+1) is in constraint set, rescale if needed
        int flag_constraint = fct_constraints( xnext );
        if ( flag_constraint > 0 )  {
            int flag_xnext = rescale_xnext( x, Jxinv, fx, xnext, fct_constraints, Nx );
            if ( flag_xnext > 0 ) {
				stoptimef = clock();
				this->time_used = (stoptimef - starttimef) / 1000.0;
				this->_fout_error( log_fout, i_nlin, loss_fx, x, xnext, fx,
					"Error! Newton_method(): Unable to find candidate in constraint set. Exiting." );
				this->Nf = NONCONVERGENCE_FLAG;
                break;
            }
        }

		/*  ----------------------------------------------------------------------------------
			step 3: (i) Compute control-flow vars and system f(k+1) at candidate x(k+1). 
			(ii) (option) Powell hybrid, check if the sum of residuals decreases from x(k) to
			x(k+1), ie |f(x(k+1))| < |f(x(k))|. If not, contract if possible.  
		-------------------------------------------------------------------------------------- */
		//	compute error at candiate x(k+1)
		fct( fxnext, xnext );
		loss_x = norm_sup( xnext, x, Nx );
        loss_x = loss_x / (1.0 + norm_sup( x, Nx ));
        loss_fx = norm_sup( fxnext, Nx );
		loss = loss_x + loss_fx;
		//	increment iteration counter
        ++i_nlin;


		//	option: Powell hybrid, break if unable to find candidate st sum of residuals decrease
		if ( (this->powell_flag) == ON ) {
			int powell_flag_tmp = this->_get_xnext_powell( xnext, x, fx, fxnext, fct );
			if ( powell_flag_tmp < 0) {
				stoptimef = clock();
				this->time_used = (stoptimef - starttimef) / 1000.0;
				this->flag = NEWTON_POWELL_FLAG;
				this->_fout_error( log_fout, i_nlin, loss_fx, x, xnext, fx,
					"Error! Newton_method(): Powell hybrid method, failed to find candidate x(k+1)" );
				break;
			}
		}
		
		/*	----------------------------------------------------------------------------------
			step 4: (i) write out optional iteration diagnostics. (ii) check stopping conditions,
			including exceeded max # iterations, x does not change, and (iii) update x and fx
		-------------------------------------------------------------------------------------- */
		//  newton iteration diagnostics
		stoptime = clock();
		time_used = ( stoptime - starttime ) / 1000.0;
		//	option : write to  outfile
		if ( (this->Nfout) != 0 ) {
			if ( ( i_nlin%(this->Nfout)==0 ) ) {
				if ( fout != 0 ) { //	write to specified logdata outfile
					this->_fout_iter( *fout, i_nlin, loss_x, loss_fx, time_used, x, xnext, 
						fx, fxnext ); 
				} else { //	write to default global logdata outfile
					this->_fout_iter( log_fout, i_nlin, loss_x, loss_fx, time_used, x, xnext, 
						fx, fxnext ); 
				}
			}
		}

		//	exit if exceeded max # iterations
		if ( i_nlin >= (this->Nmaxiter) )    {
			stoptimef = clock();
			this->time_used = (stoptimef - starttimef) / 1000.0;
			this->_fout_error( log_fout, i_nlin, loss_fx, x, xnext, fx, 
				"Error! Newton_method(): Exceeded max # iterations. Exiting." );
			this->flag = EXCEEDED_MAXITER_FLAG;
			return;
        }

		//	exit if x does not change
		if ( ( loss_x < ((this->eps_nlin)/100) ) && ( loss_fx > (this->eps_nlin) ) ) {
			stoptimef = clock();
			this->time_used = (stoptimef - starttimef) / 1000.0;
			this->_fout_error( log_fout, i_nlin, loss_fx, x, xnext, fx, 
				"Error! Newton_method(): x(k) stuck. Exiting. " );
			this->flag = X_STUCK_FLAG;
			break;
		}

        //  update candidate x and fx
        copy_array( x, xnext, Nx );
		copy_array( fx, fxnext, Nx );

	}   //  end while()

	//  Write out summary stats for function call =============================================
	stoptimef = clock();
	time_used = ( stoptimef - starttimef ) / 1000.0;
	if ( this->flag >= 0 )
		this->flag = i_nlin;
	this->loss_fx = loss_fx;
	this->time_used = time_used;
	//	write to console 
	cout << endl << "Exit: newton_method, size=" << Nx << ", iter=" << this->flag
		<< "/" << this->Nmaxiter << ", loss(x)=" << loss_x << ", loss(fx)=" << loss_fx;
	//	write to default logdata outfile
	if ( fout == 0 ) {
		this->_fout_exit( log_fout, i_nlin, loss_fx, x0, x, fx0, fx );
	} else { //	write to specified outfile
		this->_fout_exit( *fout, i_nlin, loss_fx, x0, x, fx0, fx );
	}

	//	free up memory
    delete[] x0; delete[] xnext; delete[] Jx; delete[] Jxinv; 
	delete[] fx0; delete[] fx; delete[] fxnext; delete[] step_x;
    return;
}
/*	-----------------------------------------------------------------------------------------
	Broyden_method utilites:
	(i) get_jxnext() : update Jacobian in broyden_method: 
		Set y(k) = f(x(k+1)) - f(x(k)) and set J(k+1) = J(k) + {(y(k)-J(k)s(k))s(k)'}/{s(k)'s(k)} 
		where J(k+1) = Jxnext, f(x(k+1)) = fxnext, y(k) = y_temp 
	(ii) get_xnext() : gets next candidate x(k+1) in broyden_method, also (i) gets the 
		inverse jacobian, Jinv(k), (ii) evaluates the eqns at x, f(k), and (iii) computes 
		the search direction, dir(k)
--------------------------------------------------------------------------------------------- */
void ParamNlin::_broyden_get_jxnext(double *fx, double *fxnext, double *xnext, double *y_temp, 
									double *dir, double *Jx, double *Jxnext, 
									void (*fct)(double *fx, double *x), ofstream *fout)
{
	//	evaluate system at next candidate, get f(x(k+1))
    fct( fxnext, xnext );
	//	compute y(k) = f(k+1) - f(k)
	mat_add( y_temp, fxnext, fx, -1, this->Nf );
	//	compute (y(k) + f(k))
	double *temp1; temp1 = new double [this->Nf];
	mat_add( temp1, y_temp, fx, 1, this->Nf );
	//	compute (y(k) + f(k))s(k)T  where s(k)T = transpose s(k)
	double *temp2; temp2 = new double [(this->Nf)*(this->Nf)];
	mat_mult( temp1, dir, temp2, this->Nf, 1, this->Nf);
	//	compute s(k)T s(k)
	double temp3;
	temp3 = mat_mult( dir, dir, this->Nf );
	//	compute J(k+1) = J(k) + (y(k) + f(k))s(k)T /s(k)T s(k)
	mat_add( Jxnext, Jx, temp2, 1.0/temp3, (this->Nf)*(this->Nf) );
	//	TODO write out diagnostics
	if ( fout != 0 ) {
//		fout_jxnext_broyden_diag(log_fout, Jx, Jxnext, xnext, fx, fxnext, y_temp, Nx_tmp);
	}
	//	free up memory
	delete[] temp1; delete[] temp2;
	return;

}
int ParamNlin::_broyden_get_xnext(double *Jx, double *Jxinv, double *fx, double *dir, double *x, 
								  double *xnext, void (*fct)(double *fx, double *x), ofstream *fout)
{
	const int Ntmp = this->Nf;
	//	make copy of Jacobian
	double *Jxcopy;	Jxcopy = new double [Ntmp*Ntmp];
    copy_array( Jxcopy, Jx, Ntmp*Ntmp );

	//	get inverse Jacobian, break if singular
    int mat_inv_flag = inv_mat( Jxinv, Jxcopy, Ntmp );
	if (mat_inv_flag < 0)
		return -1;

	//	evaluate equations at current candidate
    fct( fx, x );

	//	compute direction vector dir = s(k) = -J(k)^(-1)*f(k)
    mat_mult( Jxinv, fx, dir, Ntmp, Ntmp, 1 );
	mat_rescale( dir, -1.0*(this->step_nlin), Ntmp );

	//	solve for next candidate: x(k+1) = x(k) + s(k)
    mat_add( xnext, x, dir, 1, Ntmp);

	//	TODO write out function call diagnostics
	if ( fout != 0 ) {
//		fout_xnext_broyden_diag(log_fout, step_broyden_tmp, Jx, Jxinv, dir, x, xnext, fx, Nx_tmp);
	}
	//	free heap
	delete[] Jxcopy;
	return 0;
}
/*	-------------------------------------------------------------------------------------------
	Iteration utilies:
	(i) get_xnext_dpercent(): checks if x(k+1) is within step_dpercent% of x(k), i.e
		|x(k+1)-x(k)| < step_percent 
	(ii) get_xnext_powell(): implements powell hybrid method where x(k+1) must be st the SSR 
	decreases, ie 
		|f(k+1)| < |f(k)|
	In both cases (i) and (ii), contracts the candidate x(k+1) back toward x(k) if possible 
--------------------------------------------------------------------------------------------- */
int ParamNlin::_get_xnext_dpercent(double *xnext, double *x, ofstream *fout)
{
	double *dx; dx = new double [this->Nf];
	double *dxpercent; dxpercent = new double [this->Nf];
	mat_add( dx, xnext, x, -1.0, this->Nf );

	//	control-flow vars
	int process_flag = 0;
	int jiter = 0;
	int Nmaxiter_tmp = 10;
	const double rescale_factor = 0.5;
	double max_x_dpercent; // value of max percentage change
	int j_max_x_dpercent; // index of max percentage change

	//	check if change is within limits, exit immediately if within limits
	get_percentage_change( dxpercent, x, dx, this->Nf );
	mat_rescale( dxpercent, 0.01, this->Nf );
	get_max_fabs( dxpercent, this->Nf, max_x_dpercent, j_max_x_dpercent );
	if ( fabs(max_x_dpercent) <= (this->step_dpercent) ) {
		delete[] dx; delete[] dxpercent;
		return 0;
	}

	//	make copy of original next candidate
	double *xnext_tmp; xnext_tmp = new double [this->Nf];
	copy_array( xnext_tmp, xnext, this->Nf );


	while ( fabs(max_x_dpercent) > (this->step_dpercent) ) {

		//	check if exceeded max # iterations
		if ( jiter >= Nmaxiter_tmp ) {
			process_flag = -1;
			break;
		}
		//	rescale xnext by factor 1/2
		mat_rescale( dx, rescale_factor, this->Nf );
		mat_add( xnext_tmp, x, dx, 1.0, this->Nf );
		//	check if xnext is within range
		get_percentage_change( dxpercent, x, xnext_tmp, this->Nf );
		mat_rescale( dxpercent, 0.01, this->Nf );
		get_max_fabs( dxpercent, this->Nf, max_x_dpercent, j_max_x_dpercent );
		//	increment iteration flag
		++jiter;

	}

	//	option: write diagnostics to log
	if ( fout != 0 ) {
		//	TODO write fout() diag 
		this->_get_xnext_dpercent_fout( *fout, x, xnext, dx, dxpercent, max_x_dpercent, 
			j_max_x_dpercent );
	}
	delete[] dx; delete[] dxpercent; delete[] xnext_tmp;
	return process_flag;
}
void ParamNlin::_get_xnext_dpercent_fout(ofstream& fout, double *x, double *xnext, double *dx, 
										 double *dxpct, const double& max_dxpct, const int& j_max_dxpct )
{
	get_xnext_dpercent_fout( fout, x, xnext, dx, dxpct, this->Nf, max_dxpct, j_max_dxpct, 
		this->step_dpercent );
}	
int ParamNlin::_get_xnext_powell(double *xnext, double *x, double *fx, double *fxnext, 
								 void (*fct)(double *fx, double *x), ofstream *fout)
{
	double rescale_factor = 0.5;

	int process_flag = 0;

	//	compute |f(k)| and |f(k+1)| and compare; exit immediately if |f(k+1)| < |f(k)|
	double loss_fx = norm_L2( fx, this->Nf );
	double loss_fxnext = norm_L2( fxnext, this->Nf );
	if ( loss_fxnext < loss_fx ) {
		return process_flag;
	}
	//	control-flow vars
	int Nmaxiter_tmp = 10;
	int jiter = 1;
	double loss_fxnext_tmp = loss_fxnext;

	//	compute dx = x(k+1) - x(k)
	double *dx; dx = new double [this->Nf];
	mat_add( dx, xnext, x, -1, this->Nf );
	mat_rescale( dx, 0.5, this->Nf );

	/*	---------------------------------------------------------------------------------- 
		If residuals is increasing, rescale by factor half with each iter upto Nmaxiter=10; 
		Note that 2^(10) is approx 1000. 
	-------------------------------------------------------------------------------------- */
	double *fxnext_tmp; fxnext_tmp = new double [this->Nf];
	double *xnext_tmp; xnext_tmp = new double [this->Nf];

	while ( loss_fxnext_tmp > loss_fx ) {

		//	compute new candidate and eval system f at new candidate
		mat_add( xnext_tmp, x, dx, 1.0, this->Nf );
		fct( fxnext_tmp, xnext_tmp );
		loss_fxnext_tmp = norm_L2( fxnext_tmp, this->Nf );

		//	break if residual decreases
		if ( loss_fxnext_tmp < loss_fx ) {
            copy_array( xnext, xnext_tmp, this->Nf );
			copy_array( fxnext, fxnext_tmp, this->Nf );
			process_flag = jiter;
			break;
		}

		//	increment while() counter
		++jiter;
		//	break if exceeded max number of attempts
		if (jiter >= Nmaxiter_tmp) {
			process_flag = NEWTON_POWELL_FLAG;
			break;
		}
		
		//	rescale dx for next iteration
		mat_rescale( dx, 0.5, this->Nf );
	}

	//	option: write out diagnostics
	if ( fout != 0 ) {
		//	TODO write member function for _fout_powell()
		//fout_get_xnext_newton_powell(log_fout, x, xnext_old, xnext_tmp, Nx, loss_fx, loss_fxnext,
		//	loss_fxnext_tmp);
	}
	//	free heap
	delete[] fxnext_tmp; delete[] xnext_tmp; 
	return process_flag;
}
//	--------------------------------------------------------------------------------------------------
//	Iteration diagnostics - write to outfile 
void ParamNlin::_fout_iter(std::ofstream &fout, const int &i_nlin, const double &loss_x, const double &loss_fx, 
						   const double &time_used, double *x, double *xnext, double *fx, double *fxnext)
{
	int Nprec = 6;
	//	compute norms at f(k) and f(k+1)
	double norm_fx = norm_L2( fx, this->Nf );
	double norm_fxnext = norm_L2( fxnext, this->Nf );
	//	write header and summary stats
	cout << endl << "ParamNlin::_solver(): iter= " << i_nlin << ", loss(x)="
		<< loss_x << ", loss(fx)=" << loss_fx;
	fout << endl << "ParamNlin::_solver(): iter= " << i_nlin << ", loss(x)="
		<< loss_x << ", loss(fx)=" << loss_fx << ", time_used=" << time_used << endl;
	//	write out current candidate 
	fout << setw(16) << "x(k)"; 
	fout_array( fout, x, 1, this->Nf, Nprec, 0);
	//	write out equations resids at current candidate
	fout << setw(8) << "|f(k)|" << setw(8) << norm_fx; 
	fout_array( fout, fx, 1, this->Nf, Nprec, 0);
	//	write out next candidate
	fout << setw(16) << "x(k+1)"; 
	fout_array( fout, xnext, 1, this->Nf, Nprec, 0);
	//	write out equations resids at next candidate
	fout << setw(8) << "|f(k+1)|" << setw(8) << norm_fxnext;
	fout_array( fout, fxnext, 1, this->Nf, Nprec, 0);
	return;
}
//	Exit solver diagnostics - write to outfile
void ParamNlin::_fout_exit(std::ofstream &fout, const int &i_nlin, const double &loss_fx,  
						   double *x0, double *xstar, double *fx0, double *fxstar)
{
	int Nprec = 6;
	//	compute norms at intial guess and solution
	double norm_fx0 = norm_L2( fx0, this->Nf );
	double norm_fxstar = norm_L2( fxstar, this->Nf );
	//	write header 
	fout << endl << "Exit ParamNlin::_solver(): iter = " << i_nlin << "/" << this->Nmaxiter 
		<< ", time_used = " << this->time_used << endl;
	//	write out initial data: guess x0 and f(x0)
	fout << setw(14) << "x0=";
	fout_array( fout, x0, 1, this->Nf, Nprec, OFF );
	fout << setw(6) << "|fx|=" << setw(8) << norm_fx0;
	fout_array( fout, fx0, 1, this->Nf, Nprec, OFF );
	//	write out solution (or latest candidate) data 
	fout << setw(14) << "x*=";
	fout_array( fout, xstar, 1, this->Nf, Nprec, OFF );
	fout << setw(6) << "|fx*|=" << setw(8) << norm_fxstar;
	fout_array( fout, fxstar, 1, this->Nf, Nprec, OFF );
	return;
}





//
void ParamNlin::_fout_error(std::ofstream &fout, const int &i_nlin, const double &loss_fx, 
							double *x, double *xnext, double *fx, const std::string &label_tmp)
{
	int Nprec = 6;
	//	compute norms at f(k) and f(k+1)
	double norm_fx = norm_L2( fx, this->Nf );
	//	write header and summary stats (to console and to logfile)
	cout << endl << label_tmp << ", flag=" << this->flag;
	cout << endl << "ParamNlin::_solver(): iter= " << i_nlin << ", loss(fx)=" << loss_fx;
	fout << endl << label_tmp << ", flag=" << this->flag;
	fout << endl << "ParamNlin::_solver(): iter= " << i_nlin 
		<< ", loss(fx)=" << loss_fx << ", time_used=" << this->time_used << endl;
	//	write out current candidate 
	fout << setw(16) << "x(k)"; 
	fout_array( fout, x, 1, this->Nf, Nprec, OFF );
	//	write out next candidate
	fout << setw(16) << "x(k+1)"; 
	fout_array( fout, xnext, 1, this->Nf, Nprec, OFF );
	//	write out equations resids at current candidate
	fout << setw(8) << "|f(k)|" << setw(8) << norm_fx; 
	fout_array( fout, fx, 1, this->Nf, Nprec, OFF );
	return;
}
//	--------------------------------------------------------------------------------------------------

//	==================================================================================================