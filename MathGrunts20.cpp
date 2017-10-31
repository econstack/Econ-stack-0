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

/*	-----------------------------------------------------------------------------------------
	Checks whether x(k+1)=xnext is within lim_dxpct percent of x(k)=x. Rescales if needed. Option 
	to write out diagnostic to outfile. Used by nlin solvers and linemin solvers.
--------------------------------------------------------------------------------------------- */
int get_xnext_dpercent(double *xnext, double *x, const int& Nx, const double& lim_dxpct, 
						ofstream *fout)
{
	double *dx; dx = new double [Nx];
	double *dxpercent; dxpercent = new double [Nx];
	mat_add( dx, xnext, x, -1.0, Nx );

	//	control-flow vars
	int process_flag = 0;
	int jiter = 0;
	int Nmaxiter_tmp = 10;
	const double rescale_factor = 0.5;
	double max_x_dpercent; // value of max percentage change
	int j_max_x_dpercent; // index of max percentage change

	//	check if change is within limits, exit immediately if within limits
	get_percentage_change( dxpercent, x, dx, Nx );
	mat_rescale( dxpercent, 0.01, Nx );
	get_max_fabs( dxpercent, Nx, max_x_dpercent, j_max_x_dpercent );
	if ( fabs(max_x_dpercent) <= lim_dxpct ) {
		delete[] dx; delete[] dxpercent;
		return 0;
	}

	//	make copy of original next candidate
	double *xnext_tmp; xnext_tmp = new double [Nx];
	copy_array( xnext_tmp, xnext, Nx );


	while ( fabs(max_x_dpercent) > lim_dxpct ) {

		//	check if exceeded max # iterations
		if ( jiter >= Nmaxiter_tmp ) {
			process_flag = -1;
			break;
		}
		//	rescale xnext by factor 1/2
		mat_rescale( dx, rescale_factor, Nx );
		mat_add( xnext_tmp, x, dx, 1.0, Nx );
		//	check if xnext is within range
		get_percentage_change( dxpercent, x, xnext_tmp, Nx );
		mat_rescale( dxpercent, 0.01, Nx );
		get_max_fabs( dxpercent, Nx, max_x_dpercent, j_max_x_dpercent );
		//	increment iteration flag
		++jiter;

	}

	//	option: write diagnostics to log
	if ( fout != 0 ) {
		//	TODO write fout() diag 
		get_xnext_dpercent_fout( *fout, x, xnext, dx, dxpercent, Nx, max_x_dpercent, 
			j_max_x_dpercent, lim_dxpct );
	}
	delete[] dx; delete[] dxpercent; delete[] xnext_tmp;
	return process_flag;
}

void get_xnext_dpercent_fout(ofstream& fout, double *x, double *xnext, double *dx, double *dxpct, 
							 const int& Nx, const double& max_dxpct, const int& j_max_dxpct, 
							 const double& lim_dxpct)
{
	int Nprec = 6;
	//	write header and summary
	fout << endl << "get_xnext_dpercent(): max_dxpct = " << max_dxpct 
		<< ", max_dxpct(j)= " << j_max_dxpct << endl << ", lim_dxpct = " << lim_dxpct;
	//	write data
	fout << setw(12) << "x(k)";
	fout_array( fout, x, 1, Nx, Nprec, OFF );
	fout << setw(12) << "x(k+1)";
	fout_array( fout, xnext, 1, Nx, Nprec, OFF );
	fout << setw(12) << "dx";
	fout_array( fout, dx, 1, Nx, Nprec, OFF );
	fout << setw(12) << "dx_pct";
	fout_array( fout, dxpct, 1, Nx, Nprec, OFF );
	return;
}
