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
#include "fmin.h"
#include "io_diag.h"


/*	------------- GENERIC SUMMARY AND ERROR DIAGNOSTICS ------------------------------
Inverse matrix error diagnostics */
void fout_inv_mat_error(ofstream& fout, int cur_row, int Nrows,
						double *m_current, double *m_original)
{
	int Nprec = 6;
	//	write header and current row with no pivot
	fout << endl << "Matrix Inversion Failed, No pivot. Exiting..." <<
		", N=" << Nrows << ", cur_row=" << cur_row << endl;
	//	write out original matrix
	fout << "M(old)="; //fout_array(fout, m_original, Nrows, Nrows, 10);
    fout_array(fout, m_original, Nrows, Nrows, Nprec, ON, 14);

	//	write out current matrix
	fout << "M(now)="; //fout_array(fout, m_current, Nrows, Nrows, 10);
    fout_array(fout, m_current, Nrows, Nrows, Nprec, ON, 14);

	return;
}

//	generic constraint set error diagnostics
void fout_xnext_error(ofstream& fout, double *x, double *xnext, int Nx, string label_tmp)
{
	fout << endl << "ERROR: No candidate x(k+1) in constraint set. Exiting - "
		<< label_tmp << endl;
	//	write out current candidate
	fout << " x(k) ="; fout_array(fout, x, 1, Nx, 4, 0);
	//	write out next candidate
	fout << "x(k+1)="; fout_array(fout, xnext, 1, Nx, 4, 0);
	return;
}
void fout_xnext_error(ofstream& fout, double *x, double *xnext, int Nx)
{
	fout_xnext_error(fout, x, xnext, Nx, " ");
	return;
}
//	generic exceeded max # iteration diagnostics
void fout_max_iter_error(ofstream& fout, int Nmaxiter_tmp, double *xnext, double *x0,
						 int Nx_tmp, string label_tmp)
{
	fout << endl << label_tmp << ": Exceeded max # iterations, " << Nmaxiter_tmp << endl;
	//	write out initial candidate
	fout << " x0="; fout_array(fout, x0, 1, Nx_tmp, 4, 0);
	//	write out current candidate
	fout << "  x="; fout_array(fout, xnext, 1, Nx_tmp, 4, 0);

	return;
}
//	generic iterative solvers (nlin or fmin) function call summary if successful
void fout_iter_final_diag(ofstream& fout, double *x0_tmp, double *xstar_tmp, int Nx_tmp,
						  int jiter, int Nmaxiter_tmp, double time_used_tmp, double eps_tmp,
						  string label_tmp)
{
	fout << endl << "Exit " << label_tmp << ", iter " << jiter << "/" << Nmaxiter_tmp
		<< ", eps=" << eps_tmp << ", time_used=" << time_used_tmp << endl;
	//	write out initial guess
	fout << " x0="; fout_array(fout, x0_tmp, 1, Nx_tmp, 6, 0);
	//	write out solution
	fout << " x*="; fout_array(fout, xstar_tmp, 1, Nx_tmp, 6, 0);

	return;
}
void fout_iter_final_diag(ofstream& fout, double *x0_tmp, double *xstar_tmp, int Nx_tmp,
						  int jiter, int Nmaxiter_tmp, double time_used_tmp, string label_tmp)
{
	fout << endl << "Exit " << label_tmp << ", iter " << jiter << "/" << Nmaxiter_tmp
		<< ", time_used=" << time_used_tmp << endl;
	//	write out initial guess
	fout << " x0="; fout_array(fout, x0_tmp, 1, Nx_tmp, 6, 0);
	//	write out solution
	fout << " x*="; fout_array(fout, xstar_tmp, 1, Nx_tmp, 6, 0);


	return;
}
//	generic iterative solvers (nlin or fmin) function exit if unsuccessful
void fout_iter_exit(ofstream& fout, double *x0_tmp, double *x_tmp, int Nx_tmp,
					int jiter, int Nmaxiter_tmp, string label_tmp)
{
	fout << endl << "Exit " << label_tmp << ", iter " << jiter << "/"
		<< Nmaxiter_tmp << endl;
	//	write out initial guess
	fout << " x0="; fout_array(fout, x0_tmp, 1, Nx_tmp, 6, 0);
	//	write out current candidate
	fout << "  x="; fout_array(fout, x_tmp, 1, Nx_tmp, 6, 0);
	return;
}

//	generic nlin iteration statistics
void fout_nlin_iter(ofstream& fout, int Nx, int i_Nlin, double loss_x, double loss_fx,
					double time_used, double *x, double *xnext, double *fx, double *fxnext,
					string label_tmp)
{
	cout << endl << label_tmp << ", iter=" << i_Nlin << ", loss(x)="
		<< loss_x << ", loss(fx)=" << loss_fx;
	fout << endl << label_tmp << ", iter=" << i_Nlin << ", loss(x)="
		<< loss_x << ", loss(fx)=" << loss_fx << ", time_used=" << time_used << endl;
	//	write out current candidate
	fout << "   x(k)="; fout_array(fout, x, 1, Nx, 6, 0);
	//	write out equations resids at current candidate
	fout << "  fx(k)="; fout_array(fout, fx, 1, Nx, 6, 0);
	//	write out next candidate
	fout << " x(k+1)="; fout_array(fout, xnext, 1, Nx, 6, 0);
	//	write out equations resids at next candidate
	fout << "fx(k+1)="; fout_array(fout, fxnext, 1, Nx, 6, 0);

	return;	return;
}
//	generic nlin_solver function call stats summary
void fout_nlin_final(ofstream& fout, int Nx, int process_flag, double loss_x,
					 double loss_fx, double *x0, double *x, double *fx, double time_used)
{
	fout << endl << "Exit Newton's method: N=" << Nx << ", iter=" << process_flag
		<< ", loss(x)=" << loss_x << ", loss(fx)=" << loss_fx
		<< ", time_used=" << time_used << endl;
	fout << "x0=";
	fout_array(fout,x0,1,Nx,6,0);
	fout << " x=";
	fout_array(fout,x,1,Nx,6,0);
	fout << "fx=";
	fout_array(fout,fx,1,Nx,6,0);
	return;
}
/*	------------- LINEMIN NEWTON DIAGNOSTICS AND ERRORS	---------------------------			linemin_newton iteration diagnostics (note final exit summary is not fct specific). */
void fout_linemin_newton_iter(ofstream& fout, int jiter, double timeused,
							  double dfdx, double d2fdx, double dfnext,
							  double d2fnext, double *dir, double *xnext,
							  double *x, double fx, double fxnext,
							  int Nx, double loss_x)
{
	//	write header
	fout << endl << "Linemin_newton iter=" << jiter << ", time used (sec)="
		<< setiosflags(ios::fixed) << setprecision(6) << timeused << ", loss(x)=" << loss_x;
	//	write 1st and 2nd derivative and fx
	fout << endl << "df=" << dfdx << ", d2f=" << d2fdx << ", fx=" << fx
		<< ", dfnext=" << dfnext << ", d2fnext=" << d2fnext << ", fnext=" << fxnext << endl;
	//	write search direction, and current and next candidate
	fout << "   dir=";
	fout_array(fout, dir, 1, Nx, 6, 0);
	fout << "x(k+1)=";
	fout_array(fout, xnext, 1, Nx, 6, 0);
	fout << "  x(k)=";
	fout_array(fout, x, 1, Nx, 6, 0);
	return;
}
//	*** OLD-OBSOLETE VERSION ***
void fout_linemin_newton_iter_diag(ofstream& fout, int jiter, double timeused,
								   double dfdx, double d2fdx, double *dir,
								   double *xnext, double *x, double fx, int Nx,
								   double loss_x)
{
	double dfnext = -1;
	double d2fnext = -1;
	double fxnext = -1;
	fout_linemin_newton_iter(fout, jiter, timeused, dfdx, d2fdx, dfnext,
		d2fnext, dir, xnext, x, fx, fxnext, Nx, loss_x);
	return;
}
/*	-------------- LINEMIN GOLDEN DIAGNOSTICS AND ERRORS ----------------------------
Cannot find initial triple for golden section. */
void fout_linemin_golden_findtriple_error(ofstream& fout, double *x0, double *dir, int Nx)
{
	fout << endl << "Exiting linemin. Cannot find initial triple." << endl;
	fout << setw(10) << "x0=";
	fout_array(fout,x0,1,Nx,4,0);
	fout << setw(10) << "dir=";
	fout_array(fout,dir,1,Nx,4,0);
	return;
}
/*	iteration and function call summary stats. */
void fout_linemin_golden_iter(ofstream& fout, int jiter, double timeused,
							  double fa, double *xa, double fc, double *xc, int N)
{
	fout << endl << "Linemin iter diag: " << jiter << ", time used (sec)="
		<< setiosflags(ios::fixed) << setprecision(4) << timeused ;
	fout << endl << "f(x_lower)=" << fa << ", x_lower=";
	fout_array(fout,xa,1,N,4,0);
	fout << "f(x_upper)=" << fc << ", x_upper=";
	fout_array(fout,xc,1,N,4,0);
	return;
}
void fout_linemin_golden_final(ofstream& fout, int process_flag, double timeused,
							   double fa, double *xa, double fc, double *xc, int N)
{
	fout << endl << "Exiting Linemin, iter=" << process_flag << ", time used (sec)="
		<< setiosflags(ios::fixed) << setprecision(5) << timeused;
	fout << endl << "f(x_lower)=" << fa << ", x_lower=";
	fout_array(fout,xa,1,N,4,0);
	fout << "f(x_upper)=" << fc << ", x_upper=";
	fout_array(fout,xc,1,N,4,0);

	return;
}
/*	------------- GET DERIVATIVE AND GRADIENT DIAGNOSTICS ------------------------- */
void fout_get_dfdx_diag(ofstream& fout, double *dir_tmp, double fx, double fxabove,
						double fxbelow,	double *x, double *xabove, double *xbelow,
						double dfdx, int Nx, double step_jacob)
{
	int Nprecision_tmp = 6;
	//	write out header, df/dx and h=step_jacob
    fout << endl << "get_dfdx_along_s, dfdx=" << dfdx << ", step_jacob=" << step_jacob << endl;
	//	write out direction
    fout << "   dir=";
    fout_array(fout, dir_tmp, 1, Nx, Nprecision_tmp, 0);
	//	write out f(x) and x
    fout << "    fx=" << setw(10) << fx << ",    x= ";
    fout_array(fout, x, 1, Nx, Nprecision_tmp, 0);
	//	write out f(x+h) and x+h
    fout << "f(x+h)=" << setw(10) << fxabove << ", (x+h)=";
    fout_array(fout, xabove, 1, Nx, Nprecision_tmp, 0);
	//	write out f(x-h) and x-h
    fout << "f(x-h)=" << setw(10) << fxbelow << ", (x-h)=";
    fout_array(fout, xbelow, 1, Nx, Nprecision_tmp, 0);
	return;
}

void fout_get_d2fdx2_diag(ofstream& fout, double *sdir_tmp, double ffabove, double fabove,
						  double fx, double fbelow, double ffbelow, int Nx, double h,
						  double *xxabove, double *xabove, double *x,
						  double *xbelow, double *xxbelow, double d2fdx2, double dfdx)
{
	int Nprec = 6;
	fout << endl << "get_d2fdx2_diag: d2fdx2 = " << d2fdx2 << ", dfdx=" << dfdx
		<< ", h=" << h << endl;
	//	write out search dir
	fout << setw(26) << "dir";
	fout_array(fout, sdir_tmp, 1, Nx, Nprec, 0);
	//	write out function evaluations and correponding points of eval
	fout << setw(7) << "f(x+2h)" << setw(12) << ffabove << setw(7) << "x+2h";
	fout_array(fout, xxabove, 1, Nx, Nprec, 0);
	fout << setw(7) << "f(x+h)" << setw(12) << fabove << setw(7) << "x+h";
	fout_array(fout, xabove, 1, Nx, Nprec, 0);
	fout << setw(7) << "f(x)" << setw(12) << fx << setw(7) << "x";
	fout_array(fout, x, 1, Nx, Nprec, 0);
	fout << setw(7) << "f(x-h)" << setw(12) << fbelow << setw(7) << "x-h";
	fout_array(fout, xbelow, 1, Nx, Nprec, 0);
	fout << setw(7) << "f(x-2h)" << setw(12) << ffbelow << setw(7) << "x-2h";
	fout_array(fout, xxbelow, 1, Nx, Nprec, 0);

	return;
}
void fout_get_d2fdx2_diag(ofstream& fout, double *dir_tmp, double dfdx_above, double dfdx_below,
						  double *xabove, double *xbelow, int Nx, double d2fdx2,
						  double step_jacob, double dfdx)

{
	int Nprecision_tmp = 6;
	fout << endl << "get_d2fdx2_along_s diag, d2f/dx2=" << d2fdx2 << ", dfdx=" << dfdx <<
		", step_jacob=" << step_jacob << endl;
	//	write out search direction
    fout << "dir=";
    fout_array(fout, dir_tmp, 1, Nx, Nprecision_tmp, 0);
	//	write out derivative from above and evaluations point
    fout << "dfdx(+)=" << setw(10) << dfdx_above << ", x(+)=";
    fout_array(fout, xabove, 1, Nx, Nprecision_tmp, 0);
	//	write out derivative from below and evaluations point
    fout << "dfdx(-)=" << setw(10) << dfdx_below << ", x(+)=";
	fout_array(fout, xbelow, 1, Nx, Nprecision_tmp, 0);
	return;
}
void fout_get_gradient_diag(ofstream& fout, double *x, double *Jx, int Nx,
							double fx, double step_jacob)
{
	int Nprecision_tmp = 6;
	//	write header, and fx, step_jacob
	fout << setiosflags(ios::fixed) << setprecision(6);
	fout << endl << "Get gradient summary: " << " fx=" << setw(12) << fx <<
		setw(12) << ", stepjacob=" << setw(12) << step_jacob << endl;
	//	write out x and J(x)
	fout << setw(10) << "x";
	fout_array(fout, x, 1, Nx, Nprecision_tmp, 0);
	fout << setw(10) << "Jx";
	fout_array(fout, Jx, 1, Nx, Nprecision_tmp, 0);
	return;
}

void fout_get_gradient_diag(ofstream& fout, double *x, double *Jx, int Nx, double fx,
							double *fx_above, double *fx_below, double step_jacob)
{
	int Nprec = 6;
	fout_get_gradient_diag(fout, x, Jx, Nx, fx, step_jacob);
	fout << setw(10) << "fx(+)";
	fout_array(fout, fx_above, 1, Nx, Nprec, 0);
	fout << setw(10) << "fx(-)";
	fout_array(fout, fx_below, 1, Nx, Nprec, 0);
	return;
}
/*	---------------- BROYDEN METHOD DIAGNOSTICS AND ERRORS -------------------------------- get_xnext_broyden general diagnostics */
void fout_xnext_broyden_diag(ofstream& fout, double step_broyden_tmp, double *Jx, double *Jxinv,
							 double *dir, double *x, double *xnext, double *fx, int Nx_tmp)
{
	fout << endl << "get_xnext_broyden_diag: step_broydent=" << step_broyden_tmp << endl;
	//	write out current candidate
	fout << "   x(k)="; fout_array(fout, x, 1, Nx_tmp, 6, 0);
	//	write out search direction
	fout << "   s(k)=";	fout_array(fout, dir, 1, Nx_tmp, 6, 0);
	//	write out next candidate
	fout << " x(k+1)=";	fout_array(fout, xnext, 1, Nx_tmp, 6, 0);
	//	write out equations evaluated at current candidate
	fout << "   f(k)="; fout_array(fout, fx, 1, Nx_tmp, 6, 0);
	//	write out jacobian at current candidate
	fout << "   J(k)=";	fout_array(fout, Jx, Nx_tmp, Nx_tmp, 6, 1);
	//	write out jacobian inverse at current candidate
	fout << "Jinv(k)="; fout_array(fout, Jxinv, Nx_tmp, Nx_tmp, 6, 1);
	return;
}
//	get_jxnext_broyden general diagnostics
void fout_jxnext_broyden_diag(ofstream& fout, double *Jx, double *Jxnext,
							  double *xnext, double *fx, double *fxnext,
							  double *y_temp, int Nx_tmp)
{
	fout << endl << "get_jxnext_broyden_diag: " << endl;
	//	write out function at current candidate f(k)
	fout << "   f(k)=";	fout_array(fout, fx, 1, Nx_tmp, 6, 0);
	//	write out jacobian at current candidate j(k)
	fout << "   J(k)="; fout_array(fout, Jx, Nx_tmp, Nx_tmp, 6);
	//	write out y(k)
	fout << "   y(k)="; fout_array(fout, y_temp, 1, Nx_tmp, 6, 0);
	//	write out next candidate
	fout << " x(k+1)=";	fout_array(fout, xnext, 1, Nx_tmp, 6, 0);
	//	write out function at next candidate f(k+1)
	fout << " f(k+1)=";	fout_array(fout, fxnext, 1, Nx_tmp, 6, 0);
	//	write out jacobian at next candidate j(k+1)
	fout << " J(k+1)="; fout_array(fout, Jxnext, Nx_tmp, Nx_tmp, 6);
	return;
}
/*	------------------ FMIN IN R(1) DIAG AND ERRORS -------------------------------------
	iteration stats for fmin_R1 */
void fout_fminR1_iter(ofstream& fout, int jfmin, int Nmax_iter_tmp, double x,
					  double xnext, double loss_x, double loss_fx, double loss_df,
					  double df, double d2f, double time_used)
{
	fout << endl << "fmin_R1, iteration " << jfmin << "/" << Nmax_iter_tmp
		<< ", x=" << x << ", xnext=" << xnext << endl;
	fout << "loss(x)=" << loss_x << ", loss(fx)=" << loss_fx << ", loss(df)="
		<< loss_df << ", df=" << df << ", d2f=" << d2f << ", time_used=" << time_used;
	return;
}
//	fct call summary
void fout_fminR1_final(ofstream& fout, int process_flag, int Nmax_iter_tmp,
					   double loss_x, double loss_fx, double loss_df,
					   double time_used, double x0_tmp, double xmin)
{
	fout << endl << "Exit fminR1, process_flag " << process_flag << "/"
		<< Nmax_iter_tmp << ", loss(x)=" << loss_x << ", loss(fx)=" << loss_fx
		<< ", loss(df)=" << loss_df << ", time_used=" << time_used;
	fout << endl << "x=" << x0_tmp << ", x*=" << xmin;
	return;
}
/*	------------------- FMIN_BFGS DIAG AND ERRORS ----------------------------------------- 
	iteration statistics for fminBFGS */
void fout_fminBFGS_iter(ofstream& fout, int jfmin, int Nmaxiter_tmp, double timeused,
						double fx, double fxnext, double *x, double *xnext, int N,
						double norm_Jxnext, double norm_dx, double eps_fmin)
{
	int Nprecision_tmp = 6;
	fout << setiosflags(ios::fixed) << setprecision(Nprecision_tmp);
	fout << endl << "fmin_BFGS iteration: " << jfmin << "/" << Nmaxiter_tmp
		<< ", N=" << N
		<< ", time used (seconds) =" << timeused 
		<< ", |Jx(k+1)|=" << norm_Jxnext 
		<< ", |x(k+1)-x(k)|=" << norm_dx 
		<< ", eps(fmin)=" << eps_fmin << endl;
	//	write out current candidate
	fout << setw(12) << "f(x(k))" << setw(12) << fx << setw(12) << "x(k)=";
	fout_array(fout, x, 1, N, Nprecision_tmp, 0);
	//	write out next candidate
	fout << setw(12) << "f(x(k+1))" << setw(12) << fxnext << setw(12) << "x(k+1)=";
	fout_array(fout, xnext, 1, N, Nprecision_tmp, 0);
	return;
}
void fout_fminBFGS_iter(ofstream& fout, int jfmin, int Nmaxiter_tmp, double timeused,
						double fx, double fxnext, double *x, double *xnext, int N,
						double norm_Jxnext, double norm_dx)
{
	int Nprecision_tmp = 6;
	fout << endl << "fmin_BFGS iteration: " << jfmin << "/" << Nmaxiter_tmp
		<< ", time used (seconds) =" << setiosflags(ios::fixed) << setprecision(6)
		<< timeused << ", |Jx(k+1)|=" << norm_Jxnext << ", |x(k+1)-x(k)|="
		<< norm_dx << endl;
	//	write out current candidate
	fout << setw(12) << "f(x(k))" << setw(12) << fx << setw(12) << "x(k)=";
	fout_array(fout, x, 1, N, Nprecision_tmp, 0);
	//	write out next candidate
	fout << setw(12) << "f(x(k+1))" << setw(12) << fxnext << setw(12) << "x(k+1)=";
	fout_array(fout, xnext, 1, N, Nprecision_tmp, 0);
	return;
}
void fout_fminBFGS_iter(ofstream& fout, int jfmin, int Nmaxiter_tmp, double timeused,
						double fx, double fxnext, double *x, double *xnext, int N)
{
	fout_fminBFGS_iter(fout, jfmin, Nmaxiter_tmp, timeused, fx, fxnext, x, xnext, N,
		-1, -1);
	return;
}
void fout_fminBFGS_final(ofstream& fout, double *x0, double fx0, double *xmin, double fxmin,
						 int Nx, int process_flag, int Nmax_iter_fmin, double time_used,
						 double eps_fmin)
{
	fout_iter_final_diag(fout, x0, xmin, Nx, process_flag, Nmax_iter_fmin,
		time_used, eps_fmin, "fmin_BFGS");
	fout << "f(x0)= " << fx0 << ", f(x*)= " << fxmin << endl;
	return;
}
//	general debugging diagnostics for fminBFGS
void fout_fminBFGS_diag(ofstream& fout, double *x, double *xnext, double *y, double *z,
						double *Jxnext, double *Hxnext, int N, int jfmin, int Nmaxiter_tmp,
						double loss_x, double loss_Jx)
{
	int Nprec = 6;
	fout << endl << "fminBFGS diag: iter=" << jfmin << "/" << Nmaxiter_tmp
		<< ", loss_x=" << loss_x << ", loss_Jx=" << loss_Jx << endl;
	fout << "   x(k)=";
	fout_array(fout, x, 1, N, Nprec, 0);
	fout << " x(k+1)=";
	fout_array(fout, xnext, 1, N, Nprec, 0);
	fout << "   y(k)=";
	fout_array(fout, y, 1, N, Nprec, 0);
	fout << "   z(k)=";
	fout_array(fout, z, 1, N, Nprec, 0);
	fout << " J(k+1)=";
	fout_array(fout, Jxnext, 1, N, Nprec, 0);
	fout << " H(k+1)=";
	fout_array(fout, Hxnext, N, N, Nprec, 1);
	return;
}
//	suboptimal point error (non-zero jacobian at solution)
void fout_fminBFGS_badx_error(ofstream& fout, double *x, double *xnext, double *Jxnext, int Nx,
							  int jfmin, int Nmaxiter_tmp, double norm_x, double norm_Jx)
{
	fout << endl << "Exit: fmin_x_suboptimal, iter=" << jfmin << "/" << Nmaxiter_tmp
		<< ", norm(x)=" << norm_x << ", norm(Jx)=" << norm_Jx << endl;
	fout << "   x(k)=";
	fout_array(fout, x, 1, Nx, 4, 0);
	fout << " x(k+1)=";
	fout_array(fout, xnext, 1, Nx, 4, 0);
	fout << "Jx(k+1)=";
	fout_array(fout, Jxnext, 1, Nx, 4, 0);
	return;
}
//	writes out jacobian and hessian (read in from file) on resumption of fminBFGS
void fout_fminBFGS_resume(ofstream& fout, double *Jx, double *Hx, double *x,
						  double *x_linemin, int Nx)
{
	fout << endl << "fminBFGS_resume" << endl;
	fout << "Jx=";
	fout_array(fout, Jx, Nx, 1, 4);
	fout << "Hx=";
	fout_array(fout, Hx, Nx, Nx, 4);
	fout << " x=";
	fout_array(fout, x, Nx, 1, 4);
	fout << "x(line)=";
	fout_array(fout, x_linemin, Nx, 1, 4);
	return;
}

//	----------------------- POLYTOPE MINIMIZATION DIAG ------------------------------------------ iteration statistics
void fout_polytope_iter(ofstream& fout, double timeused, double *vertices,
						double *fvertices, int Nx)
{
	int Nprecision_tmp = 6;
	fout << endl << "time used (seconds) =" << setiosflags(ios::fixed) <<
			setprecision(Nprecision_tmp) << timeused << endl;
	fout << "vertices=";
	fout_array(fout, vertices, Nx+1, Nx, Nprecision_tmp, 1);
	fout << "fvertices=";
	fout_array(fout, fvertices, 1, Nx+1, Nprecision_tmp, 1);
	return;
}
//	fct call summary
void fout_polytope_final(ofstream& fout, double timeused, double *xmax,
						 double fxmax, int iter_ctr, int Nx)
{
	fout << endl << "Exiting polytope_fmax. Iterations="
		<< setiosflags(ios::fixed) << setprecision(0) << iter_ctr
		<< setiosflags(ios::fixed) << setprecision(6)
		<< " Time used (seconds) =" << timeused << endl;
	fout << "f(xmax)=" << fxmax << ", xmax=";
	fout_array(fout, xmax, Nx, 1, 6, 1);
	return;
}

/*	------------------ MISC -----------------------------------------------	*/
//	error in interpolation of 2 arrays st m(theta) = theta*m1 + (1-theta)*m2
void fout_mat_interpolation_error(ofstream& fout, double *m1, double *m2,
								  double theta_m1, int N)
{
	fout << endl << "Interpolation error, theta=" << theta_m1 << endl;
	//	write out first array
	fout << "m1="; fout_array(fout, m1, N, 1, 4);
	//	write out second array
	fout << "m2="; fout_array(fout, m2, N, 1, 4);
	return;
}
//	constraint set call <rescale_xnext> summary
void fout_xnext_diag(ofstream& fout, double *x, double *xnext, double *xnext0, int Nx)
{
	fout << endl << "Exiting rescale_xnext" << endl;
	//	current candidate
	fout << "   x(k)="; fout_array(fout, x, 1, Nx, 4, 0);
	//	original next candidate
	fout << "x0(k+1)="; fout_array(fout, xnext0, 1, Nx, 4, 0);
	//	rescaled next candidate
	fout << " x(k+1)="; fout_array(fout, xnext, 1, Nx, 4, 0);
	return;
}
void fout_vf_final(ofstream& fout, int Niters, double loss, double time_used)
{
	fout << endl << "VF iteration: No. of iters=" << Niters << ", loss=" << loss
		<< ", time_used=" << time_used;
	return;
}
/*	Writes out the diagnostics for get_eqm() */
void fout_get_eqm_diag(ofstream& fout, double *x0, double *x, int Nx,
					   void (*fct)(double *fx, double *x))
{
	fout << endl << "get_eqm_diag" << endl;

	fout << "x0= ";
	fout_array(fout, x0, 1, Nx, 6, 0);

	fout << "x=  ";
	fout_array(fout, x, 1, Nx, 6, 0);

	double *fx; fx = new double [Nx];
	fct(fx, x);
	fout << "fx= ";
	fout_array(fout, fx, 1, Nx, 6, 0);

	delete[] fx;
	return;
 }

/*	--------------------- ECON NLIN ----------------------------------------                        diagnostics for eqm2eqm recursion call */
void fout_eqm2eqm_recursion_diag(ofstream& fout, double *param_BM_tmp, double *x0_BM_tmp,
								 double *param_XP_tmp, int Nparam_tmp, int Nx_tmp)
{
	fout << endl << "eqm2eqm_solver diagnostics" << endl;
	fout << setw(12) << "params(BM)=";
	fout_array(fout, param_BM_tmp, 1, Nparam_tmp, 6, 0);
	fout << setw(12) << "params(XP)=";
	fout_array(fout, param_XP_tmp, 1, Nparam_tmp, 6, 0);
	fout << setw(12) << "x(BM)=";
	fout_array(fout, x0_BM_tmp, 1, Nx_tmp, 6, 0);
	return;
}
