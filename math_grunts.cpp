#include <math.h>
#include <fstream>
#include <string>
#include <conio.h>
#include <iostream>
#include <iomanip>

using namespace std;

#include "math_parameters.h"
#include "io_grunts.h"
#include "io_diag.h"
#include "math_grunts.h"

//  Sets all elements of an array equal to a constant ----------------------------------
void const_array(double *x, int N, double konst)
{   for (int i=0; i<N; ++i) *(x+i) = konst; return; }
void const_array(double *x, int N)
{
	const_array(x, N, 0);
}
void const_array(int *x, int N, int konst)
{   for (int i=0; i<N; ++i) *(x+i) = konst; return; }
//  Sets a column (or row) of an array equal to a constant -----------------------------
void set_array_col(double *x, int Nr, int Nc, int j_col, double konst)
{   for (int i=0; i<Nr; ++i)  *(x + i*Nc + j_col) = konst;  return; }
void set_array_col(int *x, int Nr, int Nc, int j_col, int konst)
{   for (int i=0; i<Nr; ++i)  *(x + i*Nc + j_col) = konst;  return; }
void set_array_row(double *x, int Nr, int Nc, int j_row, double konst)
{   for (int i=0; i<Nc; ++i)  *(x + j_row*Nc + i) = konst;  return; }
void set_array_row(int *x, int Nr, int Nc, int j_row, int konst)
{   for (int i=0; i<Nc; ++i)  *(x + j_row*Nc + i) = konst;  return; }

//  Copies an array --------------------------------------------------------------------
void copy_array(double *to, double *from, int N)
{   for (int i=0; i<N; ++i) *(to+i) = *(from+i); return; }
void copy_array(int *to, int *from, int N)
{   for (int i=0; i<N; ++i) *(to+i) = *(from+i); return; }
void copy_array(double *to, const double *from, int N)
{ for (int i=0; i<N; ++i) *(to+i) = *(from+i); return; }
//  Copies an array column (or row) -----------------------------------------------------
void copy_array_col(double *to, double *from, int Nr_from, int Nc_from, int j_col_from)
{   for (int i=0; i<Nr_from; ++i) *(to+i) = *(from+i*Nc_from+j_col_from); return; }
void copy_array_col(int *to, int *from, int Nr_from, int Nc_from, int j_col_from)
{   for (int i=0; i<Nr_from; ++i) *(to+i) = *(from+i*Nc_from+j_col_from); return; }
void copy_array_row(double *to, double *from, int Nr_from, int Nc_from, int j_row_from)
{   for (int i=0; i<Nc_from; ++i) *(to+i) = *(from+j_row_from*Nc_from+i); return; }
void copy_array_row(int *to, int *from, int Nr_from, int Nc_from, int j_row_from)
{   for (int i=0; i<Nc_from; ++i) *(to+i) = *(from+j_row_from*Nc_from+i); return; }

//	copies a submatrix of size (endrow-startrow+1) by (endcol-startcol+1)
void copy_submatrix(double *to, double *from, int Nr_from, int Nc_from, int startrow,
					int startcol, int Nr_sub, int Nc_sub)
{
	//	diagnostics
	if ((Nr_sub < 1) || (Nc_sub < 1)) {
		log_fout << "Error: copy_submatrix. Nonpositive submat dimensions. Exiting...";
		return;
	}
	if (((startrow+Nr_sub)>Nr_from) || ((startcol+Nc_sub)>Nc_from)) {
		log_fout << "Error: copy_submatrix. Submat dimensions bigger than original.";
		return;
	}

	//	copy submatrix row by row
	for (int i=startrow; i<startrow+Nr_sub; ++i)
		copy_array(to+i*Nc_sub, from+i*Nc_from+startcol, Nc_sub);

	return;
}
//  Puts a column (or row) into an array -----------------------------------------------
void put_array_col(double *to, double *from, int Nr_to, int Nc_to, int j_col_to)
{   for (int i=0; i<Nr_to; ++i) *(to+i*Nc_to+j_col_to) = *(from+i); return; }
void put_array_col(int *to, int *from, int Nr_to, int Nc_to, int j_col_to)
{   for (int i=0; i<Nr_to; ++i) *(to+i*Nc_to+j_col_to) = *(from+i); return; }
void put_array_row(double *to, double *from, int Nr_to, int Nc_to, int j_row_to)
{   for (int i=0; i<Nc_to; ++i) *(to+j_row_to*Nc_to+i) = *(from+i); return; }
void put_array_row(int *to, int *from, int Nr_to, int Nc_to, int j_row_to)
{   for (int i=0; i<Nc_to; ++i) *(to+j_row_to*Nc_to+i) = *(from+i); return; }

//	scalar multiplication of a matrix --------------------------------------------------
void mat_rescale(double *m_new, double *m_old, double scalar, int N)
{
	double *m_tmp;	m_tmp = new double [N];
	copy_array(m_tmp,m_old,N);
	for (int i=0; i<N; ++i)
		*(m_new+i) = *(m_tmp+i)*scalar;
	delete[] m_tmp;
	return;
}
void mat_rescale(double *m, double scalar, int N)
{
	mat_rescale(m, m, scalar, N);
	return;
}
//	Computes the norm of an array or distance between 2 arrays in R(N) ------------------
double norm_L2(double *m1, double *m2, int N)
{
    double norm = 0;
    for (int i=0; i<N; i++)   {
        double current = *(m1+i) - (*(m2+i));
        current = current*current;
        norm += current;
    }   //  end for(i)
    norm = pow(norm, 0.5) / N;
    return norm;
}   //  end norm_L2()


double norm_L2(double *x, int N)
{
    double *zeroN;      zeroN = new double [N];
    const_array(zeroN, N, 0);
    double norm = norm_L2(x, zeroN, N);
    delete[] zeroN;
    return norm;
}
double norm_sup(double *m1, double *m2, int N)
{
    double current, norm = 0;
    for (int i=0; i<N; ++i) {
        current = *(m1+i) - (*(m2+i));
        current = fabs(current);
        if (current>norm) norm = current;
    }   //  end for(i)
    return norm;
}   //  end norm_sup()


double norm_sup(int *m1, int *m2, int N)
{
    double norm = 0; double temp;
    for (int i=0; i<N; ++i) {
        temp = *(m1+i) - (*(m2+i));
        temp = fabs(temp);
        if (temp>norm) norm = temp;
    }
    return norm;
}


double norm_sup(double *m, int N)
{
    //  create zero vector
    double *zero_tmp;      zero_tmp = new double [N];
    const_array(zero_tmp, N, 0);

    double norm = norm_sup(m, zero_tmp, N);
    delete[] zero_tmp;
    return norm;
}
//	Computes the sup norm of the percentage difference between m1 and m2 using their average as the denominator
double norm_sup_percent(double *m1, double *m2, int N)
{
	double *m_aver;		m_aver = new double [N];
	double *dm_percent;	dm_percent = new double [N];

	for (int i=0; i<N; ++i) {
		//	compute average
		*(m_aver+i) = 0.5*(*(m1+i) + *(m2+i));
		//	compute change as percentage of average
		if (*(m_aver+i)!=0)
			*(dm_percent+i) = (*(m1+i) - *(m2+i))/ *(m_aver+i);
		else *(dm_percent+i) = 0;
	}

	double norm = norm_sup(dm_percent, N);

	delete[] m_aver; delete[] dm_percent;
	return norm;
}

//	------------------- MATRIX ADDITION --------------------------------------------------------         
//	MS = ML + scalar*MR
void mat_add(double *msum, double *mleft, double *mright, double scalar_mright, int N)
{
    double *mleft_tmp, *mright_tmp;
    mleft_tmp = new double [N];  mright_tmp = new double [N];
    copy_array(mleft_tmp, mleft, N);
    copy_array(mright_tmp, mright, N);

    for (int i=0; i<N; ++i)
        *(msum+i) = *(mleft_tmp+i) + scalar_mright* *(mright_tmp+i);

    delete[] mleft_tmp; delete[] mright_tmp;
    return;
}
//	MS = scalarL*ML + scalarR*MR
void mat_add(double *sum, double *mleft, double *mright,
			 double scalar_left, double scalar_right, int N)
{
	double *mleft_tmp;		mleft_tmp = new double [N];
	mat_rescale(mleft_tmp,mleft,scalar_left,N);
	mat_add(sum,mleft_tmp,mright,scalar_right,N);
	delete[] mleft_tmp;
	return;
}
//	adds a constant (scalar) to each element
void mat_add(double *sum, double *m, double scalar, int N)
{
	double *m_temp;		m_temp = new double [N];
	const_array(m_temp, N, scalar);
	mat_add(sum, m, m_temp, 1.0, N);
	delete[] m_temp;
	return;
}
//	sums the elements of an array
double mat_sum(double *m, int N)
{	double sum = 0; for (int i=0; i<N; ++i) { sum += *(m+i); } return sum; }
int mat_sum(int *m, int N)
{ int sum = 0; for (int i=0; i<N; ++i) sum += m[i]; return sum; }
//	returns the number of "non-positive" elements of a matrix
int mat_positivity(double *m, int N, double eps_tmp)
{
	int mflag = 0;
	for (int i=0; i<N; ++i) {
		if (m[i] < eps_tmp) ++mflag;
	}
	return mflag;
}
//  matrix multiplication ----------------------------------------------------------------
void mat_mult(double *left_mat, double *right_mat, double *prod_mat,
              int Nrow, int Ncommon, int Ncol)
{
	//	make copy of original matrices
    double *mleft_tmp;  mleft_tmp = new double [Nrow*Ncommon];
    double *mright_tmp; mright_tmp = new double [Ncommon*Ncol];
    copy_array(mleft_tmp, left_mat, Nrow*Ncommon);
    copy_array(mright_tmp, right_mat, Ncommon*Ncol);
	//	compute product mat
	double *col_tmp; col_tmp = new double [Ncommon];
    for (int i_row=0; i_row<Nrow; ++i_row)  {
        for (int i_col=0; i_col<Ncol; ++i_col)  {

			copy_array_col(col_tmp, right_mat, Ncommon, Ncol, i_col);
			*(prod_mat+i_row*Ncol+i_col) = mat_mult(left_mat+i_row*Ncommon,
				col_tmp, Ncommon);

        }   //  end for(i_col)
    }   //  end for(i_row)
    delete[] mleft_tmp; delete[] mright_tmp;
    return;
}   //  end mat_mult()


double mat_mult(double *left_mat, double *right_mat, int N)
{
	double dot_prod = 0;
	for (int i=0; i<N; ++i)
		dot_prod += *(left_mat+i) * *(right_mat+i);
	return dot_prod;
}
void mat_mult(double *left, double *right, double *elementwiseprod, int N)
{	for (int i=0; i<N; ++i) { *(elementwiseprod+i) = *(left+i) * *(right+i); } return; }
void mat_mult(double *left, double *right, double *outerprod, int Nleft, int Nright)
{
	for (int i=0; i<Nleft; ++i) {
		for (int j=0; j<Nright; ++j) {
			*(outerprod+i*Nright+j) = *(left+i) * *(right+j);
		}
	}
	return;
}
void mat_div(double *div_elementwise, double *num, double *den, int N)
{
	for (int i=0; i<N; ++i) {
		if (*(den+i) != 0)  *(div_elementwise+i) = *(num+i) / *(den+i);
		else *(div_elementwise+i) = 0;
	}
	return;
}

//  Computes mean, deviations from mean and variance of an array ------------------------
double get_mean(double *x, int N)
{   double mean=0; for (int i=0; i<N; ++i) mean += *(x+i); mean /= N; return mean; }
double get_mean(double *x, int N, double MISSING_FLAG)
{
    double mean = 0; int Nmean = 0;
    for (int i=0; i<N; ++i)
        if (*(x+i)!=MISSING_FLAG) { mean += *(x+i); ++Nmean; }
    mean /= Nmean;
    return mean;
}

double get_mean(double *x, int Nx, double MISSING_FLAG, int& Nmissing)
{
    double mean = 0; int Nmean = 0;
    for (int i=0; i<Nx; ++i)
        if (*(x+i)!=MISSING_FLAG) { mean += *(x+i); ++Nmean; }
	if (Nmean > 0) mean /= Nmean;
	Nmissing = Nx - Nmean;
    return mean;
}
void get_deviation_from_mean(double *xdev, double *x, int N)
{ double mean = get_mean(x,N); for (int i=0; i<N; ++i) { *(xdev+i) = *(x+i)-mean; }	return; }
double get_variance(double *x, int N)
{
	double variance;
	double *x_dev;	x_dev = new double [N];
	double x_mean = get_mean(x,N);

	mat_add(x_dev, x, x_mean, N);
	variance = mat_mult(x_dev,x_dev,N)/(N-1.0);

	delete[] x_dev;
	return variance;
}
double get_standev(double *x, int N)
{
	double std = get_variance(x, N); 
	std = pow(std, 0.5);
	return std;
}
double get_correlation(double *m1, double *m2, int N)
{
	double correl = vector_correlation(m1, m2, N);
	return correl;
}
//	computes variance and stdev of a series with (assumed) mean zero 
double get_log_variance(double *x, int N)
{
	double var_x = mat_mult(x, x, N);
	var_x /= N-1.0;
	return var_x;
}
double get_log_standev(double *x, int N)
{
	double std_x = get_log_variance(x, N);
	std_x = pow(std_x, 0.5);
	return std_x;
}
void get_percentage_change(double *dx, double *x_bench, double *x_current, int N)
{
	double *dx_raw; dx_raw = new double [N];
	int Nprec = 6;
	//	compute the absolute change from benchmark values
	mat_add(dx_raw, x_current, x_bench, -1, N);
	//	get change as fraction of benchmark values
	mat_div(dx, dx_raw, x_bench, N);
	//	rescale change into percentage
	mat_rescale(dx, dx, 100, N);
	//	free mem
	delete[] dx_raw;
	return;
}
double get_percentage_change(double xBM, double xXP)
{
	double dx = 0;
	if (xBM != 0) 
		dx = (xXP - xBM) / xBM * 100;
	return dx;
}
/*	matrix interpolation (minterp = theta*m1 + (1-theta)*m2) --------------------------- */
void mat_interpolation(double *minterp, double *m1, double *m2, double theta_m1, int N)
{
	double *m1_tmp;		m1_tmp = new double [N];
	double *m2_tmp;		m2_tmp = new double [N];

	if ((theta_m1>1)||(theta_m1<0))
		fout_mat_interpolation_error( log_fout, m1, m2, theta_m1, N);

	mat_rescale(m1_tmp,m1,theta_m1,N);
	mat_rescale(m2_tmp,m2,1-theta_m1,N);

	mat_add(minterp,m1_tmp,m2_tmp,1,N);

	delete[] m1_tmp; delete[] m2_tmp;
	return;
}
//  Transpose matrice --------------------------------------------------------------------
void transpose_mat(double *mtrans, double *m, int Nr, int Nc) // TESTED
{   //  Nrows and Ncols = rows and columns of original matrix m
    double *m_tmp;  m_tmp = new double [Nr*Nc];
    copy_array(m_tmp,m,Nr*Nc);
    for (int i=0; i<Nr; ++i)
        for (int j=0; j<Nc; ++j)
            *(mtrans+j*Nr+i) = *(m_tmp+i*Nc+j);
    delete[] m_tmp;
    return;
}   //  end transpose_mat()


void mat_transpose(double *mtrans, double *m, int Nr, int Nc)
{ transpose_mat(mtrans, m, Nr, Nc); return; }
//  creates a matrix whose rows are basis for R(N) (creates an identity matrix) ------------
void make_basis_vectors(double *u, int N)
{   const_array(u,N*N,0); for (int i=0; i<N; ++i) *(u+i*N+i) = 1; return; }
//  creates an identity matrix of dim Nx ---------------------------------------------------
void make_identity_matrix(double *mat, int Nx)
{ const_array(mat,Nx*Nx,0); for (int i=0; i<Nx; ++i) *(mat+i*Nx+i) = 1; return; }


/*  returns the max (min) of elements of an array *x of length Nx */
double get_max(double *x, int Nx)
{
    double xcurrent;
    double xmax = *(x+0);
    for (int i=0; i<Nx; ++i) {
        xcurrent = *(x+i);
        if (xcurrent>xmax) xmax = xcurrent;
    }
    return xmax;
}
double get_min(double *x, int Nx)
{
    double xcurrent;
    double xmin = *(x+0);
    for (int i=0; i<Nx; ++i) {
        xcurrent = *(x+i);
        if (xcurrent<xmin) xmin = xcurrent;
    }
    return xmin;
}
//  given an array *mat, finds the max (non-neg) element and reports its value, max_mat, and its index, jmax_mat.
void get_max(double *mat, int Nx, double& max_mat, int& jmax_mat)
{
	double maxmat = 0;
	int jmaxmat = -1;

	for (int i=0; i<Nx; ++i) {
		if (mat[i] > maxmat) {
			maxmat = mat[i];
			jmaxmat = i;
		}
	}
	max_mat = maxmat;
	jmax_mat = jmaxmat;
	return;
}
//  given an array *mat, finds the largest absolute value and reports its absolute value and its index.
void get_max_fabs(double *mat, int Nx, double& max_mat, int& jmax_mat)
{
	double maxmat = 0;
	int jmaxmat = -1;

	for (int i=0; i<Nx; ++i) {
		if (fabs(mat[i]) > maxmat) {
			maxmat = fabs(mat[i]);
			jmaxmat = i;
		}
	}
	max_mat = maxmat;
	jmax_mat = jmaxmat;
	return;
}
int get_max(int *x, int Nx)
{
    int xcurrent;
    int xmax = *(x+0);
    for (int i=0; i<Nx; ++i) {
        xcurrent = *(x+i);
        if (xcurrent>xmax) xmax = xcurrent;
    }
    return xmax;
}
int get_min(int *x, int Nx)
{
    int xcurrent;
    int xmin = *(x+0);
    for (int i=0; i<Nx; ++i) {
        xcurrent = *(x+i);
        if (xcurrent<xmin) xmin = xcurrent;
    }
    return xmin;
}
//	numerical computation of 1st and 2nd derivative in R(1)
double get_dfdx(double& fx, double x, double (*fct)(double x), double step_df)
{
	double dfdx;
	double fx_above, fx_below, x_above, x_below;
	//	set x(+) = x+h and x(-) = x-h where h = step_df
	x_above = x + step_df;
	x_below = x - step_df;
	//	evaluate f(+) = f(x+h) and f(-) = f(x-h)
	fx_above = fct(x_above);
	fx_below = fct(x_below);
	//	evaluate f = 0.5(f(+)+f(-)) and dfdx = (f(+)-f(-))/2h
	fx = 0.5*(fx_above+fx_below);
	dfdx = (fx_above-fx_below)/(2*step_df);
	return dfdx;
}
double get_d2fx2(double& fx, double x, double (*fct)(double x), double step_df)
{
	double fx_above, fx_below, x_above, x_below;

	//	set x(+) and x(-)
	x_above = x+step_df;
	x_below = x-step_df;
    //  compute derivative (from above) and (from below) at x
	double dfdx_above = get_dfdx(fx_above, x_above, fct, step_df);
    double dfdx_below = get_dfdx(fx_below, x_below, fct, step_df);
	//	compute 2nd derivative
    double d2f = (dfdx_above - dfdx_below)/(2*step_df);
	fx = 0.5*(fx_above+fx_below);

	return d2f;
}
double get_d2f(double& fx, double &df, double x, double (*fct)(double x), double h)
{
	double fabove, ffabove, ffbelow, fbelow, xxabove, xabove, xbelow, xxbelow;

	//	set x+h, x+2h, x-h and x-2h
	xabove = x+h;
	xxabove = x+2.0*h;
	xbelow = x-h;
	xxbelow = x-2.0*h;
    //  evaluate function at 5 points
	ffabove = fct(xxabove);
	fabove = fct(xabove);
	ffbelow = fct(xxbelow);
	fbelow = fct(xbelow);
	fx = fct(x);
	//	compute derivatives
	df = -ffabove + 8*fabove - 8*fabove + ffbelow;
	df /= 12*h;
	double d2f = -ffabove + 16*fabove - 30*fx + 16*fbelow - ffbelow;
	d2f /= 12*pow(h, 2.0);

	return d2f;
}
/*	<get_jacob> computes jacobian of a function <*fct> from R(N) to R(N) such that J(x) = (f(x+h) - f(x-h))/2h and f(x) = 0.5(f(x+h)+f(x-h)). The step size h can be specified as an absolute or as a percentage of x. <get_gradient> computes the gradient of a function from R(N) to R. <get_dfdx_along_s> and <get_d2fdx2_along_s> computes the 1st and 2nd derivatives of a function along the direction s. */
void get_jacob_old(double *x, double *jacob, double *fx, int N,
               void (*fct)(double *fx, double *x), double step_jacob)
{
    double *fx_above;		fx_above = new double [N];
	double *fx_below;		fx_below = new double [N];
	double *x_above;		x_above = new double [N];
	double *x_below;		x_below = new double [N];

	//	init x(+) = x_above and x(-) = x_below
	copy_array(x_above, x, N);
	copy_array(x_below, x, N);

	//	turn off updating of initial guess for calibration
	int update_x0_flag_tmp = update_x0_FLAG;
	update_x0_FLAG = OFF;

    for (int i=0; i<N; ++i)  {
        for (int j=0; j<N; ++j) {

			//	get f(x+h)
            *(x_above+j) += step_jacob;
            fct(fx_above, x_above);
			//	get f(x-h)
			*(x_below+j) -= step_jacob;
			fct(fx_below, x_below);
			//	get J(x)
            *(jacob+i*N+j)=(*(fx_above+i)-*(fx_below+i))/(2*step_jacob);
			//	reset x(+) and x(-)
            *(x_above+j) -= step_jacob;
			*(x_below+j) += step_jacob;

        }   //  end for(j)

		//	get f(x)
		*(fx+i) = 0.5* (*(fx_above+i) + *(fx_below+i));

	}    //  end for(i)

	//	restore old update_x0_flag
	update_x0_FLAG = update_x0_flag_tmp;
	delete[] fx_above; delete[] fx_below; delete[] x_above; delete[] x_below;
	return;

}
void get_jacob(double *x, double *jacob, double *fx, int N,
			   void (*fct)(double *fx, double *x), double step_jacob)
{
	double *fabove;	fabove = new double [N];
	double *fbelow;	fbelow = new double [N];
	double *df; df = new double [N];
	double *xabove;	xabove = new double [N];
	double *xbelow;	xbelow = new double [N];
	//	init x(+) = x_above and x(-) = x_below
	copy_array(xabove, x, N);
	copy_array(xbelow, x, N);
	//	turn off updating of initial guess for calibration
	int update_x0_flag_tmp = update_x0_FLAG;
	update_x0_FLAG = OFF;
	//	compute f(x)
	fct(fx, x);
	//	compute jacobian
	for (int i=0; i<N; ++i) {

		//	get f(x+h)
        *(xabove+i) += step_jacob;
        fct(fabove, xabove);
		//	get f(x-h)
		*(xbelow+i) -= step_jacob;
		fct(fbelow, xbelow);
		//	get J(x)
		mat_add(df, fabove, fbelow, -1.0, N);
		mat_rescale(df, df, 1/(2*step_jacob), N);
		put_array_col(jacob, df, N, N, i);
		//	reset x(+) and x(-)
        *(xabove+i) -= step_jacob;
		*(xbelow+i) += step_jacob;

	}
	//	restore old updating-x0 flag
	update_x0_FLAG = update_x0_flag_tmp;
	//	free mem
	delete[] fabove; delete[] fbelow; delete[] df; delete[] xabove; delete[] xbelow;
	return;
}

void get_jacob(double *x, double *jacob, double *fx,
			   void (*fct)(double *fx, double *x),
			   int N, double percent_step_jacob)
{
	double step_jacob;
	double *fx_above;		fx_above = new double [N];
	double *fx_below;		fx_below = new double [N];
	double *x_above;		x_above = new double [N];
	double *x_below;		x_below = new double [N];

	copy_array(x_above, x, N);
	copy_array(x_below, x, N);

	//	turn off updating of initial guess for calibration
	int update_x0_flag_tmp = update_x0_FLAG;
	update_x0_FLAG = OFF;

	for (int i=0; i<N; ++i)  {
		for (int j=0; j<N; ++j) {
			step_jacob = *(x+j)*percent_step_jacob;
            *(x_above+j) += step_jacob;
			*(x_below+j) -= step_jacob;
            fct(fx_above, x_above);
			fct(fx_below, x_below);
            *(jacob+i*N+j)=(*(fx_above+i)-*(fx_below+i))/(2*step_jacob);
            *(x_above+j) -= step_jacob;
			*(x_below+j) += step_jacob;
        }   //  end for(j)
		*(fx+i) = 0.5*(*(fx_above+i) + *(fx_below+i));
	}    //  end for(i)

	//	restore old update_x0_flag
	update_x0_FLAG = update_x0_flag_tmp;
	//	free mem
	delete[] x_above; delete[] x_below; delete[] fx_above; delete[] fx_below;
	return;
}
void get_gradient(double *Jx, double *x, double& fx, int Nx, double (*fct)(double *x),
				  double step_jacob)
{
	//	turn off updating of initial guess for calibration
	int update_x0_flag_tmp = update_x0_FLAG;
	update_x0_FLAG = OFF;

    //double fx_above, fx_below;
	double *fx_above;	fx_above = new double [Nx];
	double *fx_below;	fx_below = new double [Nx];
    double *x_above;    x_above = new double [Nx];
    double *x_below;	x_below = new double [Nx];
	//	init x+h and x-h
	copy_array(x_above, x, Nx);
	copy_array(x_below, x, Nx);
    //  get gradient J(x)
    for (int i=0; i<Nx; ++i) {
        *(x_above+i) += step_jacob;
		*(x_below+i) -= step_jacob;
        *(fx_above+i) = fct(x_above);
		*(fx_below+i) = fct(x_below);
        *(Jx+i) = (fx_above[i] - fx_below[i])/(2*step_jacob);
		*(x_above+i) -= step_jacob;
		*(x_below+i) += step_jacob;
    }
	//	get f(x)
	fx = fct(x);
	//	write out fct call summmary
	if (fout_get_gradient_FLAG == ON) {
		fout_get_gradient_diag(log_fout, x, Jx, Nx, fx, fx_above, fx_below, step_jacob);
    }
	//	restore old update_x0_flag
	update_x0_FLAG = update_x0_flag_tmp;
	//	free mem
    delete[] x_above; delete[] x_below;
	delete[] fx_above; delete[] fx_below;
    return;
}
void get_gradient(double *Jx, double *x, double& fx, double (*fct)(double *x),
				  int Nx, double percent_change)
{
	//	turn off updating of initial guess for calibration
	int update_x0_flag_tmp = update_x0_FLAG;
	update_x0_FLAG = OFF;

    double fx_above, fx_below;
    double *x_above;    x_above = new double [Nx];
    double *x_below;	x_below= new double [Nx];
    double step_jacob;
	//	init copies of x+h and x-h
	copy_array(x_above, x, Nx);
	copy_array(x_below, x, Nx);
	//	check for bounds error in percent change parameter
    if (percent_change>1) {
		 log_fout << endl << "get_gradient: ERROR!!! percent_change > 100% !!!!"
            << " percent change=" << percent_change;
    }
	//	get gradient J(x)
	for (int i=0; i<Nx; ++i) {
        step_jacob = *(x+i)*percent_change;
        *(x_above+i) += step_jacob;
		*(x_below+i) -= step_jacob;
        fx_above = fct(x_above);
		fx_below = fct(x_below);
        *(Jx+i) = (fx_above - fx_below)/(2*step_jacob);
        *(x_above+i) -= step_jacob;
		*(x_below+i) += step_jacob;
    }
	//	get f(x)
	fx = fct(x);
	//	write out fct call summary
	if (fout_get_gradient_FLAG ==  ON) {
		 log_fout << endl << "Percent_step_jacob=" << percent_change;
		fout_get_gradient_diag(log_fout, x, Jx, Nx, fx, step_jacob);
	}

	//	restore old update_x0_flag
	update_x0_FLAG = update_x0_flag_tmp;

    delete[] x_above; delete[] x_below;
    return;
}

double get_dfdx_along_s(double& fx, double *x, double *s_dir, int Nx,
                        double (*fct)(double *x), double step_jacob)
{
    //  compute f(x+h) where h = x + s_dir*step_Jacob
    double *xbelow;         xbelow = new double [Nx];
    double *xabove;         xabove = new double [Nx];
    double *dir_tmp;    dir_tmp = new double [Nx];
	//	make copy of search direction
    copy_array(dir_tmp, s_dir, Nx);
    //  rescale dir into unit vector
    double rescale_dir = norm_L2(dir_tmp, Nx);
	mat_rescale(dir_tmp, dir_tmp, 1/rescale_dir, Nx);
    //  compute x+h and f(x+h)
    mat_add(xabove, x, dir_tmp, step_jacob, Nx);
    double fxabove = fct(xabove);
    //  compute x-h and f(x-h)
    mat_add(xbelow, x, dir_tmp, -step_jacob, Nx);
    double fxbelow = fct(xbelow);
    //  compute dfdx(x) = [f(x+h)-f(x-h)]/2h and f(x) = 0.5*(f(x+h)+f(x-h))
    double dfdx = (fxabove-fxbelow)/(2*step_jacob);
    fx = 0.5*(fxabove+fxbelow);
	//	write out function call diagnostics
	if ( fout_get_dfdx_along_s_FLAG ==  ON) {
		fout_get_dfdx_diag( log_fout, dir_tmp, fx, fxabove, fxbelow, x, xabove, xbelow,
			dfdx, Nx, step_jacob);
    }

    delete[] xabove; delete[] xbelow; delete[] dir_tmp;
    return dfdx;
}
//	3 and 5 point approximation without calling get_dfdx_along_s (use flag to set switch, default is 3 point approx)
double get_d2fdx2_along_s(double &fx, double *x, double *sdir, int Nx,
						  double (*fct)(double *x), double h, double& dfdx, int flag)
{
	//	vars for 3pt approx
	double *xabove;		xabove = new double [Nx];
	double *xbelow;		xbelow = new double [Nx];
	//	add vars for 5pt approx
	double *xxabove;	xxabove = new double [Nx];
	double *xxbelow;	xxbelow = new double [Nx];
	const_array(xxabove, Nx, 0);
	const_array(xxbelow, Nx, 0);

	//	make copy of search direction
    double *sdir_tmp;        sdir_tmp = new double [Nx];
    copy_array(sdir_tmp, sdir, Nx);
    //  rescale search direction into unit vector
    double rescale_sdir = norm_L2(sdir_tmp, Nx);
	mat_rescale(sdir_tmp, sdir_tmp, 1/rescale_sdir, Nx);

	//	compute f(x+h)
	mat_add(xabove, x, sdir_tmp, h, Nx);
	double fabove = fct(xabove);
	//	compute f(x-h)
    mat_add(xbelow, x, sdir_tmp, -h, Nx);
	double fbelow = fct(xbelow);
	//	compute f(x)
	fx = fct(x);

	double d2fdx2, ffabove = 0, ffbelow = 0;
	//	compute 2nd derivative along sdir (default case = 3 pt approx of order O(h2) )
	if (flag == DIFF_5PT_FLAG) {
		//	compute f(x+2h)
		mat_add(xxabove, x, sdir_tmp, 2.0*h, Nx);
		ffabove = fct(xxabove);
		//	compute f(x-2h)
		mat_add(xxbelow, x, sdir_tmp, -2.0*h, Nx);
		ffbelow = fct(xxbelow);
		//	5 pt approx for 2nd derivative
		d2fdx2 = -ffabove + 16*fabove - 30*fx + 16*fbelow -ffbelow;
		d2fdx2 /= 12*pow(h, 2.0);
		//	5 pt approx for 1st derivative
		dfdx = -ffabove + 8*fabove - 8*fbelow + ffbelow;
		dfdx /= 12*h;
		//	spit out 4th derivative for diag info for computing optimal h
		double d4f = (ffabove - 4*fabove + 6*fx - 4*fbelow + ffbelow);
		d4f /= pow(h, 4.0);
		log_fout << endl << "d4f=" << d4f << endl;
	}
	else {
		d2fdx2 = (fabove - 2*fx + fbelow)/pow(h, 2.0);
		dfdx = (fabove - fbelow)/(2*h);
	}
	//	write out function call diagnostics
	if (fout_get_d2fdx2_along_s_FLAG == ON) {
		fout_get_d2fdx2_diag(log_fout, sdir_tmp, ffabove, fabove, fx, fbelow, ffbelow, Nx,
			h, xxabove, xabove, x, xbelow, xxbelow, d2fdx2, dfdx);
    }

	//	free up memory
    delete[] sdir_tmp; delete[] xabove; delete[] xbelow;
    return d2fdx2;
}
double get_d2fdx2_along_s(double& fx, double *x, double *s_dir, int Nx,
                          double (*fct)(double *x), double step_jacob)
{
    double *xabove;		xabove = new double [Nx];
	double *xbelow;		xbelow = new double [Nx];

    double fx_above, fx_below, dfdx;

	//	make copy of search direction
    double *dir_tmp;        dir_tmp = new double [Nx];
    copy_array(dir_tmp, s_dir, Nx);
    //  rescale search direction into unit vector
    double rescale_dir = norm_L2(dir_tmp, Nx);
	mat_rescale(dir_tmp, dir_tmp, 1/rescale_dir, Nx);

    //  compute derivative (from above) at x
	mat_add(xabove, x, dir_tmp, step_jacob, Nx);
    double dfdx_above = get_dfdx_along_s(fx_above, xabove, dir_tmp, Nx,fct, step_jacob);
    //  compute derivative (from below) at x
    mat_add(xbelow, x, dir_tmp, -step_jacob, Nx);
    double dfdx_below = get_dfdx_along_s(fx_below, xbelow, dir_tmp, Nx, fct, step_jacob);
    dfdx = 0.5*(dfdx_above+dfdx_below);
    //  compute 2nd derivative and function at x
    double d2fdx2 = (dfdx_above - dfdx_below)/(2*step_jacob);
	fx = 0.5*(fx_above+fx_below);
	//	write out function call diagnostics
	if (fout_get_d2fdx2_along_s_FLAG ==  ON) {
		fout_get_d2fdx2_diag( log_fout, dir_tmp, dfdx_above, dfdx_below, xabove, xbelow,
			Nx, d2fdx2, step_jacob, dfdx);
    }

	//	free up memory
    delete[] dir_tmp; delete[] xabove; delete[] xbelow;
    return d2fdx2;
}

//  matrix inverter - Given matrix and current row (pivot), finds the next candidate pivot row. */
int get_pivot(double *mat /* ptr to matrix */, int current_row, int N_rows)
{
	int pivot_row = NO_PIVOT_FLAG;
    for (int i=current_row; i<N_rows; ++i)  {
		if (fabs(*(mat+i*N_rows+current_row)) > eps_MatInv)
        { pivot_row=i; break; }
    }   //  end for(i)
    return pivot_row;
}

/*	Swaps rows pivot row and cur_row in square matrices m and minv */
void swap_row(double *m, double *minv, int pivot_row, int cur_row, int N_rows)
{
    double *row_temp;   row_temp = new double [N_rows];
	if ( fout_swap_row_FLAG ==  ON) {
		 log_fout << endl << "swap row: row1=" << pivot_row << ", row2=" << cur_row;
	}
    //  copy m old row to temp
    copy_array(row_temp,m+cur_row*N_rows,N_rows);
    //  copy m pivot row to old row
    copy_array(m+cur_row*N_rows,m+pivot_row*N_rows,N_rows);
    //  copy temp to m old row
    copy_array(m+pivot_row*N_rows,row_temp,N_rows);
    //  copy m_inv_ptr old row to temp
    copy_array(row_temp,minv+cur_row*N_rows,N_rows);
    //  copy m_inv pivot row to old row
    copy_array(minv+cur_row*N_rows,minv+pivot_row*N_rows,N_rows);
    //  copy temp to minv old row
    copy_array(minv+pivot_row*N_rows,row_temp,N_rows);
    delete[] row_temp;
    return;
}


/*	Gauss-Jordan mat inv - sets all entries below pivot to zero	*/
void zero_below(double *m, double *minv, int pivot, int N_rows)
{
    //  make current row have leading one
    double m_pivot_entry = *(m+pivot*N_rows+pivot);

    for (int j=0; j<N_rows; ++j) {
        *(m+pivot*N_rows+j) = *(m+pivot*N_rows+j)/m_pivot_entry;
        *(minv+pivot*N_rows+j) = *(minv+pivot*N_rows+j)/m_pivot_entry;
    }   //  end for(j)

    //  over loop over rows
    for (int i=pivot+1; i<N_rows; ++i)  {

        double m_factor = *(m+i*N_rows+pivot);

        //  inner loop over columns
        for (int j=0; j<N_rows; ++j)    {
            *(m+i*N_rows+j) = *(m+i*N_rows+j) - m_factor* *(m+pivot*N_rows+j);
            *(minv+i*N_rows+j) = *(minv+i*N_rows+j) - m_factor* *(minv+pivot*N_rows+j);
        }   //  end for(j)

    }   //  end for(i)
    return;
}


/*	Guass-Jordan mat inv - sets all entries above pivot to zero	*/
void zero_above(double *m, double *minv, int pivot, int N_rows)
{
    //  loop over rows
    for (int i=pivot-1; i>=0; --i)  {
        double m_factor = *(m+i*N_rows+pivot);
        //  loop over columns
        for (int j=0; j<N_rows; ++j)    {
            *(m+i*N_rows+j) = *(m+i*N_rows+j) - *(m+pivot*N_rows+j)*m_factor;
            *(minv+i*N_rows+j) = *(minv+i*N_rows+j) - *(minv+pivot*N_rows+j)*m_factor;
        }// end for(j)
    }   //  end for(i)
    return;
}


/*	Guass-Jordan matrix inversion of square matrix m_ptr with #rows=N_rows. */
int inv_mat(double *m_inv_ptr, double *m_ptr, int N_rows)
{
    double *m_original; m_original = new double [N_rows*N_rows];
    copy_array(m_original, m_ptr, N_rows*N_rows);

    int process_flag = 0;
    //  init m_inv as iden matrice
    const_array(m_inv_ptr, N_rows*N_rows, 0);
    for (int i=0; i<N_rows; ++i)    *(m_inv_ptr+i*N_rows+i) = 1;

    //  make m_ptr upper triangular - loop over rows
    for (int cur_row=0; cur_row<N_rows; ++cur_row)  {
        int pivot_row = get_pivot(m_ptr, cur_row, N_rows);
		if (pivot_row !=  NO_PIVOT_FLAG)   {
            swap_row(m_ptr, m_inv_ptr, pivot_row, cur_row, N_rows);
            zero_below(m_ptr, m_inv_ptr, cur_row, N_rows);
        }   //  end if()
		if (pivot_row ==  NO_PIVOT_FLAG)   {
			fout_inv_mat_error( log_fout, cur_row, N_rows, m_ptr, m_original);
			process_flag = SINGULAR_MATRIX_FLAG;
			break;
        }   //  end if()

    }   //  end for(cur_row)

    //  make m_ptr an iden matrix - loop over rows
    for (int j=N_rows-1; j>=0; --j)
        zero_above(m_ptr, m_inv_ptr, j, N_rows);

    copy_array(m_ptr, m_original, N_rows*N_rows);
    delete[] m_original;
    return process_flag;
}
int inv_mat(double *m, int Nrows)
{
	double *minv;	minv = new double [Nrows*Nrows];
	int inv_mat_flag = inv_mat(minv, m, Nrows);
	copy_array(m, minv, Nrows*Nrows);
	delete[] minv;
	return inv_mat_flag;
}
//  computes a smooth percentage change of variables over Nt periods given base period values and end period values where
//  xend[i] = pow(dx[i], Nt) * xbase[i] forall i. Note Nt=total # periods including base and end.
void get_growth_rate(double *dx, double *xbase, double *xend, int Nx, int Nt)
{
    double *xchange; xchange = new double [Nx];
    for (int i=0; i<Nx; ++i) {
        xchange[i] = xend[i] / xbase[i];
        dx[i] = pow(xchange[i], 1.0/(Nt-1));
    }

    return;
}
//	Cholesky decomp of a positive definite matrix A=a(i,j) into a lower triangular matrix M=m(i,j) NEED TO BE WRITTEN ********************
int mat_decomp(double *a, double *m, int N)
{
	//	local vars
	int flag = 0;
	double row_sum_m_i;

	//	init M to be zero matrix (bc upper triangular part is zero), and we write over the rest
	const_array(m, N*N, 0);

	//	loop over rows of A
	for (int i=0; i<N; ++i) {

		//	for each row i, solve 1st for m(i,i) using a(i,i)
		row_sum_m_i = 0;
		for (int k=0; k<i; ++k) row_sum_m_i += pow(*(m+i*N+k), 2.0);
		*(m+i*N+i) = *(a+i*N+i) - row_sum_m_i;
		*(m+i*N+i) = pow(*(m+i*N+i), 0.5);

		//	given m(i,i), solve for m(j,i) for j>i

	}

	return flag;
}
/*	Computes the condition number of a matrix using sup norm. Note the matrix needs to be inverted. */
double cond_sup(double *m, int N)
{
	double *minv;		minv = new double [N*N];
	double cond_number, mnorm, minvnorm;

	int invert_flag = inv_mat(minv, m, N);

	if (invert_flag>0) {
		mnorm = norm_sup(m, N*N);
		minvnorm = norm_sup(minv, N*N);
		cond_number = mnorm*minvnorm;
	}
	else {
		cond_number = invert_flag;
	}

	delete[] minv;
	return cond_number;
}
/* <get_projection> computes the projection of vector y on vector x */
void get_projection(double *yproj, double *epshat, double *y, double *x, int N)
{   //  projection of vector y on vector x, error term is yhat. Vectors' size = N
    double XX=0; double Xy=0;
    //for (int i=0; i<N; ++i) {
    //    XX += *(x_ptr+i) * *(x_ptr+i);
    //    Xy += *(x_ptr+i) * *(y_ptr+i);
    //}   //  end for(i)
	XX = mat_mult(x,x,N);
	Xy = mat_mult(x,y,N);
    double beta_temp = Xy/XX;

    //for (int i=0; i<N; ++i) {
    //    *(yhat_ptr+i) = beta_temp*(*(x_ptr+i)) - *(y_ptr+i);
    //}   //  end for(i)
	mat_add(yproj,x,y,beta_temp,-1.0,N);
	mat_add(epshat,y,yproj,-1.0,N);
    return;
}   //  end get_projection()
/* Computes the correlation of 2 arrays in R(N) */
double vector_correlation(double *v1, double *v2, int N)
{
    double cor;
    double num, den1, den2;
    double v1_mean, v2_mean;

    v1_mean = get_mean(v1,N);
    v2_mean = get_mean(v2,N);

	double *v1_tmp;	v1_tmp = new double [N];
	double *v2_tmp;	v2_tmp = new double [N];

	mat_add(v1_tmp,v1,-v1_mean,N);		//	v1_tmp = v1 - v1_mean
	mat_add(v2_tmp,v2,-v2_mean,N);		//	v2_tmp = v2 - v2_mean
	num = mat_mult(v1_tmp,v2_tmp,N);	//	dot product = v1_tmp*v2_tmp
	den1 = mat_mult(v1_tmp,v1_tmp,N);	//	dot product = v1_tmp*v1_tmp
	den2 = mat_mult(v2_tmp,v2_tmp,N);	//	dot product = v2_tmp*v2_tmp

    cor = num/pow(den1*den2,0.5);

	delete[] v1_tmp; delete[] v2_tmp;
    return cor;
}

/*	Gauss-Legendre quadrature --------------------------------------------------------------			
Given lower and upper limits of integration (x1,x2) and given n, the number of quadrature points (degree of approximation), this function returns arrays, x[] and w[] for the zeroes and weights of the Gauss-Legendre n-point quadrature formula. Adapted from Numerical Recipes for C */
void gauleg(double x1 /*lower_bound*/, double x2 /*upper_bound*/,
			double *z_ptr /*zeroes*/, double *w_ptr /*weights*/, int N)
{
    int m,j,i;
    double z1, z, xm, pp, p3, p2, p1, eps2;

    eps2=3.0e-11;         // error parameter

    m=(N+1)/2;
    xm=0.5*(x2+x1);
    x1=0.5*(x2-x1);

    for (i=1;i<=m;i++) {
        // loop over desired roots; symmetric roots so only have to find half of them
		z=cos( pi_math*(i-0.25)/(N+0.5));
        do {
            p1=1.0;
            p2=0.0;
                for (j=1;j<=N;j++)  {
                    // loop up the recurrence relation to get the Legendre polynomial at z
                    p3=p2;
                    p2=p1;
                    p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
                    }   //  end for(j)

        // p1 is now the desired Legendre polynomial. Next, compute pp, its derivative
        pp=N*(z*p1-p2)/(z*z-1.0);
        z1=z;
        z=z1-p1/pp;
        } while (fabs(z-z1)>eps2);         // end do while

        *(z_ptr+i-1)=xm-x1*z;                   // scale root to the desired interval
        *(z_ptr+N-i)=xm+x1*z;                   // put in its symmetric counterpart
        *(w_ptr+i-1)=2.0*x1/((1.0-z*z)*pp*pp);  // compute weights
        *(w_ptr+N-i)=*(w_ptr+i-1);              // and its symmetric counterpart
    } // end for(i)

    return;

}
/*	Given zeroes and weights for Gauss-Legendre quadrature defined over interval [0,1], returns rescaled zeroes and weights for interval [x_lower, x_upper] */
void gauleg_rescale(double x_lower, double x_upper, double *z_base, double *z_rescaled,
					double *w_base, double *w_rescaled, int N)
{
	for (int i=0; i<N; ++i) {
		*(z_rescaled+i) = x_lower + *(z_base+i)*(x_upper-x_lower);
		*(w_rescaled+i) = (x_upper-x_lower)* *(w_base+i);
	}
	return;
}
/*	Rescales candidate solution x - NEED TO FIX RETURN VALUE (changed from flag_soln to constraint_flag (10dec05) */
int rescale_xnext(double *x_ptr, double *Jinv_ptr, double *fx_ptr, double *xx_ptr,
                  int (*fct_constraints)(double *x), int N)
{
    int flag_soln = 0;
    int i_rescale = 0;
    int constraint_flag = 1;

	double *xnext0; xnext0 = new double [N];

	if (fout_rescale_xnext_FLAG == ON) {
		copy_array(xnext0, xx_ptr, N);
	}

    while (constraint_flag>0)   {

        for (int i=0; i<N; i++) {
            //double step_x = dotprod(Jinv_ptr+i*N,fx_ptr,N);
			double step_x = mat_mult(Jinv_ptr+i*N, fx_ptr, N);
			step_x =  stepNewton*step_x;
            step_x = step_x/pow(10.0,i_rescale);
            *(xx_ptr+i) = *(x_ptr+i) - step_x;
        }   //  end for()

        constraint_flag = fct_constraints(xx_ptr);
        ++i_rescale;

		if (i_rescale > Nmaxiter_xnext)  {
			fout_xnext_diag(log_fout, x_ptr, xx_ptr, xnext0, N);
            break;
        }   //  end if()

    }   //  end while()

	if (fout_rescale_xnext_FLAG == ON) {
		fout_xnext_diag(log_fout, x_ptr, xx_ptr, xnext0, N);
	}
	delete[] xnext0;
    return constraint_flag;
}


/*	For use with line_min_newton. */
int rescale_xnext(double *x, double *xnext, int Nx, int (*fct_constraints)(double *x),
				  int Maxiter_tmp, double rescale_factor_tmp)
{
	int constraint_flag = fct_constraints(xnext);
	int j_rescale = 1;
	double scalar_tmp;

	//	make copy of original next candidate
	double *xnext0;	xnext0 = new double [Nx];
	copy_array(xnext0, xnext, Nx);

	double *dx;		dx = new double [Nx];
	//	compute dx = (xnext - x)
	mat_add(dx, xnext, x, -1, Nx);

	while (constraint_flag>0) {
		//	set xnext = x + dx/(2^j_rescale)
		scalar_tmp = pow(rescale_factor_tmp, j_rescale);
		mat_add(xnext, x, dx, 1/scalar_tmp, Nx);
		constraint_flag = fct_constraints(xnext);
		//	break if exceeded max number of iterations
		if (j_rescale>=Maxiter_tmp) {
			fout_xnext_error( log_fout, x, xnext, Nx);
			constraint_flag = CONSTRAINT_SET_VIOLATION_FLAG;
			break;
		}
		++j_rescale;	//	increment iteration counter
	}

	if (fout_rescale_xnext_FLAG ==  ON)
		fout_xnext_diag(log_fout, x, xnext, xnext0, Nx);

	delete[] dx; delete[] xnext0;
	return constraint_flag;
}

int rescale_xnext(double x, double& xnext, int (*fct_constraints)(double x),
				  int Nmaxiter_tmp, double rescale_factor_tmp)
{
	int jrescale = 0;
	double scalar_tmp;
	int constraint_flag = fct_constraints(xnext);
	//	make copy of original xnext
	double xnext0 = xnext;
	//	compute dx = (xnext - x)
	double dx = xnext - x;
	//	rescale xnext if needed
	while (constraint_flag > 0) {
		//	set xnext = x + dx/(rescale_factor_tmp^jrescale)
		scalar_tmp = pow(rescale_factor_tmp, jrescale);
		xnext = x + dx;
		constraint_flag = fct_constraints(xnext);
		//	break if exceeded max # iterations
		if (jrescale>=Nmaxiter_tmp) {
            log_fout << endl << "Error, rescale xnext: x=" << x << ", xnext=" << xnext;
			constraint_flag =  CONSTRAINT_SET_VIOLATION_FLAG;
			break;
		}
		++jrescale;
	}
	return constraint_flag;
}

/* No constraints function ptr for the case of newton_method or broyden_method without constraints on solution vector x. */
int no_constraints(double *x)
{
	return 0;
}
int no_constraints(double x)
{ return 0; }





//	----------------------------- INTERPOLATION --------------------------------------------------------- 
//	interpolates a function f(xstar,ystar) for equally spaced points (x1,..,xn) and (y1,...,ym)
double linear_interpolate(double *x, double *y, int Nx, int Ny, double *fxy, double xstar,
						  double ystar)
{
	//	check that fxy is of right size

	//	find jx and jy where where x(jx) < xstar < x(jx+1) and y(jy) < ystar < y(jy+1)
	int jx = get_index_linear_interpolate(xstar, x, Nx);
	int jy = get_index_linear_interpolate(ystar, y, Ny);
	//	find alphax and alphay
	if ((jx == -1) || (jy == -1)) {
		log_fout << "Interpolation Error: Out of range. Exiting....";
		return -1;
	}
	double alphax = (xstar - x[jx])/(x[jx+1]-x[jx]);
	double alphay = (ystar - y[jy])/(y[jy+1]-y[jy]);
	//	interpolate
	double fstar = (1-alphax)*(1-alphay)*fxy[jx*Ny+jy]
		+ alphax*(1-alphay)*fxy[(jx+1)*Ny+jy]
		+ alphax*alphay*fxy[(jx+1)*Ny+(jy+1)]
		+ (1-alphax)*alphay*fxy[jx*Ny+(jy+1)];
	return fstar;
}
int get_index_linear_interpolate(double xstar, double *x, int Nx)
{
	int jstar  = -1;
	for (int j=0; j<Nx; ++j) {
		if (xstar >= x[j]) jstar = j;
	}
	return jstar;
}
//	generates an array of uniform deviates between 0,1 (Minimum standard for random number generator of Park and Miller, see Numerical recipes 7.1 for details)
void genr_random_array(double *mat, int N, long int seed)
{
	int IA = 16807;
	int IM = 2147483647;
	double AM = (1.0/IM);
	int IQ = 127773;
	int IR = 2836;
	long int MASK = 123459876;
	long int *seed_ptr = &seed;

	for (int i=0; i<N; ++i) {
		mat[i] = genr_random_ran0(seed_ptr, IA, IM, AM, IQ, IR, MASK);
	}
	return;
}
double genr_random_ran0(long *idum, int IA, int IM, double AM, int IQ, int IR, long int MASK)
/* “Minimal” random number generator of Park and Miller. Returns a uniform random deviate
between 0.0 and 1.0. Set or reset idum to any integer value (except the unlikely value MASK)
to initialize the sequence; idum must not be altered between calls for successive deviates in
a sequence. */
{
	long k;
	double ans;
	*idum ^= MASK;		// XORing with MASK allows use of zero and other
	k = (*idum)/IQ;		// simple bit patterns for idum.
	*idum = IA*(*idum-k*IQ)-IR*k;	// Compute idum=(IA*idum) % IM without over
	if (*idum < 0) *idum += IM;		// flows by Schrage’s method.
	ans=AM*(*idum);		// Convert idum to a floating result.
	*idum ^= MASK;		// Unmask before return.
	return ans;
}