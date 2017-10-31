//  Sets all elements of an array equal to a constant ----------------------------------
void const_array(double *x, int N, double konst);
void const_array(int *x, int N, int konst = 0);
void const_array(double *x, int N);
//  Sets a column (or row) of an array equal to a constant -----------------------------
void set_array_col(double *x, int Nr, int Nc, int j_col, double konst);
void set_array_col(int *x, int Nr, int Nc, int j_col, int konst);
void set_array_row(double *x, int Nr, int Nc, int j_row, double konst);
void set_array_row(int *x, int Nr, int Nc, int j_row, int konst);
//  Copies an array --------------------------------------------------------------------
void copy_array(double *to, double *from, int N);
void copy_array(int *to, int *from, int N);
void copy_array(double *to, const double *from, int N);
//  Copies an array column (or row) -----------------------------------------------------
void copy_array_col(double *to, double *from, int Nr_from, int Nc_from, int j_col_from);
void copy_array_col(int *to, int *from, int Nr_from, int Nc_from, int j_col_from);
void copy_array_row(double *to, double *from, int Nr_from, int Nc_from, int j_row_from);
void copy_array_row(int *to, int *from, int Nr_from, int Nc_from, int j_row_from);
void copy_submatrix(double *to, double *from, int Nr_from, int Nc_from, int startrow,
					int startcol, int endrow, int endcol);
//  Puts a column (or row) into an array
void put_array_col(double *to, double *from, int Nr_to, int Nc_to, int j_col_to);
void put_array_col(int *to, int *from, int Nr_to, int Nc_to, int j_col_to);
void put_array_row(double *to, double *from, int Nr_to, int Nc_to, int j_row_to);
void put_array_row(int *to, int *from, int Nr_to, int Nc_to, int j_row_to);
//	scalar multiplication of a matrix
void mat_rescale(double *m_new, double *m_old, double scalar, int N);
void mat_rescale(double *m, double scalar, int N);
//	Computes the norm of an array or distance between 2 arrays in R(N)
double norm_L2(double *m1, double *m2, int N);
double norm_L2(double *x, int N);
double norm_sup(double *m1, double *m2, int N);
double norm_sup(int *m1, int *m2, int N);
double norm_sup(double *m, int N);
/*	Computes the sup norm of the percentage difference between m1 and m2 using their average as the denominator. */
double norm_sup_percent(double *m1, double *m2, int N);
//	matrix addition -----------------------------------------------------------------------
void mat_add(double *sum, double *mleft, double *mright, double scalar_mright, int N);
void mat_add(double *sum, double *mleft, double *mright,
			 double scalar_left, double scalar_right, int N);
void mat_add(double *sum, double *m, double scalar, int N);
int mat_positivity(double *m, int N, double eps_tmp = 0);

//	sums the elements of an array
double mat_sum(double *m, int N);
int mat_sum(int *m, int N);
//  matrix multiplication
void mat_mult(double *left_mat, double *right_mat, double *prod_mat,
              int Nrow, int Ncommon, int Ncol);
double mat_mult(double *left_mat, double *right_mat, int N /*dot product*/);
void mat_mult(double *left, double *right, double *elementwiseprod, int N /*element_wise*/);
void mat_mult(double *left, double *right, double *outprod, int Nleft, int Nright);
//	element-wise mat division
void mat_div(double *div, double *num, double *den, int N);
//  Computes mean, deviations from mean and variance of an array
double get_mean(double *x, int N);
double get_mean(double *x, int N, double MISSING_FLAG);
double get_mean(double *x, int Nx, double MISSING_FLAG, int& Nmissing);
void get_deviation_from_mean(double *xdev, double *x, int N);
double get_variance(double *x, int N);
double get_standev(double *x, int N);
double get_correlation(double *m1, double *m2, int N);
//	computes variance and stdev of a series with (assumed) mean zero 
double get_log_variance(double *x, int N);
double get_log_standev(double *x, int N);
void get_percentage_change(double *dx, double *x_bench, double *x_current, int N);
double get_percentage_change(double xBM, double xXP);
/*	matrix interpolation (minterp = theta*m1 + (1-theta)*m2) */
void mat_interpolation(double *minterp, double *m1, double *m2, double theta_m1, int N);
//  Transpose matrice
void transpose_mat(double *mtrans, double *m, int Nr, int Nc);
void mat_transpose(double *mtrans, double *m, int Nr, int Nc);
//  creates a matrix whose rows are basis for R(N) (creates an identity matrix) ------------
void make_basis_vectors(double *u, int N);
//  creates an identity matrix of dim Nx ---------------------------------------------------
void make_identity_matrix(double *mat, int Nx);
/*  returns the max (min) of elements of an array *x of length Nx */
double get_max(double *x, int Nx);
double get_min(double *x, int Nx);
void get_max(double *mat, int Nx, double& max_mat, int& jmax_mat);
void get_max_fabs(double *mat, int Nx, double& max_mat, int& jmax_mat);
int get_max(int *x, int Nx);
int get_min(int *x, int Nx);
//	numerical computation of 1st and 2nd derivative in R(1)
double get_dfdx(double& fx, double x, double (*fct)(double x), double step_df = stepJacob);
double get_d2fx2(double& fx, double x, double (*fct)(double x), double step_df = stepJacob);
/*	<get_jacob> computes jacobian of a function <*fct> from R(N) to R(N) such that J(x) = (f(x+h) - f(x-h))/2h and f(x) = 0.5(f(x+h)+f(x-h)). The step size h can be specified as an absolute or as a percentage of x. <get_gradient> computes the gradient of a function from R(N) to R. <get_dfdx_along_s> and <get_d2fdx2_along_s> computes the 1st and 2nd derivatives of a function along the direction s. */
void get_jacob(double *x, double *jacob, double *fx, int N,
			   void (*fct)(double *fx, double *x), double step_jacob =  stepJacob);
void get_jacob_old(double *x, double *jacob, double *fx, int N,
               void (*fct)(double *fx, double *x), double step_jacob = stepJacob);
void get_jacob(double *x, double *jacob, double *fx, void (*fct)(double *fx, double *x),
			   int N, double percent_step_jacob = 0.01);
void get_gradient(double *Jx, double *x, double& fx, int Nx, double (*fct)(double *x),
				  double step_jacob =  stepJacob);
void get_gradient(double *Jx, double *x, double& fx, double (*fct)(double *x),
				  int Nx, double percent_change = 0.01);
double get_dfdx_along_s(double& fx, double *x, double *s_dir, int Nx,
						double (*fct)(double *x), double step_jacob =  linemin_stepJacob);
//	3 point approximation without calling get_dfdx_along_s
double get_d2fdx2_along_s(double &fx, double *x, double *sdir, int Nx,
						  double (*fct)(double *x), double h, double& dfdx,
						  int flag = DIFF_3PT_FLAG);
double get_d2fdx2_along_s(double& fx, double *x, double *s_dir, int Nx,
						  double (*fct)(double *x), double step_jacob = linemin_stepJacob);
//  matrix inverter - Given matrix and current row (pivot), finds the next candidate pivot row. */
int get_pivot(double *mat /* ptr to matrix */, int current_row, int N_rows);
/*	Swaps rows pivot row and cur_row in square matrices m and minv */
void swap_row(double *m, double *minv, int pivot_row, int cur_row, int N_rows);
/*	Gauss-Jordan mat inv - sets all entries below pivot to zero	*/
void zero_below(double *m, double *minv, int pivot, int N_rows);
/*	Guass-Jordan mat inv - sets all entries above pivot to zero	*/
void zero_above(double *m, double *minv, int pivot, int N_rows);
/*	Guass-Jordan matrix inversion of square matrix m_ptr with #rows=N_rows. */
int inv_mat(double *m_inv_ptr, double *m_ptr, int N_rows);
int inv_mat(double *m, int Nrows);
//  computes a smooth percentage change of variables over Nt periods given base period values and end period values where
//  xend[i] = pow(dx[i], Nt) * xbase[i] forall i
void get_growth_rate(double *dx, double *xbase, double *xend, int Nx, int Nt);
/*	Computes the condition number of a matrix using sup norm. Note the matrix needs to be inverted. */
double cond_sup(double *m, int N);
/* <get_projection> computes the projection of vector y on vector x */
void get_projection(double *yproj, double *epshat, double *y, double *x, int N);
/* Computes the correlation of 2 arrays in R(N) */
double vector_correlation(double *v1, double *v2, int N);
//	Gauss-Legendre quadrature --------------------------------------------------------------
 /*  given lower and upper limits of integration (x1,x2) and given n, the number of quadrature points (degree of approximation), this function returns arrays, x[] and w[] for the zeroes and weights of the Gauss-Legendre n-point quadrature formula. Adapted from Numerical Recipes for C */
void gauleg(double x1 /*lower_bound*/, double x2 /*upper_bound*/,
			double *z_ptr /*zeroes*/, double *w_ptr /*weights*/, int N);
/*	Given zeroes and weights for Gauss-Legendre quadrature defined over interval [0,1], returns rescaled zeroes and weights for interval [x_lower, x_upper] */
void gauleg_rescale(double x_lower, double x_upper, double *z_base, double *z_rescaled,
					double *w_base, double *w_rescaled, int N = Nquad);
/*	Rescales candidate solution x - NEED TO FIX RETURN VALUE	*/
int rescale_xnext(double *x_ptr, double *Jinv_ptr, double *fx_ptr, double *xx_ptr,
                  int (*fct_constraints)(double *x), int N);
/*	For use with line_min_newton. */
int rescale_xnext(double *x, double *xnext, int Nx, int (*fct_constraints)(double *x),
				  int Maxiter_tmp = Nmaxiter_xnext, double rescale_factor_tmp = 2.0);
int rescale_xnext(double x, double& xnext, int (*fct_constraints)(double x),
				  int Nmaxiter_tmp = Nmaxiter_xnext, double rescale_factor_tmp = 2.0);
/* No constraints function ptr for the case of newton_method or broyden_method without constraints on solution vector x. */
int no_constraints(double *x);
int no_constraints(double x);
//	------------------------------------------------------------
double linear_interpolate(double *x, double *y, int Nx, int Ny, double *fxy,
						  double xstar, double ystar);
int get_index_linear_interpolate(double xstar, double *x, int Nx);
//	generates an array of uniform deviates between 0,1 (Minimum standard for random number generator of Park and Miller, see Numerical recipes 7.1 for details)
void genr_random_array(double *mat, int N, long int seed);
double genr_random_ran0(long *idum, int IA, int IM, double AM, int IQ, int IR, long int MASK);



//  to be made obsolete --------------------------------------------------------------------
double jacob_xi(double *jj_ptr, double *fx_ptr, int i, int N);
