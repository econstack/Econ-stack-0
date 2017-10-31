
/*	-----------------------------------------------------------------------------------------
	Checks whether x(k+1)=xnext is within lim_dxpct percent of x(k)=x. Rescales if needed. Option 
	to write out diagnostic to outfile. Used by nlin solvers and linemin solvers.
--------------------------------------------------------------------------------------------- */
int get_xnext_dpercent(double *xnext, 
    double *x, 
    const int& Nx, 
    const double& lim_dxpct, 
    ofstream *fout);
void get_xnext_dpercent_fout(ofstream& fout, 
    double *x, 
    double *xnext, 
    double *dx, 
    double *dxpct, 
	const int& Nx, 
    const double& max_dxpct, 
    const int& j_max_dxpct,
    const double& limit_dxpct);
