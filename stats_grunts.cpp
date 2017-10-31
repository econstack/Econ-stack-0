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
#include "stats_grunts.h"
#include "fmin.h"

int fout_detrend_series_FLAG = ON;
const int Nt_detrend_series_max = 100;
double Gseries[Nt_detrend_series_max]; // global accessed from feval for detrend series
int Nt_detrend_series = 10; // global accessed from feval for detrend series
//	detrends a series with Nt periods, choosing a growth rate and init value to minimize the sum of square residuals, returns the series trend as well as the percent dev from trend
int detrend_series(double *series, double *series_trend, double *series_pctdev, double& grow_rate, int Nt)
{
	//	buffer overflow check - series has too many periods
	if (Nt>=Nt_detrend_series_max) {
		log_fout << endl << "Error, detrend_series()!!! Exceeded max length of series =" 
			<< Nt_detrend_series_max << "!!! Exiting...";
		return OVERFLOW_FLAG;
	}
	//	turn off (toggle) output 
	toggle_fout_nlin_FLAGS();

	//	compute init guess for growth rate
	grow_rate = pow(series[Nt-1] / series[0], 1.0/(Nt-1));
	//	store targets in global
	Nt_detrend_series = Nt;
	copy_array(&Gseries[0], series, Nt);
	//	set initial guess for detrend
	double x[2], xstar[2];
	x[0] = grow_rate;
	x[1] = series[0];
	//	choose growth rate to min sum of square residuals
	double *Hx; Hx = new double [2*2];
	ParamFmin pm_fmin;
	pm_fmin.Nfout_fmin = 0;
	pm_fmin.Nfout_LM = 0;
	pm_fmin.eps_fmin *= 0.01;
	pm_fmin.eps_LM *= 0.01;
	int Nx = 2;
	int flag = fmin_BFGS(&x[0], &xstar[0], Nx, Hx, fx_detrend_series, no_constraints, pm_fmin);

	//	return results
	grow_rate = xstar[0];
	series_trend[0] = xstar[1];
	for (int t=1; t<Nt; ++t) 
		series_trend[t] = series_trend[t-1] * grow_rate;
	for (int t=0; t<Nt; ++t)
		series_pctdev[t] = (series[t] - series_trend[t]) / series_trend[t];

	if (fout_detrend_series_FLAG == ON) {
		log_fout << endl << "series";
		fout_array(log_fout, series, 1, Nt, 6);
		log_fout << endl << "series_trend";
		fout_array(log_fout, series_trend, 1, Nt, 6);
		log_fout << endl << "series_pctdev";
		fout_array(log_fout, series_pctdev, 1, Nt, 6);
	}

	//	turn on (toggle) output 
	toggle_fout_nlin_FLAGS();

	delete[] Hx;
	return flag;
}
double fx_detrend_series(double *x)
{
	double resid, norm_fx = 0;
	double *series_trend; series_trend = new double [Nt_detrend_series];
	//	extract into labelled vars
	double grow_rate = x[0];
	series_trend[0] = x[1];
	//	compute implied trend path
	for (int t=1; t<Nt_detrend_series; ++t) 
		series_trend[t] = series_trend[t-1] * grow_rate;
	//	compute sum of square residuals
	for (int t=0; t<Nt_detrend_series; ++t) {
		resid = (Gseries[t] - series_trend[t]) / series_trend[t];
		norm_fx += pow(resid, 2.0);
	}
	delete[] series_trend;
	return norm_fx;
}
//	detrends a series with Nt periods with fixed growth rate - chooses an init value to min the sum sq resids and returns the series trend and the pct dev from trend
double Gseries_grow_rate; // glocal access from feval for detrend_series_fixedrate
int detrend_series_fixedrate(double *series, double *series_trend, double *series_pctdev, 
							 double grow_rate, int Nt)
{
	//	buffer overflow check - series has too many periods
	if (Nt>=Nt_detrend_series_max) {
		log_fout << endl << "Error, detrend_series()!!! Exceeded max length of series =" 
			<< Nt_detrend_series_max << "!!! Exiting...";
		return OVERFLOW_FLAG;
	}

	//	set init guess for init value
	double x0 = series[0];
	//	store targets in global
	Nt_detrend_series = Nt;
	copy_array(&Gseries[0], series, Nt);
	Gseries_grow_rate = grow_rate;
	//	choose init value to min sum of square residuals
	double x0_star, df=0, d2f=0, loss_fx=0;
	ParamFmin pm_fmin;
	pm_fmin.Nfout_fmin = 0;
	pm_fmin.Nfout_LM = 0;
	pm_fmin.eps_fmin *= 0.05;
	int flag = fmin_R1(x0, x0_star, df, d2f, loss_fx, fx_detrend_series_fixedrate, no_constraints, pm_fmin);
	//	return results
	series_trend[0] = x0_star;
	for (int t=1; t<Nt; ++t) 
		series_trend[t] = series_trend[t-1] * grow_rate;
	for (int t=0; t<Nt; ++t)
		series_pctdev[t] = (series[t] - series_trend[t]) / series_trend[t];

	//	option: write out diag
	if (fout_detrend_series_FLAG == ON) {
		log_fout << endl << "series";
		fout_array(log_fout, series, 1, Nt, 6);
		log_fout << endl << "series_trend";
		fout_array(log_fout, series_trend, 1, Nt, 6);
		log_fout << endl << "series_pctdev";
		fout_array(log_fout, series_pctdev, 1, Nt, 6);
	}
	return flag;
}
double fx_detrend_series_fixedrate(double x0)
{
	double resid, norm_fx = 0;
	double *series_trend; series_trend = new double [Nt_detrend_series];
	//	extract into labelled vars
	series_trend[0] = x0;
	//	compute implied trend path
	for (int t=1; t<Nt_detrend_series; ++t) 
		series_trend[t] = series_trend[t-1] * Gseries_grow_rate;
	//	compute sum of square residuals
	for (int t=0; t<Nt_detrend_series; ++t) {
		resid = (Gseries[t] - series_trend[t]) / series_trend[t];
		norm_fx += pow(resid, 2.0);
	}
	delete[] series_trend;
	return norm_fx;
}
//	HP filter: detrends a series with Nt periods with a Hodrick Prescott filter which min the weighted sum of the sum of sq resid between a series and its HP trend and the sum of sq resid of the 2nd diff of the HP trend
double lambda_HP_filter = 100; // annual = 100.0 and quarterly 1600
int HP_filter(double *series, double *series_trend, double *series_pdev, double lambda, int Nt)
{
	//	buffer overflow check - series has too many periods
	if (Nt>=Nt_detrend_series_max) {
		log_fout << endl << "Error, detrend_series()!!! Exceeded max length of series =" 
			<< Nt_detrend_series_max << "!!! Exiting...";
		return OVERFLOW_FLAG;
	}
	//	turn off (toggle) output 
	toggle_fout_nlin_FLAGS();

	//	store targets in global
	Nt_detrend_series = Nt;
	copy_array(&Gseries[0], series, Nt);
	lambda_HP_filter = lambda;
	//	set initial guess for HP filter as linear filter
	double *x; x = new double [Nt];
	double gr_factor = pow(series[Nt-1]/series[0], 1.0/(Nt-1));
	x[0] = series[0];
	for (int t=1; t<Nt; ++t) x[t] = gr_factor * x[t-1];
	//	choose growth rate to min sum of square residuals
	double *Hx; Hx = new double [Nt*Nt];
	double *xstar; xstar = new double [Nt];
	ParamFmin pm_fmin;
	pm_fmin.Nmax_iter_fmin = 100;
	pm_fmin.Nfout_fmin = 0;
	pm_fmin.Nfout_LM = 0;
	pm_fmin.eps_fmin *= 0.1;
	pm_fmin.eps_LM *= 0.1;
	int flag = fmin_BFGS(x, xstar, Nt, Hx, fx_HP_filter, no_constraints, pm_fmin);

	//	return results
	copy_array(&series_trend[0], xstar, Nt);
	for (int t=0; t<Nt; ++t)
		series_pdev[t] = (series[t] - series_trend[t]) / series_trend[t];

	if (fout_detrend_series_FLAG == ON) {
		log_fout << endl << "series";
		fout_array(log_fout, series, 1, Nt, 6);
		log_fout << endl << "series_trend";
		fout_array(log_fout, series_trend, 1, Nt, 6);
		log_fout << endl << "series_pctdev";
		fout_array(log_fout, series_pdev, 1, Nt, 6);
	}

	//	turn on (toggle) output 
	toggle_fout_nlin_FLAGS();

	delete[] Hx; delete[] xstar; delete[] x;
	return flag;
}
double fx_HP_filter(double *x)
{
	double norm_fx = 0;
	double diff1, diff2;
	double *series_trend; series_trend = new double [Nt_detrend_series];
	//	extract double* into labelled var
	for (int t=0; t<Nt_detrend_series; ++t) series_trend[t] = x[t];
	//	compute loss from sum of sq resid between series and its HP trend
	for (int t=0; t<Nt_detrend_series; ++t)
		norm_fx += pow(Gseries[t]-series_trend[t], 2.0);
	for (int t=1; t<(Nt_detrend_series-1); ++t) {
		diff1 = series_trend[t+1] - series_trend[t];
		diff2 = series_trend[t] - series_trend[t-1];
		norm_fx += lambda_HP_filter * pow(diff1 - diff2, 2.0);
	}
	delete[] series_trend;
	return norm_fx;
}