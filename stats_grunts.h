
int detrend_series(double *series, double *series_trend, double *series_pctdev, double& grow_rate, int Nt);
double fx_detrend_series(double *x);

int detrend_series_fixedrate(double *series, double *series_trend, double *series_pctdev, 
							 double grow_rate, int Nt);
double fx_detrend_series_fixedrate(double x0);

int HP_filter(double *series, double *series_trend, double *series_pdev, double lambda, int Nt);
double fx_HP_filter(double *x);
