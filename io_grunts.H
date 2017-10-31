int extern Nprec_default;
int extern ON;

//  --------- WRITES OUT ARRAYS TO FILE ------------------------------------------------
void fout_array(ofstream& , double* , int , int , int , int newline_flag = ON);
void fout_array(ofstream& , int* , int , int , int newline_flag = ON);
void fout_array(ofstream& , string *, int , int , int newline_flag = ON);
void fout_array(ofstream& fout, double *array, int Nrows, int Ncols, int Nprec, int newline_flag, int Nsetw);
void fout_submatrix(ofstream& fout, double *mat, int Nr, int Nc, int startrow,
					int startcol, int Nr_sub, int Nc_sub, int Nprec);
//	---------  READS IN ARRAYS FROM FIL ------------------------------------------------
void fin_array(ifstream& , double*, int , int );
void fin_array(ifstream& , int*, int , int );
void fin_array(ifstream& , string *, int , int );

//	writes out current system date and time to file
void fout_system_date(ofstream& fout);
