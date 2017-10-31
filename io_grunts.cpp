#include <math.h>
#include <fstream>
#include <string>
#include <conio.h>
#include <iostream>
#include <iomanip>
#include <Windows.h>

using namespace std;

#include "math_parameters.h"
#include "io_grunts.h"

int Nprec_default = 4;

//	WRITES OUT ARRAYS TO FILE
void fout_array(ofstream& fout, double *mat, int Nrows, int Ncols, int Nprecision,
				int newline_flag)
{
	if (newline_flag == ON) fout << endl;
    fout << setiosflags(ios::fixed) << setprecision(Nprecision);

    for (int i=0; i<Nrows; ++i)  {
        for (int j=0; j<Ncols; ++j)  {
			if ((fabs(*(mat+(i*Ncols)+j))>10000)||(fabs(*(mat+(i*Ncols)+j))<0.00001)) {
				fout.setf(std::ios::scientific);
				fout.precision(4);
				fout << setw(12) << *(mat+(i*Ncols)+j);
				fout.unsetf(std::ios::scientific);
				fout.precision(Nprecision);
			} else {
				fout << setw(12) << *(mat+(i*Ncols)+j);
			}
        }   //  end for(j)
        fout << endl;
    }   //  end for(i)

    return;
}
void fout_array(ofstream& fout, int *array, int Nrows, int Ncols, int newline_flag)
{
	if (newline_flag == ON) fout << endl;
     for (int i=0; i<Nrows; ++i)  {
        for (int j=0; j<Ncols; ++j)  {
            fout << setw(12) << *(array+(i*Ncols)+j);
        }   //  end for(j)
        fout << endl;
    }   //  end for(i)

    return;
}
void fout_array(ofstream& fout, std::string *array, int Nrows, int Ncols,
				int newline_flag)
{
	if (newline_flag == ON) fout << endl;
	for (int i_row=0; i_row<Nrows; ++i_row)  {
        for (int i_col=0; i_col<Ncols; ++i_col)  {
            fout << setw(12) << *(array+i_row*Ncols+i_col);
        }   //  end for(i_col)
        fout << std::endl;
    }   //  end for(i_row)
    return;
}
void fout_array(ofstream& fout, double *array, int Nrows, int Ncols, int Nprec, int newline_flag, int Nsetw)
{
	if (newline_flag == ON) fout << endl;
    fout << setiosflags(ios::fixed) << setprecision(Nprec);
    for (int i=0; i<Nrows; ++i)  {
        for (int j=0; j<Ncols; ++j)  {
            fout << setw(Nsetw) << *(array+(i*Ncols)+j);
        }   //  end for(j)
        fout << endl;
    }   //  end for(i)
    return;
}
//	writes out submatrix of a matrix
void fout_submatrix(ofstream& fout, double *mat, int Nr, int Nc, int startrow,
					int startcol, int Nr_sub, int Nc_sub, int Nprec)
{

	fout << endl;
	//	diagnostics
	if ((Nr_sub < 1) || (Nc_sub < 1)) {
		fout << "Error: fout_submatrix. Nonpositive submat dimensions. Exiting...";
		return;
	}
	if (((startrow+Nr_sub)>Nr) || ((startcol+Nc_sub)>Nc)) {
		fout << "Error: fout_submatrix. Submat dimensions bigger than original.";
		return;
	}

	//	write out submatrix
	for (int i=startrow; i<startrow+Nr_sub; ++i) {
		fout_array(fout, mat+i*Nc+startcol, 1, Nc_sub, Nprec, 0);
	}

	return;
}
//	READS IN ARRAYS FROM FILE
void fin_array(ifstream& fin, double *array_ptr, int Nrow, int Ncol)
{
	for (int i=0; i<Nrow; ++i)  {
        for (int j=0; j<Ncol; ++j)  {
            fin >> *(array_ptr+(i*Ncol)+j);
        }   //  end for(j)
    }   //  end for(i)

    return;
}

void fin_array(ifstream& fin, int *array_ptr, int Nrow, int Ncol)
{
	for (int i=0; i<Nrow; ++i)  {
        for (int j=0; j<Ncol; ++j)  {
            fin >> *(array_ptr+(i*Ncol)+j);
        }   //  end for(j)
    }   //  end for(i)

    return;
}

void fin_array(ifstream& fin, std::string *array_ptr, int Nrow, int Ncol)
{
	for (int i=0; i<Nrow; ++i)  {
        for (int j=0; j<Ncol; ++j)  {
            fin >> *(array_ptr+(i*Ncol)+j);
        }   //  end for(j)
    }   //  end for(i)

    return;
}
//	writes out current system date and time to file
void fout_system_date(ofstream& fout)
{
	SYSTEMTIME st;
    GetSystemTime(&st);
	fout << endl << "year:month:day:hour:min "
		<< st.wYear << ":"
		<< st.wMonth << ":"
		<< st.wDay << ":"
		<< st.wHour << ":"
		<< st.wMinute;
	return;
}