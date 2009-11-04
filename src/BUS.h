#include <iostream>
using namespace std;
#include <algorithm>
#include <ctime>
#include <cmath>
#include <vector>
#include <map>
#include <R.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#define NA -2000000
/* These functions are defined in the mim.cpp file */

double entropy_empirical(map< vector<double> ,int > frequencies, int nb_samples);
double mutualinfo(const double* d1, const double* d2,  int N, int n1,  int n2, int i, int j);
/* Entry points called from the R functions */
extern "C" 
{

SEXP MINempirical(SEXP data1,SEXP data2, SEXP nrows, SEXP ncols1, SEXP ncols2);


}

