#include "BUS.h"


SEXP MINempirical(SEXP Rdata1, SEXP Rdata2, SEXP Rnrows, SEXP Rncols1, SEXP Rncols2)
{
      const double *data1, *data2;
      const int *nrows,*ncols1,* ncols2;
      double *res, mi;
         SEXP Rres;
         PROTECT(Rdata1 = AS_NUMERIC(Rdata1));
         PROTECT(Rdata2 = AS_NUMERIC(Rdata2));
         PROTECT(Rnrows= AS_INTEGER(Rnrows));
         PROTECT(Rncols1= AS_INTEGER(Rncols1));
         PROTECT(Rncols2= AS_INTEGER(Rncols2));
         data1 = NUMERIC_POINTER(Rdata1);
         data2 = NUMERIC_POINTER(Rdata2);
         nrows= INTEGER_POINTER(Rnrows);
         ncols1= INTEGER_POINTER(Rncols1);
         ncols2= INTEGER_POINTER(Rncols2);       
         PROTECT(Rres = NEW_NUMERIC((*ncols1)*(*ncols2)));
         res = NUMERIC_POINTER(Rres);
      for( int i=0; i<*ncols1; ++i )
         for( int j=0; j< *ncols2; ++j ) {
                  mi = mutualinfo( data1, data2, *nrows, *ncols1, *ncols2, i, j );
                  res[j*(*ncols1)+i]  = mi;
         }
         UNPROTECT(6);
      return Rres;
}

double mutualinfo(const double *d1, const double *d2, int N, int n1, int n2, int i, int j) {
         map< vector<double> ,int > freqi;
         map< vector<double> ,int > freqj;
         map< vector<double> ,int > freqij;
         vector<double> sel;
      double Hi, Hj, Hij;
      int    ni=0, nj=0, nij=0;
      for(int s = 0; s < N; s++) 
        if(d1[s+i*N]!=NA && d2[s+j*N]!=NA ){
              sel.clear();
              sel.push_back(d1[s+i*N]);  freqi[sel]++; ni++;
              sel.push_back(d2[s+j*N]);  freqij[sel]++; nij++; 
              sel.clear();
              sel.push_back(d2[s+j*N]);  freqj[sel]++; nj++;
        }            
        else if( d1[s+i*N]!=NA ) { sel.push_back(d1[s+i*N]);  freqi[sel]++; ni++; }
        else if( d2[s+j*N!=NA] ) { sel.push_back(d2[s+j*N]);  freqj[sel]++; nj++; }                  

      //empirical
            Hi = entropy_empirical(freqi,ni);
            Hj = entropy_empirical(freqj,nj);
            Hij = entropy_empirical(freqij,nij);
      
      
      double mi = (Hi+Hj-Hij);
      if(mi<0) return 0;
      return mi;      
}

double entropy_empirical(map< vector<double> ,int > frequencies, int nb_samples) {
      double e = 0;
      for (map< vector<double> ,int>::const_iterator iter = frequencies.begin(); iter != frequencies.end(); ++iter)
            e -= iter->second * log(iter->second);
      return log(nb_samples) + e/nb_samples;
}



