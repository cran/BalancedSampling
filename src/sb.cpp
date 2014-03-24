#include <Rcpp.h>
using namespace Rcpp;

//**********************************************
// Author: Anton Grafström
// Last edit: 2014-03-18 
// Licence: GPL (>=2)
//**********************************************

// [[Rcpp::export]]
double sb(NumericVector p, NumericMatrix x, NumericVector s){
	int N = x.nrow();
	int ncol = x.ncol();
	int n = s.size();
	NumericVector incl(N);
	NumericVector nearest(n);
	int nrNearest = 0;
	double mindist, d, res = 0.0;
    for(int i=0;i<N;i++){incl[i]=0.0;}        
	for(int i=0;i<N;i++){
		nrNearest = 0;
		mindist = 1e+200; 
		for(int j=0;j<n;j++){
			d = 0.0;
			for(int k=0;k<ncol;k++){
				d += pow(x(i,k)-x(s[j]-1,k),2);
			}			
			if(d==mindist){
				nearest[nrNearest]=s[j];
				nrNearest++;
			}
			if(d<mindist){
				mindist = d;
				nearest[0] = s[j];
				nrNearest = 1;
			}
		}
		for(int k=0;k<nrNearest;k++){
			incl[nearest[k]-1] = incl[nearest[k]-1] + p[i]/nrNearest;
		}
	}
	for(int i=0;i<n;i++){
		res += pow(incl[s[i]-1]-1,2);
	}	
	return res/n;
}
