#include <Rcpp.h>
using namespace Rcpp;

//**********************************************
// Author: Anton Grafström
// Last edit: 2014-03-18 
// Licence: GPL (>=2)
//**********************************************

// [[Rcpp::export]]
double vsb(NumericVector probs, NumericVector ys, NumericMatrix xs){
	int ncol = xs.ncol();
	int n = ys.size();
	double mindist, d, dp=0.0, res = 0.0;
  NumericVector ni(n), sydp(n);
        
	for(int i=0;i<n;i++){
		mindist = 1e+200; 
		for(int j=0;j<n;j++){
			if( i != j ){
				d = 0.0;
				for(int k=0;k<ncol;k++){
					dp = xs(i,k)-xs(j,k);
					d += dp*dp;
				}
        if(d==mindist){
          ni[i] = ni[i] + 1;
          sydp[i] = sydp[i] + ys[j]/probs[j];
        }
				if(d<mindist){
					mindist = d;
          ni[i] = 2;
          sydp[i] = ys[i]/probs[i] + ys[j]/probs[j];
				}
			}
		}
    res += ni[i]/(ni[i]-1)*(ys[i]/probs[i]-1/ni[i]*sydp[i])*(ys[i]/probs[i]-1/ni[i]*sydp[i]);
	}
	return res;
}
