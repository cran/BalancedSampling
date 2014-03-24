#include <Rcpp.h>
using namespace Rcpp;

//**********************************************
// Author: Anton Grafström
// Last edit: 2014-03-18 
// Licence: GPL (>=2)
//**********************************************

// [[Rcpp::export]]
NumericVector rpm(NumericVector prob){
	int N = prob.size();
	NumericVector index(N), p(N);
	
	double wp,lp;
	int ri,rj;
	
	for(int i=0;i<N;i++){index[i]=i; p[i]=prob[i];}
		
	NumericVector r1 = runif(N);
	NumericVector r2 = runif(N);
  NumericVector r3 = runif(N);
	int move;
	int temp;
	double a;


	for(int i=0;i<N-1;i++){	
		ri = i+floor((r1[i]*(N-i)));
		rj = i+floor((r2[i]*(N-i-1)));
    if(rj>=ri){rj=rj+1;}
		a = p[index[ri]]+p[index[rj]];
		if(a > 1){
			wp = 1;
		}else{
			wp = a;
		}
		lp = p[index[ri]]+p[index[rj]]-wp;
		if( r3[i] < (wp-p[index[rj]])/(wp-lp) ){
			p[index[ri]] = wp;
			p[index[rj]] = lp;
		}else{
			p[index[ri]] = lp;
			p[index[rj]] = wp;
		}
		if(i==N-2){
			if(runif(1)[0]<p[index[ri]]){
				p[index[ri]]=1;
			}else{
				p[index[ri]]=0;
			}
			if(runif(1)[0]<p[index[rj]]){
				p[index[rj]]=1;
			}else{
				p[index[rj]]=0;
			}
		}	
		move = rj;
		if( p[index[ri]]==0 || p[index[ri]]==1 ){
			move = ri;
		}
		temp = index[i];
		index[i] = index[move];
		index[move] = temp;
	}
	int n = round(sum(p));
	NumericVector s(n);
	int count = 0;
	for(int i=0;i<N;i++){
		if(p[i]==1){
			s[count]=i+1;
			count++;
		}
	}
	return s;
}
