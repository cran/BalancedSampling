#include <Rcpp.h>
using namespace Rcpp;

//**********************************************
// Author: Anton Grafström
// Last edit: 2014-03-18 
// Licence: GPL (>=2)
//**********************************************

// [[Rcpp::export]]
NumericVector lpm(NumericVector prob, NumericMatrix x, int h){
	int N = prob.size();
	int ncol = x.ncol();
	NumericVector index(N), p(N);
	
	double wp,lp,mindist,d,dp=0.0;
	int ri,rj;
	
	for(int i=0;i<N;i++){index[i]=i; p[i]=prob[i];}
		
	NumericVector r1 = runif(N);
	NumericVector r2 = runif(N);
	int move;
	int temp;
	double a;
	int maxFrom,maxTo;
	
	for(int i=0;i<N-1;i++){	
		ri = i+floor((r1[i]*(N-i)));
		rj = 0;
		mindist = 1e+200;
		maxFrom = ri - h/2;
		if( maxFrom < i ){
			maxFrom = i;
			
		}
		maxTo = maxFrom + h;
		if( maxTo > N ){
			maxTo = N;	
		}
		for(int j=maxFrom;j<maxTo;j++){
			if(j!=ri){
				d = 0.0;
				for(int k=0;k<ncol;k++){
					dp = x(index[ri],k)-x(index[j],k);
					d += dp*dp;
				}
				if(d < mindist){
					rj = j;
					mindist = d;
				}
			}
		}
		a = p[index[ri]]+p[index[rj]];
		if(a > 1){
			wp = 1;
		}else{
			wp = a;
		}
		lp = p[index[ri]]+p[index[rj]]-wp;
		if( r2[i] < (wp-p[index[rj]])/(wp-lp) ){
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
