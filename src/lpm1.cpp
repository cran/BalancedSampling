#include <Rcpp.h>
using namespace Rcpp;

//**********************************************
// Author: Anton Grafström
// Last edit: 2014-12-16 
// Licence: GPL (>=2)
//**********************************************

// [[Rcpp::export]]
NumericVector lpm1(NumericVector prob, NumericMatrix x){
	int N = prob.size();
	int ncol = x.ncol();
	NumericVector index(N), p(N), nearest(N);
	
	double wp,lp,mindist,mindist2,d,dp=0.0;
	int ri,rj,count;
	
	for(int i=0;i<N;i++){index[i]=i; p[i]=prob[i];}
		
	NumericVector r1 = runif(N);
	NumericVector r2 = runif(N);
	int move;
	int temp;
	int found = 0;
	double a;
	
	for(int i=0;i<N-1;i++){	
	  ri = i+floor((r1[i]*(N-i)));
		found = 0;
		while(found==0){	
			count = 0;
			mindist = 1e+200;
			// find nearest
			for(int j=i;j<N;j++){
				if(j!=ri){
					d = 0.0;
					for(int k=0;k<ncol;k++){
						dp = x(index[ri],k)-x(index[j],k);
						d += dp*dp;
					}
					if(d==mindist){
						nearest[count] = j;
						count = count + 1;
						
					}
					if(d < mindist){
						nearest[0] = j;
						count = 1;
						mindist = d;
					}
				}
			}
			rj = nearest[floor(count*runif(1)[0])];
			// check that ri is also nearest to rj
			mindist2 = 1e+200;
			// find nearest
			for(int j=i;j<N;j++){
				if(j!=rj){
					d = 0.0;
					for(int k=0;k<ncol;k++){
						dp = x(index[rj],k)-x(index[j],k);
						d += dp*dp;
					}
					if(d < mindist2){
						mindist2 = d;
					}
				}
			}
			// found if equal mindist (changed from abs to fabs)
			if( fabs(mindist-mindist2) < 1e-6 ){found=1;}else{ri=rj;}
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
	count = 0;
	for(int i=0;i<N;i++){
		if(p[i]==1){
			s[count]=i+1;
			count++;
		}
	}
	return s;
}
