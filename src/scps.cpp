#include <Rcpp.h>
using namespace Rcpp;

//**********************************************
// Author: Anton Grafström
// Last edit: 2014-03-18 
// Licence: GPL (>=2)
//**********************************************

void scps_swap(IntegerVector& array, int i, int j) {double t=array[i]; array[i]=array[j]; array[j]=t; }

bool scps_compare(double d1, double d2, double w1, double w2){
  if(d1 < d2){return true;}
  if(d1 == d2){
    if(w1<w2){return true;}
  }
  return false;
}

void scps_quicksort(IntegerVector& index, NumericVector& dists, NumericVector& weights, int left, int right) {
 
    if (left < right) {
 
      int pivot = (left + right)/2;
      int left_new = left, right_new = right;
      do {
        while (scps_compare(dists[index[left_new]],dists[index[pivot]],weights[index[left_new]],weights[index[pivot]]))
          left_new++;
        while (scps_compare(dists[index[pivot]],dists[index[right_new]],weights[index[pivot]],weights[index[right_new]]))
          right_new--;
        if (left_new  <= right_new)
        
          scps_swap(index, left_new++, right_new--);
          
      } while (left_new  <= right_new);
 
      scps_quicksort(index,dists,weights, left, right_new);
      scps_quicksort(index,dists,weights, left_new, right);
 
    }
    //return index;
  }



// [[Rcpp::export]]
IntegerVector scps(NumericVector prob, NumericMatrix x){
	int N = prob.size();
	int ncol = x.ncol();
	NumericVector dists(N), weights(N), p(N);
  IntegerVector index(N), s(N);
  
	double d,dp=0.0,weight,uw,eps=1e-9;
	
	for(int i=0;i<N;i++){p[i]=prob[i];}
		
	NumericVector rand = runif(N);
	NumericVector rand2 = runif(N);
	int nr, nr_others, i,j,k, m;

	
	for(i=0;i<N;i++){
    // make decision for unit i
		if(rand[i] < p[i]){s[i]=1;}else{s[i]=0;}	
    
    // if not last unit, update others
		if(i<N-1){	
      // find distances and max weights
			for(k = i+1; k < N; k++){
				d = 0.0;
				for(m=0;m<ncol;m++){
					dp = x(i,m)-x(k,m);
					d += dp*dp;
				}
        
				dists[k-i-1] = d;
				weights[k-i-1] = std::min(p[k]/(1-p[i]),(1-p[k])/p[i]);
        index[k-i-1] = k-i-1;
			}
				
			// Sort index w.r.t. dists, weights
			scps_quicksort(index,dists,weights,0,N-i-2);
			
			// Weight to distribute
			weight = 1.0;
      
      // Check if i can distribute weight
			if(p[i] < 1 && p[i]>0){
        
        // Nr of potential recievers (at least 1)
				nr_others = N-i-1;
         
				for(k = 0; k < nr_others; k++){
          // Nr at equal distance 
          nr = 1;
          for(j=k+1;j<nr_others;j++){
            if(dists[index[k]]==dists[index[j]]){
  						nr = nr + 1;
						}else{break;}
          }	
					uw = std::min(weight/nr,weights[index[k]]);
					// Update
					p[i+index[k]+1] = p[i+index[k]+1]-(s[i]-p[i])*uw;
					// Remove used weight
					weight = weight-uw;
					if(weight<eps){break;}
				}
			}	
		}
	}
  // make sample instead of indicator
  int n = sum(s), count=0;  
  IntegerVector sa(n);
  for(i=0;i<N;i++){
    if(s[i]==1){
      sa[count]=i+1;
      count++;
    }
  }
	return sa;
}
