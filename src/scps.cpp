#include <Rcpp.h>
using namespace Rcpp;

//**********************************************
// Author: Anton Grafström
// Last edit: 2018-10-07 
// Licence: GPL (>=2)
//**********************************************

void swap(IntegerVector& ord, int a, int b) {
  int t = ord[a];
  ord[a] = ord[b];
  ord[b] = t;
  return;
}

int partition(IntegerVector& ord, NumericVector& dists, NumericVector& weights, int low, int high) {
  int ord_high = ord[high];
  double pivot_dists = dists[ord_high];
  double pivot_weights = weights[ord_high];
  int i = (low-1);
  
  for (int j = low; j <= high-1; j++) {
    int ord_j = ord[j];
    if (
        dists[ord_j] < pivot_dists ||
          (
              dists[ord_j] == pivot_dists && weights[ord_j] <= pivot_weights
          )
    ) {
      swap(ord, ++i, j);
    }
  }
  
  swap(ord, ++i, high);
  return (i);
}

void quicksort(IntegerVector& ord, NumericVector& dists, NumericVector& weights, int low, int high) {
  if (low < high) {
    int pi = partition(ord, dists, weights, low, high);
    quicksort(ord, dists, weights, low, pi-1);
    quicksort(ord, dists, weights, pi+1, high);
  }
  return;
}

// [[Rcpp::export]]
IntegerVector scps(NumericVector prob, NumericMatrix x){
  int N = prob.size();
  int ncol = x.ncol();
  NumericVector dists(N), weights(N), p(N);
  IntegerVector index(N), s(N);
  NumericVector rand = runif(N);
  
  double d,dp=0.0,weight,uw,eps=1e-9;
  
  for(int i=0;i<N;i++){p[i]=prob[i];}
  
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
      quicksort(index, dists, weights, 0, N-i-2);
      
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

// [[Rcpp::export]]
IntegerVector scps_coord(NumericVector prob, NumericMatrix x, NumericVector rand){
  int N = prob.size();
  int ncol = x.ncol();
  NumericVector dists(N), weights(N), p(N);
  IntegerVector index(N), s(N);
  
  double d,dp=0.0,weight,uw,eps=1e-9;
  
  for(int i=0;i<N;i++){p[i]=prob[i];}
  
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
      quicksort(index, dists, weights, 0, N-i-2);
      
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
