#include <Rcpp.h>
using namespace Rcpp;

//**********************************************
// Author: Anton Grafström
// Last edit: 2014-03-18 
// Licence: GPL (>=2)
//**********************************************

// reduced row echelon form
void rref(NumericMatrix& M){
  int lead = 0;
  int rowCount = M.nrow();
  int columnCount = M.ncol();
  int r,i,k;
  double temp; 
  for(r=0; r<rowCount; r++){
    if(columnCount<=lead){return;}
    i = r;
    while( M(i,lead) == 0 ){
      i = i + 1;
      if(i == rowCount){
        i = r;
        lead = lead + 1;
        if(columnCount == lead){return;}
      }
    }
    // swap rows i and r
    for(k=0;k<columnCount;k++){
      temp = M(i,k);
      M(i,k) = M(r,k);
      M(r,k) = temp;
    }  
    // If M[r, lead] is not 0 divide row r by M[r, lead]
    if( M(r,lead) != 0 ){
      temp = M(r,lead);
      for(k=0;k<columnCount;k++){
        M(r,k) = M(r,k)/temp;
      }
    }
    for(i=0;i<rowCount;i++){
      if( i != r ){
        // Subtract M[i, lead] multiplied by row r from row i
        temp = M(i,lead);
        for(k=0;k<columnCount;k++){
          M(i,k) = M(i,k) - temp * M(r,k);
        }
      }
    }
    lead = lead + 1;
  } 
  return;
}


// one step fast flight cube
NumericVector onestepfastflightcube(NumericVector prob, NumericMatrix Bm){
	int ncol = Bm.ncol();
	int nrow = Bm.nrow();
	int i, j;
  NumericVector u(ncol,0.0);
  IntegerVector uset(ncol,0);
	double lambda1 = 1e+200;
	double lambda2 = 1e+200;
	double lambda, eps = 1e-9;
  int lead;
  double v, free = -1;

	// find nonzero vector u in Ker B (null space of B, i.e. Bu = 0)
	// with both positive and negative values
	// find reduced row echelon form of B
	rref(Bm);
  
	for(i=(nrow-1);i>=0;i--){
		// find lead (first nonzero entry on row) if exists
		// if no lead, i.e lead = ncol, do nothing
		// if lead, the variables after are either set or free
		// free variables are alternately set to 1 or -1 
		lead = 0;
		for(j=0;j<ncol;j++){if(Bm(i,j)==0){lead++;}else{break;}}
		// lead found
		if(lead<ncol){
			v = 0;
			for(j=lead+1;j<ncol;j++){
				if( uset[j] == 0 ){
          uset[j] = 1;
					free *= -1;
					u[j] = free;
				}
				v -= u[j]*Bm(i,j);		
			}
			u[lead] = v/Bm(i,lead);
      uset[lead] = 1;
		}
	}
  // unset u[i] are free and are set to 1 or -1, can only exist at beginning
  for(i=0;i<ncol;i++){
		if( uset[i] == 0 ){
			free *= -1;
			u[i] = free;
		}else{break;}	
	}
	// find lambda1 and lambda2
	for(i=0;i<ncol;i++){
		if(u[i]>0){
			lambda1 = std::min(lambda1,(1-prob[i])/u[i]);
			lambda2 = std::min(lambda2,prob[i]/u[i]);
		}
		if(u[i]<0){
			lambda1 = std::min(lambda1,-prob[i]/u[i]);
			lambda2 = std::min(lambda2,(prob[i]-1)/u[i]);					
		}
	}
	// random choice of p+lambda1*u and p-lambda2*u
	if(runif(1)[0]<lambda2/(lambda1+lambda2)){
		lambda = lambda1;
	}else{
		lambda = -lambda2;
	}
	// update prob
	for(i=0;i<ncol;i++){
			prob[i] = prob[i] + lambda * u[i];
			if(prob[i] < eps){ prob[i] = 0; }
			if(prob[i] > 1-eps){ prob[i] = 1; }
	}
	return prob;	
}

// [[Rcpp::export]]
IntegerVector cube(NumericVector prob, NumericMatrix Xbal){
  // landing by supression of balancing variables, 
	// starting from the column with largest index.	  
  int N = prob.size();
  int nrow = Xbal.nrow();
  int naux = Xbal.ncol();
 
  IntegerVector index(N);
  NumericVector p(N);
  int i,j,k,howmany,first;
  for(i=0;i<N;i++){index[i]=i; p[i]=prob[i];}  
  //IntegerVector index_small;
  //NumericVector p_small, dists;
  //NumericMatrix B;
  double eps = 1e-9;
  int done = 0, tempInt, howlong;
  // randomize order of index list
  NumericVector rnd = runif(N);
	for(i=0;i<N;i++){
			k = i + floor(rnd[i] * (N-i));
			tempInt = index[i];
			index[i] = index[k];
			index[k] = tempInt;
	}
  // put finished units at beginning of list
  for(i=done;i<N;i++){
    if( p[index[i]]<eps || p[index[i]]>1-eps ){
      tempInt = index[done];
      index[done] = index[i];
      index[i] = tempInt;
      done = done + 1;
    }
  } 
  // remaining are index from done to N-1
  while( done < N ){
    // find cluster of size howmany
    howmany = std::min(naux+1,N-done);
    if( howmany > 1 ){
      NumericVector p_small(howmany);
      NumericVector dists(howmany,1e+200);
      IntegerVector index_small(howmany);
      NumericMatrix B(howmany-1,howmany);       
      for(i=0;i<howmany;i++){
        index_small[i] = index[done+i];
				for(j=0;j<howmany-1;j++){
					B(j,i) = Xbal(index_small[i],j)/prob[index_small[i]];
				}
				p_small[i] = p[index_small[i]];
			}
      p_small = onestepfastflightcube(p_small,B);
      // update prob
			for(i=0;i<howmany;i++){
				p[index_small[i]] = p_small[i];
			}
      // update done and index
      howlong = done + howmany;
      for(i=done;i<howlong;i++){
        if( p[index[i]]<eps || p[index[i]]>1-eps ){
          tempInt = index[done];
          index[done] = index[i];
          index[i] = tempInt;
          done = done + 1;
        }
      } 
      
    }else{
      // max one unit left
			if(runif(1)[0]<p[index[done]]){p[index[done]]=1;}else{p[index[done]]=0;}
      done = N;
    }
  }
  // construct sample from prob;
  int n = round(sum(p)), count = 0;
  IntegerVector s(n);
  for(i=0;i<N;i++){
    if(p[index[i]]>1-eps){
      // switch to 1-based index for sample
      s[count] = index[i] + 1;
      count = count + 1;
    }
  }
	return s;
}
