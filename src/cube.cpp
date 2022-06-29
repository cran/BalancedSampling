#include <Rcpp.h>
using namespace Rcpp;

//**********************************************
// Author: Anton Grafström
// Last edit: 2020-08-14
// Licence: GPL (>=2)
//**********************************************

// "import" print for error checking
Function print("print");

/* Is identity function
 * 
 * This function test if the matrix M is identity or not.
 * The matrix could be a n x m matrix.
 * 
 */
bool isEye(NumericMatrix M){
  bool out = true;
  int p = M.ncol();
  int N = M.nrow();
  for(int i = 0;i< N; i++){
    for(int j = 0; j < p; j++){
      if(i == j && M(i,j) != 1){
        out = false;
        break;
      }
      if(i != j && M(i,j) != 0){
        out = false;
        break;
      }
    }
  }
  return out;
}

// reduced row echelon form
void rref(NumericMatrix& M){
  int lead = 0;
  int rowCount = M.nrow();
  int columnCount = M.ncol();
  double eps = 1e-11;
  int r,i,k;
  double temp; 
  for(r=0; r<rowCount; r++){
    if(columnCount<=lead){return;}
    i = r;
    while( std::max(M(i,lead),-M(i,lead)) < eps ){
      M(i,lead) = 0.0;
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
    // If M(r, lead) is not 0 divide row r by M(r, lead)
    if( M(r,lead) != 0.0 ){
      temp = M(r,lead);
      for(k=0;k<lead;k++){M(r,k) = 0.0;}
      for(k=lead;k<columnCount;k++){
        M(r,k) = M(r,k)/temp;
      }
    }
    for(i=0;i<rowCount;i++){
      if( i != r ){
        // Subtract M(i, lead) multiplied by row r from row i
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
	double la1 = 1e+200;
	double la2 = 1e+200;
	double la, eps = 1e-9;
  int lead;
  double v, free = -1.0;
	// find nonzero vector u in Ker B (null space of B, i.e. Bu = 0)
	// with both positive and negative values
	// find reduced row echelon form of B
	// rref(Bm);
	for(i=(nrow-1);i>=0;i--){
		// find lead (first nonzero entry on row) if exists
		// if no lead, i.e lead = ncol, do nothing
		// if lead, the variables after are either set or free
		// free variables are alternately set to 1 or -1 
		lead = 0;
		for(j=0;j<ncol;j++){if(Bm(i,j)==0.0){lead++;}else{break;}}
		// lead found
		if(lead<ncol){
			v = 0.0;
			for(j=lead+1;j<ncol;j++){
				if( uset[j] == 0 ){
          uset[j] = 1;
					free *= -1.0;
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
			free *= -1.0;
			u[i] = free;
		}else{break;}	
	}
	// find lambda1 and lambda2
	for(i=0;i<ncol;i++){
		if(u[i]>0){
			la1 = std::min(la1,(1-prob[i])/u[i]);
			la2 = std::min(la2,prob[i]/u[i]);
		}
		if(u[i]<0){
			la1 = std::min(la1,-prob[i]/u[i]);
			la2 = std::min(la2,(prob[i]-1)/u[i]);					
		}
	}
	// random choice of p+lambda1*u and p-lambda2*u
	if(runif(1)[0]<la2/(la1+la2)){
		la = la1;
	}else{
		la = -la2;
	}
	// update prob
	for(i=0;i<ncol;i++){
			prob[i] = prob[i] + la * u[i];
			if(prob[i] < eps){ prob[i] = 0; }
			if(prob[i] > 1-eps){ prob[i] = 1; }
	}
	return prob;	
}


// [[Rcpp::export]]
IntegerVector cube(NumericVector prob, NumericMatrix Xbal){
  // landing by dropping balancing variables, 
	// starting from the column with largest index.	  
  int N = prob.size();
  int naux = Xbal.ncol(); 
  IntegerVector index(N);
  NumericVector p(N);
  int i,j,k,howmany;
  for(i=0;i<N;i++){index[i]=i; p[i]=prob[i];}  

  double eps = 1e-12;
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
      rref(B); 
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

// [[Rcpp::export]]
IntegerVector lcube(NumericVector prob, NumericMatrix Xspread, NumericMatrix Xbal){
  // landing by dropping balancing variables, 
  // starting from the column with largest index.    
  int N = prob.size();
  int naux = Xbal.ncol();
  int nspread = Xspread.ncol();
 
  IntegerVector index(N);
  NumericVector p(N);
  int i,j,k,howmany,first;
  for(i=0;i<N;i++){index[i]=i; p[i]=prob[i];}  
  double eps = 1e-12;
  double d;
  // put finished units at beginning of list
  int done = 0, tempInt;
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
    // choose one index randomly
    first = floor(done + runif(1)[0] * (N-done) );
    // find cluster of size howmany
    howmany = std::min(naux+1,N-done);
    if( howmany > 1 ){
      NumericVector p_small(howmany);
      NumericVector dists(howmany,1e+200);
      IntegerVector index_small(howmany);
      NumericMatrix B(howmany-1,howmany);     
      
      for(i=done;i<N;i++){
        d = 0;
        for(j=0;j<nspread;j++){
          d = d + (Xspread(index[i],j)-Xspread(index[first],j))*(Xspread(index[i],j)-Xspread(index[first],j));
        }
				// check if candidate
				if(d<dists[howmany-1]){
					for(j=0;j<howmany;j++){
						if(d<dists[j]){
              // insert into index_small, dists and shift remaining 
              // one step to the right
              for(k=howmany-1;k>j;k--){
                index_small[k] = index_small[k-1];
                dists[k] = dists[k-1];
              }
              dists[j] = d;
              index_small[j] = index[i];
							break;	
						}
					}
				}
			}
      // cluster found, build p_small and B
      for(i=0;i<howmany;i++){
				for(j=0;j<howmany-1;j++){
					B(j,i) = Xbal(index_small[i],j)/prob[index_small[i]];
				}
				p_small[i] = p[index_small[i]];
			}
      rref(B); 
      p_small = onestepfastflightcube(p_small,B);
      // update prob
			for(i=0;i<howmany;i++){
				p[index_small[i]] = p_small[i];
			}
      // update done and index
      for(i=done;i<N;i++){
        if( p[index[i]]<eps || p[index[i]]>1-eps ){
          tempInt = index[done];
          index[done] = index[i];
          index[i] = tempInt;
          done = done + 1;
        }
      }  
    }else{
      // max one unit left
			if(runif(1)[0]<p[index[first]]){p[index[first]]=1;}else{p[index[first]]=0;}
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

// [[Rcpp::export]]
NumericVector flightphase(NumericVector prob, NumericMatrix Xbal){
  int N = prob.size();
  int naux = Xbal.ncol();
  
  bool eye = false; // boolean, if matrix B is identity.

  IntegerVector index(N);
  NumericVector p(N);
  int i,j,k,howmany;
  for(i=0;i<N;i++){index[i]=i; p[i]=prob[i];}  
  double eps = 1e-12;
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
    // stop if there are less than naux units left
    // if(howmany <= naux){done=N; break;}
 
    if( howmany > 1 ){
      NumericVector p_small(howmany);
      NumericVector dists(howmany,1e+200);
      IntegerVector index_small(howmany);
      // NumericMatrix B(howmany-1,howmany); 
      NumericMatrix B(naux,howmany); // changed      
      for(i=0;i<howmany;i++){
        index_small[i] = index[done+i];
				//for(j=0;j<howmany-1;j++){
				for(j=0; j < naux;j++){
					B(j,i) = Xbal(index_small[i],j)/prob[index_small[i]];
				}
				p_small[i] = p[index_small[i]];
			}
      //p_small = onestepfastflightcube(p_small,B);
      rref(B); 
      if(howmany < naux + 1){
        eye = isEye(B);
        if(eye == true){
          break;
        }
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
  // round
  for(i=0;i<N;i++){
    if( p[index[i]] > 1-eps  ){
      p[index[i]] = 1;
    }
    if( p[index[i]] < eps  ){
      p[index[i]] = 0;
    }
  }
	return p;
}

// [[Rcpp::export]]
IntegerVector landingphase(NumericVector prob, NumericVector probflight, NumericMatrix Xbal){
  // landing by dropping balancing variables, 
  // starting from the column with largest index.	  
  int N = prob.size();
  int naux = Xbal.ncol();
 
  IntegerVector index(N);
  NumericVector p(N);
  int i,j,k,howmany;
  for(i=0;i<N;i++){index[i]=i; p[i]=probflight[i];}  
  double eps = 1e-12;
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
      rref(B); 
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
  // round
  for(i=0;i<N;i++){
    if( p[index[i]] > 1-eps  ){
      p[index[i]] = 1;
    }
    if( p[index[i]] < eps  ){
      p[index[i]] = 0;
    }
  }
  return as<IntegerVector>(p);
}

// [[Rcpp::export]]
NumericVector lcubeflightphase(NumericVector prob, NumericMatrix Xspread, NumericMatrix Xbal){    
  int N = prob.size();
  int naux = Xbal.ncol();
  int nspread = Xspread.ncol();
 
  IntegerVector index(N);
  NumericVector p(N);
  int i,j,k,howmany,first;
  for(i=0;i<N;i++){index[i]=i; p[i]=prob[i];}  
  double eps = 1e-11;
  double d;
  // put finished units at beginning of list
  int done = 0, tempInt;
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
    // choose one index randomly
    first = floor(done + runif(1)[0] * (N-done) );
    // find cluster of size howmany
    howmany = std::min(naux+1,N-done);
    // stop if there are less than naux units left
    if(howmany < naux + 1){done=N; break;}
    
    if( howmany > 1 ){
      NumericVector p_small(howmany);
      NumericVector dists(howmany,1e+200);
      IntegerVector index_small(howmany);
      NumericMatrix B(howmany-1,howmany);     
      
      for(i=done;i<N;i++){
        d = 0;
        for(j=0;j<nspread;j++){
          d = d + (Xspread(index[i],j)-Xspread(index[first],j))*(Xspread(index[i],j)-Xspread(index[first],j));
        }
  			// check if candidate
				if(d<dists[howmany-1]){
					for(j=0;j<howmany;j++){
						if(d<dists[j]){
              // insert into index_small, dists and shift remaining 
              // one step to the right
              for(k=howmany-1;k>j;k--){
                index_small[k] = index_small[k-1];
                dists[k] = dists[k-1];
              }
              dists[j] = d;
              index_small[j] = index[i];
							break;	
						}
					}
				}
			}
      // cluster found, build p_small and B
      for(i=0;i<howmany;i++){
				for(j=0;j<howmany-1;j++){
					B(j,i) = Xbal(index_small[i],j)/prob[index_small[i]];
				}
				p_small[i] = p[index_small[i]];
			}
      rref(B); 
      p_small = onestepfastflightcube(p_small,B);
      // update prob
			for(i=0;i<howmany;i++){
				p[index_small[i]] = p_small[i];
			}
      // update done and index
      for(i=done;i<N;i++){
        if( p[index[i]]<eps || p[index[i]]>1-eps ){
          tempInt = index[done];
          index[done] = index[i];
          index[i] = tempInt;
          done = done + 1;
        }
      }
    }else{
      // max one unit left
			if(runif(1)[0]<p[index[first]]){p[index[first]]=1;}else{p[index[first]]=0;}
      done = N;
    }
  }
  // round
  for(i=0;i<N;i++){
    if( p[index[i]] > 1-eps  ){
      p[index[i]] = 1;
    }
    if( p[index[i]] < eps  ){
      p[index[i]] = 0;
    }
  }
  return p;
}

// [[Rcpp::export]]
IntegerVector lcubelandingphase(NumericVector prob, NumericVector probflight, NumericMatrix Xspread, NumericMatrix Xbal){
  // landing by dropping balancing variables, 
  // starting from the column with largest index.    
  int N = prob.size();
  int naux = Xbal.ncol();
  int nspread = Xspread.ncol();
 
  IntegerVector index(N);
  NumericVector p(N);
  int i,j,k,howmany,first;
  for(i=0;i<N;i++){index[i]=i; p[i]=probflight[i];}  
  double eps = 1e-12;
  double d;
  // put finished units at beginning of list
  int done = 0, tempInt;
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
    // choose one index randomly
    first = floor(done + runif(1)[0] * (N-done) );
    // find cluster of size howmany
    howmany = std::min(naux+1,N-done);
    if( howmany > 1 ){
      NumericVector p_small(howmany);
      NumericVector dists(howmany,1e+200);
      IntegerVector index_small(howmany);
      NumericMatrix B(howmany-1,howmany);     
      
      for(i=done;i<N;i++){
        d = 0;
        for(j=0;j<nspread;j++){
          d = d + (Xspread(index[i],j)-Xspread(index[first],j))*(Xspread(index[i],j)-Xspread(index[first],j));
        }
    		// check if candidate
				if(d<dists[howmany-1]){
					for(j=0;j<howmany;j++){
						if(d<dists[j]){
              // insert into index_small, dists and shift remaining 
              // one step to the right
              for(k=howmany-1;k>j;k--){
                index_small[k] = index_small[k-1];
                dists[k] = dists[k-1];
              }
              dists[j] = d;
              index_small[j] = index[i];
							break;	
						}
					}
				}
			}
      // cluster found, build p_small and B
      for(i=0;i<howmany;i++){
				for(j=0;j<howmany-1;j++){
          B(j,i) = Xbal(index_small[i],j)/prob[index_small[i]];
				}
				p_small[i] = p[index_small[i]];
			}
      rref(B); 
      p_small = onestepfastflightcube(p_small,B);
      // update prob
			for(i=0;i<howmany;i++){
				p[index_small[i]] = p_small[i];
			}
      // update done and index
      for(i=done;i<N;i++){
        if( p[index[i]]<eps || p[index[i]]>1-eps ){
          tempInt = index[done];
          index[done] = index[i];
          index[i] = tempInt;
          done = done + 1;
        }
      }    
    }else{
      // max one unit left
			if(runif(1)[0]<p[index[first]]){p[index[first]]=1;}else{p[index[first]]=0;}
      done = N;
    }
  }
  // round
  for(i=0;i<N;i++){
    if( p[index[i]] > 1-eps  ){
      p[index[i]] = 1;
    }
    if( p[index[i]] < eps  ){
      p[index[i]] = 0;
    }
  }
  return as<IntegerVector>(p);
}

// [[Rcpp::export]]
IntegerVector cubestratified(NumericVector prob, NumericMatrix Xbal, IntegerVector integerStrata){
  
  int naux = Xbal.ncol();
  int N = prob.size();
  IntegerVector uniqueStrata = as<IntegerVector>(unique(integerStrata));
  int H = uniqueStrata.size();
  NumericVector p_star(N);
  NumericMatrix x_star(N,H+naux);
  IntegerVector counts(H);
  List strataIndex(H);  
  int i,j,h,nh,count;

  // build strata index vectors
  for(h=0;h<H;h++){
    count=0;
    for(i=0;i<N;i++){
      if(integerStrata[i]==uniqueStrata[h]){
        count++;
      }
    }
    counts[h]=count;
    IntegerVector ih(count);
    count=0;
    for(i=0;i<N;i++){
      if(integerStrata[i]==uniqueStrata[h]){
        ih[count] = i;
        count++;
      }
    }
    // TODO: randomize order of ih here
    strataIndex[h]=ih;
  }
  // 1. run flight in each stratum
  for(h=0;h<H;h++){
    IntegerVector ih = strataIndex[h];
    nh = counts[h];
    NumericVector ph(nh);
    NumericMatrix xh(nh,naux);
    for(i=0;i<nh;i++){
        xh(i,_) = Xbal(ih[i],_);
        ph[i] = prob[ih[i]];
    }
    NumericVector sh = flightphase(ph,xh);
    for(i=0;i<nh;i++){
      p_star[ih[i]] = sh[i];
      x_star(ih[i],h) = sh[i];
      for(j=H;j<H+naux;j++){
        x_star(ih[i],j) = Xbal(ih[i],j-H)*p_star[ih[i]]/prob[ih[i]];
      }
    }
  }
  // 2. run flight on full population
  NumericVector p_star_star = flightphase(p_star,x_star);
  // 3. run landing in each stratum
  for(h=0;h<H;h++){
    IntegerVector ih = strataIndex[h];
    nh = counts[h];
    NumericVector ph_flight(nh);
    NumericMatrix xh(nh,naux);
    IntegerVector sh(nh);
    for(i=0;i<nh;i++){
        ph_flight[i] = p_star_star[ih[i]];
        xh(i,_) = Xbal(ih[i],_)*(p_star_star[ih[i]]/prob[ih[i]]);
    }
    sh = landingphase(ph_flight,ph_flight,xh);
    for(i=0;i<nh;i++){
      p_star_star[ih[i]] = sh[i];
    }
  }
  return as<IntegerVector>(p_star_star);
}


// [[Rcpp::export]]
IntegerVector lcubestratified(NumericVector prob, NumericMatrix Xspread, NumericMatrix Xbal, IntegerVector integerStrata){
  
  int naux = Xbal.ncol();
  int nspread = Xspread.ncol();
  int N = prob.size();
  IntegerVector uniqueStrata = as<IntegerVector>(unique(integerStrata));
  int H = uniqueStrata.size();
  NumericVector p_star(N);
  NumericMatrix x_star(N,H+naux);
  IntegerVector counts(H);
  List strataIndex(H);  
  int i,j,h,nh,count;

  // build strata index vectors
  for(h=0;h<H;h++){
    count=0;
    for(i=0;i<N;i++){
      if(integerStrata[i]==uniqueStrata[h]){
        count++;
      }
    }
    counts[h]=count;
    IntegerVector ih(count);
    count=0;
    for(i=0;i<N;i++){
      if(integerStrata[i]==uniqueStrata[h]){
        ih[count] = i;
        count++;
      }
    }
    strataIndex[h]=ih;
  }
  // 1. run flight in each stratum
  for(h=0;h<H;h++){
    IntegerVector ih = strataIndex[h];
    nh = counts[h];
    NumericVector ph(nh);
    NumericMatrix xh(nh,naux);
    NumericMatrix xhspread(nh,nspread);
    for(i=0;i<nh;i++){
        xhspread(i,_) = Xspread(ih[i],_);
        xh(i,_) = Xbal(ih[i],_); 
        ph[i] = prob[ih[i]];
    }
    NumericVector sh = lcubeflightphase(ph,xhspread,xh);
    for(i=0;i<nh;i++){
      p_star[ih[i]] = sh[i];
      x_star(ih[i],h) = sh[i];
      for(j=H;j<H+naux;j++){
        x_star(ih[i],j) = Xbal(ih[i],j-H)*p_star[ih[i]]/prob[ih[i]];
      }
    }
  } 
  // 2. run flight on full population
  NumericVector p_star_star = lcubeflightphase(p_star,Xspread,x_star);
  // 3. run landing in each stratum
  for(h=0;h<H;h++){
    IntegerVector ih = strataIndex[h];
    nh = counts[h];
    NumericVector ph_flight(nh);
    NumericMatrix xh(nh,naux);
    IntegerVector sh(nh);
    for(i=0;i<nh;i++){
      xh(i,_) = Xbal(ih[i],_)*(p_star_star[ih[i]]/prob[ih[i]]);
      ph_flight[i] = p_star_star[ih[i]];
    }
    sh = landingphase(ph_flight,ph_flight,xh);
    for(i=0;i<nh;i++){
      p_star_star[ih[i]] = sh[i];
    }
  }
  return as<IntegerVector>(p_star_star);
}
