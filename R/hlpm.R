# Hierarcical local pivotal method (hlpm)
# selects an lpm sample and splits the sample
# into subsamples of given sizes using lpm 

# **********************************************
# Author: Anton Grafström
# Last edit: 2018-08-14 
# Licence: GPL (>=2)
# **********************************************

hlpm = function(p,X,sizes){
  
  # simple parameter check
  if(round(sum(p))!=sum(sizes)){
    stop("Probabilities do not match sizes!");
  }
  
  numberOfSamples = length(sizes);  # number of samples
  N = length(p);                    # population size
  S = SamplingBigData::lpm2_kdtree(p,X);             # initial sample
  XS = as.matrix(X[S,]);            # auxiliary for initial sample
  index = 1:length(S);              # index
  done = rep(0,length(S));          # indicator of subsampled
  sampleNumber = rep(0,length(S));  # sample number of subsampled
  
  for(t in 1:numberOfSamples){
    
    # index of remaining units
    subindex = as.vector(index[done==0]);
    
    # auxiliary of remaining units
    Xt = XS[subindex,];
    if(length(subindex)==1){
      Xt = rbind(Xt); 
    }
    
    # number of remaining units
    Nt = length(subindex); 
    
    # sample index of subsampled units
    St = subindex[SamplingBigData::lpm2_kdtree(rep(sizes[t]/Nt,Nt),Xt)]; 
    
    # set indicator for sampled units
    done[St] = 1; 
    
    # store sample number
    sampleNumber[St] = t; 
  }
  
  # return population index of sampled units and 
  # their sample number
  list(S = S,sampleNumber = sampleNumber);
} 
