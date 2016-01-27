\name{BalancedSampling-package}
\alias{BalancedSampling-package}
\alias{BalancedSampling}
\docType{package}
\title{
Balanced and Spatially Balanced Sampling
}
\description{
Select balanced and spatially balanced probability samples in multi-dimensional spaces with any prescribed inclusion probabilities. It contains fast (C++ via Rcpp) implementations of the included sampling methods. The local pivotal method and spatially correlated Poisson sampling (for spatially balanced sampling) are included. Also the cube method (for balanced sampling) and the local cube method (for doubly balanced sampling) are included. 
}
\author{
Anton Grafström, Jonathan Lisic

Maintainer: Anton Grafström <anton.grafstrom@gmail.com>

Webpage: http://www.antongrafstrom.se/balancedsampling
}
\references{
Deville, J. C. and  Tillé, Y. (2004). Efficient balanced sampling: the cube method. Biometrika, 91(4), 893-912.

Deville, J.-C. and  Tillé, Y. (1998). Unequal probability sampling without replacement through a splitting method. Biometrika 85, 89-101.

Grafström, A. (2012). Spatially correlated Poisson sampling. Journal of Statistical Planning and Inference, 142(1), 139-147.

Grafström, A. and Lundström, N.L.P. (2013). Why well spread probability samples are balanced. Open Journal of Statistics, 3(1).

Grafström, A. and Schelin, L. (2014). How to select representative samples. Scandinavian Journal of Statistics.

Grafström, A., Lundström, N.L.P. and Schelin, L. (2012). Spatially balanced sampling through the Pivotal method. Biometrics 68(2), 514-520.

Grafström, A. and Tillé, Y. (2013). Doubly balanced spatial sampling with spreading and restitution of auxiliary totals. Environmetrics, 24(2), 120-131.
}
\keyword{ sampling }
\examples{
# *********************************************************
# check inclusion probabilities
# *********************************************************
set.seed(1234567);
p = c(0.2, 0.25, 0.35, 0.4, 0.5, 0.5, 0.55, 0.65, 0.7, 0.9);
N = length(p);
X = cbind(runif(N),runif(N));
p1 = p2 = p3 = p4 = rep(0,N);
nrs = 1000; # increase for more precision 
for(i in 1:nrs){
  # lpm1
  s = lpm1(p,X);
  p1[s]=p1[s]+1;
  
  # lpm2
  s = lpm2(p,X);
  p2[s]=p2[s]+1;
  
  # scps
  s = scps(p,X);
  p3[s]=p3[s]+1;
  
  # lcube
  s = lcube(p,X,cbind(p));
  p4[s]=p4[s]+1; 
}
print(p);
print(p1/nrs);
print(p2/nrs);
print(p3/nrs);
print(p4/nrs);

# *********************************************************
# check spatial balance
# *********************************************************
set.seed(1234567);
N = 500;
n = 70;
p = rep(n/N,N);
X = cbind(runif(N),runif(N));
nrs = 10; # increase for more precision 
b1 = b2 = b3 = b4 = b5 = rep(0,nrs);

for(i in 1:nrs){
  # lpm1
  s = lpm1(p,X);
  b1[i] = sb(p,X,s);
  
  # lpm2
  s = lpm2(p,X);
  b2[i] = sb(p,X,s);
  
  # scps
  s = scps(p,X);
  b3[i] = sb(p,X,s);
  
  # lcube
  s = lcube(p,X,cbind(p));
  b4[i] = sb(p,X,s);  
  
  # srs
  s = sample(N,n);
  b5[i] = sb(p,X,s);
}
print(mean(b1));
print(mean(b2));
print(mean(b3));
print(mean(b4));
print(mean(b5));

# *********************************************************
# stratification
# *********************************************************
set.seed(1234567);
N = 10;
n = 4;
p = rep(n/N,N);
stratum1 = c(1,1,1,1,1,0,0,0,0,0); # stratum 1 indicator
stratum2 = c(0,0,0,0,0,1,1,1,1,1); # stratum 2 indicator
stratum3 = c(0,0,1,1,1,1,1,0,0,0); # overlapping 1 and 2
s = lpm1(p,cbind(stratum1,stratum2,stratum3));

# *********************************************************
# plot spatially balanced sample
# *********************************************************
set.seed(1234567);
N = 1000; # population size
n = 100; # sample size
p = rep(n/N,N); # inclusion probabilities
X = cbind(runif(N),runif(N)); # matrix of auxiliary variables
s = lpm1(p,X); # select sample 
plot(X[,1],X[,2]); # plot population
points(X[s,1],X[s,2], pch=19); # plot sample

# *********************************************************
# check cpu time (for simulation)
# *********************************************************
set.seed(1234567);
N = 2000;
n = 100;
X = cbind(runif(N),runif(N));
p = rep(n/N,N);
system.time(for(i in 1:10){lpm1(p,X)});
system.time(for(i in 1:10){lpm2(p,X)});

}
