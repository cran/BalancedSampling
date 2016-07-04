lpm2_kdtree <- function(
  prob,
  x,                      
  m=40,
  algorithm="kdtree",
  maxCheck=4,
  termDist=.1
) {

  if(!is.matrix(x)) x <- as.matrix(x)

  n <- NROW(x)
  K <- length(x) / n
  m <- min( m, n)
 

  if( algorithm=="kdtree" ) {
    algorithm <- 0 
  } else if( algorithm=="kdtree-count" ) {
    if(maxCheck < 1)
      stop(print("kdtree-count requires a positive valued integer maxCheck parameter.")) 
    algorithm <- 1
  } else if( algorithm=="kdtree-dist" ) {
    if(termDist <= 0)
      stop(print("kdtree-dist requires a positive valued floating point termDist parameter.")) 
    algorithm <- 2 
  } else {
      stop(print("Specified algorithm does not exist.")) 
  }

  # sanity check to avoid segfaults in R_lpm3 
  if( length(prob) != n ) stop(print("Length of probability vector does not match the number of rows in x."))
  
  # send our data to the C program
  r.result <- .C("R_lpm3",
                 as.double( t(x) ),                 # data we query
                 as.double( prob ),                 # probability vector
                 as.integer( n ),                   # length of prob and nrow of x 
                 as.integer( K ),                   # number of columns of x 
                 as.integer( m ),                   # max leaves per node
                 as.integer( algorithm ),           # algorithm to use 
                 as.integer( maxCheck ),            # number of leaves to check
                 as.double( termDist )              # terminal distance 
  )
  
  return( (1:n)[ r.result[[2]] > .5 ] )
}
