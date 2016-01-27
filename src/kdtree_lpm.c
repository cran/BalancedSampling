/* Copyright (c) 2015  Jonathan Lisic 
 * Last edit: 15/07/08 - 11:45:28
 * License: BSD 3-Clause
 */  

#include "kdtree_lpm.h"


/* printf fixing */
#ifdef CLI 
 #define PRINTF printf
#else 
 #include "R.h"
 #include "Rmath.h"
 #define PRINTF Rprintf
#endif

/* function to create a new Tree */
rootNodePtr createTree( size_t K, size_t leafSize, size_t n, double * data ) {
  
  rootNodePtr y = malloc( sizeof( rootNode ) );
  
  y->pointerIndex = calloc( n, sizeof( size_t * ) );

  /* setup root node */
  y->K = K;
  y->leafSize = leafSize;
  y->root = NULL;
  y->data = data;
  y->n = n;

  return(y);
}


/* function to create a new Tree */
void deleteTree( rootNodePtr r ) {
 
  if(r->pointerIndex != NULL) free( r->pointerIndex ); 
  r->pointerIndex = NULL;
 
  deleteNode( r, r->root ); 

  free(r);
  return;
}


/* add an element to the tree */
void buildIndex( rootNodePtr r, nodePtr c, nodePtr p ) {
 
  size_t i,k; 
  size_t K, n, dim;
  double * data;


  /* for the root node */
  if( c == NULL) {
    c = createNode( r, NULL);
    r->root = c;
    c->indexUsed = r->n;

    K = r->K;
    data = r->data;
    dim =  c->dim; 
    n = r->n;

    // setup the index, and also ensure that we have finite bounds
    for(i = 0; i < n; i++) {
      c->index[i] = i;
      for( k = 0; k < K; k++) {
        if( data[ i * K + k ] < c->min[k] ) c->min[k] = data[i * K + k];  
        if( data[ i * K + k ] > c->max[k] ) c->max[k] = data[i * K + k];  
      }
    }
  }
  nodePtr * children = NULL; // place holder
  
  /* do we have too many points? */
  if( c->indexUsed < r->leafSize ) {

    /* save the final pointer locations */
    for( i = 0; i < c->indexUsed; i++) 
      r->pointerIndex[ c->index[i] ] = &( c->index[i] );

    return;
  } 

  /* if we have too many points 
   * create children
   * figure out our new dim
   * figure out what the split is
   * move current contents to new children
   */
  
//PRINTF(" so you have decided to have children\n");
  children = createChildren( r, c, p);
//PRINTF(" success!, they have left the house off on their own!\n");
  c->left  = children[0];
  c->right = children[1];
  free(children);
  
  buildIndex( r, c->left, c );  
  buildIndex( r, c->right, c );   

  return;
}




/* double comparison function for qsort */
int comp ( const void * a, const void * b ) {

  double aDbl = *((double*) a);
  double bDbl = *((double*) b);

  if( aDbl > bDbl) return 1;
  if( aDbl < bDbl) return -1;

  return 0;
}


/* double * comparison function for qsort */
int compareDoublePtr ( const void * aPtr, const void * bPtr ) {
  double a = **(double **) aPtr;
  double b = **(double **) bPtr; 
  if (a < b) return -1;
  if (a > b) return  1; 
  return 0; 
}


/* function to create a simple node */
nodePtr createNode( rootNodePtr r, nodePtr p ) {

  size_t i;

  nodePtr c = malloc( sizeof( node ) );

  if( p != NULL) {
    c->index = calloc( p->indexUsed/2 +1, sizeof(size_t));
  } else {
    c->index = calloc( r->n, sizeof(size_t));
  }
  c->head  = p;
  c->left  = NULL;
  c->right = NULL;
  c->indexUsed = 0;
  c->split = 0;
  c->dim = 0;

  c->min = calloc( r->K, sizeof(double));
  c->max = calloc( r->K, sizeof(double));
  for( i = 0; i < r->K; i++) {
    (c->max)[i] = -INFINITY;
    (c->min)[i] = INFINITY;
  }
  return(c);
}


/* function to delete a node */
void deleteNode( rootNodePtr r, nodePtr c ) {

  if( c == NULL ) return;

  if( c->index != NULL) free(c->index);
  if( c->min != NULL) free(c->min);
  if( c->max != NULL) free(c->max);

  // to the left
  if( c->left != NULL ) {
    deleteNode( r, c->left );
    c->left = NULL;
  }

  // to the right
  if( c->right != NULL ) {
    deleteNode( r, c->right );
    c->right = NULL;
  }

  free(c);
  c = NULL;
  return;
}


/* split and create children from a *full* node */
nodePtr * createChildren( rootNodePtr r, nodePtr c, nodePtr p) {

  double * x = NULL;
  double ** xPtr = NULL;
  size_t i;
  size_t splitIndex;
  nodePtr * children = NULL;

  if( p != NULL) {
    c->dim = (p->dim + 1) % r->K; //increment mod K 
  } else {
    c->dim = 0;
  }



  /***********************************************************************************************************/

  /* get the median */
  x =  calloc( c->indexUsed, sizeof(double) ); //allocate some temporary space for finding the median
  xPtr =  calloc( c->indexUsed, sizeof( double * ) ); //allocate some temporary space for finding the median

  /* use quick sort to find the median */
  for( i = 0; i < c->indexUsed; i++) {
    x[i] = r->data[ c->index[i] * r->K + c->dim];
    xPtr[i] = &x[i];
  }  

  qsort( xPtr, c->indexUsed, sizeof( double * ), compareDoublePtr);  

  /* update min and max */
  c->min[c->dim] = *(xPtr[0]);
  c->max[c->dim] = *(xPtr[c->indexUsed - 1]);

  // get split
  splitIndex = c->indexUsed/2;
  if( c->indexUsed % 2 ) {  // is not even
    //splitIndex = c->indexUsed/2;
    c->split = *xPtr[splitIndex];
  } else {
    c->split = ( *xPtr[ splitIndex ] + *xPtr[ splitIndex - 1 ] )/2.0;
  }
    
  /***********************************************************************************************************/


  children = calloc( 2, sizeof( nodePtr ) );
  
  /* allocate some memory for children */ 
  children[0] = createNode( r, c );

  /* create right node from parent */
  children[1] = createNode( r, c );


  /* we are using the middle index value so we know the new size */
  children[0]->indexUsed = splitIndex;
  children[1]->indexUsed = c->indexUsed - splitIndex;

  /* now let's have some fun with pointer math */
  for( i = 0; i < splitIndex; i++) {
//    PRINTF("Left: (%p) %d: %d %d %d\n", (void *) c, (int) i,(int) c->index[xPtr[i] - x], (int) (xPtr[i]  - x), (int) c->index[i] );
    children[0]->index[i] = c->index[xPtr[i] - x];
  }
  for( i = splitIndex; i < c->indexUsed; i++) {
//    PRINTF("Right: (%p) %d: %d %d %d\n", (void *) c, (int) i,(int) c->index[xPtr[i] - x], (int) (xPtr[i]  - x), (int) c->index[i] );
    children[1]->index[i -splitIndex] = c->index[ xPtr[i] -x ];
  }

  for( i = 0; i < r->K; i ++ ) {
    if ( i == c->dim ) {
      (children[0]->max)[i] = c->split; 
      (children[1]->min)[i] = c->split; 
    } else {
      (children[0]->max)[i] = c->max[i]; 
      (children[1]->min)[i] = c->min[i]; 
    }
    (children[1]->max)[i] = c->max[i]; 
    (children[0]->min)[i] = c->min[i]; 
  }
  
  free(xPtr);
  free(x); 

  free( c->index );
  c->index = NULL;

  return( children );

}




/* a function to print the tree */
void printTree( rootNodePtr r, nodePtr c ) {

  size_t i;

  PRINTF("node: %p\n", (void *) c);
  if( c->index != NULL) {
    for( i=0; i < c->indexUsed; i++) PRINTF("%d ", (int) c->index[i]); 
  } 
  PRINTF("\n\tleft: %p right %p (split = %f) \n", (void *) c->left, (void*) c->right, c->split );
  PRINTF("\n  min= ");
  for( i = 0; i < r->K; i++) PRINTF("%f ",c->min[i] );
  PRINTF("\n  max= ");
  for( i = 0; i < r->K; i++) PRINTF("%f ",c->max[i] );
  PRINTF("\n");

  if( c->left ) {
    PRINTF("left ");
    printTree( r, c->left);
  }

  if( c->right ) {
    PRINTF("right ");
    printTree( r, c->right);
  }

}



/* function to find the minimal Euclidian distance */
size_t getClosest( rootNodePtr r, nodePtr c, size_t item, double * dist  ) {

  size_t i,j,d;
  size_t closestIndex = r->n;

  size_t K = r->K;
  double * x = r->data;
  double * y = &(r->data[item*K]);
  double currentDist;

//PRINTF("  getClosest: c->indexUsed = %d\n", (int) c->indexUsed );

  for( i = 0; i < c->indexUsed; i++) {

    currentDist = 0;
  
    j = c->index[i]; 

//PRINTF("  getClosest: Checking %d against %d", (int) j, (int) item);

    // check if it's a valid index 
    if( j  >= r->n) {
//PRINTF(" Not Valid\n");
      continue;  
    }
    // don't match what we are not looking for
    if( j == item ) {
//PRINTF(" Equal to Item\n");
      continue; 
    } 

    for( d = 0; d < K; d++) currentDist += (x[j * K + d] - y[d]) * (x[j*K + d] - y[d]);  //calculate distance
    
//PRINTF(" dist = %f", currentDist);

    if( currentDist < *dist ) {
      *dist = currentDist; 
      closestIndex = i;
//PRINTF(" newMin! ");
    }
//PRINTF("\n");


  }

  if( closestIndex < r->n ) {
    return( c->index[closestIndex] );
  }

  return( r->n );
}



/* function to find the minimal Euclidian distance */
size_t getClosestTie( rootNodePtr r, nodePtr c, size_t item, double * dist, double * tieBreak  ) {

  size_t i,j,d;
  size_t closestIndex = r->n;

  size_t K = r->K;
  double * x = r->data;
  double * y = &(r->data[item*K]);
  double currentDist;
  double newTieBreak;

//PRINTF("  getClosest: c->indexUsed = %d\n", (int) c->indexUsed );

  for( i = 0; i < c->indexUsed; i++) {

    currentDist = 0;
  
    j = c->index[i]; 

//PRINTF("  getClosest: Checking %d against %d", (int) j, (int) item);

    // check if it's a valid index 
    if( j  >= r->n) {
//PRINTF(" Not Valid\n");
      continue;  
    }
    // don't match what we are not looking for
    if( j == item ) {
//PRINTF(" Equal to Item\n");
      continue; 
    } 

    for( d = 0; d < K; d++) currentDist += (x[j * K + d] - y[d]) * (x[j*K + d] - y[d]);  //calculate distance
    
//PRINTF(" dist = %f", currentDist);

    if( currentDist < *dist ) {
      *dist = currentDist; 
      closestIndex = i;
    } else if( currentDist == *dist ) {
      newTieBreak = runif(0.0,1.0);
      if( newTieBreak > *tieBreak) *tieBreak = newTieBreak;
      closestIndex = i;
    }


  }

  if( closestIndex < r->n ) {
    return( c->index[closestIndex] );
  }

  return( r->n );
}



/* find the nearest neighbor that is not a specific index */
/* bound should be first set to the value of the node you are trying to find a neighbor for */
void find_nn_notMe( rootNodePtr r, nodePtr c, size_t item, double * dist, size_t * query, double * queryPoint, double * tieBreak ) {

//PRINTF("Entering %p\n", c);
  
  size_t i;
  double distMinLeft, distMaxLeft;
  double distMinRight, distMaxRight;
  double boundDistLeft=INFINITY, boundDistRight=INFINITY;
  size_t queryTmp = r->n;


  /* nothing to search for */
  if( item >= r->n ) {
//PRINTF("nothing to do \n");
    return;
  } 
 
  /* is there anything here ? */ 
  if( c->index != NULL ) {
    //queryTmp = getClosest( r, c, item, dist); 
    queryTmp = getClosestTie( r, c, item, dist, tieBreak); 
    if( queryTmp < r->n) *query = queryTmp;

    return;
  } 

  /* move on */
  /* get bound distance */

  for(i = 0; i < r->K; i++){


    if( c->left != NULL ) {
      distMinLeft = queryPoint[i] - (c->left)->min[i];
      distMinLeft *= distMinLeft;
      
      distMaxLeft = queryPoint[i] - (c->left)->max[i];
      distMaxLeft *= distMaxLeft;

//      Rprintf("L %d queryPoint = %f, min %f, max %f, distMin %f, distMax %f\n", 
//          (int) i, queryPoint[i], (c->left)->min[i], (c->left)->max[i], distMinLeft, distMaxLeft); 

      if( distMinLeft < distMaxLeft ) {
        if( boundDistLeft > distMinLeft ) boundDistLeft = distMinLeft;
      } else {
        if( boundDistLeft > distMaxLeft ) boundDistLeft = distMaxLeft;
      }
    }

    if( c->right != NULL ) {
      distMinRight = queryPoint[i] - (c->right)->min[i];
      distMinRight *= distMinRight;
      
      distMaxRight = queryPoint[i] - (c->right)->max[i];
      distMaxRight *= distMaxRight;
      
//      Rprintf("R %d queryPoint = %f, min %f, max %f, distMin %f, distMax %f\n", 
//          (int) i, queryPoint[i], (c->right)->min[i], (c->right)->max[i], distMinRight, distMaxRight); 

      if( distMinRight < distMaxRight ) {
        if( boundDistRight > distMinRight ) boundDistRight = distMinRight;
      } else {
        if( boundDistRight > distMaxRight ) boundDistRight = distMaxRight;
      }
    }


  }  
    
//PRINTF(" boundDist = L %f R %f ( %f )\n", boundDistLeft, boundDistRight , *dist );
    
 
  /* boundary distance */
  if( r->data[item * r->K + c->dim] <= c->split ) {
    if( boundDistLeft <= *dist ) find_nn_notMe( r, c->left, item, dist, query, queryPoint, tieBreak );  
    if( boundDistRight <= *dist ) find_nn_notMe( r, c->right, item, dist, query, queryPoint, tieBreak );   
  } else {
    if( boundDistRight <= *dist ) find_nn_notMe( r, c->right, item, dist, query, queryPoint, tieBreak );   
    if( boundDistLeft <= *dist ) find_nn_notMe( r, c->left, item, dist, query, queryPoint, tieBreak );  
  }


}


