/* Copyright (c) 2015  Jonathan Lisic 
 * Last edit: 15/07/08 - 11:45:28
 * License: GPL (>=2) 
 */  

#ifndef KDTREE_HEADER

#define KDTREE_HEADER


#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>


/* create a node */
struct node {
  size_t dim;         /* dim split across */
  size_t * index;     /* index rows for x */
  size_t indexUsed;   /* since we are using arrays, this indicates the number of elements of each array */
  double split;       /* point split on for dim */
  double * min;       /* min value along dim */
  double * max;       /* max value along dim */
  struct node * head; /* used for delete heavy trees */
  struct node * left; /* values less than or equal to split */
  struct node * right; /* value greater than or equal to split */ 
};

/* create typedef */
typedef struct node node;
typedef struct node * nodePtr;


/* tree node */
struct rootNode {
  size_t K;                 /* dims of x */
  size_t leafSize;          /* max number of elements on each leaf */
  size_t n;                 /* number of rows in x */ 
  size_t ** pointerIndex;   /* pointer index to allow easy changing of returned values */
  double * data;            /* pointer to x */
  nodePtr root;             /* root node pointer */
};

/* create typedef */
typedef struct rootNode rootNode;
typedef struct rootNode * rootNodePtr;


/* function to print a tree */
void printTree( rootNodePtr r, nodePtr c ); 


/* function to create a new Tree */
rootNodePtr createTree( size_t K, size_t leafSize, size_t n, double * data );

/* delete tree */
void deleteTree( rootNodePtr r );

/* build index for the tree */
void buildIndex( rootNodePtr r, nodePtr c, nodePtr p ); 

/* qsort comparison function for doubles */
int comp ( const void * a, const void * b ); 

/* qsort comparison function for double * */ 
int compareDoublePtr ( const void * aPtr, const void * bPtr );

/* create Node */
nodePtr createNode( rootNodePtr r, nodePtr p ); 

/* delete a node and all of it's children */
void deleteNode( rootNodePtr r, nodePtr c ); 

/* create children */
nodePtr * createChildren( rootNodePtr r, nodePtr c, nodePtr p); 

/* function to get the closest neighbor */
size_t getClosest( rootNodePtr r, nodePtr c, size_t item, double * dist  ); 

/* function to get the closest neighbor, with tie handling */
size_t getClosestTie( rootNodePtr r, nodePtr c, size_t item, double * dist, double * tieBreak  ); 

/* funciton to find neighbors */
void find_nn_notMe( rootNodePtr r, nodePtr c, size_t item, double * dist, size_t * query, double * queryPoint, double * tieBreak  ); 

#endif
