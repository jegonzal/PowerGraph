/* Copyright (c) 2003 Alexander Ihler
 * Original code from: http://www.ics.uci.edu/~ihler/code/index.html
 *
 *     All rights reserved.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing,
 *  software distributed under the License is distributed on an "AS
 *  IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 *  express or implied.  See the License for the specific language
 *  governing permissions and limitations under the License. */

 //////////////////////////////////////////////////////////////////////////////////////
// BallTreeClass  --  class definitions for a BallTree (actually KD-tree) 
//                    object, primarily for use in matlab MEX files.
//
// See BallTree.h for the class definition.
//
//////////////////////////////////////////////////////////////////////////////////////
//
// Written by Alex Ihler and Mike Mandel
// Copyright (C) 2003 Alexander Ihler
//
//////////////////////////////////////////////////////////////////////////////////////

//#define MEX
#include <math.h>
#include "BallTree.h"
#include <utility>
#include <map>
#include <assert.h>

#include "../kde.h"

const char* BallTree::FIELD_NAMES[] = {"D", "N", "centers", "ranges", "weights",
            "lower", "upper", "leftch", "rightch", "perm"};
const int BallTree::nfields = 10;

// Given the leaves, build the rest of the tree from the top down.
// Split the leaves along the most spread coordinate, build two balls
// out of those, and then build a ball around those two children.
void BallTree::buildBall(BallTree::index low, BallTree::index high, BallTree::index root)
{
  // special case for N=1 trees
  if(low == high) {
    lowest_leaf[root] = low;
    highest_leaf[root] = high;
    left_child[root] = low;

    // point right child to the same as left for calc stats, and then
    // point it to the correct NO_CHILD afterwards.  kinda kludgey
    right_child[root] = high;
    calcStats(root);
    right_child[root] = NO_CHILD;

    return;
  }

  BallTree::index coord, split, left, right;
  coord = most_spread_coord(low, high);    // find dimension of widest spread
  
  // split the current leaves into two groups, to build balls on them.
  // Chose the most spread coordinate to split them on, and make sure
  // there are the same number of points in each (+-1 for round off
  // error).
  split = (low + high) / 2;
  select(coord, split, low, high);

  // an alternative is to use partition, but that doesn't deal well
  // with repeated numbers and it doesn't split into balanced sets.
//   split = partition(coord, low, high);

  // if the left sub-tree is just one leaf, don't make a new non-leaf
  // node for it, just point left_idx directly to the leaf itself.
  if(split <= low)    left = low;
  else                left = next++;

  // same for the right
  if(split+1 >= high) right = high;
  else                right = next++;

  lowest_leaf[root]  = low;
  highest_leaf[root] = high;
  left_child[root]   = left;
  right_child[root]  = right;

  // build sub-trees if necessary
  if(left != low)    buildBall(low, split, left);
  if(right != high)  buildBall(split+1, high, right);

  calcStats(root);
}

// Find the dimension along which the leaves between low and high
// inclusive have the greatest variance
BallTree::index BallTree::most_spread_coord(BallTree::index low, BallTree::index high) const
{
  BallTree::index dimension, point, max_dim;
  double mean, variance, max_variance;

  max_variance = 0;
  max_dim = 0;
  assert(dims >0);

  for(dimension = 0; dimension<dims; dimension++) {
    mean = 0;
    for(point = dims*low + dimension; point < dims*high; point += dims)
      mean += centers[point];
    mean /= (high - low);

    variance = 0;
    for(point = dims*low + dimension; point < dims*high; point += dims)
      variance += (centers[point] - mean) * (centers[point] - mean);
    if(variance > max_variance) {
      max_variance = variance;
      max_dim = dimension;
    }
  }

  return max_dim;
}


// straight from CLR, the unrandomized partition algorithm for
// quicksort.  Partitions the leaves from low to high inclusive around
// a random pivot in the given dimension.  Does not affect non-leaf
// nodes, but does relabel the leaves from low to high.
BallTree::index BallTree::partition(unsigned int dimension, BallTree::index low, 
				  BallTree::index high) 
{
  BallTree::index pivot;

  pivot = low;  // not randomized, could set pivot to a random element

  while(low < high) {
    while(centers[dims*high + dimension] >= centers[dims*pivot + dimension])
      high--;
    while(centers[dims*low + dimension] < centers[dims*pivot + dimension])
      low++;
    
    swap(low, high);
    pivot = high;
  }

  return high;
}


// Function to partition the data into two (equal-sized or near as possible)
//   sets, one of which is uniformly greater than the other in the given
//   dimension.
void BallTree::select(unsigned int dimension, BallTree::index position,
		      BallTree::index low, BallTree::index high)
{
  BallTree::index m,r,i;
  
  while (low < high) {
    r = (low + high)/2; 
    swap(r,low);
    m = low;
    for (i=low+1; i<=high; i++) {
      if (centers[dimension+dims*i] < centers[dimension+dims*low]) {
        m++;
        swap(m,i);
      } 
    }
    swap(low,m);
    if (m <= position) low=m+1;
    if (m >= position) high=m-1;
  }    
}


// Swap the ith leaf with the jth leaf.  Actually, only swap the
// weights, permutation, and centers, so only for swapping
// leaves. Will not swap ranges correctly and will not swap children
// correctly.
void BallTree::swap(BallTree::index i, BallTree::index j) 
{
  BallTree::index k;
  double tmp;

  if (i==j) return;

  // swap weights
  tmp = weights[i];    weights[i] = weights[j];          weights[j] = tmp;

  // swap perm
  k = permutation[i];  permutation[i] = permutation[j];  permutation[j] = k;

  // swap centers
  i *= dims;   j *= dims;
  for(k=0; k<dims; i++,j++,k++) {
    tmp = centers[i];   centers[i]  = centers[j];   centers[j]  = tmp;
  }
}

//
// Calculate the statistics of level "root" based on the statistics of
//   its left and right children.
//
void BallTree::calcStats(BallTree::index root)
{
  BallTree::index Ni, NiL, NiR;
  index d;

  BallTree::index leftI = left(root), rightI=right(root);   // get children indices 
  if (!validIndex(leftI) || !validIndex(rightI)) return;    // nothing to do if this
                                                            //   isn't a parent node
  assert(dims > 0);

  // figure out the center and ranges of this ball based on it's children
  double max, min;
  for(d=0; d<dims; d++) {
    if (center(leftI)[d] + range(leftI)[d] > center(rightI)[d] + range(rightI)[d])
      max = center(leftI)[d] + range(leftI)[d];
    else
      max = center(rightI)[d] + range(rightI)[d];

    if (center(leftI)[d] - range(leftI)[d] < center(rightI)[d] - range(rightI)[d])
      min = center(leftI)[d] - range(leftI)[d];
    else
      min = center(rightI)[d] - range(rightI)[d];

    centers[root*dims+d] = (max+min) / 2;
    ranges[root*dims+d] = (max-min) / 2;
  }    
  
  // if the left ball is the same as the right ball (should only
  // happen when calling the function directly with the same argument
  // twice), don't count the weight twice
  if(leftI != rightI)
    weights[root] = weights[leftI] + weights[rightI];
  else
    weights[root] = weights[leftI];
}
  

// Public method to build the tree, just calls the private method with
// the proper starting arguments.
void BallTree::buildTree()
{
  BallTree::index i,j;
  for (j=0, i=num_points; j<num_points; i++,j++) {
    for(index k=0; k<dims; k++)
      ranges[i*dims+k] = 0;
 
    lowest_leaf[i] = highest_leaf[i] = i; 
    left_child[i] = i; 
    right_child[i] = NO_CHILD;
    permutation[i] = j;
  }
  next = 1;

  buildBall(num_points, 2*num_points - 1, 0);
}

// Figure out which of two children in this tree is closest to a given
// ball in another tree.  Returns the index in this tree of the closer
// child.
BallTree::index BallTree::closer(BallTree::index myLeft, BallTree::index myRight, const BallTree& otherTree,
			       BallTree::index otherRoot) const 
{
  if (myRight==NO_CHILD || otherRoot==NO_CHILD) return myLeft;
  double dist_sq_l = 0, dist_sq_r = 0;
  assert(dims >0);
 
  for(int i=0; i<dims; i++) {
    dist_sq_l += (otherTree.center(otherRoot)[i] - center(myLeft)[i]) * 
      (otherTree.center(otherRoot)[i] - center(myLeft)[i]);
    dist_sq_r += (otherTree.center(otherRoot)[i] - center(myRight)[i]) * 
      (otherTree.center(otherRoot)[i] - center(myRight)[i]);
  }

  if (dist_sq_l < dist_sq_r)
    return myLeft;
  else 
    return myRight;
}

//
// Perform a *slight* adjustment of the tree: move the points by delta, but
//   don't reform the whole tree; just fix up the statistics.
//
void BallTree::movePoints(double* delta)
{
  assert(dims >0);
  index i;
  for (i=leafFirst(root());i<=leafLast(root());i++)
    for (unsigned int k=0;k<dims;k++)                   // first adjust locations by delta
      centers[dims*i+k] += delta[ getIndexOf(i)*dims + k ];
  for (i=num_points-1; i != 0; i--)                     // then recompute stats of
    calcStats(i);                                       //   parent nodes
  calcStats(root());                                    //   and finally root node
}

// Assumes newWeights is the right size (num_points)
void BallTree::changeWeights(const double *newWeights) {

  assert(num_points>0);
  for(index i=num_points, j=0; i<num_points*2; i++, j++)
    weights[i] = newWeights[ getIndexOf(i) ];

  for (index i=num_points-1; i != 0; i--)
    calcStats(i);
  calcStats(root());
}



/////////////////////// k nearest neighbors functions ///////////////////////

// returns distance squared
double BallTree::minDist(index myBall, const double* point) const
{
  double dist = 0, tmp;

  assert(dims >0);
  for(index i=0; i<dims; i++) {
    tmp = fabs(center(myBall)[i] - point[i]) - range(myBall)[i];
    if(tmp >= 0)
      dist += tmp*tmp;
  }
  
  return dist;
}
double BallTree::maxDist(index myBall, const double* point) const
{
  double dist = 0, tmp;

  assert(dims >0);
  for(index i=0; i<dims; i++) {
    tmp = fabs(center(myBall)[i] - point[i]) + range(myBall)[i];
    dist += tmp*tmp;
  }
  
  return dist;
}


typedef std::multimap<double, BallTree::index> myMap;

// nns is a K x N matrix
// dists is a 1 x N vector
// points is a D x N matrix
void BallTree::kNearestNeighbors(index *nns, double *dists, const double *points, 
				 int N, int k) 
  const
{
  myMap m;
  int leavesDone = 0;
  index lastBall;
  assert(dims >0);
  assert(N>0);

  for(index target = 0; target < N*dims; target += dims, ++leavesDone) {  

    int nnsSoFar = 0;
    double leastDist = maxDist(root(), points+target);

    m.insert(myMap::value_type(minDist(root(), points+target), root()));
    // examining points in order of min dist means that when you see a
    // point, it is the closest point left
    
    while(! m.empty() && nnsSoFar < k) {
      index current = (*m.begin()).second;
      m.erase(m.begin());
      
      if(isLeaf(current)) {
	// since the nodes are sorted by minDist, a leaf at the front of
	// the pq must be the next nearest neighbor
	nns[leavesDone*k + nnsSoFar++] = current;
	lastBall = current;
      } else {  // not a leaf
	// push both children
	m.insert(myMap::value_type(minDist(left_child[current], points+target), 
				   left_child[current]));
	m.insert(myMap::value_type(minDist(right_child[current], points+target), 
				   right_child[current]));
      }
      
      bool keepInnerPruning = true;
      myMap::iterator it = m.begin(), last;
      // get rid of ineligible balls and find the min distance
      while(it != m.end()) {
	current = (*it).second;
	double max = maxDist(current, points+target);
	if(max < leastDist && Npts(current) >= k - nnsSoFar)
	  leastDist = max;
	
//	if((*it).first > leastDist) {
//	  // we see the points in order of minDist, so once we see one
//	  // that's too big, the rest will also be too big
//	  m.erase(it, m.end());
//	  break;
//	}
	
	myMap::iterator last = it++;
	if(keepInnerPruning) {
	  if(it != m.end() && max < (*it).first && Npts(current) < k - nnsSoFar) {
	    // if the closest ball doesn't have too many points and all of
	    // them are nearer than any other ball, include them all
	    //   note "<" not "<=" so that *dist will be right
	    lastBall = current;
	    for(index i=leafFirst(current); i <= leafLast(current); i++)
	      nns[leavesDone*k + nnsSoFar++] = i;
	    m.erase(last);
	  } else {
	    keepInnerPruning = false;
	  }
	}
      } // end pruning
    } // end single nearest neighbor

    // clear out the remaining points
    m.clear();

    dists[leavesDone] = sqrt(maxDist(lastBall, points+target));
    
    index i;
    for(i=leavesDone*k; i < leavesDone*k+nnsSoFar; i++)
      nns[i] = getIndexOf(nns[i]);
    for(i=leavesDone*k+nnsSoFar; i<(leavesDone+1)*k; i++)
      nns[i] = NO_CHILD;
  } // end all nearest neighbors
}

/////////////////////////////// matlab functions ////////////////////////////

// Constructor that doesn't initialize members, so that they can be
// set by the loadFromMatlab and createInMatlab functions.
//BallTree::BallTree() : next(1) {}


#ifdef MEX

// Load the arrays already allocated in matlab from the given
// structure.
BallTree::BallTree(const mxArray* structure) 
{
  dims       = (unsigned int) mxGetScalar(mxGetField(structure,0,"D")); // get the dimensions
  num_points = (BallTree::index) mxGetScalar(mxGetField(structure,0,"N")); //
  
  centers = (double*) mxGetPr(mxGetField(structure,0,"centers"));
  ranges  = (double*) mxGetPr(mxGetField(structure,0,"ranges"));
  weights = (double*) mxGetPr(mxGetField(structure,0,"weights"));

  lowest_leaf = (BallTree::index*) mxGetData(mxGetField(structure,0,"lower"));
  highest_leaf= (BallTree::index*) mxGetData(mxGetField(structure,0,"upper"));
  left_child  = (BallTree::index*) mxGetData(mxGetField(structure,0,"leftch"));
  right_child = (BallTree::index*) mxGetData(mxGetField(structure,0,"rightch"));
  permutation = (BallTree::index*) mxGetData(mxGetField(structure,0,"perm"));

  next = 1;    // unimportant
}

// Create new matlab arrays and put them in the given structure.
mxArray* BallTree::createInMatlab(const mxArray* _pointsMatrix, const mxArray* _weightsMatrix)
{
  mxArray* structure;
  structure = matlabMakeStruct(_pointsMatrix,_weightsMatrix);
  BallTree bt(structure);
  if (bt.Npts() > 0) bt.buildTree();

  return structure;
}

// Create new matlab arrays and put them in the given structure.
mxArray* BallTree::matlabMakeStruct(const mxArray* _pointsMatrix, const mxArray* _weightsMatrix)
{
  mxArray* structure;
  BallTree::index i, j;
  double *_points, *_weights;
  
  // get fields from input arguments
  unsigned int Nd = mxGetM(_pointsMatrix);
  BallTree::index Np = mxGetN(_pointsMatrix);
  _points  = (double*)mxGetData(_pointsMatrix);
  _weights = (double*)mxGetData(_weightsMatrix);

  // create structure, populate it, and get handles to the arrays
  structure = mxCreateStructMatrix(1, 1, nfields, FIELD_NAMES);
  
  mxSetField(structure, 0, "D",       mxCreateDoubleScalar((double) Nd));
  mxSetField(structure, 0, "N",       mxCreateDoubleScalar((double) Np));

  mxSetField(structure, 0, "centers", mxCreateDoubleMatrix(Nd, 2*Np, mxREAL));
  mxSetField(structure, 0, "ranges",  mxCreateDoubleMatrix(Nd, 2*Np, mxREAL));
  mxSetField(structure, 0, "weights", mxCreateDoubleMatrix(1, 2*Np, mxREAL));

  mxSetField(structure, 0, "lower",   mxCreateNumericMatrix(1, 2*Np, mxUINT32_CLASS, mxREAL));
  mxSetField(structure, 0, "upper",   mxCreateNumericMatrix(1, 2*Np, mxUINT32_CLASS, mxREAL));
  mxSetField(structure, 0, "leftch",  mxCreateNumericMatrix(1, 2*Np, mxUINT32_CLASS, mxREAL));
  mxSetField(structure, 0, "rightch", mxCreateNumericMatrix(1, 2*Np, mxUINT32_CLASS, mxREAL));
  mxSetField(structure, 0, "perm",    mxCreateNumericMatrix(1, 2*Np, mxUINT32_CLASS, mxREAL));

  // initialize arrays
  double* centers = (double *) mxGetData(mxGetField(structure, 0, "centers"));
  double* weights = (double *) mxGetData(mxGetField(structure, 0, "weights"));
  for (j=0,i=Nd*Np; j<Nd*Np; i++,j++)
    centers[i] = _points[j];
  for (j=0,i=Np; j<Np; i++,j++)
    weights[i] = _weights[j];

  return structure;
}

#else
BallTree::BallTree(const kde& structure) 
{
  //dims       = (unsigned int) mxGetScalar(mxGetField(structure,0,"D")); // get the dimensions
  dims = structure.centers.rows();
  assert(dims == 1);
  //num_points = (BallTree::index) mxGetScalar(mxGetField(structure,0,"N")); //
  num_points = structure.centers.cols();
  assert(num_points >= 1); 
 
  //centers = (double*) mxGetPr(mxGetField(structure,0,"centers"));
  itpp::vec temp = structure.centers.get_row(0);
  centers = new double[dims*2*num_points];
  for (int i=0; i< num_points; i++)
     centers[i+dims*num_points] = temp[i];
  //ranges  = (double*) mxGetPr(mxGetField(structure,0,"ranges"));

  ranges = new double[2 * dims*num_points];
  //weights = (double*) mxGetPr(mxGetField(structure,0,"weights"));
  //weights = vec2vec(&structure.weights);
  weights = new double[2*dims*num_points];
  for (int i=0; i< num_points; i++)
     weights[i+dims*num_points] = structure.weights(i);

  lowest_leaf = new unsigned int[2*num_points]; //(BallTree::index*) mxGetData(mxGetField(structure,0,"lower"));
  highest_leaf=  new unsigned int[2*num_points];//(BallTree::index*) mxGetData(mxGetField(structure,0,"upper"));
  left_child  =  new unsigned int[2*num_points];//(BallTree::index*) mxGetData(mxGetField(structure,0,"leftch"));
  right_child =  new unsigned int [2*num_points];//BallTree::index*) mxGetData(mxGetField(structure,0,"rightch"));
  permutation =  new unsigned int [2*num_points];//BallTree::index*) mxGetData(mxGetField(structure,0,"perm"));
   
  next = 1;    // unimportant

  //buildTree();
}

void BallTree::clean(){
  if (lowest_leaf != NULL)
       delete[] lowest_leaf;
  if (highest_leaf != NULL)
       delete[] highest_leaf;
  if (left_child != NULL)
       delete[] left_child;
  if (right_child != NULL)
       delete[] right_child;
  if (permutation != NULL)
       delete[] permutation;
  if (centers != NULL)
       delete[] centers;
  if (weights!= NULL)
       delete[] weights;
  dims = num_points = 0;
}

BallTree::~BallTree(){
/*  if (lowest_leaf != NULL)
       delete[] lowest_leaf;
  if (highest_leaf != NULL)
       delete[] highest_leaf;
  if (left_child != NULL)
       delete[] left_child;
  if (right_child != NULL)
       delete[] right_child;
  if (permutation != NULL)
       delete[] permutation;
  if (centers != NULL)
       delete[] centers;
  if (weights!= NULL)
       delete[] weights;
*/
//  printf("was in destructor\n");
}
#endif
