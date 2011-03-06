
// Example Call:
// Where eW is a sparse matrix of edge weights and vW is a vector of 
// vertex weights and nparts is the number parts in the cut.
// Note that eW should be a symmetric matrix.
//   bestcut = metiscut(eW, vW , nparts);
// NOTE: DO NOT INCLUDE SELF EDGES IN THE ADJACENCY MATRIX


#include "mex.h"

#include <math.h>

#include <vector>
#include <string>
#include <iostream>
#include <fstream>



#include <graphlab/extern/metis/metis.hpp>


// Specify which metis library to use
// #define KMETIS
// #define PMETIS
#define KMETIS


using namespace std;


struct metis_graph {
  size_t n; // Number of vertices
  size_t m; // Number of undirected edges
  metis::idxtype* xadj; // Offsets
  metis::idxtype* adjncy; // Neighbors
  metis::idxtype* vweight; // vertex weights
  metis::idxtype* eweight; // edge weights matches adjncy

  metis_graph(size_t n, size_t m) : 
    n(n), m(m), xadj(NULL), adjncy(NULL), 
    vweight(NULL), eweight(NULL) {
    // mexPrintf("Allocating data structures for %d vertices and %d edges\n",
    //           n, m);
    // mexEvalString("drawnow");
    xadj = new metis::idxtype[n + 1];
    adjncy = new metis::idxtype[m];
    vweight = new metis::idxtype[n];
    eweight = new metis::idxtype[m];
  }

  void destroy() {
    delete [] xadj;
    delete [] adjncy;
    delete [] vweight;
    delete [] eweight;
  }
};


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  // Do input output checking
  if(nlhs == 0) {
    mexErrMsgTxt("At least one output required.");
    return;
  }

  if(nrhs != 3) {
    mexErrMsgTxt("Three inputs required.");
    return;
  }
  
  
  //  mexPrintf("Processing Input:\n");
  // Get the edge weights
  const mxArray* matlab_eW      = prhs[0];
  // Get the vertex weights
  const mxArray* matlab_vW      = prhs[1];
  // Get the cut size k
  const mxArray* matlab_nparts  = prhs[2]; 
  
  // check that input time is correct
  if(mxGetClassID(matlab_eW) != mxDOUBLE_CLASS) {
     mexErrMsgTxt("Edge weight must be double.");
     return;
  }
  
  if(mxGetClassID(matlab_vW) != mxDOUBLE_CLASS) {
     mexErrMsgTxt("Vertex weights must be double.");
     return;
  }


  // Get the number of vertices
  metis::idxtype nverts = std::max(mxGetM(matlab_vW), 
                                   mxGetN(matlab_vW));
  
  // mexPrintf("Number of verts = %d\n", nverts);
  // mexEvalString("drawnow");

  // Get the number of edges (half the number of nonzero entries in E
  metis::idxtype nedges = mxGetNzmax(matlab_eW) ;
  //  mexPrintf("Number of edges = %d\n", nedges);

  // Get a pointer to the vertex weights
  double* vW = mxGetPr(matlab_vW);
 
  // Get pointers to sparse matrix
  mwIndex* jc = mxGetJc(matlab_eW);
  mwIndex* ir = mxGetIr(matlab_eW);
  double* pr  = mxGetPr(matlab_eW);

  mexPrintf("Running allocator\n");
  mexEvalString("drawnow");
  
  //  mexPrintf("Loading Metis Graph struct\n");
  // construct a metis data structure 
  metis_graph g(nverts, nedges);

  mexPrintf("Filling in data structures\n");
  mexEvalString("drawnow");
  // construct teh arguments to metis
  for(int i = 0; i <= nverts; ++i) {
    g.xadj[i] = static_cast<metis::idxtype>(jc[i]);
    if(g.xadj[i] < 0) mexErrMsgTxt("Negative index!");
  }

  for(int i = 0; i < nverts; ++i) {
    g.vweight[i] = static_cast<metis::idxtype>(round(vW[i]) + 1);
    if(g.vweight[i] < 0) mexErrMsgTxt("Negative weight!");
  }
  
  for(int i = 0; i < nedges; ++i) {
    g.adjncy[i] = static_cast<metis::idxtype>(ir[i]);
    if(g.adjncy[i] < 0) mexErrMsgTxt("Negative index!");
    g.eweight[i] = static_cast<metis::idxtype>(round(pr[i]) + 1);
    if(g.eweight[i] < 0) mexErrMsgTxt("Negative edge weight!");
  } 


  // Call metis
  /**
   * 0 No weights (vwgts and adjwgt are NULL) 
   * 1 Weights on the edges only (vwgts = NULL) 
   * 2 Weights on the vertices only (adjwgt = NULL) 
   * 3 Weights both on vertices and edges. 
   */
  metis::idxtype weightflag = 3;
  
  // 0 for C-style numbering starting at 0 (1 for fortran style)
  metis::idxtype numflag = 0;

  // the number of parts to cut into 
  metis::idxtype nparts = 
    static_cast<metis::idxtype>(round(mxGetPr(matlab_nparts)[0]));
  
  //  mexPrintf("Number of parts: %d\n", nparts);

  // Options array (only care about first element if first element is zero
  metis::idxtype options[5]; 
  options[0] = 1; 
  options[1] = 3; 
  options[2] = 1;
  options[3] = 1;
  options[4] = 0;

  // output argument number of edges cut
  metis::idxtype edgecut = 0;

  // output argument the array of assignments
  metis::idxtype* part = new metis::idxtype[nverts];

  // mexPrintf("Running Partitioner\n");
  // mexEvalString("drawnow");
  // Cut the graph
#ifdef KMETIS
//  mexPrintf("Calling kmetis\n");
  metis::METIS_PartGraphKway(&(nverts), 
                             g.xadj,
                             g.adjncy,
                             g.vweight,
                             g.eweight,
                             &(weightflag),
                             &(numflag),
                             &(nparts),
                             options,
                             &(edgecut),
                             part);
#else
//  mexPrintf("Calling pmetis\n");
  metis::METIS_PartGraphRecursive(&(nverts), 
                                  g.xadj,
                                  g.adjncy,
                                  g.vweight,
                                  g.eweight,
                                  &(weightflag),
                                  &(numflag),
                                  &(nparts),
                                  options,
                                  &(edgecut),
                                  part);
#endif

  // Destroy metis graph memory blocks
  g.destroy();

  //  mexPrintf("Processing Results.\n");
  //  mexEvalString("drawnow");

  // load the results
  mxArray* output = mxCreateDoubleMatrix(nverts, 1, mxREAL);
  
  // get the place to load output value
  double* outputvalues = mxGetPr(output);
  
  // write the results
  for(size_t i = 0; i < nverts; ++i) {
     outputvalues[i] = part[i] + 1;
  }
   
  // free the original results;
  delete [] part;
  
  // Save the output
  plhs[0] = output;
  
  if(nlhs == 2) {
    // load the results
    mxArray* output2 = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxGetPr(output2)[0] = edgecut;
    plhs[1] = output2;
  }
} 





  //// Indirect copy  
//  // for each column j
//   size_t metisind = 0;
//   for(size_t j = 0; j < numverts; ++j) {
//     g.vweight[j] = vW[j];
//     g.xadj[j] = metisind;
//     // Loop through all the sparse rows
//     for(size_t offset = jc[j]; offset < jc[j+1]; 
//         ++offset, ++metisind) {
//       size_t i = ir[offset];
//       double value = pr[offset];
//       g.adjncy[metisind] = i;
//       g.eweight[metisind] = value;
//     }
//   }
