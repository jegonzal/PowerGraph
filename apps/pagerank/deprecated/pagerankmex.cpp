
#include <string>
#include "mex.h"
#include "pagerankapp.hpp"


void __mexFunction__( int nlhs, mxArray *plhs[],
		int nrhs, const mxArray*prhs[] ) {
	if (nlhs != 1) {
      mexErrMsgTxt("You need one output argument");
     }

    if (nrhs != 2) {
      mexErrMsgTxt("You need two input args!!");
   }

   std::string filename = "/mnt/bigbrofs/usr5/graphlab/testdata/pagerank/p2p-Gnutella08.txt";
   pagerankapp app(filename, "", false);
  
   int argc = 0;
   char *argv[1];
   argv[0] = "pagerankapp";
   
   app.parse_args(argc, argv);
   app.start();		
}
 
void __at_exit__() {
	
}