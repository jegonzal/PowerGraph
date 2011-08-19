#ifndef _ADVANCED_CONFIG
#define _ADVANCED_CONFIG

class advanced_config{

public:
  bool debug;
  double threshold;
  int syncinterval;
  int iter;
  int ncpus;
  int algorithm;
  bool cg_resid;
  bool zero;
  bool isfloat;
  int unittest;
  bool square;
  std::string datafile;
  bool support_null_variance;
  bool round_robin;
  bool supportgraphlabcf;

  advanced_config(){
    debug = true;
    threshold = 1e-10;
    syncinterval = 10000;
    iter = 20;
    ncpus = 1;
    cg_resid = true;
    zero=true;
    isfloat = true;
    unittest = 0;
    algorithm = 0;
    square = false;
    support_null_variance = false;
    round_robin = true;
    supportgraphlabcf = false;
  }


};


#endif //_ADVANCED_CONFIG
