

#ifndef COEMAPP_HPP
#define COEMAPP_HPP


#include <graphlab.hpp>


using namespace graphlab;


class coemapp : public iapp<blob_graph> {
	
public:
  coemapp(const char * npsfile, 
          const char * contextfile, 
          const char * matrixfile, 
          const char * seedsfile, 
          const char * negseedsfile);
  ~coemapp();
  
  void start();
  

private:
  
  
  void outputTopVertices(blob_graph &g, uint8_t type);
  int loadVertices(blob_graph * g, const char * filename, uint8_t vtype);
  void loadSeeds(const  char * filename, bool positive_seed, blob_graph &g);
  void loadCoEMDataset(const char * npsfile, const char * contextfile, 
                       const char * matrixfile, const char * seedsfile, 
                       const char * negseedsfile);
  

  /* Used with analyzer_listener */
  virtual global_dumper dump_function();  
  virtual std::vector<std::string> dump_headers();  
  virtual int dump_frequency();
  
  
  const char * npsfile;
  const char * contextfile; 
  const char * matrixfile; 
  const char * seedsfile;
  const char * negseedsfile;
  
  blob_graph * g;

};

#endif

