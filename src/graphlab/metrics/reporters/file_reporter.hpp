
#ifndef GRAPHLAB_FILE_REPORTER
#define GRAPHLAB_FILE_REPORTER

#include <graphlab/metrics/imetrics_reporter.hpp>
 

/**
 * Simple metrics reporter that dumps metrics to
 * a file.
 */

namespace graphlab {

  class file_reporter : public imetrics_reporter {
  private:
    file_reporter() {}
        
    std::string filename;
    FILE * f;
  public:
            
        
    file_reporter(std::string fname) : filename(fname) {
      // Create new file
      f = fopen(fname.c_str(), "w");
      assert(f != NULL);
    }
            
    virtual void do_report(std::string name, std::string ident, 
                           std::map<std::string, metrics_entry> & entries) {
      if (ident != name) {
        fprintf(f, "[%s:%s]\n", name.c_str(), ident.c_str());
      } else {
        fprintf(f, "[%s]\n", name.c_str());
      }
      std::map<std::string, metrics_entry>::iterator it;
                    
      for(it = entries.begin(); it != entries.end(); ++it) {
        metrics_entry ent = it->second;
        switch(ent.valtype) {
        case INTEGER:
                            
          fprintf(f, "%s.%s=%ld\n", ident.c_str(), it->first.c_str(), (long int) (ent.value));
          fprintf(f, "%s.%s.count=%d\n", ident.c_str(), it->first.c_str(), ent.count);
          fprintf(f, "%s.%s.min=%ld\n", ident.c_str(), it->first.c_str(), (long int) (ent.minvalue));
          fprintf(f, "%s.%s.max=%ld\n", ident.c_str(), it->first.c_str(), (long int) (ent.maxvalue));
          fprintf(f, "%s.%s.avg=%lf\n", ident.c_str(), it->first.c_str(), ent.cumvalue/ent.count);
          break;
        case REAL:
        case TIME:
          fprintf(f, "%s.%s=%lf\n", ident.c_str(), it->first.c_str(),  (ent.value));
          fprintf(f, "%s.%s.count=%d\n", ident.c_str(), it->first.c_str(), ent.count);
          fprintf(f, "%s.%s.min=%lf\n", ident.c_str(), it->first.c_str(),  (ent.minvalue));
          fprintf(f, "%s.%s.max=%lf\n", ident.c_str(), it->first.c_str(),  (ent.maxvalue));
          fprintf(f, "%s.%s.avg=%lf\n", ident.c_str(), it->first.c_str(), ent.cumvalue/ent.count);
          break;
        case STRING:
          fprintf(f, "%s.%s=%s\n", ident.c_str(), it->first.c_str(), it->second.stringval.c_str());                                
          break;
        case VECTOR:
          fprintf(f, "%s.%s.values=",  ident.c_str(), it->first.c_str());
          for(int j=0; j<ent.v.size()-1; j++) fprintf(f, "%lf,", ent.v[j]);
          fprintf(f, "%lf\n", ent.v[ent.v.size()-1]);
          fprintf(f, "%s.%s=%lf\n", ident.c_str(), it->first.c_str(),  (ent.value));
          fprintf(f, "%s.%s.count=%d\n", ident.c_str(), it->first.c_str(), ent.count);
          fprintf(f, "%s.%s.min=%lf\n", ident.c_str(), it->first.c_str(),  (ent.minvalue));
          fprintf(f, "%s.%s.max=%lf\n", ident.c_str(), it->first.c_str(),  (ent.maxvalue));
          fprintf(f, "%s.%s.avg=%lf\n", ident.c_str(), it->first.c_str(), ent.cumvalue/ent.count);
          break;
        }
      }
                
      fflush(f);
    };
        
  };
    
};



#endif


