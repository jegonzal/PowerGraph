/**  
 * Copyright (c) 2009 Carnegie Mellon University. 
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
 *  governing permissions and limitations under the License.
 *
 * For more about this software visit:
 *
 *      http://www.graphlab.ml.cmu.edu
 *
 */



#ifndef GRAPHLAB_FILE_REPORTER
#define GRAPHLAB_FILE_REPORTER

#include <graphlab/metrics/imetrics_reporter.hpp>
#include <fstream>
#include <boost/format.hpp>

/**
 * Simple metrics reporter that dumps metrics to
 * a file.
 */

namespace graphlab {

  class file_reporter : public imetrics_reporter {
  private:
    file_reporter() {}
        
    std::string filename;
    std::ofstream fout;
  public:
            
        
    file_reporter(std::string fname) : filename(fname) {
      // Create new file
      fout.open(fname.c_str());
      ASSERT_TRUE(fout.good());
    }
            
    virtual void do_report(std::string name, std::string ident, 
                           std::map<std::string, metrics_entry> & entries) {
      if (ident != name) {
        
        fout << boost::format("[%1%:%2%]\n") % name.c_str() % ident.c_str();
      } else {
        fout << boost::format("[%1%]\n") % name.c_str();
      }
      std::map<std::string, metrics_entry>::iterator it;
                    
      for(it = entries.begin(); it != entries.end(); ++it) {
        metrics_entry ent = it->second;
        switch(ent.valtype) {
        case INTEGER:
                            
          fout << boost::format("%1%.%2%=%3%\n") % ident.c_str() % it->first.c_str() % (long int) (ent.value);
          fout << boost::format("%1%.%2%.count=%3%\n") % ident.c_str() % it->first.c_str() % ent.count;
          fout << boost::format("%1%.%2%.min=%3%\n") % ident.c_str() % it->first.c_str() % (long int) (ent.minvalue);
          fout << boost::format("%1%.%2%.max=%3%\n") % ident.c_str() % it->first.c_str() % (long int) (ent.maxvalue);
          fout << boost::format("%1%.%2%.avg=%3%\n") % ident.c_str() % it->first.c_str() % (ent.cumvalue/(double)ent.count);
          break;
        case REAL:
        case TIME:
          fout << boost::format("%1%.%2%=%3%\n") % ident.c_str() % it->first.c_str() % (ent.value);
          fout << boost::format("%1%.%2%.count=%3%\n") % ident.c_str() % it->first.c_str() % ent.count;
          fout << boost::format("%1%.%2%.min=%3%\n") % ident.c_str() % it->first.c_str() %  (ent.minvalue);
          fout << boost::format("%1%.%2%.max=%3%\n") % ident.c_str() % it->first.c_str() %  (ent.maxvalue);
          fout << boost::format("%1%.%2%.avg=%3%\n") % ident.c_str() % it->first.c_str() % (ent.cumvalue/(double)ent.count);
          break;
        case STRING:
          fout << boost::format("%1%.%2%=%3%\n") % ident.c_str() % it->first.c_str() % it->second.stringval.c_str();                                
          break;
        case VECTOR:
          fout << boost::format("%1%.%2%.values=") %  ident.c_str() % it->first.c_str();
          for(size_t j=0; j<ent.v.size()-1; j++) fout << boost::format("%1%,") % ent.v[j];
          fout << boost::format("%1%\n") % ent.v[ent.v.size()-1];
          fout << boost::format("%1%.%2%=%3%\n") % ident.c_str() % it->first.c_str() %  (ent.value);
          fout << boost::format("%1%.%2%.count=%3%\n") % ident.c_str() % it->first.c_str() % ent.count;
          fout << boost::format("%1%.%2%.min=%3%\n") % ident.c_str() % it->first.c_str() %  (ent.minvalue);
          fout << boost::format("%1%.%2%.max=%3%\n") % ident.c_str() % it->first.c_str() %  (ent.maxvalue);
          fout << boost::format("%1%.%2%.avg=%3%\n") % ident.c_str() % it->first.c_str() % (ent.cumvalue/(double)ent.count);
          break;
        }
      }
      fout.flush();
    };
        
  };
    
};



#endif


