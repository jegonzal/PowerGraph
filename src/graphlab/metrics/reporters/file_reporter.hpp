/*
This file is part of GraphLab.

GraphLab is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as 
published by the Free Software Foundation, either version 3 of 
the License, or (at your option) any later version.

GraphLab is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public 
License along with GraphLab.  If not, see <http://www.gnu.org/licenses/>.
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


