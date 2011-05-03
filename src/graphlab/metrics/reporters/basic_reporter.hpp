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


#ifndef GRAPHLAB_BASIC_REPORTER
#define GRAPHLAB_BASIC_REPORTER

#include <graphlab/metrics/imetrics_reporter.hpp>

/**
 * Simple metrics reporter that dumps metrics to
 * standard output.
 */

namespace graphlab {

  class basic_reporter : public imetrics_reporter {
    
  public:
    virtual void do_report(std::string name, std::string ident, std::map<std::string, metrics_entry> & entries) {
      // TODO: use reporters
      if (ident != name) {
        std::cout << std::endl <<  " === REPORT FOR " << name << "(" << ident << ") ===" << std::endl;
      } else {
        std::cout << std::endl << " === REPORT FOR " << name << " ===" << std::endl;
      }
                
      // First write numeral, then timings, then string entries
      for(int round=0; round<4; round++) { 
        std::map<std::string, metrics_entry>::iterator it;
        int c = 0;
                    
        for(it = entries.begin(); it != entries.end(); ++it) {
          metrics_entry ent = it->second;
          switch(ent.valtype) {
          case REAL:
          case INTEGER:
            if (round == 0) {
              if (c++ == 0) std::cout << "[Numeric]" << std::endl; 
              std::cout << it->first << ":\t\t";
              if (ent.count > 1) {
                std::cout << ent.value << "\t(count: " << ent.count <<  ", min: " << ent.minvalue <<
                  ", max: " << ent.maxvalue << ", avg: " 
                          << ent.cumvalue/(double)ent.count << ")" << std::endl;
              } else {
                std::cout << ent.value << std::endl;
              }
            }
            break;
          case TIME:
            if (round == 1) {
              if (c++ == 0) std::cout << "[Timings]" << std::endl;
              std::cout << it->first << ":\t\t";
              if (ent.count>1) {
                std::cout << ent.value << "s\t (count: " << ent.count << ", min: " << ent.minvalue <<
                  "s, " << "max: " << ent.maxvalue << ", avg: " 
                          << ent.cumvalue/(double)ent.count << "s)" << std::endl;
              } else {
                std::cout << ent.value << " s" << std::endl;
              }
            }
            break;
          case STRING:
            if (round == 2) {
              if (c++ == 0) std::cout << "[Other]" << std::endl;
              std::cout << it->first << ":\t";
              std::cout << ent.stringval << std::endl;
            }
            break;
          case VECTOR:
            if (round == 3) {
              if (c++ == 0) std::cout << "[Numeric]" << std::endl; 
              std::cout << it->first << ":\t\t";
              if (ent.count > 1) {
                std::cout << ent.value << "\t(count: " << ent.count <<  ", min: " << ent.minvalue <<
                  ", max: " << ent.maxvalue << ", avg: " 
                          << ent.cumvalue/(double)ent.count << ")" << std::endl;
              } else {
                std::cout << ent.value << std::endl;
              }
              std::cout << it->first << ".values:\t\t";
              for(size_t j=0; j<ent.v.size(); j++) std::cout << ent.v[j] << ",";
              std::cout << std::endl;
            }
            break;
          }
        }
      }
    };
        
  };
    
};



#endif

