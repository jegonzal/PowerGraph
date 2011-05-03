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


#ifndef GRAPHLAB_IMETRICS_REPORTER
#define GRAPHLAB_IMETRICS_REPORTER

#include <map>

#include <graphlab/metrics/metrics.hpp>

namespace graphlab {
   class imetrics_reporter {
   
        public:
            virtual void do_report(std::string name, std::string id, std::map<std::string, metrics_entry> &  entries) = 0;
    
   };

};

#endif


