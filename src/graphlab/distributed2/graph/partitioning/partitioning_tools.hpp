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

#ifndef GRAPHLAB_PARTITIONING_TOOLS
#define GRAPHLAB_PARTITIONING_TOOLS


#include <string>

namespace graphlab {


  namespace partitioning_tools {
    
    void construct_partitioning(int argc, char** argv, 
                                int numparts,
                                const std::string& path,
                                const double acceptance_probability);
  };
  



}; // end graphlab namspace


#endif

