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

#ifndef GRAPHLAB_BINARY_PARSER_HPP
#define GRAPHLAB_BINARY_PARSER_HPP

#include <iostream>
#include <fstream>

namespace graphlab {

  /**
   * \ingroup util_internal
   * A thin wrapper around ifstream to provide simplicity of reading
   * of binary data.
   * \see binary_output_stream
   */
  class binary_input_stream : public std::ifstream {
    typedef std::ifstream base_type;
    using base_type::bad;
  public:
    binary_input_stream(const char* fname) :
      base_type(fname, std::ios::binary | std::ios::in) {
      assert(bad() == false);
    }
    
    /**
     * Read an arbitrary type.
     */
    template<typename T> T read() {
      T t;
      base_type::read(reinterpret_cast<char*>(&t), sizeof(T));
      if(bad()) {
        std::cout << "Error reading file!" << std::endl;
        assert(false);
      }
      return t;
    }
  };



  /**
   * \ingroup util_internal
   * A thin wrapper around ifstream to provide simplicity of writing
   * of binary data.
   * \see binary_input_stream
   */
  class binary_output_stream : public std::ofstream {
  typedef std::ofstream base_type;
    using std::ofstream::bad;
  public:
    binary_output_stream(const char* fname) : 
    std::ofstream(fname, std::ios::binary | std::ios::out) {
      assert(bad() == false);
    }
    
    //! Write the arbitrary data type to file
    template<typename T> void write(T t) {
      base_type::write(reinterpret_cast<char*>(&t), sizeof(T));
      if(bad()) {
        std::cout << "Error writing file!" << std::endl;
        assert(false);
      }
    }
  }; // end of binary_output_stream

  

}




#endif




