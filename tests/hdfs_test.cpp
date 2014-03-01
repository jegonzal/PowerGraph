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


#include <iostream>
#include <vector>


#include <graphlab/util/hdfs.hpp>

int main(int argc, char **argv) {
  {
    graphlab::hdfs hdfs;
    const bool write = true;
    graphlab::hdfs::fstream file(hdfs, "/tmp/joeytest.txt", write);
    file.good();
    file << "Hello World\n";
    file.close();
    std::vector<std::string> files = hdfs.list_files("/tmp/");
    for(size_t i = 0; i < files.size(); ++i) 
      std::cout << files[i] << std::endl;
  }

  {
    graphlab::hdfs hdfs;
    graphlab::hdfs::fstream file(hdfs, "/tmp/joeytest.txt");
    file.good();
    std::string answer;
    std::getline(file, answer);
    std::cout << "contents: " << std::endl;
    std::cout << answer << std::endl;
    file.close();
  }
    std::cout << "Done!" << std::endl;
}
