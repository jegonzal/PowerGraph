#include <cassert>
#include <cstdlib>
#include <cmath>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>




// MAIN
// ============================================================================>

inline std::string trim(const std::string& str) {
  std::string::size_type pos1 = str.find_first_not_of(" \t\r");
  std::string::size_type pos2 = str.find_last_not_of(" \t\r");
  return str.substr(pos1 == std::string::npos ? 0 : pos1,
                    pos2 == std::string::npos ? str.size()-1 : pos2-pos1+1);
}


int main(int argc, char** argv) {
  std::cout << "This program computes the difference between two belief files."
            << std::endl;

  if(argc != 3) {
    std::cout << "Usage example: " << std::endl
              << " %> " << argv[0] << " bp_beliefs.csv gibbs_beliefs.csv "
              << std::endl;
  }


  std::ifstream fin1(argv[1]);
  if(!fin1.good()) {
    std::cout << "Unable to open: " << argv[1] << std::endl;
    exit(EXIT_FAILURE);
  }
  std::ifstream fin2(argv[2]);
  if(!fin2.good()) {
    std::cout << "Unable to open: " << argv[2] << std::endl;
    exit(EXIT_FAILURE);
  }

  double L1_L1_error    = 0;
  double L1_L_inf_error  = 0;
  double L_inf_L1_error    = 0;
  size_t disagree = 0;
  size_t vertices = 0;
  while(fin1.good() && fin2.good()) {
    // Read the line for both files
    std::string line1;
    std::getline(fin1, line1);
    std::string line2;
    std::getline(fin2, line2);
    line1 = trim(line1);
    line2 = trim(line2);
    if(line1.size() == 0 || line2.size() == 0) continue;


    // Trim the varialble name
    size_t index1 = line1.find_first_of('/');
    assert(index1 != std::string::npos);
    size_t index2 = line2.find_first_of('/');
    assert(index2 != std::string::npos);

    std::string var1 = trim(line1.substr(0, index1));
    std::string var2 = trim(line1.substr(0, index2));
    assert(var1 == var2);

    line1 = line1.substr(index1);
    line2 = line2.substr(index2);
    

    //  compute the rest of the string
    std::stringstream strm1(line1), strm2(line2);
    
    size_t index = 0;
    double sum = 0, max = 0;
    size_t mapInd1 = 0;
    double mapVal1 = -1;
    size_t mapInd2 = 0;
    double mapVal2 = -1;

    while(strm1.good() && strm2.good()) {
      
      double value1=-1, value2=-1;      
      strm1 >> value1;
      strm2 >> value2;
      strm1.ignore(1); strm2.ignore(1);
      
      assert(value1 >= 0 && value1 <= 1);
      assert(value2 >= 0 && value2 <= 1);
      // Compute the map for each line
      if(value1 > mapVal1) {
        mapVal1 = value1;
        mapInd1 = index;
      }
      if(value2 > mapVal2) {
        mapVal2 = value2;
        mapInd2 = index;
      }
      // Compute the difference in the probabilities 
      double diff = std::abs(value1 - value2);
      sum += diff;
      max = std::max(max, diff);
      // Increment index
      ++index;
    }

    // Assert both lines end at the same state
    assert(!strm1.good() && !strm2.good());

    // update the global counters
    double l1_error = sum / index;
    L1_L1_error += l1_error;
    L1_L_inf_error += max;
    L_inf_L1_error = std::max(L_inf_L1_error, l1_error);
    disagree += (mapInd1 != mapInd2? 1 : 0);
    ++vertices;
  } // end of while loop
  assert(!fin1.good() && !fin2.good());
  fin1.close(); fin2.close();

  std::cout << "Read " << vertices << " beliefs."
            << std::endl
            << "L1 L1 error:       " << L1_L1_error / vertices << std::endl
            << "L1 Linf error:     " << L1_L_inf_error / vertices << std::endl
            << "Linf L1 error:     " << L_inf_L1_error << std::endl
            << "%Map Disagree:      " << double(disagree) / vertices << std::endl;  
  return EXIT_SUCCESS;
}
