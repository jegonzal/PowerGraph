#include <iostream>
#include <vector>
namespace graphlab {
  class pds {
   public:
   pds() {}
   std::vector<size_t> get_pds(p) {
        // verify pdsness
        int pdslength = p *p + p + 1;
        std::vector<int> count(pdslength, 0);
        for (int i = 0;i < result.size(); ++i) {
          for (int j = 0;j < result.size(); ++j) {
            if (i == j) continue;
            count[(result[i] - result[j] + pdslength) % pdslength]++;
          }
        }
        bool ispds = true;
        for (int i = 1;i < count.size(); ++i) {
          if (count[i] != 1) ispds = false;
        }
        if (ispds) {
          return result;
        } else {
          logstream(LOG_ERROR) << "Fail to generate pds for p = " << p << std::endl;
          return std::vector<size_t>();
        }
   }

   private:
    bool test_seq(int a, int b, int c, int p, std::vector<int>& result) {
      std::vector<int> seq;
      int pdslength = p*p + p + 1;
      seq.resize(pdslength + 3);
      seq[0] = 0; seq[1] = 0; seq[2] = 1;
      int ctr = 2;
      for (int i = 3; i < seq.size(); ++i) {
        seq[i] = a * seq[i - 1] + b * seq[i - 2] + c * seq[i - 3];
        seq[i] = seq[i] % p;
        ctr += (seq[i] == 0);
        // PDS must be of length p + 1
        // and are the 0's of seq.
        if (i < pdslength && ctr > p + 1) return false;
      }
      if (seq[pdslength] == 0 && seq[pdslength + 1] == 0){ 
        // we are good to go
        // now find the 0s
        for (int i = 0; i < pdslength; ++i) {
          if (seq[i] == 0) {
            result.push_back(i);
          }
        }
        // probably not necessary. but verify that the result has length p + 1
        if (result.size() != p + 1) {
          result.clear();
          return false;
        }
        return true;
      }
      else {
        return false;
      } 
    }


    std::vector<int> find_pds(int p) {
      std::vector<int> result;
      for (int a = 0; a < p; ++a) {
        for (int b = 0; b < p; ++b) {
          if (b == 0 && a == 0) continue;
          for (int c = 1; c < p; ++c) {
            if (test_seq(a,b,c,p,result)) {
              return result;
            }
          }
        }
      } 
      return result;
    }


  }
}
