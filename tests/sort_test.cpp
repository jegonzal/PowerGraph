#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/sample_sort.hpp>
#include <algorithm>

using namespace graphlab;

int main(int argc, char** argv) {
  mpi_tools::init(argc, argv);
  distributed_control dc;
  std::vector<size_t> keys;
  std::vector<size_t> values;
  for (size_t i = 0;i < 1000000; ++i) {
    size_t s = rand();
    keys.push_back(s); values.push_back(s);
  }

  sample_sort<size_t, size_t> sorter(dc);
  sorter.sort(keys.begin(), keys.end(),
              values.begin(), values.end());

  std::vector<std::vector<std::pair<size_t, size_t> > > result(dc.numprocs());
  
  std::swap(result[dc.procid()], sorter.result());
  dc.gather(result, 0);
  if (dc.procid() == 0) {
    // test that it is sorted and the values are correct
    size_t last = 0;
    for (size_t i = 0;i < result.size(); ++i) {
      dc.cout() << result[i].size() << ",";
      for (size_t j = 0; j < result[i].size(); ++j) {
        ASSERT_EQ(result[i][j].first, result[i][j].second);
        ASSERT_GE(result[i][j].first, last);
        last = result[i][j].first;
      }
    }
    dc.cout() << std::endl;
  }
  mpi_tools::finalize();
}
