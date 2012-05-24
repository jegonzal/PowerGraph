#include <graphlab/util/hdfs.hpp>

namespace graphlab {


hdfs& hdfs::get_hdfs() {
  static hdfs fs;
  return fs;
}

}