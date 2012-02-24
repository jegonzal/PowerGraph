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
