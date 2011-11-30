#include <stdint.h>
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <graphlab/serialization/serialization_includes.hpp>
#include <graphlab/util/timer.hpp>
#include <graphlab/util/resizing_array_sink.hpp>
#include <graphlab/logger/assertions.hpp>

using namespace graphlab;

struct a{
  size_t i;
  a(): i(rand()) { }
  void save(oarchive& oarc) const {
    oarc << i;
  }
  void load(iarchive& iarc)  {
    iarc >> i;
  }
};

struct b{
  size_t i;
  b(): i(rand()) { }
  void save(oarchive& oarc) const {
    oarc << i;
  }
  void load(iarchive& iarc)  {
    iarc >> i;
  }
};

SERIALIZABLE_POD(b);

int main(int argc, char** argv) {
  timer ti;
  ti.start();
  char c[10];
  int64_t u2;
  

  ti.start();
  size_t totallen = 0;
  int64_t rangelow = -100000000;
  int64_t rangehigh = 100000000;
  std::cout << "compressing and decompressing range [" << rangelow << ", "<< rangehigh << "]" << std::endl;
  for (int64_t u = rangelow;u <= rangehigh; ++u) {
    unsigned char len = compress_int(u, c);
    totallen += len;
    decompress_int(c + 10 - len, u2);
  }
  std::cout << "runtime: " << ti.current_time() << std::endl;
  std::cout << "Range compression rate: " << (8 * (rangehigh - rangelow) / ti.current_time())/(1024 * 1024) << " MBps" << std::endl;
  std::cout << "Compression Ratio: " << (double)totallen / (8 * (rangehigh - rangelow)) << std::endl;


  
    
  size_t testlen = 20000000;
  std::vector<a> testvec(testlen);
  std::vector<b> testvecpod(testlen);  
  for (size_t i =0;i < testvecpod.size(); ++i) testvecpod[i].i = testvec[i].i;
  std::stringstream strm;
  charstream chstream(0);
  
  {  
    ti.start();
    oarchive oarc(strm);
    oarc << testvec;
    std::cout << "Write to strstream at " << testlen * sizeof(a) / ti.current_time() / (1024 * 1024) << " MBps" << std::endl;
  }


  {  
    ti.start();
    iarchive iarc(strm);
    iarc >> testvec;
    std::cout << "Read from strstream at " << testlen * sizeof(a) / ti.current_time()/ (1024 * 1024) << " MBps" << std::endl;
  }

  
  strm.clear();

  
  {  
    ti.start();
    oarchive oarc(chstream);
    oarc << testvec;
    chstream.flush();
    std::cout << "Write to charstream at " << testlen * sizeof(a) / ti.current_time() / (1024 * 1024) << " MBps" << std::endl;
  }


  { 
    boost::iostreams::stream<boost::iostreams::array_source> as(chstream->c_str(), chstream->size());
    ti.start();
    
    iarchive iarc(as);
    iarc >> testvec;
    std::cout << "Read from charstream at " << testlen * sizeof(a) / ti.current_time()/ (1024 * 1024) << " MBps" << std::endl;
  }
  
  {
    ti.start();
    oarchive oarc(strm);
    oarc << testvecpod;
    std::cout << "Thunked write to strstream at " << testlen * sizeof(size_t) / ti.current_time() / (1024 * 1024) << " MBps" << std::endl;
  }


  {
    ti.start();
    iarchive iarc(strm);
    iarc >> testvecpod;
    std::cout << "Thunked read from strstream at " << testlen * sizeof(size_t) / ti.current_time() /(1024 * 1024) << " MBps" << std::endl;
  }
  
  {  
    ti.start();
    std::ofstream fout("test.bin", std::ios::binary);
    oarchive oarc(fout);
    oarc << testvec;
    fout.close();
    std::cout << "Write to file at " << testlen * sizeof(a) / ti.current_time() / (1024 * 1024) << " MBps" << std::endl;
  }


  {  
    ti.start();
    std::ifstream fin("test.bin", std::ios::binary);
    iarchive iarc(fin);
    iarc >> testvec;
    fin.close();
    std::cout << "Read from file at " << testlen * sizeof(a) / ti.current_time()/ (1024 * 1024) << " MBps" << std::endl;
  }
  

  {
    ti.start();
    std::ofstream fout("test2.bin", std::ios::binary);
    oarchive oarc(fout);
    oarc << testvecpod;
    fout.close();
    std::cout << "Thunked write to file at " << testlen * sizeof(size_t) / ti.current_time() / (1024 * 1024) << " MBps" << std::endl;
  }


  {
    ti.start();
    std::ifstream fin("test2.bin", std::ios::binary);
    iarchive iarc(fin);
    iarc >> testvecpod;
    fin.close();
    std::cout << "Thunked read from file at " << testlen * sizeof(size_t) / ti.current_time() /(1024 * 1024) << " MBps" << std::endl;
  }

}
