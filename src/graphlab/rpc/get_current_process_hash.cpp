#include <graphlab/rpc/get_current_process_hash.hpp>
#include <graphlab/ui/mongoose/mongoose.h>

#ifdef __APPLE__
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <libproc.h>
#endif
namespace graphlab {
namespace dc_impl {


#ifdef __linux
std::string get_current_process_hash() {
  char buf[33];
  mg_md5_file(buf, "/proc/self/exe");
  buf[32] = '\0';
  std::string ret = buf;
  if (ret.length() != 32) {
    ret = std::string(32, '0');
  }
  return ret;
}
#elif __APPLE_
std::string get_current_process_hash() {
  std::string ret;

  pid_t pid = getpid();
  char pathbuf[PROC_PIDPATHINFO_MAXSIZE];
  int pidsuccess = proc_pidpath (pid, pathbuf, sizeof(pathbuf));
  if (pidsuccess > 0) {
    char buf[33];
    mg_md5_file(buf,  pathbuf);
    buf[32] = '\0';
    ret = buf;
  }
  if (ret.length() != 32) {
    ret = std::string(32, '0');
  }
  return ret;
}
#endif

} // dc_impl
} // graphlab
