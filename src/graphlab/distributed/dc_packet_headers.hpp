#ifndef DISTRIBUTED_CONTROL_PACKET_HEADERS_HPP
#define DISTRIBUTED_CONTROL_PACKET_HEADERS_HPP

#include <graphlab/distributed/distributed_control.hpp>
namespace graphlab {

enum packet_ids{
  REMOTECALL_ID, REMOTECALLX_ID, REMOTECALLXS_ID, REMOTECALL_CONTROL_ID
};

struct remotecall_packdata{
  uint16_t packtype;
  procid_t srcnodeid;
  void* fnptr;
  uint32_t len;
  uint32_t numargs;
  handlerarg_t args[0];
};

struct remotecallx_packdata{
  uint16_t packtype;
  procid_t srcnodeid;
  void* fnptr;
  uint32_t len;
  uint32_t stacklen;
};

struct remotecallxs_packdata{
  uint16_t packtype;
  procid_t srcnodeid;
  void* fnptr;
  uint32_t len;
  uint32_t stacklen;
};

struct fence_packdata{
  uint16_t packtype;
};

typedef uint32_t pheaderlen_t;
}


#endif