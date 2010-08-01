#ifndef DC_INTERNAL_HPP
#define DC_INTERNAL_HPP

#include <graphlab/distributed/dc_packet_headers.hpp>

namespace graphlab {
struct dc_thread_local_struct{
  char* sendbuffer;
  size_t sendbufferlen;
};

dc_thread_local_struct& get_thread_dc_buffer();
}

#endif