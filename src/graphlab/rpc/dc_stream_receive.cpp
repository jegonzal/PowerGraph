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

#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_internal_types.hpp>
#include <graphlab/rpc/dc_stream_receive.hpp>

//#define DC_RECEIVE_DEBUG
namespace graphlab {
namespace dc_impl {



char* dc_stream_receive::get_buffer(size_t& retbuflength) {
  header_read = 0;
  retbuflength = sizeof(block_header_type);
  return reinterpret_cast<char*>(&cur_chunk_header);
}


char* dc_stream_receive::advance_buffer(char* c, size_t wrotelength, 
                            size_t& retbuflength) {
  
  if (header_read != sizeof(block_header_type)) {
    // tcp is still writing into cur`writelen
    header_read += wrotelength;
    ASSERT_LE(header_read, sizeof(block_header_type));
    // are we done reading the header?
    if (header_read < sizeof(block_header_type)) {
      // nope!
      retbuflength = sizeof(block_header_type) - header_read;
      return (reinterpret_cast<char*>(&cur_chunk_header) + header_read);
    }
    else {
      // ok now lets switch it to the write buffer
      ASSERT_TRUE(writebuffer == NULL);
      writebuffer = (char*)malloc(cur_chunk_header);
      retbuflength = cur_chunk_header;
      write_buffer_written = 0;
      return writebuffer;
    }
  }
  else {
    // we read the entire header and is reading buffers now
    // try to store the buffer and see if we are full yet.
    ASSERT_EQ(header_read, sizeof(block_header_type));
    write_buffer_written += wrotelength;
    if (write_buffer_written < cur_chunk_header) {
      retbuflength = cur_chunk_header - write_buffer_written;
      return writebuffer + write_buffer_written;
    }
  }

  // if we reach here, we have an available block
  // give away the buffer to dc
  dc->deferred_function_call_chunk(writebuffer, cur_chunk_header, associated_proc);
  writebuffer = NULL;
  write_buffer_written = 0;
  header_read = 0;
  return get_buffer(retbuflength);
}


  
void dc_stream_receive::shutdown() { }

} // namespace dc_impl
} // namespace graphlab

