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

#ifndef GRAPHLAB_HDFS_HPP
#define GRAPHLAB_HDFS_HPP

// Requires the hdfs library
#ifdef HAS_HADOOP
extern "C" {
  #include <hdfs.h>
}
#endif

#include <vector>
#include <boost/iostreams/stream.hpp>


#include <graphlab/logger/assertions.hpp>


namespace graphlab {

#ifdef HAS_HADOOP
  class hdfs {
  private:
    /** the primary filesystem object */
    hdfsFS filesystem;
  public:
    /** hdfs file source is used to construct boost iostreams */
    class hdfs_device {
    public: // boost iostream concepts
      typedef char                                          char_type;
      struct category : 
        public boost::iostreams::bidirectional_device_tag, 
        public boost::iostreams::multichar_tag,
        public boost::iostreams::closable_tag { };
    private:
      hdfsFS filesystem;
      
      hdfsFile file;
     
    public:
      hdfs_device() : filesystem(NULL), file(NULL) { }
      hdfs_device(const hdfs& hdfs_fs, const std::string& filename,
                  const bool write = false) :
        filesystem(hdfs_fs.filesystem) {
        ASSERT_TRUE(filesystem != NULL);
        // open the file
        const int flags = write? O_WRONLY : O_RDONLY;
        const int buffer_size = 0; // use default
        const short replication = 0; // use default
        const tSize block_size = 0; // use default;
        file = hdfsOpenFile(filesystem, filename.c_str(), flags, buffer_size,
                            replication, block_size);
      }
      //      ~hdfs_device() { if(file != NULL) close(); }

      void close(std::ios_base::openmode mode = std::ios_base::openmode() ) { 
        if(file == NULL) return;
        if(file->type == OUTPUT) {
          const int flush_error = hdfsFlush(filesystem, file);
          ASSERT_EQ(flush_error, 0);
        }
        const int close_error = hdfsCloseFile(filesystem, file);
        ASSERT_EQ(close_error, 0);
        file = NULL;
      }

      /** the optimal buffer size is 0. */
      inline std::streamsize optimal_buffer_size() const { return 0; }

      std::streamsize read(char* strm_ptr, std::streamsize n) {
        return hdfsRead(filesystem, file, strm_ptr, n);
      } // end of read
      std::streamsize write(const char* strm_ptr, std::streamsize n) {
         return hdfsWrite(filesystem, file, strm_ptr, n);
      }
      bool good() const { return file != NULL; }
    }; // end of hdfs device
    
    /**
     * The basic file type has constructor matching the hdfs device.
     */
    typedef boost::iostreams::stream<hdfs_device> fstream;

    /**
     * Open a connection to the filesystem. The default arguments
     * should be sufficient for most uses 
     */
    hdfs(const std::string& host = "default", tPort port = 0) {
      filesystem =  hdfsConnect(host.c_str(), port);
      ASSERT_TRUE(filesystem != NULL); 
    } // end of constructor

    ~hdfs() { 
      const int error = hdfsDisconnect(filesystem);
      ASSERT_EQ(error, 0);
    } // end of ~hdfs
    
    inline std::vector<std::string> list_files(const std::string& path) {
      int num_files = 0;
      hdfsFileInfo* hdfs_file_list_ptr = 
        hdfsListDirectory(filesystem, path.c_str(), &num_files);
      // copy the file list to the string array
      std::vector<std::string> files(num_files);
      for(int i = 0; i < num_files; ++i) 
        files[i] = std::string(hdfs_file_list_ptr[i].mName);
      // free the file list pointer
      hdfsFreeFileInfo(hdfs_file_list_ptr, num_files);
      return files;
    } // end of list_files

    inline static bool has_hadoop() { return true; }
    
    static hdfs& get_hdfs();
  }; // end of class hdfs
#else



  class hdfs {
  public:
    /** hdfs file source is used to construct boost iostreams */
    class hdfs_device {
    public: // boost iostream concepts
      typedef char                                          char_type;
      typedef boost::iostreams::bidirectional_device_tag    category;
    public:
      hdfs_device(const hdfs& hdfs_fs, const std::string& filename,
                  const bool write = false) { 
        logstream(LOG_FATAL) << "Libhdfs is not installed on this system." 
                             << std::endl;
      }
      void close() { }
      std::streamsize read(char* strm_ptr, std::streamsize n) {
        logstream(LOG_FATAL) << "Libhdfs is not installed on this system." 
                             << std::endl;
        return 0;
      } // end of read
      std::streamsize write(const char* strm_ptr, std::streamsize n) {
        logstream(LOG_FATAL) << "Libhdfs is not installed on this system." 
                             << std::endl;
        return 0;
      }
      bool good() const { return false; }
    }; // end of hdfs device
    
    /**
     * The basic file type has constructor matching the hdfs device.
     */
    typedef boost::iostreams::stream<hdfs_device> fstream;

    /**
     * Open a connection to the filesystem. The default arguments
     * should be sufficient for most uses 
     */
    hdfs(const std::string& host = "default", int port = 0) {
      logstream(LOG_FATAL) << "Libhdfs is not installed on this system." 
                           << std::endl;
    } // end of constructor


    
    inline std::vector<std::string> list_files(const std::string& path) {
      logstream(LOG_FATAL) << "Libhdfs is not installed on this system." 
                           << std::endl;
      return std::vector<std::string>();;
    } // end of list_files

    // No hadoop available
    inline static bool has_hadoop() { return false; }
    
    static hdfs& get_hdfs();
  }; // end of class hdfs


#endif

}; // end of namespace graphlab
#endif
