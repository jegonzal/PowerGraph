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
 
/**
 * @file process.hpp
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 */
 
#ifndef GRAPHLAB_PROCESS_HPP
#define GRAPHLAB_PROCESS_HPP

#include <iostream>
#include <boost/asio.hpp>
#include <graphlab.hpp>

namespace graphlab {
//////////////////////////////// PROCESS CLASS ////////////////////////////////

  class json_message;

  /**
   * Process
   * Handles the spawning and killing of children.
   */
  class process {
  
  public:
    typedef boost::asio::posix::stream_descriptor stream_descriptor;
    typedef boost::asio::io_service io_service;
    typedef boost::asio::streambuf streambuf;
  private:
  
    /** path to the executable */
    static std::string executable;
    
    io_service ios;
    
    /** output stream to write to process */
    stream_descriptor pout;
    
    /** input stream to read from process */
    stream_descriptor pin;
    
    /** private constructor */
    process();
    
    /** Redirects i/o to pipes - used only in constructor */
    void redirect_io();
    
    /** Writes output to process */
    std::size_t write(json_message& message);
    
    /**
     * Reads input from process, delimited by JSON braces.
     */
    json_message& read(json_message& message);
    
  public:
    
    virtual ~process();
    
    std::size_t send(json_message& message);
    json_message& receive(json_message& message);
    
    /** ID of process in thread local store */
    static const size_t PROC_ID;
  
    /** Sets the path to the executable */
    static void set_executable(const std::string path);
  
    /**
     * Retrieves the associated process for the current thread.
     * If the process can be found in the thread-local store, returns
     * immediately; otherwise, that means that the current thread does not have
     * an associated process yet. In that case, this function will create a new
     * process and store it in the thread-local store.
     * @return process associated with the current thread.
     */
    static process& get_process();
    
    /**
     * Detaches the process from the thread.
     * If a pointer to the process cannot be found in the thread-local
     * store, that means that this thread has already been detached, and the
     * function will return immediately. Otherwise, the thread is detached and
     * the pointer to the process is removed from the thread-local store.
     */
    static void detach_process();
  
  }; // end of process

////////////////////////////////////////////////////////////////////////////////
}

#endif /* #ifndef GRAPHLAB_PROCESS_HPP */
