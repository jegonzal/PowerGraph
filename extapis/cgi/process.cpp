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
 * @file process.cpp
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 */
 
#include "process.hpp"

using namespace graphlab;
namespace io = boost::asio;

/////////////////////////////// INSTANCE MEMBERS ///////////////////////////////

process::
  process() : ios(), pout(ios) {
  
  if (executable.empty()) return;
  
  int pipefd[2];
  CHECK (!::pipe(pipefd));    // assert pipe created
  
  pid_t pid = ::fork();
  CHECK (0 <= pid);           // assert child created
  
  if (0 == pid){              // child process
    CHECK (!::close(pipefd[1]));
    CHECK (0 <= ::dup2(pipefd[0], 0));
    CHECK (!::close(pipefd[0]));
    CHECK (0 <= execve(executable.c_str(), NULL, NULL));
  } /* child process goes no further than here */

  // parent process
  CHECK (!::close(pipefd[0]));
  pout.assign(pipefd[1]);
  return;
  
}

process::~process(){
  // TODO: tell child process to terminate
  // close if child is still alive
  if (pout.is_open()) pout.close();
}

std::size_t process::write(const std::string str){
  
  if (!pout.is_open()){
    std::cerr << "Pipe closed unexpectedly." << std::endl;
    return 0;
  }

  boost::system::error_code ec;
  std::size_t bytes = io::write(pout, io::buffer(str), ec);
  if (ec) std::cerr << boost::system::system_error(ec).what() << std::endl;
  
  return bytes;
  
}

///////////////////////////////// CLASS MEMBERS ////////////////////////////////

const size_t process::PROC_ID = 2;
std::string process::executable = "";

void process::set_executable(const std::string path){
  executable = path;
}

process& process::get_process(){
       
  if (!thread::contains(PROC_ID)) {
    // store process in thread-local storage
    process *p = new process;
    thread::get_local(PROC_ID) = p;
    thread::set_thread_destroy_callback(detach_process);   
  }
  
  // return the process associated with the current thread
  process *p =  thread::get_local(PROC_ID).as<process *>();
  return *p;
  
}

void process::detach_process(){
  if (!thread::contains(PROC_ID)) return; // nothing to do
  process *p =  thread::get_local(PROC_ID).as<process *>();
  if (p) delete p;
  thread::get_local(PROC_ID) = NULL;
}
