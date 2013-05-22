/*  
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
 *
 * \brief This application performs MAP inference on Markov Nets 
 * provided in standard UAI file format via Dual-Decomposition. 
 *
 *
 *  \author Dhruv Batra
 */


#ifndef __DD_OPTS_HPP__
#define __DD_OPTS_HPP__

#include <string>

/////////////////////////////////////////////////////////////////////////
// Option Struct
struct Options 
{
    // graphlab options
    std::string exec_type;
    
    // input output dirs
    std::string graph_file;
    std::string output_dir;
    std::string history_file;
    
    int verbose;
    
    double dualimprovthres;
    double pdgapthres;
    int maxiter;
    
    int stepsize_type;    
    
    
    // Default values
    Options(): 
    exec_type("sync"),
    output_dir("pred"),
    history_file("\0"),
    verbose(0), stepsize_type(1),
    dualimprovthres(1e-5), pdgapthres(1e-1), 
    maxiter(10000)
    {}
};

extern Options opts;

#endif
