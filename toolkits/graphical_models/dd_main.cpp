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
 *  \authors Dhruv Batra, André Martins, Aroma Mahendru
 */

#include "dd_main.hpp"
#include "dd_grlab.hpp"
#include "ad3_qp.hpp"

Options opts;

/////////////////////////////////////////////////////////////////////////
// Main function
int main(int argc, char** argv) 
{
    global_logger().set_log_level(LOG_INFO);
    global_logger().set_log_to_console(true);
    
    ///! Initialize control plain using mpi
    graphlab::mpi_tools::init(argc, argv);
    graphlab::distributed_control dc;
    
    
    // Parse command line options -----------------------------------------------
    const std::string description = "Dual-Decomposition for MAP-MRF Inference";
    graphlab::command_line_options clopts(description);
    
    clopts.attach_option("graph", opts.graph_file,
                         "The path to UAI file containing the factors");
    clopts.add_positional("graph");
    clopts.attach_option("output", opts.output_dir,
                         "The directory in which to save the predictions");
    clopts.add_positional("output");
    clopts.attach_option("history", opts.history_file, " for saving objective values");
    clopts.add_positional("history");
    clopts.attach_option("step_size", opts.step_size, "initial/fixed stepsize");
    clopts.add_positional("step_size");
    clopts.attach_option("debug_mode", opts.debug, "to activate debug mode set it to true");
    clopts.add_positional("debug_mode");
    clopts.attach_option("dualimprovthres", opts.dualimprovthres,
                         "The tolerance level for Dual Convergence.");
     clopts.add_positional("dualimprovthres");                    
    //clopts.attach_option("pdgapthres", opts.pdgapthres,
    //                     "The tolerance level for Primal-Dual Gap.");
    clopts.attach_option("maxiter", opts.maxiter,
                         "The maximum no. of DD iterations.");
    clopts.attach_option("verbose", opts.verbose,
                         "Verbosity of Printing: 0 (default, no printing) or 1 (lots).");
    clopts.attach_option("engine", opts.exec_type,
                         "The type of engine to use {async, sync}.");
    clopts.attach_option("algorithm", opts.algorithm, 
                         "specify type of algorithm: 0 for dd_symmetric, 1 for dd_projected, 2 for admm");
    clopts.add_positional("algorithm");
    clopts.attach_option("format", opts.file_format, 
                         "specify file format : uai or distr_uai");
    clopts.add_positional("format");
    clopts.attach_option("agg_time", opts.agg_time, 
                         "specify the time period after aggregator works");
    clopts.add_positional("agg_time");
    clopts.attach_option("alpha", opts.alpha, 
                         "specify the value of parameter alpha for bethe admm");
    clopts.add_positional("alpha");
    
    if(!clopts.parse(argc, argv)) 
    {
        graphlab::mpi_tools::finalize();
        return clopts.is_set("help")? EXIT_SUCCESS : EXIT_FAILURE;
    }
    
    if(opts.graph_file.empty()) 
    {
        logstream(LOG_ERROR) << "No adjacency file provided." << std::endl;
        return EXIT_FAILURE;
    }

   
    ///! display settings  
    if(dc.procid() == 0) 
    {
        std::cout 
        << "ncpus:          " << clopts.get_ncpus() << std::endl
        << "engine:         " << opts.exec_type << std::endl
        << "scheduler:      " << clopts.get_scheduler_type() << std::endl
        << "graph_file:     " << opts.graph_file << std::endl
        << "verbose:        " << opts.verbose << std::endl
        ;
    }
    
       
    // Instantiate graph object
    graph_type graph(dc, clopts);  
    
    
    // load the graph
    //    graph.load(prior_dir, vertex_loader);
    //    graph.load(graph_dir, edge_loader);
    int nnodes  = 0;
    if(opts.file_format == "uai")
    loadUAIfile(dc, graph, opts.graph_file, nnodes);
    else
    graph.load(opts.graph_file.c_str(), line_parser);

    graph.finalize();
  
    // run dual decomposition    
    switch (opts.algorithm){
    case 1: run_dd_projected(dc, graph, opts.exec_type, clopts);
            break;
    case 2: run_ad3(dc, graph, opts.exec_type, clopts);
            break;
    case 3: run_bethe_admm(dc, graph, opts.exec_type, clopts);
            break;
    default : run_dd_symmetric(dc, graph, opts.exec_type, clopts);
    }
    
    // save predictions
    if(opts.output_dir.find("/", opts.output_dir.size()-1) == std::string::npos){
    opts.output_dir.append("/"); }
    opts.output_dir.append("output.txt");
    graph.save(opts.output_dir.c_str(), graph_writer(), false, true, false,1); 

    // save history 
    if ( opts.history_file.size() < 4 || opts.history_file.find(".txt",opts.history_file.size()-4) == std::string::npos)
    {opts.history_file.append(".txt");} 
    char *filename = (char*)opts.history_file.c_str();
    ofstream file;
    file.open(filename);
    int i = 0;
    while(i< global_vars.history[0].size())
    { file<<global_vars.history[0][i]<<" "<<global_vars.history[1][i]<<" "
          <<global_vars.history[2][i]<<" "<<global_vars.history[3][i]<<endl;
     i++;}
    file.close();
    
    graphlab::mpi_tools::finalize();
    return EXIT_SUCCESS;
    
    
} // end of main
