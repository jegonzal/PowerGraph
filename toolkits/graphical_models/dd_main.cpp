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

#include "dd_main.hpp"
#include "dd_grlab.hpp"



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
    
    std::string graph_file;
    std::string output_dir = "pred";
    std::string exec_type = "async";
    Options opts;
    
    clopts.attach_option("graph", graph_file,
                         "The path to UAI file containing the factors");
    clopts.add_positional("graph");
    clopts.attach_option("output", output_dir,
                         "The directory in which to save the predictions");
    clopts.add_positional("output");
    clopts.attach_option("dualimprovthres", opts.dualimprovthres,
                         "The tolerance level for Dual Convergence.");
    clopts.attach_option("pdgapthres", opts.pdgapthres,
                         "The tolerance level for Primal-Dual Gap.");
    clopts.attach_option("maxiter", opts.pdgapthres,
                         "The maximum no. of DD iterations.");
    clopts.attach_option("verbose", opts.verbose,
                         "Verbosity of Printing: 0 (default, no printing) or 1 (lots).");
    clopts.attach_option("engine", exec_type,
                         "The type of engine to use {async, sync}.");
    
    if(!clopts.parse(argc, argv)) 
    {
        graphlab::mpi_tools::finalize();
        return clopts.is_set("help")? EXIT_SUCCESS : EXIT_FAILURE;
    }
    
    ///! display settings  
    if(dc.procid() == 0) 
    {
        std::cout 
        << "ncpus:          " << clopts.get_ncpus() << std::endl
        << "scheduler:      " << clopts.get_scheduler_type() << std::endl
        ;
    }
    
        
    ///! load the graph
    graph_type graph(dc, clopts);  
    
    
    ///! load the graph
    //    graph.load(prior_dir, vertex_loader);
    //    graph.load(graph_dir, edge_loader);
    loadUAIfile(dc, graph, graph_file, opts);
    graph.finalize();
    
#if 0    
    // Test subproblem_map.
    cout << graph.num_vertices() << endl;
    for (size_t i = 0; i < graph.num_vertices(); ++i)
    {
        subproblem_map(graph.vertex(i));
    }
#endif
    
    /* Nothing execute right now. So no engine.
    typedef graphlab::omni_engine<bp_vertex_program> engine_type;
    engine_type engine(dc, graph, exec_type, clopts);
    engine.signal_all();
    graphlab::timer timer;
    
    // The main command. Run graphlab
    engine.start();  
    
    const double runtime = timer.current_time();    
    dc.cout() 
    << "----------------------------------------------------------" << std::endl
    << "Final Runtime (seconds):   " << runtime 
    << std::endl
    << "Updates executed: " << engine.num_updates() << std::endl
    << "Update Rate (updates/second): " 
    << engine.num_updates() / runtime << std::endl;
    */
    
    graphlab::mpi_tools::finalize();
    return EXIT_SUCCESS;
    
    
} // end of main





















