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
    clopts.attach_option("stepsize_type", opts.stepsize_type, "to specify type of stepsize, 0 for no variation, 1 for 1/t type, 2 for polyak");
    clopts.add_positional("stepsize_type");
    clopts.attach_option("dualimprovthres", opts.dualimprovthres,
                         "The tolerance level for Dual Convergence.");
    clopts.attach_option("pdgapthres", opts.pdgapthres,
                         "The tolerance level for Primal-Dual Gap.");
    clopts.attach_option("maxiter", opts.maxiter,
                         "The maximum no. of DD iterations.");
    clopts.attach_option("verbose", opts.verbose,
                         "Verbosity of Printing: 0 (default, no printing) or 1 (lots).");
    clopts.attach_option("engine", opts.exec_type,
                         "The type of engine to use {async, sync}.");
    
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
    loadUAIfile(dc, graph, opts.graph_file);
    graph.finalize();
    
    // Define the engine.
    typedef graphlab::omni_engine<dd_vertex_program_symmetric> engine_type;
    //typedef graphlab::omni_engine<dd_vertex_program_projected> engine_type;

    // Instantiate the engine object
    engine_type engine(dc, graph, opts.exec_type, clopts);
    engine.signal_all();
    graphlab::timer timer;
    
    // Attach an aggregator to compute primal/dual objective, periodic with 0.5 s intervals.
    
     engine.add_vertex_aggregator<objective>("pd_obj",sum, print_obj); 
     engine.aggregate_periodic("pd_obj",0.0001);

    // The main command. Run graphlab
    engine.start();  
    // engine.aggregate_now("pd_obj");
    
    const double runtime = timer.current_time();    
    dc.cout() 
    << "----------------------------------------------------------" << std::endl
    << "Final Runtime (seconds):   " << runtime 
    << std::endl
    << "Updates executed: " << engine.num_updates() << std::endl
    << "Update Rate (updates/second): " 
    << engine.num_updates() / runtime << std::endl;
    
    if ( opts.history_file.size() < 4 || opts.history_file.find(".txt",opts.history_file.size()-4) == std::string::npos)
    {opts.history_file.append(".txt");} 
    char *filename = (char*)opts.history_file.c_str();
    ofstream file;
    file.open(filename);
    int i = 0;
    while(i< history[0].size())
    { file<<history[0][i]<<" "<<history[1][i]<<" "<<history[2][i]<<endl;
      i++;}
     file.close();
    graphlab::mpi_tools::finalize();
    return EXIT_SUCCESS;
    
    
} // end of main
