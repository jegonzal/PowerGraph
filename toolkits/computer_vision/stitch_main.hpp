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
 *
 * \brief This file contains an example of graphlab used for stitching
 * multiple images into a panorama. The code is based on a example
 * stiching application in OpenCV.
 *
 *  \author Dhruv Batra
 */


#ifndef __STITCH_MAIN_HPP__
#define __STITCH_MAIN_HPP__

#include "utils.hpp"
#include "stitch_grlab.hpp"


/////////////////////////////////////////////////////////////////////////
// Vertex Loader (used to read images and load the vertex data of the graph)
//bool vertex_loader(graph_type& graph, const std::string& fname, 
//                   const std::string& line) 
bool vertex_loader(graphlab::distributed_control& dc, graph_type& graph, string img_path)
{ 
    // force a "/" at the end of the path
    // make sure to check that the path is non-empty. (you do not
    // want to make the empty path "" the root path "/" )
    string path = img_path;
    if (path.length() > 0 && path[path.length() - 1] != '/') path = path + "/";
    
    vector<string> graph_files;
    string search_prefix;
    graphlab::fs_util::list_files_with_prefix(path, search_prefix, graph_files);
    
    if (graph_files.size() == 0) 
        logstream(LOG_WARNING) << "No files found in " << path << std::endl;
    
    // vertex data & id
    graphlab::vertex_id_type vid(-1);
    
    ///////////////////////////////////////////////////////
    // Loop over files
    for(size_t i = 0; i < graph_files.size(); ++i) 
    {
        // Each machine loads corresponding file
        if (i % dc.numprocs() == dc.procid()) 
        {
            if (opts.verbose > 0)
                logstream(LOG_EMPH) 
                << "Process: " << dc.procid() << "/" << dc.numprocs() << " "
                << "picked image: " << graph_files[i] << "\n";
            
            
            vid = i;
            vertex_data vdata;
            vdata.empty = false;
            vdata.img_path = graph_files[i];
            vdata.features.img_idx = i;
            
            graph.add_vertex(vid, vdata);
            
        }
    }
    
    return true;
}


/////////////////////////////////////////////////////////////////////////
// Edge Loader (used to read the adjacency list and add edges to the graph)
bool edge_loader(graph_type& graph, const std::string& fname, 
                   const std::string& textline) 
{
    if ( textline.length() == 0 || textline[0] == '#' )
        return true; // empty or comment line, return
    
    std::stringstream strm(textline);
    graphlab::vertex_id_type vid;
        
    // first entry in the line is a vertex ID
    strm >> vid;
    
    if (opts.verbose > 0)
        logstream(LOG_EMPH) << "Here's the input: "
        << textline << "\n"
        << vid << "\n";
    
    // Line should contain at least 1 more number (degree of node)
    if (!strm.good())
    {
        logstream(LOG_ERROR) << "The following ajacency list line is incomplete(check adj_list standard):\n" 
        << textline << std::endl;
        return EXIT_FAILURE;
    }
        
    // second entry is the out-degree
    int outdeg;
    strm >> outdeg;
    
    graphlab::vertex_id_type other_vid;
    for (int i=0; i!=outdeg; ++i) 
    {
        // Line should contain more numbers (id of neighbours)
        if (!strm.good())
        {
            logstream(LOG_ERROR) << "The following ajacency list line is incomplete(check adj_list standard):\n" 
            << textline << std::endl;
            return EXIT_FAILURE;
        }
        
        strm >> other_vid;
        
        // only add edges in one direction
        if (other_vid < vid)
            continue;
        
        if (opts.verbose > 0)
            logstream(LOG_EMPH) << "Adding edge: (" << vid << "," << other_vid << ")\n";
        
        edge_data edata; edata.empty = false;
        graph.add_edge(vid,other_vid,edata);
    }
    
    return true;
}

/////////////////////////////////////////////////////////////////////////
// Function to extract features in parallel
void compute_features(graph_type::vertex_type vertex)
{
    // Get vertex data
    vertex_data &vdata = vertex.data();
    
    // Load image
    //            // open the stream
    //            std::ifstream in_file(graph_files[i].c_str(), 
    //                                  std::ios_base::in | std::ios_base::binary);
    //
    //            boost::iostreams::filtering_stream<boost::iostreams::input> fin;  
    //            fin.push(in_file);
    //            
    //            // Get data from stream into a buffer
    //            fin.pop();
    
    // Ignore the above hdfs-setup for now. Just read from file directly.
    Mat full_img = imread(vdata.img_path);
    Mat img;
    
    // Scale image if necessary
    double work_scale = min(1.0, sqrt(opts.work_megapix * 1e6 / full_img.size().area()));
    
    if ( abs(work_scale-1) > 1e-3 )
        resize(full_img, img, Size(), work_scale, work_scale);
    else
        img = full_img;
    
    if (img.empty())
    {
        logstream(LOG_EMPH) << "Could not imread image: " << vdata.img_path << "\n";
        //exit();
        //return EXIT_FAILURE;
    }
    
    // compute features
    SurfFeaturesFinder finder;
    finder(img, vdata.features);
    finder.collectGarbage();
    
    if (opts.verbose > 0)
        logstream(LOG_EMPH) << "Features in image #" << vertex.id() << ": " << vdata.features.keypoints.size() << "\n";
    
}

/////////////////////////////////////////////////////////////////////////
// Function to compute feature-matches in parallel on edges
void match_features(graph_type::edge_type edge)
{
    // Get edge data
    edge_data &edata = edge.data();
    
    // Get vertex ids of two vertices involved
    vertex_data &vdata1 = edge.source().data();
    vertex_data &vdata2 = edge.target().data();
    
    // Match features 
    BestOf2NearestMatcher matcher;
    matcher(vdata1.features, vdata2.features, edata.matchinfo);
    matcher.collectGarbage();
    
    if (opts.verbose > 0)
        logstream(LOG_EMPH) << "#Matches in Image Pair "
        "(" << edge.source().id() << "," << edge.target().id() << ")" 
        << ": (" << edata.matchinfo.matches.size() 
        << "," << edata.matchinfo.num_inliers << ")"
        << "\n";
    
//    // Estimate focal length from Homography
//    double f0, f1; bool f0ok, f1ok;
//    focalsFromHomography(edata.matchinfo.H, f0, f1, f0ok, f1ok);
}

/////////////////////////////////////////////////////////////////////////
// Map Function to compile a list of features
//vector<vertex_data> compile_features(const graph_type::vertex_type& vertex)
vector<vertex_data> compile_features(engine_type::icontext_type& context, 
                         const graph_type::vertex_type& vertex)
{
    vector<vertex_data> temp(context.num_vertices());
    
    temp[vertex.id()] = vertex.data();
    return temp;
}

/////////////////////////////////////////////////////////////////////////
// Map Function to compile a list of matches
//vector<vertex_data> compile_features(const graph_type::vertex_type& vertex)
vector<edge_data> compile_matches(engine_type::icontext_type& context, 
                                     const graph_type::edge_type& edge)
{
    
    int edlist_len = context.num_vertices() * context.num_vertices();
    vector<edge_data> temp(edlist_len);
    
    int pair_idx = edge.source().id() * context.num_vertices() + edge.target().id();
    temp[pair_idx] = edge.data();
    return temp;
}

#endif
