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
 * This file contains an example of graphlab used for stitching
 * multiple images into a panorama. The code is based on a example
 * stiching application in OpenCV.
 *
 *  \author Dhruv Batra
 */


#include <vector>
#include <string>
#include <fstream>


#include <Eigen/Dense>

#include <graphlab.hpp>
#include <graphlab/util/fs_util.hpp>


#include "eigen_serialization.hpp"

#include <graphlab/macros_def.hpp>

#include "opencv2/opencv_modules.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/stitching/detail/autocalib.hpp"
#include "opencv2/stitching/detail/blenders.hpp"
#include "opencv2/stitching/detail/camera.hpp"
#include "opencv2/stitching/detail/exposure_compensate.hpp"
#include "opencv2/stitching/detail/matchers.hpp"
#include "opencv2/stitching/detail/motion_estimators.hpp"
#include "opencv2/stitching/detail/seam_finders.hpp"
#include "opencv2/stitching/detail/util.hpp"
#include "opencv2/stitching/detail/warpers.hpp"
#include "opencv2/stitching/warpers.hpp"

using namespace std;
using namespace cv;
using namespace cv::detail;

typedef Eigen::VectorXd vec;
typedef Eigen::MatrixXd mat;


/////////////////////////////////////////////////////////////////////////
// Edge and Vertex data and Graph Type
struct vertex_data 
{
    // path to image
    std::string img_path;
    
    //cv::Mat img;
    cv::detail::ImageFeatures features;
    
    // constructor
    vertex_data()
    { }
    
    void save(graphlab::oarchive& arc) const 
    {
        arc << img_path;
        //arc << features;
        //arc << feautres.descriptors;
    }
    void load(graphlab::iarchive& arc) 
    {
        arc >> img_path;
        //arc >> features;
    }
}; // End of vertex data


typedef graphlab::empty edge_data;

/**
 * The graph type
 */
typedef graphlab::distributed_graph<vertex_data, edge_data> graph_type;


/////////////////////////////////////////////////////////////////////////
// GraphLab Vertex Program (and Gather Type)
/**
 * The type passed around during the gather phase
 */
//struct gather_type 
//{
//    
//    gather_type& operator+=(const gather_type& other) 
//    {
//    } // end of operator +=
//    void save(graphlab::oarchive& arc) const 
//    {
//        arc << delf_i << delf_j;
//    }
//    void load(graphlab::iarchive& arc) 
//    {
//        arc >> delf_i >> delf_j;
//    }
//}; // end of gather type
typedef graphlab::empty gather_type;

/** 
 * The core stitching update function.  
 */
class stitch_vertex_program : 
public graphlab::ivertex_program<graph_type, gather_type, 
graphlab::messages::sum_priority>,
public graphlab::IS_POD_TYPE 
{
private:
    
public:
    
    stitch_vertex_program()  { }
    
    edge_dir_type gather_edges(icontext_type& context,
                               const vertex_type& vertex) const 
    { 
        return graphlab::ALL_EDGES; 
    }; // end of gather_edges 

    // Run the gather operation over all in edges
    gather_type gather(icontext_type& context, const vertex_type& target_vertex, 
                       edge_type& edge) const 
    {
    } // end of gather

    void apply(icontext_type& context, vertex_type& vertex, 
               const gather_type& sum) 
    {
        
    } // end of apply

    edge_dir_type scatter_edges(icontext_type& context,
                                const vertex_type& vertex) const 
    { 
        //return graphlab::ALL_EDGES;
        return graphlab::NO_EDGES;
    }; // end of gather_edges 
    
};


/**
 * Define the engine type
 */
//typedef graphlab::synchronous_engine<mplp_vertex_program> engine_type;
typedef graphlab::async_consistent_engine<stitch_vertex_program> engine_type;


/////////////////////////////////////////////////////////////////////////
// Vertex Loader (used to read load the vertex data of the graph)
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
            logstream(LOG_EMPH) 
            //<< "Process: " << dc.procid() << "/" << dc.numprocs()
            << "picked image: " << graph_files[i] << "\n";

            
            vid = i;
            vertex_data vdata;            
            vdata.img_path = graph_files[i];
            vdata.features.img_idx = i;
            
            graph.add_vertex(vid, vdata);

        }
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
    Mat img = imread(vdata.img_path);

    if (img.empty())
    {
        logstream(LOG_EMPH) << "Could not imread image: " << vdata.img_path << "\n";
        //exit();
        //return EXIT_FAILURE;
    }
    else
        logstream(LOG_EMPH) << "Done\n";
    
    // compute features
    SurfFeaturesFinder finder;
    finder(img, vdata.features);

    
    logstream(LOG_EMPH) << "Features in image #" << vertex.id() << ": " << vdata.features.keypoints.size() << "\n";
    
}