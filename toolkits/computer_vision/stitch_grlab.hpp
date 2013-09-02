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


#ifndef __STITCH_GRLAB_HPP__
#define __STITCH_GRLAB_HPP__


#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <math.h>
#include <unistd.h>


#include <Eigen/Dense>

#include <graphlab.hpp>
#include <graphlab/util/fs_util.hpp>


#include "eigen_serialization.hpp"
#include "opencv_serialization.hpp"
#include "stitch_opts.hpp"
#include "utils.hpp"
//#include "seam_finders_gr.hpp"
#include "precomp.hpp"

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
    bool empty; // used to quickly check if this is a dummy vertex.
   
    // path to image
    std::string img_path;
   
    cv::Mat full_img;       // Original image
    cv::Mat img;            // Used for feature computation
    cv::Mat img_warped;     // Used by gain compensator
    cv::Mat img_warped_f;   // Used by seam_finder
   
    cv::Size full_img_size;

    cv:: Size warp_size;
        
    cv::detail::ImageFeatures features;
   
    cv::detail::CameraParams camera;
   
    cv::Point2f corner;
    //cv::Mat mask;
    cv::Mat mask_warped;
   
    // constructor
    vertex_data() : empty(true)
    { }
   
    void save(graphlab::oarchive& arc) const
    {
        arc << empty << img_path
        << full_img << img << img_warped << img_warped_f
        << full_img_size << warp_size << features << camera
        << corner
        << mask_warped;
    }
    void load(graphlab::iarchive& arc)
    {
        arc >> empty >> img_path
        >> full_img >> img >> img_warped >> img_warped_f
        >> full_img_size >> warp_size >> features >> camera
        >> corner 
        >> mask_warped;
    }
   //>> mask
    vertex_data operator+ (vertex_data& othervertex)
    {
        vertex_data sum;
       
        if (!empty && !othervertex.empty)
        {
            logstream(LOG_ERROR) << "Don't know about to merge two non-empty vertex-data structures" << std::endl;
            //return EXIT_FAILURE;
        }
        else if (!empty && othervertex.empty)
            sum = *this;
        else if (empty && !othervertex.empty)
            sum = othervertex;        
        // Nothing to do if both empty.
       
        return sum;
    }

    vertex_data& operator+= (const vertex_data& othervertex)
    {
        if (!empty && !othervertex.empty)
        {
            logstream(LOG_ERROR) << "Don't know about to merge two non-empty vertex-data structures" << std::endl;
            //return EXIT_FAILURE;
        }
        else if (empty && !othervertex.empty)
            *this = othervertex;        
        // Nothing to do if both empty or othervertex empty.
       
        return *this;
    }
}; // End of vertex data


//typedef graphlab::empty edge_data;
struct edge_data
{
    bool empty; // used to quickly check if this is a dummy edge.

    cv::detail::MatchesInfo matchinfo;
   
    // constructor
    edge_data() : empty(true)
    { }
   
    void save(graphlab::oarchive& arc) const
    {
        arc << empty << matchinfo;
    }
    void load(graphlab::iarchive& arc)
    {
        arc >> empty >> matchinfo;
    }

    edge_data operator+ (edge_data& otheredge)
    {
        edge_data sum;
       
        if (!empty && !otheredge.empty)
        {
            logstream(LOG_ERROR) << "Don't know about to merge two non-empty edge-data structures" << std::endl;
            //return EXIT_FAILURE;
        }
        else if (!empty && otheredge.empty)
            sum = *this;
        else if (empty && !otheredge.empty)
            sum = otheredge;        
        // Nothing to do if both empty.
       
        return sum;
    }

    edge_data& operator+= (const edge_data& otheredge)
    {
        if (!empty && !otheredge.empty)
        {
            logstream(LOG_ERROR) << "Don't know about to merge two non-empty edge-data structures" << std::endl;
            //return EXIT_FAILURE;
        }
        else if (empty && !otheredge.empty)
            *this = otheredge;        
        // Nothing to do if both empty or othervertex empty.
       
        return *this;
    }
}; // End of edge data

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
        return gather_type();
    } // end of gather

    void apply(icontext_type& context, vertex_type& vertex,
               const gather_type& sum)
    {
       
        // Get vertex data
        vertex_data &vdata = vertex.data();

        logstream(LOG_EMPH) << "Features in image #" << vertex.id() << ": " << vdata.features.keypoints.size() << "\n";
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
//typedef graphlab::async_consistent_engine<stitch_vertex_program> engine_type;
typedef graphlab::omni_engine<stitch_vertex_program> engine_type;

#endif

