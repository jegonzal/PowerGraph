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


#ifndef __STITCH_OPTS_HPP__
#define __STITCH_OPTS_HPP__

#include <string>

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


/////////////////////////////////////////////////////////////////////////
// Option Struct
struct Options
{
    // graphlab options
    std::string exec_type;
   
    // input output dirs
    std::string output_dir;

    int verbose;
   
    bool try_gpu;
    // size of images
    double work_megapix;
    double seam_megapix;
    double compose_megapix;
    double output_megapix;

    double work_scale;
    double seam_scale;
    double compose_scale;
    double output_scale;
   
    double seam_work_aspect;
    double compose_seam_aspect;
    double compose_work_aspect;
 //   double output_work_aspect;
   
    double warped_image_scale;
    std::string warp_type;

    // match options
    double conf_thresh;
    float match_conf;

    // seam options
    std::string seam_find_type;
    float terminal_cost;
    float bad_region_penalty; 

    //wave correction options
    std::string wave_correct_type;

    //bundle adjustment options
    std::string ba_cost_func;
    std::string ba_refine_mask;

    //gain compensation options
    std::string expose_comp_type;

    //blending options
    std::string blending_type;
    float blend_strength;

    //saving output
    std::string result_name;
        
   //saving the adjacency list for creating the graph
   //std::string graph_name;
   
    // Default values
    Options():
    exec_type("async"),
    output_dir("./"),
    verbose(0),
    try_gpu(false),
    work_megapix(0.6), seam_megapix(0.1), compose_megapix(-1), output_megapix(0.6),
    work_scale(1), seam_scale(1), compose_scale(1), output_scale(1),
    seam_work_aspect(1/6), compose_seam_aspect(1), compose_work_aspect(1),
    warped_image_scale(-1), warp_type("spherical"),
    conf_thresh(1.f), match_conf(0.3f),
    seam_find_type("gc_color"), terminal_cost(10000.f), bad_region_penalty(1000.f),
    wave_correct_type("horiz"),
    ba_cost_func("ray"),
    ba_refine_mask("xxxxx"),
    expose_comp_type("gain_blocks"),
    blending_type("multiband"), blend_strength(5),
    result_name("result_stitch.jpg")
    {}
};

// output_megapix(1), output_scale(1), 
extern Options opts;

#endif

