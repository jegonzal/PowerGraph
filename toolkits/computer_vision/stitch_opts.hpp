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

/////////////////////////////////////////////////////////////////////////
// Option Struct
struct Options 
{
    // graphlab options
    std::string exec_type;
    
    // input output dirs
    std::string output_dir;

    int verbose;
    
    // size of images
    double work_megapix;
    double seam_megapix;
    double compose_megapix;

    double work_scale;
    double seam_scale;
    double compose_scale;
    
    double seam_work_aspect;
    double compose_seam_aspect;
    double compose_work_aspect;
    
    double warped_image_scale;

    // match options
    double conf_thresh;

    // seam options
    std::string seam_find_type;
    float terminal_cost;
    float bad_region_penalty;    
    
    // Default values
    Options(): 
    exec_type("async"),
    output_dir("./"),
    verbose(0), 
    work_megapix(0.6), seam_megapix(0.1), compose_megapix(-1),
    work_scale(1), seam_scale(1), compose_scale(1),
    seam_work_aspect(1/6), compose_seam_aspect(1), compose_work_aspect(1),
    warped_image_scale(-1),
    conf_thresh(1.0),
    seam_find_type("gc_color"), terminal_cost(10000.f), bad_region_penalty(1000.f)
    {}
};

extern Options opts;

#endif
