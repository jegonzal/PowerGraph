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

#ifndef OPENCV_SERIALIZATION_HPP
#define OPENCV_SERIALIZATION_HPP

#include <graphlab.hpp>

#include "opencv2/opencv_modules.hpp"
#include "opencv2/opencv.hpp"


//////////////////////////////////////////////////
// For Size
BEGIN_OUT_OF_PLACE_SAVE(arc, cv::Size, img_size) 
{
    arc << img_size.width << img_size.height;
} END_OUT_OF_PLACE_SAVE()


BEGIN_OUT_OF_PLACE_LOAD(arc, cv::Size, img_size) 
{
    arc >> img_size.width >> img_size.height;
} END_OUT_OF_PLACE_LOAD()


//////////////////////////////////////////////////
// For Point2f
BEGIN_OUT_OF_PLACE_SAVE(arc, cv::Point2f, pt) 
{
    arc << pt.x << pt.y;
} END_OUT_OF_PLACE_SAVE()


BEGIN_OUT_OF_PLACE_LOAD(arc, cv::Point2f, pt) 
{
    arc >> pt.x >> pt.y;
} END_OUT_OF_PLACE_LOAD()

//////////////////////////////////////////////////
// For KeyPoint
BEGIN_OUT_OF_PLACE_SAVE(arc, cv::KeyPoint, keypoint) 
{
    arc << keypoint.pt 
    << keypoint.size << keypoint.angle << keypoint.response
    << keypoint.octave << keypoint.class_id;
} END_OUT_OF_PLACE_SAVE()


BEGIN_OUT_OF_PLACE_LOAD(arc, cv::Size, img_size) 
{
    arc >> keypoint.pt 
    >> keypoint.size >> keypoint.angle >> keypoint.response
    >> keypoint.octave >> keypoint.class_id;
} END_OUT_OF_PLACE_LOAD()



#endif
