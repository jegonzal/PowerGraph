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

#include "opencv_serialization.hpp"



//////////////////////////////////////////////////
// For Size
graphlab::oarchive& operator<<(graphlab::oarchive& arc, const cv::Size& img_size) 
{
    arc << img_size.width << img_size.height;  
    return arc;
} // end of save vector

graphlab::iarchive& operator>>(graphlab::iarchive& arc, cv::Size& img_size) 
{
    arc >> img_size.width >> img_size.height;  
    return arc;
} // end of save vector


//////////////////////////////////////////////////
// For Point2f
graphlab::oarchive& operator<<(graphlab::oarchive& arc, const cv::Point2f& pt) 
{
    arc << pt.x << pt.y;
    return arc;
} // end of save vector

graphlab::iarchive& operator>>(graphlab::iarchive& arc, cv::Point2f& pt) 
{
    arc >> pt.x >> pt.y;
    return arc;
} // end of save vector


//////////////////////////////////////////////////
// For KeyPoint
graphlab::oarchive& operator<<(graphlab::oarchive& arc, const cv::KeyPoint& keypoint) 
{
    arc << keypoint.pt 
    << keypoint.size << keypoint.angle << keypoint.response
    << keypoint.octave << keypoint.class_id;
    return arc;
} // end of save vector

graphlab::iarchive& operator>>(graphlab::iarchive& arc, cv::KeyPoint& keypoint) 
{
    arc >> keypoint.pt 
    >> keypoint.size >> keypoint.angle >> keypoint.response
    >> keypoint.octave >> keypoint.class_id;
    return arc;
} // end of save vector


//////////////////////////////////////////////////
// For Mat
graphlab::oarchive& operator<<(graphlab::oarchive& arc, const cv::Mat& mat) 
{
    size_t elem_size = mat.elemSize();
    size_t elem_type = mat.type();
    
    arc << mat.cols << mat.rows 
    << elem_size << elem_type;
    
    const size_t data_size = mat.cols * mat.rows * elem_size;
    graphlab::serialize(arc, mat.ptr(), data_size);
    return arc;
} // end of save vector

graphlab::iarchive& operator>>(graphlab::oarchive& arc, const cv::Mat& mat) 
{
    int cols, rows; size_t elem_size, elem_type;
    
    arc >> cols >> rows 
    >> elem_size >> elem_type;
    
    mat.create(rows, cols, elem_type);
    
    size_t data_size = mat.cols * mat.rows * elem_size;
    graphlab::deserialize(arc, mat.ptr(), data_size);
    return arc;
} // end of save vector
