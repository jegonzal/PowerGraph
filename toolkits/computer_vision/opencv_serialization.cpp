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
} 

graphlab::iarchive& operator>>(graphlab::iarchive& arc, cv::Size& img_size) 
{
    arc >> img_size.width >> img_size.height;  
    return arc;
} 


//////////////////////////////////////////////////
// For Point2f
graphlab::oarchive& operator<<(graphlab::oarchive& arc, const cv::Point2f& pt) 
{
    arc << pt.x << pt.y;
    return arc;
} 

graphlab::iarchive& operator>>(graphlab::iarchive& arc, cv::Point2f& pt) 
{
    arc >> pt.x >> pt.y;
    return arc;
} 


//////////////////////////////////////////////////
// For KeyPoint
graphlab::oarchive& operator<<(graphlab::oarchive& arc, const cv::KeyPoint& keypoint) 
{
    arc << keypoint.pt 
    << keypoint.size << keypoint.angle << keypoint.response
    << keypoint.octave << keypoint.class_id;
    return arc;
} 

graphlab::iarchive& operator>>(graphlab::iarchive& arc, cv::KeyPoint& keypoint) 
{
    arc >> keypoint.pt 
    >> keypoint.size >> keypoint.angle >> keypoint.response
    >> keypoint.octave >> keypoint.class_id;
    return arc;
} 


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
} 

graphlab::iarchive& operator>>(graphlab::oarchive& arc, const cv::Mat& mat) 
{
    int cols, rows; size_t elem_size, elem_type;
    
    arc >> cols >> rows 
    >> elem_size >> elem_type;
    
    mat.create(rows, cols, elem_type);
    
    size_t data_size = mat.cols * mat.rows * elem_size;
    graphlab::deserialize(arc, mat.ptr(), data_size);
    return arc;
} 


//////////////////////////////////////////////////
// For ImageFeatures
graphlab::oarchive& operator<<(graphlab::oarchive& arc, const cv::detail::ImageFeatures& features) 
{
    arc << features.img_idx << features.img_size
    << features.keypoints
    << features.descriptors;    
    return arc;
} 


graphlab::iarchive& operator>>(graphlab::oarchive& arc, const cv::detail::ImageFeatures& features) 
{
    arc >> features.img_idx >> features.img_size
    >> features.keypoints
    >> features.descriptors;
    return arc;
} 


//////////////////////////////////////////////////
// For DMatch
graphlab::oarchive& operator<<(graphlab::oarchive& arc, const cv::DMatch& match) 
{
    arc << match.queryIdx << match.trainIdx << match.imgIdx
    << match.distance;     
    return arc;
}


graphlab::iarchive& operator>>(graphlab::oarchive& arc, const cv::DMatch& match) 
{
    arc >> match.queryIdx >> match.trainIdx >> match.imgIdx
    >> match.distance;     
    return arc;
}


//////////////////////////////////////////////////
// For MatchesInfo 
graphlab::oarchive& operator<<(graphlab::oarchive& arc, const cv::detail::MatchesInfo& matchesinfo) 
{
    arc << matchesinfo.src_img_idx << matchesinfo.dst_img_idx
    << matchesinfo.matches 
    << matchesinfo.inliers_mask << matchesinfo.num_inliers
    << matchesinfo.H 
    << matchesinfo.confidence;
    return arc;
}


graphlab::iarchive& operator>>(graphlab::oarchive& arc, const cv::detail::MatchesInfo& matchesinfo) 
{
    arc >> matchesinfo.src_img_idx >> matchesinfo.dst_img_idx
    >> matchesinfo.matches 
    >> matchesinfo.inliers_mask >> matchesinfo.num_inliers
    >> matchesinfo.H 
    >> matchesinfo.confidence;
    return arc;
}


//////////////////////////////////////////////////
// For CameraParams
graphlab::oarchive& operator<<(graphlab::oarchive& arc, const cv::detail::CameraParams& camera) 
{
    arc << camera.focal << camera.aspect 
    << camera.ppx << camera.ppy
    << camera.R << camera.t;
    return arc;
} 


graphlab::iarchive& operator>>(graphlab::oarchive& arc, const cv::detail::CameraParams& camera) 
{
    arc >> camera.focal >> camera.aspect 
    >> camera.ppx >> camera.ppy
    >> camera.R >> camera.t;
    return arc;
} 
