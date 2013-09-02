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
#include "opencv2/stitching/stitcher.hpp"


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


BEGIN_OUT_OF_PLACE_LOAD(arc, cv::KeyPoint, keypoint) 
{
    arc >> keypoint.pt 
    >> keypoint.size >> keypoint.angle >> keypoint.response
    >> keypoint.octave >> keypoint.class_id;
} END_OUT_OF_PLACE_LOAD()


//////////////////////////////////////////////////
// For Mat
BEGIN_OUT_OF_PLACE_SAVE(arc, cv::Mat, mat) 
{
    size_t elem_size = mat.elemSize();
    size_t elem_type = mat.type();
    
    arc << mat.cols << mat.rows 
    << elem_size << elem_type;
    
    const size_t data_size = mat.cols * mat.rows * elem_size;
    graphlab::serialize(arc, mat.ptr(), data_size);
} END_OUT_OF_PLACE_SAVE()


BEGIN_OUT_OF_PLACE_LOAD(arc, cv::Mat, mat) 
{
    int cols, rows; size_t elem_size, elem_type;

    arc >> cols >> rows 
    >> elem_size >> elem_type;
    
    mat.create(rows, cols, elem_type);
    
    size_t data_size = mat.cols * mat.rows * elem_size;
    graphlab::deserialize(arc, mat.ptr(), data_size);    
} END_OUT_OF_PLACE_LOAD()


//////////////////////////////////////////////////
// For ImageFeatures
BEGIN_OUT_OF_PLACE_SAVE(arc, cv::detail::ImageFeatures, features) 
{
    arc << features.img_idx << features.img_size
    << features.keypoints
    << features.descriptors;
} END_OUT_OF_PLACE_SAVE()


BEGIN_OUT_OF_PLACE_LOAD(arc, cv::detail::ImageFeatures, features) 
{
    arc >> features.img_idx >> features.img_size
    >> features.keypoints
    >> features.descriptors;
} END_OUT_OF_PLACE_LOAD()


//////////////////////////////////////////////////
// For DMatch
BEGIN_OUT_OF_PLACE_SAVE(arc, cv::DMatch, match) 
{
    arc << match.queryIdx << match.trainIdx << match.imgIdx
    << match.distance;     
} END_OUT_OF_PLACE_SAVE()


BEGIN_OUT_OF_PLACE_LOAD(arc, cv::DMatch, match) 
{
    arc >> match.queryIdx >> match.trainIdx >> match.imgIdx
    >> match.distance;     
} END_OUT_OF_PLACE_LOAD()


//////////////////////////////////////////////////
// For MatchesInfo 
BEGIN_OUT_OF_PLACE_SAVE(arc, cv::detail::MatchesInfo, matchesinfo) 
{
    arc << matchesinfo.src_img_idx << matchesinfo.dst_img_idx
    << matchesinfo.matches 
    << matchesinfo.inliers_mask << matchesinfo.num_inliers
    << matchesinfo.H 
    << matchesinfo.confidence;    
} END_OUT_OF_PLACE_SAVE()


BEGIN_OUT_OF_PLACE_LOAD(arc, cv::detail::MatchesInfo, matchesinfo) 
{
    arc >> matchesinfo.src_img_idx >> matchesinfo.dst_img_idx
    >> matchesinfo.matches 
    >> matchesinfo.inliers_mask >> matchesinfo.num_inliers
    >> matchesinfo.H 
    >> matchesinfo.confidence;
} END_OUT_OF_PLACE_LOAD()


//////////////////////////////////////////////////
// For CameraParams
BEGIN_OUT_OF_PLACE_SAVE(arc, cv::detail::CameraParams, camera) 
{
    arc << camera.focal << camera.aspect 
    << camera.ppx << camera.ppy
    << camera.R << camera.t;
} END_OUT_OF_PLACE_SAVE()


BEGIN_OUT_OF_PLACE_LOAD(arc, cv::detail::CameraParams, camera) 
{
    arc >> camera.focal >> camera.aspect 
    >> camera.ppx >> camera.ppy
    >> camera.R >> camera.t;
} END_OUT_OF_PLACE_LOAD()

#endif
