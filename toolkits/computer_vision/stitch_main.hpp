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
#include "iostream"

/////////////////////////////////////////////////////////////////////////
// Graph Loader (used to read images and load the vertex and edge data of
// the graph. There is an edge between every pair of vertices in the graph,  
// i.e. the graph is fully connected. (Used in stitch_full_main)

bool graph_loader(graphlab::distributed_control& dc, graph_type& graph, string img_dir)
{
    // force a "/" at the end of the path
    // make sure to check that the path is non-empty. (you do not
    // want to make the empty path "" the root path "/" )
    string path = img_dir;
    if (path.length() > 0 && path[path.length() - 1] != '/') path = path + "/";
    
    vector<string> graph_files;
    string search_prefix;
    graphlab::fs_util::list_files_with_prefix(path, search_prefix, graph_files);
    
    if (graph_files.size() == 0)
        logstream(LOG_WARNING) << "No files found in " << path << std::endl;
    
    
    if (opts.verbose > 2)
        logstream(LOG_EMPH)
        << "Total number of images: " << graph_files.size() << "\n";
    
    // vertex data & id
    graphlab::vertex_id_type vid(-1);
    graphlab::vertex_id_type other_vid;
    
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
            
            if (opts.verbose > 2)
                logstream(LOG_EMPH)
                << "Vertex " << i << " Image: " << vdata.img_path << "\n"; 
            
        }
    }
    
    // Adding edges between every pair of vertices to create a fully connected graph
    // no duplicate edges are added
    for(size_t i = 0; i < graph_files.size()-1; ++i)
    {
        vid = i;
        for(size_t j = i+1; j < graph_files.size(); ++j)
        {
            other_vid = j;
            if (opts.verbose > 0)
                logstream(LOG_EMPH) << "Adding edge: (" << vid << "," << other_vid << ")\n";
            
            edge_data edata; edata.empty = false;
            graph.add_edge(vid,other_vid,edata);
        }
    }
    return true;
}

// Second Graph loader that only a single machine calls and pre-loads cameras.
// It adds only the selected edges to the subgraph of the fully connected graph.
// (Used in stitch_full_main)
bool graph_loader(graph_type& graph, string img_dir, vector<CameraParams>& cameras, 
                  vector<string>& img_path, vector<int> indices, vector<MatchesInfo>& pairwise_matches)
{
    // force a "/" at the end of the path
    // make sure to check that the path is non-empty. (you do not
    // want to make the empty path "" the root path "/" )
    string path = img_dir;
    if (path.length() > 0 && path[path.length() - 1] != '/') path = path + "/";
    
    // vertex data & id   
    graphlab::vertex_id_type vid(-1);
    graphlab::vertex_id_type other_vid;
    
    if (opts.verbose > 2)
        logstream(LOG_EMPH)
        << "Total number of vertices in the second graph: " << indices.size() << "\n";
    
    ///////////////////////////////////////////////////////
    // Loop over files
    for(size_t i = 0; i < indices.size(); ++i)
    {
        vid = i;
        vertex_data vdata;
        vdata.empty = false;
        vdata.img_path = img_path[i];
        vdata.features.img_idx = i;
        
        vdata.camera = cameras[i]; // addition to above function.
        
        graph.add_vertex(vid, vdata);
    }
    
    // Adding edges between selected pair of vertices to create a subgraph of the fully connected graph
    // no duplicate edges are added
    
    if (opts.verbose > 2)
        logstream(LOG_EMPH)
        << "Pairwise_matches size : " << pairwise_matches.size() << "\n";
    
    for(size_t i = 0; i < pairwise_matches.size(); ++i)
    {
        if (opts.verbose > 2)
            logstream(LOG_EMPH)
            << "Pairwise_matches : (" << pairwise_matches[i].src_img_idx << "," << pairwise_matches[i].dst_img_idx << ")\n";
        
        if (pairwise_matches[i].src_img_idx >= 0 && pairwise_matches[i].dst_img_idx >= 0)
        {
            if (pairwise_matches[i].src_img_idx < pairwise_matches[i].dst_img_idx)	// no duplicate edges are allowed
            {
                vid = pairwise_matches[i].src_img_idx;
                other_vid = pairwise_matches[i].dst_img_idx;
                if (opts.verbose > 0)
            	    logstream(LOG_EMPH) << "Adding edge: (" << vid << "," << other_vid << ")\n";
                
                edge_data edata; edata.empty = false;
                graph.add_edge(vid,other_vid,edata);
            }
        }
    }
    
    return true;
}

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

// Second loader that only a single machine calls and pre-loads cameras.
bool vertex_loader(graph_type& graph, string img_path, vector<CameraParams>& cameras)
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
        vid = i;
        vertex_data vdata;
        vdata.empty = false;
        vdata.img_path = graph_files[i];
        vdata.features.img_idx = i;
        
        vdata.camera = cameras[i]; // addition to above function.
        
        graph.add_vertex(vid, vdata);
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
// Map-Aggregator function to find the largest image size
struct ImgArea
{
    double full_img_area;
    
    void save(graphlab::oarchive& arc) const
    {
        arc << full_img_area;
    }
    void load(graphlab::iarchive& arc)
    {
        arc >> full_img_area;
    }
    
    ImgArea operator+= (ImgArea other) // computes max
    {
        ImgArea max;
        max.full_img_area =
        (full_img_area > other.full_img_area) ? full_img_area : other.full_img_area;
        return max;
    }
};

ImgArea find_largest_img(engine_type::icontext_type& context,
                         graph_type::vertex_type& vertex)
{
    // Get vertex data
    vertex_data &vdata = vertex.data();
    
    // Ignore hdfs-setup for now. Just read from file directly.
    Mat full_img = imread(vdata.img_path);
    
    ImgArea imgarea;
    imgarea.full_img_area = full_img.size().area();
    
    if (opts.verbose > 2)
        LOGLN("largest image area: " << imgarea.full_img_area << "\n");
    
    return imgarea;
}

void set_scales(engine_type::icontext_type& context, ImgArea& largestimg)
{
    
    if (opts.work_megapix > 0)
        opts.work_scale = min(1.0, sqrt(opts.work_megapix * 1e6 / largestimg.full_img_area));
    
    opts.seam_scale = min(1.0, sqrt(opts.seam_megapix * 1e6 / largestimg.full_img_area));
    
    if (opts.compose_megapix > 0)
        opts.compose_scale = min(1.0, sqrt(opts.compose_megapix * 1e6 / largestimg.full_img_area));
    
    opts.seam_work_aspect = opts.seam_scale / opts.work_scale;
    opts.compose_seam_aspect = opts.compose_scale / opts.seam_scale;
    opts.compose_work_aspect = opts.compose_scale / opts.work_scale;
}


/////////////////////////////////////////////////////////////////////////
// Function to extract features in vertex parallel
void compute_features(graph_type::vertex_type& vertex)
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
    Mat &full_img = vdata.full_img;
    Mat &img = vdata.img;
    full_img = imread(vdata.img_path);
    
    if ( abs(opts.work_scale-1) > 1e-3 )
        resize(full_img, img, Size(), opts.work_scale, opts.work_scale);
    else
        img = full_img;
    
    if (img.empty())
        logstream(LOG_ERROR) << "Could not imread image: " << vdata.img_path << "\n";
    
    // compute features
    SurfFeaturesFinder finder;
    finder(img, vdata.features);
    finder.collectGarbage();
    
    if (opts.verbose > 0)
    {
        logstream(LOG_EMPH) << "Features in image #" << vertex.id() << ": " << vdata.features.keypoints.size() << "\n";
        LOGLN("Size of feature image #" << img.cols << "  " << img.rows); //
    }
    
}

/////////////////////////////////////////////////////////////////////////
// Function to compute feature-matches in parallel on edges
void match_features(graph_type::edge_type& edge)
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
}

/////////////////////////////////////////////////////////////////////////
// Function to warp images in parallel
void warp_images(graph_type::vertex_type& vertex)
{
    // Get vertex data
    vertex_data &vdata = vertex.data();
    
    Mat full_img = imread(vdata.img_path);
    
    Mat &img = vdata.img;
    Mat &img_warped = vdata.img_warped;
    Mat &img_warped_f = vdata.img_warped_f;
    //Mat &mask = vdata.mask;
    Mat mask;
    Mat &mask_warped = vdata.mask_warped;
    CameraParams &camera = vdata.camera;
    Point2f &corner = vdata.corner;
    Size &size = vdata.warp_size;
    
    if (full_img.empty())
        logstream(LOG_ERROR) << "Could not imread image: " << vdata.img_path << "\n";
    
    vdata.full_img_size = full_img.size();
    
    // Scale image if necessary
    double seam_scale = min(1.0, sqrt(opts.seam_megapix * 1e6 / full_img.size().area()));
    
    if ( abs(seam_scale-1) > 1e-3 )
        resize(full_img, img, Size(), seam_scale, seam_scale);
    else
        img = full_img.clone();
    
    // Prepare images mask
    mask.create(img.size(), CV_8U);
    mask.setTo(Scalar::all(255));
    
    if (opts.verbose > 2)
        LOGLN("Size of mask image #" << img.cols << "  " << img.rows); //
    
    // Warp images and their masks
    Ptr<WarperCreator> warper_creator;
    
    if (opts.warp_type == "plane") warper_creator = new cv::PlaneWarper();
    else if (opts.warp_type == "cylindrical") warper_creator = new cv::CylindricalWarper();
    else if (opts.warp_type == "spherical") warper_creator = new cv::SphericalWarper();
    else if (opts.warp_type == "fisheye") warper_creator = new cv::FisheyeWarper();
    else if (opts.warp_type == "stereographic") warper_creator = new cv::StereographicWarper();
    else if (opts.warp_type == "compressedPlaneA2B1") warper_creator = new cv::CompressedRectilinearWarper(2, 1);
    else if (opts.warp_type == "compressedPlaneA1.5B1") warper_creator = new cv::CompressedRectilinearWarper(1.5, 1);
    else if (opts.warp_type == "compressedPlanePortraitA2B1") warper_creator = new cv::CompressedRectilinearPortraitWarper(2, 1);
    else if (opts.warp_type == "compressedPlanePortraitA1.5B1") warper_creator = new cv::CompressedRectilinearPortraitWarper(1.5, 1);
    else if (opts.warp_type == "paniniA2B1") warper_creator = new cv::PaniniWarper(2, 1);
    else if (opts.warp_type == "paniniA1.5B1") warper_creator = new cv::PaniniWarper(1.5, 1);
    else if (opts.warp_type == "paniniPortraitA2B1") warper_creator = new cv::PaniniPortraitWarper(2, 1);
    else if (opts.warp_type == "paniniPortraitA1.5B1") warper_creator = new cv::PaniniPortraitWarper(1.5, 1);
    else if (opts.warp_type == "mercator") warper_creator = new cv::MercatorWarper();
    else if (opts.warp_type == "transverseMercator") warper_creator = new cv::TransverseMercatorWarper();
    
    
    if (warper_creator.empty())
        logstream(LOG_ERROR) << "Can't create the following warper '" << opts.warp_type << "'\n";
    
    Ptr<RotationWarper> warper = warper_creator->create(static_cast<float>(opts.warped_image_scale * opts.seam_work_aspect));
    Mat_<float> K;
    camera.K().convertTo(K, CV_32F);
    float swa = (float)opts.seam_work_aspect;
    K(0,0) *= swa; K(0,2) *= swa;
    K(1,1) *= swa; K(1,2) *= swa;
    
    corner = warper->warp(img, K, camera.R, INTER_LINEAR, BORDER_REFLECT, img_warped);
    if (opts.verbose > 2)
        LOGLN("Warp corners x : " << corner.x << "   y : " << corner.y << "\n");
    
    size = img_warped.size();
    if (opts.verbose > 2)
        LOGLN("Warp sizes height : " << size.height << "   width : " << size.width << "\n");
    
    warper->warp(mask, K, camera.R, INTER_NEAREST, BORDER_CONSTANT, mask_warped);
    
    img_warped.convertTo(img_warped_f, CV_32F);
    
    // If no gain compensator, then clear.
    img_warped.release();
    
}

/////////////////////////////////////////////////////////////////////////
// Function to composite images in parallel
void composite_images(graph_type::vertex_type& vertex)
{
    // Get vertex data
    vertex_data &vdata = vertex.data();
    CameraParams &camera = vdata.camera;
    Point2f &corner = vdata.corner;
    Mat full_img = imread(vdata.img_path);	//we have to check it later for speed
    
    Mat &img_warped = vdata.img_warped;		//added by me
    Mat &mask_warped = vdata.mask_warped;	//added by me
    Size &size = vdata.warp_size;		//added by me
    Mat mask, dilated_mask, seam_mask, masks_warped;		//added by me
    
    
    if (full_img.empty())
        logstream(LOG_ERROR) << "Could not imread image: " << vdata.img_path << "\n";
    
    Mat &img = vdata.img;
    
    // Update warped image scale
    Ptr<WarperCreator> warper_creator;
    
    if (opts.warp_type == "plane") warper_creator = new cv::PlaneWarper();
    else if (opts.warp_type == "cylindrical") warper_creator = new cv::CylindricalWarper();
    else if (opts.warp_type == "spherical") warper_creator = new cv::SphericalWarper();
    else if (opts.warp_type == "fisheye") warper_creator = new cv::FisheyeWarper();
    else if (opts.warp_type == "stereographic") warper_creator = new cv::StereographicWarper();
    else if (opts.warp_type == "compressedPlaneA2B1") warper_creator = new cv::CompressedRectilinearWarper(2, 1);
    else if (opts.warp_type == "compressedPlaneA1.5B1") warper_creator = new cv::CompressedRectilinearWarper(1.5, 1);
    else if (opts.warp_type == "compressedPlanePortraitA2B1") warper_creator = new cv::CompressedRectilinearPortraitWarper(2, 1);
    else if (opts.warp_type == "compressedPlanePortraitA1.5B1") warper_creator = new cv::CompressedRectilinearPortraitWarper(1.5, 1);
    else if (opts.warp_type == "paniniA2B1") warper_creator = new cv::PaniniWarper(2, 1);
    else if (opts.warp_type == "paniniA1.5B1") warper_creator = new cv::PaniniWarper(1.5, 1);
    else if (opts.warp_type == "paniniPortraitA2B1") warper_creator = new cv::PaniniPortraitWarper(2, 1);
    else if (opts.warp_type == "paniniPortraitA1.5B1") warper_creator = new cv::PaniniPortraitWarper(1.5, 1);
    else if (opts.warp_type == "mercator") warper_creator = new cv::MercatorWarper();
    else if (opts.warp_type == "transverseMercator") warper_creator = new cv::TransverseMercatorWarper();
    
    if (warper_creator.empty())
        logstream(LOG_ERROR) << "Can't create the following warper '" << opts.warp_type << "'\n";
    
    Ptr<RotationWarper> warper = warper_creator->create(static_cast<float>(opts.warped_image_scale * opts.compose_work_aspect));
    
    // Update intrinsics
    camera.focal *= opts.compose_work_aspect;
    camera.ppx *= opts.compose_work_aspect;
    camera.ppy *= opts.compose_work_aspect;
    
    // Update corner and size
    //vdata.full_img_size = full_img.size();
    Size sz = vdata.full_img_size;
    if (std::abs(opts.compose_scale - 1) > 1e-1)
    {
        sz.width = cvRound(vdata.full_img_size.width * opts.compose_scale);
        sz.height = cvRound(vdata.full_img_size.height * opts.compose_scale);
    }
    
    Mat K;
    camera.K().convertTo(K, CV_32F);
    Rect roi = warper->warpRoi(sz, K, camera.R);
    corner = roi.tl();
    if (opts.verbose > 2)
        LOGLN("Compose corner x : " << corner.x << "   y : " << corner.y << "\n");
    
    size = roi.size();
    if (opts.verbose > 2)
        LOGLN("Compose size height : " << size.height << "   width : " << size.width << "\n");
    
    if (abs(opts.compose_scale - 1) > 1e-1)
        resize(full_img, img, Size(), opts.compose_scale, opts.compose_scale);
    else
        img = full_img;
    Size img_size = img.size();
    
    // Warp the current image
    warper->warp(img, K, camera.R, INTER_LINEAR, BORDER_REFLECT, img_warped);
    
    // Warp the current image mask
    mask.create(img_size, CV_8U);
    mask.setTo(Scalar::all(255));
    warper->warp(mask, K, camera.R, INTER_NEAREST, BORDER_CONSTANT, masks_warped);
    
    // Compensate exposure
    //compensator->apply(img_idx, corner[img_idx], img_warped, mask_warped);
    
    img.release();
    mask.release();
    
    dilate(mask_warped, dilated_mask, Mat());
    resize(dilated_mask, seam_mask, masks_warped.size());
    mask_warped = seam_mask & masks_warped;
    
}

/////////////////////////////////////////////////////////////////////////
// Function to compute feature-matches in parallel on edges
void find_seams(graph_type::edge_type& edge)
{
    // Get edge data
    //edge_data &edata = edge.data(); //commented by me as it was unused
	
    // Get vertex ids of two vertices involved
    vertex_data &vdata1 = edge.source().data();
    vertex_data &vdata2 = edge.target().data();
    
    // Not sure why this is needed anymore?
    //Ptr<SeamFinder> seam_finder;
    //seam_finder = new detail::GraphCutSeamFinder(GraphCutSeamFinderBase::COST_COLOR);
    
    
    // Code from PairwiseSeamFinder::Impl::findInPair()
    //Mat img1 = images_[first], img2 = images_[second];
    Mat &img1 = vdata1.img_warped_f; Mat &img2 = vdata2.img_warped_f;
    vector<Mat> src; src.push_back(img1); src.push_back(img2);
    
    vector<Mat> dx_(2), dy_(2);  
    Mat dx, dy;
    for (size_t i = 0; i < src.size(); ++i)
    {
        CV_Assert(src[i].channels() == 3);
        Sobel(src[i], dx, CV_32F, 1, 0);
        Sobel(src[i], dy, CV_32F, 0, 1);
        dx_[i].create(src[i].size(), CV_32F);
        dy_[i].create(src[i].size(), CV_32F);
        for (int y = 0; y < src[i].rows; ++y)
        {
            const Point3f* dx_row = dx.ptr<Point3f>(y);
            const Point3f* dy_row = dy.ptr<Point3f>(y);
            float* dx_row_ = dx_[i].ptr<float>(y);
            float* dy_row_ = dy_[i].ptr<float>(y);
            for (int x = 0; x < src[i].cols; ++x)
            {
                dx_row_[x] = normL2(dx_row[x]);
                dy_row_[x] = normL2(dy_row[x]);
            }
        }
    }
    
    //Mat dx1 = dx_[first], dx2 = dx_[second];
    //Mat dy1 = dy_[first], dy2 = dy_[second];
    Mat &dx1 = dx_[0]; Mat &dx2 = dx_[1];
    Mat &dy1 = dy_[0]; Mat &dy2 = dy_[1];
    
    //Mat mask1 = masks_[first], mask2 = masks_[second];
    Mat &mask1 = vdata1.mask_warped; Mat &mask2 = vdata2.mask_warped;
    //Point tl1 = corners_[first], tl2 = corners_[second];
    Point2f &tl1 = vdata1.corner; Point2f &tl2 = vdata2.corner;
    
    Rect roi;
    overlapRoi(tl1, tl2, img1.size(), img2.size(), roi);
    
    const int gap = 10;
    Mat subimg1(roi.height + 2 * gap, roi.width + 2 * gap, CV_32FC3);
    Mat subimg2(roi.height + 2 * gap, roi.width + 2 * gap, CV_32FC3);
    Mat submask1(roi.height + 2 * gap, roi.width + 2 * gap, CV_8U);
    Mat submask2(roi.height + 2 * gap, roi.width + 2 * gap, CV_8U);
    Mat subdx1(roi.height + 2 * gap, roi.width + 2 * gap, CV_32F);
    Mat subdy1(roi.height + 2 * gap, roi.width + 2 * gap, CV_32F);
    Mat subdx2(roi.height + 2 * gap, roi.width + 2 * gap, CV_32F);
    Mat subdy2(roi.height + 2 * gap, roi.width + 2 * gap, CV_32F);
    
    // Cut subimages and submasks with some gap
    for (int y = -gap; y < roi.height + gap; ++y)
    {
        for (int x = -gap; x < roi.width + gap; ++x)
        {
            int y1 = roi.y - tl1.y + y;
            int x1 = roi.x - tl1.x + x;
            if (y1 >= 0 && x1 >= 0 && y1 < img1.rows && x1 < img1.cols)
            {
                subimg1.at<Point3f>(y + gap, x + gap) = img1.at<Point3f>(y1, x1);
                submask1.at<uchar>(y + gap, x + gap) = mask1.at<uchar>(y1, x1);
                subdx1.at<float>(y + gap, x + gap) = dx1.at<float>(y1, x1);
                subdy1.at<float>(y + gap, x + gap) = dy1.at<float>(y1, x1);
            }
            else
            {
                subimg1.at<Point3f>(y + gap, x + gap) = Point3f(0, 0, 0);
                submask1.at<uchar>(y + gap, x + gap) = 0;
                subdx1.at<float>(y + gap, x + gap) = 0.f;
                subdy1.at<float>(y + gap, x + gap) = 0.f;
            }
            
            int y2 = roi.y - tl2.y + y;
            int x2 = roi.x - tl2.x + x;
            if (y2 >= 0 && x2 >= 0 && y2 < img2.rows && x2 < img2.cols)
            {
                subimg2.at<Point3f>(y + gap, x + gap) = img2.at<Point3f>(y2, x2);
                submask2.at<uchar>(y + gap, x + gap) = mask2.at<uchar>(y2, x2);
                subdx2.at<float>(y + gap, x + gap) = dx2.at<float>(y2, x2);
                subdy2.at<float>(y + gap, x + gap) = dy2.at<float>(y2, x2);
            }
            else
            {
                subimg2.at<Point3f>(y + gap, x + gap) = Point3f(0, 0, 0);
                submask2.at<uchar>(y + gap, x + gap) = 0;
                subdx2.at<float>(y + gap, x + gap) = 0.f;
                subdy2.at<float>(y + gap, x + gap) = 0.f;
            }
        }
    }
    
    const int vertex_count = (roi.height + 2 * gap) * (roi.width + 2 * gap);
    const int edge_count = (roi.height - 1 + 2 * gap) * (roi.width + 2 * gap) + (roi.width - 1 + 2 * gap) * (roi.height + 2 * gap);
    GCGraph<float> graph(vertex_count, edge_count);	
    
    
    const Size img_size = subimg1.size();
    
    if (opts.seam_find_type.compare("gc_color") ==0)
    {
    	// Set terminal weights
    	for (int y = 0; y < img_size.height; ++y)
    	{
            for (int x = 0; x < img_size.width; ++x)
            {
                int v = graph.addVtx();
                graph.addTermWeights(v, submask1.at<uchar>(y, x) ? opts.terminal_cost : 0.f, 
                                     submask2.at<uchar>(y, x) ? opts.terminal_cost : 0.f);
            }
    	}
        
    	// Set regular edge weights
    	const float weight_eps = 1.f;
        
        for (int y = 0; y < img_size.height; ++y)
    	{
            for (int x = 0; x < img_size.width; ++x)
            {
                int v = y * img_size.width + x;
                
                if (x < img_size.width - 1)
                {
                    float weight = normL2(subimg1.at<Point3f>(y, x), subimg2.at<Point3f>(y, x)) +
                    normL2(subimg1.at<Point3f>(y, x + 1), subimg2.at<Point3f>(y, x + 1)) + weight_eps;
                    
                    if (!submask1.at<uchar>(y, x) || !submask1.at<uchar>(y, x + 1) || 
                        !submask2.at<uchar>(y, x) || !submask2.at<uchar>(y, x + 1))
                        weight += opts.bad_region_penalty;
                    
                    graph.addEdges(v, v + 1, weight, weight);
            	}
            	if (y < img_size.height - 1)
            	{
                    float weight = normL2(subimg1.at<Point3f>(y, x), subimg2.at<Point3f>(y, x)) +
                    normL2(subimg1.at<Point3f>(y + 1, x), subimg2.at<Point3f>(y + 1, x)) + weight_eps;
                    
                    if (!submask1.at<uchar>(y, x) || !submask1.at<uchar>(y + 1, x) ||
                    	!submask2.at<uchar>(y, x) || !submask2.at<uchar>(y + 1, x))
                        weight += opts.bad_region_penalty;
                    
                    graph.addEdges(v, v + img_size.width, weight, weight);
            	}
       	    }
        }
    }
    
    
    else if (opts.seam_find_type.compare("gc_colorgrad") ==0)
    {
        // Set terminal weights
    	for (int y = 0; y < img_size.height; ++y)
    	{
            for (int x = 0; x < img_size.width; ++x)
            {
                int v = graph.addVtx();
                graph.addTermWeights(v, submask1.at<uchar>(y, x) ? opts.terminal_cost : 0.f, 
                                     submask2.at<uchar>(y, x) ? opts.terminal_cost : 0.f);
            }
        }
        
        // Set regular edge weights
        const float weight_eps = 1.f;
        for (int y = 0; y < img_size.height; ++y)
        {
            for (int x = 0; x < img_size.width; ++x)
            {
                int v = y * img_size.width + x;
                if (x < img_size.width - 1)
                {
                    float grad = subdx1.at<float>(y, x) + subdx1.at<float>(y, x + 1) +
                    subdx2.at<float>(y, x) + subdx2.at<float>(y, x + 1) + weight_eps;
                    float weight = (normL2(subimg1.at<Point3f>(y, x), subimg2.at<Point3f>(y, x)) +
                                    normL2(subimg1.at<Point3f>(y, x + 1), subimg2.at<Point3f>(y, x + 1))) / grad + weight_eps;
                    if (!submask1.at<uchar>(y, x) || !submask1.at<uchar>(y, x + 1) ||
                        !submask2.at<uchar>(y, x) || !submask2.at<uchar>(y, x + 1))
                        weight += opts.bad_region_penalty;
                    
                    graph.addEdges(v, v + 1, weight, weight);
                }
                if (y < img_size.height - 1)
                {
                    float grad = subdy1.at<float>(y, x) + subdy1.at<float>(y + 1, x) + 
                    subdy2.at<float>(y, x) + subdy2.at<float>(y + 1, x) + weight_eps;
                    float weight = (normL2(subimg1.at<Point3f>(y, x), subimg2.at<Point3f>(y, x)) + 
                                    normL2(subimg1.at<Point3f>(y + 1, x), subimg2.at<Point3f>(y + 1, x))) / grad + weight_eps;
                    
                    if (!submask1.at<uchar>(y, x) || !submask1.at<uchar>(y + 1, x) ||
                        !submask2.at<uchar>(y, x) || !submask2.at<uchar>(y + 1, x))
                        weight += opts.bad_region_penalty;
                    
                    graph.addEdges(v, v + img_size.width, weight, weight);
                }
            }
        }
    }
    
    else
        CV_Error(CV_StsBadArg, "unsupported pixel similarity measure");
    
    graph.maxFlow();
    
    for (int y = 0; y < roi.height; ++y)
    {
        for (int x = 0; x < roi.width; ++x)
        {
            if (graph.inSourceSegment((y + gap) * (roi.width + 2 * gap) + x + gap))
            {
                if (mask1.at<uchar>(roi.y - tl1.y + y, roi.x - tl1.x + x))
                    mask2.at<uchar>(roi.y - tl2.y + y, roi.x - tl2.x + x) = 0;
            }
            else
            {
                if (mask2.at<uchar>(roi.y - tl2.y + y, roi.x - tl2.x + x))
                    mask1.at<uchar>(roi.y - tl1.y + y, roi.x - tl1.x + x) = 0;
            }
        }
    }
}

/////////////////////////////////////////////////////////////////////////
// Map Function to compile a list of features
//vector<vertex_data> compile_features(const graph_type::vertex_type& vertex)
vector<vertex_data> compile_vertices(engine_type::icontext_type& context,
                                     const graph_type::vertex_type& vertex)
{
    vector<vertex_data> temp(context.num_vertices());
    
    temp[vertex.id()] = vertex.data();
    return temp;
}

/////////////////////////////////////////////////////////////////////////
// Map Function to compile a list of matches
//vector<vertex_data> compile_features(const graph_type::vertex_type& vertex)
vector<edge_data> compile_edges(engine_type::icontext_type& context,
                                const graph_type::edge_type& edge)
{
    
    int edlist_len = context.num_vertices() * context.num_vertices();
    vector<edge_data> temp(edlist_len);
    
    int pair_idx = edge.source().id() * context.num_vertices() + edge.target().id();
    temp[pair_idx] = edge.data();
    return temp;
}

#endif

