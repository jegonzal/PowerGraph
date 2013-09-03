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


#include "stitch_main.hpp"

Options opts;

/////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
   
    int64 app_start_time = getTickCount();
    ///////////////////////////////////////////////////////
    // Set up Graphlab
    global_logger().set_log_level(LOG_INFO);
    global_logger().set_log_to_console(true);

    ///! Initialize control plain using mpi
    graphlab::mpi_tools::init(argc, argv);
    graphlab::distributed_control dc;

    ///////////////////////////////////////////////////////
    // Set up OpenCV
    cv::setBreakOnError(true);

    ///////////////////////////////////////////////////////
    // Graphlab parse input
    const std::string description = "Image Stitching";
    graphlab::command_line_options clopts(description);

    string img_dir;
    string graph_path;
    
    clopts.attach_option("img", img_dir,
                         "The directory containing the images");
    clopts.add_positional("img");
    clopts.attach_option("graph", graph_path,
                         "The path to the adjacency list file (could be the prefix in case of multiple files)");
    clopts.add_positional("graph");
    clopts.attach_option("output", opts.output_dir,
                         "The directory in which to save the output");
    clopts.attach_option("verbose", opts.verbose,
                         "Verbosity of Printing: 0 (default, no printing) or 1 (lots).");
    clopts.attach_option("work_megapix", opts.work_megapix,
                         "Resolution for image registration step. The default is 0.6 Mpx.");
    clopts.attach_option("engine", opts.exec_type,
                         "The type of engine to use {async, sync}.");

    if(!clopts.parse(argc, argv))
    {
        graphlab::mpi_tools::finalize();
        return clopts.is_set("help")? EXIT_SUCCESS : EXIT_FAILURE;
    }
   
    if(img_dir.empty())
    {
        logstream(LOG_ERROR) << "No image directory was provided." << std::endl;
        return EXIT_FAILURE;
    }
   
    if(graph_path.empty())
    {
        logstream(LOG_ERROR) << "No adjacency file provided." << std::endl;
        return EXIT_FAILURE;
    }
   
    if (opts.work_megapix > 10)
    {
        logstream(LOG_ERROR) << "Inappropriate value for work_megapix." << std::endl;
        return EXIT_FAILURE;
    }
   
   
    // display settings  
    dc.cout()
    << "ncpus:          " << clopts.get_ncpus() << std::endl
    << "engine:         " << opts.exec_type << std::endl
    << "scheduler:      " << clopts.get_scheduler_type() << std::endl
    << "img_dir:        " << img_dir << std::endl
    << "graph_path:     " << graph_path << std::endl
    << "work_megapix:   " << opts.work_megapix << std::endl
    << "verbose:        " << opts.verbose << std::endl;
   
   
    ///////////////////////////////////////////////////////
    // Feature Graph
    graph_type graph_feat(dc, clopts);
       
    // load the graph
    //graph.load(img_dir, vertex_loader);
    vertex_loader(dc, graph_feat, img_dir);
    graph_feat.load(graph_path, edge_loader);
    graph_feat.finalize();
   
    ///////////////////////////////////////////////////////
    // Graphlab Engine
    engine_type engine_feat(dc, graph_feat, opts.exec_type, clopts);
   
    ///////////////////////////////////////////////////////
    // Run Aggregator to find size of largest image
    engine_feat.add_vertex_aggregator<ImgArea>("find_largest_img", find_largest_img, set_scales);
    engine_feat.aggregate_now("find_largest_img");

   
    ///////////////////////////////////////////////////////
    // Computer features in parallel on vertices
    graph_feat.transform_vertices(compute_features);

    ///////////////////////////////////////////////////////
    // Match features in parallel on edges
    graph_feat.transform_edges(match_features);

   
    //if (dc.procid()==0) {
    ///////////////////////////////////////////////////////
    // Compile features
    typedef vector<vertex_data> VecVD;
    VecVD vdlist = engine_feat.map_reduce_vertices<VecVD>(compile_vertices);
   
    vector<ImageFeatures> features(vdlist.size());
    for (size_t i=0; i!=vdlist.size(); ++i)
    {
        features[i] = vdlist[i].features;
    }
    vdlist.clear();
   
    int num_images = features.size();
   
    ///////////////////////////////////////////////////////
    // Compile matches
    typedef vector<edge_data> VecED;
    VecED edlist = engine_feat.map_reduce_edges<VecED>(compile_edges);
   
    if ((opts.verbose > 0) & (dc.procid()==0))
        logstream(LOG_EMPH) << "edlist.size() =  " << edlist.size()
        << "\n";


    vector<MatchesInfo> pairwise_matches(edlist.size());
    int r,c; int pair_idx;
    for (size_t i=0; i!=edlist.size(); ++i)
    {
        IND2SUB_RM(i,r,c,num_images)
       
        if (r==c)
            continue;
       
        if (r<c)
            pair_idx = i;
        else
            pair_idx = SUB2IND_RM(c,r,num_images);

        pairwise_matches[i] = edlist[pair_idx].matchinfo;
        pairwise_matches[i].src_img_idx = r;
        pairwise_matches[i].dst_img_idx = c;
       
        if (r>c) // Swap & invert a few things in the match
        {
            if (!pairwise_matches[i].H.empty())
                pairwise_matches[i].H = pairwise_matches[i].H.inv();
           
            for (size_t j = 0; j < pairwise_matches[i].matches.size(); ++j)
                std::swap(pairwise_matches[i].matches[j].queryIdx,
                          pairwise_matches[i].matches[j].trainIdx);
        }
       
        if ((opts.verbose > 0) & (dc.procid()==0))
            logstream(LOG_EMPH) << "#Matches in Pair "
            "(" << pairwise_matches[i].src_img_idx
            << "," << pairwise_matches[i].dst_img_idx << ")"
            << ": (" << pairwise_matches[i].matches.size()
            << "," << pairwise_matches[i].num_inliers
            << "," << pairwise_matches[i].confidence << ")"
            << "\n";

    }
    edlist.clear();

    ///////////////////////////////////////////////////////
    // Homography-Based Initialization
    int64 t;
    t = getTickCount();
    HomographyBasedEstimator estimator;
    vector<CameraParams> cameras;
    estimator(features, pairwise_matches, cameras);
    logstream(LOG_EMPH) << "Homography-based init, time: " << ((getTickCount() - t) / getTickFrequency()) << " sec\n";
   
    for (size_t i = 0; i < cameras.size(); ++i)
    {
        Mat R;
        cameras[i].R.convertTo(R, CV_32F);
        cameras[i].R = R;
        if (dc.procid() == 0)
            logstream(LOG_EMPH) << "Initial intrinsics #" << i << ":\n" << cameras[i].K() << "\n\n";
    }

   
    ///////////////////////////////////////////////////////
    // Bundle Adjustment
    t = getTickCount();
    Ptr<detail::BundleAdjusterBase> adjuster;
    //adjuster = new detail::BundleAdjusterRay();
    if (opts.ba_cost_func == "reproj") adjuster = new detail::BundleAdjusterReproj();
    else if (opts.ba_cost_func == "ray") adjuster = new detail::BundleAdjusterRay();
    else
    {
        cout << "Unknown bundle adjustment cost function: '" << opts.ba_cost_func << "'.\n";
        return -1;
    }
    adjuster->setConfThresh(opts.conf_thresh);
    Mat_<uchar> refine_mask = Mat::zeros(3, 3, CV_8U);
    if (opts.ba_refine_mask[0] == 'x') refine_mask(0,0) = 1;
    if (opts.ba_refine_mask[1] == 'x') refine_mask(0,1) = 1;
    if (opts.ba_refine_mask[2] == 'x') refine_mask(0,2) = 1;
    if (opts.ba_refine_mask[3] == 'x') refine_mask(1,1) = 1;
    if (opts.ba_refine_mask[4] == 'x') refine_mask(1,2) = 1;
    adjuster->setRefinementMask(refine_mask);
    (*adjuster)(features, pairwise_matches, cameras);
    if (dc.procid() == 0)
        logstream(LOG_EMPH) << "Bundle Adjustment, time: " << ((getTickCount() - t) / getTickFrequency()) << " sec\n";

    ///////////////////////////////////////////////////////
    // Find median focal length    
    vector<double> focals;
    for (size_t i = 0; i < cameras.size(); ++i)
    {
        if (dc.procid() == 0)
            logstream(LOG_EMPH) << "Camera #" << i << ":\n" << cameras[i].K() << "\n\n";
        focals.push_back(cameras[i].focal);
    }
   
    sort(focals.begin(), focals.end());

    if (opts.verbose > 2)
    	LOGLN("Focals size: " << focals.size() << "	focals: " << focals << "\n");    

    if (focals.size() % 2 == 1)
        opts.warped_image_scale = static_cast<float>(focals[focals.size() / 2]);
    else
        opts.warped_image_scale = static_cast<float>(focals[focals.size() / 2 - 1] + focals[focals.size() / 2]) * 0.5f;
   
    ///////////////////////////////////////////////////////
    // Wave-Correction
    vector<Mat> rmats;
    WaveCorrectKind wave_correct; // in future it needs to be added in clopts to change it in runtime
    if (opts.wave_correct_type == "horiz")
	wave_correct = detail::WAVE_CORRECT_HORIZ;
    if (opts.wave_correct_type == "vert")
	wave_correct = detail::WAVE_CORRECT_VERT;
    for (size_t i = 0; i < cameras.size(); ++i)
        rmats.push_back(cameras[i].R);
    waveCorrect(rmats, wave_correct);
    for (size_t i = 0; i < cameras.size(); ++i)
        cameras[i].R = rmats[i];
   
    //} // End of if procid=0

   
    ///////////////////////////////////////////////////////
    // Create a second graph with cameras
    graph_type graph_cam(dc, clopts);

    // load the graph
    if (dc.procid()==0) {
    vertex_loader(graph_cam, img_dir, cameras);
    graph_cam.load(graph_path, edge_loader);
    }
    graph_cam.finalize();

    ///////////////////////////////////////////////////////
    // Warp Images in parallel on vertices
    graph_cam.transform_vertices(warp_images);

    ///////////////////////////////////////////////////////
    // Gain Normalize

    ///////////////////////////////////////////////////////
    // Find seams in parallel on edges
    graph_cam.transform_edges(find_seams);

    ///////////////////////////////////////////////////////
    // Composite Images in parallel on vertices
    graph_cam.transform_vertices(composite_images);
   

    ///////////////////////////////////////////////////////
    // Second Graphlab Engine
    engine_type engine_cam(dc, graph_cam, opts.exec_type, clopts);

    ///////////////////////////////////////////////////////
    // blend images, gather vertices
    VecVD veclist = engine_cam.map_reduce_vertices<VecVD>(compile_vertices);
    vector<Point> corner(veclist.size());
    vector<Mat> img_warped(veclist.size());
    vector<Mat> mask_warped(veclist.size());
    vector<Size> size(veclist.size());
    Mat img_warped_s;

    for (size_t i=0; i!=veclist.size(); ++i)
    {
        corner[i] = veclist[i].corner;

	if (opts.verbose > 2)
            LOGLN("Blending corners x : " << veclist[i].corner.x << "   y : " << veclist[i].corner.y << "\n");

        img_warped[i] = veclist[i].img_warped;
        mask_warped[i] = veclist[i].mask_warped;
        size[i] = veclist[i].warp_size;

        if (opts.verbose > 2)
            LOGLN("Blending sizes height : " << size[i].height << "   width : " << size[i].width << "\n");
    }
    veclist.clear();
   
    num_images = corner.size();
    
    Ptr<Blender> blender;
    
    int blend_type;
    if (opts.blending_type == "no")
	blend_type = Blender::NO;
    if (opts.blending_type == "feather")
	blend_type = Blender::FEATHER;
    if (opts.blending_type == "multiband")
        blend_type = Blender::MULTI_BAND;

    bool try_gpu = false;
            
    if (blender.empty())
    {
        blender = Blender::createDefault(blend_type, try_gpu);
        Size dst_sz = resultRoi(corner, size).size();
        float blend_width = sqrt(static_cast<float>(dst_sz.area())) * opts.blend_strength / 100.f;
        if (blend_width < 1.f)
            blender = Blender::createDefault(Blender::NO, try_gpu);
        else if (blend_type == Blender::MULTI_BAND)
        {
            MultiBandBlender* mb = dynamic_cast<MultiBandBlender*>(static_cast<Blender*>(blender));
            mb->setNumBands(static_cast<int>(ceil(log(blend_width)/log(2.)) - 1.));
            LOGLN("Multi-band blender, number of bands: " << mb->numBands());
        }
        else if (blend_type == Blender::FEATHER)
        {
            FeatherBlender* fb = dynamic_cast<FeatherBlender*>(static_cast<Blender*>(blender));
            fb->setSharpness(1.f/blend_width);
            LOGLN("Feather blender, sharpness: " << fb->sharpness());
        }
        blender->prepare(corner, size);
    }

    // Blend all images 
    for (int j=0; j!=num_images; ++j)
    {
	img_warped[j].convertTo(img_warped_s, CV_16S);
        blender->feed(img_warped_s, mask_warped[j], corner[j]);
        img_warped_s.release();
    }
    
    Mat result, result_mask;
    blender->blend(result, result_mask);

    imwrite(opts.result_name, result);
    
    LOGLN("Finished, total time: " << ((getTickCount() - app_start_time) / getTickFrequency()) << " sec");

    ///////////////////////////////////////////////////////
    // Run everything
//    engine.signal_all();
//    graphlab::timer timer;
//    engine.start();  
//    const double runtime = timer.current_time();
//    dc.cout()
//    << "----------------------------------------------------------" << std::endl
//    << "Final Runtime (seconds):   " << runtime
//    << std::endl
//    << "Updates executed: " << engine.num_updates() << std::endl
//    << "Update Rate (updates/second): "
//    << engine.num_updates() / runtime << std::endl;
       
}


