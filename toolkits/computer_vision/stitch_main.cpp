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
 * multiple images into a panorama. 
 *
 *  \author Dhruv Batra
 */


#include "stitch_grlab.hpp"

void printUsage()
{
    cout <<
    "Rotation model images stitcher.\n\n"
    "stitching_detailed img1 img2 [...imgN] [flags]\n\n"
    "Flags:\n"
    "  --preview\n"
    "      Run stitching in the preview mode. Works faster than usual mode,\n"
    "      but output image will have lower resolution.\n"
    "  --try_gpu (yes|no)\n"
    "      Try to use GPU. The default value is 'no'. All default values\n"
    "      are for CPU mode.\n"
    "\nMotion Estimation Flags:\n"
    "  --work_megapix <float>\n"
    "      Resolution for image registration step. The default is 0.6 Mpx.\n"
    "  --features (surf|orb)\n"
    "      Type of features used for images matching. The default is surf.\n"
    "  --match_conf <float>\n"
    "      Confidence for feature matching step. The default is 0.65 for surf and 0.3 for orb.\n"
    "  --conf_thresh <float>\n"
    "      Threshold for two images are from the same panorama confidence.\n"
    "      The default is 1.0.\n"
    "  --ba (reproj|ray)\n"
    "      Bundle adjustment cost function. The default is ray.\n"
    "  --ba_refine_mask (mask)\n"
    "      Set refinement mask for bundle adjustment. It looks like 'x_xxx',\n"
    "      where 'x' means refine respective parameter and '_' means don't\n"
    "      refine one, and has the following format:\n"
    "      <fx><skew><ppx><aspect><ppy>. The default mask is 'xxxxx'. If bundle\n"
    "      adjustment doesn't support estimation of selected parameter then\n"
    "      the respective flag is ignored.\n"
    "  --wave_correct (no|horiz|vert)\n"
    "      Perform wave effect correction. The default is 'horiz'.\n"
    "  --save_graph <file_name>\n"
    "      Save matches graph represented in DOT language to <file_name> file.\n"
    "      Labels description: Nm is number of matches, Ni is number of inliers,\n"
    "      C is confidence.\n"
    "\nCompositing Flags:\n"
    "  --warp (plane|cylindrical|spherical|fisheye|stereographic|compressedPlaneA2B1|compressedPlaneA1.5B1|compressedPlanePortraitA2B1|compressedPlanePortraitA1.5B1|paniniA2B1|paniniA1.5B1|paniniPortraitA2B1|paniniPortraitA1.5B1|mercator|transverseMercator)\n"
    "      Warp surface type. The default is 'spherical'.\n"
    "  --seam_megapix <float>\n"
    "      Resolution for seam estimation step. The default is 0.1 Mpx.\n"
    "  --seam (no|voronoi|gc_color|gc_colorgrad)\n"
    "      Seam estimation method. The default is 'gc_color'.\n"
    "  --compose_megapix <float>\n"
    "      Resolution for compositing step. Use -1 for original resolution.\n"
    "      The default is -1.\n"
    "  --expos_comp (no|gain|gain_blocks)\n"
    "      Exposure compensation method. The default is 'gain_blocks'.\n"
    "  --blend (no|feather|multiband)\n"
    "      Blending method. The default is 'multiband'.\n"
    "  --blend_strength <float>\n"
    "      Blending strength from [0,100] range. The default is 5.\n"
    "  --output <result_img>\n"
    "      The default is 'result.jpg'.\n";
}


// Default command line args
vector<string> img_names;
bool preview = false;
bool try_gpu = false;
double work_megapix = 0.6;
double seam_megapix = 0.1;
double compose_megapix = -1;
float conf_thresh = 1.f;
string features = "surf";
string ba_cost_func = "ray";
string ba_refine_mask = "xxxxx";
bool do_wave_correct = true;
WaveCorrectKind wave_correct = detail::WAVE_CORRECT_HORIZ;
bool save_graph = false;
std::string save_graph_to;
string warp_type = "spherical";
int expos_comp_type = ExposureCompensator::GAIN_BLOCKS;
float match_conf = 0.3f;
string seam_find_type = "gc_color";
int blend_type = Blender::MULTI_BAND;
float blend_strength = 5;
string result_name = "result.jpg";

int parseCmdArgs(int argc, char** argv)
{
    if (argc == 1)
    {
        printUsage();
        return -1;
    }
    for (int i = 1; i < argc; ++i)
    {
        if (string(argv[i]) == "--help" || string(argv[i]) == "/?")
        {
            printUsage();
            return -1;
        }
        else if (string(argv[i]) == "--preview")
        {
            preview = true;
        }
        else if (string(argv[i]) == "--try_gpu")
        {
            if (string(argv[i + 1]) == "no")
                try_gpu = false;
            else if (string(argv[i + 1]) == "yes")
                try_gpu = true;
            else
            {
                cout << "Bad --try_gpu flag value\n";
                return -1;
            }
            i++;
        }
        else if (string(argv[i]) == "--work_megapix")
        {
            work_megapix = atof(argv[i + 1]);
            i++;
        }
        else if (string(argv[i]) == "--seam_megapix")
        {
            seam_megapix = atof(argv[i + 1]);
            i++;
        }
        else if (string(argv[i]) == "--compose_megapix")
        {
            compose_megapix = atof(argv[i + 1]);
            i++;
        }
        else if (string(argv[i]) == "--result")
        {
            result_name = argv[i + 1];
            i++;
        }
        else if (string(argv[i]) == "--features")
        {
            features = argv[i + 1];
            if (features == "orb")
                match_conf = 0.3f;
            i++;
        }
        else if (string(argv[i]) == "--match_conf")
        {
            match_conf = static_cast<float>(atof(argv[i + 1]));
            i++;
        }
        else if (string(argv[i]) == "--conf_thresh")
        {
            conf_thresh = static_cast<float>(atof(argv[i + 1]));
            i++;
        }
        else if (string(argv[i]) == "--ba")
        {
            ba_cost_func = argv[i + 1];
            i++;
        }
        else if (string(argv[i]) == "--ba_refine_mask")
        {
            ba_refine_mask = argv[i + 1];
            if (ba_refine_mask.size() != 5)
            {
                cout << "Incorrect refinement mask length.\n";
                return -1;
            }
            i++;
        }
        else if (string(argv[i]) == "--wave_correct")
        {
            if (string(argv[i + 1]) == "no")
                do_wave_correct = false;
            else if (string(argv[i + 1]) == "horiz")
            {
                do_wave_correct = true;
                wave_correct = detail::WAVE_CORRECT_HORIZ;
            }
            else if (string(argv[i + 1]) == "vert")
            {
                do_wave_correct = true;
                wave_correct = detail::WAVE_CORRECT_VERT;
            }
            else
            {
                cout << "Bad --wave_correct flag value\n";
                return -1;
            }
            i++;
        }
        else if (string(argv[i]) == "--save_graph")
        {
            save_graph = true;
            save_graph_to = argv[i + 1];
            i++;
        }
        else if (string(argv[i]) == "--warp")
        {
            warp_type = string(argv[i + 1]);
            i++;
        }
        else if (string(argv[i]) == "--expos_comp")
        {
            if (string(argv[i + 1]) == "no")
                expos_comp_type = ExposureCompensator::NO;
            else if (string(argv[i + 1]) == "gain")
                expos_comp_type = ExposureCompensator::GAIN;
            else if (string(argv[i + 1]) == "gain_blocks")
                expos_comp_type = ExposureCompensator::GAIN_BLOCKS;
            else
            {
                cout << "Bad exposure compensation method\n";
                return -1;
            }
            i++;
        }
        else if (string(argv[i]) == "--seam")
        {
            if (string(argv[i + 1]) == "no" ||
                string(argv[i + 1]) == "voronoi" ||
                string(argv[i + 1]) == "gc_color" ||
                string(argv[i + 1]) == "gc_colorgrad")
                seam_find_type = argv[i + 1];
            else
            {
                cout << "Bad seam finding method\n";
                return -1;
            }
            i++;
        }
        else if (string(argv[i]) == "--blend")
        {
            if (string(argv[i + 1]) == "no")
                blend_type = Blender::NO;
            else if (string(argv[i + 1]) == "feather")
                blend_type = Blender::FEATHER;
            else if (string(argv[i + 1]) == "multiband")
                blend_type = Blender::MULTI_BAND;
            else
            {
                cout << "Bad blending method\n";
                return -1;
            }
            i++;
        }
        else if (string(argv[i]) == "--blend_strength")
        {
            blend_strength = static_cast<float>(atof(argv[i + 1]));
            i++;
        }
        else if (string(argv[i]) == "--output")
        {
            result_name = argv[i + 1];
            i++;
        }
        else
            img_names.push_back(argv[i]);
    }
    if (preview)
    {
        compose_megapix = 0.6;
    }
    return 0;
}


/////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv) 
{
    
    string img_dir; 
    string graph_dir;
    string output_dir = "stitch_output";

    ///////////////////////////////////////////////////////
    // Set up Graphlab
    global_logger().set_log_level(LOG_INFO);
    global_logger().set_log_to_console(true);
    ///! Initialize control plain using mpi
    graphlab::mpi_tools::init(argc, argv);
    graphlab::distributed_control dc;
    const std::string description = "Stitching Application";

    ///////////////////////////////////////////////////////
    // Set up OpenCV
    cv::setBreakOnError(true);

    ///////////////////////////////////////////////////////
    // Graphlab parse input
    graphlab::command_line_options clopts(description);
    clopts.attach_option("img", img_dir,
                         "The directory containing the images");
    clopts.add_positional("img");
    clopts.attach_option("graph", graph_dir,
                         "The directory containing the adjacency graph");
    clopts.add_positional("graph");
    clopts.attach_option("output", output_dir,
                         "The directory in which to save the output");

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
    
//    if(graph_dir.empty()) 
//    {
//        logstream(LOG_ERROR) << "No graph was provided." << std::endl;
//        return EXIT_FAILURE;
//    }
    
    // display settings  
    if(dc.procid() == 0) 
    {
        std::cout 
        << "ncpus:          " << clopts.get_ncpus() << std::endl
        << "img_dir:        " << img_dir << std::endl
        << "graph_dir:      " << graph_dir << std::endl;
    }
    

    
    ///////////////////////////////////////////////////////
    // Graphlab Graph
    graph_type graph(dc, clopts);  
        
    // load the graph
    //graph.load(img_dir, vertex_loader);
    //graph.load(graph_dir, edge_loader);
    vertex_loader(dc, graph, img_dir);
    graph.finalize();
    
    ///////////////////////////////////////////////////////
    // Computer features in parallel
    graph.transform_vertices(compute_features);

    
    ///////////////////////////////////////////////////////
    // Graphlab Engine
    engine_type engine(dc, graph, clopts);

    
    ///////////////////////////////////////////////////////
    // Run everything
    engine.signal_all();
    graphlab::timer timer;
    engine.start();  
    const double runtime = timer.current_time();
    dc.cout() 
    << "----------------------------------------------------------" << std::endl
    << "Final Runtime (seconds):   " << runtime 
    << std::endl
    << "Updates executed: " << engine.num_updates() << std::endl
    << "Update Rate (updates/second): " 
    << engine.num_updates() / runtime << std::endl;
    
//    int retval = parseCmdArgs(argc, argv);
//    if (retval)
//        return retval;
//    
//    // Check if have enough images
//    int num_images = static_cast<int>(img_names.size());
//    if (num_images < 2)
//    {
//        LOGLN("Need more images");
//        return -1;
//    }
//    
//    double work_scale = 1, seam_scale = 1, compose_scale = 1;
//    bool is_work_scale_set = false, is_seam_scale_set = false, is_compose_scale_set = false;
    
//    LOGLN("Finding features...");
    
}

