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
 * This file contains an example of graphlab used for MAP inference 
 * in a discrete graphical model (pairwise MRF). The algorithm
 * implemented is the MPLP LP-Relaxation scheme of Globerson & Jaakkola. 
 *
 *  \author Dhruv Batra
 */


#include <vector>
#include <string>
#include <fstream>


#include <Eigen/Dense>


#include <Magick++.h> 
#undef restrict

#include <graphlab.hpp>

#include "eigen_serialization.hpp"

#include <graphlab/macros_def.hpp>

typedef Eigen::VectorXd vector;
typedef Eigen::MatrixXd matrix;

template <typename T>
inline std::ostream& operator<<(std::ostream& os, std::vector<T>& x)
{
    typename std::vector<T>::const_iterator i(x.begin());
    while(i != x.end()) os << *i++ << ' ';
    return os;
}

// Global variables
size_t NCOLORS;
double SIGMA;
double BOUND;

// LP-based upper-bound on MAP
double MAPLPval;

// Shared base edge potential
matrix THETA_ij; 

// STRUCTS (Edge and Vertex data) =============================================>

/**
 * Each GraphLab vertex is a (pairwise) factor from the MRF
 */
struct vertex_data {
    /** variable ids */
    int i, j; 
    
    // degree of these nodes in the MRF
    int deg_i, deg_j; 
    
    /** observed color for each variable */
    float obs_color_i, obs_color_j;
    /** predicted color for each variable */
    float pred_color_i, pred_color_j;
    /** dual variables being optimized (or messages) */
    vector delf_i, delf_j;  
    // constructor
    vertex_data(): i(-1), j(-1), deg_i(0), deg_j(0), obs_color_i(-1), obs_color_j(-1) { }
    void save(graphlab::oarchive& arc) const 
    {
        arc << i << j 
        << deg_i << deg_j 
        << obs_color_i << obs_color_j 
        << pred_color_i << pred_color_j 
        << delf_i << delf_j;
    }
    void load(graphlab::iarchive& arc) {
        arc >> i >> j 
        >> deg_i >> deg_j
        >> obs_color_i >> obs_color_j 
        >> pred_color_i >> pred_color_j 
        >> delf_i >> delf_j;
    }
}; // End of vertex data


// /**
//  * The data associated with a pair of factors in a pairwise MRF
//  */
// struct edge_data : public graphlab::IS_POD_TYPE {
//   // primal labelling; We assume pairwise factors, so intersection has
//   // a single node
//   int pred_label; 
// }; // End of edge data
typedef graphlab::empty edge_data;

/**
 * The graph type
 */
typedef graphlab::distributed_graph<vertex_data, edge_data> graph_type;


// GraphLab Vertex Program ====================================================
/**
 * The type passed around during the gather phase
 */
struct gather_type 
{
    vector delf_i, delf_j;
    
    gather_type& operator+=(const gather_type& other) 
    {
        if(!other.delf_i.size() == 0) 
        {
            if(delf_i.size() == 0) delf_i = other.delf_i;
            else delf_i += other.delf_i;
        }
        if(!other.delf_j.size() == 0) 
        {
            if(delf_j.size() == 0) delf_j = other.delf_j;
            else delf_j += other.delf_j;
        }
        return *this;
    } // end of operator +=
    void save(graphlab::oarchive& arc) const 
    {
        arc << delf_i << delf_j;
    }
    void load(graphlab::iarchive& arc) 
    {
        arc >> delf_i >> delf_j;
    }
}; // end of gather type



/** 
 * The core belief propagation update function.  This update satisfies
 * the graphlab update_function interface.  
 */
class mplp_vertex_program : 
public graphlab::ivertex_program<graph_type, gather_type, 
graphlab::messages::sum_priority>,
public graphlab::IS_POD_TYPE {
private:
    double priority;
public:
    
    mplp_vertex_program() : priority(0) { }
    
    // void save(graphlab::oarchive& arc) const { /** save members */ }
    // void load(graphlab::iarchive& arc) { /** load members */ }
    
    /**
     * The init function is called once for each vertex before the
     * start of the GraphLab program.  If the vertex program does not
     * implement this function then the default implementation (NOP)
     * is used.
     */
    void init(icontext_type& context, vertex_type& vertex)
    { 
        vertex_data& vdata = vertex.data();
        
        // Create zero messages
        vdata.delf_i = vector::Zero(NCOLORS);
        vdata.delf_j = vector::Zero(NCOLORS);
        
        // Initialize predicted values and add to global bound
        vdata.pred_color_i = vdata.pred_color_j = 0; 
        vector theta_i = make_unary_potential(vertex, 'i');
        vector theta_j = make_unary_potential(vertex, 'j');

        MAPLPval += theta_i[vdata.pred_color_i]/vdata.deg_i; 
        MAPLPval += theta_j[vdata.pred_color_j]/vdata.deg_j;
        
        MAPLPval += THETA_ij(vdata.pred_color_i, vdata.pred_color_j);
    }
    
    /**
     * Recv message is called by the engine to receive a message to this
     * vertex program.  The vertex program can use this to initialize
     * any state before entering the gather phase.  If the vertex
     * program does not implement this function then the default
     * implementation (NOP) is used.
     */
    // void recv_message(icontext_type& context, const vertex_type& vertex, 
    //                   const message_type& msg) { /** NOP */ }
    
    /**
     * Since the MRF is undirected we will use all edges for gather and
     * scatter
     */
    edge_dir_type gather_edges(icontext_type& context,
                               const vertex_type& vertex) const { 
        return graphlab::ALL_EDGES; 
    }; // end of gather_edges 
    
    
    // Run the gather operation over all in edges
    gather_type gather(icontext_type& context, const vertex_type& target_vertex, 
                       edge_type& edge) const 
    {
        const vertex_type source_vertex = get_other_vertex(edge, target_vertex);   
        const vertex_data& source_vdata = source_vertex.data();
        const vertex_data& target_vdata = target_vertex.data();
        
        // Accumulate message
        gather_type ret_value;
        if (target_vdata.i == source_vdata.i)
            ret_value.delf_i = source_vdata.delf_i;
        else if (target_vdata.j == source_vdata.i)
            ret_value.delf_j = source_vdata.delf_i;
        else if (target_vdata.i == source_vdata.j)
            ret_value.delf_i = source_vdata.delf_j;
        else if (target_vdata.j == source_vdata.j)
            ret_value.delf_j = source_vdata.delf_j;
        else assert(false); // invalid state
        
        return ret_value;
    } // end of gather
    
    /** Update the dual parameters */
    void apply(icontext_type& context, vertex_type& vertex, 
               const gather_type& sum) {
        // Make sure this vertex has neighbours. Everyone should have neighbours
        ASSERT_GT(vertex.num_in_edges() + vertex.num_out_edges(), 0);
        vertex_data& vdata = vertex.data();  
        vector theta_i = make_unary_potential(vertex, 'i');
        vector theta_j = make_unary_potential(vertex, 'j');
        
        
        ASSERT_EQ(THETA_ij.rows(), theta_i.size());
        ASSERT_EQ(THETA_ij.rows(), sum.delf_i.size());
        ASSERT_EQ(THETA_ij.cols(), theta_j.size());
        ASSERT_EQ(THETA_ij.cols(), sum.delf_j.size());   
        
        std::cout << MAPLPval << "\t" ;
        // Remove contribution of old labels from MAPLPval
        MAPLPval -= theta_i[vdata.pred_color_i]/vdata.deg_i; 
        MAPLPval -= theta_j[vdata.pred_color_j]/vdata.deg_j;        
        MAPLPval -= THETA_ij(vdata.pred_color_i, vdata.pred_color_j);        
        
        // Compute the prediction
        Eigen::MatrixXf::Index maxI = 0, maxJ = 0;
        (theta_i + sum.delf_i + vdata.delf_i).maxCoeff(&maxI);
        (theta_j + sum.delf_j + vdata.delf_j).maxCoeff(&maxJ);
        vdata.pred_color_i = maxI;
        vdata.pred_color_j = maxJ;
        
        // Add new contributions to MAPLPval
        MAPLPval += theta_i[vdata.pred_color_i]/vdata.deg_i; 
        MAPLPval += theta_j[vdata.pred_color_j]/vdata.deg_j;        
        MAPLPval += THETA_ij(vdata.pred_color_i, vdata.pred_color_j);        

        
        // Backup the old prediction
        const vector old_delf_i = vdata.delf_i;
        const vector old_delf_j = vdata.delf_j;
        
        // Update del fi
        vdata.delf_i = -(theta_i + sum.delf_i)/2 + 
        (THETA_ij + sum.delf_j.rowwise().replicate(THETA_ij.rows())).
        rowwise().maxCoeff()/2;
        // Update del fj
        vdata.delf_j = -(theta_j + sum.delf_j)/2 + 
        ((THETA_ij + sum.delf_i.rowwise().replicate(THETA_ij.cols()).transpose()).
         colwise().maxCoeff()).transpose()/2;
        
        const double residual = (vdata.delf_i - old_delf_i).cwiseAbs().sum() +
        (vdata.delf_j - old_delf_j).cwiseAbs().sum();
        
        priority = residual;
        //std::cout << "priority: " << priority << std::endl;
        //std::cout << MAPLPval << std::endl;
        //test code; for now, only run 1 iteration
        priority = 0;
    } // end of apply
    
    /**
     * Since the MRF is undirected we will use all edges for gather and
     * scatter
     */
    edge_dir_type scatter_edges(icontext_type& context,
                                const vertex_type& vertex) const { 
        return priority < BOUND? graphlab::NO_EDGES : graphlab::ALL_EDGES; 
    }; // end of gather_edges 
    
    
    /** reschedule neighbors with a given priority and updated
     predictions on each edge*/
    void scatter(icontext_type& context, const vertex_type& vertex, 
                 edge_type& edge) const {  
        context.signal(get_other_vertex(edge, vertex), priority);
    } // end of scatter
    
private:
    
    /**
     * Construct the unary evidence potential
     */
    vector make_unary_potential(const vertex_type& vertex, 
                                const char varid) const {    
        vector potential(NCOLORS);
        const double obs = varid == 'i'? 
        vertex.data().obs_color_i : vertex.data().obs_color_j;
        const double sigmaSq = SIGMA*SIGMA;
        for(int pred = 0; pred < potential.size(); ++pred) {
            potential(pred) = -(obs - pred)*(obs - pred) / (2.0 * sigmaSq);
        }
        potential /= potential.sum();
        return potential;
    } // end of make_potentail
    
    /**
     * Return the other vertex
     */
    vertex_type get_other_vertex(edge_type& edge, 
                                 const vertex_type& vertex) const {
        return vertex.id() == edge.source().id()? edge.target() : edge.source();
    } // end of other_vertex
    
}; // end of MPLP vertex program




/**
 * Define the engine type
 */
//typedef graphlab::synchronous_engine<mplp_vertex_program> engine_type;
typedef graphlab::async_consistent_engine<mplp_vertex_program> engine_type;
//typedef graphlab::asynchronous_consistent_engine<mplp_vertex_program> engine_type;











// Helper functions ===========================================================>
graphlab::vertex_id_type pixel_ind(size_t rows, size_t cols,
                                   size_t r, size_t c) {
    return r * cols + c;
}; // end of pixel_ind


graphlab::vertex_id_type factor_ind(size_t rows, size_t cols,
                                    size_t i, size_t j) {
    if(i > j) std::swap(i,j);
    return i * (rows * cols) + j;
}; // end of factor_ind


void create_synthetic_cluster_graph(graphlab::distributed_control& dc,
                                    graph_type& graph,
                                    const size_t rows, const size_t cols) {
    dc.barrier();
    // Generate the image on all machines --------------------------------------->
    // Need to ensure that all machines generate the same noisy image
    graphlab::random::generator gen; gen.seed(314);
    std::vector<float>    obs_pixels(rows * cols);
    std::vector<uint16_t> true_pixels(rows * cols);
    const double center_r = rows / 2.0;
    const double center_c = cols / 2.0;
    const double max_radius = std::min(rows, cols) / 2.0;
    for(size_t r = 0; r < rows; ++r) 
    {
        for(size_t c = 0; c < cols; ++c) 
        {
            // Compute the true pixel value
            const double distance = sqrt((r-center_r)*(r-center_r) + 
                                         (c-center_c)*(c-center_c));
            // Compute ring of sunset
            const uint16_t ring_color =  
            std::floor(std::min(1.0, distance/max_radius) * (NCOLORS - 1) );
            // Compute the true pixel color by masking with the horizon
            const uint16_t true_color = r < rows/2 ? ring_color : 0;
            // compute the predicted color
            const float obs_color = true_color + gen.normal(0, SIGMA);
            // determine the true pixel id
            const size_t pixel = pixel_ind(rows,cols,r,c);
            true_pixels[pixel] = true_color; obs_pixels[pixel] = obs_color;
        } // end of loop over cols
    } // end of loop over rows
    
    if(dc.procid() == 0) 
    {
        // temp code 
        std::ofstream ne, ee; 
        ne.open("./node_en.txt");
        //ee.open("./edge_en.txt");
        
        ne << rows*cols << " " << NCOLORS << " "
        << 0 << " " << rows << " " << cols << std::endl;
        // end temp
        
        std::vector<graphlab::vertex_id_type> nbrs;
        // load the graph
        for(size_t r = 0; r < rows; ++r) 
        {
            for(size_t c = 0; c < cols; ++c) 
            {
                // temp code to write out potential to file:
                std::vector<double> potential(NCOLORS);
                const double obs = obs_pixels[pixel_ind(rows,cols,r,c)];
                const double sigmaSq = SIGMA*SIGMA;
                double sum = 0;
                for(int pred = 0; pred < potential.size(); ++pred) 
                {
                    potential[pred] = +(obs - pred)*(obs - pred) / (2.0 * sigmaSq);
                    sum += potential[pred];
                }
                for(int pred = 0; pred < potential.size(); ++pred)                 
                    potential[pred] /= sum;
                //ne << true_pixels[pixel_ind(rows,cols,r,c)] << " " <<  obs << " " << potential << std::endl;
                ne << 4 - int((r==0)||(r==(rows-1))) - int((c==0)||(c==(cols-1))) << " " << potential << std::endl;
                // end temp
                
                // Add the two vertices (factors to the right and below this
                // pixel)
                if(r + 1 < rows) 
                {
                    vertex_data vdata;
                    vdata.i = pixel_ind(rows,cols,r,c);
                    vdata.j = pixel_ind(rows,cols,r+1,c);
                    vdata.deg_i = 4 - int((r==0)||(r==(rows-1))) - int((c==0)||(c==(cols-1)));
                    vdata.deg_j = 4 - int(((r+1)==0)||((r+1)==(rows-1))) - int((c==0)||(c==(cols-1)));
                    vdata.obs_color_i = obs_pixels[vdata.i];
                    vdata.obs_color_j = obs_pixels[vdata.j];
                    graph.add_vertex(factor_ind(rows,cols,vdata.i,vdata.j), vdata);
                }
                if(c + 1 < cols) 
                {
                    vertex_data vdata;
                    vdata.i = pixel_ind(rows,cols,r,c);
                    vdata.j = pixel_ind(rows,cols,r,c+1);
                    vdata.deg_i = 4 - int((r==0)||(r==(rows-1))) - int((c==0)||(c==(cols-1)));
                    vdata.deg_j = 4 - int((r==0)||(r==(rows-1))) - int(((c+1)==0)||((c+1)==(cols-1)));
                    vdata.obs_color_i = obs_pixels[vdata.i];
                    vdata.obs_color_j = obs_pixels[vdata.j];
                    graph.add_vertex(factor_ind(rows,cols,vdata.i,vdata.j), vdata);
                }
                // Compute all the factors that contain this pixel
                nbrs.clear();
                if(r+1 < rows)
                    nbrs.push_back(factor_ind(rows,cols,
                                              pixel_ind(rows,cols,r,c),
                                              pixel_ind(rows,cols,r+1,c)));
                if(r-1 < rows)
                    //if(r-1 >= 0)
                    nbrs.push_back(factor_ind(rows,cols,
                                              pixel_ind(rows,cols,r-1,c),
                                              pixel_ind(rows,cols,r,c)));
                if(c+1 < cols)
                    nbrs.push_back(factor_ind(rows,cols,
                                              pixel_ind(rows,cols,r,c),
                                              pixel_ind(rows,cols,r,c+1)));
                if(c-1 < cols)
                    //if(c-1 >= 0)
                    nbrs.push_back(factor_ind(rows,cols,
                                              pixel_ind(rows,cols,r,c-1),
                                              pixel_ind(rows,cols,r,c)));
                // construct the clique over the factors
                for(size_t i = 0; i < nbrs.size(); ++i) 
                {
                    for(size_t j = i+1; j < nbrs.size(); ++j) 
                    {
                        graph.add_edge(nbrs[i], nbrs[j]);
                    }
                }
            } // end of for cols
        } // end of for rows
        
        // temp code
        ne.close(); //ee.close();
    } // end of if proc 0
    dc.barrier();
} // end of create synthetic cluster graph





void initialize_theta_ij(const std::string& smoothing,
                         const double lambda) 
{
    THETA_ij.resize(NCOLORS, NCOLORS);
    // Set the smoothing type
    if(smoothing == "laplace") 
    {
        for(int i = 0; i < THETA_ij.rows(); ++i) 
            for(int j = 0; j < THETA_ij.cols(); ++j) 
                THETA_ij(i,j) = -std::abs(double(i) - double(j)) * lambda;
    } 
    else 
    {   
        for(int i = 0; i < THETA_ij.rows(); ++i) 
            for(int j = 0; j < THETA_ij.cols(); ++j) 
                THETA_ij(i,j) = -(i == j? 0 : lambda);
    } 
} // end of initialize_theta_ij



template<typename T>
struct merge_reduce {
    std::set<T> values;
    void save(graphlab::oarchive& arc) const { arc << values; }
    void load(graphlab::iarchive& arc) { arc >> values; }
    merge_reduce& operator+=(const merge_reduce& other) {
        values.insert(other.values.begin(), other.values.end());
        return *this;
    }
}; // end of merge_reduce

typedef std::pair<graphlab::vertex_id_type, float> pred_pair_type; 
typedef merge_reduce<pred_pair_type> merge_reduce_type;

merge_reduce_type pred_map_function(graph_type::vertex_type vertex) {
    merge_reduce<pred_pair_type> ret;
    ret.values.insert(pred_pair_type(vertex.data().i, vertex.data().pred_color_i));
    ret.values.insert(pred_pair_type(vertex.data().j, vertex.data().pred_color_j));
    return ret;
} // end of pred_map_function

merge_reduce_type obs_map_function(graph_type::vertex_type vertex) {
    merge_reduce<pred_pair_type> ret;
    ret.values.insert(pred_pair_type(vertex.data().i, vertex.data().obs_color_i));
    ret.values.insert(pred_pair_type(vertex.data().j, vertex.data().obs_color_j));
    return ret;
} // end of obs_map_function




std::pair<int,int> ind2sub(size_t rows, size_t cols,
                           size_t ind) {
    return std::make_pair(ind / cols, ind % cols);
}; // end of sub2ind


// /**
//  * Saving an image as a pgm file.
//  */
// void save_image(const size_t rows, const size_t cols,
//                 const std::set<pred_pair_type>& values,
//                 const std::string& fname) {
//     std::cout << "NPixels: " << values.size() << std::endl;
//     image img(rows, cols);
//     foreach(pred_pair_type pair, values) 
//     img.pixel(pair.first) = pair.second;
//     img.save(fname);
// } // end of save_image


/**
 * Saving an image as a pgm file.
 */
void save_image(const size_t rows, const size_t cols,
                const std::set<pred_pair_type>& values,
                const std::string& fname) {
    using namespace Magick;
    std::cout << "NPixels: " << values.size() << std::endl;
    // determine the max and min colors
    float max_color = -std::numeric_limits<float>::max();
    float min_color =  std::numeric_limits<float>::max();
    foreach(pred_pair_type pair, values) {
        max_color = std::max(max_color, pair.second);
        min_color = std::min(min_color, pair.second);
    }
    
    Image img(Magick::Geometry(rows, cols), "white");
    // img.modifyImage();
    // Pixels img_cache(img);
    // PixelPackets* pixels = img_cache.
    foreach(pred_pair_type pair, values) {
        std::pair<int,int> coords = ind2sub(rows,cols, pair.first);
        float value = (pair.second - min_color) / (max_color - min_color);
        Color color(MaxRGB * value, MaxRGB * value, MaxRGB * value, 0);
        img.pixelColor(coords.second, coords.first, color);
    }
    img.write(fname);  
} // end of save_image



// MAIN =======================================================================>
int main(int argc, char** argv) {
    std::cout << "This program creates and denoises a synthetic " << std::endl
    << "image using loopy belief propagation inside " << std::endl
    << "the graphlab framework." << std::endl;
    
    // // set the global logger
    // global_logger().set_log_level(LOG_WARNING);
    // global_logger().set_log_to_console(true);
    
    // Set initial values for members ------------------------------------------->
    NCOLORS = 5;
    SIGMA = 2;
    BOUND = 1E-4;
    
    
    size_t nrows = 200;
    size_t ncols = 200;
    double lambda = 0.2;
    
    std::string smoothing = "square";
    
    std::string orig_fn =  "source_img.jpeg";
    std::string noisy_fn = "noisy_img.jpeg";
    std::string pred_fn = "pred_img.jpeg";
    
    // std::string orig_fn =  "source_img.pgm";
    // std::string noisy_fn = "noisy_img.pgm";
    // std::string pred_fn = "pred_img.pgm";
    
    
    
    // Parse command line arguments --------------------------------------------->
    graphlab::command_line_options clopts("Loopy BP image denoising");
    clopts.attach_option("bound",
                         &BOUND, BOUND,
                         "Residual termination bound");
    clopts.attach_option("ncolors",
                         &NCOLORS, NCOLORS,
                         "The number of colors in the noisy image");
    clopts.attach_option("sigma",
                         &SIGMA, SIGMA,
                         "Standard deviation of noise.");
    clopts.attach_option("nrows",
                         &nrows, nrows,
                         "The number of rows in the noisy image");
    clopts.attach_option("ncols",
                         &ncols, ncols,
                         "The number of columns in the noisy image");
    clopts.attach_option("lambda",
                         &lambda, lambda,
                         "Smoothness parameter (larger => smoother).");
    clopts.attach_option("smoothing",
                         &smoothing, smoothing,
                         "Options are {square, laplace}");
    clopts.attach_option("orig",
                         &orig_fn, orig_fn,
                         "Original image file name.");
    clopts.attach_option("noisy",
                         &noisy_fn, noisy_fn,
                         "Noisy image file name.");
    clopts.attach_option("pred",
                         &pred_fn, pred_fn,
                         "Predicted image file name.");
    
    ///! Initialize control plain using mpi
    graphlab::mpi_tools::init(argc, argv);
    const bool success = clopts.parse(argc, argv);
    if(!success) {
        graphlab::mpi_tools::finalize();
        return EXIT_FAILURE;
    }
    
    ///! Create a distributed control object 
    graphlab::distributed_control dc;
    ///! display settings  
    if(dc.procid() == 0) {
        std::cout << "ncpus:          " << clopts.get_ncpus() << std::endl
        << "bound:          " << BOUND << std::endl
        << "colors:         " << NCOLORS << std::endl
        << "nrows:           " << nrows << std::endl
        << "ncols:           " << ncols << std::endl
        << "sigma:          " << SIGMA << std::endl
        << "lambda:         " << lambda << std::endl
        << "smoothing:      " << smoothing << std::endl
        << "scheduler:      " << clopts.get_scheduler_type() << std::endl
        << "orig_fn:        " << orig_fn << std::endl
        << "noisy_fn:       " << noisy_fn << std::endl
        << "pred_fn:        " << pred_fn << std::endl;
    }
    
    
    
    
    // Create synthetic images -------------------------------------------------->
    std::cout << "Creating a synthetic noisy image." << std::endl;
    graph_type graph(dc, clopts);
    create_synthetic_cluster_graph(dc, graph, nrows, ncols);
    std::cout << "Finalizing the graph." << std::endl;
    graph.finalize();
    
    std::cout << "Collect the noisy image. " << std::endl;
    merge_reduce_type obs_image = 
    graph.map_reduce_vertices<merge_reduce_type>(obs_map_function);
    std::cout << "saving the noisy image." << std::endl;
    if(dc.procid() == 0) {
        save_image(nrows, ncols, obs_image.values, noisy_fn);
    }
    
    // Initialze the edge factor ----------------------------------------------->
    std::cout << "Initializing shared edge factor. " << std::endl;
    // dummy variables 0 and 1 and num_rings by num_rings
    initialize_theta_ij(smoothing, lambda);
    if(dc.procid() == 0) std::cout << THETA_ij << std::endl;
    
    // Create the engine -------------------------------------------------------->
    std::cout << "Creating the engine. " << std::endl;
    engine_type engine(dc, graph, clopts);
    
    std::cout << "Scheduling all vertices" << std::endl;
    engine.signal_all();
    std::cout << "Starting the engine" << std::endl;
    engine.start();
    const float runtime = engine.elapsed_time();
    size_t update_count = engine.num_updates();
    std::cout << "Finished Running engine in " << runtime 
    << " seconds." << std::endl
    << "Total updates: " << update_count << std::endl
    << "Efficiency: " << (double(update_count) / runtime)
    << " updates per second "
    << std::endl;  
    
    
    // Saving the output -------------------------------------------------------->
    std::cout << "Saving the predicted image" << std::endl;
    std::cout << "Collect the noisy image. " << std::endl;
    merge_reduce_type pred_image = 
    graph.map_reduce_vertices<merge_reduce_type>(pred_map_function);
    std::cout << "saving the pred image." << std::endl;
    if(dc.procid() == 0) {
        save_image(nrows, ncols, pred_image.values, pred_fn);
    }
    
    std::cout << "Done!" << std::endl;
    graphlab::mpi_tools::finalize();
    return EXIT_SUCCESS;
} // End of main



