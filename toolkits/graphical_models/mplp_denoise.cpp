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

#include <cv.h>
#include <highgui.h>  



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
graphlab::mutex mutex;
//mutex.lock();
double LPval = 0;
double MAPval = 0;
double MAPrepval = 0;
//mutex.unlock();

// Shared base edge potential
matrix THETA_ij; 

// keep track of predictions at each node
vector PRED_COLOR;

// check if all nodes are visited
//Eigen::Matrix<graphlab::atomic<int>, Eigen::Dynamic,1> vinit;
//Eigen::Matrix<graphlab::atomic<int>, Eigen::Dynamic,1> vapply;
//vector vinit;
//vector vapply;
//std::vector<graphlab::atomic<int> > vinit;
//std::vector<graphlab::atomic<int> > vapply;
std::vector<int> vinit;
std::vector<int> vapply;



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
    
    // current maximizers of reparameterized theta_i, theta_j and theta_IJ
    int maxI, maxJ, maxIJ_i, maxIJ_j;
    
    // current contribution to LP dual value
    double vali, valj, valij;
    // current contribution to MAP value
    double pvali, pvalj, pvalij;
    // current contribution to MAPrep value
    double prvali, prvalj, prvalij;
    
    // since variables i and j are present in multiple factors, this determines who owns them
    bool iowner, jowner; 
    
    /** dual variables being optimized (or messages) */
    vector delf_i, delf_j;  

    // constructor
    vertex_data(): i(-1), j(-1), deg_i(0), deg_j(0), 
    obs_color_i(-1), obs_color_j(-1), 
    pred_color_i(0), pred_color_j(0),
    vali(0), valj(0), valij(0),
    pvali(0), pvalj(0), pvalij(0), 
    prvali(0), prvalj(0), prvalij(0), 
    iowner(false), jowner(false)
    { }
    
    void save(graphlab::oarchive& arc) const 
    {
        arc << i << j 
        << deg_i << deg_j 
        << obs_color_i << obs_color_j 
        << pred_color_i << pred_color_j 
        << maxI << maxJ << maxIJ_i << maxIJ_j
        << vali << valj << valij 
        << pvali << pvalj << pvalij 
        << prvali << prvalj << prvalij 
        << iowner << jowner
        << delf_i << delf_j;
    }
    void load(graphlab::iarchive& arc) 
    {
        arc >> i >> j 
        >> deg_i >> deg_j
        >> obs_color_i >> obs_color_j 
        >> pred_color_i >> pred_color_j 
        >> maxI >> maxJ >> maxIJ_i >> maxIJ_j
        >> vali >> valj >> valij
        >> pvali >> pvalj >> pvalij
        >> prvali >> prvalj >> prvalij
        >> iowner >> jowner
        >> delf_i >> delf_j;
    }
}; // End of vertex data


// /**
//  * The data associated with a pair of factors in a pairwise MRF
//  */
//struct edge_data : public graphlab::IS_POD_TYPE 
//{
//    // primal labelling; We assume pairwise factors, so intersection has
//    // a single node
//    int pred_color;
//    
//    // current contribution to LP dual value
//    double dval;
//    // current contribution to MAP value
//    double pval;
//    // current contribution to MAPrep value
//    double prval;
//    
//    edge_data():
//    pred_color(0),
//    dval(0), pval(0), prval(0)
//    {}
//     
//    void save(graphlab::oarchive& arc) const 
//    {
//        arc << pred_color
//        << dval << pval << prval; 
//    }
//    void load(graphlab::iarchive& arc) 
//    {
//        arc >> pred_color
//        >> dval >> pval >> prval;
//    }
//}; // End of edge data
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
     * This function is now called in the main by invoking:
     * engine.transform_vertices(mplp_vertex_program::init)
     */
    static void init_vertex_data(icontext_type& context, vertex_type& vertex)
    { 
        vertex_data& vdata = vertex.data();
        
        // Create zero messages
        vdata.delf_i = vector::Zero(NCOLORS);
        vdata.delf_j = vector::Zero(NCOLORS);
        
        // create temporary node potentials
        vector theta_i = make_unary_potential(vertex, 'i');
        vector theta_j = make_unary_potential(vertex, 'j');

        // if we own i
        if (vdata.iowner) 
        {
            // get dual contribution
            vdata.vali = theta_i.maxCoeff(&vdata.maxI); 

            // get primal contribution
            vdata.pred_color_i = vdata.maxI;
            vdata.pvali = vdata.vali;
            
            // also update the global copy
            PRED_COLOR[vdata.i] = vdata.pred_color_i;

            // get rep primal contribution
            vdata.prvali = vdata.pvali;
        }
        else // if we don't own then just copy over the global predicted color
            vdata.pred_color_i = PRED_COLOR[vdata.i];

        // if we own j
        if (vdata.jowner) 
        {
            // get dual contribution
            vdata.valj = theta_j.maxCoeff(&vdata.maxJ); 
            
            // get primal contribution
            vdata.pred_color_j = vdata.maxJ;
            vdata.pvalj = vdata.valj;
            
            // also update the global copy
            PRED_COLOR[vdata.j] = vdata.pred_color_j;

            // get rep primal contribution
            vdata.prvalj = vdata.pvalj;
        }
        else 
            vdata.pred_color_j = PRED_COLOR[vdata.j];
        
        // we always own edge i,j
        vdata.valij = THETA_ij.maxCoeff(&vdata.maxIJ_i,&vdata.maxIJ_j); 
        vdata.pvalij = THETA_ij(vdata.pred_color_i, vdata.pred_color_j);
        vdata.prvalij = vdata.pvalij;
             
        mutex.lock();
        LPval += vdata.vali; LPval += vdata.valj; LPval += vdata.valij;
        MAPval += vdata.pvali; MAPval += vdata.pvalj; MAPval += vdata.pvalij;
        MAPrepval += vdata.prvali; MAPrepval += vdata.prvalj; MAPrepval += vdata.prvalij;
        mutex.unlock();

        // debug code to check in all nodes are inited
        if (vinit[vertex.id()] == 0)
            vinit[vertex.id()] = 1;
    }
    
    /**
     * Recv message is called by the engine to receive a message to this
     * vertex program.  The vertex program can use this to initialize
     * any state before entering the gather phase.  If the vertex
     * program does not implement this function then the default
     * implementation (NOP) is used.
     */
    // void init(icontext_type& context, const vertex_type& vertex, 
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
               const gather_type& sum) 
    {
        // Make sure this vertex has neighbours. Everyone should have neighbours
        ASSERT_GT(vertex.num_in_edges() + vertex.num_out_edges(), 0);
        
        vertex_data& vdata = vertex.data();  
        vector theta_i = make_unary_potential(vertex, 'i');
        vector theta_j = make_unary_potential(vertex, 'j');
                
        ASSERT_EQ(THETA_ij.rows(), theta_i.size());
        ASSERT_EQ(THETA_ij.rows(), sum.delf_i.size());
        ASSERT_EQ(THETA_ij.cols(), theta_j.size());
        ASSERT_EQ(THETA_ij.cols(), sum.delf_j.size());   
        
        // debug code to check in all nodes are applied
        if (vapply[vertex.id()] == 0)
            vapply[vertex.id()] = 1;
        
        
        ////////////////////////////////////////////
        // Update outgoing messages (coordinate descent)
        
        // Backup the old prediction
        const vector old_delf_i = vdata.delf_i;
        const vector old_delf_j = vdata.delf_j;
        
        // Update del fi
        vdata.delf_i = -(theta_i + sum.delf_i)/2 + 
        (THETA_ij + (theta_j + sum.delf_j).transpose().replicate(THETA_ij.rows(),1)).
        rowwise().maxCoeff()/2;
        // Update del fj
        vdata.delf_j = -(theta_j + sum.delf_j)/2 + 
        ((THETA_ij + (theta_i + sum.delf_i).replicate(1,THETA_ij.cols())).
         colwise().maxCoeff()).transpose()/2;
        
        ////////////////////////////////////////////
        // Compute contributions to dual, primal and rep primal
        
        // Remove contribution of old labels from LPval
        double LPremove=0, MAPremove=0, MAPrepremove=0;
        LPremove += vdata.vali; LPremove += vdata.valj; LPremove += vdata.valij;
        MAPremove += vdata.pvali; MAPremove += vdata.pvalj; MAPremove += vdata.pvalij;
        MAPrepremove += vdata.prvali; MAPrepremove += vdata.prvalj; MAPrepremove += vdata.prvalij;

        // Update dual, primal and rep primal contributions. 
        // TODO: if primal labelling changes at a node we own, update it's edge potential too
        if (vdata.iowner)
        {
            // reparameterized node potential
            vector thetarep_i = theta_i + sum.delf_i + vdata.delf_i;
            
            vdata.vali = thetarep_i.maxCoeff(&vdata.maxI);
            
            vdata.pred_color_i = vdata.maxI;
            vdata.pvali = theta_i[vdata.pred_color_i];
            
            PRED_COLOR[vdata.i] = vdata.pred_color_i;
            
            vdata.prvali = thetarep_i[vdata.pred_color_i];
        } 
        else
            vdata.pred_color_i = PRED_COLOR[vdata.i];
        if (vdata.jowner)
        {
            // reparameterized node potential
            vector thetarep_j = theta_j + sum.delf_j + vdata.delf_j;
            
            vdata.valj = thetarep_j.maxCoeff(&vdata.maxJ);
            
            vdata.pred_color_j = vdata.maxJ;
            vdata.pvalj = theta_j[vdata.pred_color_j];
            
            PRED_COLOR[vdata.j] = vdata.pred_color_j;

            vdata.prvalj = thetarep_j[vdata.pred_color_j];
        }
        else
            vdata.pred_color_j = PRED_COLOR[vdata.j];
        
        // We always own edge i,j
        matrix thetarep_ij = THETA_ij - (vdata.delf_i.replicate(1,THETA_ij.cols()))
                            - (vdata.delf_j.transpose().replicate(THETA_ij.rows(),1));
        
        vdata.valij = thetarep_ij.maxCoeff(&vdata.maxIJ_i, &vdata.maxIJ_j);
        vdata.pvalij = THETA_ij(vdata.pred_color_i, vdata.pred_color_j);
        vdata.prvalij = thetarep_ij(vdata.pred_color_i, vdata.pred_color_j);
        
        mutex.lock();
        LPval -= LPremove; MAPval -= MAPremove; MAPrepval -= MAPrepremove;
        LPval += vdata.vali; LPval += vdata.valj; LPval += vdata.valij;
        MAPval += vdata.pvali; MAPval += vdata.pvalj; MAPval += vdata.pvalij;
        MAPrepval += vdata.prvali; MAPrepval += vdata.prvalj; MAPrepval += vdata.prvalij;
        mutex.unlock();
        
        ////////////////////////////////////////////
        // Debugging printing and residuals
        
        //std::cout << vertex.id() << ": " << vdata.i << "," << vdata.j << "\n";
        if (0) // (vdata.i == 0 )
        {
            mutex.lock();
            std::cout << "Applying at vertex: " << vertex.id() << "(" << vdata.i << "," << vdata.j << ")\n";
            
            std::cout << LPval << "," << MAPval << "," << MAPrepval << "\t" ;        
            int vinitsum = 0, vapplysum = 0;
            //std::vector<graphlab::atomic<int> >::iterator it;
            std::vector<int>::iterator it;
            for (it = vinit.begin(); it < vinit.end(); ++it)
                vinitsum += *it;
            for (it = vapply.begin(); it < vapply.end(); ++it)
                vapplysum += *it;
            
            // if all vertices have been visited start counting again
            if (vapply.size() == vapplysum)
                for (int i=0; i!=vapply.size(); ++i)
                    vapply[i] = 0;
            
            std::cout << "Verted Id: " << vertex.id() << " "  << vdata.i << "," << vdata.j << " "; 
            std::cout << "Inited: " << vinitsum << " Applied: " << vapplysum << "\n";
            if (vinit.size() == vinitsum)
                std::cout << "Restarting counting of apply\n";
            std::cout.flush();
            mutex.unlock();
        }
        
        if (0)//vdata.i == 1) 
        {
            std::cout << "\n\n";
            
            std::cout << "Pairwise Potential\n" << THETA_ij << "\n\n";
            std::cout << "theta_ij reparameterized: \n" <<         
            (THETA_ij - vdata.delf_i.replicate(1,THETA_ij.cols()) 
             - vdata.delf_j.transpose().replicate(THETA_ij.rows(),1)   )  << "\n\n";
            std::cout << "maxIJ_i: " << vdata.maxIJ_i << " maxIJ_j: " << vdata.maxIJ_j << "\n\n";
            
            std::cout << "thetai \n" << theta_i << "\n\n";
            std::cout << "sum of incomming messages into i\n" << sum.delf_i << "\n\n";
            std::cout << "outgoing message to i\n" << vdata.delf_i << "\n\n";
            std::cout << " Reparamterized thetai\n" << (theta_i + sum.delf_i + vdata.delf_i) << "\n\n";
            std::cout << "maxI: " << vdata.maxI << "\n\n";
            
            std::cout << "thetaj \n" << theta_j << "\n\n";
            std::cout << "sum of incomming messages into j\n" << sum.delf_j << "\n\n";
            std::cout << "outgoing message to j\n" << vdata.delf_j << "\n\n";
            std::cout << " Reparamterized thetaj\n" << (theta_j + sum.delf_j + vdata.delf_j) << "\n\n";
            std::cout << "maxJ: " << vdata.maxJ << "\n\n";
            
            std::cout << "thetaij + j message\n" << (THETA_ij + sum.delf_j.transpose().replicate(THETA_ij.rows(),1))/2 << "\n\n";
            
            std::cout << (THETA_ij + sum.delf_j.transpose().replicate(THETA_ij.rows(),1)).
            rowwise().maxCoeff()/2 << std::endl << std::endl;
            
            std::cout << "thetaij + i message\n" << (THETA_ij + sum.delf_i.replicate(1,THETA_ij.cols()))/2 << "\n\n";
            std::cout << 	  ((THETA_ij + sum.delf_i.replicate(1,THETA_ij.cols())).
                               colwise().maxCoeff()).transpose()/2 << "\n\n";
            
            std::cout << "Old del_fi\n" << vdata.delf_i << "\n\n";
            std::cout << "Old del_fj\n" << vdata.delf_j << "\n\n";
            std::cout << "New del_fi\n" << -(theta_i + sum.delf_i)/2 + 
            (THETA_ij + (theta_j + sum.delf_j).transpose().replicate(THETA_ij.rows(),1)).
            rowwise().maxCoeff()/2 << "\n\n";
            
            std::cout << "New del_fj\n" << -(theta_j + sum.delf_j)/2 + 
            ((THETA_ij + (theta_i + sum.delf_i).replicate(1,THETA_ij.cols())).
             colwise().maxCoeff()).transpose()/2 << "\n\n";
            
            getchar();
        }
        
        // const double residual = (vdata.delf_i - old_delf_i).cwiseAbs().sum() +
        // (vdata.delf_j - old_delf_j).cwiseAbs().sum();
        
        //priority = residual;
        priority = LPval - MAPval;
        //std::cout << "priority: " << priority << std::endl;
        //std::cout << LPval << std::endl;
        //test code; for now, only run 1 iteration
        //priority = 0;
    } // end of apply
    
    /**
     * Since the MRF is undirected we will use all edges for gather and
     * scatter
     */
    edge_dir_type scatter_edges(icontext_type& context,
                                const vertex_type& vertex) const { 
        //return priority < BOUND? graphlab::NO_EDGES : graphlab::ALL_EDGES; 
        return graphlab::ALL_EDGES;
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
    static vector make_unary_potential(const vertex_type& vertex, 
                                       const char varid) {
        vector potential(NCOLORS);
        const double obs = varid == 'i'? 
        vertex.data().obs_color_i : vertex.data().obs_color_j;
        const double sigmaSq = SIGMA*SIGMA;
        for(int pred = 0; pred < potential.size(); ++pred) {
            potential(pred) = -(obs - pred)*(obs - pred) / (2.0 * sigmaSq);
        }
        //potential /= std::abs(potential.sum());
        
        //float tmp = potential.minCoeff();
        //potential.array() -= tmp; // (float) potential.minCoeff();
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




/////////////////////////////////////////////////////////////////////////////////////
// Aggregator functions to compute primal & dual values
double get_energy_fun(mplp_vertex_program::icontext_type& context, const mplp_vertex_program::vertex_type& vertex) 
{
    double tmp = 0;
    
    const vertex_data &vdata = vertex.data();
    
    if (vdata.iowner)
        tmp += vdata.vali;
    if (vdata.jowner)
        tmp += vdata.valj;
    tmp += vdata.valij;

    return tmp;
}

void finalize_fun(mplp_vertex_program::icontext_type& context, double total) 
{
    if(context.procid() == 0) 
        std::cout << "Dual value: " << total << std::endl;
}



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

graphlab::vertex_id_type factor_ind2(size_t rows, size_t cols,
                                    size_t i, size_t j) {
    if(i > j) std::swap(i,j);
    
    if (j == (i+1)) // horizontal edge
        return i - std::floor(i/cols);
    else if (j == (i+cols))
        return rows*(cols-1) + i;
    else
      std::cout << "Problem ";
      //ASSERT_TRUE(false);
    return 0;
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
        //int nedges = 2*rows*cols -rows-cols;
        int nedges = factor_ind2(rows,cols,rows*cols-cols,rows*cols);
        //vinit = vector::Zero(2*rows*cols -rows-cols);
        //vapply = vector::Zero(2*rows*cols -rows-cols);
        vinit.clear(); vinit.resize(nedges, 0);
        vapply.clear(); vapply.resize(nedges,0);

        int max_vid = 0; 
        
        PRED_COLOR = vector::Zero(rows*cols);
        int ownercount = 0; 
        
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
                //                for(int pred = 0; pred < potential.size(); ++pred)                 
                //                    potential[pred] /= sum;
                //                ne << true_pixels[pixel_ind(rows,cols,r,c)] << " " <<  obs << " " << potential << std::endl;
                //                ne << 4 - int((r==0)||(r==(rows-1))) - int((c==0)||(c==(cols-1))) << " " << potential << std::endl;
                ne << potential << std::endl;
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
                    graph.add_vertex(factor_ind2(rows,cols,vdata.i,vdata.j), vdata);
                    
                    // temp code
                    max_vid = std::max((int)factor_ind2(rows,cols,vdata.i,vdata.j), max_vid); 
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

                    vdata.iowner = true; // give i-ownership to horizontal edges
                    ++ownercount;
                    if ((c+1)==(cols-1)) // and j-ownership too if last node in this row
                    {
                        vdata.jowner = true;
                        ++ownercount;
                    }
                    graph.add_vertex(factor_ind2(rows,cols,vdata.i,vdata.j), vdata);
                    
                    // temp code
                    max_vid = std::max((int)factor_ind2(rows,cols,vdata.i,vdata.j), max_vid); 
                }
                // Compute all the factors that contain this pixel
                nbrs.clear();
                if(r+1 < rows)
                    nbrs.push_back(factor_ind2(rows,cols,
                                              pixel_ind(rows,cols,r,c),
                                              pixel_ind(rows,cols,r+1,c)));
                if(r-1 < rows)
                    //if(r-1 >= 0)
                    nbrs.push_back(factor_ind2(rows,cols,
                                              pixel_ind(rows,cols,r-1,c),
                                              pixel_ind(rows,cols,r,c)));
                if(c+1 < cols)
                    nbrs.push_back(factor_ind2(rows,cols,
                                              pixel_ind(rows,cols,r,c),
                                              pixel_ind(rows,cols,r,c+1)));
                if(c-1 < cols)
                    //if(c-1 >= 0)
                    nbrs.push_back(factor_ind2(rows,cols,
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
        std::cout << "Max vid fed into graphlab: " << max_vid << "\n";
        std::cout << "No. of owners: " << ownercount << "\n";
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
  std::cout << "NPixels: " << values.size() << std::endl;
  // determine the max and min colors
  float max_color = -std::numeric_limits<float>::max();
  float min_color =  std::numeric_limits<float>::max();
  foreach(pred_pair_type pair, values) {
    max_color = std::max(max_color, pair.second);
    min_color = std::min(min_color, pair.second);
  }

  cv::Mat img(cols, rows, CV_8UC1);
  foreach(pred_pair_type pair, values) {
    std::pair<int,int> coords = ind2sub(rows,cols, pair.first);
    float value = (pair.second - min_color) / (max_color - min_color);
    int color = 255 * value > 255 ? 255 : 255 * value;
    img.at<unsigned char>(coords.first, coords.second) = color;
  }
  cv::imwrite(fname, img);
}




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
    
    
//    size_t nrows = 200;
//    size_t ncols = 200;
    size_t nrows = 20;
    size_t ncols = 20;
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
    clopts.attach_option("bound", BOUND,
                         "Residual termination bound");
    clopts.attach_option("ncolors", NCOLORS,
                         "The number of colors in the noisy image");
    clopts.attach_option("sigma", SIGMA,
                         "Standard deviation of noise.");
    clopts.attach_option("nrows", nrows,
                         "The number of rows in the noisy image");
    clopts.attach_option("ncols", ncols,
                         "The number of columns in the noisy image");
    clopts.attach_option("lambda", lambda,
                         "Smoothness parameter (larger => smoother).");
    clopts.attach_option("smoothing", smoothing,
                         "Options are {square, laplace}");
    clopts.attach_option("orig", orig_fn,
                         "Original image file name.");
    clopts.attach_option("noisy", noisy_fn,
                         "Noisy image file name.");
    clopts.attach_option("pred", pred_fn,
                         "Predicted image file name.");
    
    ///! Initialize control plain using mpi
    graphlab::mpi_tools::init(argc, argv);
    const bool success = clopts.parse(argc, argv);
    if(!success) {
        clopts.print_description();
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

    engine.add_vertex_aggregator<double>("energy", get_energy_fun, finalize_fun);
    engine.aggregate_periodic("energy", 3); // run every 3 seconds

    engine.transform_vertices(mplp_vertex_program::init_vertex_data);

    std::cout << "Scheduling all vertices" << std::endl;
    engine.signal_all();
    std::cout << "Starting the engine" << std::endl;
    engine.start();
    const float runtime = engine.elapsed_seconds();
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



