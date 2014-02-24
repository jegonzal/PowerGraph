/**  
 * Copyright (c) 2013 GraphLab Inc.
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
 * Implementation of the adpredictor algorithm as given in the paper:
 * Web-Scale Bayesian Click-Through Rate Prediction for Sponsored Search Advertising in Microsoft’s Bing Search Engine
 * Thore Graepel, Joaquin Quinonero Candela, Thomas Borchert, and Ralf Herbrich
 * ICML 2010
 * Implemented by Danny Bickson, GraphLab Inc.
 *
 */


#include <graphlab.hpp>
#include <graphlab/util/timer.hpp>
#include <graphlab/warp.hpp>
#include "stats.hpp"
#include "cdf.hpp"

//when using negative node id range, we are not allowed to use
//0 and 1 so we add 2.
const static int SAFE_NEG_OFFSET=2;
const double pi = 3.14159265;
const double gaussian_normalization = 1/sqrt(2 * pi);
double beta = 1;
bool debug = false;


enum data_role_type { TRAIN = 0, VALIDATE = 1, PREDICT =2 };

data_role_type mode;
/** 
 */
struct vertex_data : graphlab::IS_POD_TYPE{
	int y;
	float xT_mu; 
	float sigma;
	float predict;
	float err;
	float likelihood;
	float weights;
	data_role_type type;

	vertex_data() {
		xT_mu = 0;
		y = 0;
		sigma  = 1;
		predict = 0;
		err = 0;
		likelihood = 0;
		weights = 0;
		type = TRAIN;
	}


}; // end of vertex data


/**
 * \brief The edge data stores the entry in the matrix.
 *
 * In addition the edge data adpredictoro stores the most recent error estimate.
 */
struct edge_data : public graphlab::IS_POD_TYPE {
	/**
	 * \brief The type of data on the edge;
	 *
	 * \li *Train:* the observed value is correct and used in training
	 * \li *Validate:* the observed value is correct but not used in training
	 * \li *Predict:* The observed value is not correct and should not be
	 *        used in training.
	 */

	/** \brief the observed value for the edge */
	float x_ij;

	/** \brief The train/validation/test designation of the edge */
	data_role_type role;

	/** \brief basic initialization */
	edge_data(float x_ij = 1, data_role_type role = TRAIN) :
		x_ij(x_ij), role(role) { }

}; // end of edge data



/**
 * \brief The graph type is defined in terms of the vertex and edge
 * data.
 */ 
typedef graphlab::distributed_graph<vertex_data, edge_data> graph_type;






/* compute v(t) according to equation (9) left */
double v(double t){
	return gaussian_normalization * exp(-t*t/2) / phi(t);
}

/* compute w(t) according to equation (9) right */
double w(double t){
	double vt = v(t);
	return vt * (vt+t);
}


struct gather_type: public graphlab::IS_POD_TYPE{
	float sigma;
	float mu;
	float mult_sigma;

	gather_type(){ sigma = 0; mu = 0; mult_sigma = 1; }
	gather_type& operator+=(const gather_type& other) {
		sigma += other.sigma;
		mu += other.mu;
		mult_sigma *= other.mult_sigma;
		return *this;
	}
};

/** compute probability for click as given in equation (2) */
float ctr_predict( const vertex_data& data, 
		const float rating, 
		double & prediction, 
		void * extra = NULL){

	assert(beta > 0);
	prediction = data.xT_mu;
	double prob = phi(data.xT_mu * data.y / beta);
	if (debug)
		//std::cout<<"prediction: " << prediction << " y: " << data.y << std::endl;
		printf("prediction %12.8lf y: %d \n", prediction, data.y);
	return prob; 
}



gather_type adpredictor_map(graph_type::edge_type edge, graph_type::vertex_type other) {
	gather_type ret;
	assert(edge.data().x_ij  == 1);
	/* compute equation (6) */
	ret.sigma = edge.data().x_ij * other.data().sigma;
	ret.mu = edge.data().x_ij * other.data().xT_mu;			
	return ret;
}

// the function arguments of the combiner must match the return type of the
// map function.
void adpredictor_combine2(gather_type &a, const gather_type & b, const vertex_data unused) {
	a.mu += b.mu;
	a.mult_sigma *= b.mult_sigma;
}



gather_type adpredictor_map2(graph_type::edge_type edge, graph_type::vertex_type other, vertex_data  vertex){
	gather_type ret;
	assert(vertex.sigma > 0);
	assert(other.data().y == -1 || other.data().y == 1);
	double product = other.data().y * other.data().xT_mu / sqrt(other.data().sigma);
	//assert(product > 0);
	ret.mu = (other.data().y * edge.data().x_ij * vertex.sigma / sqrt(other.data().sigma))  * v(product);
	double factor = 1.0 - (edge.data().x_ij * vertex.sigma / other.data().sigma)*w(product);
	assert(factor > 0);
	ret.sigma = factor;
	return ret;
}

void adpredictor_update(graph_type::vertex_type vertex) {
	//go over all row nodes
	if ( vertex.num_out_edges() > 0){
                if (debug) printf("Entered vertex %llu role %d \n", vertex.id(), vertex.data().type);
		if (vertex.data().type == TRAIN){
			vertex_data & row = vertex.data(); 
			row.likelihood = 0;
			row.err = 0;
			assert(row.y == -1 || row.y == 1);
			assert(beta > 0);

			if (debug)
				std::cout<<"Entered item " << vertex.id() << " y: " << row.y << std::endl;
			row.sigma = beta*beta;
			row.xT_mu = 0;

			gather_type sum = graphlab::warp::map_reduce_neighborhood<gather_type>(vertex, graphlab::OUT_EDGES, adpredictor_map);
			row.sigma = sum.sigma;
			row.xT_mu = sum.mu;

			double prediction;
			double ret = ctr_predict(row, row.y, prediction);
			double predicted_target = prediction < 0 ? -1: 1;
			if ((predicted_target == -1  && row.y == 1) || (predicted_target == 1 && row.y == -1))
				row.err += 1.0;  
			if (debug)
				std::cout<<"Prediction was: " << prediction << " real value: " << row.y << std::endl;
			row.likelihood += ret;

			assert(row.sigma > 0);
		}
		else if (vertex.data().type == VALIDATE){
			vertex_data & row = vertex.data(); 
			row.likelihood = 0;
			row.err = 0;
			assert(row.y == -1 || row.y == 1);
			gather_type sum = graphlab::warp::map_reduce_neighborhood<gather_type>(vertex, graphlab::OUT_EDGES, adpredictor_map);
			double predict = sum.mu > 0 ? 1 : -1;                       
			if (predict != row.y)
				row.err++;
		}
                else assert(false);

	}
}

void adpredictor_update2(graph_type::vertex_type vertex) {
	if (vertex.num_in_edges() > 0){
		gather_type sum = graphlab::warp::map_reduce_neighborhood<gather_type>(vertex, graphlab::IN_EDGES, vertex.data(), adpredictor_map2,adpredictor_combine2);
		vertex.data().sigma *= sum.mult_sigma;
		vertex.data().xT_mu += sum.mu;
	}
}
gather_type count_vertices(const graph_type::vertex_type& vertex) {
	gather_type ret;
	if (vertex.data().type == TRAIN){
		ret.mu = 1;
	}
	else if (vertex.data().type == VALIDATE){
		ret.sigma = 1;
	}
	return ret;
}

struct model_saver {
	typedef graph_type::vertex_type vertex_type;
	typedef graph_type::edge_type   edge_type;
	/* save the linear model, using the format:
	 */
	std::string save_vertex(const vertex_type& vertex) const {
		std::stringstream strm;
		strm << vertex.id() << " " << vertex.data().xT_mu << " " << vertex.data().sigma << std::endl;
		return strm.str();
	}
	std::string save_edge(const edge_type& edge) const {
		return "";
	}
}; // end of prediction_saver


/**
 * \brief The graph loader function is a line parser used for
 * distributed graph construction.
 */
inline bool graph_loader(graph_type& graph, 
		const std::string& filename,
		const std::string& line) {
	ASSERT_FALSE(line.empty()); 

	// Parse the line
	std::stringstream strm(line);
	float weight = 0;
	float label = 0;

	strm >> label;
	if (label != -1 && label != 1)
		logstream(LOG_FATAL)<<"Each line must have label -1 or 1 as the first item in the row. Row was : " << line << " label: " << label << std::endl;

	// Determine the role of the data
	data_role_type role = TRAIN;
	if(boost::ends_with(filename,".validate")) role = VALIDATE;
	else if(boost::ends_with(filename, ".predict")) role = PREDICT;

	int myid = rand();
	int num_vals = 0;
	while (strm.good()) {
		graphlab::vertex_id_type target;
		strm >> target;
		if (strm.fail()) break;
		char col;
		strm >> col;
		if (strm.fail()) break;
		strm >> weight;
		if (strm.fail()) break;
		if (weight != 1)
			logstream(LOG_FATAL)<<"Currently we support only binary edges. Line was: " << line << " in file: " << filename << std::endl;
		num_vals++;
		target = -(graphlab::vertex_id_type(target + SAFE_NEG_OFFSET));
		graph.add_edge(myid, target, edge_data(weight, role));
	}  

	if (num_vals == 0)
		logstream(LOG_FATAL)<<"Failed to load line: " << line << " in file: " << filename << std::endl;

	vertex_data data;
	data.y = label;
	data.type = role;
        if (debug) printf("Adding vertex %u with role %d\n", myid, role);
	graph.add_vertex(myid, data);
	return true; // successful load
} // end of graph_loader






int MAX_ITER = 5;

gather_type calc_error(const graph_type::vertex_type& vertex) {
	gather_type ret;
	if (mode == vertex.data().type){
		ret.mu = vertex.data().err;
		ret.sigma = vertex.data().likelihood;
	}
	return ret;
}


int main(int argc, char** argv) {
	global_logger().set_log_level(LOG_INFO);
	global_logger().set_log_to_console(true);

	// Parse command line options -----------------------------------------------
	const std::string description = 
		"adPredictor algorithm";
	graphlab::command_line_options clopts(description);
	std::string input_dir;
	std::string save_model;
	std::string exec_type = "synchronous";
	clopts.attach_option("matrix", input_dir,
			"The directory containing the matrix file");
	clopts.add_positional("matrix");
	clopts.attach_option("max_iter", MAX_ITER,
			"The maxumum number of udpates allowed for a vertex");
	clopts.attach_option("debug", debug, 
			"debug - additional verbose info"); 
	clopts.attach_option("save_model", save_model,
			"The prefix (folder and filename) to save predictions.");
	clopts.attach_option("beta", beta, "gaussian bandwidth");

	if(!clopts.parse(argc, argv) || input_dir == "") {
		std::cout << "Error in parsing command line arguments." << std::endl;
		clopts.print_description();
		return EXIT_FAILURE;
	}

	graphlab::mpi_tools::init(argc, argv);
	graphlab::distributed_control dc;

	dc.cout() << "Loading graph." << std::endl;
	graphlab::timer timer; 
	graph_type graph(dc, clopts);  
	graph.load(input_dir, graph_loader); 
	dc.cout() << "Loading graph. Finished in " 
		<< timer.current_time() << std::endl;
	dc.cout() << "Finalizing graph." << std::endl;
	timer.start();
	graph.finalize();
	dc.cout() << "Finalizing graph. Finished in " 
		<< timer.current_time() << std::endl;


	dc.cout() 
		<< "========== Graph statistics on proc " << dc.procid() 
		<< " ==============="
		<< "\n Num vertices: " << graph.num_vertices()
		<< "\n Num edges: " << graph.num_edges()
		<< "\n Num replica: " << graph.num_replicas()
		<< "\n Replica to vertex ratio: " 
		<< float(graph.num_replicas())/graph.num_vertices()
		<< "\n --------------------------------------------" 
		<< "\n Num local own vertices: " << graph.num_local_own_vertices()
		<< "\n Num local vertices: " << graph.num_local_vertices()
		<< "\n Replica to own ratio: " 
		<< (float)graph.num_local_vertices()/graph.num_local_own_vertices()
		<< "\n Num local edges: " << graph.num_local_edges()
		//<< "\n Begin edge id: " << graph.global_eid(0)
		<< "\n Edge balance ratio: " 
		<< float(graph.num_local_edges())/graph.num_edges()
		<< std::endl;

	dc.cout() << "Running adPredictor" << std::endl;
	dc.cout() << "(C) Code by Danny Bickson, GraphLab Inc. " << std::endl;
	dc.cout() << "Please send bug reports to danny.bickson@gmail.com" << std::endl;
	timer.start();

	gather_type edge_count = graph.map_reduce_vertices<gather_type>(count_vertices);
	dc.cout()<<"Training rows: " << edge_count.mu << " validation rows: " << edge_count.sigma << std::endl;
        if (edge_count.mu <= 0)
          logstream(LOG_FATAL)<< "Failed to read training data. Aborting" << std::endl;

	graphlab::timer mytimer; mytimer.start();

	for (int i = 0; i < MAX_ITER; ++i) {
		graphlab::warp::parfor_all_vertices(graph, adpredictor_update); 
		graphlab::warp::parfor_all_vertices(graph, adpredictor_update2); 
		mode = TRAIN;
		gather_type ret = graph.map_reduce_vertices<gather_type>(calc_error);
		dc.cout() << i << ") Log likelihood: " << std::setw(10) << ret.sigma << " Avg error: " << std::setw(10) << ret.mu/edge_count.mu << std::endl; 
		mode = VALIDATE;
		ret = graph.map_reduce_vertices<gather_type>(calc_error);
		dc.cout() << i << " Avg validation error: " << std::setw(10) << ret.mu/edge_count.sigma << std::endl; 

	}

	const double runtime = timer.current_time();
	dc.cout() << "----------------------------------------------------------"
		<< std::endl
		<< "Final Runtime (seconds):   " << runtime << std::endl;


	// Make predictions ---------------------------------------------------------
	if(!save_model.empty()) {
		std::cout << "Saving predictions" << std::endl;
		const bool gzip_output = false;
		const bool save_vertices = true;
		const bool save_edges = false;
		const size_t threads_per_machine = 1;
		//save the predictions
		graph.save(save_model, model_saver(), gzip_output, save_vertices, save_edges, threads_per_machine);
	}

	graphlab::mpi_tools::finalize();
	return EXIT_SUCCESS;
} // end of main



