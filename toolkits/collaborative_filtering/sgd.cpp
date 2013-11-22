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
 * \file
 * 
 * \brief The main file for the SGD matrix factorization algorithm.
 *
 * This file contains the main body of the SGD matrix factorization
 * algorithm. 
 */

#include <graphlab/util/stl_util.hpp>
#include <graphlab.hpp>

#include <Eigen/Dense>
#include "eigen_serialization.hpp"
#include <graphlab/macros_def.hpp>


typedef Eigen::VectorXd vec_type;
typedef Eigen::MatrixXd mat_type;

//when using negative node id range, we are not allowed to use
//0 and 1 so we add 2.
const static int SAFE_NEG_OFFSET=2;
static bool debug;
int iter = 0;

bool isuser(uint node){
	return ((int)node) >= 0;
}

/** 
 * \ingroup toolkit_matrix_pvecization
 *
 * \brief the vertex data type which contains the latent pvec.
 *
 * Each row and each column in the matrix corresponds to a different
 * vertex in the SGD graph.  Associated with each vertex is a pvec
 * (vector) of latent parameters that represent that vertex.  The goal
 * of the SGD algorithm is to find the values for these latent
 * parameters such that the non-zero entries in the matrix can be
 * predicted by taking the dot product of the row and column pvecs.
 */
struct vertex_data {
	/**
	 * \brief A shared "constant" that specifies the number of latent
	 * values to use.
	 */
	static size_t NLATENT;
	/** \brief The latent pvec for this vertex */
	vec_type pvec;

	int nupdates;
	/** 
	 * \brief Simple default constructor which randomizes the vertex
	 *  data 
	 */
	vertex_data() : nupdates(0) { if (debug) pvec = vec_type::Ones(NLATENT); else randomize(); } 
	/** \brief Randomizes the latent pvec */
	void randomize() { pvec.resize(NLATENT); pvec.setRandom(); }
	/** \brief Save the vertex data to a binary archive */
	void save(graphlab::oarchive& arc) const { 
		arc << nupdates << pvec;
	}
	/** \brief Load the vertex data from a binary archive */
	void load(graphlab::iarchive& arc) { 
		arc >> nupdates >> pvec;
	}
}; // end of vertex data


/**
 * \brief The edge data stores the entry in the matrix.
 *
 * In addition the edge data sgdo stores the most recent error estimate.
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
	enum data_role_type { TRAIN, VALIDATE, PREDICT  };

	/** \brief the observed value for the edge */
	float obs;

	/** \brief The train/validation/test designation of the edge */
	data_role_type role;

	/** \brief basic initialization */
	edge_data(float obs = 0, data_role_type role = PREDICT) :
		obs(obs), role(role) { }

}; // end of edge data


/**
 * \brief The graph type is defined in terms of the vertex and edge
 * data.
 */ 
typedef graphlab::distributed_graph<vertex_data, edge_data> graph_type;

#include "implicit.hpp"

stats_info count_edges(const graph_type::edge_type & edge){
	stats_info ret;

	if (edge.data().role == edge_data::TRAIN)
		ret.training_edges = 1;
	else if (edge.data().role == edge_data::VALIDATE)
		ret.validation_edges = 1;
	ret.max_user = (size_t)edge.source().id();
	ret.max_item = (-edge.target().id()-SAFE_NEG_OFFSET);
	return ret;
}


double extract_l2_error(const graph_type::edge_type & edge);


/**
 * \brief Given a vertex and an edge return the other vertex in the
 * edge.
 */
inline graph_type::vertex_type
get_other_vertex(graph_type::edge_type& edge, 
		const graph_type::vertex_type& vertex) {
	return vertex.id() == edge.source().id()? edge.target() : edge.source();
}; // end of get_other_vertex





/**
 *
 */
class gather_type {
	public:
		/**
		 * \brief Stores the current sum of nbr.pvec.transpose() *
		 * nbr.pvec
		 */

		/**
		 * \brief Stores the current sum of nbr.pvec * edge.obs
		 */
		vec_type pvec;
		/** \brief basic default constructor */
		gather_type() { }

		/**
		 * \brief This constructor computes XtX and Xy and stores the result
		 * in XtX and Xy
		 */
		gather_type(const vec_type& X) {
			pvec = X;
		} // end of constructor for gather type

		/** \brief Save the values to a binary archive */
		void save(graphlab::oarchive& arc) const { arc << pvec; }

		/** \brief Read the values from a binary archive */
		void load(graphlab::iarchive& arc) { arc >> pvec; }  

		/** 
		 */
		gather_type& operator+=(const gather_type& other) {
			if (pvec.size() == 0){
				pvec = other.pvec;
				return *this;
			}
			else if (other.pvec.size() == 0)
				return *this;
			pvec += other.pvec;
			return *this;
		} // end of operator+=

}; // end of gather type

typedef vec_type message_type;

bool isuser_node(const graph_type::vertex_type& vertex){
	return isuser(vertex.id());
}

/**
 * SGD vertex program type
 */ 
class sgd_vertex_program : 
	public graphlab::ivertex_program<graph_type, gather_type,
	message_type> {
		public:
			/** The convergence tolerance */
			static double TOLERANCE;
			static double LAMBDA;
			static double GAMMA;
			static double MAXVAL;
			static double MINVAL;
			static double STEP_DEC;
			static bool debug;
			static size_t MAX_UPDATES;
			vec_type pmsg;

			void save(graphlab::oarchive& arc) const { 
				arc << pmsg;
			}
			/** \brief Load the vertex data from a binary archive */
			void load(graphlab::iarchive& arc) { 
				arc >> pmsg;
			}

			/** The set of edges to gather along */
			edge_dir_type gather_edges(icontext_type& context, 
					const vertex_type& vertex) const { 
				return graphlab::ALL_EDGES; 
			}; // end of gather_edges 


			gather_type gather(icontext_type& context, const vertex_type& vertex, 
					edge_type& edge) const {

				vec_type delta, other_delta;
				//this is user node
				if (vertex.num_in_edges() == 0){
					vertex_type my_vertex(vertex);
					//get a copy of the item node
					vertex_type other_vertex(get_other_vertex(edge, vertex));
					//compute the current prediction by computing a dot production of user and item nodes
					double pred = vertex.data().pvec.dot(other_vertex.data().pvec);
					//truncte predictions into allowed range
					pred = std::min(pred, sgd_vertex_program::MAXVAL);
					pred = std::max(pred, sgd_vertex_program::MINVAL); 
					//compute the prediction error 
					const float err = edge.data().obs - pred;
					if (debug)
						std::cout<<"entering edge " << (int)edge.source().id() << ":" << (int)edge.target().id() << " err: " << err << " rmse: " << err*err <<std::endl;
					if (std::isnan(err))
						logstream(LOG_FATAL)<<"Got into numeric errors.. try to tune step size and regularization using --lambda and --gamma flags" << std::endl;

					//for training edges, update the linear model
					if (edge.data().role == edge_data::TRAIN){
						//compute the change in gradient for this user node
						delta = GAMMA*(err*other_vertex.data().pvec - LAMBDA*vertex.data().pvec);
						//compute the change in gradient for this item node
						other_delta = GAMMA*(err*vertex.data().pvec - LAMBDA*other_vertex.data().pvec);

						//heuristic: update the current gradient with the change (this change is discarded when this function exists)
						//my_vertex.data().pvec += delta;
						//other_vertex.data().pvec += other_delta;
						if (debug)
							std::cout<<"new val:" << (int)edge.source().id() << ":" << (int)edge.target().id() << " U " << my_vertex.data().pvec.transpose() << " V " << other_vertex.data().pvec.transpose() << std::endl;
						//send the delta gradient for the item node to be updated in the next iteration
						if(std::fabs(err) > TOLERANCE && other_vertex.data().nupdates < MAX_UPDATES) 
							context.signal(other_vertex, other_delta);
					}
				} 
				return gather_type(delta);
			} // end of gather function

			void init(icontext_type& context,
					const vertex_type& vertex,
					const message_type& msg) {
				//if this is an item node, store the change in the gradient (sum of changes) to be
				//applied in the apply() function
				if (vertex.num_in_edges() > 0){
					pmsg = msg;
				}
			}

			void apply(icontext_type& context, vertex_type& vertex,
					const gather_type& sum) {

				vertex_data& vdata = vertex.data(); 
				//this is a user node, update the gradient using the comulative sum of gradient updates computed in gather
				if (sum.pvec.size() > 0){
					vdata.pvec += sum.pvec; 
					assert(vertex.num_in_edges() == 0);
				}
				//if this is an item node, update the gradient using the received sum from the init() function
				else if (pmsg.size() > 0){
					vdata.pvec += pmsg;
					assert(vertex.num_out_edges() == 0); 
				}
				++vdata.nupdates;
			} // end of apply

			/** The edges to scatter along */
			edge_dir_type scatter_edges(icontext_type& context,
					const vertex_type& vertex) const { 
				return graphlab::ALL_EDGES; 
			}; // end of scatter edges

			/** Scatter reschedules neighbors */  
			void scatter(icontext_type& context, const vertex_type& vertex, 
					edge_type& edge) const {
				edge_data& edata = edge.data();
				if(edata.role == edge_data::TRAIN) {
					const vertex_type other_vertex = get_other_vertex(edge, vertex);
					// Reschedule neighbors ------------------------------------------------
					if(other_vertex.data().nupdates < MAX_UPDATES) 
						context.signal(other_vertex, vec_type::Zero(vertex_data::NLATENT));
				}
			} // end of scatter function


			/**
			 * \brief Signal all vertices on one side of the bipartite graph
			 */
			static graphlab::empty signal_left(icontext_type& context,
					vertex_type& vertex) {
				if(vertex.num_out_edges() > 0) context.signal(vertex, vec_type::Zero(vertex_data::NLATENT));
				return graphlab::empty();
			} // end of signal_left 

	}; // end of sgd vertex program


struct error_aggregator : public graphlab::IS_POD_TYPE {
	typedef sgd_vertex_program::icontext_type icontext_type;
	typedef graph_type::edge_type edge_type;
	double train_error, validation_error;

	error_aggregator() : 
		train_error(0), validation_error(0){ }
	error_aggregator& operator+=(const error_aggregator& other) {
		train_error += other.train_error;
		assert(!std::isnan(train_error));
		validation_error += other.validation_error;
		return *this;
	}
	static error_aggregator map(icontext_type& context, const graph_type::edge_type& edge) {
		error_aggregator agg;
		if (edge.data().role == edge_data::TRAIN){
			if (isuser_node(edge.source())) 
				agg.train_error = extract_l2_error(edge); 
			assert(!std::isnan(agg.train_error));
		}
		else if (edge.data().role == edge_data::VALIDATE){
			if (isuser_node(edge.source())) 
				agg.validation_error = extract_l2_error(edge); 
		}
		return agg;
	}


	static void finalize(icontext_type& context, const error_aggregator& agg) {
		iter++;
		if (iter%2 == 0)
			return; 
		const double train_error = std::sqrt(agg.train_error / info.training_edges);
		assert(!std::isnan(train_error));
		context.cout() << std::setw(8) << context.elapsed_seconds()  << "  " << std::setw(8) << train_error;
		if(info.validation_edges > 0) {
			const double validation_error = 
				std::sqrt(agg.validation_error / info.validation_edges);
			context.cout() << "   " << std::setw(8) << validation_error; 
		}
		context.cout() << std::endl;
		sgd_vertex_program::GAMMA *= sgd_vertex_program::STEP_DEC;
	}
}; // end of error aggregator

/**
 * \brief Given an edge compute the error associated with that edge
 */
double extract_l2_error(const graph_type::edge_type & edge) {
	double pred = 
		edge.source().data().pvec.dot(edge.target().data().pvec);
	pred = std::min(sgd_vertex_program::MAXVAL, pred);
	pred = std::max(sgd_vertex_program::MINVAL, pred);
	double rmse = (edge.data().obs - pred) * (edge.data().obs - pred);
	assert(rmse <= pow(sgd_vertex_program::MAXVAL-sgd_vertex_program::MINVAL,2));
	return rmse;
} // end of extract_l2_error


struct prediction_saver {
	typedef graph_type::vertex_type vertex_type;
	typedef graph_type::edge_type   edge_type;
	/* save the linear model, using the format:
	   nodeid) factor1 factor2 ... factorNLATENT \n
	 */
	std::string save_vertex(const vertex_type& vertex) const {
		return "";
	}
	std::string save_edge(const edge_type& edge) const {
		if (edge.data().role != edge_data::PREDICT)
			return "";

		std::stringstream strm;
		const double prediction = 
			edge.source().data().pvec.dot(edge.target().data().pvec);
		strm << edge.source().id() << '\t' 
			<< -edge.target().id()-SAFE_NEG_OFFSET << '\t'
			<< prediction << '\n';
		return strm.str();
	}
}; // end of prediction_saver

struct linear_model_saver_U {
	typedef graph_type::vertex_type vertex_type;
	typedef graph_type::edge_type   edge_type;
	/* save the linear model, using the format:
	   nodeid) factor1 factor2 ... factorNLATENT \n
	 */
	std::string save_vertex(const vertex_type& vertex) const {
		if (vertex.num_out_edges() > 0){
			std::string ret = boost::lexical_cast<std::string>(vertex.id()) + ") ";
			for (uint i=0; i< vertex_data::NLATENT; i++)
				ret += boost::lexical_cast<std::string>(vertex.data().pvec[i]) + " ";
			ret += "\n";
			return ret;
		}
		else return "";
	}
	std::string save_edge(const edge_type& edge) const {
		return "";
	}
}; 

struct linear_model_saver_V {
	typedef graph_type::vertex_type vertex_type;
	typedef graph_type::edge_type   edge_type;
	/* save the linear model, using the format:
	   nodeid) factor1 factor2 ... factorNLATENT \n
	 */
	std::string save_vertex(const vertex_type& vertex) const {
		if (vertex.num_out_edges() == 0){
			std::string ret = boost::lexical_cast<std::string>(-vertex.id()-SAFE_NEG_OFFSET) + ") ";
			for (uint i=0; i< vertex_data::NLATENT; i++)
				ret += boost::lexical_cast<std::string>(vertex.data().pvec[i]) + " ";
			ret += "\n";
			return ret;
		}
		else return "";
	}
	std::string save_edge(const edge_type& edge) const {
		return "";
	}
}; 



/**
 * \brief The graph loader function is a line parser used for
 * distributed graph construction.
 */
inline bool graph_loader(graph_type& graph, 
		const std::string& filename,
		const std::string& line) {

	// Parse the line
	std::stringstream strm(line);
	graph_type::vertex_id_type source_id(-1), target_id(-1);
	float obs(0);
	strm >> source_id >> target_id;

	if (source_id == graph_type::vertex_id_type(-1) || target_id == graph_type::vertex_id_type(-1)){
		logstream(LOG_WARNING)<<"Failed to read input line: "<< line << " in file: "  << filename << " (or node id is -1). " << std::endl;
		return true;
	}

	// Determine the role of the data
	edge_data::data_role_type role = edge_data::TRAIN;
	if(boost::ends_with(filename,".validate")) role = edge_data::VALIDATE;
	else if(boost::ends_with(filename, ".predict")) role = edge_data::PREDICT;

	if(role == edge_data::TRAIN || role == edge_data::VALIDATE){
		strm >> obs;
		if (obs < sgd_vertex_program::MINVAL || obs > sgd_vertex_program::MAXVAL){
			logstream(LOG_WARNING)<<"Rating values should be between " << sgd_vertex_program::MINVAL << " and " << sgd_vertex_program::MAXVAL << ". Got value: " << obs << " [ user: " << source_id << " to item: " <<target_id << " ] " << std::endl; 
			assert(false); 
		}
	}
	target_id = -(graphlab::vertex_id_type(target_id + SAFE_NEG_OFFSET));
	// Create an edge and add it to the graph
	graph.add_edge(source_id, target_id, edge_data(obs, role)); 
	return true; // successful load
} // end of graph_loader








size_t vertex_data::NLATENT = 20;
double sgd_vertex_program::TOLERANCE = 1e-3;
double sgd_vertex_program::LAMBDA = 0.001;
double sgd_vertex_program::GAMMA = 0.001;
size_t sgd_vertex_program::MAX_UPDATES = -1;
double sgd_vertex_program::MAXVAL = 1e+100;
double sgd_vertex_program::MINVAL = -1e+100;
double sgd_vertex_program::STEP_DEC = 0.9;
bool sgd_vertex_program::debug = false;


/**
 * \brief The engine type used by the SGD matrix factorization
 * algorithm.
 *
 * The SGD matrix factorization algorithm currently uses the
 * synchronous engine.  However we plan to add support for alternative
 * engines in the future.
 */
typedef graphlab::omni_engine<sgd_vertex_program> engine_type;

int main(int argc, char** argv) {
	global_logger().set_log_level(LOG_INFO);
	global_logger().set_log_to_console(true);

	// Parse command line options -----------------------------------------------
	const std::string description = 
		"Compute the SGD factorization of a matrix.";
	graphlab::command_line_options clopts(description);
	std::string input_dir;
	std::string predictions;
	size_t interval = 0;
	std::string exec_type = "synchronous";
	clopts.attach_option("matrix", input_dir,
			"The directory containing the matrix file");
	clopts.add_positional("matrix");
	clopts.attach_option("D", vertex_data::NLATENT,
			"Number of latent parameters to use.");
	clopts.attach_option("engine", exec_type, 
			"The engine type synchronous or asynchronous");
	clopts.attach_option("max_iter", sgd_vertex_program::MAX_UPDATES,
			"The maxumum number of udpates allowed for a vertex");
	clopts.attach_option("lambda", sgd_vertex_program::LAMBDA, 
			"SGD regularization weight"); 
	clopts.attach_option("gamma", sgd_vertex_program::GAMMA, 
			"SGD step size"); 
	clopts.attach_option("debug", sgd_vertex_program::debug, 
			"debug - additional verbose info"); 
	clopts.attach_option("tol", sgd_vertex_program::TOLERANCE,
			"residual termination threshold");
	clopts.attach_option("maxval", sgd_vertex_program::MAXVAL, "max allowed value");
	clopts.attach_option("minval", sgd_vertex_program::MINVAL, "min allowed value");
	clopts.attach_option("step_dec", sgd_vertex_program::STEP_DEC, "multiplicative step decrement");
	clopts.attach_option("interval", interval, 
			"The time in seconds between error reports");
	clopts.attach_option("predictions", predictions,
			"The prefix (folder and filename) to save predictions.");

	parse_implicit_command_line(clopts);

	if(!clopts.parse(argc, argv) || input_dir == "") {
		std::cout << "Error in parsing command line arguments." << std::endl;
		clopts.print_description();
		return EXIT_FAILURE;
	}
	debug = sgd_vertex_program::debug;
	//  omp_set_num_threads(clopts.get_ncpus());
	///! Initialize control plain using mpi
	graphlab::mpi_tools::init(argc, argv);
	graphlab::distributed_control dc;

	dc.cout() << "Loading graph." << std::endl;
	graphlab::timer timer; 
	graph_type graph(dc, clopts);  
	graph.load(input_dir, graph_loader); 
	dc.cout() << "Loading graph. Finished in " 
		<< timer.current_time() << std::endl;

	if (dc.procid() == 0) 
		add_implicit_edges<edge_data>(implicitratingtype, graph, dc);


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

	dc.cout() << "Creating engine" << std::endl;
	engine_type engine(dc, graph, exec_type, clopts);

	// Add error reporting to the engine
	const bool success = engine.add_edge_aggregator<error_aggregator>
		("error", error_aggregator::map, error_aggregator::finalize) &&
		engine.aggregate_periodic("error", interval);
	ASSERT_TRUE(success);


	// Signal all vertices on the vertices on the left (libersgd) 
	engine.map_reduce_vertices<graphlab::empty>(sgd_vertex_program::signal_left);
	info = graph.map_reduce_edges<stats_info>(count_edges);
	dc.cout()<<"Training edges: " << info.training_edges << " validation edges: " << info.validation_edges << std::endl;


	// Run the PageRank ---------------------------------------------------------
	dc.cout() << "Running SGD" << std::endl;
	dc.cout() << "(C) Code by Danny Bickson, CMU " << std::endl;
	dc.cout() << "Please send bug reports to danny.bickson@gmail.com" << std::endl;
	dc.cout() << "Time   Training    Validation" <<std::endl;
	dc.cout() << "       RMSE        RMSE " <<std::endl;
	timer.start();
	engine.start();  

	const double runtime = timer.current_time();
	dc.cout() << "----------------------------------------------------------"
		<< std::endl
		<< "Final Runtime (seconds):   " << runtime 
						    << std::endl
								<< "Updates executed: " << engine.num_updates() << std::endl
											      << "Update Rate (updates/second): " 
												      << engine.num_updates() / runtime << std::endl;

	// Compute the final training error -----------------------------------------
	dc.cout() << "Final error: " << std::endl;
	engine.aggregate_now("error");

	// Make predictions ---------------------------------------------------------
	if(!predictions.empty()) {
		std::cout << "Saving predictions" << std::endl;
		const bool gzip_output = false;
		const bool save_vertices = false;
		const bool save_edges = true;
		const size_t threads_per_machine = 1;
		//save the predictions
		graph.save(predictions, prediction_saver(),
				gzip_output, save_vertices, 
				save_edges, threads_per_machine);
		//save the linear model
		graph.save(predictions + ".U", linear_model_saver_U(),
				gzip_output, save_edges, save_vertices, threads_per_machine);
		graph.save(predictions + ".V", linear_model_saver_V(),
				gzip_output, save_edges, save_vertices, threads_per_machine);

	}



	graphlab::mpi_tools::finalize();
	return EXIT_SUCCESS;
} // end of main



