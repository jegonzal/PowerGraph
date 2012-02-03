#include "linear.h"
#include "../../libs/matrixmarket/mmio.h"
#include "advanced_config.h"


template<typename vertex_data>
struct node_cache{
   int from;
   int to;
   vertex_data *pfrom;
   vertex_data *pto;

   node_cache(){
      from = to = -1;
      pfrom = pto = NULL;
   }

};
node_cache<vertex_data_shotgun> vnode_cache;

extern problem_setup ps;
extern advanced_config config;
extern const char * runmodesnames[];

void save_matrix_market_format();
void write_output_itpp();
void write_vec(FILE * f, int len, double * array);
FILE * load_matrix_metadata(const char * filename);
template <typename graph_type, typename vertex_data, typename edge_data>
void load_square_matrix(FILE * f, graph_type& graph, advanced_config & config);
template<typename graph_type, typename vertex_data, typename edge_data>
void load_non_square_matrix(FILE * f, graph_type& graph, advanced_config & config);


extern graphlab::glshared<double> RELATIVE_NORM_KEY;

FILE * open_file(const char * name, const char * mode){
  FILE * f = fopen(name, mode);
  if (f == NULL){
      perror("fopen failed");
      logstream(LOG_FATAL) <<" Failed to open file" << name << std::endl;
   }
  return f;
}


template<typename graph_type, typename vertex_data>
void add_vertices(graph_type& graph){
     for (int i=0; i< (int)(config.square ? ps.m : ps.m+ps.n); i++){
       vertex_data data;
       init_vertex_data(data, config.supportgraphlabcf ? 1 : 0, GABP_PRIOR_MEAN_OFFSET, ps.last_node, i); 
       graph.add_vertex(data);
     }
}



void fill_output(graph_type * g){

  double diff = 0;
  int start = ps.m;
  if (config.square) 
    start = 0;
  for (size_t i = start; i < g->num_vertices(); i++){
    const vertex_data& vdata = g->vertex_data(i);
     if (config.algorithm == JACOBI || config.algorithm == GaBP){
       diff += pow(vdata.real - vdata.cur_mean,2);
       ps.means.push_back(vdata.cur_mean);
       ps.prec.push_back(vdata.cur_prec);
     }
     //TODO: else if (algorithm == CONJUGATE_GRADIENT)
     //   diff += pow(means[i-m] - vdata.real,2);
  }
  assert(ps.means.size() > 0);
   if (config.algorithm == JACOBI || config.algorithm == GaBP){
  std::cout << "Assuming the linear system is Ax=y, and the correct solution is x*," << runmodesnames[config.algorithm] << " converged to an accuracy norm(x-x*) of " << diff
            << " msg norm is: " << RELATIVE_NORM_KEY.get_val() << std::endl;
   }
}
void fill_output(graph_type_inv * g){

  for (size_t i = 0; i < g->num_vertices(); i++){
    const vertex_data_inv& vdata = g->vertex_data(i);
       for (size_t j=0; j< g->num_vertices(); j++){
	 ps.means.push_back(vdata.cur_mean[j]);
       }  
       ps.prec.push_back(vdata.cur_prec);
  }
}

void fill_output(graph_type_shotgun * g){
  for (size_t i = ps.m; i < g->num_vertices(); i++){
    const vertex_data_shotgun& vdata = g->vertex_data(i);
       ps.means.push_back(vdata.x);
    }
  
 
}

void write_output(){

   if (config.matrixmarket)
      save_matrix_market_format();
   else write_output_itpp();
}

void write_output_itpp(){

   FILE * f = fopen((config.datafile+".out").c_str(), "w");
   assert(f!= NULL);

   std::cout<<"Writing result to file: "<<config.datafile<<".out"<<std::endl;
   std::cout<<"You can read the file in Matlab using the load_c_gl.m matlab script"<<std::endl;
   write_vec(f, ps.means.size(), &ps.means[0]);
   if (config.algorithm == GaBP)
     write_vec(f, ps.prec.size(), &ps.prec[0]);

   fclose(f);




}


/* handle reading diagonal entries of square matrix, which are usually not edges! */
void init_prior_prec(vertex_data & data, double val){
   data.prior_prec = val;
}
void init_prior_prec(vertex_data_inv & data, double val){
   data.prior_prec = val;
}
void init_prior_prec(vertex_data_shotgun & data, double val){
  //TODOw
}



void init_vertex_data(vertex_data_inv & data, double val, int offset, int nodes, int i){
      assert(config.square);
      data.prior_mean = new sdouble[nodes];
      memset(data.prior_mean, 0, nodes*sizeof(sdouble));
      data.prior_mean[i] = 1;
      data.cur_mean = new sdouble[nodes];
      memset(data.cur_mean, 0, nodes*sizeof(sdouble));
      data.prev_mean = new sdouble[nodes];
      memset(data.prev_mean, 1, nodes*sizeof(sdouble));
}

void init_vertex_data(vertex_data_shotgun & data, double val, int offset, int nodes, int i){
     if (i < (int)ps.m){
        if (config.supportgraphlabcf)
           data.y = (drand48() < 0.5 ? -1: 1);
        else data.y = val;
     } 
     data.x = 0;
     data.expAx = 1;
     set_size(data.features, (i < (int)ps.m ? ps.n: ps.m ));
}

void init_vertex_data(vertex_data & data, double val, int offset, int nodes, int i){
   if (offset == GABP_PRIOR_MEAN_OFFSET)
          data.prior_mean = val;
   else if (offset == GABP_REAL_OFFSET)
          data.real = val;
   else assert(false);
}

template <typename graph_type, typename vertex_data>
void load_matrix_market_vector(graph_type * g)
{
    int ret_code;
    MM_typecode matcode;
    int M, N, nz;   
    int i;

    FILE * f = open_file((config.datafile + "v").c_str(), "r");
    if (mm_read_banner(f, &matcode) != 0)
        logstream(LOG_FATAL) << "Could not process Matrix Market banner in file " << config.datafile << "v " << std::endl;

    /*  This is how one can screen matrix types if their application */
    /*  only supports a subset of the Matrix Market data types.      */

    if (mm_is_complex(matcode) && mm_is_matrix(matcode) && 
            mm_is_sparse(matcode) )
        logstream(LOG_FATAL) << "sorry, this application does not support " << std::endl << 
          "Market Market type: " << mm_typecode_to_str(matcode) << std::endl;

    /* find out size of sparse matrix .... */
    if (mm_is_sparse(matcode)){
      if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0){
       logstream(LOG_FATAL) << "failed to read matrix market cardinality size " << std::endl; 
      } 
    } 
    else {
      if ((ret_code = mm_read_mtx_array_size(f, &M, &N)) != 0){
       logstream(LOG_FATAL) << "failed to read matrix market cardinality size " << std::endl; 
      }
      nz = M*N;
    }
    
    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

    int I,J; 
    double val;


    for (i=0; i<nz; i++)
    {
        if (mm_is_sparse(matcode)){
          int rc = fscanf(f, "%d %d %lg\n", &I, &J, &val);
          if (rc != 3)
            logstream(LOG_FATAL) << "Failed to read line: " << i << " in file: " << config.datafile << std::endl;
          I--;  /* adjust from 1-based to 0-based */
          J--;
        }
        else {
          int rc = fscanf(f, "%lg", &val);
          if (rc != 1)
            logstream(LOG_FATAL) << "Failed to read line: " << i << " in file: " << config.datafile << std::endl;
           J = 0;
          I = i;
        }
         if (config.scalerating != 1.0)
	     val /= config.scalerating;
         if (!config.zero)
	   assert(val!=0 );
        
        assert(I< M);
        assert(J< N);
        vertex_data & vdata = g->vertex_data(I);
        int offset = GABP_PRIOR_MEAN_OFFSET;
        if (!config.square && I >= (int)ps.m)
	  offset = GABP_REAL_OFFSET;
        init_vertex_data(vdata, val, offset, config.square ? ps.n : ps.m+ps.n, I);
    }
    fclose(f);

}

template<typename graph_type, typename vertex_data, typename edge_data>
void load_matrix_market_matrix(graph_type * g)
{
    int ret_code;
    MM_typecode matcode;
    int M, N, nz;   
    int i;

    FILE * f = open_file(config.datafile.c_str(), "r");
    if (mm_read_banner(f, &matcode) != 0)
        logstream(LOG_FATAL) << "Could not process Matrix Market banner: " << config.datafile << std::endl;

    /*  This is how one can screen matrix types if their application */
    /*  only supports a subset of the Matrix Market data types.      */

    if (mm_is_complex(matcode) && mm_is_matrix(matcode) && 
            mm_is_sparse(matcode) )
        logstream(LOG_FATAL) << "sorry, this application does not support " << std::endl << 
          "Market Market type: " << mm_typecode_to_str(matcode) << std::endl;

    /* find out size of sparse matrix .... */
    if (mm_is_sparse(matcode)){
      if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)
        logstream(LOG_FATAL) << "failed to read matrix market cardinality size " << std::endl; 
    } else {
      if ((ret_code = mm_read_mtx_array_size(f, &M, &N)) != 0){
        logstream(LOG_FATAL) << "failed to read matrix market cardinality size " << std::endl; 
      }
      nz = M*N;
    }
    ps.m = M; ps.n = N; 
    if (ps.m == ps.n){
       config.square = true;
       ps.last_node = ps.n;
    }
    else {
       config.square = false;
       ps.last_node = ps.n+ps.m;
    }

    add_vertices<graph_type,vertex_data>(*g);

    int I,J; 
    double val;

    int offset = 0;
    if (!config.square)
	offset = ps.m;

    for (i=0; i<nz; i++)
    {

        if (mm_is_sparse(matcode)){
          int rc= fscanf(f, "%d %d %lg\n", &I, &J, &val);
          if (rc != 3)
            logstream(LOG_FATAL) << "Failed to read line: " << i << " in file: " << config.datafile << std::endl;
           I--;  /* adjust from 1-based to 0-based */
           J--;
         } else {
	  int rc = fscanf(f, "%lg", &val);
          if (rc != 1)
            logstream(LOG_FATAL) << "Failed to read line: " << i << " in file: " << config.datafile << std::endl;
          I = i / N;
          J = i % N;  
        }
        if (config.scalerating != 1.0)
	     val /= config.scalerating;
         if (!config.zero)
	   assert(val!=0 );
        
        assert(I< M);
        assert(J< N);
     
        if (config.square && I == J){
	   vertex_data & data = g->vertex_data(I);
           init_prior_prec(data, val);
        }
	else {
          add_edge(val, I, J, g);
          if (mm_is_symmetric(matcode))
            add_edge(val, J, I, g);
        }
        
   }
    ps.e = nz;
    fclose(f);

}


void save_matrix_market_vector(const std::string& filename, const std::vector<double> & a, std::string comment, bool integer,bool issparse){
    MM_typecode matcode;                        
    int i;

    mm_initialize_typecode(&matcode);
    mm_set_matrix(&matcode);
    if (issparse){
      mm_set_sparse(&matcode);
      mm_set_coordinate(&matcode);
    }
    else {
      mm_set_dense(&matcode);
      mm_set_array(&matcode);
    }
    if (!integer)
      mm_set_real(&matcode);
    else
      mm_set_integer(&matcode);

    FILE * f = open_file(filename.c_str(),"w");
    mm_write_banner(f, matcode); 
    if (comment.size() > 0)
      fprintf(f, "%s%s", "%", comment.c_str());
    
    if (issparse)
      mm_write_mtx_crd_size(f, a.size(), 1, a.size());
    else 
      mm_write_mtx_array_size(f, a.size(), 1);

    for (i=0; i<a.size(); i++){
      if (issparse){
        if (integer)
          fprintf(f, "%d %d %d\n", i+1, 1, (int)a[i]);
        else fprintf(f, "%d %d %10.13g\n", i+1, 1, a[i]);
      }
      else {//dense
        if (integer)
          fprintf(f,"%d ", (int)a[i]);
        else fprintf(f, "%10.13g\n", a[i]);
      }
    }

    logstream(LOG_INFO) << "Saved output vector to file: " << filename << std::endl;
    logstream(LOG_INFO) << "You can read it with Matlab/Octave using the script mmread.m found on http://graphlab.org/mmread.m" << std::endl;
}


void save_matrix_market_format()
{
    if (ps.means.size() > 0){
      save_matrix_market_vector(config.datafile + ".x.mtx", 
 			        ps.means, ps.x_comment, 
                                false, false);
    }

    if (ps.prec.size() > 0){
      save_matrix_market_vector(config.datafile + ".prec.mtx", 
			        ps.prec, ps.prec_comment, 
	 			false, false);
    }

}


template<typename graph_type, typename vertex_data, typename edge_data>
void load_data(graph_type * g){

  if (config.matrixmarket){
     load_matrix_market_matrix<graph_type, vertex_data, edge_data>(g);
     load_matrix_market_vector<graph_type, vertex_data>(g);
     return; 
 }

  // Load the graph --------------------------------------------------
  FILE * f = load_matrix_metadata(config.datafile.c_str());
  if (ps.m == ps.n){ //square matrix
     config.square = true;
     ps.last_node = ps.n;
     load_square_matrix<graph_type, vertex_data, edge_data>(f, *g, config);
  }
  else {
     config.square = false;
     ps.last_node = ps.m+ps.n;
     if (config.algorithm == JACOBI || config.algorithm == GaBP_INV){
        logstream(LOG_ERROR)<<" Jacobi/GaBP-INV can not run with non-square mastrix." << std::endl;
        exit(1);
     }

     load_non_square_matrix<graph_type, vertex_data, edge_data>(f, *g, config);
  }
}


void load_data_gamp(graph_type_gamp *g);



template<>
void load_data<graph_type_gamp, vertex_data_gamp, edge_data_gamp>(graph_type_gamp* g){
    load_data_gamp(g);
 }


#define BUFSIZE 500000
template <typename graph_type, typename vertex_data>
void read_nodes(FILE * f, int len, int offset, int nodes,
                graph_type * g){

  assert(offset>=0 && offset < len);
  assert(nodes>0);

  int toread = nodes;
  if (nodes > BUFSIZE)
    toread = BUFSIZE;
  int remain = nodes;

  while(remain > 0){
    double * temp = new double[remain];
    toread = (remain < BUFSIZE)?remain:toread;
    int rc = (int)fread(temp, sizeof(double), toread, f);
    //assert(rc == toread);
    remain -= rc;
    for (int i=0; i< rc; i++){
      vertex_data data;
      init_vertex_data(data, temp[i], offset, nodes, i);     
      g->add_vertex(data);
    }
    delete [] temp;
  }
}



//read a vector from file and return an array
double * read_vec(FILE * f, size_t len){
  double * vec = new double[len];
  assert(vec != NULL);
  fread(vec, len, sizeof(double), f);
  return vec;
}
//write an output vector to file
void write_vec(FILE * f, int len, double * array){
  assert(f != NULL && array != NULL);
  fwrite(array, len, sizeof(double), f);
}


template<typename graph>
void dispatch_vec(int start_pos, int end_pos, int offset,
                  graph * g, double * vec, int len,
                  bool free_vec = false, int d = 1){
  assert(end_pos - start_pos > 1);
  assert(len == end_pos - start_pos);
  assert(vec != NULL);
  assert(offset>=0 && offset < len);
  for (int i=start_pos; i < end_pos; i++){
    sdouble * data = (sdouble*)&g->vertex_data(i);
    for (int j=0; j<d; j++){
      data[offset+d*j] = vec[d*(i-start_pos)+j];
    }
  }
  if (free_vec)
    delete [] vec;
}

template<typename graph>
void dispatch_vec(int start_pos, int end_pos, int offset,
                  graph * g,sdouble val,int d = 1){
  assert(end_pos - start_pos > 1);
  for (int i=start_pos; i < end_pos; i++){
    sdouble * data = (sdouble*)&g->vertex_data(i);
    for (int j=0; j<d; j++)
      data[offset+d*j] = val;
  }
}



//struct for holding edge data in file
struct edata{
  float from;
  float to;
  float weight;
};

struct graphlab_cf_format{
  float from;
  float to;
  float time;
  float weight;
};
struct graphlab_cf_old_format{
  int from;
  int to;
  double time;
  double weight;
};


struct graphlab_old_format{
  int from;
  int to; 
  double weight;
};

void add_edge(double val, int from, int to, graph_type *g){
   edge_data edge;
   edge.weight = val;
   if (!config.zero)
      assert(val != 0);
     
   assert(from >= 0 && to >= 0);
   assert(from < (int)ps.m && to < ps.last_node);
   assert(from != to);
    
    // Matlab export has ids starting from 1, ours start from 0
   g->add_edge(from, to, edge);
   if (!config.square) { //add the reverse edge as well
        // Matlab export has ids starting from 1, ours start from 0
        g->add_edge(to, from, edge);
   }
}
void add_edge(double val, int from, int to, graph_type_inv *g){
      edge_data_inv edge;
      edge.mean = new sdouble[ps.m];
      memset(edge.mean, 0, sizeof(sdouble)*ps.m);
      edge.weight =  val;
      if (!config.zero)
        assert(val != 0);
      assert(from >= 0 && to >= 0);
      assert(from < (int)ps.m && to < ps.last_node);
      assert(from != to);
  
      // Matlab export has ids starting from 1, ours start from 0
      g->add_edge(from, to, edge);
  }

void add_edge(double val, int from, int to, graph_type_shotgun * g){

      assert(from >= 0 && from < (int)ps.m);
      assert(to >= (int)ps.m && to < (int)ps.last_node);
      assert(from != to);
      if (vnode_cache.from == from){
        set_new(vnode_cache.pfrom->features, to-ps.m, val); 
      }
      else {
       vnode_cache.pfrom = &g->vertex_data(from);
       vnode_cache.from = from;  
       set_new(vnode_cache.pfrom->features, to-ps.m, val); 
      }
      if (vnode_cache.to == to){
        set_new(vnode_cache.pto->features, from, val);
      }
      else {
        vnode_cache.to = to;
        vnode_cache.pto = &g->vertex_data(to);   
        set_new(vnode_cache.pto->features, from, val); 
      }
}

//read edges from file into the graph
template<typename edatatype, typename graph_type, typename vertex_data, typename edge_data>
int read_edges(FILE * f, int nodes, graph_type * g, advanced_config & config){
  assert(nodes > 0);

  unsigned int e,g0;
  int rc = fread(&e,1,4,f); //read the number of edges
  assert(rc == 4);
  if (!config.supportgraphlabcf){
    rc = fread(&g0,1,4,f); //zero pad
    assert(rc == 4);
    assert(g0 == 0);
  }
  printf("Creating %d edges...\n", e);
  assert(e>0);
  int total = 0;
  edatatype* ed = new edatatype[200000];
  int edgecount_in_file = e;

  //add unioque edges
  while(true){
    memset(ed, 0, 200000*sizeof(edatatype));
    rc = (int)fread(ed, sizeof(edatatype),
                    std::min(200000, edgecount_in_file - total), f);
    total += rc;
   
    for (int i=0; i<rc; i++){
      add_edge(ed[i].weight, ed[i].from-1, ed[i].to-1, g);
    }
    printf(".");
    fflush(0);
    if (rc == 0 || total >= edgecount_in_file)
      break;
  }
  delete [] ed; ed = NULL;
  return e;
}


FILE * load_matrix_metadata(const char * filename){
   printf("Loading %s\n", filename);
   FILE * f = fopen(filename, "r");
   assert(f!= NULL);

   fread(&ps.m, 1, 4, f);
   fread(&ps.n, 1, 4, f);
   if (ps.n == 0) 
      ps.n=ps.m; //compatability with older file format, will be removed later
   return f;
}

/*
 *  READ A SQUARE MATRIX A of size nxn
 *  For Gapbe: the main digonal is the precision vector
 * */

template <typename graph_type, typename vertex_data, typename edge_data>
void load_square_matrix(FILE * f, graph_type& graph, advanced_config & config) {

  assert(ps.m == ps.n);
  assert(ps.n > 0);
  //read the observation vector y of size n
  read_nodes<graph_type,vertex_data>(f, sizeof(vertex_data)/sizeof(sdouble), GABP_PRIOR_MEAN_OFFSET,
             ps.n, &graph);

  //read the real solution of size n (if given, otherwise it is zero)
  double * real = read_vec(f, ps.n);
  if (config.algorithm == GaBP || config.algorithm == JACOBI || config.algorithm == CONJUGATE_GRADIENT)
    dispatch_vec(0,ps.n,GABP_REAL_OFFSET, &graph, real, ps.n, true);
  else delete [] real;

  //read the precition of size n (the main diagonal of the matrix A)
  double * prec = read_vec(f, ps.n);
  dispatch_vec(0,ps.n,GABP_PRIOR_PREC_OFFSET, &graph, prec, ps.n, true);

  if (config.oldformat)
    ps.e = read_edges<graphlab_old_format,graph_type, vertex_data, edge_data>(f, ps.n, &graph, config);
  else
    ps.e = read_edges<edata,graph_type, vertex_data, edge_data>(f, ps.n, &graph, config);
  fclose(f);
}



/**
 * READ A NON-SQUARE matrix of size m rows x n cols 
 * Where the observation y is a vector of size m, the solution vector x=A\y is 
 * a vector of size n.
 */
template<typename graph_type, typename vertex_data, typename edge_data>
void load_non_square_matrix(FILE * f, graph_type& graph, advanced_config & config) {
  
  assert( ps.n > 0);
  assert(!config.square);  

  printf("Loading a non-square matrix A of size %d x %d\n", ps.m,ps.n);

  if (config.supportgraphlabcf){ //read matrix factorization file (GraphLab collabrative filtering format)
     int tmp;
     fread(&tmp, 1, 4, f); //skip over time bin number 
     add_vertices<graph_type, vertex_data>(graph);

     if (config.oldformat)
       ps.e = read_edges<graphlab_cf_old_format, graph_type, vertex_data, edge_data>(f, ps.n+ps.m, &graph, config);
     else
       ps.e = read_edges<graphlab_cf_format, graph_type, vertex_data, edge_data>(f, ps.n+ps.m, &graph, config);
  }
  else { //read A, x, y, from file
  
    //read y (the observation) of size m
    read_nodes<graph_type, vertex_data>(f, sizeof(vertex_data)/sizeof(sdouble),
             GABP_PRIOR_MEAN_OFFSET,ps.m,&graph);
    
    //read x (the real solution, if given) of size n
    read_nodes<graph_type, vertex_data>(f, sizeof(vertex_data)/sizeof(sdouble),
             GABP_REAL_OFFSET,ps.n,&graph);
    
    assert((int)graph.num_vertices() == ps.last_node);


    //read the precision vector of size m+n (of both solution and observation)
    double * prec = read_vec(f, ps.n+ps.m);
    if (config.algorithm == JACOBI || config.algorithm == GaBP || config.algorithm == CONJUGATE_GRADIENT || config.algorithm == GaBP_INV){
      dispatch_vec(0,ps.n+ps.m,GABP_PRIOR_PREC_OFFSET, &graph, prec, ps.n+ps.m, true);
      dispatch_vec(0,ps.n+ps.m,GABP_PREV_MEAN_OFFSET, &graph, 1);
      dispatch_vec(0,ps.n+ps.m,GABP_PREV_PREC_OFFSET, &graph, 1);
    }

    if (config.oldformat)
       ps.e = read_edges<graphlab_old_format,graph_type, vertex_data, edge_data>(f, ps.n+ps.m, &graph,config);
    else
       ps.e = read_edges<edata,graph_type, vertex_data, edge_data>(f, ps.n+ps.m, &graph,config);
  } 
  fclose(f);
}

