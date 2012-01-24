/* 
*   Matrix Market I/O example program
*
*   Read a real (non-complex) sparse matrix from a Matrix Market (v. 2.0) file.
*   and copies it to stdout.  This porgram does nothing useful, but
*   illustrates common usage of the Matrix Matrix I/O routines.
*   (See http://math.nist.gov/MatrixMarket for details.)
*
*   Usage:  a.out [filename] > output
*
*       
*   NOTES:
*
*   1) Matrix Market files are always 1-based, i.e. the index of the first
*      element of a matrix is (1,1), not (0,0) as in C.  ADJUST THESE
*      OFFSETS ACCORDINGLY offsets accordingly when reading and writing 
*      to files.
*
*   2) ANSI C requires one to use the "l" format modifier when reading
*      double precision floating point numbers in scanf() and
*      its variants.  For example, use "%lf", "%lg", or "%le"
*      when reading doubles, otherwise errors will occur.
*/

#include <stdio.h>
#include <stdlib.h>
#include "clustering.h"
#include "../../libs/matrixmarket/mmio.h"
#include "../gabp/advanced_config.h"
#include <assert.h>

extern advanced_config ac;
extern problem_setup ps;

void init();
void compact(graph_type *g);
FILE * open_file(const char * filename, const char * mode);

void load_matrix_market(const char * filename, graph_type_kcores *_g, testtype type)
{
    int ret_code;
    MM_typecode matcode;
    int M, N, nz;   
    int i;

    FILE * f = fopen(filename, "r");
    if (f== NULL){
        if (type ==VALIDATION && ac.algorithm != ITEM_KNN && ac.algorithm != USER_KNN)
	   return;
        if (type == TEST)
           return;
	logstream(LOG_ERROR) << " can not find input file: " << filename << ". aborting " << std::endl;
	exit(1);
    }

    if (mm_read_banner(f, &matcode) != 0)
    {
        logstream(LOG_ERROR) << "Could not process Matrix Market banner." << std::endl;
        exit(1);
    }

    /*  This is how one can screen matrix types if their application */
    /*  only supports a subset of the Matrix Market data types.      */

    if (mm_is_complex(matcode) && mm_is_matrix(matcode) && 
            mm_is_sparse(matcode) )
    {
        logstream(LOG_ERROR) << "sorry, this application does not support " << std::endl << 
          "Market Market type: " << mm_typecode_to_str(matcode) << std::endl;
        exit(1);
    }

    /* find out size of sparse matrix .... */

    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0){
       logstream(LOG_ERROR) << "failed to read matrix market cardinality size " << std::endl; 
       exit(1);
    }

    if (type ==TRAINING){
    ps.M = M; ps.N = N; ps.K = ac.K;
    }
    else if (type ==VALIDATION){
      ps.M_validation = M; ps.N_validation = N;
    }
    else if (type == TEST){
      ps.M_test = M; ps.N_test = N;
    }
    init();
    add_vertices(_g, type); 

    /* reseve memory for matrices */

    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

    int I,J; 
    double val;
    int step = nz/10;
    if (step == 0)
      step = nz - 1;

    for (i=0; i<nz; i++)
    {
        int rc = fscanf(f, "%d %d %lg\n", &I, &J, &val);
        if (rc != 3)
          logstream(LOG_FATAL) << "Failed to read matrix market file " << filename << " line: " << i << std::endl;
  
        //printf("Found row %d %d %g\n", I, J, val);        
        I--;  /* adjust from 1-based to 0-based */
        J--;
        if (I<=0 || J<= 0){
          logstream(LOG_ERROR) << "Matrix market values should be >= 1, observed values: " << I << " " << J << " In item number " << nz << std::endl;
          exit(1);
        }
 if (ac.scalerating != 1.0)
	     val /= ac.scalerating;
         if (!ac.zero)
	   assert(val!=0 );
        
        assert(I >= 0 && I< M);
        assert(J >=0 && J< N);
        kcores_edge edge;
        edge.weight = val;
        _g->add_edge(I,J, edge);
        _g->add_edge(J,I, edge);
        if (i % step == 0)
          logstream(LOG_INFO) << "Matrix market read " << i << " edges" << std::endl;
    }
    ps.L = nz;
    fclose(f);
}

void load_matrix_market_clusters(const std::string & filename, graph_type_kcores *_g){
  assert(false);
}
void load_matrix_market_clusters(const std::string & filename, graph_type *_g)
{
    assert(ac.algorithm == K_MEANS);
    int ret_code;
    MM_typecode matcode;
    int M, N, nz;   
    int i;

    FILE * f = open_file(filename.c_str() , "r");

    if (mm_read_banner(f, &matcode) != 0)
        logstream(LOG_FATAL) << "Could not process Matrix Market banner." << std::endl;

    if (mm_is_complex(matcode) && mm_is_matrix(matcode) && 
            mm_is_sparse(matcode) )
        logstream(LOG_FATAL) << "sorry, this application does not support " << std::endl << 
          "Market Market type: " << mm_typecode_to_str(matcode) << std::endl;

    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)
       logstream(LOG_FATAL) << "failed to read matrix market cardinality size " << std::endl; 


    int I,J; 
    double val;
    int step=nz/10;
  
    if (mm_is_sparse(matcode)){
    for (i=0; i<nz; i++)
    {
        int rc = fscanf(f, "%d %d %lg\n", &I, &J, &val);
        if (rc != 3)
          logstream(LOG_FATAL) << "Failed to read matrix market file " << filename << " line: " << i << std::endl;
        //printf("%d) Found row %d %d %g\n", i, I, J, val);        
        I--;  /* adjust from 1-based to 0-based */
        J--;
        
        assert(I>=0 && I< ac.K);
        assert(J >=0 && J< N);
        set_val(ps.clusts.cluster_vec[I].location, J, val);

        if (i % step == 0)
          logstream(LOG_INFO) << "Matrix market read " << i << " edges" << std::endl;
     }
     }
     else {
	for (I=0; I<M; I++){
           for (J=0; J< N; J++){
              int rc = fscanf(f, "%lg", &val);
              if (rc != 1)
                 logstream(LOG_FATAL) << " Failed to read matrix market file: " << filename << " on line: " << I << " column: " << J << std::endl;
              set_val(ps.clusts.cluster_vec[I].location, J, val);
           }
        }
      }
    
    
    fclose(f);
}
void load_matrix_market_assignments(const std::string & filename, graph_type_kcores *_g){
   assert(false);
}
void load_matrix_market_assignments(const std::string & filename, graph_type *_g)
{
    assert(ac.algorithm == K_MEANS);
    int ret_code;
    MM_typecode matcode;
    int M, N, nz;   
    int i;

    FILE * f = open_file(filename.c_str() , "r");

    if (mm_read_banner(f, &matcode) != 0)
        logstream(LOG_FATAL) << "Could not process Matrix Market banner." << std::endl;

    if (mm_is_complex(matcode) && mm_is_matrix(matcode) && 
            mm_is_sparse(matcode) )
        logstream(LOG_FATAL) << "sorry, this application does not support " << std::endl << 
          "Market Market type: " << mm_typecode_to_str(matcode) << std::endl;

    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)
       logstream(LOG_FATAL) << "failed to read matrix market cardinality size " << std::endl; 


    int I,J; 
    double val;
    int step=nz/10;
    if (step == 0)
      step = nz - 1;
    ps.M = nz;
    assert(ps.N > 0);

    if (_g->num_vertices() == 0){
      init();
      add_vertices(_g, TRAINING); 
    }


    for (i=0; i<nz; i++)
    {
	if (mm_is_sparse(matcode)){
          int rc = fscanf(f, "%d %d %lg\n", &I, &J, &val);
          if (rc != 3)
            logstream(LOG_FATAL) << "Failed to read matrix market file " << filename << " line: " << i << std::endl;
	  I--;  /* adjust from 1-based to 0-based */
          J--;
        
          assert(I>=0 && I< nz);
          assert(J == 0);
          assert(val >= 0 && val < ac.K);
	  _g->vertex_data(i).current_cluster = val;
        }
        else {
          int rc = fscanf(f, "%lg\n", &val);
          if (rc != 1)
            logstream(LOG_FATAL) << "Failed to read matrix market file: " << filename << " line: " << i << std::endl;
          _g->vertex_data(i).current_cluster = val; 
        }
        if (i % step == 0)
          logstream(LOG_INFO) << "Matrix market read " << i << " edges" << std::endl;
     }
    
    fclose(f);
}





void load_matrix_market(const char * filename, graph_type *_g, testtype type)
{
    int ret_code;
    MM_typecode matcode;
    int M, N, nz;   
    int i;

    FILE * f = fopen(filename, "r");
    if (f== NULL){
        if (type == VALIDATION && ac.algorithm != ITEM_KNN && ac.algorithm != USER_KNN)
           return;
        if (type == TEST)
           return;
	logstream(LOG_ERROR) << " can not find input file. aborting " << std::endl;
	exit(1);
    }

    if (mm_read_banner(f, &matcode) != 0)
    {
        logstream(LOG_ERROR) << "Could not process Matrix Market banner." << std::endl;
        exit(1);
    }

    /*  This is how one can screen matrix types if their application */
    /*  only supports a subset of the Matrix Market data types.      */

    if (mm_is_complex(matcode) && mm_is_matrix(matcode) && 
            mm_is_sparse(matcode) )
    {
        logstream(LOG_ERROR) << "sorry, this application does not support " << std::endl << 
          "Market Market type: " << mm_typecode_to_str(matcode) << std::endl;
        exit(1);
    }

    /* find out size of sparse matrix .... */

    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0){
       logstream(LOG_ERROR) << "failed to read matrix market cardinality size " << std::endl; 
       exit(1);
    }

    if (type ==  TRAINING){
      ps.M = M; ps.N = N; ps.K = ac.K;
    }
    else if (type ==VALIDATION){
      ps.M_validation = M; ps.N_validation = N;
    }
    else if (type == TEST){
      ps.M_test = M; ps.N_test = N;
    }
 
    if (ps.algorithm == SVD_EXPERIMENTAL && ac.reduce_mem_consumption && ac.svd_compile_eigenvectors)
      return;


    if (_g->num_vertices() == 0){
      init();
      add_vertices(_g, type); 
    }

    /* reseve memory for matrices */

    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

    int I,J; 
    double val;
    int step=nz/10;

    for (i=0; i<nz; i++)
    {
        int rc = fscanf(f, "%d %d %lg\n", &I, &J, &val);
        if (rc != 3)
	  logstream(LOG_FATAL)<<"Failed to read matrix market file: " << filename << " on line: " << i << std::endl;

        //printf("%d) Found row %d %d %g\n", i, I, J, val);        
        I--;  /* adjust from 1-based to 0-based */
        J--;
         if (ac.scalerating != 1.0)
	     val /= ac.scalerating;
         if (!ac.zero)
	   assert(val!=0 );
        
        assert(I>=0 && I< M);
        assert(J >=0 && J< N);
        vertex_data & vdata = _g->vertex_data(I);
        set_new(vdata.datapoint,J, val);   

        if (ps.algorithm == SVD_EXPERIMENTAL){
          vertex_data & other = _g->vertex_data(J + ps.M);
          set_new(other.datapoint, I, val);
	  other.reported = true;
        }
	if (ac.algorithm == K_MEANS){ //compute mean for each cluster by summing assigned points
           ps.clusts.cluster_vec[vdata.current_cluster].cur_sum_of_points[J] += val;  
        }
        if (!vdata.reported){
	      vdata.reported = true;
              if (ac.algorithm == K_MEANS)
	        ps.clusts.cluster_vec[vdata.current_cluster].num_assigned_points++;
              ps.total_assigned++;//count the total number of non-zero rows we encountered
        }
        if (i % step == 0)
          logstream(LOG_INFO) << "Matrix market read " << i << " edges" << std::endl;
     }
    ps.L = nz;
    fclose(f);


    if (ac.reduce_mem_consumption)
      compact(_g);
}

void save_matrix_market_matrix(const char * filename, const flt_dbl_mat & a, std::string comment, bool integer, bool issparse){
    MM_typecode matcode;                        
    int i,j;

    mm_initialize_typecode(&matcode);
    mm_set_matrix(&matcode);
    mm_set_coordinate(&matcode);
    if (issparse)
      mm_set_sparse(&matcode);
    else mm_set_dense(&matcode);  

    if (!integer)
       mm_set_real(&matcode);
    else
       mm_set_integer(&matcode);

    FILE * f = open_file(filename,"w");
    assert(f != NULL);
    mm_write_banner(f, matcode); 
    if (comment.size() > 0)
      fprintf(f, "%s%s", "%", comment.c_str());

    mm_write_mtx_crd_size(f, a.rows(), a.cols(), a.size());
           

    for (i=0; i<a.rows(); i++){
       for (j=0; j<a.cols(); j++){
          if (issparse){
            if (get_val(a,i,j) != 0){
               if (integer)
                 fprintf(f, "%d %d %d\n", i+1, j+1, (int)get_val(a,i,j));
               else  fprintf(f, "%d %d %10.13g\n", i+1, j+1, (double)get_val(a,i,j)); 
            }
          } else { //dense
             if (integer)
                fprintf(f, "%d ", (int)get_val(a,i,j));
             else fprintf(f, "%10.13g ", (double)get_val(a,i,j));
	     if (j == a.cols() -1 )
                fprintf(f, "\n");
          }
       }
    }
    logstream(LOG_INFO) << "Saved output matrix to file: " << filename << std::endl;
    logstream(LOG_INFO) << "You can read it with Matlab/Octave using the script mmread.m found on http://graphlab.org/mmread.m" << std::endl;

}



void save_matrix_market_vector(const char * filename, const flt_dbl_vec & a, std::string comment, bool integer,bool issparse){
    MM_typecode matcode;                        
    int i;

    mm_initialize_typecode(&matcode);
    mm_set_matrix(&matcode);
    mm_set_coordinate(&matcode);
    if (issparse)
      mm_set_sparse(&matcode);
    else 
      mm_set_dense(&matcode);

    if (!integer)
      mm_set_real(&matcode);
    else
      mm_set_integer(&matcode);

    FILE * f = open_file(filename,"w");
    assert(f != NULL);
    mm_write_banner(f, matcode); 
    if (comment.size() > 0)
      fprintf(f, "%s%s", "%", comment.c_str());
     mm_write_mtx_crd_size(f, a.size(), 1, a.size());

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



void save_matrix_market_format(const char * filename)
{
    if (ps.algorithm != SVD_EXPERIMENTAL){
      save_matrix_market_matrix((std::string(filename) + ".clusters.mtx").c_str(),ps.output_clusters,ps.output_clusters_comment, false,false);
      save_matrix_market_matrix((std::string(filename) + ".assignments.mtx").c_str(),ps.output_assignements, ps.output_assignements_comment, ps.output_assignements_integer,false);
    }
    else {
     
      if (!ac.reduce_mem_consumption){
        save_matrix_market_matrix((std::string(filename) + ".V").c_str(),ps.V, ps.V_comment,false,false); /* for conforming to wikipedia convention, I swap U and V*/
        save_matrix_market_matrix((std::string(filename) + ".U").c_str(),ps.U, ps.U_comment,false,false);
      
        save_matrix_market_vector((std::string(filename) + ".EigenValues_AAT").c_str(),get_col(ps.T,0),ps.output_comment3, false,false);
        save_matrix_market_vector((std::string(filename) + ".EigenValues_ATA").c_str(),get_col(ps.T,1),ps.output_comment4, false,false);
      }
    }
}


