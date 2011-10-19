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

    ps.M = M; ps.N = N; ps.K = ac.K;

    init();
    add_vertices(_g, type); 

    /* reseve memory for matrices */

    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

    int I,J; 
    double val;


    for (i=0; i<nz; i++)
    {
        fscanf(f, "%d %d %lg\n", &I, &J, &val);
        //printf("Found row %d %d %g\n", I, J, val);        
        I--;  /* adjust from 1-based to 0-based */
        J--;
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

    }
    ps.L = nz;
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

    ps.M = M; ps.N = N; ps.K = ac.K;

    init();
    add_vertices(_g, type); 

    /* reseve memory for matrices */

    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

    int I,J; 
    double val;


    for (i=0; i<nz; i++)
    {
        fscanf(f, "%d %d %lg\n", &I, &J, &val);
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
    }
    ps.L = nz;
    fclose(f);


    if (ac.reduce_mem_consumption)
      compact(_g);
}

void save_matrix_market_format(const char * filename)
{
    MM_typecode matcode;                        
    int i,j;

    mm_initialize_typecode(&matcode);
    mm_set_matrix(&matcode);
    mm_set_coordinate(&matcode);
    mm_set_real(&matcode);

    FILE * f = fopen((std::string(filename) + ".clusters.mtx").c_str(),"w");
    assert(f != NULL);
    mm_write_banner(f, matcode); 
    mm_write_mtx_crd_size(f, ps.K, ps.N, ps.K*ps.N);

    for (i=0; i<ps.K; i++)
       for (j=0; j<ps.N; j++)
        fprintf(f, "%d %d %10.3g\n", i+1, j+1, get_val(ps.output_clusters,i,j));

    fclose(f);
    f = fopen((std::string(filename) + ".assignments.mtx").c_str(),"w");
    assert(f != NULL);
    mm_write_banner(f, matcode); 
    int rows = ps.output_assignements.rows();
    int cols = ps.output_assignements.cols();

    mm_write_mtx_crd_size(f, rows, cols, rows*cols);

    for (i=0; i< rows; i++)
    for (j=0; j< cols; j++)
        if (get_val(ps.output_assignements,i,j) > 0)
          fprintf(f, "%d %d %10.3g\n", i+1, j+1, get_val(ps.output_assignements,i,j));

    fclose(f);

}

