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
#include "pmf.h"
#include "../../libs/matrixmarket/mmio.h"
#include "../gabp/advanced_config.h"
#include <assert.h>

extern advanced_config ac;
extern problem_setup ps;
extern const char * testtypename[];

void load_matrix_market(const char * filename, graph_type *_g, testtype data_type)
{
#ifdef GL_NO_MULT_EDGES
    int ret_code;
    MM_typecode matcode;
    int M, N, nz;   
    int i;

    printf("Loading %s %s\n", filename, testtypename[data_type]);
    FILE * f = fopen(filename, "r");
    if (data_type!=TRAINING && f == NULL){//skip optional files, if the file is missing
      printf("skipping file\n");
      return;
    }

    if(data_type==TRAINING && f== NULL){
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
        logstream(LOG_ERROR) << "sOrry, this application does not support " << std::endl << 
          "Market Market type: " << mm_typecode_to_str(matcode) << std::endl;
        exit(1);
    }

    /* find out size of sparse matrix .... */

    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0){
       logstream(LOG_ERROR) << "failed to read matrix market cardinality size " << std::endl; 
       exit(1);
    }

    ps.M = M; ps.N = N; ps.K = 1;
    verify_size(data_type, M, N, 1);
    add_vertices(_g, data_type); 

    /* reseve memory for matrices */

    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

    int I,J; 
    double val;


    edge_data edge;
    edge.time = 0;
    for (i=0; i<nz; i++)
    {
        fscanf(f, "%d %d %lg\n", &I, &J, &val);
        I--;  /* adjust from 1-based to 0-based */
        J--;
        edge.weight = val;
        if (!ac.zero)
	   assert(val!=0 );
        assert(I< M);
        assert(J< N);
        _g->add_edge(I,J+ps.M,edge);
    }
    set_num_edges(nz, data_type);
    verify_edges(_g, data_type);

    count_all_edges(ps.g);

    fclose(f);
#else
   logstream(LOG_ERROR) << "loading matrix market format is not supported with multiple edges between user and movie mode. Uncomment the flag GL_NO_MULT_EDGES and recompile" << std::endl;
#endif

}

void save_matrix_market_format(const char * filename, mat U, mat V)
{
    MM_typecode matcode;                        
    int i,j;

    mm_initialize_typecode(&matcode);
    mm_set_matrix(&matcode);
    mm_set_coordinate(&matcode);
    mm_set_real(&matcode);

    FILE * f = fopen((std::string(filename) + ".U").c_str(),"w");
    assert(f != NULL);
    mm_write_banner(f, matcode); 
    mm_write_mtx_crd_size(f, ps.M, ac.D, ps.M*ac.D);

    for (i=0; i<ps.M; i++)
       for (j=0; j<ac.D; j++)
        fprintf(f, "%d %d %10.3g\n", i+1, j+1, get_val(U,i,j));

    fclose(f);
    f = fopen((std::string(filename) + ".V").c_str(),"w");
    assert(f != NULL);
    mm_write_banner(f, matcode); 
    mm_write_mtx_crd_size(f, ps.N, ac.D, ps.N*ac.D);

    for (i=0; i<ps.N; i++)
       for (j=0; j<ac.D; j++)
        fprintf(f, "%d %d %10.3g\n", i+1, j+1, get_val(V,i,j));

    fclose(f);

}

