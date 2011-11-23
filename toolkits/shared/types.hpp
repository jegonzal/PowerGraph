#ifndef TYPES_COMMON
#define TYPES_COMMON

typedef double real_type;

/*
 * store the size of the matrix
 */
struct matrix_descriptor {
  int rows, cols, nonzeros;

  matrix_descriptor() : rows(0), cols(0), nonzeros(0) { }
   // is the matrix square?
  bool is_square(){ return rows == cols; }
  // get the position of the starting row/col node
  int get_start_node(bool rows){ if (is_square()) return 0; else return rows?0:rows; }
  // get the position of the ending row/col node 
  int get_end_node(bool rows){ if (is_square()) return rows; else return rows?rows:rows+cols; }
  // get howmany row/column nodes
  int howmany(bool rows){ if (is_square()) return rows; else return rows?rows:cols; }
  // how many total nodes
  int total(){ if (is_square()) return rows; else return rows+cols; }

}; // end of matrix descriptor



#endif
