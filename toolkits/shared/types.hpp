#ifndef TYPES_COMMON
#define TYPES_COMMON

typedef double real_type;

/*
 * store a matrix is a bipartite graph. One side is the rows and the other is the column.
 */
struct bipartite_graph_descriptor {
  int rows, cols, nonzeros;

  bipartite_graph_descriptor() : rows(0), cols(0), nonzeros(0) { }
   // is the matrix square?
  bool is_square(){ return rows == cols; }
  // get the position of the starting row/col node
  int get_start_node(bool rows){ if (is_square()) return 0; else return rows?0:rows; }
  // get the position of the ending row/col node 
  int get_end_node(bool rows){ if (is_square()) return rows; else return rows?rows:rows+cols; }
  // get howmany row/column nodes
  int num_nodes(bool _rows){ if (_rows) return rows; else return cols; }
  // how many total nodes
  int total(){ if (is_square()) return rows; else return rows+cols; }
  //is this a row node
  bool is_row_node(int id){ return id < rows; }
  //debug print?
  bool toprint(int id){ return (id == 0) || (id == rows - 1) || (id == rows) || (id == rows+cols-1); }
  
}; // end of bipartite graph descriptor



#endif
