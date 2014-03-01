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


#ifndef _MATH_HPP
#define _MATH_HPP

#include "types.hpp"
#include "graphlab.hpp"
#include "graphlab/util/tracepoint.hpp"


DECLARE_TRACER(Axbtrace);
DECLARE_TRACER(Axbtrace2);
DECLARE_TRACER(vecequals);
DECLARE_TRACER(orthogonalize_vs_alltrace);
DECLARE_TRACER(als_lapack_trace);
DECLARE_TRACER(orth1);
DECLARE_TRACER(orth2);
DECLARE_TRACER(orth3);

double regularization = 0;
bool debug;
bool regnormal;

void print_vec(const char * name, const vec & pvec, bool high);

struct math_info{
  //for Axb operation
  int increment;
  double  c;
  double  d;
  int x_offset, b_offset , y_offset, r_offset, div_offset, prev_offset, mat_offset, vec_offset;
  int orthogonalization;
  bool A_offset, A_transpose;
  std::vector<std::string> names;
  bool use_diag;
  int ortho_repeats;
  int start, end;
  bool update_function;

  //for backslash operation
  bool dist_sliced_mat_backslash;
  mat eDT;
  double maxval, minval;

  math_info(){
    reset_offsets();
  }

  void reset_offsets(){
    increment = 2;
    c=1.0; d=0.0;
    x_offset = b_offset = y_offset = r_offset = div_offset = prev_offset = mat_offset = vec_offset = -1;
    A_offset = false;
    A_transpose = false;
    use_diag = true;
    start = end = -1;
    update_function = false;
    dist_sliced_mat_backslash = false;
    orthogonalization = 0;
  }
  int increment_offset(){
    return increment++;
  }


};


bipartite_graph_descriptor info;
math_info mi;
class DistMat; 
class DistDouble;
class DistSlicedMat;
DistSlicedMat * curMat = NULL;
gather_type alphas;
gather_type sum_alpha;



#define MAX_PRINT_ITEMS 25
double runtime = 0;

using namespace graphlab;

vec curvec;
/***
 * UPDATE FUNCTION (ROWS)
 */
class Axb :
  public graphlab::ivertex_program<graph_type, double>,
  public graphlab::IS_POD_TYPE {
    float last_change;
    public:
    /* Gather the weighted rank of the adjacent page   */
    double gather(icontext_type& context, const vertex_type& vertex,
        edge_type& edge) const {

      if (edge.data().role == edge_data::PREDICT)
         return 0;

      bool brows = vertex.id() < (uint)info.get_start_node(false);
      if (info.is_square()) 
        brows = !mi.A_transpose;
      if (mi.A_offset  && mi.x_offset >= 0){
        double val = edge.data().obs * (brows ? edge.target().data().pvec[mi.x_offset] :
            edge.source().data().pvec[mi.x_offset]);
        //printf("gather edge on vertex %d val %lg obs %lg\n", vertex.id(), val, edge.data().obs);
        return val;
      }
      //printf("edge on vertex %d val %lg\n", vertex.id(), 0.0);
      return 0;
    }

    /* Use the total rank of adjacent pages to update this page */
    void apply(icontext_type& context, vertex_type& vertex,
        const double& total) {

      //printf("Entered apply on node %d value %lg\n", vertex.id(), total);
      vertex_data & user = vertex.data();
      assert(mi.x_offset >=0 || mi.y_offset >= 0);
      assert(mi.r_offset >=0);

      /* perform orthogonalization of current vector */
      if (mi.orthogonalization){
         for (int i=mi.mat_offset; i< mi.vec_offset; i++){
            vertex.data().pvec[mi.vec_offset] -= alphas.pvec[i-mi.mat_offset] * vertex.data().pvec[i]; 
         }
         return;
      }

      double val = total;
      //assert(total != 0 || mi.y_offset >= 0);

      //store previous value for convergence detection
      if (mi.prev_offset >= 0)
        user.pvec[mi.prev_offset ] = user.pvec[mi.r_offset];

      assert(mi.x_offset >=0 || mi.y_offset>=0);
      if (mi.A_offset  && mi.x_offset >= 0){
        if  (info.is_square() && mi.use_diag)// add the diagonal term
          val += (/*mi.c**/ (user.A_ii+ regularization) * user.pvec[mi.x_offset]);
        //printf("node %d added diag term: %lg\n", vertex.id(), user.A_ii);
        val *= mi.c;
      }
      /***** COMPUTE r = c*I*x  *****/
      else if (!mi.A_offset && mi.x_offset >= 0){
        val = mi.c*user.pvec[mi.x_offset];
      }

      /**** COMPUTE r+= d*y (optional) ***/
      if (mi.y_offset>= 0){
        val += mi.d*user.pvec[mi.y_offset]; 
      }

      /***** compute r = (... ) / div */
      if (mi.div_offset >= 0){
        val /= user.pvec[mi.div_offset];
      }

      user.pvec[mi.r_offset] = val;
      //printf("Exit apply on node %d value %lg\n", vertex.id(), val);
    }

    edge_dir_type gather_edges(icontext_type& context,
        const vertex_type& vertex) const {
      if (vertex.id() < rows)
        return OUT_EDGES;
      else return IN_EDGES;
    }


    edge_dir_type scatter_edges(icontext_type& context,
        const vertex_type& vertex) const {
      return NO_EDGES;
    }

    /* The scatter function just signal adjacent pages */
    //void scatter(icontext_type& context, const vertex_type& vertex,
    //    edge_type& edge) const {
    //}

  }; 

void init_lanczos_mapr( graph_type::vertex_type& vertex) {
  assert(actual_vector_len > 0);
  vertex.data().pvec = zeros(actual_vector_len);
} 




void init_math(graph_type * _pgraph, bipartite_graph_descriptor & _info, double ortho_repeats = 3, 
    bool update_function = false){
  pgraph = _pgraph;
  info = _info;
  mi.reset_offsets();
  mi.update_function = update_function;
  mi.ortho_repeats = ortho_repeats;
}


class DistVec{
  public:
    int offset; //real location in memory
    int display_offset; //offset to print out
    int prev_offset;
    std::string name; //optional
    bool transpose;
    bipartite_graph_descriptor info;
    int start; 
    int end;

    void init(){
      start = info.get_start_node(!transpose);
      end = info.get_end_node(!transpose);
      assert(start < end && start >= 0 && end >= 1);
      //debug_print(name);
    };

    int size(){ return end-start; }

    DistVec(const bipartite_graph_descriptor &_info, int _offset, bool _transpose, const std::string & _name){
      offset = _offset;
      display_offset = _offset;
      name = _name;
      info = _info;
      transpose = _transpose;
      prev_offset = -1;
      init();
    }
    DistVec(const bipartite_graph_descriptor &_info, int _offset, bool _transpose, const std::string & _name, int _prev_offset){
      offset = _offset;
      display_offset = _offset;
      name = _name;
      info = _info;
      transpose = _transpose;
      assert(_prev_offset < data_size);
      prev_offset = _prev_offset;
      init();
    }


    DistVec& operator-(){
      mi.d=-1.0;
      return *this; 
    }
    DistVec& operator-(const DistVec & other){
      mi.x_offset = offset;
      mi.y_offset = other.offset;
      transpose = other.transpose;
      if (mi.d == 0)
        mi.d = -1.0;
      else 
        mi.d*=-1.0;
      return *this;
    }
    DistVec& operator+(){
      if (mi.d == 0)
        mi.d=1.0;
      return *this;
    }
    DistVec& orthogonalize(){
      mi.orthogonalization = 1;
      return *this;
    }

    DistVec& operator+(const DistVec &other){
      mi.x_offset =offset;
      mi.y_offset = other.offset;
      transpose = other.transpose;
      return *this; 
    }
    DistVec& operator+(const DistMat &other);

    DistVec& operator-(const DistMat &other);

    DistVec& operator/(const DistVec &other){
      mi.div_offset = other.offset;
      return *this;
    }
    DistVec& operator/(const DistDouble & other);

    DistVec& operator/(double val){
      assert(val != 0);
      assert(mi.d == 0);
      mi.d = 1/val;
      return *this;
    }

    DistVec& operator=(const DistVec & vec);

    DistVec& operator=(const vec & pvec);

    vec to_vec(int dmax = -1, int doffset = -1);


    void debug_print(const char * name){
      if (debug){
        std::cout<<name<<"["<<display_offset<<"]" << std::endl;
        vec pvec = this->to_vec(MAX_PRINT_ITEMS, mi.r_offset == -1? offset:mi.r_offset);
        for (int i=0; i< pvec.size(); i++){  
          //TODO printf("%.5lg ", fabs(pgraph->vertex_data(i).pvec[(mi.r_offset==-1)?offset:mi.r_offset]));
          printf("%.5lg ", fabs(pvec[i]));
        }
        printf("\n");
      }
    }
    void debug_print(std::string name){ return debug_print(name.c_str());}

    double operator[](int i){
      assert(i < end - start);
      assert(false);
      // TODO   return pgraph->vertex_data(i+start).pvec[offset];
    }

    DistDouble operator*(const DistVec & other);

    DistVec& operator*(const double val){
      assert(val!= 0);
      mi.d=val;
      return *this;
    }
    DistVec& operator*(const DistDouble &dval);

    DistMat &operator*(DistMat & v);

    DistVec& _transpose() { 
      /*if (!config.square){
        start = n; end = m+n;
        }*/
      return *this;
    }

    DistVec& operator=(DistMat &mat);

};

DistVec * pcurrent = NULL;


class DistSlicedMat{
  public:
    bipartite_graph_descriptor info;
    int start_offset;
    int end_offset; 
    std::string name; //optional
    int start;
    int end;
    bool transpose;

    DistSlicedMat(int _start_offset, int _end_offset, bool _transpose, const bipartite_graph_descriptor &_info, std::string _name){
      assert(_start_offset < _end_offset);
      assert(_start_offset >= 0);
      assert(_info.total() > 0);
      transpose = _transpose;
      info = _info;
      init();
      start_offset = _start_offset;
      end_offset = _end_offset;
      name = _name;
    }

    DistSlicedMat& operator=(DistMat & other);

    void init(){
      start = info.get_start_node(!transpose);
      end = info.get_end_node(!transpose);
      assert(start < end && start >= 0 && end >= 1);
      //debug_print(name);
    };

    int size(int dim){ return (dim == 1) ? (end-start) : (end_offset - start_offset) ; }

    void set_cols(int start_col, int end_col, const mat& pmat){
      assert(start_col >= 0);
      assert(end_col <= end_offset - start_offset);
      assert(pmat.rows() == end-start);
      assert(pmat.cols() >= end_col - start_col);
      for (int i=start_col; i< end_col; i++)
        this->operator[](i) = get_col(pmat, i-start_col);
    }
    mat get_cols(int start_col, int end_col){
      assert(start_col < end_offset - start_offset);
      assert(start_offset + end_col <= end_offset);
      mat retmat = zeros(end-start, end_col - start_col);
      for (int i=start_col; i< end_col; i++)
        set_col(retmat, i-start_col, this->operator[](i-start_col).to_vec());
      return retmat;
    }

    void operator=(mat & pmat){
      assert(end_offset-start_offset <= pmat.cols());
      assert(end-start == pmat.rows());
      set_cols(0, pmat.cols(), pmat);
    }

    std::string get_name(int pos){
      assert(pos < end_offset - start_offset);
      assert(pos >= 0);
      return name;
    }

    DistVec operator[](int pos){
      assert(pos < end_offset-start_offset);
      assert(pos >= 0);
      DistVec ret(info, start_offset + pos, transpose, get_name(pos));
      ret.display_offset = pos;
      return ret;
    }

};



void assign_vec(graph_type::vertex_type & vertex){
  if (!info.is_square())
    assert(vertex.id() - pcurrent->start >= 0 && vertex.id() - pcurrent->start < curvec.size());
  vertex.data().pvec[pcurrent->offset] = curvec[vertex.id() - pcurrent->start];
}  

gather_type output_vector(const graph_type::vertex_type & vertex){
   assert(pcurrent && pcurrent->offset >= 0 && pcurrent->offset < vertex.data().pvec.size());
   gather_type ret;
   assert(pcurrent->end - pcurrent->start > 0);
   assert(vertex.id() - pcurrent->start >= 0);
   ret.pvec = vec::Zero(pcurrent->end - pcurrent->start);
   ret.pvec[vertex.id() - pcurrent->start] = vertex.data().pvec[pcurrent->offset];
   return ret;
}
bool select_in_range(const graph_type::vertex_type & vertex){
   return vertex.id() >= (uint)pcurrent->start && vertex.id() < (uint)pcurrent->end;
}
DistVec& DistVec::operator=(const DistVec & vec){
      assert(offset < (info.is_square() ? 2*data_size: data_size));
      if (mi.x_offset == -1 && mi.y_offset == -1){
        mi.y_offset = vec.offset;
      }  
      mi.r_offset = offset;
      assert(prev_offset < data_size);
      mi.prev_offset = prev_offset;
      if (mi.d == 0.0)
        mi.d=1.0;
      transpose = vec.transpose;
      end = vec.end; 
      start = vec.start;
      mi.start = start;
      mi.end = end;
      INITIALIZE_TRACER(Axbtrace2, "Update function Axb");
      BEGIN_TRACEPOINT(Axbtrace2);
      pcurrent = (DistVec*)&vec;
      start_engine();
      debug_print(name);
      mi.reset_offsets();
      return *this;
    }

DistVec& DistVec::operator=(const vec & pvec){
  assert(offset >= 0);
  assert(pvec.size() == info.num_nodes(true) || pvec.size() == info.num_nodes(false));
  assert(start < end);
  if (!info.is_square() && pvec.size() == info.num_nodes(false)){
    transpose = true;
  }
  else {
    transpose = false;
  }
  //#pragma omp parallel for    
  INITIALIZE_TRACER(vecequals, "vector assignment");
  BEGIN_TRACEPOINT(vecequals);
  //for (int i=start; i< end; i++){  
  //  pgraph->vertex_data(i).pvec[offset] = pvec[i-start];
  //}
  pcurrent = this;
  curvec = pvec;
  graphlab::vertex_set nodes = pgraph->select(select_in_range);
  pgraph->transform_vertices(assign_vec, nodes);
  END_TRACEPOINT(vecequals);
  debug_print(name);
  return *this;       
}

vec DistVec::to_vec(int dmax, int doffset){
  pcurrent = this;
  if (doffset >= 0)
    pcurrent->offset = doffset;
  if (dmax >= 0)
    pcurrent->end = std::min(pcurrent->start + dmax, pcurrent->end);
  graphlab::vertex_set nodes = pgraph->select(select_in_range);
  //    for (int i=start; i< end; i++){
  //      //TODO ret[i-start] = pgraph->vertex_data(i).pvec[offset];
  //    }
  gather_type curvec = pgraph->map_reduce_vertices<gather_type>(output_vector, nodes);
  return curvec.pvec;
}



/*
 * wrapper for computing r = c*A*x+d*b*y
 */
class DistMat{
  public:
    bool transpose;
    bipartite_graph_descriptor info;

    DistMat(const bipartite_graph_descriptor& _info) { 
      info = _info;
      transpose = false;
    };


    DistMat &operator*(const DistVec & v){
      mi.x_offset = v.offset;
      mi.A_offset = true;
      //v.transpose = transpose;
      //r_offset = A_offset;
      return *this;
    }
    DistMat &operator*(const DistDouble &d);

    DistMat &operator-(){
      mi.c=-1.0;
      return *this;
    }

    DistMat &operator/(const DistVec & v){
      mi.div_offset = v.offset;
      return *this;
    }

    DistMat &operator+(){
      mi.c=1.0;
      return *this;
    }
    DistMat &operator+(const DistVec &v){
      mi.y_offset = v.offset;
      if (mi.d == 0.0)
        mi.d=1.0;
      return *this;
    }
    DistMat &operator-(const DistVec &v){
      mi.y_offset = v.offset;
      if (mi.d == 0.0)
        mi.d=-1.0;
      else 
        mi.d*=-1.0;
      return *this;
    }
    DistMat & _transpose(){
      transpose = true;
      mi.A_transpose = true;
      return *this;
    }
    DistMat & operator~(){
      return _transpose();
    }
    DistMat & backslash(DistSlicedMat & U){
      mi.dist_sliced_mat_backslash = true;
      transpose = U.transpose;
      return *this;
    }
    void set_use_diag(bool use){
      mi.use_diag = use;
    }   
};


DistVec& DistVec::operator=(DistMat &mat){
  mi.r_offset = offset;
  assert(prev_offset < data_size);
  mi.prev_offset = prev_offset;
  transpose = mat.transpose;
  mi.start = info.get_start_node(!transpose);
  mi.end = info.get_end_node(!transpose);
  INITIALIZE_TRACER(Axbtrace, "Axb update function");
  BEGIN_TRACEPOINT(Axbtrace);
  pcurrent = this;
  int old_start = start; int old_end = end;
  start = mi.start; end = mi.end;
  start_engine();
  start = old_start; end = old_end;
  END_TRACEPOINT(Axbtrace);
  debug_print(name);
  mi.reset_offsets();
  mat.transpose = false;
  return *this;
}
DistVec& DistVec::operator+(const DistMat &other){
  mi.y_offset = offset;
  transpose = other.transpose;
  return *this; 
}
DistVec& DistVec::operator-(const DistMat & other){
  mi.y_offset = offset;
  transpose = other.transpose;
  if (mi.c == 0)
    mi.c = -1;
  else mi.c *= -1;
  return *this;
}

DistMat& DistVec::operator*(DistMat & v){
  mi.x_offset = offset;
  mi.A_offset = true;
  return v;
}


class DistDouble{
  public:
    double val;
    std::string name;

    DistDouble() {};
    DistDouble(double _val) : val(_val) {};


    DistVec& operator*(DistVec & dval){
      mi.d=val;
      return dval;
    }
    DistMat& operator*(DistMat & mat){
      mi.c = val;
      return mat;
    }

    DistDouble  operator/(const DistDouble dval){
      DistDouble mval;
      mval.val = val / dval.val;
      return mval;
    }
    bool operator<(const double other){
      return val < other;
    }
    DistDouble & operator=(const DistDouble & other){
      val = other.val;
      debug_print(name);
      return *this;
    }
    bool operator==(const double _val){
      return val == _val;
    }
    void debug_print(const char * name){
      std::cout<<name<<" "<<val<<std::endl;
    }
    double toDouble(){
      return val;
    }
    void debug_print(std::string name){ return debug_print(name.c_str()); }


};

DistDouble DistVec::operator*(const DistVec & vec){
  mi.y_offset = offset;
  mi.b_offset = vec.offset;
  if (mi.d == 0) 
    mi.d = 1.0;
  assert(mi.y_offset >=0 && mi.b_offset >= 0);

  double val = 0;
  for (int i=start; i< end; i++){  
    assert(false);//not yet
    //TODO const vertex_data * data = &pgraph->vertex_data(i);
    //TODO double * pv = (double*)&data->pvec[0];
    //TODO val += mi.d* pv[mi.y_offset] * pv[mi.b_offset];
  }
  mi.reset_offsets();
  DistDouble mval;
  mval.val = val;
  return mval;
}
DistVec& DistVec::operator*(const DistDouble &dval){
  mi.d = dval.val;
  return *this;
}


int size(DistMat & A, int pos){
  assert(pos == 1 || pos == 2);
  return A.info.num_nodes(!A.transpose);
}

DistMat &DistMat::operator*(const DistDouble &d){
  mi.c = d.val;
  return *this;
}

DistDouble sqrt(DistDouble & dval){
  DistDouble mval;
  mval.val=sqrt(dval.val);
  return mval;
}

gather_type calc_norm(const graph_type::vertex_type & vertex){
  gather_type ret;
  assert(pcurrent && pcurrent->offset < vertex.data().pvec.size());
  ret.training_rmse = pow(vertex.data().pvec[pcurrent->offset], 2);
  return ret;
}

DistDouble norm(const DistVec &vec){
  assert(vec.offset>=0);
  assert(vec.start < vec.end);

  DistDouble mval;
  mval.val = 0;
  pcurrent = (DistVec*)&vec;
  vertex_set nodes = pgraph->select(select_in_range);
  //for (int i=vec.start; i < vec.end; i++){
    // TODO const vertex_data * data = &pgraph->vertex_data(i);
    //double * px = (double*)&data->pvec[0];
    // mval.val += px[vec.offset]*px[vec.offset];
  gather_type ret = pgraph->map_reduce_vertices<gather_type>(calc_norm);
  //}
  mval.val = sqrt(ret.training_rmse);
  return mval;
}


DistDouble norm(DistMat & mat){
  DistVec vec(info, 0, mat.transpose, "norm");
  vec = mat;
  return norm((const DistVec&)vec);
}

vec diag(DistMat & mat){
  assert(info.is_square());
  vec ret = zeros(info.total());
  for (int i=0; i< info.total(); i++){
    //TODO ret[i] = pgraph->vertex_data(i).A_ii;
    assert(false);
  }
  return ret;
}

int curoffset = -1;
gather_type map_reduce_ortho(const graph_type::vertex_type & vertex){
  gather_type ret;
  assert(curoffset >= 0);
  assert(curMat && curMat->start_offset - pcurrent->offset);
  ret.pvec = vec::Zero(curoffset);
  assert(curMat != NULL && curMat->start_offset < pcurrent->offset);
  //for (int i=mat.start_offset; i< current.offset; i++){
  for (int i=curMat->start_offset; i< pcurrent->offset; i++){
    ret.pvec[i - curMat->start_offset] = vertex.data().pvec[i] * vertex.data().pvec[pcurrent->offset];
  }
  //printf("map_Reduce_ortho: node %d\n", vertex.id());
  //std::cout<<ret.pvec<<std::endl;
  return ret;
}
  gather_type map_reduce_sum_power(const graph_type::vertex_type & vertex){
    gather_type ret;
    assert(pcurrent->offset >= 0 && pcurrent->offset < vertex.data().pvec.size());
    ret.training_rmse = pow(vertex.data().pvec[pcurrent->offset], 2);
    return ret;
  }
  void divide_by_sum(graph_type::vertex_type& vertex){
    assert(pcurrent->offset >= 0 && pcurrent->offset < vertex.data().pvec.size());
    vertex.data().pvec[pcurrent->offset] /= sum_alpha.training_rmse;
  }



  void transform_ortho(graph_type::vertex_type & vertex){
    assert(curMat != NULL && curMat->start_offset < pcurrent->offset);
    for (int i=curMat->start_offset; i< pcurrent->offset; i++){
      //assert(alphas.pvec[i-curMat->start_offset] != 0);
      vertex.data().pvec[pcurrent->offset] -= alphas.pvec[i-curMat->start_offset] * vertex.data().pvec[i]; 
    }
  }

  bool selected_node(const graph_type::vertex_type& vertex){
    if (info.is_square())
      return true;
    else return ((vertex.id() >= (uint)info.get_start_node(!pcurrent->transpose)) &&
        (vertex.id() < (uint)info.get_end_node(!pcurrent->transpose)));
  }


  double orthogonalize_vs_all(DistSlicedMat & mat, int _curoffset, double &alpha){
    assert(mi.ortho_repeats >=1 && mi.ortho_repeats <= 3);
    curoffset = _curoffset;
    curMat = &mat;
    mi.mat_offset = mat.start_offset;
    INITIALIZE_TRACER(orthogonalize_vs_alltrace, "orthogonalization step - optimized");
    BEGIN_TRACEPOINT(orthogonalize_vs_alltrace);
    bool old_debug = debug;
    debug = false;
    DistVec current = mat[curoffset];
    pcurrent =&current;
    mi.vec_offset = pcurrent->offset;
    assert(mat.start_offset <= current.offset); 
    vertex_set nodes = pgraph->select(selected_node);
    if (curoffset > 0){
      for (int j=0; j < mi.ortho_repeats; j++){
        INITIALIZE_TRACER(orth1, "map reduce in ortho");
        BEGIN_TRACEPOINT(orth1);
        alphas = pgraph->map_reduce_vertices<gather_type>(map_reduce_ortho, nodes);
        END_TRACEPOINT(orth1);
        //pgraph->transform_vertices(transform_ortho, nodes);
        mat[_curoffset] = mat[_curoffset].orthogonalize(); 
      } //for ortho_repeast 
    }

    debug = old_debug;
    current.debug_print(current.name);
    INITIALIZE_TRACER(orth2, "map reduce in ortho2");
    BEGIN_TRACEPOINT(orth2);
    sum_alpha = pgraph->map_reduce_vertices<gather_type>(map_reduce_sum_power, nodes);
    END_TRACEPOINT(orth2);
    sum_alpha.training_rmse = sqrt(sum_alpha.training_rmse);
    alpha = sum_alpha.training_rmse;
    if (alpha >= 1e-10 ){
       INITIALIZE_TRACER(orth3, "transform_vertices in ortho3");
       BEGIN_TRACEPOINT(orth3);
       //pgraph->transform_vertices(divide_by_sum, nodes);    
       mat[_curoffset] = mat[_curoffset] / alpha;
       END_TRACEPOINT(orth3);
    }
    END_TRACEPOINT(orthogonalize_vs_alltrace);
    return alpha;
  }


  DistVec& DistVec::operator/(const DistDouble & other){
    assert(other.val != 0);
    assert(mi.d == 0);
    mi.d = 1/other.val;
    return *this;
  }

#endif //_MATH_HPP
