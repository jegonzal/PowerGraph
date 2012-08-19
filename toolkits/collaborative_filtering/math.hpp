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
//#include "printouts.hpp"


DECLARE_TRACER(Axbtrace);
DECLARE_TRACER(Axbtrace2);
DECLARE_TRACER(vecequals);
DECLARE_TRACER(orthogonalize_vs_alltrace);
DECLARE_TRACER(als_lapack_trace);

double regularization;
bool debug;
bool regnormal;

void print_vec(const char * name, const vec & pvec, bool high);

struct math_info{
  //for Axb operation
  int increment;
  double  c;
  double  d;
  int x_offset, b_offset , y_offset, r_offset, div_offset, prev_offset, div_const;
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
    x_offset = b_offset = y_offset = r_offset = div_offset = prev_offset = -1;
    div_const = 0;
    A_offset = false;
    A_transpose = false;
    use_diag = true;
    start = end = -1;
    update_function = false;
    dist_sliced_mat_backslash = false;
  }
  int increment_offset(){
    return increment++;
  }


};


bipartite_graph_descriptor info;
math_info mi;

#define MAX_PRINT_ITEMS 25
double runtime = 0;

using namespace graphlab;

//typedef vertex_data vertex_data_type;
//typedef edge_data edge_data_type;

graph_type * pgraph = NULL;

/***
 * UPDATE FUNCTION (ROWS)
 */
class Axb :
  public graphlab::ivertex_program<graph_type, float>,
  public graphlab::IS_POD_TYPE {
  float last_change;
public:
  /* Gather the weighted rank of the adjacent page   */
  float gather(icontext_type& context, const vertex_type& vertex,
               edge_type& edge) const {
       bool rows = vertex.id() < (uint)info.get_start_node(false);
       if (info.is_square()) 
         rows = mi.A_transpose;
       if (mi.A_offset  && mi.x_offset >= 0){
         return edge.data().obs * (rows ? edge.target().data().pvec[mi.x_offset] :
				      edge.source().data().pvec[mi.x_offset]);
       }
  }

  /* Use the total rank of adjacent pages to update this page */
  void apply(icontext_type& context, vertex_type& vertex,
             const gather_type& total) {
      
     vertex_data & user = vertex.data();
     assert(mi.r_offset >=0);

     double val = total;

      //store previous value for convergence detection
      if (mi.prev_offset >= 0)
        user.pvec[mi.prev_offset ] = user.pvec[mi.r_offset];

      assert(mi.x_offset >=0 || mi.y_offset>=0);
      if (mi.A_offset  && mi.x_offset >= 0){
          if  (info.is_square() && mi.use_diag)// add the diagonal term
            val += (/*mi.c**/ (user.A_ii+ regularization) * user.pvec[mi.x_offset]);

          val *= mi.c;
      }
      /***** COMPUTE r = c*I*x  *****/
      else if (!mi.A_offset && mi.x_offset >= 0){
          val = mi.c*user.pvec[mi.x_offset];
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
  
  }

  edge_dir_type scatter_edges(icontext_type& context,
                              const vertex_type& vertex) const {
    bool rows = mi.end < info.get_start_node(false);
    if (info.is_square()) 
      rows = mi.A_transpose;
    return rows ? graphlab::OUT_EDGES : graphlab::IN_EDGES;
  }

  /* The scatter function just signal adjacent pages */
  void scatter(icontext_type& context, const vertex_type& vertex,
               edge_type& edge) const {
    //context.signal(edge.target());
  }
}; 



 
#if 0
struct Axb:
  public iupdate_functor<graph_type, Axb> {
 
  void operator()(icontext_type &context){
      if (context.vertex_id() < (uint)mi.start || context.vertex_id() >= (uint)mi.end)
        return;

      vertex_data& user = context.vertex_data();
      bool rows = context.vertex_id() < (uint)info.get_start_node(false);
      if (info.is_square()) 
        rows = mi.A_transpose;
      assert(mi.r_offset >=0);
      //store previous value for convergence detection
      if (mi.prev_offset >= 0)
        user.pvec[mi.prev_offset ] = user.pvec[mi.r_offset];

      double val = 0;
      assert(mi.x_offset >=0 || mi.y_offset>=0);
      timer t; t.start();
 
      /*** COMPUTE r = c*A*x  ********/
      if (mi.A_offset  && mi.x_offset >= 0){
        edge_list edges = rows?
          context.out_edges() : context.in_edges(); 
        for (size_t i = 0; i < edges.size(); i++){
            const edge_data & edge = context.edge_data(edges[i]);
            const vertex_data  & movie = context.const_vertex_data(rows ? edges[i].target() : edges[i].source());
            val += (edge.weight * movie.pvec[mi.x_offset]);
          }
  
          if  (info.is_square() && mi.use_diag)// add the diagonal term
            val += (/*mi.c**/ (user.A_ii+ regularization) * user.pvec[mi.x_offset]);

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
      }
      void operator+=(const Axb& other) { 
      }

      void finalize(iglobal_context_type& context) {
      } 
    };
#endif
    

    void init_math(graph_type * _pgraph, bipartite_graph_descriptor & _info, double ortho_repeats = 3, 
                   bool update_function = false){
      pgraph = _pgraph;
      info = _info;
      mi.reset_offsets();
      mi.update_function = update_function;
      mi.ortho_repeats = ortho_repeats;
    }


    class DistMat; 
    class DistDouble;

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


      DistVec& operator=(const DistVec & vec){
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
        pengine->signal_all();
        pengine->start();
        debug_print(name);
        mi.reset_offsets();
        return *this;
      }

      DistVec& operator=(const vec & pvec){
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
        for (int i=start; i< end; i++){  
          //TODO pgraph->vertex_data(i).pvec[offset] = pvec[i-start];
        }
        END_TRACEPOINT(vecequals);
        debug_print(name);
        return *this;       
      }


      vec to_vec(){
        vec ret = vec::Zero(end-start);
        for (int i=start; i< end; i++){
          //TODO ret[i-start] = pgraph->vertex_data(i).pvec[offset];
        }
        return ret;
      }

      double get_pos(int i){
        //TODO return pgraph->vertex_data(i).pvec[offset];
        return 0;
      }

      void debug_print(const char * name){
        if (debug){
          std::cout<<name<<"["<<display_offset<<"]" << std::endl;
          for (int i=start; i< std::min(end, start+MAX_PRINT_ITEMS); i++){  
            //TODO printf("%.5lg ", fabs(pgraph->vertex_data(i).pvec[(mi.r_offset==-1)?offset:mi.r_offset]));
          }
          printf("\n");
        }
      }
      void debug_print(std::string name){ return debug_print(name.c_str());}

      double operator[](int i){
        assert(i < end - start);
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
      //TODO glcore->aggregate_now("Axb");
      pengine->signal_all();
      pengine->start();
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
    DistDouble norm(const DistVec &vec){
      assert(vec.offset>=0);
      assert(vec.start < vec.end);

      DistDouble mval;
      mval.val = 0;
      for (int i=vec.start; i < vec.end; i++){
        // TODO const vertex_data * data = &pgraph->vertex_data(i);
        //double * px = (double*)&data->pvec[0];
        // mval.val += px[vec.offset]*px[vec.offset];
      }
      mval.val = sqrt(mval.val);
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
      }
      return ret;
    }
    void orthogonalize_vs_all(DistSlicedMat & mat, int curoffset, double &alpha){
      assert(mi.ortho_repeats >=1 && mi.ortho_repeats <= 3);
      INITIALIZE_TRACER(orthogonalize_vs_alltrace, "orthogonalization step - optimized");
      BEGIN_TRACEPOINT(orthogonalize_vs_alltrace);
      bool old_debug = debug;
      debug = false;
      DistVec current = mat[curoffset];
      assert(mat.start_offset <= current.offset); 
      double * alphas = new double[curoffset];
      //DistDouble * alphas = new DistDouble[curoffset];
      //cout<<current.to_vec().transpose() << endl;
      if (curoffset > 0){
        for (int j=0; j < mi.ortho_repeats; j++){
          memset(alphas, 0, sizeof(double)*curoffset);
#pragma omp parallel for
          for (int i=mat.start_offset; i< current.offset; i++){
            for (int k=info.get_start_node(!current.transpose); k< info.get_end_node(!current.transpose); k++){
              //TODO alphas[i-mat.start_offset] += pgraph->vertex_data(k).pvec[i] * pgraph->vertex_data(k).pvec[current.offset];
            }
          }
          for (int i=mat.start_offset; i< current.offset; i++){
#pragma omp parallel for
            for (int k=info.get_start_node(!current.transpose); k< info.get_end_node(!current.transpose); k++){
              //TODO pgraph->vertex_data(k).pvec[current.offset] -= alphas[i-mat.start_offset]  * pgraph->vertex_data(k).pvec[i];
            }
          }
        } //for ortho_repeast 
      }

      delete [] alphas; 
      debug = old_debug;
      current.debug_print(current.name);
      //    alpha = 0;
      double sum = 0;
      int k;
      //#pragma omp parallel for private(k) reduction(+: sum)
      for (k=info.get_start_node(!current.transpose); k< info.get_end_node(!current.transpose); k++){
        //TODO sum = sum + pow(pgraph->vertex_data(k).pvec[current.offset],2);
      }    
      alpha = sqrt(sum);
      if (alpha >= 1e-10 ){
#pragma omp parallel for
        for (int k=info.get_start_node(!current.transpose); k< info.get_end_node(!current.transpose); k++){
          //TODO pgraph->vertex_data(k).pvec[current.offset]/=alpha;
        }    
      }
      END_TRACEPOINT(orthogonalize_vs_alltrace);
    }
    void multiply(DistSlicedMat & mat, int curoffset, double a){

      INITIALIZE_TRACER(multiply, "multiply sliced mat");
      BEGIN_TRACEPOINT(multiply);
  
      assert(a>0);
      DistVec current = mat[curoffset];
      assert(mat.start_offset <= current.offset); 
      vec result = zeros(curoffset);
 
      if (curoffset > 0){

#pragma omp parallel for
        for (int i=mat.start_offset; i< current.offset; i++){
          for (int k=info.get_start_node(!current.transpose); k< info.get_end_node(!current.transpose); k++){
            //TODO result[i-mat.start_offset] += pgraph->vertex_data(k).pvec[i] * pgraph->vertex_data(k).pvec[current.offset];
          }
        }
#pragma omp parallel for
        for (int k=info.get_start_node(!current.transpose); k< info.get_end_node(!current.transpose); k++){
          //TODO pgraph->vertex_data(k).pvec[curoffset] /= a;
        }
     
        for (int i=mat.start_offset; i< current.offset; i++){
#pragma omp parallel for
          for (int k=info.get_start_node(!current.transpose); k< info.get_end_node(!current.transpose); k++){
            //TODO pgraph->vertex_data(k).pvec[current.offset] -= result[i-mat.start_offset]/a  * pgraph->vertex_data(k).pvec[i];
          }
        }
      }

      current.debug_print(current.name);
      END_TRACEPOINT(multiply);
    }



    DistVec& DistVec::operator/(const DistDouble & other){
      assert(other.val != 0);
      assert(mi.d == 0);
      mi.d = 1/other.val;
      return *this;
    }
 
#endif //_MATH_HPP
