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
#include "mathlayer.hpp"
#include "graphlab.hpp"


extern double regularization;
extern bool debug;

struct math_info{
  int increment;
  double  c;
  double  d;
  int x_offset, b_offset , y_offset, r_offset;
  bool A_offset;
  std::vector<std::string> names;

  math_info(){
    reset_offsets();
  }

  void reset_offsets(){
    increment = 2;
    c=1.0; d=0.0;
    x_offset = b_offset = y_offset = r_offset = -1;
    A_offset = false;
  }
  int increment_offset(){
    return increment++;
  }


};


bipartite_graph_descriptor info;
math_info mi;

#define MAX_PRINT_ITEMS 21
#define MAX_OFFSET 7
double runtime = 0;

using namespace graphlab;
//std::vector<vertex_id_t> rows,cols;

typedef graph<vertex_data,edge_data>::edge_list_type edge_list_type;
typedef graph<vertex_data,edge_data>::edge_type edge_type;

graph_type * pgraph;

/***
 * UPDATE FUNCTION (ROWS)
 */
struct Axb:
 public iupdate_functor<graph_type, Axb> {
 
  void operator()(icontext_type &context){

  vertex_data& user = context.vertex_data();
  bool rows = context.vertex_id() < (uint)info.get_start_node(false);
  assert(mi.r_offset >=0);
  double val = 0;
  assert(mi.x_offset >=0 || mi.y_offset>=0);
  timer t; t.start();
  
  /*** COMPUTE r = c*A*x  ********/
  if (mi.A_offset  && mi.x_offset >= 0){
   edge_list_type edges = rows?
	context.out_edges() : context.in_edges(); 
   for (size_t i = 0; i < edges.size(); i++){
      const edge_data & edge = context.edge_data(edges[i]);
      const vertex_data  & movie = context.const_vertex_data(rows ? edges[i].target() : edges[i].source());
      val += (mi.c * edge.weight * movie.pvec[mi.x_offset]);
    }
  
    if (info.is_square())// add the diagonal term
      val += (mi.c* (user.A_ii+ regularization) * user.pvec[mi.x_offset]);
  }
 /***** COMPUTE r = c*I*x  *****/
  else if (!mi.A_offset && mi.x_offset >= 0){
     val = mi.c*user.pvec[mi.x_offset];
  }
  
  /**** COMPUTE r+= d*y (optional) ***/
  if (mi.y_offset>= 0){
    val += mi.d*user.pvec[mi.y_offset]; 
  }

  user.pvec[mi.r_offset] = val;
}
};

core<graph_type, Axb> * glcore = NULL;
void init_math(graph_type * _pgraph, core<graph_type, Axb> * _glcore){
  pgraph = _pgraph;
  glcore = _glcore;
}


class DistMat; 
class DistDouble;

class DistVec{
   public:
   int offset;
   std::string name; //optional
   bool transpose;
   bipartite_graph_descriptor info;
   int start; 
   int end;

   void init(){
     start = info.get_start_node(!transpose);
     end = info.get_end_node(!transpose);
     assert(start < end && start >= 0 && end >= 1);
     debug_print(name);
   };

   DistVec(bipartite_graph_descriptor &_info, int _offset, bool _transpose, const std::string & _name){
     offset = _offset;
     name = _name;
     info = _info;
     transpose = _transpose;
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

   DistVec& operator=(const DistVec & vec){
     assert(offset < MAX_OFFSET);
     if (mi.x_offset == -1 && mi.y_offset == -1){
         mi.y_offset = vec.offset;
       }  
      mi.r_offset = offset;
      if (mi.d == 0.0)
        mi.d=1.0;
      transpose = vec.transpose;
      end = vec.end; 
      start = vec.start;
      for (vertex_id_type i = start; i < (vertex_id_type)end; i++)
        glcore->schedule(i, Axb()); 
      runtime += glcore->start();
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
    for (int i=start; i< end; i++){  
         pgraph->vertex_data(i).pvec[offset] = pvec[i-start];
    }
    debug_print(name);
    return *this;       
  }


    vec to_vec(){
      vec ret = zeros(end-start);
      for (int i=start; i< end; i++){
         ret[i-start] = pgraph->vertex_data(i).pvec[offset];
      }
      return ret;
    }

  void debug_print(const char * name){
     if (debug){
       std::cout<<name<<" ("<<name<<" "<<offset<<" [ " << (end-start) << "] ";
       for (int i=start; i< std::min(end, start+MAX_PRINT_ITEMS); i++){  
         std::cout<<pgraph->vertex_data(i).pvec[mi.r_offset==-1?offset:mi.r_offset]<<" ";
       }
       std::cout<<std::endl;
     }
  }
  void debug_print(std::string name){ return debug_print(name.c_str());}

  DistDouble operator*(const DistVec & other);
  
  DistVec& operator*(const double val){
     mi.d=val;
     return *this;
  }
  DistVec& _transpose() { 
     /*if (!config.square){
       start = n; end = m+n;
     }*/
     return *this;
  }

  DistVec& operator=(DistMat &mat);
 
 };



/*
 * wrapper for computing r = c*A*x+d*b*y
 */
class DistMat{
  public:
    bool transpose;
    bipartite_graph_descriptor info;

    DistMat(bipartite_graph_descriptor& _info) { 
      info = _info;
      transpose = false;
    };


    DistMat &operator*(DistVec & v){
      	mi.x_offset = v.offset;
        mi.A_offset = true;
        //v.transpose = transpose;
        //r_offset = A_offset;
        return *this;
    }
    DistMat &operator-(){
        mi.c=-1.0;
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
       return *this;
    }
   
};

DistVec& DistVec::operator=(DistMat &mat){
  mi.r_offset = offset;
  transpose = mat.transpose;
  for (vertex_id_type start = info.get_start_node(!transpose); start< (vertex_id_type)info.get_end_node(!transpose); start++)
    glcore->schedule(start, Axb());
  runtime += glcore->start();
  debug_print(name);
  mi.reset_offsets();
  mat.transpose = false;
  return *this;
}


class DistDouble{
  public:
     double val;
     std::string name;

     DistDouble() {};
     DistDouble(double _val) : val(_val) {};
   
   
     const DistVec& operator*(const DistVec & dval){
        mi.d=val;
        return dval;
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
         const vertex_data * data = &pgraph->vertex_data(i);
         double * pv = (double*)&data->pvec[0];
        // if (y_offset >= 0 && b_offset == -1)
	     //val += pv[y_offset] * pv[y_offset];
         val += mi.d* pv[mi.y_offset] * pv[mi.b_offset];
      }
      mi.reset_offsets();
      DistDouble mval;
      mval.val = val;
      return mval;
 }


int size(DistMat & A, int pos){
   assert(pos == 1 || pos == 2);
   return A.info.num_nodes(!A.transpose);
}

DistDouble sqrt(DistDouble & dval){
    DistDouble mval;
    mval.val=sqrt(dval.val);
    return mval;
}

DistDouble norm(DistVec & vec){
    assert(vec.offset>=0);
    assert(vec.start < vec.end);

    DistDouble mval;
    mval.val = 0;
    for (int i=vec.start; i < vec.end; i++){
       const vertex_data * data = &pgraph->vertex_data(i);
       double * px = (double*)&data->pvec[0];
       mval.val += px[vec.offset]*px[vec.offset];
    }
    mval.val = sqrt(mval.val);
    return mval;
}



#endif //_MATH_HPP
