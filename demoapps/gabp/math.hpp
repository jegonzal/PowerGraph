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

#include <graphlab/macros_def.hpp>
int increment=2;
double  c=1.0;
double  d=0.0;
extern bool debug;
extern bool square;
int x_offset = -1, b_offset = -1, y_offset = -1, r_offset = -1;
bool A_offset = false;
extern uint n; 
extern uint m;
gl_types::core * glcore;
#define MAX_PRINT_ITEMS 21
#define MAX_OFFSET 7
double runtime = 0;
const char * names[]={"b","","r","p","x","Ap","t"};
extern bool cg_noop;

using namespace graphlab;
std::vector<vertex_id_t> rows,cols;
void Axb(gl_types::iscope &scope, 
	 gl_types::icallback &scheduler);
void fast_Axb(graph_type * g, std::vector<vertex_id_t> nodes);

void reset_offsets(){
  c=1.0; d=0.0;
  x_offset = b_offset = y_offset = r_offset = -1;
  A_offset = false;
  if (debug)
	  std::cout<<"reset offsets"<<std::endl;
}

int increment_offset(){
  return increment++;
}


class DistMat; 
class DistDouble;

class DistVec{
   public:
   int offset;
   std::string name; //optional
   bool transpose;
   int start;
   int end;

   void init(){
     if (square || !transpose){
       start = 0;
       end = n;
     }
     else {
       start = n;
       end = m+n;
     }
     debug_print(name);
   };

   DistVec(){
     offset = increment_offset();
     transpose = false;
     init();
   }

   DistVec(bool _transpose){
     transpose = _transpose;
     offset = increment_offset();
     init();
   }

   DistVec(int _offset){
     offset = _offset;
     transpose = false;
     init();
   }   

   DistVec(int _offset, bool _transpose){
     offset = _offset;
     transpose = _transpose;
     init();
   }

   DistVec& operator-(){
     d=-1.0;
     return *this; 
   }
   DistVec& operator-(const DistVec & other){
     x_offset = offset;
     y_offset = other.offset;
     transpose = other.transpose;
     start = other.start;
     end = other.end;
     if (d == 0)
       d = -1.0;
     else 
       d*=-1.0;
     return *this;
   }
   DistVec& operator+(){
     if (d == 0)
       d=1.0;
     return *this;
   }
   DistVec& operator+(const DistVec &other){
      x_offset =offset;
      y_offset = other.offset;
      start = other.start;
      end = other.end;
      return *this; 
   }

   DistVec& operator=(const DistVec & vec){
     assert(offset < MAX_OFFSET);
     if (x_offset == -1 && y_offset == -1){
         y_offset = vec.offset;
       }  
      r_offset = offset;
      if (d == 0.0)
        d=1.0;
      transpose = vec.transpose;
      start = vec.start;
      end = vec.end;
      if (cg_noop && !A_offset)
         fast_Axb(
               &glcore->graph(),
               (!transpose)?rows:cols); 
      else   
         glcore->add_tasks((!transpose)?rows:cols, Axb, 1); 
      
      runtime += glcore->start();
      debug_print(name);
      reset_offsets();
      return *this;
  }

  DistVec& operator=(const std::vector<double> & pvec){
    assert(offset >= 0);
    assert(pvec.size() == m || pvec.size() == n);
    if (pvec.size() == m){
      transpose = true;
      assert(!square);
      start = n; end = m+n;
    }
    
    for (int i=start; i< (int)end; i++){  
         const vertex_data * data = &glcore->graph().vertex_data(i);
         double * pv = (double*)data;
         pv[offset] = pvec[i-start];  
    }
    debug_print(name);
    return *this;       
  }

  void to_vec(std::vector<double> & v){
     for (int i=start; i<end; i++){
         const vertex_data * data = &glcore->graph().vertex_data(i);        double * pv = (double*)data;
         v[i-start] = pv[offset];
     }
  }

  void debug_print(const char * name){
     if (debug){
       std::cout<<name<<" ("<<names[offset]<<" "<<offset<<" [ " << end-start << "] ";
       for (int i=start; i< std::min(end, start+MAX_PRINT_ITEMS); i++){  //TODO
         const vertex_data * data = &glcore->graph().vertex_data(i);
         double * pv = (double*)data;
         std::cout<<pv[r_offset==-1?offset:r_offset]<<" ";
       }
       std::cout<<std::endl;
     }
  }
  void debug_print(std::string name){ return debug_print(name.c_str());}

  DistDouble operator*(const DistVec & other);
  
  DistVec& operator*(const double val){
     d=val;
     return *this;
  }
  DistVec& _transpose() { 
     /*if (!square){
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
    int start; 
    int end;
    bool transpose;

    DistMat() { 
      start = 0; 
      end = n; //@TODO
      transpose = false;
    };


    DistMat &operator*(DistVec & v){
      	x_offset = v.offset;
        A_offset = true;
        //v.transpose = transpose;
        //r_offset = A_offset;
        return *this;
    }
    DistMat &operator-(){
        c=-1.0;
        return *this;
    }
    
    DistMat &operator+(){
        c=1.0;
        return *this;
    }
    DistMat &operator+(const DistVec &v){
        y_offset = v.offset;
        if (d == 0.0)
            d=1.0;
        return *this;
    }
    DistMat &operator-(const DistVec &v){
        y_offset = v.offset;
        if (d == 0.0)
          d=-1.0;
        else 
          d*=-1.0;
        return *this;
    }
    DistMat & _transpose(){
       transpose = true;
       if (!square){
         start=n ; end = m+n;
       }
       else {
         start =0; end = n;
       }
       return *this;
    }
   
};

DistVec& DistVec::operator=(DistMat &mat){
  r_offset = offset;
  transpose = mat.transpose;
  if (!transpose || square){
    start = 0; end = n;
  }
  else {
    start = n; end = m+n;
  }
  glcore->add_tasks((!mat.transpose)?rows:cols, Axb, 1);
  runtime += glcore->start();
  debug_print(name);
  reset_offsets();
  mat.transpose = false;
  return *this;
}


class DistDouble{
  public:
     double val;
     std::string name;
     DistDouble(){};
   
     const DistVec& operator*(const DistVec & dval){
        d=val;
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
      y_offset = offset;
      b_offset = vec.offset;
      if (d == 0) 
        d = 1.0;
      assert(y_offset >=0 && b_offset >= 0);

      double val = 0;
      start = vec.start; end = vec.end;
      for (int i=start; i<  end; i++){  
         const vertex_data * data = &glcore->graph().vertex_data(i);
         double * pv = (double*)data;
        // if (y_offset >= 0 && b_offset == -1)
	     //val += pv[y_offset] * pv[y_offset];
         val += d* pv[y_offset] * pv[b_offset];
      }
      reset_offsets();
      DistDouble mval;
      mval.val = val;
      return mval;
 
     }


int size(DistMat & A, int pos){
   assert(pos == 1 || pos == 2);
   if (square || pos == 1)
     return n;
   if (pos == 2){
     assert(m!= 0);
     return m;
   }
   return -1;
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
       const vertex_data * data = &glcore->graph().vertex_data(i);
       double * px = (double*)data;
       mval.val += px[vec.offset]*px[vec.offset];
    }
    mval.val = sqrt(mval.val);
    return mval;
}

edge_list get_edges(gl_types::iscope & scope){
     return (scope.vertex() < n) ? scope.out_edge_ids(): scope.in_edge_ids();
}
edge_list get_edges(graph_type *g, vertex_id_t i){
     return (i < n) ? g->out_edge_ids(i) : g->in_edge_ids(i);
}
const vertex_data& get_neighbor(gl_types::iscope & scope, edge_id_t oedgeid){
     return (scope.vertex() < n) ? scope.neighbor_vertex_data(scope.target(oedgeid)) :  scope.neighbor_vertex_data(scope.source(oedgeid));
}
const vertex_data& get_neighbor(graph_type *g, vertex_id_t i, edge_id_t oedgeid){
     return (i < n) ? g->vertex_data(g->target(oedgeid)) : g->vertex_data(g->source(oedgeid));
}
 
/***
 * UPDATE FUNCTION (ROWS)
 */
void Axb(gl_types::iscope &scope, 
	 gl_types::icallback &scheduler) {
    

  /* GET current vertex data */
  vertex_data& user = scope.vertex_data();
  double * pr = (double*)&user;
  assert(r_offset >=0);
  double val = 0;
  assert(x_offset >=0 || y_offset>=0);

  /*** COMPUTE r = c*A*x  ********/
  if (A_offset  && x_offset >= 0){
    edge_list outs = get_edges(scope);
    foreach(graphlab::edge_id_t oedgeid, outs) {
      edge_data & edge = scope.edge_data(oedgeid);
      const vertex_data  & movie = get_neighbor(scope, oedgeid);
      double * px = (double*)&movie;
      val += (c * edge.weight * px[x_offset]);
    }
  
    if (square)// add the diagonal term
      val += (c* user.prior_prec * pr[x_offset]);
  }
 /***** COMPUTE r = c*I*x  *****/
  else if (!A_offset && x_offset >= 0){
     val = c*pr[x_offset];
  }
  
  /**** COMPUTE r+= d*y (optional) ***/
  if (y_offset>= 0){
    val += d*pr[y_offset]; 
  }

  pr[r_offset] = val;
}

void fast_Axb(graph_type * g, std::vector<vertex_id_t> nodes){
   for (int j=0; j< (int)nodes.size(); j++){
       vertex_id_t i =  nodes[j];
       vertex_data & user = g->vertex_data(i);
       double * pr = (double*)&user;
       
       assert(r_offset >=0);
       double val = 0;
       assert(x_offset >=0 || y_offset>=0);

      /*** COMPUTE r = c*A*x  ********/
     if (A_offset  && x_offset >= 0){
     edge_list outs = get_edges(g, i);
    foreach(graphlab::edge_id_t oedgeid, outs) {
      edge_data & edge = g->edge_data(oedgeid);
      const vertex_data  & movie = get_neighbor(g, i, oedgeid);
      double * px = (double*)&movie;
      val += (c * edge.weight * px[x_offset]);
    }
  
    if (square)// add the diagonal term
      val += (c* user.prior_prec * pr[x_offset]);
    }
   /***** COMPUTE r = c*I*x  *****/
    else if (!A_offset && x_offset >= 0){
      val = c*pr[x_offset];
    }
  
  /**** COMPUTE r+= d*y (optional) ***/
    if (y_offset>= 0){
      val += d*pr[y_offset]; 
    }
     pr[r_offset] = val;
   }
}   

void init_row_cols(){
      
    for (int i=0; i< (int)n; i++)
      rows.push_back(i);

    if (square){
      cols = rows;
    }
    else {
      assert(m!= 0);
      for (int i=n; i < (int)(m+n); i++)
        cols.push_back(i);
    }
}

std::vector<double> ones(double val, int len){
    std::vector<double> ret;
    for (int i=0; i< len; i++)
      ret.push_back(val);
    return ret;
}


#include <graphlab/macros_undef.hpp>
#endif //_MATH_HPP
