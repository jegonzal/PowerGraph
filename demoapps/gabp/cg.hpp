#ifndef _CG_HPP
#define _CG_HPP

#include <itpp/itbase.h>
#include <itpp/itstat.h>
#include <itpp/stat/misc_stat.h>
#include <graphlab.hpp>
#include <graphlab/macros_def.hpp>

extern uint32_t m,n;
double  c=1.0;
double  d=0.0;
int x_offset = -1, b_offset = -1, y_offset = -1, r_offset = -1;
bool A_offset = false;
int increment=2;
gl_types::core * glcore;
double runtime = 0;
extern bool debug;
extern bool square;
extern int cg_maxiter;
#define MAX_PRINT_ITEMS 21

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
using namespace graphlab;

/*function [x] = conjgrad(A,b,x)
    r=b-A*x; ///DIST
    p=r;     //SER
    rsold=r'*r;  //SER

    for i=1:size(A,1)
        Ap=A*p;               //DIST
        alpha=rsold/(p'*Ap);  //SER
        x=x+alpha*p;          //SER
        r=r-alpha*Ap;         //SER
        rsnew=r'*r;           //SER
        if sqrt(rsnew)<1e-10  //SER
              break;
        end
        p=r+rsnew/rsold*p;    //SER
        rsold=rsnew;          //SER
    end
end
*/
std::vector<vertex_id_t> rows,cols;
void Axb(gl_types::iscope &scope, gl_types::icallback &scheduler);
 

class DistMat;
class DistDouble;



class DistVec{
   public:
   int offset;
   std::string name; //optional
   bool transpose;
   int start;
   int end;

   DistVec(){ 
     offset = increment_offset();
     transpose = false;
     start = 0;
     end = n;
     debug_print(name);
   };

   DistVec(int _offset){
     offset = _offset;
     start = 0;
     end = n;
     debug_print(name);
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
       if (x_offset == -1 && y_offset == -1){
         y_offset = vec.offset;
       }  
      r_offset = offset;
      if (d == 0.0)
        d=1.0;
      transpose = vec.transpose;
      start = vec.start;
      end = vec.end;
      glcore->add_tasks((!transpose)?rows:cols, Axb, 1); 
      
      runtime += glcore->start();
      debug_print(name);
      reset_offsets();
      return *this;
  }

  DistVec& operator=(const itpp::vec & pvec){
    assert(offset >= 0);
    assert(pvec.size() == (int)m || pvec.size() == (int)n);
    if (pvec.size() == (int)m){
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


  void debug_print(const char * name){
     if (debug){
       std::cout<<name<<" ("<<offset<<" [ " << end-start << "] ";
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
        v.transpose = transpose;
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
   if (pos == 2)
     return m;
}
DistDouble sqrt(DistDouble & dval){
    DistDouble mval;
    mval.val=sqrt(dval.val);
    return mval;
}

DistDouble norm(DistVec & vec){
    assert(vec.offset>=0);
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
     return ((int)scope.vertex() < (int)n) ? scope.out_edge_ids(): scope.in_edge_ids();
}
const vertex_data& get_neighbor(gl_types::iscope & scope, edge_id_t oedgeid){
     return ((int)scope.vertex() < (int)n) ? scope.neighbor_vertex_data(scope.target(oedgeid)) :  scope.neighbor_vertex_data(scope.source(oedgeid));
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



double cg(gl_types::core * _glcore){

    glcore = _glcore;
    init_row_cols();

    DistMat A;
    DistVec b(0), r, p, x, Ap, t;
    b.transpose=true;
    b.start=n; b.end=m+n;
    if (!square)
      x = itpp::ones(m)/2;
    else
      x = itpp::ones(n)/2;
    DistDouble rsold, rnew, alpha, tmpdiv;
 

    /* r = -A*x+b;
       p = r;
       rsold = r'*r;
    */
    if (square){
      r=-A*x+b; 
    }
    else {
      r=-A*x;
      r=A._transpose()*r;
      r=r+b;
    }
    p = r;
    rsold = r._transpose()*r;

     /*
     for i=1:size(A,1)
        Ap=A*p;               
        alpha=rsold/(p'*Ap); 
        x=x+alpha*p;        
        r=r-alpha*Ap;      
        rsnew=r'*r;       
        if sqrt(rsnew)<1e-10 
              break;
        end
        p=r+rsnew/rsold*p;  
        rsold=rsnew;       
    end
    */

    for (int i=1; i <= std::min(cg_maxiter,size(A,1)); i++){
        Ap=A*p;
        if (!square)
          Ap= A._transpose()*Ap;
        tmpdiv = p._transpose()*Ap;
        alpha=rsold/tmpdiv;
        x=x+alpha*p;
    
        t=A*x;
        if (!square)
           t=A._transpose()*t-b;
        else
          t=t-b;
        logstream(LOG_INFO)<<"Iteration " << i << " approximated solution redidual is " << norm(t).toDouble() << std::endl;
        
        r=r-alpha*Ap;
        rnew=r._transpose()*r;
        if (sqrt(rnew)<1e-10){
          logstream(LOG_INFO)<<" Conjugate gradient converged in iteration "<<i<<" to an accuracy of "  << sqrt(rnew).toDouble() << std::endl; 
          break;
        }
        tmpdiv = rnew/rsold;
        p=r+tmpdiv*p;
        rsold=rnew;
    }

    return runtime;

}
#include <graphlab/macros_undef.hpp>
#endif
