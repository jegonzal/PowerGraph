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
int x_offset = -1, b_offset = -1, y_offset = -1, r_offset = -1, A_offset = -1;
int increment=0;
gl_types::core * glcore;
double runtime = 0;

void reset_offsets(){
  c=1.0; d=0.0;
  x_offset = b_offset = y_offset = r_offset = A_offset = -1;
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
void ATxb(gl_types::iscope &scope, gl_types::icallback &scheduler);
 

class DistMat;
class DistDouble;

class DistVec{
   public:
   int offset;

   DistVec(){ 
     offset = increment_offset();
   };

   /*DistVec& operator*(DistVec & v){
     x_offset = v.offset;
     return *this;
   }*/
   DistVec& operator-(){
     d=-1.0;
     return *this; 
   }
   DistVec& operator-(const DistVec & other){
     x_offset = offset;
     y_offset = other.offset;
     d=-1.0;
     return *this;
   }
   DistVec& operator+(){
     d=1.0;
     return *this;
   }
   DistVec& operator+(const DistVec &other){
      x_offset =offset;
      y_offset = other.offset;
      return *this; //TODO
   }

   /*DistVec& operator=(const DistVec & vec){
    output = itpp::zeros(end-start);
    for (int i=start; i< end; i++){ 
       const vertex_data * data = &glcore->graph().vertex_data(i);
       output[i-start] = ((double*)data)[offset];
     }
    return *this;  
  }*/

   DistVec& operator=(const DistVec & vec){
        //if ( == 0)
          glcore->add_tasks(rows, Axb, 1); //TODO
        //else if (mat.start == (int)m)
	//  glcore->add_tasks(cols, ATxb, 1);
        //else assert(false);
        runtime += glcore->start();
        reset_offsets();
      return *this;
  }

  DistDouble operator*(const DistVec & other);
  
  DistVec& operator*(const double val){
     d=val;
     return *this;
  }
  DistVec& transpose() { 
    return *this;
  }

  DistVec& operator=(const DistMat &mat);
  //DistDouble operator=(const DistVec & vec);
 

 };



/*
 * wrapper for computing r = c*A*x+d*b*y
 */
class DistMat{
  public:
    int start; 
    int end;
 
    DistMat() { 
      start = 0; 
      end = m+n; //@TODO
    };


    DistMat &operator*(const DistVec & v){
	x_offset = v.offset;
        A_offset = increment_offset();
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
        d=1.0;
        return *this;
    }
    DistMat &operator-(const DistVec &v){
        d=-1.0;
        return *this;
    }
     
};

  DistVec& DistVec::operator=(const DistMat &mat){
        if (mat.start == 0)
          glcore->add_tasks(rows, Axb, 1);
        else if (mat.start == (int)m)
	  glcore->add_tasks(cols, ATxb, 1);
        else assert(false);
        runtime += glcore->start();
        reset_offsets();
        return *this;
  }


class DistDouble{
  public:
     double val;
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
     DistDouble operator=(const DistDouble & other){
         DistDouble mval;
         mval.val = other.val;
         return mval;
     }

 };

 DistDouble DistVec::operator*(const DistVec & vec){
      double val = 0;
      for (int i=0; i< (int)(m+n); i++){  //TODO
         const vertex_data * data = &glcore->graph().vertex_data(i);
         double * pv = (double*)data;
         if (y_offset >= 0 && b_offset == -1)
	   val += pv[y_offset] * pv[y_offset];
         else if (y_offset >=0 && b_offset >= 0)
         val += pv[y_offset] * pv[b_offset];
      }
      reset_offsets();
      DistDouble mval;
      mval.val = val;
      return mval;
 
     }


int size(DistMat & A, int pos){
   assert(pos == 1 || pos == 2);
   return pos==1 ? m+n: n; //TODO
}
DistDouble sqrt(DistDouble & dval){
    DistDouble mval;
    mval.val=sqrt(dval.val);
    return mval;
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
  pr[r_offset] = 0;

  assert(x_offset >=0 || y_offset>=0);

  /*** COMPUTE r = c*A*x  ********/
  if (A_offset >= 0 && x_offset >= 0){
    edge_list outs = scope.out_edge_ids();
    foreach(graphlab::edge_id_t oedgeid, outs) {
      edge_data & edge = scope.edge_data(oedgeid);
      vertex_data  & movie = scope.neighbor_vertex_data(scope.target(oedgeid));
      double * px = (double*)&movie;
      pr[r_offset] += (c * edge.weight * px[x_offset]);
    }
  }
 /***** COMPUTE r = c*I*x  *****/
  else if (A_offset== -1 && x_offset >= 0 && y_offset == -1){
     pr[r_offset] = d*pr[x_offset];
  }
  
  /**** COMPUTE r+= d*y (optional) ***/
  if (y_offset>= 0){
    pr[r_offset] += d*pr[y_offset]; 
  }
}
/***
 * UPDATE FUNCTION (COLS)
 */
void ATxb(gl_types::iscope &scope, 
	 gl_types::icallback &scheduler) {
    


  /* GET current vertex data */
  vertex_data& user = scope.vertex_data();
  double * pr = (double*)&user;
  assert(r_offset>= 0);
  pr[r_offset] = 0;

  edge_list ins = scope.in_edge_ids();
  timer t;
  t.start(); 

   foreach(graphlab::edge_id_t iedgeid, ins) {
      edge_data & edge = scope.edge_data(iedgeid);
      vertex_data  & movie = scope.neighbor_vertex_data(scope.source(iedgeid));
      double * px = (double*)&movie;
      pr[r_offset] += edge.weight * px[x_offset];
   }

   pr[r_offset] += d*pr[y_offset];

}



double cg(gl_types::core * _glcore){

    glcore = _glcore;

    for (int i=0; i< (int)(m==0?m+n:n); i++)
      rows.push_back(i);
 
    for (int i=m; i< (int)(m+n); i++)
      cols.push_back(i);
 
    DistMat A;
    DistVec b, r, p, x, Ap;
    DistDouble rsold, rnew, alpha;
  

    /* r = -A*x+b;
       p = r;
       rsold = r'*r;
    */

    r=-A*x+b; //r = -A*x+b; @TODO
    p = r;
    rsold = r.transpose()*r;

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

    for (int i=1; i < size(A,1); i++){
        Ap=A*p;
        alpha=rsold/(p.transpose()*Ap);
        x=x+alpha*p;
        r=r-alpha*Ap;
        rnew=r.transpose()*r;
        if (sqrt(rnew)<1e-10)
            break;
        p=r+(rnew/rsold)*p;
        rsold=rnew;
    }

    return runtime;

}
#include <graphlab/macros_undef.hpp>
#endif
