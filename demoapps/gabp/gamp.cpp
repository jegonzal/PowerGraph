
/*
   Copyright [2011] [Aapo Kyrola, Joseph Bradley, Danny Bickson, Carlos Guestrin / Carnegie Mellon University]

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF Aps.m KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.


    Implementation of the GAMP algorithm.
    Based on Matlab code by Phil Schniter, Ohio State University
    Adapted to GraphLab by Danny Bickson, CMU
*/

 
//#include "common.h"
#include "linear.h"
#include "cas_array.h"
#include "advanced_config.h"
#include <assert.h>

// Parameters
extern problem_setup ps;
extern advanced_config config;
graph_type_gamp *g_gamp;




extern int offset, offset2, offset3;

struct gamp_struct{
   //observation node
   double * Z_var, * Z_mean, * P_var, *P_mean, *S_var, *S_mean, * A_i, *Y;
   //variable node
   double * R_var, * R_mean, * gam, * nu, * pi, * X_mean, *X_var, * AT_i, *X_true;
  
   unsigned int id;

   gamp_struct(vertex_data_gamp &vdata, unsigned int _id){
      assert(config.D > 0);

      id = _id;
      double* pdata = (double*)data(vdata.data);
      assert(pdata);
      //observation nodes
      if (id >= ps.m){
         Z_var = pdata;
         Z_mean = pdata+config.D;
         P_var = pdata+2*config.D;
         P_mean = pdata+3*config.D;
         S_var = pdata+4*config.D;
         S_mean = pdata+5*config.D;
	 Y = pdata+6*config.D;
	 A_i = pdata+7*config.D;

	 R_var = R_mean = gam = nu = pi = X_mean = X_var = AT_i = X_true = NULL;
      } //iable nodes
      else {
         R_var = pdata;
	 R_mean = pdata+config.D;
         gam = pdata+2*config.D;
         nu = pdata+3*config.D;
         pi = pdata+4*config.D;
         X_mean = pdata+5*config.D;
         X_var = pdata+6*config.D;
         X_true = pdata+7*config.D;

 	 AT_i = pdata+8*config.D;

 	 Z_var = Z_mean = P_var = P_mean = S_var = S_mean  = Y = A_i = NULL;
      }
   }

   vec get_X(){
      return init_vec(X_mean, config.D);
   }

   void print(){
      if (id >= ps.m) { //obs
         printf("Printg obs node %d\n", id);
         debug_print_vec("Z_var", Z_var, config.D);
         debug_print_vec("Z_mean", Z_mean, config.D);
         debug_print_vec("P_var", P_var, config.D);
         debug_print_vec("P_mean", P_mean, config.D);
         debug_print_vec("S_var", S_var, config.D);
         debug_print_vec("S_mean", S_mean, config.D);
         debug_print_vec("Y", Y, config.D);
	 debug_print_vec("A_i", A_i, ps.m);
      } 
      else {
         printf("Printg var node %d\n", id);
         debug_print_vec("R_var", R_var, config.D);
         debug_print_vec("R_mean", R_mean, config.D);
         debug_print_vec("gam", gam, config.D);
         debug_print_vec("nu", nu, config.D);
         debug_print_vec("pi", pi, config.D);
         debug_print_vec("X_mean", X_mean, config.D);
         debug_print_vec("X_var", X_var, config.D);
	 debug_print_vec("X_true", X_true, config.D);
         debug_print_vec("AT_i", AT_i, ps.n);
      }
   }

};

int calc_size(unsigned int id){
 if (id < ps.m)
   return ps.n+config.D * 8;
  else return ps.m+config.D*9;
}


/**
 *
 * [n,k] = size(A);
   V = zeros(k,m+1);
   V(:,2) = rand(k,1);
   V(:,2)=V(:,2)/norm(V(:,2),2);
   beta(2)=0;
 *
 * */
/*
void init_gamp(){
   int m = config.iter;
   double sum = 0;

  const graph_type_gamp *g = ps.g<graph_type_gamp>(TRAINING);

  for (int i=ps.m; i< ps.m+ps.n; i++){ 
    vertex_data_gamp * data = (vertex_data_gamp*)&g->vertex_data(i); 
     gamp_struct obs(*data);
 }

  for (int i=0; i< ps.m; i++){ 
    vertex_data * data = (vertex_data*)&g->vertex_data(i);
    gamp_struct (*data);
    
  }

}*/

/***
 * UPDATE FUNCTION (ROWS)
 */
void gamp_obs_update_function(gl_types_gamp::iscope &scope, 
	 gl_types_gamp::icallback &scheduler) {
    


  /* GET current vertex data */
  vertex_data_gamp& vdata = scope.vertex_data();
  int id = scope.vertex(); 
  gamp_struct obs(vdata, id);
  /* print statistics */
  if (config.debug && (id == (int)ps.m || (id == (int)(ps.m+ps.n-1)))){ 
    printf("gamp: entering observation node  %u \n",  id);   
  }

  graphlab::timer t;
  t.start(); 
/*
 * observation nodes
 */
 
   memset(obs.Z_var, 0, sizeof(sdouble)*config.D);
   memset(obs.Z_mean, 0, sizeof(sdouble)*config.D);
   memset(obs.P_var, 0, sizeof(sdouble)*config.D);
   for (unsigned int k=0; k< ps.m; k++){
      vertex_data_gamp  & other= g_gamp->vertex_data(k);
      gamp_struct var(other, k);

      for (int i=0; i< config.D; i++){
        // Z_var = A_mean2*X_var;
        obs.Z_var[i] += pow(obs.A_i[k],2) * var.X_var[i];
        // Z_mean = A_mean*X_mean;
        obs.Z_mean[i] += obs.A_i[k] * var.X_mean[i];
        // P_var = Z_var + A_var*(X_var + X_mean.^2);
        obs.P_var[i] +=  pow(fabs(obs.A_i[k]),2) * (var.X_var[i] + pow(var.X_mean[i],2)) ;
      }
  }


  for (int i=0; i< config.D; i++){
     // P_var = Z_var + A_var*(X_var + X_mean.^2);
     obs.P_var[i] += obs.Z_var[i];
     // P_mean = Z_mean - Z_var.*S_mean;
     obs.P_mean[i] = obs.Z_mean[i] - obs.Z_var[i] * obs.S_mean[i];
     // S_ = 1./(P_var+psi);
     obs.S_var[i] = 1.0 / (obs.P_var[i] + ps.gamp_psi);
     // S_mean = S_var.*(Y-P_mean);
     obs.S_mean[i] = obs.S_var[i] *(obs.Y[i] - obs.P_mean[i]);
  }
 
  ps.counter[GAMP_MULT_A] += t.current_time();

  if (config.debug && (id == (int)ps.m || id == (int)(ps.m+ps.n-1))){
    printf("gamp: computed obs node  %u\n", id);   
    obs.print(); 
 }


}



/***
 * UPDATE FUNCTION (COLS)
 */
void gamp__update_function(gl_types_gamp::iscope &scope, 
	 gl_types_gamp::icallback &scheduler) {
    


  /* GET current vertex data */
  vertex_data_gamp& vdata = scope.vertex_data();
  unsigned int id = scope.vertex(); 
  gamp_struct var(vdata, id);
  /* print statistics */
  if (config.debug && (id == ps.m || (id == ps.m+ps.n-1))){ 
    printf("gamp: entering  var node  %u \n",  id);   
  }

  graphlab::timer t;
  t.start(); 


/*
 *  % computations at variable nodes
*/
   memset(var.R_var, 0, config.D*sizeof(sdouble)); 
   memset(var.R_mean, 0, config.D*sizeof(sdouble));
   for (unsigned int k=ps.m; k< ps.m+ps.n; k++){
      vertex_data_gamp  & other= g_gamp->vertex_data(k);
      gamp_struct obs(other,k);

      for (int i=0; i< config.D; i++){
        var.R_var[i] += pow(var.AT_i[k-ps.m],2) * obs.S_var[i];
        var.R_mean[i] += var.AT_i[k-ps.m] * obs.S_mean[i];
      }
  }


  for (int i=0; i< config.D; i++){
     //R_var = 1./(A_
     //mean2'*S_var + eps);
     var.R_var[i] = 1.0 / (var.R_var[i] + ps.gamp_eps);
     //R_mean = X_mean + R_var.*(A_mean'*S_mean);
     var.R_mean[i] = var.X_mean[i] + var.R_var[i] * var.R_mean[i];
     //gam = (R_.var*theta + phi.*R_mean)./(R_var+phi);
     var.gam[i] = (var.R_var[i]*ps.gamp_theta + ps.gamp_phi*var.R_mean[i]) / (var.R_var[i] + ps.gamp_phi);
     //nu = phi.*R_./(R_var+phi);
     var.nu[i] = ps.gamp_phi*var.R_var[i] / (var.R_var[i] + ps.gamp_phi);
     //pi = 1./(1+(1-lam)./lam .*sqrt((R_+phi)./R_var) .* exp( ...
     var.pi[i] = 1.0 / (1.0 + (1.0 - ps.gamp_lam) / ps.gamp_lam * sqrt((var.R_var[i]+ps.gamp_phi)/var.R_var[i]) * exp( \
      //  -(R_mean.^2)./(2*R_) + ((R_mean-theta).^2)./(2*(phi+R_var)) ));
        -(powf(var.R_mean[i],2)) / (2*var.R_var[i]) + (powf(var.R_mean[i] - ps.gamp_theta,2)) / (2*(ps.gamp_phi+var.R_var[i])) ));
     //X_mean = pi.*gam;
     var.X_mean[i] = var.pi[i] * var.gam[i];
      //X_ = pi.*(nu+gam.^2) - (pi.*gam).^2;
     var.X_var[i] = var.pi[i] * (var.nu[i]+ pow(var.gam[i],2)) - pow(var.pi[i] * var.gam[i],2);
  }
 
  ps.counter[GAMP_MULT_AT] += t.current_time();

  if (config.debug && (id == 0 || id == ps.m-1)){
    var.print();
  }


}


void calc_normalized_mse(){

  double diff_norm = 0;
  double X_norm = 0;


  for (unsigned int i=0; i< ps.m; i++){
    vertex_data_gamp & vdata = g_gamp->vertex_data(i);    
    gamp_struct var(vdata,i);
    for (int j=0; j< config.D; j++){
      diff_norm += pow(var.X_true[j] - var.X_mean[j],2);
      X_norm += pow(var.X_mean[j],2);
    }
  }
  diff_norm = sqrt(diff_norm);
  X_norm = sqrt(X_norm);
  if (config.debug)
     std::cout<<"Diff norm is: " << diff_norm << " X_NORM IS: " << X_norm << std::endl;
  logstream(LOG_INFO) << "MMSE (db) is: " << 20*log10(diff_norm / X_norm) << " X_norm is: " << X_norm << std::endl;
}

void add_vertices_gamp(mat &A, mat &X_true, mat &X_mean, mat &X_var, mat &Y, graph_type_gamp * g){
   for (unsigned int i=0; i<ps.m+ps.n; i++){
       vertex_data_gamp node;
       node.data=zeros(calc_size(i)); 
       if (i < ps.m){
          gamp_struct var(node, i);
          memcpy(var.X_true, data(X_true.get_row(i)), config.D*sizeof(double)); 
          memcpy(var.X_mean, data(X_mean.get_row(i)), config.D*sizeof(double)); 
          memcpy(var.X_var, data(X_var.get_row(i)), config.D*sizeof(double)); 
	  memcpy(var.AT_i, data(A.get_row(i)), ps.n*sizeof(double));
          if (config.debug && (i == 0 || i == ps.m-1))
		var.print();
      }
       else {
          gamp_struct obs(node, i);
          memcpy(obs.Y, data(Y.get_row(i-ps.m)), config.D*sizeof(double)); 
	  memcpy(obs.A_i, data(A.get_col(i-ps.m)), ps.m*sizeof(double));
          if (config.debug && (i == ps.m || i == ps.m+ps.n-1))
                obs.print();
       }
       g->add_vertex(node);
    }
}

void load_data_gamp(graph_type_gamp * g){
    mat A,Y, X_true, X_mean, X_var;
    vec values;
    it_file input(config.datafile);
    input >> Name("A_mean");
    input >> A;
    input >> Name("Y");
    input >> Y;
    input >> Name("X_true");
    input >> X_true;
    //input >> Name("X_mean");
    //input >> X_mean;
  
    input >> Name("X_var");
    input >> X_var; 
    input >> Name("values");
    input >> values;

    assert(values.size() == 8);
    ps.gamp_phi = values[0];
    ps.gamp_lam = values[1];
    ps.gamp_theta = values[2];
    ps.gamp_psi = values[3];
    ps.gamp_eps = values[4];
    ps.m = (uint)values[6];
    ps.n = (uint)values[5];
    config.D = (int)values[7];
    X_mean = zeros(ps.m,config.D);

    assert(A.rows() == (int)ps.m && A.cols() == (int)ps.n);
    assert(X_true.rows() == (int)ps.m && X_true.cols() == (int)config.D);
    assert(X_mean.rows() == (int)ps.m && X_true.cols() == (int)config.D);
    assert(X_var.rows() == (int)ps.m && X_true.cols() == (int)config.D);
    assert(Y.rows() == (int)ps.n && X_true.cols() == (int)config.D);

    add_vertices_gamp(A, X_true, X_mean, X_var, Y, g);
}

mat fill_output_gamp(graph_type_gamp * g){
  mat out(ps.m, config.D);
  for (uint i=0; i<ps.m; i++){
     vertex_data_gamp &varnode = g->vertex_data(i);
     gamp_struct var(varnode, i);
     set_row(out, i, var.get_X()); 
  }
  return out;
}

void write_vec(FILE * f, int len, const double * parray);

void write_output_gamp(graph_type_gamp *g, mat & X){

   
   it_file output((config.datafile+".out").c_str());
   std::cout<<"Writing result to file: "<<config.datafile<<".out"<<std::endl;
   output<<X;
   output.close();
}



void gamp_main(gl_types_gamp::core & glcore){
  
  g_gamp=&glcore.graph();
 
  std::vector<graphlab::vertex_id_t> rows,cols;
   for (unsigned int i=0; i< ps.m; i++)
      rows.push_back(i);
   for (unsigned int i=ps.m; i< ps.m+ps.n; i++)
      cols.push_back(i);
 
   //for j=2:m+2
   for (int j=1; j<= config.iter+1; j++){
        //w = A*V(:,j) 
	glcore.add_tasks(cols, gamp_obs_update_function,1);
        glcore.start();
        

        glcore.add_tasks(rows, gamp__update_function, 1);
	glcore.start();

        logstream(LOG_INFO) << "GAMP: Finished iteration " << j << " in time: " << ps.gt.current_time() << std::endl;
   	calc_normalized_mse();
   } 



 

}


