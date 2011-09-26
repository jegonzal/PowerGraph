#include "clustering.h"
#include "../gabp/advanced_config.h"

extern problem_setup ps;
extern advanced_config ac;

void print_point(FILE * f, vertex_data & data){
   fprintf(f, "\t[");
   bool first = true;
   for (int i=0; i< data.datapoint.nnz(); i++){
      if (!first){
	fprintf(f, " ");
      }
      else first = false;
      fprintf(f, "%d:%g", data.datapoint.get_nz_index(i), data.datapoint.get_nz_data(i));
   }
   fprintf(f, "]\n");
}



void print_points_for_cluster(FILE * f, int cluster_id){
  for (int i=0; i< ps.M; i++){
     vertex_data & data = ps.g<graph_type>()->vertex_data(i);
     if (data.current_cluster == cluster_id)
        print_point(f, data);

  }
}



void find_radius(double * radius){
   
  memset(radius, 0, sizeof(double)*ac.K);
  for (int i=0; i<ps.M; i++){
     vertex_data & data = ps.g<graph_type>()->vertex_data(i);
     if (radius[data.current_cluster] < data.min_distance)
       radius[data.current_cluster] = data.min_distance;
   }

}


void dumpcluster(){


   char outfile[256];
   sprintf(outfile, "%s%d.clusters.txt", ac.datafile.c_str(), ac.K);
   FILE * f = fopen(outfile, "w");
   logstream(LOG_INFO)<< "Dumping clusters info into file " << std::string(outfile) << std::endl;

   double * radius = new double[ac.K];
   find_radius(radius);

    for (int i=0; i< ac.K; i++){
      fprintf(f, "CL-%d { n=%d c=[", i, ps.clusts.cluster_vec[i].num_assigned_points);
       bool first = true;
        for (int j=0; j< ps.N; j++){
	 //print cluster location
         if (ps.clusts.cluster_vec[i].location[j] != 0){
             if (!first){
		fprintf(f, " ");
             }
	     else
                first = false;

             fprintf(f, "%d:%g", j, ps.clusts.cluster_vec[i].location[j]);
         }
      } 
     //print cluster radius
     fprintf(f, "] r=[%g]\n", radius[i]);

      //print all points belonging to this cluster
      print_points_for_cluster(f,i);

    }
      

    fclose(f);
    delete []radius;







}

