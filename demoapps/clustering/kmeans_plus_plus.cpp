// parallel implementation of the k-means++ algorithm
// see: http://en.wikipedia.org/wiki/K-means%2B%2B
// written by Danny Bickson, CMU


#include "clustering.h"
#include "../gabp/advanced_config.h"

extern advanced_config ac;
extern problem_setup ps;

void add_tasks(gl_types::core &glcore);
void run_graphlab(gl_types::core & glcore);

using namespace std;

void initialize_clusters(gl_types::core &glcore){
   int first = -1;
   while(true){
     first = ::randi(0, ps.M-1);
     if (!ps.g<graph_type>()->vertex_data(first).reported){
        logstream(LOG_WARNING) << " node " << first << " has no edges - can not be selcted to cluster head " << endl;
     }
     else break;
   }
   if (ac.debug)
     cout<<"****Seed node for kmeans++ is: " << first << endl;
   cluster first_clust;
   assign(first_clust.location, ps.g<graph_type>()->vertex_data(first).datapoint, ps.N);
   first_clust.sum_sqr =sum_sqr(ps.g<graph_type>()->vertex_data(first).datapoint); 
   ps.g<graph_type>()->vertex_data(first).clusterhead = true;
   ps.g<graph_type>()->vertex_data(first).current_cluster = ac.K-1; //start from last to first
   first_clust.num_assigned_points = 1;
   assign(first_clust.cur_sum_of_points, ps.g<graph_type>()->vertex_data(first).datapoint, ps.N);
   std::vector<cluster>::iterator it = ps.clusts.cluster_vec.begin();
   ps.clusts.cluster_vec.insert(it,first_clust);

   flt_dbl_vec distances = zeros(ps.M+1);
   for (int i=0; i<ac.K-1; i++){
     glcore.start();
     add_tasks(glcore);

      int max_node = -1; double max_dist = 0;
    
      for (int j=0; j< ps.M; j++){
        const vertex_data & data = ps.g<graph_type>()->vertex_data(j);
        if (max_dist < data.min_distance){
           max_dist = data.min_distance;
           max_node = j;
        }
        if (data.clusterhead)
          distances[j+1] = 0;//do not use this point as clusterhead in this stage if this is already a clusterhead from previous stage
        else if (data.reported)              
          distances[j+1] = data.min_distance;
        else
	  distances[j+1] = 0;//if this matrix row is all zeros, don't take this node as cluster head
      }

      int thenode = -1; double thenode_dist;
      distances = pow(distances,2);
      assert(sum(distances) != 0);
      distances = distances / sum(distances); //normalize to one
      distances = cumsum(distances);
      double loc = randu();
      for (int j=0; j< ps.M; j++){
        if (loc >= distances[j] && loc < distances[j+1]){
          thenode = j; 
          thenode_dist = ps.g<graph_type>()->vertex_data(j).min_distance; 
          assert(!ps.g<graph_type>()->vertex_data(j).clusterhead);
          break;
        }
      }
      assert(thenode != -1);
      cluster cur_clust;
      assign(cur_clust.location,ps.g<graph_type>()->vertex_data(thenode).datapoint, ps.N);
      cur_clust.sum_sqr = sum_sqr(ps.g<graph_type>()->vertex_data(thenode).datapoint);
      cur_clust.num_assigned_points = 1;
      assign(cur_clust.cur_sum_of_points, ps.g<graph_type>()->vertex_data(thenode).datapoint, ps.N);
      ps.g<graph_type>()->vertex_data(thenode).clusterhead = true;
      ps.g<graph_type>()->vertex_data(thenode).current_cluster = ac.K-2-i;
      it = ps.clusts.cluster_vec.begin();
      ps.clusts.cluster_vec.insert(it,cur_clust);     
 
      if (ac.debug){
         cout<<"Found maximal distance: " << max_dist << " node: " << max_node << endl;
         cout<<"*****Selected node for cluster " << i+1 << " is: " << thenode << " with dist: " << thenode_dist << endl;
      }
   }
   
   logstream(LOG_INFO)<<"Finished with k-means++ initialization, going to run K-means with " <<ac.K << " clusters " << endl;

   assert((int)ps.clusts.cluster_vec.size() == ac.K);
  
   ac.algorithm = K_MEANS; ps.algorithm = K_MEANS; 
   run_graphlab(glcore);
}










