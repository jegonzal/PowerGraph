#ifndef _BLOCK_UTILS_H
#define _BLOCK_UTILS_H

#include <graphlab.hpp>

#include <graphlab/macros_def.hpp>

//#include <itpp/itbase.h>
//#include <itpp/stat/misc_stat.h>
//#include <itpp/base/smat.h>

#define GABP_PRIOR_PREC_OFFSET 0

using namespace graphlab;
using namespace itpp;

typedef Sparse_Mat<double> ssmat;
typedef Sparse_Vec<double> ssvec;



template<typename graph>
double  collect_diff(int start_pos, int end_pos, graph * g, int first, int second, int d){
    typedef typename graph::vertex_data_type vdata;

        assert(end_pos >= start_pos );
        assert(first>=0 && second >=0 && first != second);
        sdouble diff = 0;
        for (int i=start_pos; i < end_pos; i++){
           vdata * data = &g->vertex_data(i);
           for (int j=0; j<d; j++){
           	diff += pow((*data->data[first])(j) - (*data->data[second])(j),2); 
           }
        }
	return sqrt(diff);
     }

 template<typename graph> 
 edge_id_t **  parse_edge(edata * ed, int i, int d, graph * g, int num, edge_id_t **& edges){
        //printf("adding an edge from %d to %d val %e\n", ed[i].to-1, ed[i].from-1, ed[i].weight);
       typedef typename graph::edge_data_type gabp_block_edge_data;
       typedef typename graph::vertex_data_type gabp_block_vertex_data;
       typedef std::pair<bool, edge_id_t> edgeflag;
       
         int pto =( ed[i].to-1) / d;
        int pfrom = (ed[i].from-1) /d;
       
	if (edges==NULL){
        	edges = new edge_id_t*[num];
	for (int t=0; t<num; t++){
           edges[t] = new edge_id_t[num];
           for (int j=0; j< num; j++)
             edges[t][j] = -1;
         }
        }

        edgeflag ef;
        const gabp_block_edge_data * edge = NULL;
         if (pto == pfrom){
		 gabp_block_vertex_data * data = &g->vertex_data(pto);
	 	 assert(data);
         	 data->pmat[GABP_PRIOR_PREC_OFFSET]->set((ed[i].from-1)%d,(ed[i].to-1)%d, ed[i].weight);
         	 data->pmat[GABP_PRIOR_PREC_OFFSET]->set((ed[i].to-1)%d,(ed[i].from-1)%d, ed[i].weight);
         }
         else {
	        edge_id_t e = edges[pto][ pfrom];
                if (e==0xFFFFFFFF){
         	edge = new gabp_block_edge_data();
                edge->weight->set((ed[i].from-1) % d, (ed[i].to-1) % d, ed[i].weight);
                assert(pfrom != pto);
                edges[pfrom][ pto] = g->add_edge(pfrom, pto, *edge); // Matlab export has ids starting from 1, ours start from 0
                
        	//printf("adding an edge %d %d \n", pfrom, pto);
                edge = new gabp_block_edge_data();
                edge->weight->set((ed[i].to-1) % d, (ed[i].from-1) % d, ed[i].weight);
                edges[pto][ pfrom] = g->add_edge(pto, pfrom, *edge); // Matlab export has ids starting from 1, ours start from 0
                
                // assert((*(edge->prec))(0,0) == 0);                
        }
        else {
                assert(edges[pto][pfrom] != 0xFFFFFFFF);
                assert(edges[pfrom][pto] != 0xFFFFFFFF);
                edge = &g->edge_data(edges[pfrom][pto]);
                edge->weight->set((ed[i].from-1)%d, (ed[i].to-1) % d, ed[i].weight);
                edge = &g->edge_data(edges[pto][pfrom]);
		edge->weight->set((ed[i].to-1)%d, (ed[i].from-1) % d , ed[i].weight);
        }
        }
        return edges;
    }

 template<typename graph> 
 edge_id_t **  read_block_edges(FILE * f, int len, int offset, int nodes, graph * g, bool symmetry, int d, ssmat * sm, int M, bool doparse, int num){
     assert(offset>=0 && offset < len);
     

     unsigned int e;
     int rc = fread(&e,1,4,f);
     assert(rc == 4);
     printf("Creating %d edges...\n", e);
     assert(e>0);
     int total = 0;
     edata* ed = new edata[200000];
     printf("symmetry: %d\n", symmetry);
     int edgecount_in_file = e;
     if (symmetry) edgecount_in_file /= 2;
     int ec = 0;   
 
     edge_id_t ** edges = NULL;

      while(true){
       memset(ed, 0, 200000*sizeof(edata));
       rc = (int)fread(ed, sizeof(edata), _min(200000, edgecount_in_file - total), f);
       total += rc;

       for (int i=0; i<rc; i++){
        if (symmetry && ec >= (int)(e/2))
             break;

    	assert(ed[i].weight != 0);
    	if (M == 0){//adj matrix
	  assert(ed[i].from >= 1 && ed[i].from <= nodes);
    	  assert(ed[i].to >=1 && ed[i].to <= nodes);
         }
        else {//a bipartite graph
          if (symmetry && (ec % 2) == 1)
		continue;

          assert(ed[i].from >= M && ed[i].from <= nodes);
    	  assert(ed[i].to >=1 && ed[i].to <= M);
        }
    	assert(ed[i].to != ed[i].from);
       
        if (sm){
	     sm->set(ed[i].to-1, ed[i].from-M-1,ed[i].weight);

        }
        if (doparse)
		parse_edge(ed, i, d, g, num, edges);
       ec++;
	}   
        printf(".");
        fflush(0);
        if (rc == 0 || total >= edgecount_in_file)
	   break;
    }
    delete [] ed; ed = NULL;
    return edges;
  }

 void debug_print_matrix(const char * name, mat** pmat, int m, int n){
      for (int i=0; i<m; i++)
         for(int j=0; j<n; j++)
            std::cout<<i<<":"<<j<<" " << pmat[i][j] << std::endl;
  }
 

 template <typename graph> 
 void debug_print_graph(const char * name, graph * g, int start_pos, int end_pos, int offset, int m, int n, int d){
       typedef typename graph::vertex_data_type gabp_block_vertex_data;   
       typedef typename graph::edge_data_type gabp_block_edge_data;   
 
        assert(start_pos >= 0 && start_pos < end_pos && offset>=0);
        assert(m>0 && n>0);
        mat ** pmat = new mat*[m];
      	for (int i=0; i < m; i++){
            pmat[i] = new mat[n];
            for (int j=0; j<n; j++){
                pmat[i][j] = zeros(d,d);
            }
           
            gabp_block_vertex_data * data = &g->vertex_data(i);
            pmat[i][i] = *data->pmat[GABP_PRIOR_PREC_OFFSET];
    	    foreach(edge_id_t eid, g->out_edge_ids(i)){          
                gabp_block_edge_data * data= &g->edge_data(eid);
	        int to = g->target(eid);
                int from = g->source(eid);
                pmat[from][to] = (*data->weight);
          }
        }  
        debug_print_matrix(name, pmat,m,n);
        for (int i=0; i<m; i++)
		delete [] pmat[i];
        delete [] pmat;
  } 

template<typename graph> 
  void dispatch_vecb(int start_pos, int end_pos, int offset, graph * g, double * _vec, bool free_vec, int d ){
        
      typedef typename graph::vertex_data_type gabp_block_vertex_data;
       assert(end_pos >=start_pos);
        assert(_vec != NULL);
	for (int i=start_pos; i < end_pos; i++){
              gabp_block_vertex_data * data = &g->vertex_data(i);
              vec * pvec = data->data[offset];
              assert(pvec->size() == d);
              for (int j=0; j<d; j++){
                   (*pvec)[j] = _vec[d*(i-start_pos)+j];
	   }
	}
        if (free_vec)
           delete [] _vec;
     }
  template<typename graph> 
 void dispatch_vecb(int start_pos, int end_pos, int offset, graph * g, double _vec, int d ){
      typedef typename graph::vertex_data_type gabp_block_vertex_data;
      
       assert(end_pos >=start_pos);
	for (int i=start_pos; i < end_pos; i++){
              gabp_block_vertex_data * data = &g.vertex_data(i);
              vec * pvec = data->data[offset];
              assert(pvec->size() == d);
              for (int j=0; j<d; j++){
                   (*pvec)[j] = _vec;
	   }
	}
 }


  template<typename graph> 
  void dispatch_mat(int start_pos, int end_pos, int offset, graph * g, double * _vec, bool free_vec, int d ){
        typedef typename graph::vertex_data_type gabp_block_vertex_data;

        assert(end_pos >= start_pos);
	for (int i=start_pos; i < end_pos; i++){
              gabp_block_vertex_data * data = &g->vertex_data(i);
              mat * pmat = data->pmat[offset];
              //assert(pvec->size() == D);
              for (int j=0; j<d; j++){
                   pmat->set(j,j, _vec[d*(i-start_pos)+j]);
	   }
	}
        if (free_vec)
           delete [] _vec;
     }


  template<typename graph> 
double* collect_vec(int start_pos, int end_pos, int offset, graph * g, int d ){
        typedef typename graph::vertex_data_type gabp_block_vertex_data;
        
        assert(end_pos >= start_pos );
        assert(d>0);
        double * ret = new double[(end_pos-start_pos)*d];
	for (int i=start_pos; i < end_pos; i++){
              gabp_block_vertex_data * data = &g->vertex_data(i);
              vec * pvec = data->data[offset];
              assert(pvec->size() == d);
              for (int j=0; j<d; j++){
                   ret[d*(i-start_pos)+j] = (*pvec)[j];
	   }
	}
       return ret;   
  }



  template<typename graph> 
  void read_nodesb(FILE * f, int len, int nodes, graph * g,int offset, int d){
  
    typedef typename graph::vertex_data_type gabp_block_vertex_data;

    assert(nodes>0);

    int toread = nodes;
    if (nodes > BUFSIZE)
       toread = BUFSIZE;
    int remain = nodes;

    double * temp = new double[remain*d];
    while(remain > 0){
        toread = (remain < BUFSIZE)?remain:toread;
        int rc = (int)fread(temp, d*sizeof(double), toread, f);
    	//assert(rc == toread);
        remain -= rc;	
        
	for (int i=0; i< rc; i++){
          gabp_block_vertex_data * data = new gabp_block_vertex_data();
           for (int j=0; j< d; j++){
	   	(*data->data[offset])[j]= temp[i*d+j];
           }
	   g->add_vertex(*data);
        }
   }

  
   delete [] temp;
  
  }




#include <graphlab/macros_undef.hpp>
#endif //_BLOCK_UTILS_H
