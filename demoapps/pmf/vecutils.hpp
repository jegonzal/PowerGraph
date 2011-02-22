#ifndef _VECUTILS_HPP
#define _VECUTILS_HPP

#include <graphlab.hpp>
#include <graphlab/macros_def.hpp>


/**

 Vector utilizies written by Danny Bickson, CMU


*/


#ifndef _min
#define _min(a,b) (a>b)?b:a
#endif


using namespace graphlab;

struct edata{
  int from;
  int to;
  double weight;
};
struct edata2{
  int from;
  int to;
  double time;
  double weight;
};
struct edata3{
  float from;
  float to;
  float time;
  float weight;
};
  
sdouble * clean_vec(sdouble * vec){
  delete [] vec;
  return NULL;
}


template<typename graph>
sdouble * collect_vec(int start_pos, int end_pos, int offset,
                      graph * g, sdouble *& ret);



template<typename graph>
sdouble * collect_vec(int start_pos, int end_pos, int offset, graph * g){
  sdouble * ret = new sdouble[end_pos - start_pos];
  return collect_vec(start_pos, end_pos, offset, g, ret);

}
  
template<typename graph>
sdouble * collect_vec(int start_pos, int end_pos, int offset, graph* g,
                      sdouble *& ret){
  typedef typename graph::vertex_data_type vdata_type;
  typedef typename graph::edge_data_type edata_type;

  assert(end_pos - start_pos > 1);
  assert(offset>=0);
  assert(ret != NULL);
  for (int i=start_pos; i < end_pos; i++){
    double* data = (double*)&(g->vertex_data(i));
    ret[i-start_pos] = data[offset]; 
  }
  return ret;

}


     
template<typename graph>
double  collect_diff(int start_pos, int end_pos, graph * g,
                     int first, int second){

  assert(end_pos - start_pos > 1);
  assert(first>=0 && second >=0 && first != second);
  sdouble diff = 0;
  for (int i=start_pos; i < end_pos; i++){
    sdouble * data = (sdouble*)&g->vertex_data(i);
    diff += powf(data[first] - data[second],2); 
  }
  return sqrt(diff);
}



template<typename graph>
void dispatch_vec(int start_pos, int end_pos, int offset, graph * g,
                  float * vec, int len, bool free_vec = false, int d = 1){
  assert(end_pos - start_pos > 1);
  assert(len == end_pos - start_pos);
  assert(vec != NULL);
  assert(offset>=0 && offset < len);
  for (int i=start_pos; i < end_pos; i++){
    sdouble * data = (sdouble*)&g->vertex_data(i);
    for (int j=0; j<d; j++){
      data[offset+d*j] = vec[d*(i-start_pos)+j];
    }
  }
  if (free_vec)
    delete [] vec;
}
  
template<typename graph>
void dispatch_vec(int start_pos, int end_pos, int offset, graph * g,
                  double * vec, int len, bool free_vec = false, int d = 1){
  assert(end_pos - start_pos > 1);
  assert(len == end_pos - start_pos);
  assert(vec != NULL);
  assert(offset>=0 && offset < len);
  for (int i=start_pos; i < end_pos; i++){
    sdouble * data = (sdouble*)&g->vertex_data(i);
    for (int j=0; j<d; j++){
      data[offset+d*j] = vec[d*(i-start_pos)+j];
    }
  }
  if (free_vec)
    delete [] vec;
}

template<typename graph>
void dispatch_vec(int start_pos, int end_pos, int offset, graph * g,
                  sdouble val,int d = 1){
  assert(end_pos - start_pos > 1);
  for (int i=start_pos; i < end_pos; i++){
    sdouble * data = (sdouble*)&g->vertex_data(i);
    for (int j=0; j<d; j++)
      data[offset+d*j] = val;
  }
}


template<typename graph>
void dispatch_increase(int start_pos, int end_pos, int offset, graph * g,
                       sdouble val){
  assert(end_pos - start_pos > 1);
  for (int i=start_pos; i < end_pos; i++){
    sdouble * data = (sdouble*)&g->vertex_data(i);
    data[offset] += val;
  }
}



template<typename graph>
void update_all_edges(int start_pos, int end_pos, int offset, graph * g,
                      sdouble * _vec, int len, bool free_vec = false){
  assert(end_pos - start_pos > 1);
  assert(len == end_pos - start_pos);
  assert(_vec != NULL);
  assert(offset>=0 && offset < len);
  for (int i=start_pos; i < end_pos; i++){
    foreach(edge_id_t eid, g->out_edge_ids(i)){          
      sdouble * data= (sdouble*)&g->edge_data(eid);
      memcpy(&data[offset], &_vec[i-start_pos], sizeof(sdouble));
      assert(_vec[i-start_pos] != 0);
    }
    foreach(edge_id_t eid, g->in_edge_ids(i)){          
      sdouble * data= (sdouble*)&g->edge_data(eid);
      memcpy(&data[offset], &_vec[i-start_pos], sizeof(sdouble));
      assert(_vec[i-start_pos] != 0);
    }

  }  
  if (free_vec)
    _vec = clean_vec(_vec);
}


template<typename graph>
sdouble * multiply_by_out_edge(int start_pos, int end_pos, int offset,
                               graph * g,sdouble * vec, int len){
  // assert(len == end_pos - start_pos); 
  sdouble * ret = new sdouble[len];
  memset(ret, 0, sizeof(sdouble)*len);
  assert(ret != NULL);
  for (int i=start_pos; i<end_pos; i++){
    foreach(edge_id_t eid, g->out_edge_ids(i)){          
      sdouble * data= (sdouble*)&g->edge_data(eid);
      sdouble weight;
      memcpy(&weight, &data[offset], sizeof(sdouble));
      assert(weight != 0);
      int to = g->target(eid);
      assert(to != i);
      assert(to>=0 && to < (int)g->num_vertices());
      ret[i-start_pos] += (vec[to-end_pos] * weight);
    }  
  }
  return ret; 
}

template<typename graph>
sdouble * multiply_by_in_edge(int start_pos, int end_pos, int offset,
                              graph * g,sdouble * vec, int len, int outlen){
  assert(len == end_pos - start_pos); 
  sdouble * ret = new sdouble[outlen];
  memset(ret, 0, sizeof(sdouble)*outlen);
  assert(ret != NULL);
  for (int i=start_pos; i<end_pos; i++){
    foreach(edge_id_t eid, g->in_edge_ids(i)){          
      sdouble * data= (sdouble*)&g->edge_data(eid);
      sdouble weight;
      memcpy(&weight, &data[offset], sizeof(sdouble));
      assert(weight != 0);
      int to = g->target(eid);
      int from = g->source(eid);
      assert(to == i);
      assert(from != i);
      ret[from-end_pos] += (vec[to] * weight);
    }  
  }
  return ret; 
}
template<typename graph>
void normalize(int start_pos, int end_pos, graph * g, int offset){
  assert(end_pos > start_pos && start_pos >= 0);
  assert(offset >= 0 && offset < 10);
  int wrong = 0;
  for (int i=start_pos; i<end_pos; i++){
    sdouble sum =0;
    foreach(edge_id_t eid, g->out_edge_ids(i)){          
      sdouble * data= (sdouble*)&g->edge_data(eid);
      sum += data[offset] * data[offset];
    }
    if (sum == 0){
      printf("Warning: node %d has no edges!!\n", i);
      wrong++;
    }
    //assert(sum > 0);
    sum = sqrt(sum);
    if (i % 2 == 0)
      sum = -sum;

    foreach(edge_id_t eid, g->out_edge_ids(i)){          
      sdouble * data= (sdouble*)&g->edge_data(eid);
      data[offset] /= sum;
    }
    foreach(edge_id_t eid, g->in_edge_ids(i)){          
      sdouble * data= (sdouble*)&g->edge_data(eid);
      data[offset] /= sum;
    }

    if (wrong >0)
      printf("total wrong nodes are %d\n", wrong); 
	 
  }
}


template<typename graph>
int read_edges(FILE * f, int len, int offset, int nodes, graph * g,
               bool symmetry = false){
  assert(offset>=0 && offset < len/(int)sizeof(sdouble));

  typedef typename graph::edge_data_type etype;

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
  while(true){
    memset(ed, 0, 200000*sizeof(edata));
    rc = (int)fread(ed, sizeof(edata),
                    _min(200000, edgecount_in_file - total), f);
    total += rc;

    sdouble tmp[len/sizeof(sdouble)];
    for (int i=0; i<rc; i++){
      memset(tmp, 0, len/sizeof(sdouble));
      tmp[offset] =  ed[i].weight;
      assert(ed[i].weight != 0);
      assert(ed[i].from >= 1 && ed[i].from <= nodes);
      assert(ed[i].to >= 1 && ed[i].to <= nodes);
      assert(ed[i].to != ed[i].from);
      // Matlab export has ids starting from 1, ours start from 0
      g->add_edge(ed[i].from-1, ed[i].to-1, *(etype*)&tmp);
      if (symmetry) { //add the reverse edge as well
        // Matlab export has ids starting from 1, ours start from 0
        g->add_edge(ed[i].to-1, ed[i].from-1, *(etype*)&tmp); 
        // printf("adding an edge %d -> %d (%e) \n", ed[i].from, ed[i].to, ed[i].weight);
      }
    }
    printf(".");
    fflush(0);
    if (rc == 0 || total >= edgecount_in_file)
      break;
  }
  delete [] ed; ed = NULL;
  return e;
}

template<typename graph>
int read_edges2(FILE * f, int len, int offset, int nodes, graph * g,
                bool symmetry = false){
  assert(offset>=0 && offset < len);

  typedef typename graph::edge_data_type edge_data;
     
  unsigned int e;
  int rc = fread(&e,1,4,f);
  assert(rc == 4);
  printf("Creating %d edges...\n", e);
  assert(e>0);
  int total = 0;
  edata2* ed = new edata2[200000];
  int edgecount_in_file = e;
  while(true){
    memset(ed, 0, 200000*sizeof(edata2));
    rc = (int)fread(ed, sizeof(edata2),
                    _min(200000, edgecount_in_file - total), f);
    total += rc;

    edge_data edge;
    for (int i=0; i<rc; i++){
      assert(ed[i].weight != 0); // && ed[i].weight <= 5);
      assert(ed[i].from >= 1 && ed[i].from <= nodes);
      assert(ed[i].to >= 1 && ed[i].to <= nodes);
      assert(ed[i].to != ed[i].from);
      edge.weight = ed[i].weight;
      edge.time = ed[i].time - 1;
      // Matlab export has ids starting from 1, ours start from 0
      g->add_edge(ed[i].from-1, ed[i].to-1, edge); 
    }
    printf(".");
    fflush(0);
    if (rc == 0 || total >= edgecount_in_file)
      break;
  }
  assert(total == (int)e);
  delete [] ed; ed = NULL;
  return e;
}


#define BUFSIZE 500000
template<typename graph>
void read_nodes(FILE * f, int len, int offset, int nodes, graph * g){
  typedef typename graph::vertex_data_type vtype;

  assert(offset>=0 && offset < len);
  assert(nodes>0);

  int toread = nodes;
  if (nodes > BUFSIZE)
    toread = BUFSIZE;
  int remain = nodes;

  while(remain > 0){
    double * temp = new double[remain];
    toread = (remain < BUFSIZE)?remain:toread;
    int rc = (int)fread(temp, sizeof(double), toread, f);
    //assert(rc == toread);
    remain -= rc;
    sdouble ndata[len];
    for (int i=0; i< rc; i++){
      memset(ndata,0,len*sizeof(sdouble));
      ndata[offset] = temp[i];
      g->add_vertex(*(vtype*)&ndata);
    }
    delete [] temp;
  }

   
   
}
  
  

double * read_vec(FILE * f, int len){
  double * vec = new double[len];
  assert(vec != NULL);
  fread(vec, len, sizeof(sdouble), f);
  return vec;
}


void debug_print_vec(const char * name,const sdouble * vec, int len){
  printf("%s ) ", name);
  for (int i=0; i< len; i++)
    if (vec[i] == 0)
      printf("      0    ");
    else printf("%12.4g    ", vec[i]);
  printf("\n");
}

void debug_print_matrix(const char * name,sdouble ** mat, int x, int y){
  for (int i=0; i< x; i++)
    debug_print_vec(name, mat[i], y);
} 


double current_time()  {
  timeval current_time;
  gettimeofday(&current_time, NULL);
  double answer =
    (current_time.tv_sec + ((double)current_time.tv_usec)/1.0E6);
  return answer;
}

template<typename graph>
void debug_print_graph(const char * name, graph * g, int start_pos,
                       int end_pos, int offset, int m, int n){
  assert(start_pos >= 0 && start_pos < end_pos && offset>=0);
  assert(m>0 && n>0);
  sdouble ** mat = new sdouble*[m];
  for (int i=0; i < m; i++){
    mat[i] = new sdouble[n];
    memset(mat[i], 0, n*sizeof(sdouble));
    foreach(edge_id_t eid, g->out_edge_ids(i)){          
      sdouble * data= (sdouble*)&g->edge_data(eid);
      int to = g->target(eid);
      int col = to-m;
      assert(i>=0 && col>=0 && col<n);
      mat[i][to-m] = data[offset];
    }
  }  
  debug_print_matrix(name, mat,m,n);
  for (int i=0; i<m; i++)
    delete [] mat[i];
  delete [] mat;
} 

#include <graphlab/macros_undef.hpp>
#endif


