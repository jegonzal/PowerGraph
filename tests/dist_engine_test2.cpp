#include <stdint.h>
#include <string>

#include <graphlab/distributed/graph/cloned_graph.hpp>
#include <graphlab/distributed/distributed_shared_data.hpp>
#include <graphlab/shared_data/shared_data_ops.hpp>
#include <graphlab.hpp>


#include <graphlab/macros_def.hpp>

using namespace graphlab;

const size_t TOTALSUM = 0;
typedef types<blob_cloned_graph> gl_types;
struct MyVertexData {
  double DBLVAL;
};

struct MyEdgeData {
  int64_t INT64VAL;
};

void identity_print(size_t index,
                      const ishared_data<blob_cloned_graph>& shared_data,
                      any& current_data,
                      const any& new_data) {
  current_data.as<double>() = new_data.as<double>();
  std::cout << "We got : " << current_data.as<double>() << std::endl;
}

void merge(size_t index,
            const ishared_data<blob_cloned_graph>& shared_data,
            any& a,
            const any& b) {
  a.as<double>() += b.as<double>();
  std::cout << "Merging\n";
}

double get_number(const blob &m) {
  const MyVertexData *vdata = m.as_ptr<MyVertexData>();
  double curval = vdata->DBLVAL;
  std::cout << ".";
  std::cout.flush();
  return curval;
}

void MultiplyFunction(gl_types::iscope& scope,
                      gl_types::icallback& scheduler,
                      gl_types::ishared_data* data_manager) {
  MyVertexData * vdata = scope.vertex_data().as_ptr<MyVertexData>();
  double curval = vdata->DBLVAL;
  logger(LOG_INFO, 
         "Run scope vertex %d with vertex data: %lf", scope.vertex(), curval);

  vdata->DBLVAL = curval * 100.0;

  curval = vdata->DBLVAL;
  logger(LOG_INFO, "new vertex data: %lf", curval);

  foreach(edge_id_t edge, scope.out_edge_ids()) {
    if (scope.target(edge) > scope.vertex()) {
      scheduler.add_task(update_task<blob_cloned_graph>(scope.target(edge), 
                                                 MultiplyFunction), 0.0);
    }
  }
  timer ti;
  ti.start();
  while (ti.current_time() < 1 ) sched_yield();
}


void CreateGraph(blob_cloned_graph &g, size_t &vid1, size_t &vid2, size_t &vid3) {
  
  MyVertexData vd;
  vd.DBLVAL = 0.0;
  vid1 = g.add_vertex(blob(sizeof(MyVertexData), &vd));
  
  vd.DBLVAL = 1.0;
  vid2 = g.add_vertex(blob(sizeof(MyVertexData), &vd));
  
  vd.DBLVAL = 2.0;
  vid3 =g.add_vertex(blob(sizeof(MyVertexData), &vd));
  
  MyEdgeData e;
  e.INT64VAL = 5;
  g.add_edge(vid1, vid2, blob(sizeof(MyEdgeData), &e));
  g.add_edge(vid2, vid3, blob(sizeof(MyEdgeData), &e));
  g.add_edge(vid3, vid1, blob(sizeof(MyEdgeData), &e));
  g.finalize();
}
  
  
int main(int argc, char** argv) {
  global_logger().set_log_level(LOG_INFO);

    distributed_control dc(&argc, &argv);
    dc.init_message_processing(2);
    
    blob_cloned_graph g;
    size_t vid1, vid2, vid3;
    CreateGraph(g, vid1, vid2, vid3);
    distributed_shared_data<blob_cloned_graph> dsdm(dc);
    
    g.distributed_partition(dc, partition_method::PARTITION_BFS,1);
    g.distribute(dc);
    dc.barrier();
  
    std::cout << "Set graphlab args" << std::endl;
    distributed_engine<blob_cloned_graph,
      distributed_scheduler_wrapper<blob_cloned_graph, fifo_scheduler<blob_cloned_graph> > > graphlab(dc, g);
    dc.barrier();
    graphlab.set_default_scope(scope_range::EDGE_CONSISTENCY);
    graphlab.set_shared_data_manager(&dsdm);
    dc.barrier();
    // create the sync task
    if (dc.procid() == 0) {
      dsdm.set_tree_sync(TOTALSUM,
                    sync_ops<blob_cloned_graph>::sum<double, get_number>,
                    identity_print,
                    merge,
                    double(0),
                    100);
    }
    
    std::cout << "Add tasks" << std::endl;
    // Add only vertex 1, it will spawn vertex 2 which will spawn vertex 3
    if (dc.procid() == 0) {
      graphlab.get_scheduler().
        add_task(update_task<blob_cloned_graph>(vid1, MultiplyFunction), 1.0);
    }
    dc.barrier();
    std::cout << "Start" << std::endl;
    /* Start */
    graphlab.start();

    std::cout << "Finished" << std::endl;
    dc.barrier();
    /* Check values changed as expected */
    MyVertexData * vd1 = (MyVertexData *) g.vertex_data(vid1).data();
    MyVertexData * vd2 = (MyVertexData *) g.vertex_data(vid2).data();
    MyVertexData * vd3 = (MyVertexData *) g.vertex_data(vid3).data();
    std::cout << "v1: " << vd1->DBLVAL << std::endl;
    std::cout << "v2: " << vd2->DBLVAL << std::endl;
    std::cout << "v3: " << vd3->DBLVAL << std::endl;
  }




