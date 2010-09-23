#include <stdint.h>
#include <string>
#include <cxxtest/TestSuite.h>
#include <tasks/function_types.hpp>
#include <tasks/update_task.hpp>
#include <graph/graph.hpp>
#include <graph/scope/unsync_scope_factory.hpp>
#include <graph/scope/locked_scope_factory.hpp>
#include <schedulers/set_scheduler/set_scheduler.hpp>
#include <schedulers/set_scheduler/selector.hpp>
#include <schedulers/set_scheduler/scoped_selector.hpp>
#include <schedulers/set_scheduler/multinomial_vertex_set.hpp>
#include <schedulers/set_scheduler/assign_priority.hpp>
#include <schedulers/set_scheduler/edge_change_selector.hpp>
#include <engine/single_thread_engine.hpp>
#include <engine/multi_thread_engine.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <graphlab/logger/logger.hpp>

#include <macros_def.hpp>

using namespace graphlab;

struct MyVertexData {
  double val;
};


void MultiplyFunction(iscope &scope,
                      ischeduler_callback &scheduler,
                      const shared_data_manager &data_manager) {
  MyVertexData* d = (MyVertexData*)scope.vertex_blob().mutable_data();
  logger(LOG_INFO, "Run scope vertex %d with vertex data: %f", scope.vertex(),
                                                               d->val);
//  TS_ASSERT_LESS_THAN(size_t(0),d->val);
  foreach(edge_id_t ei, scope.in_edge_ids()) {
    MyVertexData* inedata = (MyVertexData*)scope.in_edge_blob(ei).data();
    d->val = inedata->val;
  }
  
  foreach(edge_id_t ei, scope.out_edge_ids()) {
    MyVertexData* edata = (MyVertexData*)scope.out_edge_blob(ei).mutable_data();
    edata->val = d->val;
  }
}


double get_priority(vertex_id_t v, const blob &b) {
  MyVertexData* d = (MyVertexData*)b.data;
  logger(LOG_INFO, "Get Priority at %d with %f", v,d->val);
  if (d->val <= 0.01) return 0;
  else return d->val;
}

bool edge_selector(vertex_id_t cur,
                     blob_wrapper data,
                     vertex_id_t src,
                     read_only_blob_wrapper edgedata,
                     bool& repropagate) {
  MyVertexData* vd = (MyVertexData*)data.data();
  MyVertexData* ed = (MyVertexData*)edgedata.data();
  logger(LOG_INFO, "Run eselect %d with %f, %f", cur,vd->val, ed->val);
  repropagate = true;
  return ed->val != vd->val;
}

void test_schedule(set_scheduler &sched) {
  // All sets must be created before scheduling calls
  ivertex_set &vs1 = sched.attach(edge_change_selector(edge_selector),sched.root_set());
  multinomial_vertex_set &mvs = sched.attach(assign_priority_to_multinomial(get_priority), vs1);
  sched.init();
  sched.execute_rep(mvs, MultiplyFunction);
}
/********************************************************************/

class SetScheduleTestSuite: public CxxTest::TestSuite {
public:

  void testEdgeSelector() {
    global_logger().set_log_to_console(true);
    global_logger().set_log_level(LOG_INFO);
    graph g;
    // make a chain graph
    for (size_t i = 0; i < 100; ++i) {
      MyVertexData vd;
      vd.val = i;
      g.add_vertex(&vd, sizeof(MyVertexData));
      vd.val = i-1;
      if (i > 0) g.add_edge(i-1,i, &vd, sizeof(MyVertexData));
    }
    g.finalize();
    multi_thread_engine<set_scheduler, unsync_scope_factory> engine(g, 1);
    ((set_scheduler&)engine.get_scheduler()).begin_schedule(test_schedule);
    engine.start();
    for (size_t i = 0; i < g.num_vertices(); ++i) {
        TS_ASSERT_EQUALS(((MyVertexData*)(g.vertex_blob(i).data))->val, size_t(0));
    }
  }





};





