#include <stdint.h>
#include <string>
#include <cxxtest/TestSuite.h>
#include <tasks/function_types.hpp>
#include <tasks/update_task.hpp>
#include <graph/graph.hpp>
#include <graph/scope/unsync_scope_factory.hpp>
#include <graph/scope/locked_scope_factory.hpp>
#include <schedulers/set_scheduler/set_scheduler.hpp>
#include <schedulers/set_scheduler/vertex_set.hpp>
#include <schedulers/set_scheduler/restricted_vertex_set.hpp>
#include <schedulers/set_scheduler/multinomial_vertex_set.hpp>
#include <engine/single_thread_engine.hpp>
#include <engine/multi_thread_engine.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <graphlab/logger/logger.hpp>

#include <macros_def.hpp>

using namespace graphlab;

struct MyVertexData {
  double vdata;
};


void MultiplyFunction(iscope &scope,
                      ischeduler_callback &scheduler) {
  MyVertexData* d = (MyVertexData*)scope.vertex_blob().mutable_data();
  logger(LOG_INFO, "Run scope vertex %d with vertex data: %f", scope.vertex(),
                                                               d->vdata);
//  TS_ASSERT_LESS_THAN(size_t(0),d->vdata);
  d->vdata = rand() % 2;
}

/****************  Test 1: Selector with rep execution **************************/
bool nonzero(vertex_id_t v, const blob &b) {
  MyVertexData* d = (MyVertexData*)b.data;
  logger(LOG_INFO, "Run nonzero on vertex %d, %f ", v, d->vdata);
  return d->vdata > 0;
}


void test_schedule(set_scheduler &sched) {
  // All sets must be created before scheduling calls
  ivertex_set &vs = sched.attach(restricted_vertex_set(nonzero), sched.root_set());
  
  sched.init();
//  while(vs.size() > 0) {
  sched.execute_rep(vs, MultiplyFunction);
//    sched.wait();
//  }
}



/**********  Test 2: Scoped Selector with phased execution **********************/

void MultiplyFunctionWithAssert(iscope &scope,
                      ischeduler_callback &scheduler,
                      const shared_data_manager &data_manager) {
  MyVertexData* d = (MyVertexData*)scope.vertex_blob().mutable_data();
  logger(LOG_INFO, "Run scope vertex %d with vertex data: %f", scope.vertex(),
                                                               d->vdata);
  TS_ASSERT_LESS_THAN(size_t(0),d->vdata);
  d->vdata = rand() % 2;
}


bool nonzero_scope(vertex_id_t v, iscope &b, bool &repropagatechanges) {
  MyVertexData* d = (MyVertexData*)b.vertex_blob().data();
  logger(LOG_INFO, "Run nonzero on vertex %d, %f ", v, d->vdata);
  return d->vdata > 0;
}



void test_schedule_scope(set_scheduler &sched) {
  // All sets must be created before scheduling calls
  ivertex_set &vs = sched.attach(restricted_vertex_set(nonzero_scope), sched.root_set());
  sched.init();
  while(vs.size() > 0) {
    sched.execute(vs, MultiplyFunctionWithAssert);
    sched.wait();
  }
}



/****************  Test 3: Priority Set with phased execution***********************/


double priority(vertex_id_t v, const blob &b) {
  MyVertexData* d = (MyVertexData*)(b.data);
  logger(LOG_INFO, "Run get_priority on vertex %d, %f ", v, d->vdata);
  return d->vdata;
}

void test_schedule_priority(set_scheduler &sched) {
  // All sets must be created before scheduling calls
  multinomial_vertex_set &vs = sched.attach(multinomial_vertex_set(priority), sched.root_set());
  sched.init();
  while(vs.size() > 0) {
    sched.execute(vs, MultiplyFunctionWithAssert);
    sched.wait();
  }
}


/****************  Test 4: Selected Priority Set with execution********************/


bool less_than_selector(vertex_id_t v, const blob &b) {
  MyVertexData* d = (MyVertexData*)b.data;
  logger(LOG_INFO, "Run nonzero on vertex %d, %f ", v, d->vdata);
  return d->vdata > 1.0;
}

void test_selected_priority(set_scheduler &sched) {
  // All sets must be created before scheduling calls
  ivertex_set &vs1 = sched.attach(restricted_vertex_set(less_than_selector),sched.root_set());
  multinomial_vertex_set &mvs = sched.attach(multinomial_vertex_set(priority), vs1);
  sched.init();
  sched.execute_rep(mvs, MultiplyFunction);
}
/********************************************************************/

class SetScheduleTestSuite: public CxxTest::TestSuite {
public:

  void testBasicSetSchedule() {
    global_logger().set_log_to_console(true);
    global_logger().set_log_level(LOG_INFO);
    graph g;
    for (size_t i = 0; i < 100; ++i) {
      MyVertexData vd;
      vd.vdata = (i % 2 == 0);
      g.add_vertex(&vd, sizeof(MyVertexData));
      if (i > 0) g.add_edge(i-1,i, &vd, sizeof(MyVertexData));
    }
    
    multi_thread_engine<set_scheduler, unsync_scope_factory> engine(g, 2);
    ((set_scheduler&)engine.get_scheduler()).begin_schedule(test_schedule);
    engine.start();
    for (size_t i = 0; i < g.num_vertices(); ++i) {
        TS_ASSERT_EQUALS(((MyVertexData*)(g.vertex_blob(i).data))->vdata, size_t(0));
    }
  }

  void testBasicSetSchedule2() {
    global_logger().set_log_to_console(true);
    global_logger().set_log_level(LOG_INFO);
    graph g;
    for (size_t i = 0; i < 100; ++i) {
      MyVertexData vd;
      vd.vdata = (i % 2 == 0);
      g.add_vertex(&vd, sizeof(MyVertexData));
      if (i > 0) g.add_edge(i-1,i, &vd, sizeof(MyVertexData));
    }

    multi_thread_engine<set_scheduler, locked_scope_factory> engine(g, 2);
    ((set_scheduler&)engine.get_scheduler()).begin_schedule(test_schedule_scope);
    engine.start();
    for (size_t i = 0; i < g.num_vertices(); ++i) {
        TS_ASSERT_EQUALS(((MyVertexData*)(g.vertex_blob(i).data))->vdata, size_t(0));
    }
  }


  void testPrioritySet() {
    global_logger().set_log_to_console(true);
    global_logger().set_log_level(LOG_INFO);
    graph g;
    for (size_t i = 0; i < 100; ++i) {
      MyVertexData vd;
      vd.vdata = (i % 2 == 0);
      g.add_vertex(&vd, sizeof(MyVertexData));
      if (i > 0) g.add_edge(i-1,i, &vd, sizeof(MyVertexData));
    }

    multi_thread_engine<set_scheduler, locked_scope_factory> engine(g, 2);
    ((set_scheduler&)engine.get_scheduler()).begin_schedule(test_schedule_priority);
    engine.start();
    for (size_t i = 0; i < g.num_vertices(); ++i) {
        TS_ASSERT_EQUALS(((MyVertexData*)(g.vertex_blob(i).data))->vdata, size_t(0));
    }
  }


  void testSelectedPrioritySet() {
    global_logger().set_log_to_console(true);
    global_logger().set_log_level(LOG_INFO);
    graph g;
    for (size_t i = 0; i < 100; ++i) {
      MyVertexData vd;
      vd.vdata = double(i)/10;
      g.add_vertex(&vd, sizeof(MyVertexData));
      if (i > 0) g.add_edge(i-1,i, &vd, sizeof(MyVertexData));
    }

    multi_thread_engine<set_scheduler, locked_scope_factory> engine(g, 2);
    ((set_scheduler&)engine.get_scheduler()).begin_schedule(test_selected_priority);
    engine.start();
    for (size_t i = 0; i < g.num_vertices(); ++i) {
        TS_ASSERT_LESS_THAN_EQUALS(((MyVertexData*)(g.vertex_blob(i).data))->vdata, 1.0);
    }
  }


};





