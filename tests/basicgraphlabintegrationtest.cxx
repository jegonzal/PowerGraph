#include <stdint.h>
#include <string>
#include <cxxtest/TestSuite.h>


#include <graphlab.hpp>

#include <graphlab/macros_def.hpp>


struct MyVertexData {
  double DBLVAL;
};

struct MyEdgeData {
  int64_t INT64VAL;
};


typedef graphlab::types<MyVertexData, MyEdgeData> gl_types;

void MultiplyFunction(gl_types::iscope& scope,
                      gl_types::icallback& scheduler) {
  MyVertexData& vdata = scope.vertex_data();
  double curval = vdata.DBLVAL;
  logger(LOG_INFO, 
         "Run scope vertex %d with vertex data: %lf", scope.vertex(), curval);

  vdata.DBLVAL = curval * 100.0;

  curval = vdata.DBLVAL;
  logger(LOG_INFO, "new vertex data: %lf", curval);

  foreach(edge_id_t edge, scope.out_edge_ids()) {
    if (scope.target(edge) > scope.vertex()) {
      scheduler.add_task(update_task<blob_graph>(scope.target(edge), 
                                                 MultiplyFunction), 0.0);
    }
  }
}


class BasicGraphlabIntegrationTestSuite: public CxxTest::TestSuite {
public:
  
  void CreateGraph(blob_graph &g, size_t &vid1, size_t &vid2, size_t &vid3) {
    
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
  }
  
  
  void TestUnsyncScope(void) {

    blob_graph g;
    size_t vid1, vid2, vid3;
    CreateGraph(g, vid1, vid2, vid3);
    global_logger().set_log_to_console(true);
    global_logger().set_log_file("test.logger");
    global_logger().set_log_level(LOG_INFO);
    
    logger(LOG_INFO, "Test logger");
    TS_TRACE("Set graphlab args");
    single_thread_engine<blob_graph, 
      fifo_scheduler<blob_graph>, 
      general_scope_factory<blob_graph> > graphlab(g);

    graphlab.set_default_scope(scope_range::VERTEX_CONSISTENCY);
    
    TS_TRACE("Create func");
    
    TS_TRACE("Add tasks");
    // Add only vertex 1, it will spawn vertex 2 which will spawn vertex 3
    graphlab.get_scheduler().
      add_task(update_task<blob_graph>(vid1, MultiplyFunction), 1.0);
    
    TS_TRACE("Going to start");
    /* Start */
    graphlab.start();
    
    TS_TRACE("Finished");
    
    /* Check values changed as expected */
    MyVertexData * vd1 = g.vertex_data(vid1).as_ptr<MyVertexData>();
    MyVertexData * vd2 = g.vertex_data(vid2).as_ptr<MyVertexData>();
    MyVertexData * vd3 = g.vertex_data(vid3).as_ptr<MyVertexData>();
    TS_ASSERT_EQUALS(0.0,  vd1->DBLVAL);
    TS_ASSERT_EQUALS(100.0, vd2->DBLVAL);
    TS_ASSERT_EQUALS(200.0, vd3->DBLVAL);
  }

  void TestLockedScope(void) {

    blob_graph g;
    size_t vid1, vid2, vid3;
    CreateGraph(g, vid1, vid2, vid3);
    global_logger().set_log_to_console(true);
    global_logger().set_log_file("test.logger");
    global_logger().set_log_level(LOG_INFO);

    logger(LOG_INFO, "Test logger");
    TS_TRACE("Set graphlab args");
    single_thread_engine<blob_graph,
      fifo_scheduler<blob_graph>, 
      general_scope_factory<blob_graph> > graphlab(g);

    graphlab.set_default_scope(scope_range::EDGE_CONSISTENCY);
    
    TS_TRACE("Create func");

    TS_TRACE("Add tasks");
    // Add only vertex 1, it will spawn vertex 2 which will spawn vertex 3
    graphlab.get_scheduler().
      add_task(update_task<blob_graph>(vid1, MultiplyFunction), 1.0);

    TS_TRACE("Going to start");
    /* Start */
    graphlab.start();

    TS_TRACE("Finished");

    /* Check values changed as expected */
    MyVertexData * vd1 = (MyVertexData *) g.vertex_data(vid1).data();
    MyVertexData * vd2 = (MyVertexData *) g.vertex_data(vid2).data();
    MyVertexData * vd3 = (MyVertexData *) g.vertex_data(vid3).data();
    TS_ASSERT_EQUALS(0.0,  vd1->DBLVAL);
    TS_ASSERT_EQUALS(100.0, vd2->DBLVAL);
    TS_ASSERT_EQUALS(200.0, vd3->DBLVAL);
  }
};



