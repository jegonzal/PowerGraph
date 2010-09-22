#include <stdint.h>
#include <string>
#include <cxxtest/TestSuite.h>

#include <graphlab/graph/graph.hpp>
#include <graphlab/schedulers/set_scheduler/vertex_set.hpp>
#include <graphlab/schedulers/set_scheduler/restricted_vertex_set.hpp>
#include <logger/logger.hpp>

using namespace graphlab;

struct MyVertexData {
  double vdata;
};

bool is_even(vertex_id_t v, const blob &b) {
  MyVertexData* d = (MyVertexData*)b.data();
  logger(LOG_INFO, "Run is_even on vertex %d, %f ", v, d->vdata);
  bool even = (size_t(d->vdata) % 2 == 0);
  return even;
}


class VertexSetTestSuite: public CxxTest::TestSuite {
public:

  void testVertexSet() {
    global_logger().set_log_to_console(true);
    global_logger().set_log_level(LOG_INFO);
    blob_graph g;
    ss_set_type allv;
    for (size_t i = 0; i < 100; ++i) {
      MyVertexData vd;
      vd.vdata = (i % 2 == 0);
      g.add_vertex(blob( sizeof(MyVertexData), &vd) );
      ss_insert(allv, i);
    }
    // attach v2 to v1
    vertex_set<blob_graph> v1,v2;
    v1.init(&g, NULL, NULL, 1);
    v2.init(&g, &v1, NULL, 1);
    v1.resolve_event_handlers();
    // set v1 to everything
    v1.rebuild(NULL, allv);
    
    TS_ASSERT_EQUALS(v2.size(), size_t(100));

    v1.erase(NULL, 5);
    TS_ASSERT_EQUALS(v2.size(), size_t(99));

    v1.erase(NULL, 5);
    TS_ASSERT_EQUALS(v2.size(), size_t(99));

    v1.insert(NULL, 5);
    TS_ASSERT_EQUALS(v2.size(), size_t(100));
  }

  void testSelectorVertexSet() {
    global_logger().set_log_to_console(true);
    global_logger().set_log_level(LOG_INFO);
    blob_graph g;
    ss_set_type allv;
    for (size_t i = 0; i < 100; ++i) {
      MyVertexData vd;
      vd.vdata = i;
      g.add_vertex(blob(sizeof(MyVertexData), &vd) );
      ss_insert(allv, i);
    }
    // attach v2 to v1
    vertex_set<blob_graph> v1;
    restricted_vertex_set<blob_graph> v2(is_even);
    v1.init(&g, NULL, NULL, 1);
    v2.init(&g, &v1, NULL, 1);
    v1.resolve_event_handlers();
    // set v1 to everything
    v1.rebuild(NULL, allv);

    TS_ASSERT_EQUALS(v2.size(), size_t(50));

    v1.erase(NULL, 5);
    TS_ASSERT_EQUALS(v2.size(), size_t(50));

    v1.erase(NULL, 4);
    TS_ASSERT_EQUALS(v2.size(), size_t(49));

    v1.modify_vertex(NULL, 4);
    TS_ASSERT_EQUALS(v2.size(), size_t(49));

    v1.insert(NULL, 4);
    TS_ASSERT_EQUALS(v2.size(), size_t(50));
    
    
  }
};
