// standard C++ headers
#include <iostream>

// includes the entire graphlab framework
#include <graphlab.hpp>
#include <graphlab/macros_def.hpp>

#define NUM_VERTICES 10000

struct vertex_data{
  int val;
  int ucount;
};


typedef char edge_data;
typedef graphlab::graph<vertex_data, edge_data> graph_type;
typedef graphlab::types<graph_type> gl;


void update_function(gl::iscope& scope,
                     gl::icallback& scheduler) {
  vertex_data& curvdata = scope.vertex_data();
  curvdata.ucount += 1;
  curvdata.val += 1;
  foreach(graphlab::edge_id_t eid, scope.in_edge_ids()) {
    graphlab::vertex_id_t sourcev = scope.source(eid);
    const vertex_data& nbrvertex =
      scope.const_neighbor_vertex_data(sourcev);
    curvdata.val += nbrvertex.val;
  }
}



void init_graph(graph_type& g,
                  size_t length) {
  vertex_data vd;
  vd.val = 0;
  vd.ucount = 0 ;
  for (graphlab::vertex_id_t i = 0;i < (graphlab::vertex_id_t)length; ++i) {
    g.add_vertex(vd);
    if (i > 0) {
      g.add_edge(i-1, i, char(0));
      g.add_edge(i, i-1, char(0));
    }
  }
  g.finalize();
}


bool test_graphlab_static(gl::core &glcore,
                          bool sequentializationcheck){
  init_graph(glcore.graph(),NUM_VERTICES);
  glcore.add_task_to_all(update_function, 1.0);

  glcore.start();

  TS_ASSERT_EQUALS(glcore.engine().last_update_count(),
                   size_t(NUM_VERTICES));
  
  for (graphlab::vertex_id_t i = 0;i < NUM_VERTICES; ++i) {
    if (glcore.graph().vertex_data(i).ucount < 1 ||
        glcore.graph().vertex_data(i).val < 1) {
      return false;
    }
  }

  if (sequentializationcheck) {
    for (graphlab::vertex_id_t i = 1;i < NUM_VERTICES; ++i) {
      if (glcore.graph().vertex_data(i).val ==
            glcore.graph().vertex_data(i-1).val) {
        return false;
      }
    }
  }
  return true;
}

bool test_graphlab_round_robin(gl::core &glcore,
                            bool sequentializationcheck,
                            int iterations,
                            int skiptov){

  init_graph(glcore.graph(),NUM_VERTICES);

  glcore.sched_options().add_option("max_iterations",iterations);
  size_t startv = rand() % 10000;
  glcore.sched_options().add_option("start_vertex",startv);

  for (graphlab::vertex_id_t v = skiptov;v < NUM_VERTICES; ++v) {
    glcore.add_task(gl::update_task(v, update_function), 1.0);
  }
  
  glcore.start();

  size_t numupdates = iterations *  (NUM_VERTICES - skiptov);
  TS_ASSERT_EQUALS(glcore.engine().last_update_count(), numupdates);
  
  for (graphlab::vertex_id_t i = 0;i < NUM_VERTICES; ++i) {
    if (i < size_t(skiptov)) {
      TS_ASSERT_EQUALS(glcore.graph().vertex_data(i).ucount, 0);
    }
    else {
      TS_ASSERT_EQUALS(glcore.graph().vertex_data(i).ucount, iterations);
    }
  }
  
  if (sequentializationcheck) {
    for (graphlab::vertex_id_t i = skiptov;i < NUM_VERTICES; ++i) {
      if (glcore.graph().vertex_data(i).val ==
                glcore.graph().vertex_data(i-1).val) {
        return false;
      }
    }
  }
  return true;
}



bool test_graphlab_colored(gl::core &glcore,
                          bool sequentializationcheck,
                          size_t numiterations){
  init_graph(glcore.graph(),NUM_VERTICES);
  glcore.graph().compute_coloring();
  glcore.add_task_to_all(update_function, 1.0);
  glcore.sched_options().add_option("max_iterations",numiterations);
  glcore.start();
  
  TS_ASSERT_EQUALS(glcore.engine().last_update_count(),
                   size_t(NUM_VERTICES * numiterations));
  
  for (graphlab::vertex_id_t i = 0;i < NUM_VERTICES; ++i) {
    if (glcore.graph().vertex_data(i).ucount < 1 ||
      glcore.graph().vertex_data(i).val < 1) {
      return false;
      }
  }
  
  if (sequentializationcheck) {
    for (graphlab::vertex_id_t i = 1;i < NUM_VERTICES; ++i) {
      if (glcore.graph().vertex_data(i).val ==
        glcore.graph().vertex_data(i-1).val) {
        return false;
        }
    }
  }
  return true;
}


class GraphlabTestSuite: public CxxTest::TestSuite {
public:


  void test_static(void) {
    global_logger().set_log_level(LOG_WARNING);
    global_logger().set_log_to_console(true);
    
    const char* engine_types[] = {"async", "async_sim"};
    const char* scope_types[] = {"vertex", "edge", "full"};
    const char* schedulers[]  = {"fifo", "multiqueue_fifo", "priority", "multiqueue_priority", "sweep", "clustered_priority"};
    std::cout << "\n\n\n";
    std::cout << "engine\tscheduler\tscope\tncpus" << std::endl;
    for (size_t e = 0;e < 2; ++e) {
      for (size_t c = 0; c < 3; ++c) {
        for (size_t s = 0;s < 6; ++s) {
          // skip async_sim + sweep/clustered. doesn't work
          if (e == 1 && s >= 4) continue; 
          for (size_t n =1; n <= 4; ++n) {
            gl::core glcore;
            glcore.set_engine_type(engine_types[e]);
            glcore.set_scheduler_type(schedulers[s]);
            glcore.set_scope_type(scope_types[c]);
            glcore.set_ncpus(n);
            std::cout << engine_types[e] << "\t" << schedulers[s] << "\t" << scope_types[c] << "\t" << n << std::endl;
            
            bool ret = test_graphlab_static(glcore,
                                            c != 0);
            // sequentialization check only if not vertex scope
            TS_ASSERT_EQUALS(ret, true);
          }
        }
      }
    }

  }


  void test_round_robin(void) {
    global_logger().set_log_level(LOG_WARNING);
    global_logger().set_log_to_console(true);
    std::cout << "\n\n\n";
    std::cout << "engine\tscheduler\tscope\tncpus\titerations" << std::endl;
    const char* engine_types[] = {"async", "async_sim"};
    const char* scope_types[] = {"vertex", "edge", "full"};
    for (size_t e = 0;e < 2; ++e) {
      for (size_t c = 0; c < 3; ++c) {
        for (size_t n =1; n <= 4; ++n) {
          for (size_t iter = 1;iter < 4; ++iter) {
            gl::core glcore;
            glcore.set_engine_type(engine_types[e]);
            glcore.set_scheduler_type("round_robin");
            glcore.set_scope_type(scope_types[c]);
            glcore.set_ncpus(n);
            std::cout << engine_types[e] << "\t" << "round_robin" << "\t" << scope_types[c] << "\t" << n << "\t" << iter << std::endl;
            bool ret = test_graphlab_round_robin(glcore,
                                                c != 0,
                                                (int)iter,
                                                (graphlab::vertex_id_t)(rand() % 1000));
                        // sequentialization check only if not vertex scope
            TS_ASSERT_EQUALS(ret, true);
          }
        }
      }
    }

  }



  void test_colored(void) {
    global_logger().set_log_level(LOG_WARNING);
    global_logger().set_log_to_console(true);
    std::cout << "\n\n\n";
    std::cout << "engine\tscheduler\tscope\tncpus\titerations" << std::endl;
    const char* engine_types[] = {"async", "async_sim"};
    const char* scope_types[] = {"vertex", "edge", "full"};
    for (size_t e = 0;e < 2; ++e) {
      for (size_t c = 0; c < 3; ++c) {
        for (size_t n =1; n <= 4; ++n) {
          for (size_t iter = 1;iter < 4; ++iter) {
            gl::core glcore;
            glcore.set_engine_type(engine_types[e]);
            glcore.set_scheduler_type("chromatic");
            glcore.set_scope_type(scope_types[c]);
            glcore.set_ncpus(n);
            
            
            std::cout << engine_types[e] << "\t" << "colored" << "\t" << scope_types[c] << "\t" << n << "\t" << iter << std::endl;
            bool ret = test_graphlab_colored(glcore,
                                            c != 0,
                                            iter);
            // sequentialization check only if not vertex scope
            TS_ASSERT_EQUALS(ret, true);
          }
        }
      }
    }

  }
};
