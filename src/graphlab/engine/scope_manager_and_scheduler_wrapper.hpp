#ifndef GRAPHLAB_SCOPE_MANAGER_AND_SCHEDULER_WRAPPER_HPP
#define GRAPHLAB_SCOPE_MANAGER_AND_SCHEDULER_WRAPPER_HPP

#include <graphlab/graph/graph.hpp>
#include <graphlab/engine/iengine.hpp>
#include <graphlab/schedulers/scheduler_options.hpp>

#include <graphlab/scope/iscope.hpp>
#include <graphlab/tasks/update_task.hpp>

namespace graphlab {

  /** \internal
   * This is the implementation for the scope manager and scheduler
   * "dynamic reconstruction" capabilities of the engine.
   */
  template <typename Graph, typename Scheduler, typename ScopeFactory>
  class scope_manager_and_scheduler_wrapper : public iengine<Graph>{
  public:
    typedef typename Scheduler::update_task_type update_task_type;
    typedef typename Scheduler::update_function_type update_function_type;

  private:

    Graph& graph;

    size_t ncpus;

    /** Responsible for managing the update of scopes */
    ScopeFactory *scope_manager;

    /** Responsible for maintaining the schedule over tasks
     * The scheduler is instantiated as late as possible.
     * Its first instantiation only happens when add task is called
     */
    Scheduler *scheduler;

    /// lazy deletion.
    bool deletionmark; // = has_engine_run
  
    scheduler_options schedopts;

    struct {
      size_t nvertices;
      size_t nedges;
      size_t changeid;
    }  graphtracker ;
  
    void construct_members() {
      // if deletion mark is set
      // We only recreate the members if the graph changed
      // if the members did not change, clear the deletionmark
      // (members are now fine and useable!) and return;
      if (deletionmark) {
        if (graph_changed()) {
          // delete the current members
          delete scheduler;
          delete scope_manager;
          scheduler = NULL;
          scope_manager = NULL;
          deletionmark = false;
          // re-query the graph
        } else {          
          assert(scheduler != NULL);
          assert(scope_manager != NULL);
          deletionmark = false;
          return;
        }
      }

      // really, both should be NULL
      if (scheduler == NULL && scope_manager == NULL) {
        scheduler = new Scheduler(this, graph, std::max(ncpus, size_t(1)));
        scheduler->set_options(schedopts);
        scope_manager = new ScopeFactory(graph, std::max(ncpus, size_t(1)));
        update_graph_tracker();
      }
      assert(scheduler != NULL);
      assert(scope_manager != NULL);
    }


    void update_graph_tracker() {
      graphtracker.nvertices = graph.num_vertices();
      graphtracker.nedges = graph.num_edges();
      graphtracker.changeid = graph.get_changeid();
    }

    bool graph_changed() {
      return (graph.num_vertices() != graphtracker.nvertices ||
              graph.num_edges() != graphtracker.nedges ||
              graph.get_changeid() != graphtracker.changeid);
    }

    // The graph should not be modified once
    // tasks are being inserted since there is no gaurantee
    // that the tasks will still be valid
    void check_if_add_task_is_safe() {
      // if scheduler was already insteantiated
      if (scheduler != NULL && deletionmark == false) {
        // and the graph has changed
        ASSERT_MSG(!graph_changed(),
                   "Graph modifications not allowed once tasks have been added!");
      } else {
        construct_members();
        update_graph_tracker();
      }
    }
  protected:

    Scheduler* get_scheduler() {
      construct_members();
      return scheduler;
    }

    void apply_scheduler_options() {
      construct_members();
      scheduler->set_options(schedopts);
    }

    ScopeFactory* get_scope_manager() {
      construct_members();
      return scope_manager;
    }

    void release_scheduler_and_scope_manager() {
      // lazy deletion
      deletionmark = true;
    }
  
  public:
    scope_manager_and_scheduler_wrapper(Graph &graph, size_t ncpus) :
      graph(graph),
      ncpus(ncpus),
      scope_manager(NULL),
      scheduler(NULL),
      deletionmark(false) {
      update_graph_tracker();
    }

    virtual ~scope_manager_and_scheduler_wrapper() {
      if (scheduler != NULL) {
        delete scheduler;
        scheduler = NULL;
      }
      if (scope_manager != NULL) {
        delete scope_manager;
        scope_manager = NULL;
      }
    }

    void set_scheduler_options(const scheduler_options& opts) {
      schedopts = opts;
      apply_scheduler_options();
    }
    

    /**
     * Adds an update task with a particular priority.
     * This function is forwarded to the scheduler.
     */
    void add_task(update_task_type task, double priority) {
      check_if_add_task_is_safe();
      scheduler->add_task(task, priority);
    }

    /**
     * Creates a collection of tasks on all the vertices in
     * 'vertices', and all with the same update function and priority
     * This function is forwarded to the scheduler.
     */
    void add_tasks(const std::vector<vertex_id_t>& vertices,
                   update_function_type func, double priority) {
      check_if_add_task_is_safe();
      scheduler->add_tasks(vertices, func, priority);
    }

    /**
     * Creates a collection of tasks on all the vertices in the graph,
     * with the same update function and priority
     * This function is forwarded to the scheduler.
     */
    void add_task_to_all(update_function_type func,
                         double priority) {
      check_if_add_task_is_safe();
      scheduler->add_task_to_all(func, priority);
    }
  };

}
#endif
