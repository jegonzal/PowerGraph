#ifndef GRAPHLAB_CORE_HPP
#define GRAPHLAB_CORE_HPP


#include <graphlab.hpp>
#include <graphlab/metrics/metrics.hpp>
#include <graphlab/metrics/reporters/basic_reporter.hpp>
#include <graphlab/metrics/reporters/file_reporter.hpp>
#include <graphlab/metrics/reporters/html_reporter.hpp>

#include <graphlab/macros_def.hpp>
namespace graphlab {
  template<typename Graph> struct types;
  
  /**
   * A GraphLab core is the base (or core) data structure in GraphLab.
   * The core contains the 

   <ol>
   <li> Data Graph: which represents the structured data dependencies. </li>
   <li> Shared Data: which represents the global constants and mutable data. </li>
   <li> Engine: The computational structure which contains the
        scheduling and execution statistics for the GraphLab program. </li>
   </ol>
   
   * 
   */
  template <typename VertexType, typename EdgeType>
  class core {
  public:
    typedef graphlab::types<graphlab::graph<VertexType, EdgeType> > types;

  public:
    //! default constructor does nothing
    core() : mengine(NULL) { }

    /**
     * Destroy the core clearing any state associated with the graph,
     * shared data or engine
     */
    ~core() { 
        destroy_engine(); 
        
        // Write options to metrics
        metrics & coremetrics = metrics::create_metrics_instance("core", true);
        coremetrics.set("ncpus", meopts.ncpus);
        coremetrics.set("engine", meopts.engine_type);
        coremetrics.set("scope", meopts.scope_type);
        coremetrics.set("scheduler", meopts.scheduler_type);
        coremetrics.set("affinities", meopts.enable_cpu_affinities ? "true" : "false");
        coremetrics.set("schedyield", meopts.enable_sched_yield ? "true" : "false");
        coremetrics.set("compile_flags", meopts.compile_flags);
        
        // Metrics dump: basic 
        if (meopts.metrics_type == "basic") { 
            basic_reporter reporter = basic_reporter();
            metrics::report_all(reporter); 
        }
        // Metrics dump: file
        if (meopts.metrics_type == "file") {
            file_reporter freporter = file_reporter("graphlab_metrics.txt");
            metrics::report_all(freporter);
        }
        if (meopts.metrics_type == "html") {
            html_reporter hreporter = html_reporter("graphlab_metrics.html");
            metrics::report_all(hreporter);
        }
     } 
       
    //! Get a reference to the graph associated with this core
    typename types::graph& graph() { return mgraph; }

    //! Get a reference to the graph associated with this core
    const typename types::graph& graph() const { return mgraph; }

    /**
     * Set the type of scheduler.  This will destroy the current
     * engine and any tasks currently associated with the scheduler.
     * See schedulers for the list of supported scheduler (strings)
     */
    void set_scheduler_type(const std::string& scheduler_type) {
      meopts.scheduler_type = scheduler_type;
      destroy_engine();
    }

    /**
     * Set the scope consistency model used in this engine.  This will
     * destroy the current engine and any tasks associated with the
     * current scheduler.  The available scopes are: 
     * <ol>
     *  <li> full : This ensures full data consistency within the scope </li>
     *  <li> edge : This ensures data consistency with just the 
                    vertex and edges </li>
     *  <li> vertex : This ensures that a vertex cannot be updated 
     *                by two processors simultaneously </li>
     *  <li> none : This eliminates all locking </li>
     * </ol>
     */
    void set_scope_type(const std::string& scope_type) {
      meopts.scope_type = scope_type;
      destroy_engine();
    }


    /**
     * Set the engine type { threaded, sequential, sim, synchronous}.
     * This will destroy the current engine and any tasks assocaited
     * with the current scheduler. 
     *
     */
    void set_engine_type(const std::string& engine_type) {
      meopts.engine_type = engine_type;
      destroy_engine();
    }

    
    /**
     * Set the number of cpus that the core will use
     */
    void set_ncpus(size_t ncpus) {
      meopts.ncpus = ncpus;
      destroy_engine();
    }


    /**
     * Get a reference to the active engine.  If no enge exists one is
     * created.
     */
    typename types::iengine& engine() { 
      bool success = auto_build_engine();
      assert(success);
      return *mengine; 
    }



    /**
     * Get a reference to the shared data associated with this core.
     */
    typename types::ishared_data_manager& shared_data() {
      return mshared_data;
    }

    
    /**
     * Get a const reference to the shared data associated with this
     * core.
     */
    const typename types::ishared_data_manager& shared_data() const {
      return mshared_data;
    }


    /**
     * Finalize the core by clearing the current engine and
     * reprocessing the engine arguments.  
     */
    bool rebuild_engine() {
      destroy_engine();
      assert(mengine == NULL);
      return auto_build_engine();
    }

    /**
     * Set the engine options by passing in an engine options object.
     */
    void set_engine_options(const engine_options& opts) {
      meopts = opts;
    }

    /**
     * View the engine options
     */
    const engine_options& get_engine_options() const { 
      return meopts;
    }

    scheduler_options& sched_options() {
      return meopts.sched_options();
    }

    const scheduler_options& sched_options() const{
      return meopts.sched_options();
    }


    /**
     * Set the engien options by simply parsing the command line
     * arguments.
     */
    bool parse_engine_options(int argc, char **argv) {
      command_line_options clopts;
      bool success = clopts.parse(argc, argv);
      assert(success);
      return set_engine_options(clopts);
    }


    /**
     * Run the engine until a termination condition is reached or
     * there are no more tasks remaining to execute.
     */
    double start() {
      bool success = auto_build_engine();
      assert(success);
      assert(mengine != NULL);
      graphlab::timer ti;
      ti.start();
      mengine->start();
      return ti.current_time();
    }
  

    /**
     * Add a single update function to a single vertex.
     */
    void add_task(vertex_id_t vertex,
                  typename types::update_function func,
                  double priority) {
      typename types::update_task task(vertex, func);
      add_task(task, priority);
    }


    /**
     * Add a single task with a fixed priority.
     */
    void add_task(typename types::update_task task, double priority) {
      engine().add_task(task, priority);
    }

    /**
     * Add the update function to all the veritces in the provided
     * vector with the given priority.
     */
    void add_tasks(const std::vector<vertex_id_t>& vertices, 
                   typename types::update_function func, double priority) {
      engine().add_tasks(vertices, func, priority);
    }


    /**
     * Add the given function to all vertices using the given priority
     */
    void add_task_to_all(typename types::update_function func, 
                         double priority) {
      engine().add_task_to_all(func, priority);
    }
    
    /**
     * Get the number of updates executed by the engine
     */
    size_t last_update_count() {
      if(mengine == NULL) return 0;
      else return mengine->last_update_count();
    }

  private:

    /**
     * Build the engine if it has not already been built.
     */
    bool auto_build_engine() {
      if(mengine == NULL) {
        // create the engine
        mengine = meopts.create_engine(mgraph);
        if(mengine == NULL) return false;
        else mengine->set_shared_data_manager(&mshared_data);
      }
      else {
        // scheduler options is one parameter that is allowed
        // to change without rebuilding the engine
        mengine->sched_options() = sched_options();
      }
      return true;
    }

    /**
     * Destroy the engine if one exists.
     */
    void destroy_engine() {
      if(mengine != NULL) {
        delete mengine;
        mengine = NULL;
      }
    }



    /** Save the core to a file */
    void save(const std::string& filename) const {
      std::ofstream fout(filename.c_str());
      assert(fout.good());
      oarchive oarc(fout);
      oarc << *this;
      fout.close();
    } // end of save
    
    /** Save the core to an archive */
    void save(oarchive& arc) const {
      arc << mgraph
          << mshared_data
          << meopts;
    } // end of save


    /** Load the core from a file. */
    void load(const std::string& filename) {
      std::ifstream fin(filename.c_str());
      assert(fin.good());
      iarchive iarc(fin);
      iarc >> *this;
      fin.close();
    } // end of load


    /** Load the core from an archive. */
    void load(iarchive& arc) {
      arc >> mgraph
          >> mshared_data
          >> meopts;
    } // end of load


    
    // graph and data objects
    typename types::graph mgraph;
    typename types::thread_shared_data mshared_data;    
    engine_options meopts;
    typename types::iengine *mengine;
      
  };

}
#include <graphlab/macros_undef.hpp>
#endif
