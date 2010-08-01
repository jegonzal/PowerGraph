
// #include <queue>
// #include <cmath>
// #include <cassert>

// #include <graphlab/graph/graph.hpp>
// #include <graphlab/scope/iscope.hpp>
// #include <graphlab/util/synchronized_circular_queue.hpp>
// #include <graphlab/tasks/update_task.hpp>
// #include <graphlab/schedulers/ischeduler.hpp>
// #include <graphlab/schedulers/set_scheduler/ivertex_set.hpp>
// #include <graphlab/schedulers/set_scheduler/vertex_set.hpp>
// #include <graphlab/schedulers/set_scheduler/set_scheduler.hpp>
// #include <graphlab/parallel/pthread_tools.hpp>
// #include <graphlab/parallel/atomic.hpp>
// #include <graphlab/schedulers/support/unused_scheduler_callback.hpp>
// #include <graphlab/schedulers/support/vertex_task_set.hpp>
// #include <graphlab/schedulers/set_scheduler/set_generic_typedefs.hpp>
// #include <graphlab/logger/logger.hpp>


// namespace graphlab {
//   bool set_scheduler::get_task_from_active_set(size_t cpuid, update_task &ret_task) {
//     if (activeset == -1) return false;

//     int localactiveset = activeset;

//     if (localactiveset == -1) {
//       return false;
//     }
//     if (vertexsets[localactiveset]->next(nextschedpos[cpuid], cpuid) == false) {
//       return false;
//     }
//     size_t tmp = nextschedpos[cpuid];
//     ret_task = update_task(tmp, schedupdate);
//     return true;
//   }

//   void set_scheduler::init() {
//     ss_set_type allvertices;
//     for (size_t i = 0; i < g.num_vertices(); ++i) {
//       ss_insert(allvertices, i);
//     }
//     static ss_set_type unused;
//     rootset.rebuild(NULL, allvertices);
//     rootset.resolve_event_handlers();
//     rootset.set_as_root_set();
//     rootsetevents = rootset.all_dependent_events();
//     logger(LOG_INFO, "Set Scheduler Initialization Complete");
//     logger(LOG_INFO, "Set Scheduler Active Event Set: %d", rootsetevents);
//     startschedule->wait();
//   }

  
//   void set_scheduler::wait() {
  
//     schedposlock.lock();
//     while(1) {  
//       if ((activeset == -1 || vertexsets[activeset]->size() == 0) && 
//           pendingtaskctr.value == 0 && 
//           executedtaskctr.value == finishedtaskctr.value &&
//           (cureplan == NULL || cureplan->done())) {
//         break;
//       }
//       schedposcond.wait(schedposlock);
//     }
//     schedposlock.unlock();

//   }

//   void set_scheduler::begin_schedule(schedule_function s) {
//     executing = false;
//     complete = false;
//     size_t nschedulethreads = 1;
//     //  startbarrier = new graphlab::barrier(nschedulethreads);
//     worker.resize(nschedulethreads);
//     for (size_t i = 0; i < nschedulethreads; ++i) {
//       worker[i].func = s;
//       worker[i].parent = this;
//       scheduling_thread.launch(&worker[i]);
//     }
//     startschedule->wait();
//   }


// }; //end of namespace
