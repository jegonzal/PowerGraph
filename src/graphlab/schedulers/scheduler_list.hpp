/*
This file is part of GraphLab.

GraphLab is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as 
published by the Free Software Foundation, either version 3 of 
the License, or (at your option) any later version.

GraphLab is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public 
License along with GraphLab.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef GRAPHLAB_SCHEDULER_LIST_HPP
#define GRAPHLAB_SCHEDULER_LIST_HPP
#include <string>
#include <vector>
#include <iostream>
#include <boost/preprocessor.hpp>

#define __SCHEDULER_LIST__                                              \
  (("chromatic", chromatic_scheduler,                                   \
    "a scheduler which performs #iterations sweeps of the graph using " \
    "a graph color ordering."))                                         \
  (("sweep", sweep_scheduler,                                           \
    "very fast dynamic scheduler. Scans all vertices in sequence, "     \
    "running all update tasks on each vertex evaluated."))              \
  (("fifo", fifo_scheduler,                                             \
    "Standard FIFO task queue, poor parallelism, but task evaluation "  \
    "sequence is highly predictable. Useful for debugging and testing.")) \
  (("priority", priority_scheduler,                                     \
    "Standard Priority queue, poor parallelism, but task evaluation "   \
    "sequence is highly predictable. Useful for debugging"))            \
  (("multiqueue_fifo", multiqueue_fifo_scheduler,                       \
    "One or more FIFO task queues is assigned to each processor, "      \
    "where the queues are stochastically load balanced. Like the "      \
    "fifo scheduler, but less predictable, and much faster."))          \
  (("multiqueue_priority", multiqueue_priority_scheduler,               \
    "One or more Priority task queues is assigned to each processor, "  \
    "where the queues are stochastically load balanced. Like the "      \
    "priority scheduler, but less predictable, and much faster."))      \
  (("splash", splash_scheduler,                                         \
    "Similar to the priority queue scheduler, but allows for only one " \
    "update function. Updates are evaluted in a \"splash\" ordering"))  \
  (("round_robin", round_robin_scheduler,                               \
    "Loops over a sequence of tasks repeatedly for # iterations."))     \
  (("clustered_priority", clustered_priority_scheduler,                 \
    "Like the priority scheduler, but groups vertices into clusters "   \
    "where the entire cluster has a single priority"))                  \
  (("sampling", sampling_scheduler,                                     \
    "A scheduler which samples vertices to update based on a "          \
    "multinomial probability which can be updated dynamically."))


#include <graphlab/schedulers/fifo_scheduler.hpp>
#include <graphlab/schedulers/priority_scheduler.hpp>
#include <graphlab/schedulers/sampling_scheduler.hpp>
#include <graphlab/schedulers/round_robin_scheduler.hpp>
#include <graphlab/schedulers/chromatic_scheduler.hpp>
#include <graphlab/schedulers/sweep_scheduler.hpp>
#include <graphlab/schedulers/splash_scheduler.hpp>
#include <graphlab/schedulers/multiqueue_fifo_scheduler.hpp>
#include <graphlab/schedulers/multiqueue_priority_scheduler.hpp>
#include <graphlab/schedulers/clustered_priority_scheduler.hpp>
#include <graphlab/graph/graph.hpp>

namespace graphlab {
  /// get all the scheduler names
  std::vector<std::string> get_scheduler_names();

  /// get all the scheduler names concated into a string
  std::string get_scheduler_names_str();

  /// Display the scheduler options for a particular scheduler
  void print_scheduler_info(std::string s, std::ostream &out);
}

#endif

