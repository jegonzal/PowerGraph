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

#ifndef GRAPHLAB_DISTRIBUTED_SCHEDULER_LIST_HPP
#define GRAPHLAB_DISTRIBUTED_SCHEDULER_LIST_HPP
#include <string>
#include <vector>
#include <iostream>
#include <boost/preprocessor.hpp>

#define __DISTRIBUTED_SCHEDULER_LIST__                                              \
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
    "priority scheduler, but less predictable, and much faster."))      


#include <graphlab/schedulers/fifo_scheduler.hpp>
#include <graphlab/schedulers/priority_scheduler.hpp>
#include <graphlab/schedulers/sweep_scheduler.hpp>
#include <graphlab/schedulers/multiqueue_fifo_scheduler.hpp>
#include <graphlab/schedulers/multiqueue_priority_scheduler.hpp>

namespace graphlab {
  /// get all the scheduler names
  std::vector<std::string> get_distributed_scheduler_names();

  /// get all the scheduler names concated into a string
  std::string get_distributed_scheduler_names_str();

  /// Display the scheduler options for a particular scheduler
  void print_distributed_scheduler_info(std::string s, std::ostream &out);
}

#endif

