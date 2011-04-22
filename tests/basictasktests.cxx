#include <stdint.h>
#include <string>
#include <cxxtest/TestSuite.h>

#include <graph/graph.hpp>
#include <graph/iscope.hpp>
#include <util/synchronized_queue.hpp>
#include <tasks/update_task.hpp>
#include <schedulers/ischeduler.hpp>
#include <parallel/pthread_tools.hpp>
#include <schedulers/support/buffered_scheduler_callback.hpp>
#include <schedulers/support/vertex_task_set.hpp>
#include <graphlab/logger/logger.hpp>

#include <util/mutable_queue.hpp>

using namespace graphlab;

int count = 0;

void very_dummy_update(iscope &scope,
                       ischeduler_callback &scheduler) {
    count++;
  }





class BasicTaskTestSuite: public CxxTest::TestSuite {
   public:
        
  void test_synchronized_fifo_taskqueue(void) {
    synchronized_queue<update_task> * queue = new synchronized_queue<update_task>();
    TS_TRACE("Created queue");
    
    int i=0;
  
    /* Fill queue */
    for(i=0; i<100; i++) {
      update_task task(i, very_dummy_update);
      queue->push(task);
      TS_ASSERT_EQUALS(i+1, queue->size());
    }
    
    /* Take stuff out from queue */
    for(i=0; i<100; i++) {
      update_task  task = queue->pop();
      TS_ASSERT_EQUALS(task.vertex_id(), i);
      TS_ASSERT_EQUALS(100-i-1, queue->size());
    }
    TS_ASSERT_EQUALS(0, queue->size());
    TS_ASSERT_EQUALS(0, queue->size());

    delete queue;
  }
  
  
  /*** PRIORITY QUEUE ***/
  void test_mutable_queue_basics() {
    mutable_queue<update_task, float, std::less<update_task>, update_task_hash> priority_queue;
    
    /* Insert some elements and check they come out in expected order. 
       Delibrately use arbitary vertex indices so a queue that would sort
       tasks by vertex id would not pass the test. */
    priority_queue.push(update_task(10, very_dummy_update), 0.2);
    priority_queue.push(update_task(12, very_dummy_update), 0.1);
    priority_queue.push(update_task(13, very_dummy_update), 0.15);
    priority_queue.push(update_task(9, very_dummy_update), 1.0);
    
    TS_TRACE("Created priority queue, now popping");

    TS_ASSERT_EQUALS(priority_queue.pop().first.vertex_id(), 9);
    TS_ASSERT_EQUALS(priority_queue.pop().first.vertex_id(), 10);
    TS_ASSERT_EQUALS(priority_queue.pop().first.vertex_id(), 13);
    TS_ASSERT_EQUALS(priority_queue.pop().first.vertex_id(), 12);
    
    TS_ASSERT_EQUALS(priority_queue.empty(), true);
  }
  
  void test_mutable_queue_task_update_max() {
    mutable_queue<update_task, float, std::less<update_task>, update_task_hash> priority_queue;

    update_task t1(10, very_dummy_update);
    update_task t2(12, very_dummy_update);
    update_task t3(13, very_dummy_update);
    update_task t4(9, very_dummy_update);
    
    priority_queue.push(t1, 0.2f);
    priority_queue.push(t2, 0.1f);
    priority_queue.push(t3, 0.15f);
    priority_queue.push(t4, 1.0f);
    
    /* Check top is right */
    TS_ASSERT_EQUALS(priority_queue.top().first.vertex_id(), 9);
    
    /* Promote one */
    priority_queue.updating_insert_max(t3, 1.2);
    TS_ASSERT_EQUALS(priority_queue.pop().first, t3);
    
    /* Second best */
    TS_ASSERT_EQUALS(priority_queue.top().first, t4);
    
    /* Check that if we have lower new priority, it does not affect anything */
    priority_queue.updating_insert_max(t1, 0.05);
    TS_ASSERT_EQUALS(priority_queue.get(t1) , 0.2f);

    /* Check rests */
    TS_ASSERT_EQUALS(priority_queue.pop().first.vertex_id(), 9);
    TS_ASSERT_EQUALS(priority_queue.pop().first.vertex_id(), 10);
    TS_ASSERT_EQUALS(priority_queue.pop().first.vertex_id(), 12);
    
    TS_ASSERT_EQUALS(priority_queue.empty(), true);

  }
  
  void test_cumulative_priority() {
    mutable_queue<update_task, float, std::less<update_task>, update_task_hash> priority_queue;
    
    update_task t1(10, very_dummy_update);
    update_task t2(12, very_dummy_update);
    update_task t3(13, very_dummy_update);
    update_task t4(9, very_dummy_update);
    priority_queue.push(t1, 0.2f);
    priority_queue.push(t2, 0.1f);
    priority_queue.push(t3, 0.15f);
    priority_queue.push(t4, 1.0f);
      
    priority_queue.updating_insert_cumulative(t4, 0.5f);
    TS_ASSERT_EQUALS(priority_queue.get(t4), 1.5f);
  }
  
  
};



