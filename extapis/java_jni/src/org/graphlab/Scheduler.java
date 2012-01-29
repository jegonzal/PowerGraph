package org.graphlab;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import org.graphlab.util.StringUtils;

/**
 * Mirrors scheduler_list.hpp
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 * @see <a href="http://graphlab.org/doxygen/html/Schedulers.html">Schedulers</a>
 */
public enum Scheduler {

  // TODO: chromatic scheduler

  /**
   * The fifo scheduler executes tasks in the classical first in first out
   * order. When update functions generate tasks they are added to the back of
   * the fifo queue. Because there is a single central queue the fifo scheduler
   * can become a synchronizing bottleneck for algorithms with relatively light
   * update functions.
   */
  FIFO("fifo"),
  
  /**
   * This scheduler maintains a shared FIFO queue of FIFO queues. 
   * Each thread maintains its own smaller in and out queues. When a 
   * threads out queue is too large (greater than "queuesize") then 
   * the thread puts its out queue at the end of the master queue.
   */
  QUEUED_FIFO("queued_fifo"),

  /**
   * The priority scheduler maintains a single priority scheduling queue. The
   * task with highest priority is always executed next. If add_task() is
   * invoked on an already present task the existing task's priority is set to
   * the max of the two tasks.
   */
  PRIORITY("priority"),

  /**
   * This scheduler loops over vertices executing a task if one is associated
   * with the vertex. Each vertex maintains its own local queue. This scheduler
   * has the least possible overhead.
   */
  SWEEP("sweep");

  private final String mName;
  private Map<String, String> mOptions;

  private Scheduler(String name) {
    mName = name;
    mOptions = new HashMap<String, String>();
  }
  
  public void setOption(String key, String value){
    mOptions.put(key, value);
  }
  
  /**
   * @return just the name of the scheduler.
   */
  public String type(){
    return mName;
  }

  /*
   * Creates the command-line options string by concatenating the
   * scheduler options.
   * (non-Javadoc)
   * @see java.lang.Enum#toString()
   */
  @Override
  public String toString() {
    
    // if no options, just return name
    if (mOptions.size() <= 0) return mName;
    
    // otherwise, prepare options string
    List<String> options = new LinkedList<String>();
    for (Entry<String, String> option : mOptions.entrySet()){
      options.add(option.getKey() + "=" + option.getValue());
    }
    
    return mName + "[" + StringUtils.join(options, ",") + "]";
    
  }
  
}
