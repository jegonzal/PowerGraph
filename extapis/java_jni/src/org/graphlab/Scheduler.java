package org.graphlab;

/**
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 * @see <a href="http://graphlab.org/doxygen/html/Schedulers.html">Schedulers</a>
 */
public enum Scheduler {

  // TODO: chromatic scheduler

  /** This scheduler executes repeated cycles through all vertices. */
  ROUND_ROBIN("round_robin"),

  /**
   * The fifo scheduler executes tasks in the classical first in first out
   * order. When update functions generate tasks they are added to the back of
   * the fifo queue. Because there is a single central queue the fifo scheduler
   * can become a synchronizing bottleneck for algorithms with relatively light
   * update functions.
   */
  FIFO("fifo"),

  /**
   * The Multiqueue fifo scheduler is like the fifo scheduler but instead of
   * using a single queue multiple queues (2 x ncpus) are used and tasks are
   * added to queues using a randomized balancing strategy. Each processor only
   * draws from its own pair of queue.
   */
  MULTI_QUEUE_FIFO("multiqueue_fifo"),

  /**
   * The priority scheduler maintains a single priority scheduling queue. The
   * task with highest priority is always executed next. If add_task() is
   * invoked on an already present task the existing task's priority is set to
   * the max of the two tasks.
   */
  PRIORITY("priority"),

  /**
   * Same as the priority scheduler except multiple queues are maintained.
   */
  MULTI_QUEUE_PRIORITY("multiqueue_priority"),

  /**
   * The clustered priority schedule maintains a priority queue over blocks of
   * vertices. This schedule begins by partitioning the graph using the
   * partitioning method into blocks which each contain vert_per_part vertices.
   * The priorities are maintained over blocks but within each block a sweep
   * schedule is run.
   */
  CLUSTERED_PRIORITY("clustered_priority"),

  /**
   * This scheduler loops over vertices executing a task if one is associated
   * with the vertex. Each vertex maintains its own local queue. This scheduler
   * has the least possible overhead.
   */
  SWEEP("sweep"),

  /**
   * We currently only provide the Splash scheduler which grows small spanning
   * trees for each cpu and then sequential executes the update function in a
   * forward backward fashion. This scheduler is loosely based on the Splash BP
   * algorithm by Gonzalez et al. AISTATS 2009.
   */
  SPLASH("splash");

  private final String mStr;

  private Scheduler(String str) {
    mStr = str;
  }

  @Override
  public String toString() {
    return mStr;
  }

}
