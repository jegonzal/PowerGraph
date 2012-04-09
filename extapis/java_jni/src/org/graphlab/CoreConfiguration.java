package org.graphlab;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import org.graphlab.util.StringUtils;

/**
 * Set of configuration options for {@link org.graphlab.Core}.
 * 
 * <h3>Usage</h3>
 * <p>To use, create an instance of this class and call the necessary <code>set*</code>
 * methods to specify your desired configuration. Then, pass the instance to
 * {@link Core#Core(CoreConfiguration)}.</p>
 * 
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 */
public class CoreConfiguration {

  /** Options Map - options keys to options values */
  private Map<String, String> mOptions = new HashMap<String, String>();
  
  /**
   * Creates a configuration with the default options.
   * <dl>
   *  <dt>Scheduler<dt>
   *  <dd>{@link org.graphlab.Scheduler#SWEEP}</dd>
   *  <dt>Scope</dt>
   *  <dd>{@link org.graphlab.Scope#EDGE}</dd>
   *  <dt>Number of CPUs</dt>
   *  <dd>Number of CPUs available on the system</dd>
   * </dl>
   */
  public CoreConfiguration(){}
  
  /**
   * Creates a configuration with the given scheduler, scope, and number of
   * CPUs.
   * 
   * @param scheduler
   *          one of {@link org.graphlab.Scheduler}
   * @param scope
   *          one of {@link org.graphlab.Scope}
   * @param ncpus
   *          number of CPUs 
   */
  public CoreConfiguration(Scheduler scheduler, Scope scope, int ncpus){
    setScheduler(scheduler);
    setScopeType(scope);
    setNCpus(ncpus);
  }
  
  /**
   * Specifies the scheduler to use.
   * 
   * @param scheduler
   *          one of {@link org.graphlab.Scheduler}
   * @throws NullPointerException
   *          if <tt>scheduler</tt> was null.
   */
  public void setScheduler(Scheduler scheduler){
    if (null == scheduler)
      throw new NullPointerException ("scheduler must not be null.");
    mOptions.put("--scheduler", scheduler.toString());
  }
  
  /**
   * Specifies the scope to use.
   * 
   * @param scope
   *          one of {@link org.graphlab.Scope}
   * @throws NullPointerException
   *          if <tt>scope</tt> was null.
   */
  public void setScopeType(Scope scope){
    if (null == scope)
      throw new NullPointerException ("scope must not be null.");
    mOptions.put("--scope", scope.toString());
  }
  
  /**
   * Specifies the number of CPUs GraphLab should use. This roughly corresponds to
   * the number of threads that will be spawned by the scheduler.
   * 
   * @param n
   *      number of CPUs
   * @throws IllegalArgumentException
   *          if <tt>n</tt> was negative.
   */
  public void setNCpus(int n){
    if (0 >= n)
      throw new IllegalArgumentException ("n must be positive.");
    mOptions.put("--ncpus", n + "");
  }
  
  /*
   * Creates the command line string by concatenating option keys
   * and values.
   * (non-Javadoc)
   * @see java.lang.Object#toString()
   */
  @Override
  public String toString(){
    
    // concatenate options into a string of command line args
    List<String> options = new LinkedList<String>();
    for (Entry<String, String> option : mOptions.entrySet()){
      options.add(option.getKey() + "=" + option.getValue());
    }
    
    return StringUtils.join(options, " ");
    
  }
  
}
