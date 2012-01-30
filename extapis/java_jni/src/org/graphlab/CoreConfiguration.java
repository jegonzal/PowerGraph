package org.graphlab;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import org.graphlab.util.StringUtils;

/**
 * Set of configuration options for {@link org.graphlab.Core}.
 * Pass an instance of this class to {@link Core#Core(CoreConfiguration)}.
 * 
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 */
public class CoreConfiguration {

  /** Options Map - options keys to options values */
  private Map<String, String> mOptions = new HashMap<String, String>();
  
  /** Default options */
  public CoreConfiguration(){}
  
  /**
   * Creates a configuration with the given scheduler, scope, and number of
   * CPUs.
   * 
   * @param scheduler
   * @param scope
   * @param ncpus
   */
  public CoreConfiguration(Scheduler scheduler, Scope scope, int ncpus){
    setScheduler(scheduler);
    setScopeType(scope);
    setNCpus(ncpus);
  }
  
  /**
   * @param scheduler
   * @throws NullPointerException
   *          if <tt>scheduler</tt> was null.
   */
  public void setScheduler(Scheduler scheduler){
    if (null == scheduler)
      throw new NullPointerException ("scheduler must not be null.");
    mOptions.put("--scheduler", scheduler.toString());
  }
  
  /**
   * @param scope
   * @throws NullPointerException
   *          if <tt>scope</tt> was null.
   */
  public void setScopeType(Scope scope){
    if (null == scope)
      throw new NullPointerException ("scope must not be null.");
    mOptions.put("--scope", scope.toString());
  }
  
  /**
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
   * Creates the command line string by concatentating option keys
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
