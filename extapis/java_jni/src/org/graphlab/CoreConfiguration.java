package org.graphlab;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import org.graphlab.util.StringUtils;

public class CoreConfiguration {

  private Map<String, String> mOptions = new HashMap<String, String>();
  
  public CoreConfiguration(){}
  
  public CoreConfiguration(Scheduler scheduler, Scope scope, int ncpus){
    setScheduler(scheduler);
    setScopeType(scope);
    setNCpus(ncpus);
  }
  
  public void setScheduler(Scheduler scheduler){
    if (null == scheduler)
      throw new IllegalArgumentException ("scheduler must not be null.");
    mOptions.put("--scheduler", scheduler.toString());
  }
  
  public void setScopeType(Scope scope){
    if (null == scope)
      throw new IllegalArgumentException ("scope must not be null.");
    mOptions.put("--scope", scope.toString());
  }
  
  public void setNCpus(int n){
    if (0 >= n)
      throw new IllegalArgumentException ("n must be positive.");
    mOptions.put("--ncpus", n + "");
  }
  
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
