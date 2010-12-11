package graphlab;

import java.util.List;

/**
 * Interface to Scheduler, passed to UpdateFunction. Allows scheduling
 * neighboring vertices.
 * @author akyrola
 */
public interface Scheduler {
    
    /** 
      * Add task to vertex with given id. <b>Note:</b> UpdateFunction instance
      * has to have been provided to GraphLabJNIWrapper in constructor. 
      * Use "this" to schedule the current function.
      * @param vertex vertex id
      * @param function Update function instance. 
      */
    void addTask(int vertex, UpdateFunction function);
    void addTask(int vertex, UpdateFunction function, float priority);

    void addTask(UpdateFunction function, int ...v);
    
    /**
      * Schedule current update function to list of vertices.
      * @param v list of vertex ids
      */
    public void addTask(int... v);
    public void addTask(List<Integer> vv);
    
    /**
      * Schedule all outbound neighbors.
      * @param scope Scope passed to the update function.
      */
    public void addTaskToOutbound(Scope scope);    
    
}
