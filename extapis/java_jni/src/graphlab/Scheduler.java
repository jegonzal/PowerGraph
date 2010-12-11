package graphlab;

import java.util.List;

/**
 * @author akyrola
 */
public interface Scheduler {

    void addTask(int vertex, UpdateFunction function);
    void addTask(int vertex, UpdateFunction function, float priority);

    void addTask(UpdateFunction function, int ...v);

    public void addTask(int... v);
    public void addTask(List<Integer> vv);
    public void addTaskToOutbound(Scope scope);    
    
}
