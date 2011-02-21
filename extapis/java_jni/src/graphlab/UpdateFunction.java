package graphlab;

/**
 * Update function class. Override to implement your own.
 * @author akyrola
 */
public abstract class UpdateFunction {

   public int functionId;

   abstract public void update(Scope scope, Scheduler scheduler);

}
