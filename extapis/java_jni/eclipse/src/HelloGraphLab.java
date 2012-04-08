import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.graphlab.Context;
import org.graphlab.Core;
import org.graphlab.Updater;
import org.graphlab.data.ScalarVertex;
import org.jgrapht.DirectedGraph;
import org.jgrapht.graph.DefaultDirectedWeightedGraph;
import org.jgrapht.graph.DefaultWeightedEdge;

/**
 * Say HELLO! to GraphLab.
 * 
 * Before running, add "-Djava.library.path=${resource_loc:/GraphLab-App/lib}" to the
 * VM arguments in the Run Configuration.
 * 
 * @author Jiunn Haur Lim
 */
public class HelloGraphLab {
  
  private static Core core;
  private static DirectedGraph<ScalarVertex, DefaultWeightedEdge> graph;

  public static void main(String[] args) throws Core.CoreException {
    
    BasicConfigurator.configure();
    Logger.getLogger(Core.class).setLevel(Level.OFF);
    
    core = new Core();
    // JGraphT
    graph = new DefaultDirectedWeightedGraph<ScalarVertex, DefaultWeightedEdge>(DefaultWeightedEdge.class);
    int[] hello = {'G', 'D', 'K', 'K', 'N', ' '};
    
    for (int i=0; i<hello.length; i++)
      graph.addVertex(new ScalarVertex(i, hello[i]));
    
    core.setGraph(graph);
    core.scheduleAll(new HelloUpdater());
    core.start();
    
    // don't forget to destroy!
    core.destroy();
    
    for(ScalarVertex v : graph.vertexSet())
      System.out.print((char) v.value());
    
  }
  
  private static class HelloUpdater
    extends Updater<ScalarVertex, DefaultWeightedEdge, HelloUpdater>{
    
    @Override
    public void update(Context context, ScalarVertex vertex) {
      // increment value stored in vertex
      int c = (int) vertex.value();
      vertex.setValue(c+1);
    }
    
    @Override
    protected HelloUpdater clone() {
      return new HelloUpdater();
    }
    
  }
  
}
