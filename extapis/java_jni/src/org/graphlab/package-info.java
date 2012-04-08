/**
 * Provides basic GraphLab classes.
 * <p>More utilities may be found in {@link org.graphlab.data}
 * and {@link org.graphlab.toolkits}, but this package has the minimum classes you need to
 * get started.</p>
 * 
 * <h2>Introduction</h2>
 * <p>Welcome to GraphLab Java 2. This page will give you a brief introduction; if you would like
 * to have a more detailed tutorial, please refer to
 * <a href='http://graphlab.org/' title='GraphLab'>the official website</a>.</p>
 * 
 * <h3>GraphLab</h3>
 * <p>GraphLab is an abstraction for graph algorithms. Much akin to functional programming, it
 * allows you to define functions (which we
 * refer to as {@link org.graphlab.Updater updaters}) to be executed on vertices and edges. These
 * functions are scheduled using a scheduler, which takes the task of handling parallelism out of
 * your hands. What you just downloaded (and hopefully compiled!) is a framework that seeks to
 * exemplify and motivate this abstraction.</p>
 * 
 * <p>For example, suppose you wish to write a program to increment the value stored at every vertex
 * in a graph. In GraphLab, the pseudo-code will be written as:</p>
<pre>
  def update(vertex)
    read value from vertex
    increment value
    store value in vertex
  end
  
  schedule_all(graph, update)
</pre>
 * <p>GraphLab allows you to choose between several {@link org.graphlab.Scheduler schedulers} and specify
 * different consistencies. For example, if you wish to ensure safe read/write access to the current vertex
 * as well as adjacent edges, you may use the {@link org.graphlab.Scope#EDGE edge consistency}. If you only
 * need read/write access to the current vertex, just opt for {@link org.graphlab.Scope#VERTEX vertex consistency}.
 * In general, it's recommended to use the least consistency you need, so that the scheduler may exploit more
 * opportunities for parallelism.</p>
 * 
 * <p>In general, applications written using the GraphLab abstraction look like:</p>
<pre>
  def update(vertex)
    gather data
    compute data
    save data
    reschedule neighbors if necessary
  end
  
  schedule(update)
</pre>
 * <h3>Installation</h3>
 * <p>First, obtain a copy of GraphLab from <a href='http://graphlab.org/download.html'>the official site</a>.
 * Ensure that you have Java 1.6 and above before running the compilation script. The script will compile
 * the GraphLab shared library as well as the JNI interface.</p>
 * 
 * <p>To see an example, navigate to <code>graphlabapi/extapis/java_jni</code> and execute <code>ant &lt;demo-name&gt; -Dfile=&lt;input-file&gt;</code>.
 * Currently available demos include:</p>
 * <ul>
 *  <li>PageRank</li>
 *  <li>Alternating Least Squares</li>
 *  <li>Graph Coloring</li>
 *  <li>Shortest Path</li>
 * </ul>
 * 
 * <h3>Making an App</h3>
 * <p>To create a Java app, execute the following commands:</p>
<pre>
 $> cd /path/to/graphlabapi/release
 $> make java_app
</pre>
 * <p>An Eclipse project with all the necessary libraries will be created for you at
 * <code>graphlabapi/release/extapis/java_jni/my-app</code>. If you use Eclipse, follow the following instructions:</p>
 * <ol>
 *  <li>Go to File > Import ...</li>
 *  <li>In the dialog box that appears, choose "Existing Projects into Workspace"</li>
 *  <li>Select <code>graphlabapi/release/extapis/java_jni/my-app</code> as the root directory</li>
 * </ol>
 * <p>You may use any IDE you desire, but remember to add the jar files in the <code>lib</code>
 * directory to your build path.</p>
 * <p>Once that is done, you should see the following directory structure in the <code>GraphLab-App</code> project:</p>
 * <p>
 *  <img src='doc-files/graphlab-1.png' alt='GraphLab-App Directory Structure' title='GraphLab-App Directory Structure'/>
 * </p>
 * <ul>
 *  <li><code>src</code> contains an example program</li>
 *  <li><code>docs</code> contains this documentation</li>
 *  <li><code>lib</code> contains the necessary libraries</li>
 * </ul>
 * <p>Every time you create a new Run Configuration, remember to add "<code>-Djava.library.path=${resource_loc:/GraphLab-App/lib}</code>"
 * to the VM arguments.</p>
 */
package org.graphlab;

