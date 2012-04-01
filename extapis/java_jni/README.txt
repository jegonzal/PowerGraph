
GraphLab Java API

== BUILDING ==

GraphLab Java requires Ant, which you may obtain from http://ant.apache.org/

To build the examples, use the following commands:
$> ./configure
$> cd release
$ make graphlabjni

To use GraphLab in a project, use the following commands:
$> ./configure
$> cd release
$> make java_app

== RUNNING PROGRAMS ==

Simple demos may be executed using:
$> ant <example name>

For example, the following executes the ShortestPath demo:
$> ant ShortestPath

Demos may be explored in org.graphlab.demo package. To run the demos
without ant, you will need to specify java.library.path. E.g.
$> java \
    -classpath bin:lib/* \
    -Djava.library.path=../../release/src/graphlab/jni \
    org.graphlab.demo.Coloring test-graphs/toy.tsv
	 
== LOGGING ==

The library uses log4j for logging. Refer to
http://logging.apache.org/log4j/1.2/manual.html
for instructions on how to turn off or control the log output.


Updated December 10, 2010.
Aapo Kyrola, akyrola@cs.cmu.edu

Updated (see mercurial timestamp)
Jiunn Haur Lim, jiunnhal@cmu.edu

