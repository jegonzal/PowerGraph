
GraphLab Java and Jython API

== BUILDING ==

To build with Ant, use:
$>  cd extapis/java_jni
$>  ant
   
To create a distribution and generate javadocs, use:
$>  ant dist

You can also simply use your favorite IDE, and 
link the src/ and lib/ directories.


== JAVA NATIVE INTERFACE LIBRARY ==

Java API interfaces with GraphLab using Java Native Interface (JNI).
Precompiled versions for some architectures will be made available (later)
in lib/ directory.
Currently we support x86_64 Linux and Mac OS X libraries.

To build your own version of the JNI shared library from C++ GraphLab
distribution:
  ./configure
	cd release
	make

This will create the library (.so suffix for Linux, .dylib for Mac) in
release/src/graphlab/jni. You may also choose to build debug or profile
versions by cd'ing into the corresponding directory before calling make.

== RUNNING PROGRAMS ==

Simple demos may be executed using:
$> ant <example name>

For example, the following executes the ShortestPath demo:
$> ant ShortestPath

Demos may be explored in org.graphlab.demo package. To run the demos
without ant, you will need to specify java.library.path. E.g.
$> java \
    -classpath bin:lib/log4j-1.2.16.jar \
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

