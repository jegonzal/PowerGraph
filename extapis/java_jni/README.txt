
GraphLab Java and Jython API

== BUILDING ==

To build with Ant, just type
   ant

You can also simply use your favorite IDE, and 
link the src/ and lib/ directories.


== JAVA NATIVE INTERFACE LIBRARY ==

Java API interfaces with GraphLab using Java Native Interface (JNI).
Precompiled versions for some architectures will be made available (later)
in native/ -directory.
Currently we support x86_64 Linux and Mac OS X libraries.

To build your own version of the JNI shared library from C++ GraphLab
distribution:
  ./configure
	cd release
  cmake <... run the cmake command from the ./configure output ...>
	make

This will create the library (.so suffix for Linux, .dylib for Mac) in
release/src/graphlab/jni. You may also choose to build debug or profile
versions by cd'ing into the corresponding directory before calling cmake.

== RUNNING PROGRAMS ==

See directory /examples for some demo apps.

Remember to set the Java library path:
	 java -Djava.library.path=<path to graphlab>/release/src/graphlab/jni
	 
	 
	 


Updated December 10, 2010.
Aapo Kyrola, akyrola@cs.cmu.edu

Updated (see mercurial timestamp)
Jiunn Haur Lim, jiunnhal@cmu.edu

