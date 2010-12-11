
GraphLab Java and Jython API

== BUILDING ==

To build with Ant, just type
   ant

You can also simply use your favorite IDE, and 
include the src/ and lib/ directories.


== JAVA NATIVE INTERFACE LIBRARY ==

Java API interfaces with GraphLab using Java Native Interface (JNI).
Precompiled versions for some architectures are in native/ -directory.
Currently we support x86_64 Linux and Mac OS X libraries.

To build your own version of the JNI shared library from C++ GraphLab
distribution:
	cd release/src/graphlab/jni
	make

This will create the library (.so suffix for Linux, .jnilib for Mac).
Copy this to directory /native.

== RUNNING PROGRAMS ==

See directory /examples for some demo apps.

Remember to set the Java library path to directory /native:
	 java -Djava.library.path=<DIR>/native [...]




Updated, December 10, 2010.
Aapo Kyrola, akyrola@cs.cmu.edu

