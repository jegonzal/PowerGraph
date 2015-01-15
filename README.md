Graphlab 
==========


License
-------

GraphLab is free software licensed under the Apache 2.0 License. See
license/LICENSE.txt for details.


Introduction
------------

GraphLab is a graph-based, high performance, distributed computation framework
written in C++. The GraphLab project started in 2009 to develop a new parallel computation abstraction 
tailored to machine learning. GraphLab 1.0 represents our first shared memory design,
and in GraphLab 2.1, we completely redesigned the framework to target the distributed environment addressing
the difficulties with real world power-law graphs, achieving unparalleled performance.
In GraphLab 2.2, we introduce the Warp System which provides a new
flexible, distributed architecture around fine-grained user-mode threading (fibers).
The new Warp system will allow us to easily extend the abstraction to cover new ground
while improving useability, and will also allow us to realize new system optimizations that 
are not available in the past.

GraphLab Features:

* **Unified multicore/distributed API:**       write once run anywhere 

* **Tuned for performance:** optimized C++ execution engine leverages extensive multi-threading and asynchronous IO 

* **Scalable:**              Run on large cluster deployments by intelligently placing data and computation 

* **HDFS Integration:**      Access your data directly from HDFS 

* **Powerful Machine Learning Toolkits:**     Tackle challenging machine learning problems with ease


For more details on the GraphLab see http://graphlab.org, including
documentation, tutorial, etc.



Dependencies
------------

GraphLab now automatically satisfied most dependencies. 
There are however, a few dependencies which we cannot easily satisfy:

* On OS X: g++ (>= 4.2) or clang (>= 3.0) [Required]
  +  Required for compiling GraphLab.

* On Linux: g++ (>= 4.3) or clang (>= 3.0) [Required]
  +  Required for compiling GraphLab.

* *nix build tools: patch, make [Required]
   +  Should come with most Mac/Linux systems by default. Recent Ubuntu version will require to install the build-essential package.

* zlib [Required]
   +   Comes with most Mac/Linux systems by default. Recent Ubuntu version will require the zlib1g-dev package.

* Open MPI or MPICH2 [Strongly Recommended]
   + Required for running GraphLab distributed. 

* JDK 6 or greater [Optional]
   + Required for HDFS support 


    
### Satisfying Dependencies on Mac OS X

Installing XCode with the command line tools (in XCode 4.3 you have to do this
manually in the XCode Preferences -> Download pane), satisfies all of these
dependencies.  



### Satisfying Dependencies on Ubuntu

All the dependencies can be satisfied from the repository:

    apt-get install gcc g++ build-essential libopenmpi-dev default-jdk cmake zlib1g-dev



Compiling
---------

    ./configure

In the graphlab directory, will create two sub-directories, release/ and
debug/ . cd into either of these directories and running make will build the
release or the debug versions respectively. Note that this will compile all of
GraphLab, including all toolkits. Since some toolkits require additional
dependencies (for instance, the Computer Vision toolkit needs OpenCV), this
will also download and build all optional dependencies.

We recommend using make’s parallel build feature to accelerate the compilation
process. For instance:

    make -j 4

will perform up to 4 build tasks in parallel. When building in release/ mode,
GraphLab does require a large amount of memory to compile with the
heaviest toolkit requiring 1GB of RAM. Where K is the amount of memory you
have on your machine in GB, we recommend not exceeding make -j K

Alternatively, if you know exactly which toolkit you want to build, cd into the
toolkit’s sub-directory and running make, will be significantly faster as it
will only download the minimal set of dependencies for that toolkit. For
instance:

    cd release/toolkits/graph_analytics
    make -j4

will build only the Graph Analytics toolkit and will not need to obtain OpenCV,
Eigen, etc used by the other toolkits.




Writing Your Own Apps
---------------------

There are two ways to write your own apps.

1: To work in the GraphLab source tree,    (recommended)
2: Install and link against Graphlab       (not recommended)



### 1:  Working in the GraphLab Source Tree

This is the best option if you just want to try using GraphLab quickly. GraphLab
uses the CMake build system which enables you to quickly create
a c++ project without having to write complicated Makefiles. 

1: Create your own sub-directory in the apps/ directory. for example apps/my_app
   
2: Create a CMakeLists.txt in apps/my_app containing the following lines:

    project(GraphLab) 
    add_graphlab_executable(my_app [List of cpp files space seperated]) 

  Substituting the right values into the square brackets. For instance:

    project(GraphLab) 
    add_graphlab_executable(my_app my_app.cpp) 

4: Running "make" in the apps/ directory of any of the build directories 
should compile your app. If your app does not show up, try running

    cd [the GraphLab API directory]
    touch apps/CMakeLists.txt

and try again.



### 2: Installing and Linking Against GraphLab

To install graphlab and use GraphLab this way will require your system
to completely satisfy all remaining dependencies, which GraphLab normally 
builds automatically. This path is not extensively tested and is 
**not recommended**

You will require the following additional dependencies
 - libevent (>=2.0.18)
 - libjson (>=7.6.0)
 - libboost (>=1.53)
 - libhdfs (required for HDFS support)
 - tcmalloc (optional)

Follow the instructions in the [Compiling] section to build the release/ 
version of the library. Then cd into the release/ build directory and 
run make install . This will install the following:

* include/graphlab.hpp
 +   The primary GraphLab header 
*  include/graphlab/...
 +   The folder containing the headers for the rest of the GraphLab library 
*  lib/libgraphlab.a
 +   The GraphLab static library.
    
Once you have installed GraphLab you can compile your program by running:

    g++ -O3 -pthread -lzookeeper_mt -lzookeeper_st -lboost_context -lz -ltcmalloc -levent -levent_pthreads -ljson -lboost_filesystem -lboost_program_options -lboost_system -lboost_iostreams -lboost_date_time -lhdfs -lgraphlab hello_world.cpp 
    
If you have compiled with MPI support, you will also need

   -lmpi -lmpi++ 
  
  

