/**
\page RPC GraphLab RPC

GraphLab RPC primary design goal was to provide a convenient and easy to use
asynchronous communication system between \b identical binaries running
on different machines over a distributed network. It therefore provides 
MPI-like capabilities together with RPC functionality. The GraphLab distributed 
implementation is built on top of this RPC library.

GraphLab RPC uses extensive template meta-programming techniques to provide
an \b IDL-free (http://en.wikipedia.org/wiki/Interface_description_language) 
RPC system, allowing arbitrary functions to be called on program running on 
remote machines (Note that all machines must the running the same binary).

For instance, this is a particularly interesting example:
\code
#include <iostream>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_init_from_mpi.hpp>
using namespace graphlab;

int main(int argc, char ** argv) {
  mpi_tools::init(argc, argv);

  dc_init_param param;
  if (init_param_from_mpi(param) == false) {
    return 0;
  }
  distributed_control dc(param);
  
  if (dc.procid() == 0 && dc.numprocs() >= 2) {
    dc.remote_call(1, printf, "%d + %f = %s\n", 1, 2.0, "three");
  }
  dc.barrier();
}
\endcode

GraphLab RPC uses MPI for initialization. <tt> mpi_tools::init(argc, argv) </tt>
initializes the MPI library, and <tt> init_param_from_mpi </tt> performs the 
initial RPC negotiation.

\remarks
MPI provides a convenient platform for spawning distributed processes,
and also allows GraphLab RPC to perform negotiation of initial port numbers.
An alternate spawning is available (\ref init_param_from_env ), but it is less
convenient as the default ports may be used, or blocked.


Once the distributed_control object is created, \ref distributed_control::procid "dc.procid()"
provides the current machine number, while \ref distributed_control::numprocs "dc.numprocs()"
provide the total number of machines.

The if-condition is therefore entered by only the first machine, which
performs a remote call to the second machine (the first argument of remote_call
is the target machine ID). The second machine will then
execute the equivalent of
\code
  printf("%d + %f = %s\n", 1, 2.0, "three");
\endcode



*/