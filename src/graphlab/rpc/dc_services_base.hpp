#ifndef DC_SERVICES_BASE_HPP
#define DC_SERVICES_BASE_HPP
/**
A base class describing an implemention of 
distributed communication services. These include things such
as "barrier" and "reduce" or "ping".
This allows different interfaces to provide different 
implementions of such services. (for instance, if MPI is available,
we would just use MPI_Barrier()
*/
class dc_services_base {
 public:
  /** Performs a distributed barrier. When the process
  hits a barrier, it will pause until all machines hits
  the barrier. Do note that barriers are not named. 
  Therefore it is possible for different machines to hit the
  barrier in different lines of code.*/
  virtual void barrier() = 0;
  
  /**
  This function allows one machine to broadcasts a string of data 
  to all machines.
  
  The originator calls broadcast with data provided in 
  in 'data' and length len. All other machines must call
  broadcast with data = NULL. 
  
  The originator will then return 'data'. All other machines
  will return a new pointer, with the length of the string
  returned in 'len'. The returned pointer must be freed
  by the caller (with the exception of the originator). 
  
  This function is not guaranteed to have barrier-like behavior.
  That is, broadcast could be implemented in a buffered fashion.
  */
  virtual char* broadcast(char* data, size_t &len) = 0;
  virtual ~dc_services_base() { }
};
#endif
