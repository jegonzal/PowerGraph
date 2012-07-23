#ifndef CHANDY_MISRA_INTERFACE_HPP
#define CHANDY_MISRA_INTERFACE_HPP

template <typename GraphType>
class chandy_misra_interface {
 public:
  typedef typename GraphType::lvid_type lvid_type;

  virtual size_t num_clean_forks() const = 0;
  virtual void make_philosopher_hungry_per_replica(lvid_type p_id) = 0; 
  virtual void philosopher_stops_eating_per_replica(lvid_type p_id) = 0;
 
};
#endif
