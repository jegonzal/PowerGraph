/* 
 * File:   iaggregate.hpp
 * Author: jegonzal
 *
 * Created on July 14, 2011, 5:56 PM
 */

#ifndef GRAPHLAB_IAGGREGATE_HPP
#define	GRAPHLAB_IAGGREGATE_HPP



namespace graphlab {

  /**
   * \brief  This class is the base type of an aggregation (sync)
   * operation. 
   *
   * Each thread is assumed to have a a single iaggregate which it
   * uses to assemble the the partial sum. 
   */
  template<typename Graph> 
  class iaggregate {
  public:
    typedef Graph graph_type;
    typedef iengine<graph_type> iengine_type;
    typedef typename iengine_type::iscope_type;
    virtual ~iaggregate();
    virtual void add(iscope_type& scope) = 0;
    virtual void add(const iaggregate* other_ptr) = 0;
    virtual bool add_delta(const void* accumulator_ptr) = 0;
    virtual void apply(glshared_base& base) = 0;
    virtual void clear() = 0;
  };






  template<typename Graph, typename Accum>
  class aggregate_task : public iaggregate<Graph> {
  public:
    typedef Accum accumulator_type;
    typedef Accum(*map_function_type)(iscope_type& scope);
    typedef Accum(*reduce_function_type)(const Accum& lvalue, 
                                         const Accum& rvalue);
    typedef void(*apply_function_type)(const Accum& accum, glshared_base& target);

  private:
    accumulator_type partial_sum;
    bool cleared;
    map_function_type map_function;
    reduce_function_type reduce_function;
  public:
    void add(iscope_type& scope) {
      if(cleared) partial_sum = map_function(scope);
      else partial_sum = reduce_function(partial_sum, map_function(scope));
    }
    
    void add(const iaggregate* other_ptr) {
      ASSERT_TRUE(other_ptr != NULL);
      const aggregate_task* agg_task = *((const aggregate_task*)other_ptr);
      if(!agg_task->cleared) 
        partial_sum = reduce_function(partial_sum, agg_task->partial_sum); 
    }
    
    
    

  };
  
  
}; // end of namespace graphlab


#endif	/* IAGGREGATE_HPP */

