#ifndef FAKE_SHARED_DATA_HPP
#define FAKE_SHARED_DATA_HPP
#include <graphlab/util/generics/any.hpp>
namespace graphlab {
/** \internal
 * Temporary backward compatibility hack for the shared data manager
 * Do not use!
 */
template<typename Engine>
class fake_shared_data:public ishared_data<typename Engine::graph_type> {
 private:
  Engine* engine;
  any trash;
 public:
  fake_shared_data(Engine* engine):engine(engine) {}
  const any& get_constant(size_t index) const { ASSERT_TRUE(false); return trash;}
  any get(size_t index) const { ASSERT_TRUE(false); return any(); }
  any atomic_get(size_t index) const { ASSERT_TRUE(false); return any(); }
  void atomic_set(size_t index, const any& data) { ASSERT_TRUE(false); }
  any atomic_exchange(size_t index, const any& data) { ASSERT_TRUE(false); return any(); }
  void atomic_apply(size_t index,
                    typename ishared_data<typename Engine::graph_type>::apply_function_type fun,     
                    const any& closure) { ASSERT_TRUE(false); }

  void trigger_sync(size_t index) { ASSERT_TRUE(false); }
  void trigger_sync_all() {
    engine->sync_all_soon();
  }
  void trigger_sync(glshared_base &var) {
    engine->sync_soon(var);
  }
};
}
#endif