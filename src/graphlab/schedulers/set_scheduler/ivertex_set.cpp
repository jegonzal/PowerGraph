
// #include <graphlab/schedulers/set_scheduler/ivertex_set.hpp>
// #include <graphlab/schedulers/set_scheduler/set_generic_typedefs.hpp>

// namespace graphlab {

//   void ivertex_set::add_event_target(ivertex_set *ev) {
//     if (ev != NULL) {
//       all_children.push_back(ev);
//     }
//   }


//   size_t ivertex_set::all_dependent_events() {
//     unsigned char dependentevs = supported_events();
//     for (size_t i = 0; i < all_children.size(); ++i) {
//       dependentevs |= all_children[i]->all_dependent_events();
//     }
//     return dependentevs;
//   }

//   void ivertex_set::resolve_event_handlers() {
//     for (size_t i = 0; i < all_children.size(); ++i) {
//       ivertex_set *ev = all_children[i];
//       size_t supportevs = ev->all_dependent_events();
//       if (supportevs & INSERT_EVENT) insert_events.push_back(ev);
//       if (supportevs & ERASE_EVENT) erase_events.push_back(ev);
//       if (supportevs & MODIFY_VERTEX_EVENT) modify_vertex_events.push_back(ev);
//       if (supportevs & MODIFY_EDGE_EVENT) modify_edge_events.push_back(ev);
//       ev->resolve_event_handlers();
//     }
//   }

// }
