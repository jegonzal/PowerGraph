#ifdef __GNUC__
#if (GRAPHLAB_MACROS_INC_LEVEL != __INCLUDE_LEVEL__)
  #error "A <macros_def.hpp> was not paired with a <macros_undef.hpp>"  
#endif
#undef GRAPHLAB_MACROS_INC_LEVEL
#endif


#undef GRAPHLAB_MACROS
#undef DISALLOW_COPY_AND_ASSIGN
#undef foreach
#undef rev_foreach
