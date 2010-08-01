#include <boost/foreach.hpp>
#include <stdint.h>

// if GNUC is available, this checks if the file which included
// macros_def.hpp is the same file which included macros_undef.hpp
#ifdef __GNUC__
#define GRAPHLAB_MACROS_INC_LEVEL __INCLUDE_LEVEL__
#endif


// prevent this file from being included before other graphlab headers
#ifdef GRAPHLAB_MACROS
#error "Repeated include of <macros_def.hpp>. This probably means that macros_def.hpp was not the last include, or some header file failed to include <macros_undef.hpp>"
#endif

#define GRAPHLAB_MACROS

/** A macro to disallow the copy constructor and operator= functions
    This should be used in the private: declarations for a class */
#define DISALLOW_COPY_AND_ASSIGN(TypeName) \
TypeName(const TypeName&);               \
void operator=(const TypeName&);



// Shortcut macro definitions
//! see http://www.boost.org/doc/html/foreach.html 
#define foreach BOOST_FOREACH

#define rev_foreach BOOST_REVERSE_FOREACH




