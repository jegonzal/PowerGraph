/*
This file is part of GraphLab.

GraphLab is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as 
published by the Free Software Foundation, either version 3 of 
the License, or (at your option) any later version.

GraphLab is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public 
License along with GraphLab.  If not, see <http://www.gnu.org/licenses/>.
*/

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




