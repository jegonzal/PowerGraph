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

