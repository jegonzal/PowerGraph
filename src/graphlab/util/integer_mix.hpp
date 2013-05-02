#ifndef GRAPHLAB_INTEGER_MIX_HPP
#define GRAPHLAB_INTEGER_MIX_HPP
#include <stdint.h>
namespace graphlab {
// Jenkin's 32 bit integer mix from
// http://burtleburtle.net/bob/hash/integer.html
uint32_t integer_mix(uint32_t a) {
  a -= (a<<6);
  a ^= (a>>17);
  a -= (a<<9);
  a ^= (a<<4);
  a -= (a<<3);
  a ^= (a<<10);
  a ^= (a>>15);
  return a;
}

}
#endif

