#ifndef GRAPHLAB_UTIL_BRANCH_HINTS_HPP
#define GRAPHLAB_UTIL_BRANCH_HINTS_HPP

#define __likely__(x)       __builtin_expect((x),1)
#define __unlikely__(x)     __builtin_expect((x),0)

#endif //GRAPHLAB_UTIL_BRANCH_HINTS_HPP

