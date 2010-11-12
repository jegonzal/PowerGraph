#include <mex.h>
#include "b2/rtwtypes.h"
#include "b2/test_types.h"
#include "b2/generator.hpp"
#include "b2/test_initialize.h"
#include "b2/test.h"



#include <boost/type_traits/decay.hpp>
#include <boost/type_traits/remove_pointer.hpp>
#include <boost/type_traits/remove_const.hpp>
#include <boost/type_traits/function_traits.hpp>

#include <boost/function.hpp>

#define REMOVE_CONST_REF(REF) typename boost::remove_const<typename boost::remove_reference<REF>::type>::type
#define NIF0 REMOVE_CONST_REF(typename boost::function<typename boost::remove_pointer<F>::type>::arg1_type)
#define NIF1 REMOVE_CONST_REF(typename boost::function<typename boost::remove_pointer<F>::type>::arg2_type)
#define FRESULT REMOVE_CONST_REF(typename boost::function<typename boost::remove_pointer<F>::type>::result_type)



template <typename F>
void exec2arg(F f, const mxArray *rhs, mxArray *&lhs) {
  REMOVE_CONST_REF(typename boost::remove_pointer<NIF0>::type) b;
  //REMOVE_CONST_REF(typename boost::remove_pointer<NIF1>::type) a;
  REMOVE_CONST_REF(typename boost::remove_pointer<FRESULT>::type) a;
  clearemx(a);
  mxarray2emx(rhs, b);
  
  a = f(b);

  emx2mxarray(a, lhs);
  freeemx(a);
  freeemx(b);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  if (nlhs != 1 || nrhs != 1) {
    mexPrintf("args invalid\n");
    return;
  }
  test_initialize();
  exec2arg(test, prhs[0], plhs[0]);
}
