#ifndef SCALAR_GENERATORS_HPP
#define SCALAR_GENERATORS_HPP
#include <boost/preprocessor.hpp>

#define TYPES (char)(int8_T)(uint8_T)(int16_T)(uint16_T)(int32_T)\
              (uint32_T)(int64_T)(uint64_T)(real32_T)(real64_T)
#define COMPLEX_TYPES (cint8_T)(cuint8_T)(cint16_T)(cuint16_T)(cint32_T)\
                      (cuint32_T)(cint64_T)(cuint64_T)(creal_T)(creal32_T)(creal64_T)


#define GEN_CONVERTERS(r, _, TYPENAME) \
template <>                                                        \
bool mxarray2emx<TYPENAME>(const mxArray* mx, TYPENAME &emxdata) { \
  emxdata = mxGetScalar(mx);                                       \
}                                                                  \
template <>                                                        \
bool emx2mxarray<TYPENAME>(TYPENAME &emxdata, mxArray* &mx) {      \
  mx = mxCreateDoubleScalar(emxdata);                              \
}                                                                  \
template <>                                                        \
void clearemx<TYPENAME>(TYPENAME &emxdata) {                       \
  emxdata = (TYPENAME)0;                                           \
}                               

BOOST_PP_SEQ_FOR_EACH(GEN_CONVERTERS, _, TYPES);

#undef GEN_CONVERTERS
#undef TYPES
#undef COMPLEX_TYPES

#endif