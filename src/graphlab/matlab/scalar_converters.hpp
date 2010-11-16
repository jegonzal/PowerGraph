#ifndef SCALAR_GENERATORS_HPP
#define SCALAR_GENERATORS_HPP
#include <boost/preprocessor.hpp>

#define TYPES (char)(int8_T)(uint8_T)(int16_T)(uint16_T)(int32_T)\
              (uint32_T)(int64_T)(uint64_T)(real32_T)(real64_T)
#define COMPLEX_TYPES (cint8_T)(cuint8_T)(cint16_T)(cuint16_T)(cint32_T)\
                      (cuint32_T)(cint64_T)(cuint64_T)(creal_T)(creal32_T)(creal64_T)


#define GEN_CONVERTERS(r, _, TYPENAME) \
template <>                                                        \
struct converter<TYPENAME>{                                                  \
static bool mxarray2emx(const mxArray* mx, TYPENAME &emxdata) {           \
  emxdata = mxGetScalar(mx);                                       \
}                                                                  \
static bool emx2mxarray(TYPENAME &emxdata, mxArray* &mx) {                \
  mx = mxCreateDoubleScalar(emxdata);                              \
}                                                                  \
static void clearemx(TYPENAME &emxdata) {                                 \
  emxdata = (TYPENAME)0;                                           \
}                                                                  \
static void freeemx(TYPENAME &emxdata) {                                  \
}                                                                  \
static void emxcopy(TYPENAME &dest, const TYPENAME &src) {        \
  dest = src;                                                     \
}                                                                  \
};
BOOST_PP_SEQ_FOR_EACH(GEN_CONVERTERS, _, TYPES);

#undef GEN_CONVERTERS
#undef TYPES
#undef COMPLEX_TYPES

#endif