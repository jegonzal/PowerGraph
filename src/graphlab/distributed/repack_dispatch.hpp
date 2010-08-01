#ifndef REPACK_DISPATCH_HPP
#define REPACK_DISPATCH_HPP
#include <boost/type_traits/decay.hpp>
#include <boost/type_traits/remove_pointer.hpp>
#include <boost/type_traits/function_traits.hpp>
namespace graphlab {
class distributed_control;

#define F0 typename boost::function_traits<typename boost::remove_pointer<F>::type>::arg5_type
#define F1 typename boost::function_traits<typename boost::remove_pointer<F>::type>::arg6_type
#define F2 typename boost::function_traits<typename boost::remove_pointer<F>::type>::arg7_type
#define F3 typename boost::function_traits<typename boost::remove_pointer<F>::type>::arg8_type
#define F4 typename boost::function_traits<typename boost::remove_pointer<F>::type>::arg9_type
#define F5 typename boost::function_traits<typename boost::remove_pointer<F>::type>::arg10_type

/*
This set of functions REPACKSTRUCT(F, ...) and DISPATCHn(...)
Allows for a set of arguments to be repacked and transmitted to remote machines.
The basic idea is to generate a struct which stores the arguments, and
transmit a byte for byte copy of the struct. 

We also transmit a pointer to DISPATCH<structtype> which is a generated template 
function which knows how to understand the struct.
*/
template <typename F>
struct PACKSTRUCT0{ F targetfn;};

template <typename F, typename T0>
struct PACKSTRUCT1{ F targetfn; T0 i0; };

template <typename F, typename T0, typename T1>
struct PACKSTRUCT2{ F targetfn; T0 i0; T1 i1; };

template <typename F, typename T0, typename T1, typename T2>
struct PACKSTRUCT3{ F targetfn; T0 i0; T1 i1; T2 i2; };

template<typename F, typename T0, typename T1, typename T2, typename T3>
struct PACKSTRUCT4{ F targetfn; T0 i0; T1 i1; T2 i2; T3 i3; };

template<typename F, typename T0, typename T1, typename T2, typename T3, 
                    typename T4>
struct PACKSTRUCT5{ F targetfn; T0 i0; T1 i1; T2 i2; T3 i3; T4 i4; };

template<typename F, typename T0, typename T1, typename T2, typename T3, 
                    typename T4, typename T5>
struct PACKSTRUCT6{ F targetfn; T0 i0; T1 i1; T2 i2; T3 i3; T4 i4; T5 i5; };


template<typename F>
inline PACKSTRUCT0<F> REPACKSTRUCT(F fn) {
  PACKSTRUCT0<F> p;
  p.targetfn = fn;
  return p;
}

template<typename F, typename T0>
inline PACKSTRUCT1<F, F0> REPACKSTRUCT(F fn, T0 i0) {
  PACKSTRUCT1<F, F0> p;
  p.targetfn = fn;
  p.i0 = i0;
  return p;
}

template<typename F, typename T0, typename T1>
PACKSTRUCT2<F, F0, F1> REPACKSTRUCT(F fn, T0 i0, T1 i1) {
  PACKSTRUCT2<F, F0, F1> p;
  p.targetfn = fn;
  p.i0 = i0; p.i1 = i1;
  return p;
}

template<typename F, typename T0, typename T1, typename T2>
inline PACKSTRUCT3<F, F0, F1, F2> REPACKSTRUCT(F fn, T0 i0, T1 i1, T2 i2) {
  PACKSTRUCT3<F, F0, F1, F2> p;
  p.targetfn = fn;
  p.i0 = i0; p.i1 = i1; p.i2 = i2;
  return p;
}
template<typename F, typename T0, typename T1, typename T2, typename T3>
inline PACKSTRUCT4<F, F0, F1, F2, F3> REPACKSTRUCT(F fn, T0 i0, T1 i1, T2 i2, T3 i3) {
  PACKSTRUCT4<F, F0, F1, F2, F3> p;
  p.targetfn = fn;
  p.i0 = i0; p.i1 = i1; p.i2 = i2;
  p.i3 = i3;
  return p;
}
template<typename F, typename T0, typename T1, typename T2, typename T3,
         typename T4>
inline PACKSTRUCT5<F, F0, F1, F2, F3, F4> 
       REPACKSTRUCT(F fn, T0 i0, T1 i1, T2 i2, T3 i3, T4 i4) {
  PACKSTRUCT5<F, F0, F1, F2, F3, F4> p;
  p.targetfn = fn;
  p.i0 = i0; p.i1 = i1; p.i2 = i2;
  p.i3 = i3; p.i4 = i4;
  return p;
}

template<typename F, typename T0, typename T1, typename T2, typename T3,
         typename T4, typename T5>
inline PACKSTRUCT6<F, F0, F1, F2, F3, F4, F5> 
       REPACKSTRUCT(F fn, T0 i0, T1 i1, T2 i2, T3 i3, T4 i4, T5 i5) {
  PACKSTRUCT6<F, F0, F1, F2, F3, F4, F5> p;
  p.targetfn = fn;
  p.i0 = i0; p.i1 = i1; p.i2 = i2;
  p.i3 = i3; p.i4 = i4; p.i5 = i5;
  return p;
}



template<typename PSTRUCT>
void DISPATCH0(distributed_control& dc, size_t source, 
         void* ptr, size_t len, void* stack) {
  PSTRUCT* ps = (PSTRUCT*)stack;
  (ps->targetfn)(dc, source, ptr, len);
}

template<typename PSTRUCT>
void DISPATCH1(distributed_control& dc, size_t source, 
         void* ptr, size_t len, void* stack) {
  PSTRUCT* ps = (PSTRUCT*)stack;
  (ps->targetfn)(dc, source, ptr, len, ps->i0);
}

template<typename PSTRUCT>
void DISPATCH2(distributed_control& dc, size_t source, 
         void* ptr, size_t len, void* stack) {
  PSTRUCT* ps = (PSTRUCT*)stack;
  (ps->targetfn)(dc, source, ptr, len, ps->i0, ps->i1);
}

template<typename PSTRUCT>
void DISPATCH3(distributed_control& dc, size_t source, 
         void* ptr, size_t len, void* stack) {
  PSTRUCT* ps = (PSTRUCT*)stack;
  (ps->targetfn)(dc, source, ptr, len, ps->i0, ps->i1, ps->i2);
}

template<typename PSTRUCT>
void DISPATCH4(distributed_control& dc, size_t source, 
         void* ptr, size_t len, void *stack) {
  PSTRUCT* ps = (PSTRUCT*)stack;
  (ps->targetfn)(dc, source, ptr, len, ps->i0, ps->i1, ps->i2, ps->i3);
}

template<typename PSTRUCT>
void DISPATCH5(distributed_control& dc, size_t source, 
         void* ptr, size_t len, void* stack) {
  PSTRUCT* ps = (PSTRUCT*)stack;
  (ps->targetfn)(dc, source, ptr, len, ps->i0, ps->i1, ps->i2, ps->i3, ps->i4);
}

template<typename PSTRUCT>
void DISPATCH6(distributed_control& dc, size_t source, 
         void* ptr, size_t len, void* stack) {
  PSTRUCT* ps = (PSTRUCT*)stack;
  (ps->targetfn)(dc, source, ptr, len, ps->i0, ps->i1, ps->i2, ps->i3, ps->i4, 
                 ps->i5);
}


}
#undef F0
#undef F1
#undef F2
#undef F3
#undef F4
#undef F5
#undef F6
#undef F7
#endif

