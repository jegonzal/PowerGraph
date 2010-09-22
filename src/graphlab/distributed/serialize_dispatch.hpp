#ifndef SERIALIZE_DISPATCH_HPP
#define SERIALIZE_DISPATCH_HPP
#include <string>
#include <boost/type_traits/decay.hpp>
#include <boost/type_traits/remove_pointer.hpp>
#include <boost/type_traits/function_traits.hpp>
#include <serialization/serialization_includes.hpp>
namespace graphlab {
class distributed_control;

#define F0 typename boost::remove_reference<typename boost::function_traits<typename boost::remove_pointer<F>::type>::arg5_type>::type
#define F1 typename boost::remove_reference<typename boost::function_traits<typename boost::remove_pointer<F>::type>::arg6_type>::type
#define F2 typename boost::remove_reference<typename boost::function_traits<typename boost::remove_pointer<F>::type>::arg7_type>::type
#define F3 typename boost::remove_reference<typename boost::function_traits<typename boost::remove_pointer<F>::type>::arg8_type>::type
#define F4 typename boost::remove_reference<typename boost::function_traits<typename boost::remove_pointer<F>::type>::arg9_type>::type
#define F5 typename boost::remove_reference<typename boost::function_traits<typename boost::remove_pointer<F>::type>::arg10_type>::type

/*
This set of functions REPACKSTRUCT(F, ...) and DISPATCHn(...)
Allows for a set of arguments to be repacked and transmitted to remote machines.
The basic idea is to generate a struct which stores the arguments, and
transmit a byte for byte copy of the struct. 

We also transmit a pointer to DISPATCH<structtype> which is a generated template 
function which knows how to understand the struct.
*/


template <typename F>
struct SERPACKSTRUCT0{ 
  F targetfn;
  // for efficiency reasons, and to avoid an additional copy, the save
  // is not actually done here, but is done by the REPACKSTRUCT functions
  void save(oarchive &arc) const {
    arc << reinterpret_cast<size_t>(targetfn);
  }
  void load(iarchive &arc) {
    size_t i;
    arc >> i;
    targetfn = reinterpret_cast<F>(i);
  }
};

template <typename F, typename T0>
struct SERPACKSTRUCT1{ 
  F targetfn; T0 i0; 
  void save(oarchive &arc) const {
    arc << reinterpret_cast<size_t>(targetfn) << i0;
  }
  void load(iarchive &arc) {
    size_t i;
    arc >> i >> i0;
    targetfn = reinterpret_cast<F>(i);
  }
};

template <typename F, typename T0, typename T1>
struct SERPACKSTRUCT2{ 
  F targetfn; T0 i0; T1 i1; 
  void save(oarchive &arc) const {
    arc << reinterpret_cast<size_t>(targetfn) << i0 << i1;
  }
  void load(iarchive &arc) {
    size_t i;
    arc >> i >> i0 >> i1;
    targetfn = reinterpret_cast<F>(i);
  }
};

template <typename F, typename T0, typename T1, typename T2>
struct SERPACKSTRUCT3{ 
  F targetfn; T0 i0; T1 i1; T2 i2; 
  void save(oarchive &arc) const {
    arc << reinterpret_cast<size_t>(targetfn) << i0 << i1 << i2;
  }
  void load(iarchive &arc) {
    size_t i;
    arc >> i >> i0 >> i1 >> i2;
    targetfn = reinterpret_cast<F>(i);
  }
};

template<typename F, typename T0, typename T1, typename T2, typename T3>
struct SERPACKSTRUCT4{ 
  F targetfn; T0 i0; T1 i1; T2 i2; T3 i3; 
  void save(oarchive &arc) const {
    arc << reinterpret_cast<size_t>(targetfn) << i0 << i1 << i2 << i3;
  }
  void load(iarchive &arc) {
    size_t i;
    arc >> i >> i0 >> i1 >> i2 >> i3;
    targetfn = reinterpret_cast<F>(i);
  }
};

template<typename F, typename T0, typename T1, typename T2, typename T3, 
                    typename T4>
struct SERPACKSTRUCT5{ 
  F targetfn; T0 i0; T1 i1; T2 i2; T3 i3; T4 i4; 
  void save(oarchive &arc) const {
    arc << reinterpret_cast<size_t>(targetfn) << i0 << i1 << i2 << i3 << i4;
  }
  void load(iarchive &arc) {
    size_t i;
    arc >> i >> i0 >> i1 >> i2 >> i3 >> i4;
    targetfn = reinterpret_cast<F>(i);
  }
};

template<typename F, typename T0, typename T1, typename T2, typename T3, 
                    typename T4, typename T5>
struct SERPACKSTRUCT6{ 
  F targetfn; T0 i0; T1 i1; T2 i2; T3 i3; T4 i4; T5 i5; 
  void save(oarchive &arc) const {
    arc << reinterpret_cast<size_t>(targetfn) << i0 << i1 << i2 << i3 << i4 << i5;
  }
  void load(iarchive &arc) {
    size_t i;
    arc >> i >> i0 >> i1 >> i2 >> i3 >> i4 >> i5;
    targetfn = reinterpret_cast<F>(i);
  }
};

template <typename F>
class SER_GET_STRUCT_TYPE0 { public: typedef SERPACKSTRUCT0<F> struct_type; };

template <typename F>
class SER_GET_STRUCT_TYPE1 { public: typedef SERPACKSTRUCT1<F, F0> struct_type; };

template <typename F>
class SER_GET_STRUCT_TYPE2 { public: typedef SERPACKSTRUCT2<F, F0, F1> struct_type; };

template <typename F>
class SER_GET_STRUCT_TYPE3 { public: typedef SERPACKSTRUCT3<F, F0, F1, F2> struct_type; };

template <typename F>
class SER_GET_STRUCT_TYPE4 { public: typedef SERPACKSTRUCT4<F, F0, F1, F2, F3> struct_type; };

template <typename F>
class SER_GET_STRUCT_TYPE5 { public: typedef SERPACKSTRUCT5<F, F0, F1, F2, F3, F4> struct_type; };

template <typename F>
class SER_GET_STRUCT_TYPE6 { public: typedef SERPACKSTRUCT6<F, F0, F1, F2, F3, F4, F5> struct_type; };


template<typename F>
void SERREPACKSTRUCT(F fn, oarchive &arc) {
  arc << reinterpret_cast<size_t>(fn);
}

template<typename F, typename T0>
void SERREPACKSTRUCT(F fn, oarchive &arc, T0 i0) {
  arc << reinterpret_cast<size_t>(fn) << (F0)(i0);
}

template<typename F, typename T0, typename T1>
void SERREPACKSTRUCT(F fn, oarchive &arc, T0 i0, T1 i1) {
  arc << reinterpret_cast<size_t>(fn) << (F0)(i0) << (F1)(i1);
}

template<typename F, typename T0, typename T1, typename T2>
inline void SERREPACKSTRUCT(F fn, oarchive &arc, T0 i0, T1 i1, T2 i2) {
  arc << reinterpret_cast<size_t>(fn) << (F0)(i0) << (F1)(i1) << (F2)(i2);
}

template<typename F, typename T0, typename T1, typename T2, typename T3>
inline void SERREPACKSTRUCT(F fn, oarchive &arc, T0 i0, T1 i1, T2 i2, T3 i3) {
  arc << reinterpret_cast<size_t>(fn) << (F0)(i0) << (F1)(i1) << (F2)(i2) << (F3)(i3);
}

template<typename F, typename T0, typename T1, typename T2, typename T3,
         typename T4>
inline void SERREPACKSTRUCT(F fn, oarchive &arc, T0 i0, T1 i1, T2 i2, T3 i3, T4 i4) {
  arc << reinterpret_cast<size_t>(fn) << (F0)(i0) << (F1)(i1) << (F2)(i2) << (F3)(i3) << F4(i4);
}

template<typename F, typename T0, typename T1, typename T2, typename T3,
         typename T4, typename T5>
inline void SERREPACKSTRUCT(F fn, oarchive &arc, T0 i0, T1 i1, T2 i2, T3 i3, T4 i4, T5 i5) {
  arc << reinterpret_cast<size_t>(fn) << (F0)(i0) << (F1)(i1) << (F2)(i2) << (F3)(i3) << F4(i4) << F5(i5);
}



template<typename PSTRUCT>
void SERDISPATCH0(distributed_control& dc, size_t source, 
         void* ptr, size_t len, iarchive &arc) {
  PSTRUCT ps;
  ps.load(arc);
  (ps.targetfn)(dc, source, ptr, len);
}

template<typename PSTRUCT>
void SERDISPATCH1(distributed_control& dc, size_t source, 
         void* ptr, size_t len, iarchive &arc) {
  PSTRUCT ps;
  ps.load(arc);
  (ps.targetfn)(dc, source, ptr, len, ps.i0);
}

template<typename PSTRUCT>
void SERDISPATCH2(distributed_control& dc, size_t source, 
         void* ptr, size_t len, iarchive &arc) {
  PSTRUCT ps;
  ps.load(arc);
  (ps.targetfn)(dc, source, ptr, len, ps.i0, ps.i1);
}

template<typename PSTRUCT>
void SERDISPATCH3(distributed_control& dc, size_t source, 
         void* ptr, size_t len, iarchive &arc) {
  PSTRUCT ps;
  ps.load(arc);
  (ps.targetfn)(dc, source, ptr, len, ps.i0, ps.i1, ps.i2);
}

template<typename PSTRUCT>
void SERDISPATCH4(distributed_control& dc, size_t source, 
         void* ptr, size_t len, iarchive &arc) {
  PSTRUCT ps;
  ps.load(arc);
  (ps.targetfn)(dc, source, ptr, len, ps.i0, ps.i1, ps.i2, ps.i3);
}

template<typename PSTRUCT>
void SERDISPATCH5(distributed_control& dc, size_t source, 
         void* ptr, size_t len, iarchive &arc) {
  PSTRUCT ps;
  ps.load(arc);
  (ps.targetfn)(dc, source, ptr, len, ps.i0, ps.i1, ps.i2, ps.i3, ps.i4);
}

template<typename PSTRUCT>
void SERDISPATCH6(distributed_control& dc, size_t source, 
         void* ptr, size_t len, iarchive &arc) {
  PSTRUCT ps;
  ps.load(arc);
  (ps.targetfn)(dc, source, ptr, len, ps.i0, ps.i1, ps.i2, ps.i3, ps.i4, 
                 ps.i5);
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

