#ifdef F0
#error "multiple includes of function arg types" 
#endif

#include <boost/type_traits/decay.hpp>
#include <boost/type_traits/remove_pointer.hpp>
#include <boost/type_traits/remove_const.hpp>
#include <boost/type_traits/function_traits.hpp>

#include <boost/function.hpp>
/*
A huge collection of useful typedefs.
F0... F5: identifies the arguments for an RPC aware function F. Dropping const and dropping references 
          (therefore allowing you to use F0....F5 to do casting.

NIF0... NIF5: identifies the arguments for an RPC unaware function F

R0.... R7: Identifies the actual arguments of the function F, without de-consting and de-reffing

FRESULT: de-const and de-refed type of F's return type

FARITY: the number of arguments F takes
*/
#define REMOVE_CONST_REF(REF) typename boost::remove_const<typename boost::remove_reference<REF>::type>::type
// non-intrusive calls
#define NIF0 REMOVE_CONST_REF(typename boost::function<typename boost::remove_pointer<F>::type>::arg1_type)
#define NIF1 REMOVE_CONST_REF(typename boost::function<typename boost::remove_pointer<F>::type>::arg2_type)
#define NIF2 REMOVE_CONST_REF(typename boost::function<typename boost::remove_pointer<F>::type>::arg3_type)
#define NIF3 REMOVE_CONST_REF(typename boost::function<typename boost::remove_pointer<F>::type>::arg4_type)
#define NIF4 REMOVE_CONST_REF(typename boost::function<typename boost::remove_pointer<F>::type>::arg5_type)
#define NIF5 REMOVE_CONST_REF(typename boost::function<typename boost::remove_pointer<F>::type>::arg6_type)
#define NIF6 REMOVE_CONST_REF(typename boost::function<typename boost::remove_pointer<F>::type>::arg7_type)
#define NIF7 REMOVE_CONST_REF(typename boost::function<typename boost::remove_pointer<F>::type>::arg8_type)

// for rpc aware functions
#define F0 NIF2
#define F1 NIF3
#define F2 NIF4
#define F3 NIF5
#define F4 NIF6
#define F5 NIF7



#define R0 typename boost::function<typename boost::remove_pointer<F>::type>::arg1_type
#define R1 typename boost::function<typename boost::remove_pointer<F>::type>::arg2_type
#define R2 typename boost::function<typename boost::remove_pointer<F>::type>::arg3_type
#define R3 typename boost::function<typename boost::remove_pointer<F>::type>::arg4_type
#define R4 typename boost::function<typename boost::remove_pointer<F>::type>::arg5_type
#define R5 typename boost::function<typename boost::remove_pointer<F>::type>::arg6_type
#define R6 typename boost::function<typename boost::remove_pointer<F>::type>::arg7_type
#define R7 typename boost::function<typename boost::remove_pointer<F>::type>::arg8_type

#define FRESULT REMOVE_CONST_REF(typename boost::function<typename boost::remove_pointer<F>::type>::result_type)

#define FARITY boost::function<typename boost::remove_pointer<F>::type>::arity

