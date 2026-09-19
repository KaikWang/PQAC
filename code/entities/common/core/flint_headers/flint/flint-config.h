/* Define if the compiler supports and should use thread-local storage */
#define FLINT_USES_TLS 1

/* Define if the library should be thread-safe, no matter whether FLINT_USES_TLS is used */
#define FLINT_REENTRANT 1

/* Define if -DCMAKE_BUILD_TYPE=Debug was given, to enable some ASSERT()s */
/* #undef FLINT_WANT_ASSERT */

/* Define if you cpu_set_t in sched.h */
#define FLINT_USES_CPUSET 0

#define FLINT_USES_PTHREAD 1

#define FLINT_USES_POPCNT

#define FLINT_USES_BLAS 0

#define FLINT_USES_FENV 1

#ifdef _MSC_VER
#define access _access
#define strcasecmp _stricmp
#define strncasecmp	_strnicmp
#define alloca _alloca
#define MSC_C_(x) #x  
#define MSC_CC_(x)  MSC_C_(x)
#define MSC_VERSION "Microsoft C++ (Version " MSC_CC_(_MSC_FULL_VER) ")"
#endif

#if defined (FLINT_BUILD_DLL)
#define FLINT_DLL __declspec(dllexport)
#elif defined(MSC_USE_DLL)
#define FLINT_DLL __declspec(dllimport)
#else
#define FLINT_DLL
#endif
