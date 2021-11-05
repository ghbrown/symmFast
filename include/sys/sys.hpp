#ifndef SYMMFAST_SYS_SYSHPP
#define SYMMFAST_SYS_SYSHPP

#ifdef __GNUC__
#  define SF_FUNCTION_NAME __PRETTY_FUNCTION__
#else
#  define SF_FUNCTION_NAME __func__
#endif

#if __cpluslpus >= 201703L // c++17
#  define SF_NODISCARD [[nodiscard]]
#else
#  define SF_NODISCARD
#endif

#define SF_ATTRIBUTE(attr) __attribute__((attr))

#define SF_UNUSED SF_ATTRIBUTE(unused)

#define SF_VISIBILITY_EXTERNAL SF_ATTRIBUTE(visibility ("default"))
#define SF_VISIBILITY_INTERNAL SF_ATTRIBUTE(visibility ("hidden"))

#define	SF_INTERN extern SF_NODISCARD SF_VISIBILITY_INTERNAL
#define SF_EXTERN extern SF_NODISCARD SF_VISIBILITY_EXTERNAL

#define SF_INLINE        inline
#define	SF_STATIC_INLINE static SF_INLINE

#define sfunlikely(pred) __builtin_expect(!!(pred),0)
#define sflikely(pred)   __builtin_expect(!!(pred),1)

// todo implement the detail:: parts
#ifdef SF_DEBUG
#  define SFCHECK(expr) do {						\
    const auto errc = expr;						\
    if (sfunlikely(detail::sf_check_error(errc))) {			\
      return detail::error_handler(__FILE__,SF_FUNCTION_NAME,		\
				   __LINE__,SF_COMM_WORLD,errc);	\
    }									\
  } while (0)
#else
#  define SFCHECK(expr)	do { const auto SF_UNUSED errc = expr; } while (0)
#endif

#include <mpi.h>

namespace sf
{

using sf_error_t = int;

SF_EXTERN MPI_Comm SF_COMM_WORLD;
SF_EXTERN MPI_Comm SF_COMM_SELF;

SF_EXTERN sf_error_t initialize(int,char*[]);
SF_EXTERN sf_error_t finalize();

} // namespace sf

#endif // SYMMFAST_SYS_SYSHPP
