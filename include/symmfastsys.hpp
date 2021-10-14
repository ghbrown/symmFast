#ifndef SYMMFASTSYS_HPP
#define SYMMFASTSYS_HPP

#ifdef __GNUC__
#  define SF_FUNCTION_NAME __PRETTY_FUNCTION__
#else
#  define SF_FUNCTION_NAME __func__
#endif

#define SF_VISIBILITY_EXTERNAL __attribute__((visibility ("default")))
#define SF_VISIBILITY_INTERNAL __attribute__((visibility ("hidden")))

#define	SF_INTERN extern SF_VISIBILITY_INTERNAL
#define SF_EXTERN extern SF_VISIBILITY_EXTERNAL

#define SF_INLINE        inline
#define	SF_STATIC_INLINE static SF_INLINE

#define SFUnlikely(pred) __builtin_expect(!!(pred),0)
#define SFLikely(pred)   __builtin_expect(!!(pred),1)

#if 0
// todo implement the detail:: parts
#ifdef SF_DEBUG
#  define SFCHKERR(expr) do {						\
    const auto errc = expr;						\
    if (SFUnlikely(detail::sf_check_error(errc))) {			\
      return detail::error_handler(__FILE__,SF_FUNCTION_NAME,		\
				   __LINE__,SF_COMM_WORLD,errc);	\
    }									\
  } while (0)
#else
#  define SFCHKERR(expr) expr
#endif
#endif

#include <mpi.h>

namespace sf
{

using sf_error_t = int;

SF_EXTERN MPI_Comm SF_COMM_WORLD;

} // namespace sf

#endif // SYMMFASTSYS_HPP
