#ifndef SYMMFAST_LINALG_LINEAR_OPERATOR_HPP
#define SYMMFAST_LINALG_LINEAR_OPERATOR_HPP

#include <sys/sys.hpp>
#include <type_traits>

namespace sf
{

// abstract base class from which all other linear algebra types will derive from. a
// linear_operator has both a width (dimension of the input space) and a height (dimension of
// the output space)
template <typename T>
class SF_VISIBILITY_EXTERNAL linear_operator
{
  static_assert(std::is_floating_point<T>::value,"");
public:
  using value_type	     = T;
  using reference_type	     = value_type&;
  using const_reference_type = const value_type&;
  using size_type	     = std::size_t;

  linear_operator(MPI_Comm comm = SF_COMM_SELF, size_type hl = 0, size_type wl = 0, size_type hg = -1, size_type wg = -1) noexcept
    : comm_(comm), heightl_(hl), heightg_(hg), widthl_(wl), widthg_(wg)
  { }

  virtual ~linear_operator() noexcept = default;

  // assemble this linear operator and make it ready for use
  virtual sf_error_t assemble() noexcept = 0;

  // apply this linear operator, resulting in another linear operator
  virtual linear_operator& apply(const linear_operator&) noexcept = 0;

protected:
  MPI_Comm  comm_;
  size_type heightl_;
  size_type heightg_;
  size_type widthl_;
  size_type widthg_;
};

} // namespace sf

#define SF_LINEAR_OPERATOR_HEADER(base_name,base_class_T)	\
  using base_name = base_class_T;				\
  /* values */							\
  using base_name::comm_;					\
  using base_name::heightl_;					\
  using base_name::heightg_;					\
  using base_name::widthl_;					\
  using base_name::widthg_;					\
  /* typedefs */						\
  using typename base_name::value_type;				\
  using typename base_name::reference_type;			\
  using typename base_name::const_reference_type;		\
  using typename base_name::size_type

#endif // SYMMFAST_LINALG_LINEAR_OPERATOR_HPP
