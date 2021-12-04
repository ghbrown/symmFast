#ifndef SYMMFAST_LINALG_MATRIX_HPP
#define SYMMFAST_LINALG_MATRIX_HPP

#include <sys/sys.hpp>
#include <linalg/linear_operator.hpp>

namespace sf
{

// a simple matrix abstract base class
template <typename T>
class SF_VISIBILITY_EXTERNAL matrix : public linear_operator<T>
{
public:
  SF_LINEAR_OPERATOR_HEADER(base_type,linear_operator<T>);

  // constructor
  matrix(MPI_Comm comm = SF_COMM_SELF, size_type hl = 0, size_type wl = 0, size_type hg = -1, size_type wg = -1) noexcept
    : base_type(comm,hl,wl,hg,wg)
  { }

  virtual ~matrix() noexcept = default;
  virtual sf_error_t assemble() noexcept = 0;
  virtual linear_operator<T>& apply(const linear_operator<T>&) noexcept = 0;
};

} // namespace sf

#endif // SYMMFAST_LINALG_MATRIX_HPP
