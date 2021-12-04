#ifndef SF_LINALG_SPARSE_MATRIX_HPP
#define SF_LINALG_SPARSE_MATRIX_HPP

#include <linalg/matrix.hpp>

namespace sf
{

template <typename T>
class SF_VISIBILITY_EXTERNAL sparse_matrix : public matrix<T>
{
public:
  SF_LINEAR_OPERATOR_HEADER(base_type,matrix<T>);

  sparse_matrix(MPI_Comm comm = SF_COMM_SELF, size_type hl = 0, size_type wl = 0, size_type hg = -1, size_type wg = -1) noexcept
    : base_type(comm,hl,wl,hg,wg)
  { }

  virtual ~sparse_matrix() noexcept = default;
  // construct the matrix
  virtual sf_error_t assemble() noexcept;
  virtual linear_operator<T>& apply(const linear_operator<T>&) noexcept;
};

template <typename T>
sf_error_t sparse_matrix<T>::assemble() noexcept
{
  return 0;
}

template <typename T>
linear_operator<T>& sparse_matrix<T>::apply(const linear_operator<T>&) noexcept
{
  return *this;
}

}  // namespace sf

#endif // SF_LINALG_SPARSE_MATRIX_HPP