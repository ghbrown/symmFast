#ifndef SYMMFAST_LINALG_MATRIX_PETSC_HPP
#define SYMMFAST_LINALG_MATRIX_PETSC_HPP

#include <linalg/matrix/matrix_core.hpp>
#include <petscmat.h>

namespace sf
{

extern template class detail::matrix_base<PetscScalar>;

class SF_VISIBILITY_EXTERNAL petsc_matrix;

class petsc_matrix : public detail::matrix_base<PetscScalar>
{
public:
  SF_LINEAR_OPERATOR_HEADER(base_type,detail::matrix_base<PetscScalar>);

  petsc_matrix(MPI_Comm comm = SF_COMM_SELF, MatType type = MATAIJ, size_type hl = 0, size_type wl = 0, size_type hg = -1, size_type wg = -1) noexcept
    : base_type(comm,hl,wl,hg,wg), mat_(nullptr), type_(type)
  { }

  virtual ~petsc_matrix() noexcept;
  virtual sf_error_t setrandom() noexcept;
  virtual sf_error_t assemble() noexcept;

private:
  Mat     mat_;
  MatType type_;

  sf_error_t initialize_() noexcept;
};

} // namespace sf

#endif // SYMMFAST_LINALG_MAT_PETSC_HPP
