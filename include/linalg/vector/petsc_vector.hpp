#ifndef SYMMFAST_LINALG_VECTOR_PETSC_VECTOR_HPP
#define SYMMFAST_LINALG_VECTOR_PETSC_VECTOR_HPP

#include <linalg/vector/vector_core.hpp>
#include <petscvec.h>

namespace sf
{

class SF_VISIBILITY_EXTERNAL petsc_vector;

class petsc_vector : public detail::vector_base<PetscScalar>
{
public:
  SF_LINEAR_OPERATOR_HEADER(base_type,detail::vector_base<PetscScalar>);

  petsc_vector(MPI_Comm comm = SF_COMM_SELF, VecType type = VECSTANDARD, size_type ll = 0, size_type lg = -1) noexcept
    : base_type(comm,ll,lg),type_(type)
  { }

  ~petsc_vector() noexcept { if (vec_) SFCHECKABORT(VecDestroy(&vec_)); }

  sf_error_t set_random() noexcept;
  sf_error_t assemble() noexcept;

private:
  Vec     vec_;
  VecType type_;

  sf_error_t initialize_() noexcept;
};

} // namespace sf

#endif // SYMMFAST_LINALG_VECTOR_PETSC_VECTOR_HPP
