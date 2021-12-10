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

  petsc_vector(MPI_Comm comm = SF_COMM_SELF, size_type ll = 0, size_type lg = -1) noexcept
    : base_type(comm,ll,lg)
  { }

  ~petsc_vector() noexcept { if (v_) SFCHECKABORT(VecDestroy(&v_)); }

  sf_error_t assemble() noexcept;

private:
  Vec v_;

  sf_error_t initialize_() noexcept;
};

} // namespace sf

#endif // SYMMFAST_LINALG_VECTOR_PETSC_VECTOR_HPP
