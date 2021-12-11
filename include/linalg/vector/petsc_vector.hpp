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

  petsc_vector(MPI_Comm comm = SF_COMM_SELF, VecType type = VECSTANDARD, size_type len_local = 0, size_type len_global = -1) noexcept
    : base_type(comm,len_local,len_global),type_(type)
  { }

  ~petsc_vector() noexcept { if (vec_) SFCHECKABORT(VecDestroy(&vec_)); }

  sf_error_t set_random() noexcept;
  sf_error_t assemble() noexcept;

  sf_error_t resize(size_type) noexcept;

private:
  Vec     vec_ = nullptr;
  VecType type_;

  sf_error_t initialize_() noexcept;
};

} // namespace sf

#endif // SYMMFAST_LINALG_VECTOR_PETSC_VECTOR_HPP
