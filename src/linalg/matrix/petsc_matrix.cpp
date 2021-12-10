#include <linalg/matrix/petsc_matrix.hpp>

namespace sf
{

petsc_matrix::~petsc_matrix() noexcept
{
  if (mat_) SFCHECKABORT(MatDestroy(&mat_));
}

sf_error_t petsc_matrix::initialize_() noexcept
{
  if (mat_) return 0;
  SFCHECK(MatCreate(comm_,&mat_));
  SFCHECK(MatSetSizes(mat_,heightl_,widthl_,heightg_,widthg_));
  SFCHECK(MatSetType(mat_,type_));
  SFCHECK(MatSetFromOptions(mat_));
  SFCHECK(MatSetUp(mat_));
  return 0;
}

sf_error_t petsc_matrix::setrandom() noexcept
{
  SFCHECK(initialize_());
  SFCHECK(MatSetRandom(mat_,nullptr));
  return 0;
}

sf_error_t petsc_matrix::assemble() noexcept
{
  SFCHECK(initialize_());
  SFCHECK(MatAssemblyBegin(mat_,MAT_FINAL_ASSEMBLY));
  SFCHECK(MatAssemblyEnd(mat_,MAT_FINAL_ASSEMBLY));
  return 0;
}


template class detail::matrix_base<PetscScalar>;

} // namespace sf
