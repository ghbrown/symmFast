#include <linalg/petsc.hpp>

namespace sf
{

petsc_matrix::~petsc_matrix() noexcept
{
  SFCHECKABORT(MatDestroy(&mat_));
}

sf_error_t petsc_matrix::initialize_() noexcept
{
  if (!mat_) SFCHECK(MatCreate(comm_,&mat_));
  SFCHECK(MatSetSizes(mat_,heightl_,widthl_,heightg_,widthg_));
  return 0;
}

sf_error_t petsc_matrix::assemble() noexcept
{
  SFCHECK(initialize_());
  SFCHECK(MatSetType(mat_,type_));
  SFCHECK(MatSetFromOptions(mat_));
  SFCHECK(MatAssemblyBegin(mat_,MAT_FINAL_ASSEMBLY));
  SFCHECK(MatAssemblyEnd(mat_,MAT_FINAL_ASSEMBLY));
  return 0;
}

petsc_matrix::linear_operator& petsc_matrix::apply(const linear_operator &other) noexcept
{
  return *this;
}

} // namespace sf
