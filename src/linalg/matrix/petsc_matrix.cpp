#include <linalg/matrix/petsc_matrix.hpp>

namespace sf
{

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

sf_error_t petsc_matrix::set_random() noexcept
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
  assembled_ = true;
  return 0;
}

sf_error_t petsc_matrix::resize(size_type new_height, size_type new_width) noexcept
{
  assembled_ = false;
  heightl_   = new_height;
  widthl_    = new_width;
  return 0;
}

sf_error_t petsc_matrix::apply(const petsc_vector& vec, petsc_vector& result) const noexcept
{
  PetscBool assembled;

  SFCHECK(MatAssembled(mat_,&assembled));
  SFASSERT(assembled == PETSC_TRUE,SF_COMM_SELF);
  SFASSERT(vec.assembled(),SF_COMM_WORLD);
  SFCHECK(result.resize(heightl_));

  return 0;
}


template class detail::matrix_base<PetscScalar>;

} // namespace sf
