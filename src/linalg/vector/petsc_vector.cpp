#include <linalg/vector/petsc_vector.hpp>

#include <iostream>

namespace sf
{

sf_error_t petsc_vector::initialize_() noexcept
{
  if (vec_) return 0;
  SFCHECK(VecCreate(comm_,&vec_));
  SFCHECK(VecSetType(vec_,type_));
  SFCHECK(VecSetSizes(vec_,heightl_,heightg_));
  SFCHECK(VecSetFromOptions(vec_));
  return 0;
}

sf_error_t petsc_vector::set_random() noexcept
{
  SFCHECK(initialize_());
  SFCHECK(VecSetRandom(vec_,nullptr));
  return 0;
}

sf_error_t petsc_vector::assemble() noexcept
{
  SFCHECK(initialize_());
  SFCHECK(VecAssemblyBegin(vec_));
  SFCHECK(VecAssemblyEnd(vec_));
  assembled_ = true;
  return 0;
}

sf_error_t petsc_vector::resize(size_type new_length) noexcept
{
  heightl_   = new_length;
  assembled_ = false;
  return 0;
}

} // namespace sf
