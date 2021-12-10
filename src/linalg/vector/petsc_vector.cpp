#include <linalg/vector/petsc_vector.hpp>

namespace sf
{

sf_error_t petsc_vector::initialize_() noexcept
{
  if (vec_) return 0;
  SFCHECK(VecCreate(comm_,&vec_));
  SFCHECK(VecSetType(vec_,type_));
  return 0;
}

sf_error_t petsc_vector::assemble() noexcept
{
  SFCHECK(initialize_());
  return 0;
}

} // namespace sf
