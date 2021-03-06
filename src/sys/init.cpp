#include <sys/sys.hpp>
#include <petscsys.h>

namespace sf
{

MPI_Comm SF_COMM_WORLD = MPI_COMM_NULL;
MPI_Comm SF_COMM_SELF  = MPI_COMM_NULL;

sf_error_t initialize(int argc, char *argv[]) noexcept
{
  SFCHECK(PetscInitialize(&argc,&argv,nullptr,nullptr));
  SF_COMM_SELF  = PETSC_COMM_SELF;
  SF_COMM_WORLD = PETSC_COMM_WORLD;
  return 0;
}

sf_error_t finalize() noexcept
{
  SFCHECK(PetscFinalize());
  SF_COMM_SELF  = MPI_COMM_NULL;
  SF_COMM_WORLD = MPI_COMM_NULL;
  return 0;
}

} // namespace sf
