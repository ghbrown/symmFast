#include <symmfast.hpp>

using namespace sf;

int main(int argc, char *argv[])
{
  SFCHECK(initialize(argc,argv));
  auto comm = SF_COMM_WORLD;

  auto pmat = petsc_matrix(comm,10,10,10,10);
  SFCHECK(finalize());
  return 0;
}
