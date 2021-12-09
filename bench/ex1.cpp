#include <symmfast.hpp>

using namespace sf;

int main(int argc, char *argv[])
{
  SFCHECK(initialize(argc,argv));
  auto comm = SF_COMM_WORLD;

  {
    auto pmat = petsc_matrix(comm,MATDENSE);
    pmat.setrandom();
    SFCHECK(pmat.assemble());
  }

  SFCHECK(finalize());
  return 0;
}
