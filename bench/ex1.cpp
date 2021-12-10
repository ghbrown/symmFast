#include <symmfast.hpp>
#include <iostream>

using namespace sf;

int main(int argc, char *argv[])
{
  SFCHECK(initialize(argc,argv));
  auto comm = SF_COMM_WORLD;

  {
    auto pmat = petsc_matrix(comm,MATDENSE);
    SFCHECK(pmat.set_random());
    SFCHECK(pmat.assemble());

    auto pvec = petsc_vector(comm);
    SFCHECK(pvec.set_random());
    SFCHECK(pvec.assemble());
  }

  SFCHECK(finalize());
  std::cout<<"success\n";
  return 0;
}
