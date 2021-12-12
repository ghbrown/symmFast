#include <symmfast.hpp>
#include <iostream>

using namespace sf;

int main(int argc, char *argv[])
{
  SFCHECK(initialize(argc,argv));
  auto comm = SF_COMM_WORLD;

  {
    auto pmat = petsc_matrix(comm,MATDENSE,10,10);
    SFCHECK(pmat.set_random());
    SFCHECK(pmat.assemble());

    auto pvec = petsc_vector(comm,VECSTANDARD,pmat.local_height());
    SFCHECK(pvec.set_random());
    SFCHECK(pvec.assemble());

    auto pvec_out = petsc_vector(comm,VECSTANDARD);

    pmat.apply(pvec,pvec_out);
  }

  SFCHECK(finalize());
  std::cout<<"success\n";
  return 0;
}
