#include <symmfast.hpp>
#include <iostream>

using namespace sf;

int main(int argc, char *argv[])
{
  SFCHECK(initialize(argc,argv));
  auto comm = SF_COMM_WORLD;

  {
    auto pmat = petsc_matrix(comm,MATDENSE);
    SFCHECK(pmat.setrandom());
    SFCHECK(pmat.assemble());
  }

  SFCHECK(finalize());
  std::cout<<"success\n";
  return 0;
}
