#include <symmfast.hpp>

using namespace sf;

int main(int argc, char *argv[])
{
  SFCHECK(initialize(argc,argv));

  SFCHECK(finalize());
  return 0;
}
