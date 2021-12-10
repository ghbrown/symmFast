#ifndef SYMMFAST_LINALG_VECTOR_CORE_HPP
#define SYMMFAST_LINALG_VECTOR_CORE_HPP

#include <sys/sys.hpp>
#include <linalg/linear_operator.hpp>

namespace sf
{

template <typename T> class SF_VISIBILITY_EXTERNAL vector;

template <typename T>
class vector : public linear_operator<T>
{
public:

};

}

#endif // SYMMFAST_LINALG_VECTOR_CORE_HPP
