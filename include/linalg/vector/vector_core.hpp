#ifndef SYMMFAST_LINALG_VECTOR_CORE_HPP
#define SYMMFAST_LINALG_VECTOR_CORE_HPP

#include <sys/sys.hpp>
#include <linalg/linear_operator.hpp>

#include <vector>

namespace sf
{

namespace detail
{

template <typename T> class SF_VISIBILITY_EXTERNAL vector_base;

template <typename T>
class vector_base : public linear_operator<T>
{
public:
  SF_LINEAR_OPERATOR_HEADER(base_type,linear_operator<T>);

  vector_base(MPI_Comm comm = SF_COMM_SELF, size_type len_local = 0, size_type len_global = -1) noexcept
    : base_type(comm,len_local,1,len_global,1)
  { }

  virtual ~vector_base() noexcept = default;
  virtual sf_error_t assemble() noexcept = 0;
};

} // namespace detail

} // namespace sf

#endif // SYMMFAST_LINALG_VECTOR_CORE_HPP
