#ifndef SYMMFAST_LINALG_MATRIX_HPP
#define SYMMFAST_LINALG_MATRIX_HPP

#include <sys/sys.hpp>
#include <linalg/linear_operator.hpp>

#include <vector>

namespace sf
{

// a simple matrix abstract base class
template <typename T>
class matrix : public linear_operator<T>
{
public:
  SF_LINEAR_OPERATOR_HEADER(linop,T);

  // constructor
  matrix(MPI_Comm comm = SF_COMM_SELF, size_type hl = 0, size_type wl = 0, size_type hg = -1, size_type wg = -1) noexcept
    : linop(comm,hl,wl,hg,wg)
  { }

  virtual ~matrix() noexcept = default;

  // construct the matrix
  virtual sf_error_t assemble() noexcept = 0;
};

// template <typename T>
// class dense_matrix : public matrix<T>
// {
// public:
//   SF_LINEAR_OPERATOR_HEADER(linop,T);
//   using base_type = matrix<T>;

//   constexpr dense_matrix(MPI_Comm comm = SF_COMM_SELF, size_type h = 0, size_type w = 0) noexcept
//     : base_type(comm,h,w), data_(h*w)
//   { }

//   // getter
//   const_reference_type operator()(size_type i, size_type j) const noexcept final
//   {
//     return data_[i*this->width_+j];
//   }

//   reference_type operator()(size_type i, size_type j) noexcept final
//   {
//     return data_[i*this->width_+j];
//   }

//   virtual ~dense_matrix() noexcept { }

// private:
//   std::vector<value_type> data_;
// };

template <typename T>
class sparse_matrix : public matrix<T>
{
public:
  SF_LINEAR_OPERATOR_HEADER(linop,T);
  using base_type = matrix<T>;

  sparse_matrix(MPI_Comm comm = SF_COMM_SELF, size_type hl = 0, size_type wl = 0, size_type hg = -1, size_type wg = -1) noexcept
    : base_type(comm,hl,wl,hg,wg)
  { }

  virtual ~sparse_matrix() noexcept;

  // construct the matrix
  virtual sf_error_t assemble() noexcept override;

  virtual linop& apply(const linop&) noexcept override;
};

template <typename T>
inline sparse_matrix<T>::~sparse_matrix() noexcept
{

}

template <typename T>
inline sf_error_t sparse_matrix<T>::assemble() noexcept
{
  return 0;
}

template <typename T>
inline linear_operator<T>& sparse_matrix<T>::apply(const linear_operator<T>&) noexcept
{
  return *this;
}

}

#endif // SYMMFAST_LINALG_MATRIX_HPP
