#ifndef SYMMFAST_LINALG_MATRIX_HPP
#define SYMMFAST_LINALG_MATRIX_HPP

#include <type_traits>

namespace sf
{

// abstract base class from which all other linear algebra types will derive from. a
// linear_operator has both a width (dimension of the input space) and a height (dimension of
// the output space)
class linear_operator
{
  using size_type = std::size_t;

  constexpr linear_operator(size_type h, size_type w = 1) noexcept : _height(h), _width(w) { }

  // apply this linear operator, resulting in another linear operator
  ~virtual linear_operator& apply(const linear_operator&) noexcept const = 0;

  ~virtual linear_operator() { }

protected:
  size_type _height;
  size_type _width;
};

// a simple matrix abstract base class
template <typename T>
class matrix : public linear_operator
{
  static_assert(std::is_floating_point<T>::value,"");
public:
  using value_type = T;
  using linear_operator::size_type;

  // square matrix constructor
  constexpr matrix(size_type h) noexcept : linear_operator(h,h) { }

  // generic matrix constructor
  constexpr matrix(size_type h, size_type w) noexcept : linear_operator(h,w) { }

  virtual ~matrix() { }

  // getter
  virtual const value_type& operator()(size_type,size_type) const noexcept = 0;
  // setter
  virtual value_type& operator()(size_type,size_type) noexcept = 0;

  // construct the matrix
  virtual sf_error_t assemble() noexcept = 0;
};

template <typename T>
class dense_matrix : public matrix<T>
{
  // TODO
};

template <typename T>
class sparse_matrix : public matrix<T>
{
  // TODO
};

}

#endif // SYMMFAST_LINALG_MATRIX_HPP
