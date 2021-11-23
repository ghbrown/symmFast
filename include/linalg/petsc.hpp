#ifndef SYMMFAST_LINALG_PETSC_HPP
#define SYMMFAST_LINALG_PETSC_HPP

#include <linalg/matrix.hpp>
#include <petscmat.h>

namespace sf
{

class petsc_matrix;

class petsc_matrix : public sparse_matrix<PetscScalar>
{
public:
  SF_LINEAR_OPERATOR_HEADER(linop,PetscScalar);
  using base_type = sparse_matrix<PetscScalar>;

  petsc_matrix(MPI_Comm comm = SF_COMM_SELF, size_type hl = 0, size_type wl = 0, size_type hg = -1, size_type wg = -1)
    : base_type(comm,hl,wl,hg,wg), mat_(nullptr), type_(MATAIJ)
  { }

  virtual ~petsc_matrix() noexcept;

  sf_error_t assemble() noexcept override;

  linop& apply(const linop&) noexcept override;

private:
  Mat     mat_;
  MatType type_;

  sf_error_t initialize_() noexcept;
};

} // namespace sf

#endif // SYMMFAST_LINALG_PETSC_HPP
