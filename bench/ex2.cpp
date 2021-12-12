#include <petscvec.h>
#include <petsc/private/matimpl.h>

#include <type_traits>
#include <vector>
#include <iostream>

#define MATSYMMFAST "symmfast"

template <typename T>
static constexpr auto PetscObjectCast(T& obj)
{
  static_assert(std::is_pointer_v<T>);
  return reinterpret_cast<PetscObject>(obj);
}

template <typename T>
static constexpr auto PetscObjComm(T obj) { return PetscObjectComm(PetscObjectCast(obj)); }

class MatSymmFast
{
private:
  using array_type = std::vector<PetscScalar>;

  MPI_Comm   comm_;
  PetscInt   rbegin_ = PETSC_DECIDE;
  PetscInt   rend_   = PETSC_DECIDE;
  PetscInt   cbegin_ = PETSC_DECIDE;
  PetscInt   cend_   = PETSC_DECIDE;
  array_type data_   = {};

  static constexpr auto cast_(void* obj) noexcept { return static_cast<MatSymmFast*>(obj); }

  static constexpr auto impl_cast_(Mat m) noexcept { return cast_(m->data); }

  template <typename T>
  static auto destroy_(T *&obj) noexcept
  {
    PetscFunctionBegin;
    CHKERRCXX(delete cast_(obj));
    obj = nullptr;
    PetscFunctionReturn(0);
  }

  auto nrows_() const noexcept { return rend_-rbegin_; }
  auto ncols_() const noexcept { return cend_-cbegin_; }

  auto check_setup_() const noexcept
  {
    PetscFunctionBegin;
    if (PetscUnlikelyDebug(data_.size() && (nrows_()*ncols_() != data_.size()))) SETERRQ(comm_,PETSC_ERR_ARG_WRONGSTATE,"Mat was not setup");
    PetscFunctionReturn(0);
  }

  auto set_ownership_(PetscInt rbegin, PetscInt rend, PetscInt cbegin, PetscInt cend) noexcept
  {
    PetscFunctionBegin;
    rbegin_ = rbegin;
    rend_   = rend;
    cbegin_ = cbegin;
    cend_   = cend_;
    CHKERRCXX(data_.resize(nrows_()*ncols_()));
    PetscFunctionReturn(0);
  }

public:
  MatSymmFast(MPI_Comm comm) noexcept : comm_(comm) { }

  static PetscErrorCode create(Mat)                 noexcept;
  static PetscErrorCode destroy(Mat)                noexcept;
  static PetscErrorCode setup(Mat)                  noexcept;
  static PetscErrorCode set_random(Mat,PetscRandom) noexcept;
  static PetscErrorCode mat_mult(Mat,Vec,Vec)       noexcept;
  static PetscErrorCode view(Mat,PetscViewer)       noexcept;
};

template auto MatSymmFast::destroy_(MatSymmFast*&) noexcept;
template auto MatSymmFast::destroy_(void*&)        noexcept;

PetscErrorCode MatSymmFast::setup(Mat m) noexcept
{
  PetscFunctionBegin;
  const auto rmap = m->rmap;
  CHKERRQ(PetscLayoutSetUp(rmap));
  CHKERRQ(PetscIntView(rmap->size+1,rmap->range,PETSC_VIEWER_STDOUT_(PetscObjComm(m))));
  const auto cmap = m->cmap;
  CHKERRQ(PetscLayoutSetUp(cmap));
  PetscMPIInt size,rank;
  CHKERRMPI(MPI_Comm_size(PetscObjComm(m),&size));
  CHKERRMPI(MPI_Comm_rank(PetscObjComm(m),&rank));
  const auto total_size = rmap->N*cmap->N;
  const auto half_size  = rmap->N*cmap->N/2;
  const auto num_row_blocks = int(sqrt(rmap->N)+1);
  const auto num_col_blocks = int(sqrt(cmap->N)+1);
  const auto total_num_blocks = num_row_blocks*num_col_blocks;
  CHKERRQ(impl_cast_(m)->set_ownership_(rmap->rend-rmap->rstart,cmap->rend-cmap->rstart));
  CHKERRQ(PetscSynchronizedPrintf(PetscObjComm(m),"rend: %" PetscInt_FMT "\n",rmap->rend));
  CHKERRQ(PetscSynchronizedFlush(PetscObjComm(m),PETSC_STDOUT));
  PetscFunctionReturn(0);
}

PetscErrorCode MatSymmFast::destroy(Mat m) noexcept
{
  PetscFunctionBegin;
  CHKERRQ(destroy_(m->data));
  CHKERRQ(PetscObjectChangeTypeName(PetscObjectCast(m),nullptr));
  PetscFunctionReturn(0);
}

PetscErrorCode MatSymmFast::mat_mult(Mat m, Vec vin, Vec vout) noexcept
{
  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

PetscErrorCode MatSymmFast::set_random(Mat m, PetscRandom rand) noexcept
{
  PetscFunctionBegin;
  const auto msf = impl_cast_(m);
  CHKERRQ(msf->check_setup_());
  std::generate(msf->data_.begin(),msf->data_.end(),[=](){
    PetscScalar val;

    CHKERRABORT(PETSC_COMM_SELF,PetscRandomGetValue(rand,&val));
    return val;
  });
  PetscFunctionReturn(0);
}

PetscErrorCode MatSymmFast::create(Mat m) noexcept
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  m->ops->mult	    = mat_mult;
  m->ops->setrandom = set_random;
  m->ops->setup	    = setup;
  m->ops->destroy   = destroy;
  m->ops->view      = view;

  m->symmetric	       = PETSC_TRUE;
  m->hermitian	       = PETSC_TRUE;
  m->symmetric_eternal = PETSC_TRUE;
  CHKERRQ(MatSetOption(m,MAT_SYMMETRIC,PETSC_TRUE));
  CHKERRQ(MatSetOption(m,MAT_SYMMETRY_ETERNAL,PETSC_TRUE));
  CHKERRQ(PetscFree(m->data));
  m->data = new MatSymmFast(PetscObjComm(m));
  CHKERRQ(PetscObjectChangeTypeName(PetscObjectCast(m),MATSYMMFAST));
  PetscFunctionReturn(0);
}

PetscErrorCode MatSymmFast::view(Mat m, PetscViewer vwr) noexcept
{
  PetscFunctionBegin;
  const auto values = impl_cast_(m)->data_;
  CHKERRQ(PetscScalarView(values.size(),values.data(),PETSC_VIEWER_STDOUT_(PetscObjComm(m))));
  PetscFunctionReturn(0);
}

int main(int argc, char*argv[])
{
  PetscErrorCode ierr;
  MPI_Comm       comm;
  Mat            mat;
  Vec            vin,vout;

  ierr = PetscInitialize(&argc,&argv,nullptr,nullptr);if (PetscUnlikely(ierr)) return ierr;
  ierr = MatRegister(MATSYMMFAST,MatSymmFast::create);CHKERRQ(ierr);
  comm = PETSC_COMM_WORLD;

  ierr = MatCreate(comm,&mat);CHKERRQ(ierr);
  ierr = MatSetSizes(mat,10,10,PETSC_DECIDE,PETSC_DECIDE);CHKERRQ(ierr);
  ierr = MatSetType(mat,MATSYMMFAST);CHKERRQ(ierr);
  ierr = MatSetFromOptions(mat);CHKERRQ(ierr);
  ierr = MatSetUp(mat);CHKERRQ(ierr);
  ierr = MatSetRandom(mat,nullptr);CHKERRQ(ierr);

  ierr = MatCreateVecs(mat,&vin,&vout);CHKERRQ(ierr);
  ierr = VecSetFromOptions(vin);CHKERRQ(ierr);
  ierr = VecSetFromOptions(vout);CHKERRQ(ierr);
  ierr = VecSetRandom(vin,nullptr);CHKERRQ(ierr);
  ierr = VecSetRandom(vout,nullptr);CHKERRQ(ierr);

  ierr = MatMult(mat,vin,vout);CHKERRQ(ierr);

  ierr = VecDestroy(&vin);CHKERRQ(ierr);
  ierr = VecDestroy(&vout);CHKERRQ(ierr);
  ierr = MatDestroy(&mat);CHKERRQ(ierr);
  return PetscFinalize();
}
