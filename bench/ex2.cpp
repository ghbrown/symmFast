#include <petscvec.h>
#include <petsc/private/matimpl.h>

#include <type_traits>
#include <vector>
#include <numeric>
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
  enum class DataOrder {ROW_MAJOR,COLUMN_MAJOR};

  static const auto order_ = DataOrder::ROW_MAJOR;

  const Mat  m_;
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

  auto ncols_(PetscInt i) const noexcept { return cend_ - (i < m_->rmap->rstart ? cbegin_ : i); }
  auto nrows_(PetscInt i) const noexcept { return i; }

  auto check_setup_() const noexcept
  {
    PetscFunctionBegin;
    if (PetscUnlikelyDebug(!data_.size())) SETERRQ(PetscObjComm(m_),PETSC_ERR_ARG_WRONGSTATE,"Mat was not setup");
    PetscFunctionReturn(0);
  }

  auto set_ownership_(PetscInt rbegin, PetscInt rend, PetscInt cbegin, PetscInt cend) noexcept
  {
    PetscFunctionBegin;
    rbegin_ = rbegin;
    rend_   = rend;
    cbegin_ = cbegin;
    cend_   = cend;
    PetscFunctionReturn(0);
  }

  auto& operator()(PetscInt i, PetscInt j) noexcept
  {
    return data_[order_ == DataOrder::ROW_MAJOR ? (i*ncols_(i))+j : (j*nrows_(j))+i];
  }

  const auto& operator()(PetscInt i, PetscInt j) const noexcept
  {
    return data_[order_ == DataOrder::ROW_MAJOR ? (i*ncols_(i))+j : (j*nrows_(j))+i];
  }

public:
  MatSymmFast(Mat m) noexcept : m_(m) { }

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
  const auto rmap = m->rmap;
  const auto cmap = m->cmap;
  const auto msf  = impl_cast_(m);

  PetscFunctionBegin;
  CHKERRQ(PetscLayoutSetUp(rmap));
  CHKERRQ(PetscLayoutSetUp(cmap));
  CHKERRQ(PetscIntView(rmap->size+1,rmap->range,PETSC_VIEWER_STDOUT_(PetscObjComm(m))));
  CHKERRQ(msf->set_ownership_(0,rmap->N,cmap->rstart,cmap->rend));
  // own a whole column each
  const auto begin = cmap->rstart+1;
  const auto end   = cmap->rend+1;
  CHKERRCXX(msf->data_.resize((begin+end)*(end-begin)/2));
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
  const auto&       msf = *impl_cast_(m);
  const auto        N = m->cmap->N;
  const auto        nrow = msf.rend_-msf.rbegin_;
  const PetscScalar *array_in;
  PetscScalar       *array_out;

  PetscFunctionBegin;
  CHKERRQ(VecGetArrayRead(vin,&array_in));
  CHKERRQ(VecGetArrayWrite(vout,&array_out));

  // compute local row-sums
  auto rowsums = std::vector<PetscScalar>(m->cmap->N,0);
  std::generate_n(std::back_inserter(rowsums),nrow,[&,it=msf.data_.cbegin(),i=0]() mutable {
    const auto begin = it;
    it += msf.ncols_(i++);
    return std::accumulate(begin,it,0);
  });

  // all reduce for global row sums
  CHKERRMPI(
    MPI_Allreduce(rowsums.data(),rowsums.data(),rowsums.size(),MPIU_SCALAR,MPI_SUM,PetscObjComm(m))
  );

  for (auto i = 0; i < nrow; ++i) {
    auto       lhs_sum = PetscScalar(0);
    auto       rhs_sum = PetscScalar(0);
    const auto bi      = array_in[i];

    // compute local row-sums of A

    for (auto k = i; k < N; ++k) {
      const auto aik = msf(i,k);
      const auto bk  = array_in[k];

      lhs_sum += aik*(bi+bk);
      rhs_sum += aik;
    }
    array_out[i] = lhs_sum-(rhs_sum*bi);
  }
  CHKERRQ(VecRestoreArrayWrite(vout,&array_out));
  CHKERRQ(VecRestoreArrayRead(vin,&array_in));
  PetscFunctionReturn(0);
}

PetscErrorCode MatSymmFast::set_random(Mat m, PetscRandom rand) noexcept
{
  const auto msf = impl_cast_(m);

  PetscFunctionBegin;
  CHKERRQ(msf->check_setup_());
  CHKERRCXX(std::generate(msf->data_.begin(),msf->data_.end(),[=](){
    PetscScalar val;

    CHKERRABORT(PETSC_COMM_SELF,PetscRandomGetValue(rand,&val));
    return val;
  }));
  PetscFunctionReturn(0);
}

PetscErrorCode MatSymmFast::create(Mat m) noexcept
{
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
  m->data = new MatSymmFast{m};
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
  Mat            mat;
  Vec            vin,vout;

  ierr = PetscInitialize(&argc,&argv,nullptr,nullptr);if (PetscUnlikely(ierr)) return ierr;
  CHKERRQ(MatRegister(MATSYMMFAST,MatSymmFast::create));

  auto comm = PETSC_COMM_WORLD;
  CHKERRQ(MatCreate(comm,&mat));
  CHKERRQ(MatSetSizes(mat,10,10,PETSC_DECIDE,PETSC_DECIDE));
  CHKERRQ(MatSetType(mat,MATSYMMFAST));
  CHKERRQ(MatSetFromOptions(mat));
  CHKERRQ(MatSetUp(mat));
  CHKERRQ(MatSetRandom(mat,nullptr));

  CHKERRQ(MatCreateVecs(mat,&vin,&vout));
  CHKERRQ(VecSetFromOptions(vin));
  CHKERRQ(VecSetFromOptions(vout));
  CHKERRQ(VecSetRandom(vin,nullptr));
  CHKERRQ(VecZeroEntries(vout));

  CHKERRQ(VecViewFromOptions(vout,nullptr,"-result_view"));
  CHKERRQ(MatMult(mat,vin,vout));
  CHKERRQ(VecViewFromOptions(vout,nullptr,"-result_view"));

  CHKERRQ(VecDestroy(&vin));
  CHKERRQ(VecDestroy(&vout));
  CHKERRQ(MatDestroy(&mat));
  ierr = PetscFinalize();
  return ierr;
}
