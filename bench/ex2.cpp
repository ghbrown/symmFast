#include <petscvec.h>
#include <petsc/private/matimpl.h>

#include <type_traits>
#include <vector>
#include <numeric>
#include <iostream>
#include <array>

#define MATSYMMFAST "symmfast"

template <typename T>
static constexpr auto PetscObjectCast(T& obj)
{
  static_assert(std::is_pointer_v<T>);
  return reinterpret_cast<PetscObject>(obj);
}

template <typename T>
static constexpr auto PetscObjComm(T obj) { return PetscObjectComm(PetscObjectCast(obj)); }

template <typename T>
static void vec_view(const std::vector<T>& v)
{
  for (auto&& x : v) std::cout<<x<<std::endl;
}

template <typename T, std::size_t N>
static void vec_view(const std::vector<std::array<T,N>>& v)
{
  for (auto&& x : v) std::cout<<std::get<0>(x)<<", "<<std::get<1>(x)<<std::endl;
}

template <typename F>
static auto generate_vector(std::size_t size, F&& fn) -> std::vector<decltype(fn())>
{
  std::vector<decltype(fn())> vec;
  vec.reserve(size);

  std::generate_n(std::back_inserter(vec),size,fn);
  vec.shrink_to_fit();
  return vec;
}

namespace impl
{

template <typename T, typename...> struct type_impl { using type = T; };
template <typename... T> struct type_impl<void,T...> : std::common_type<T...> { };

template <typename T, typename... Rest>
using array_type = std::array<typename type_impl<T,Rest...>::type,sizeof...(Rest)>;

} // namespace impl

template <typename R = void, typename... T>
static impl::array_type<R,T...> make_array(T&&... args) noexcept
{
  return {std::forward<T>(args)...};
}

class MatSymmFast
{
private:
  using array_type = std::vector<PetscScalar>;
  enum class DataOrder {ROW_MAJOR,COLUMN_MAJOR};

  static const auto order_ = DataOrder::ROW_MAJOR;

  const Mat  m_;
  MPI_Comm   row_comm_;
  MPI_Comm   col_comm_;
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

  auto on_diagonal_() const noexcept { return rbegin_ == cbegin_ && rend_ == cend_; }

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
  const auto  msf  = impl_cast_(m);
  const auto  cmap = m->cmap;
  PetscMPIInt size,rank;
  PetscInt    col_color,row_color;
  auto        r = PetscInt(3);

  PetscFunctionBegin;
  CHKERRMPI(MPI_Comm_size(PetscObjComm(m),&size));
  CHKERRMPI(MPI_Comm_rank(PetscObjComm(m),&rank));
  CHKERRQ(PetscLayoutSetUp(m->rmap));
  CHKERRQ(PetscLayoutSetUp(cmap));

  PetscObjectOptionsBegin(PetscObjectCast(m));
  CHKERRQ(PetscOptionsInt("-mat_symmfast_r","num communicators per dim","Mat",r,&r,nullptr));
  PetscOptionsEnd();
  if (PetscUnlikely(size != r*(r+1)/2)) SETERRQ2(PetscObjComm(m),PETSC_ERR_ARG_SIZ,"Communicator of size %d not %d",size,r*(r+1)/2);

  if (rank) {
    auto recv_buf = make_array(-1,-1,-1,-1,-1,-1);

    CHKERRMPI(
      MPI_Recv(recv_buf.data(),recv_buf.size(),MPIU_INT,0,0,PetscObjComm(m),MPI_STATUS_IGNORE)
    );
    for (auto&& x: recv_buf) if (PetscUnlikely(x == -1)) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_MPI,"MPI_Recv from root failed");

    std::tie(msf->rbegin_,msf->rend_,msf->cbegin_,msf->cend_,col_color,row_color) = recv_buf;
  } else {
    const auto block_size = std::max(cmap->N/r,1);
    auto       bounds     = std::vector<PetscInt>();

    CHKERRCXX(bounds.reserve(r));
    for (auto i = 0,sum = 0; i < r; ++i) {
      CHKERRCXX(bounds.insert(bounds.cend(),{sum,(i == r-1 ? cmap->N : sum+=block_size)-1}));
    }
    vec_view(bounds);
    for (auto i = 1,cur_row = 0,cur_col = 1; i < size; ++i) {
      const auto row_idx  = 2*i+cur_row;
      const auto col_idx  = 2*i+cur_row;
      const auto send_buf = make_array(
        bounds[row_idx],bounds[row_idx+1],bounds[col_idx],bounds[col_idx+1],cur_row,cur_col
      );

      CHKERRMPI(MPI_Send(send_buf.data(),send_buf.size(),MPIU_INT,i,0,PetscObjComm(m)));
      if (++cur_col > r-1) cur_col = ++cur_row;
    }

    // rank 0 always know its ownership
    CHKERRQ(msf->set_ownership_(0,bounds[1],0,bounds[1]));
    col_color = row_color = 0;
  }
  CHKERRMPI(MPI_Comm_split(PetscObjComm(m),col_color,0,&msf->col_comm_));
  CHKERRMPI(MPI_Comm_split(PetscObjComm(m),row_color,0,&msf->row_comm_));

  const auto n = msf->cend_-msf->cbegin_;
  CHKERRCXX(msf->data_.resize(msf->on_diagonal_() ? n*(n+1)/2 : n*(msf->rend_-msf->rbegin_)));
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
  // A \in N x M
  // A[n:n+bs,m:m+bs]
  // _______
  // |xxxxxx| ""
  // | xxxxx| ncols -> 5
  // |  xxxx| ncols -> 4
  // |   xxx|
  // |    xx|
  // |     x|

  // for row in [0 ... n]
  //   for col in ncols(row)
  //
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
  // only diagonal blocks own data
  if (!msf->on_diagonal_()) PetscFunctionReturn(0);
  CHKERRCXX(std::generate(msf->data_.begin(),msf->data_.end(),[=]{
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
