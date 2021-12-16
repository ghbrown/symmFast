#include <petsc/private/vecimpl.h>
#include <petsc/private/matimpl.h>

#include <type_traits>
#include <vector>
#include <numeric>
#include <iostream>
#include <array>

#define MATSYMMFAST "symmfast"

template <typename T>
static constexpr auto PetscObjectCast(T& obj) noexcept
{
  static_assert(std::is_pointer_v<T>);
  return reinterpret_cast<PetscObject>(obj);
}

template <typename T>
static constexpr auto PetscObjComm(T o) noexcept { return PetscObjectComm(PetscObjectCast(o)); }

template <typename T>
static void vec_view(const std::vector<T>& v) noexcept
{
  for (auto&& x : v) std::cout<<x<<std::endl;
}

template <typename T, std::size_t N>
static void vec_view(const std::vector<std::array<T,N>>& v) noexcept
{
  for (auto&& x : v) std::cout<<std::get<0>(x)<<", "<<std::get<1>(x)<<std::endl;
}

template <typename F>
static auto generate_vector(std::size_t size, F&& fn) noexcept -> std::vector<decltype(fn())>
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
  MPI_Comm   row_comm_ = MPI_COMM_NULL;
  MPI_Comm   col_comm_ = MPI_COMM_NULL;
  PetscInt   rbegin_   = PETSC_DECIDE;
  PetscInt   rend_     = PETSC_DECIDE;
  PetscInt   cbegin_   = PETSC_DECIDE;
  PetscInt   cend_     = PETSC_DECIDE;
  array_type data_     = {};

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

  auto ncols_(PetscInt row) const noexcept { return cend_-cbegin_-on_diagonal_()*row; }
  auto nrows_(PetscInt col) const noexcept { return rend_-rbegin_-on_diagonal_()*col; }
  auto ncols_sum_(PetscInt row) const noexcept
  {
    return row*(on_diagonal_() ? (ncols_(row)+(cend_-cbegin_))/2 : ncols_(row));
  }

  auto nrows_sum_(PetscInt col) const noexcept
  {
    return col*(on_diagonal_() ? (1+nrows_(col))/2 : nrows_(col));
  }

  auto check_setup_() const noexcept
  {
    PetscFunctionBegin;
    if (PetscUnlikelyDebug(row_comm_ != MPI_COMM_NULL)) SETERRQ(PetscObjComm(m_),PETSC_ERR_ARG_WRONGSTATE,"Mat was not setup");
    PetscFunctionReturn(0);
  }

  auto set_ownership_(PetscInt rbegin, PetscInt rend, PetscInt cbegin, PetscInt cend) noexcept
  {
    PetscFunctionBegin;
    rbegin_ = rbegin;
    rend_   = rend;
    cbegin_ = cbegin;
    cend_   = cend;
    data_.resize((rend_-rbegin_)*(cend_-cbegin_));
    PetscFunctionReturn(0);
  }

  auto& operator()(PetscInt i, PetscInt j) noexcept
  {
    return data_[i*(cend_-cbegin_)+j];
  }

  const auto& operator()(PetscInt i, PetscInt j) const noexcept
  {
    return data_[i*(cend_-cbegin_)+j];
  }

public:
  MatSymmFast(Mat m) noexcept : m_(m) { }

  static PetscErrorCode create(Mat)                 noexcept;
  static PetscErrorCode destroy(Mat)                noexcept;
  static PetscErrorCode setup(Mat)                  noexcept;
  static PetscErrorCode set_random(Mat,PetscRandom) noexcept;
  static PetscErrorCode mat_mult(Mat,Vec,Vec)       noexcept;
  static PetscErrorCode view(Mat,PetscViewer)       noexcept;
  static PetscErrorCode get_vecs(Mat,Vec*,Vec*)     noexcept;
};

template auto MatSymmFast::destroy_(MatSymmFast*&) noexcept;
template auto MatSymmFast::destroy_(void*&)        noexcept;

PetscErrorCode MatSymmFast::setup(Mat m) noexcept
{
  const auto  msf  = impl_cast_(m);
  const auto  cmap = m->cmap;
  PetscMPIInt size,rank;
  PetscInt    col_color,row_color;
  auto        r = 3;

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
      MPI_Recv(recv_buf.data(),recv_buf.size(),MPI_INT,0,0,PetscObjComm(m),MPI_STATUS_IGNORE)
    );
    for (auto&& x: recv_buf) if (PetscUnlikely(x == -1)) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_MPI,"MPI_Recv from root failed");

    const auto [rb,re,cb,ce,rowc,colc] = recv_buf;
    row_color = rowc;
    col_color = colc;
    CHKERRQ(msf->set_ownership_(rb,re,cb,ce));
  } else {
    const auto N          = cmap->N;
    const auto block_size = std::max(N/r,1);
    auto       bounds     = std::vector<PetscInt>{};

    CHKERRCXX(bounds.reserve(2*r));
    for (auto i = 0,sum = 0; i < r; ++i) CHKERRCXX(
      bounds.insert(bounds.cend(),{sum,i == r-1 ? N : sum+=block_size})
    );

    for (auto i = 1,cur_row = 0,cur_col = 1; i < size; ++i) {
      const auto row_idx  = 2*cur_row;
      const auto col_idx  = 2*cur_col;
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
  CHKERRMPI(MPI_Comm_split(PetscObjComm(m),row_color,rank,&msf->row_comm_));
  CHKERRMPI(MPI_Comm_split(PetscObjComm(m),col_color,rank%r,&msf->col_comm_));
  PetscFunctionReturn(0);
}

PetscErrorCode MatSymmFast::destroy(Mat m) noexcept
{
  const auto msf = impl_cast_(m);

  PetscFunctionBegin;
  if (msf->col_comm_ != MPI_COMM_NULL) CHKERRMPI(MPI_Comm_free(&msf->col_comm_));
  if (msf->row_comm_ != MPI_COMM_NULL) CHKERRMPI(MPI_Comm_free(&msf->row_comm_));
  CHKERRQ(destroy_(m->data));
  CHKERRQ(PetscObjectChangeTypeName(PetscObjectCast(m),nullptr));
  PetscFunctionReturn(0);
}

PetscErrorCode MatSymmFast::mat_mult(Mat m, Vec vin, Vec vout) noexcept
{
  //computes matrix vector product m @ vin = vout
  const auto&  msf  = *impl_cast_(m);
  const auto   nrow = msf.rend_ - msf.rbegin_; //local row block size
  const auto   ncol = msf.cend_ - msf.cbegin_; //local column block size
  PetscScalar *array_out = nullptr;
  PetscScalar *b_row,*b_col;
  PetscMPIInt rsize,csize,rrank,crank;

  PetscFunctionBegin;
  if (msf.on_diagonal_()) {
    CHKERRQ(VecGetArrayRead(vin,const_cast<const PetscScalar**>(&b_row)));
    CHKERRQ(VecGetArrayWrite(vout,&array_out));
    b_col = b_row;
  } else {
    b_row = new PetscScalar[nrow+ncol];
    b_col = b_row+nrow;
  }

  // get proper slices of b on local rank via broadcast from diagonal to column communicators
  // and row communicators
  CHKERRMPI(MPI_Comm_size(msf.row_comm_,&rsize));
  CHKERRMPI(MPI_Comm_size(msf.row_comm_,&csize));
  if (rsize) CHKERRMPI(MPI_Bcast(b_row,nrow,MPIU_SCALAR,0,msf.row_comm_));
  if (csize) CHKERRMPI(MPI_Bcast(b_col,ncol,MPIU_SCALAR,0,msf.col_comm_));

  // local updates:
  // locals[0:nrow] -> updates for A row sum row comm
  // locals[nrow:ncol] -> updates for Z row sum row comm
  // locals[nrow+ncol:nrow+ncol+nrow] -> updates for A row sum col comm
  // locals[nrow+ncol+nrow:-1] -> updates for Z row sum col comm
  auto local = std::vector<PetscScalar>(nrow+nrow+ncol+ncol,0);
  const auto rc_size = nrow+nrow;

  for (auto i = 0; i < nrow; ++i) {
    const auto bi = b_row[i];

    for (auto k = 0; k < ncol; ++k) {
      const auto aik     = msf(i,k);
      const auto aikbibk = aik*(bi+b_col[k]);

      // updates to row sum of A - row comm
      local[i] += aik;

      // updates to row sum of Z - row comm
      local[nrow+i] += aikbibk;

      //symmetric updates if not on true diagonal of global matrix
      if (!(msf.on_diagonal_() && (i == k))) {
        // updates to row sum of A - col comm
        local[rc_size+k] += aik;
        // updates to row sum of Z - col comm
        local[rc_size+ncol+k] += aikbibk;
      }
    }
  }
  CHKERRMPI(MPI_Comm_rank(msf.row_comm_,&rrank));
  CHKERRMPI(MPI_Comm_rank(msf.col_comm_,&crank));
  {
    auto rowdat = local.data();
    // reduce the row comms
    if (rsize) CHKERRMPI(MPI_Reduce(rrank ? rowdat : MPI_IN_PLACE,rowdat,rc_size,MPIU_SCALAR,MPI_SUM,0,msf.row_comm_));

    auto coldat = local.data()+rc_size;
    // reduce the col comms
    if (csize) CHKERRMPI(MPI_Reduce(crank ? coldat : MPI_IN_PLACE,coldat,ncol+ncol,MPIU_SCALAR,MPI_SUM,0,msf.col_comm_));
  }
  if (msf.on_diagonal_()) {
    // final reduction here
    CHKERRQ(VecRestoreArrayRead(vin,const_cast<const PetscScalar**>(&b_row)));
    CHKERRQ(VecRestoreArrayWrite(vout,&array_out));
  } else delete[] b_row;
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
  m->ops->getvecs   = get_vecs;

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
  PetscMPIInt rank;
  const auto  msf = impl_cast_(m);

  PetscFunctionBegin;
  CHKERRMPI(MPI_Comm_rank(PetscObjComm(m),&rank));
  CHKERRQ(PetscViewerASCIIPushSynchronized(vwr));
  CHKERRQ(PetscViewerASCIISynchronizedPrintf(vwr,"[%d] rows [%" PetscInt_FMT ":%" PetscInt_FMT "], cols [%" PetscInt_FMT ":%" PetscInt_FMT "]\n",rank,msf->rbegin_,msf->rend_,msf->cbegin_,msf->cend_));
  CHKERRQ(PetscViewerFlush(vwr));
  const auto& values = msf->data_;
  CHKERRQ(PetscScalarView(values.size(),values.data(),vwr));
  CHKERRQ(PetscViewerFlush(vwr));
  CHKERRQ(PetscViewerASCIIPopSynchronized(vwr));
  PetscFunctionReturn(0);
}

PetscErrorCode MatSymmFast::get_vecs(Mat m, Vec *right, Vec *left) noexcept
{
  const auto  msf     = impl_cast_(m);
  const auto  on_diag = msf->on_diagonal_();
  PetscMPIInt size;

  PetscFunctionBegin;
  CHKERRMPI(MPI_Comm_size(PetscObjComm(m),&size));
  if (right) {
    Vec v;

    CHKERRQ(VecCreate(PetscObjComm(m),&v));
    CHKERRQ(VecSetType(v,m->defaultvectype));
    if (on_diag) {
      CHKERRQ(VecSetSizes(v,msf->ncols_(0),m->cmap->N));
      v->ops->setrandom = [](Vec v, PetscRandom r) {
        const auto  n = v->map->n;
        PetscScalar *array;

        PetscFunctionBegin;
        CHKERRQ(VecGetArrayWrite(v,&array));
        for (auto i = 0; i < n; ++i) CHKERRQ(PetscRandomGetValue(r,array+i));
        CHKERRQ(VecRestoreArrayWrite(v,&array));
        PetscFunctionReturn(0);
      };
    } else {
      CHKERRQ(VecSetSizes(v,0,m->cmap->N));
      v->ops->setrandom = [](Vec,PetscRandom) { return 0; };
    }
    *right = v;
  }
  if (left) {
    Vec v;

    CHKERRQ(VecCreate(PetscObjComm(m),&v));
    CHKERRQ(VecSetType(v,m->defaultvectype));
    if (on_diag) {
      CHKERRQ(VecSetSizes(v,msf->nrows_(0),m->rmap->N));
      v->ops->setrandom = [](Vec v, PetscRandom r) {
        const auto  n = v->map->n;
        PetscScalar *array;

        PetscFunctionBegin;
        CHKERRQ(VecGetArrayWrite(v,&array));
        for (auto i = 0; i < n; ++i) CHKERRQ(PetscRandomGetValue(r,array+i));
        CHKERRQ(VecRestoreArrayWrite(v,&array));
        PetscFunctionReturn(0);
      };
    } else {
      CHKERRQ(VecSetSizes(v,0,m->rmap->N));
      v->ops->setrandom = [](Vec,PetscRandom) { return 0; };
    }
    *left = v;
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode MatMult_Private(Mat m, Vec x, Vec y) noexcept
{
  PetscFunctionBegin;
  CHKERRQ(VecLockReadPush(x));
  if (!m->ops->mult) SETERRQ1(PetscObjectComm(PetscObjectCast(m)),PETSC_ERR_SUP,"Matrix type %s does not have a multiply defined",(PetscObjectCast(m))->type_name);
  CHKERRQ(PetscLogEventBegin(MAT_Mult,m,x,y,0));
  CHKERRQ(m->ops->mult(m,x,y));
  CHKERRQ(PetscLogEventEnd(MAT_Mult,m,x,y,0));
  if (m->erroriffailure) CHKERRQ(VecValidValues(y,3,PETSC_FALSE));
  CHKERRQ(VecLockReadPop(x));
  PetscFunctionReturn(0);
}

static const auto default_lambda = [](Mat,Vec,Vec){ return 0; };

template <typename F = decltype(default_lambda)>
static PetscErrorCode MatMultTime(MPI_Comm comm, MatType type, PetscInt rows, PetscInt cols, F pre_process_fn = default_lambda, std::size_t its = 1000) noexcept
{
  PetscMPIInt size;
  Mat         mat;
  Vec         vin,vout;
  double      root_elapsed = 0;

  PetscFunctionBegin;
  CHKERRQ(MatCreate(comm,&mat));
  CHKERRQ(PetscObjectSetOptionsPrefix(PetscObjectCast(mat),type));
  CHKERRQ(MatSetSizes(mat,PETSC_DECIDE,PETSC_DECIDE,rows,cols));
  CHKERRQ(MatSetType(mat,type));
  CHKERRQ(MatSetOption(mat,MAT_SYMMETRIC,PETSC_TRUE));
  CHKERRQ(MatSetOption(mat,MAT_SYMMETRY_ETERNAL,PETSC_TRUE));
  CHKERRQ(MatSetFromOptions(mat));
  CHKERRQ(MatSetUp(mat));
  CHKERRQ(MatSetRandom(mat,nullptr));

  CHKERRQ(MatCreateVecs(mat,&vin,&vout));
  CHKERRQ(PetscObjectSetOptionsPrefix(PetscObjectCast(vin),type));
  CHKERRQ(PetscObjectSetOptionsPrefix(PetscObjectCast(vout),type));
  CHKERRQ(VecSetFromOptions(vin));
  CHKERRQ(VecSetFromOptions(vout));
  CHKERRQ(VecSetRandom(vin,nullptr));
  CHKERRQ(VecAssemblyBegin(vin));
  CHKERRQ(VecAssemblyEnd(vin));
  CHKERRQ(VecAssemblyBegin(vout));
  CHKERRQ(VecAssemblyEnd(vout));

  // warmup
  CHKERRQ(MatMult_Private(mat,vin,vout));
  CHKERRQ(PetscBarrier(PetscObjectCast(mat)));

  // for realsies
  const auto start_time = MPI_Wtime();
  for (auto i = 0; i < its; ++i) CHKERRQ(MatMult_Private(mat,vin,vout));
  const auto avg_time = (MPI_Wtime()-start_time)/static_cast<decltype(start_time)>(its);

  CHKERRMPI(MPI_Reduce(&avg_time,&root_elapsed,1,MPI_DOUBLE,MPI_SUM,0,comm));

  CHKERRMPI(MPI_Comm_size(comm,&size));
  CHKERRQ(PetscPrintf(comm,"%s %g\n",type,root_elapsed/size));

  // cleanup
  CHKERRQ(VecDestroy(&vin));
  CHKERRQ(VecDestroy(&vout));
  CHKERRQ(MatDestroy(&mat));
  PetscFunctionReturn(0);
}

int main(int argc, char*argv[])
{
  PetscInt       rows = 10, cols = 10;
  PetscErrorCode ierr;

  ierr = PetscInitialize(&argc,&argv,nullptr,nullptr);if (PetscUnlikely(ierr)) return ierr;
  CHKERRQ(MatRegister(MATSYMMFAST,MatSymmFast::create));
  auto comm = PETSC_COMM_WORLD;

  CHKERRQ(MatMultTime(comm,MATDENSE,rows,cols));
  CHKERRQ(MatMultTime(comm,MATSYMMFAST,rows,cols));

  ierr = PetscFinalize();
  return ierr;
}
