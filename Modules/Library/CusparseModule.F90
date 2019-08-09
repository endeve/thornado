!***************************************************************************************************
! CusparseModule.f90 10/18/17
! this file contains the module defining fortran interfaces for the cublas library
!***************************************************************************************************

module CusparseModule
  !-------------------------------------------------------------------------------------------------
  ! Interface to cuSPARSE routines
  !-------------------------------------------------------------------------------------------------
  use, intrinsic :: iso_c_binding

  type(c_ptr) :: cusparse_handle
  !!$omp threadprivate(cusparse_handle)

  enum, bind(c) !:: cusparseIndexBase_t
    enumerator :: CUSPARSE_INDEX_BASE_ZERO = 0
    enumerator :: CUSPARSE_INDEX_BASE_ONE = 1
  end enum !cusparseIndexBase_t

  enum, bind(c) !:: cusparseDiagType_t
    enumerator :: CUSPARSE_DIAG_TYPE_NON_UNIT = 0
    enumerator :: CUSPARSE_DIAG_TYPE_UNIT = 1
  end enum !cusparseDiagType_t

  enum, bind(c) !:: cusparseFillMode_t
    enumerator :: CUSPARSE_FILL_MODE_LOWER = 0
    enumerator :: CUSPARSE_FILL_MODE_UPPER = 1
  end enum !cusparseFillMode_t

  enum, bind(c) !:: cusparseOperation_t
    enumerator :: CUSPARSE_OPERATION_NON_TRANSPOSE = 0
    enumerator :: CUSPARSE_OPERATION_TRANSPOSE = 1
    enumerator :: CUSPARSE_OERATIONP_CONJUGATE_TRANSPOSE = 2
  end enum !cusparseOperation_t

  enum, bind(c) !:: cusparseStatus_t
    enumerator :: CUSPARSE_STATUS_SUCCESS = 0
    enumerator :: CUSPARSE_STATUS_NOT_INITIALIZED  = 1
    enumerator :: CUSPARSE_STATUS_ALLOC_FAILED = 3
    enumerator :: CUSPARSE_STATUS_INVALID_VALUE = 7
    enumerator :: CUSPARSE_STATUS_ARCH_MISMATCH = 8
    enumerator :: CUSPARSE_STATUS_MAPPING_ERROR = 11
    enumerator :: CUSPARSE_STATUS_EXECUTION_FAILED = 13
    enumerator :: CUSPARSE_STATUS_INTERNAL_ERROR = 14
    enumerator :: CUSPARSE_MATRIX_TYPE_NOT_SUPPORTED = 15
  end enum !cusparseStatus_t

  interface

    function cusparseCreate(handle) &
        & bind(c, name="cusparseCreate")
      use, intrinsic :: iso_c_binding
      integer(c_int) :: cusparseCreate
      type(c_ptr) :: handle
    end function cusparseCreate

    function cusparseCreateSolveAnalysisInfo(info) &
        & bind(c, name="cusparseCreateSolveAnalysisInfo")
      use, intrinsic :: iso_c_binding
      integer(c_int) :: cusparseCreateSolveAnalysisInfo
      type(c_ptr) :: info
    end function cusparseCreateSolveAnalysisInfo

    function cusparseDestroy(handle) &
        & bind(c, name="cusparseDestroy")
      use, intrinsic :: iso_c_binding
      integer(c_int) :: cusparseDestroy
      type(c_ptr), value :: handle
    end function cusparseDestroy

    function cusparseDestroySolveAnalysisInfo(info) &
        & bind(c, name="cusparseDestroySolveAnalysisInfo")
      use, intrinsic :: iso_c_binding
      integer(c_int) :: cusparseDestroySolveAnalysisInfo
      type(c_ptr), value :: info
    end function cusparseDestroySolveAnalysisInfo

    function cusparseSetStream(handle, stream) &
        & bind(c, name="cusparseSetStream")
      use, intrinsic :: iso_c_binding
      integer(c_int) :: cusparseSetStream
      type(c_ptr), value :: handle
      type(c_ptr), value :: stream
    end function cusparseSetStream

    function cusparseDgtsvStridedBatch(handle, m, dl, d, du, x, batchcount, batchstride) &
        & bind(c, name="cusparseDgtsvStridedBatch")
      use, intrinsic :: iso_c_binding
      integer(c_int) :: cusparseDgtsvStridedBatch
      type(c_ptr), value :: handle
      integer(c_int), value :: m
      type(c_ptr), value :: dl
      type(c_ptr), value :: d
      type(c_ptr), value :: du
      type(c_ptr), value :: x
      integer(c_int), value :: batchcount
      integer(c_int), value :: batchstride
    end function cusparseDgtsvStridedBatch

    function cusparseDgthr(handle, nnz, y, xval, xind, idxbase ) &
        & bind(c, name="cusparseDgthr")
      use, intrinsic :: iso_c_binding
      integer(c_int) :: cusparseDgthr
      type(c_ptr), value :: handle
      integer(c_int), value :: nnz
      type(c_ptr), value :: y
      type(c_ptr), value :: xval
      type(c_ptr), value :: xind
      integer(c_int), value :: idxbase
    end function cusparseDgthr

  end interface

end module CusparseModule
