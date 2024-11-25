!***************************************************************************************************
! CublasModule.f90 10/18/17
! This file contains the module defining Fortran interfaces for the cuBLAS library
!***************************************************************************************************

module CublasModule
  !-------------------------------------------------------------------------------------------------
  ! Interface to cuBLAS routines
  !-------------------------------------------------------------------------------------------------
  use, intrinsic :: iso_c_binding

  type(c_ptr) :: cublas_handle
  !!$omp threadprivate(cublas_handle)

  enum, bind(c) !:: cublasStatus_t
    enumerator :: CUBLAS_STATUS_SUCCESS = 0
    enumerator :: CUBLAS_STATUS_NOT_INITIALIZED = 1
    enumerator :: CUBLAS_STATUS_ALLOC_FAILED = 3
    enumerator :: CUBLAS_STATUS_INVALID_VALUE = 7
    enumerator :: CUBLAS_STATUS_ARCH_MISMATCH = 8
    enumerator :: CUBLAS_STATUS_MAPPING_ERROR = 11
    enumerator :: CUBLAS_STATUS_EXECUTION_FAILED = 13
    enumerator :: CUBLAS_STATUS_INTERNAL_ERROR = 14
  end enum !cublasStatus_t

  enum, bind(c) !:: cublasFillMode_t
    enumerator :: CUBLAS_FILL_MODE_LOWER = 0
    enumerator :: CUBLAS_FILL_MODE_UPPER = 1
  end enum !cublasFillMode_t

  enum, bind(c) !:: cublasDiag type_t
    enumerator :: CUBLAS_DIAG_NON_UNIT = 0
    enumerator :: CUBLAS_DIAG_UNIT = 1
  end enum !cublasDiag    type_t

  enum, bind(c) !:: cublasSideMode_t
    enumerator :: CUBLAS_SIDE_LEFT = 0
    enumerator :: CUBLAS_SIDE_RIGHT = 1
  end enum !cublasSideMode_t

  enum, bind(c) !:: cublasOperation_t
    enumerator :: CUBLAS_OP_N = 0
    enumerator :: CUBLAS_OP_T = 1
    enumerator :: CUBLAS_OP_C = 2
  end enum !cublasOperation_t

  interface

    integer(c_int) function &
        & cublasInit() &
        & bind(c, name="cublasInit")
      use, intrinsic :: iso_c_binding
    end function cublasInit

    integer(c_int) function &
        & cublasShutdown() &
        & bind(c, name="cublasShutdown")
      use, intrinsic :: iso_c_binding
    end function cublasShutdown

    integer(c_int) function &
        & cublasCreate_v2(handle) &
        & bind(c, name="cublasCreate_v2")
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: handle
    end function cublasCreate_v2

    integer(c_int) function &
        & cublasDestroy_v2(handle) &
        & bind(c, name="cublasDestroy_v2")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: handle
    end function cublasDestroy_v2

    integer(c_int) function &
        & cublasGetStream_v2(handle, stream) &
        & bind(c, name="cublasGetStream_v2")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: handle
      type(c_ptr) :: stream
    end function cublasGetStream_v2

    integer(c_int) function &
        & cublasSetStream_v2(handle, stream) &
        & bind(c, name="cublasSetStream_v2")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: handle
      type(c_ptr), value :: stream
    end function cublasSetStream_v2

    integer(c_int) function &
        & cublasGetVector(n, elemSize, dx_src, incx, hy_dst, incy) &
        & bind(c, name="cublasGetVector")
      use, intrinsic :: iso_c_binding
      integer(c_int), value :: n
      integer(c_size_t), value :: elemSize
      type(c_ptr), value :: dx_src
      integer(c_int), value :: incx
      type(c_ptr), value :: hy_dst
      integer(c_int), value :: incy
    end function cublasGetVector

    integer(c_int) function &
        & cublasGetVectorAsync(n, elemSize, dx_src, incx, hy_dst, incy, stream) &
        & bind(c, name="cublasGetVectorAsync")
      use, intrinsic :: iso_c_binding
      integer(c_int), value :: n
      integer(c_size_t), value :: elemSize
      type(c_ptr), value :: dx_src
      integer(c_int), value :: incx
      type(c_ptr), value :: hy_dst
      integer(c_int), value :: incy
      type(c_ptr), value :: stream
    end function cublasGetVectorAsync

    integer(c_int) function &
        & cublasSetVector(n, elemSize, hx_src, incx, dy_dst, incy) &
        & bind(c, name="cublasSetVector")
      use, intrinsic :: iso_c_binding
      integer(c_int), value :: n
      integer(c_size_t), value :: elemSize
      type(c_ptr), value :: hx_src
      integer(c_int), value :: incx
      type(c_ptr), value :: dy_dst
      integer(c_int), value :: incy
    end function cublasSetVector

    integer(c_int) function &
        & cublasSetVectorAsync(n, elemSize, hx_src, incx, dy_dst, incy, stream) &
        & bind(c, name="cublasSetVectorAsync")
      use, intrinsic :: iso_c_binding
      integer(c_int), value :: n
      integer(c_size_t), value :: elemSize
      type(c_ptr), value :: hx_src
      integer(c_int), value :: incx
      type(c_ptr), value :: dy_dst
      integer(c_int), value :: incy
      type(c_ptr), value :: stream
    end function cublasSetVectorAsync

    integer(c_int) function &
        & cublasSetMatrix(rows, cols, elemSize, hA_src, lda, dB_dst, lddb) &
        & bind(c, name="cublasSetMatrix")
      use, intrinsic :: iso_c_binding
      integer(c_int), value :: rows
      integer(c_int), value :: cols
      integer(c_size_t), value :: elemSize
      type(c_ptr), value :: hA_src
      integer(c_int), value :: lda
      type(c_ptr), value :: dB_dst
      integer(c_int), value :: lddb
    end function cublasSetMatrix

    integer(c_int) function &
        & cublasSetMatrixAsync(rows, cols, elemSize, hA_src, lda, dB_dst, lddb, stream) &
        & bind(c, name="cublasSetMatrixAsync")
      use, intrinsic :: iso_c_binding
      integer(c_int), value :: rows
      integer(c_int), value :: cols
      integer(c_size_t), value :: elemSize
      type(c_ptr), value :: hA_src
      integer(c_int), value :: lda
      type(c_ptr), value :: dB_dst
      integer(c_int), value :: lddb
      type(c_ptr), value :: stream
    end function cublasSetMatrixAsync

    integer(c_int) function &
        & cublasSetBatchMatrixAsync(rows, cols, batch, elemSize, hA_src, lda, dB_dst, lddb, stream) &
        & bind(c, name="cublasSetBatchMatrixAsync")
      use, intrinsic :: iso_c_binding
      integer(c_int), value :: rows
      integer(c_int), value :: cols
      integer(c_int), value :: batch
      integer(c_size_t), value :: elemSize
      type(c_ptr), value :: hA_src
      integer(c_int), value :: lda
      type(c_ptr), value :: dB_dst
      integer(c_int), value :: lddb
      type(c_ptr), value :: stream
    end function cublasSetBatchMatrixAsync

    integer(c_int) function &
        & cublasGetMatrix(rows, cols, elemSize, dA_src, ldda, hB_dst, ldb) &
        & bind(c, name="cublasGetMatrix")
      use, intrinsic :: iso_c_binding
      integer(c_int), value :: rows
      integer(c_int), value :: cols
      integer(c_size_t), value :: elemSize
      type(c_ptr), value :: dA_src
      integer(c_int), value :: ldda
      type(c_ptr), value :: hB_dst
      integer(c_int), value :: ldb
    end function cublasGetMatrix

    integer(c_int) function &
        & cublasGetMatrixAsync(rows, cols, elemSize, dA_src, ldda, hB_dst, ldb, stream) &
        & bind(c, name="cublasGetMatrixAsync")
      use, intrinsic :: iso_c_binding
      integer(c_int), value :: rows
      integer(c_int), value :: cols
      integer(c_size_t), value :: elemSize
      type(c_ptr), value :: dA_src
      integer(c_int), value :: ldda
      type(c_ptr), value :: hB_dst
      integer(c_int), value :: ldb
      type(c_ptr), value :: stream
    end function cublasGetMatrixAsync

    integer(c_int) function &
        & cublasDnrm2_v2(handle, n, dx, incx, xnorm) &
        & bind(c, name="cublasDnrm2_v2")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: handle
      integer(c_int), value :: n
      type(c_ptr), value :: dx
      integer(c_int), value :: incx
      real(c_double) :: xnorm
    end function cublasDnrm2_v2

    integer(c_int) function &
        & cublasDaxpy_v2(handle, n, alpha, dx, incx, dy, incy) &
        & bind(c, name="cublasDaxpy_v2")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: handle
      integer(c_int), value :: n
      real(c_double), value :: alpha
      type(c_ptr), value :: dx
      integer(c_int), value :: incx
      type(c_ptr), value :: dy
      integer(c_int), value :: incy
    end function cublasDaxpy_v2

    integer(c_int) function &
        & cublasDgemv_v2(handle, trans, m, n, alpha, dA, ldda, dx, incx, beta, dy, incy) &
        & bind(c, name="cublasDgemv_v2")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: handle
      integer(c_int), value :: trans
      integer(c_int), value :: m
      integer(c_int), value :: n
      real(c_double) :: alpha
      type(c_ptr), value :: dA
      integer(c_int), value :: ldda
      type(c_ptr), value :: dx
      integer(c_int), value :: incx
      real(c_double) :: beta
      type(c_ptr), value :: dy
      integer(c_int), value :: incy
    end function cublasDgemv_v2

    integer(c_int) function &
        & cublasDgetrfBatched(handle, n, dA, ldda, dP, dInfo, nbatch) &
        & bind(c, name="cublasDgetrfBatched")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: handle
      integer(c_int), value :: n
      type(c_ptr), value :: dA
      integer(c_int), value :: ldda
      type(c_ptr), value :: dP
      type(c_ptr), value :: dInfo
      integer(c_int), value :: nbatch
    end function cublasDgetrfBatched

    integer(c_int) function &
        & cublasDgetriBatched(handle, n, dA, ldda, dP, dC, lddc, dInfo, nbatch) &
        & bind(c, name="cublasDgetriBatched")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: handle
      integer(c_int), value :: n
      type(c_ptr), value :: dA
      integer(c_int), value :: ldda
      type(c_ptr), value :: dP
      type(c_ptr), value :: dC
      integer(c_int), value :: lddc
      type(c_ptr), value :: dInfo
      integer(c_int), value :: nbatch
    end function cublasDgetriBatched

    integer(c_int) function &
        & cublasDtrsmBatched(handle, side, uplo, trans, diag, m, n, alpha, dA, ldda, dB, lddb, nbatch) &
        & bind(c, name="cublasDtrsmBatched")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: handle
      integer(c_int), value :: side
      integer(c_int), value :: uplo
      integer(c_int), value :: trans
      integer(c_int), value :: diag
      integer(c_int), value :: m
      integer(c_int), value :: n
      real(c_double) :: alpha
      type(c_ptr), value :: dA
      integer(c_int), value :: ldda
      type(c_ptr), value :: dB
      integer(c_int), value :: lddb
      integer(c_int), value :: nbatch
    end function cublasDtrsmBatched

    integer(c_int) function &
        & cublasDgemmBatched(handle, transa, transb, m, n, k, alpha, dA, ldda, dB, lddb, beta, dC, lddc, nbatch) &
        & bind(c, name="cublasDgemmBatched")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: handle
      integer(c_int), value :: transa
      integer(c_int), value :: transb
      integer(c_int), value :: m
      integer(c_int), value :: n
      integer(c_int), value :: k
      real(c_double) :: alpha
      type(c_ptr), value :: dA
      integer(c_int), value :: ldda
      type(c_ptr), value :: dB
      integer(c_int), value :: lddb
      real(c_double) :: beta
      type(c_ptr), value :: dC
      integer(c_int), value :: lddc
      integer(c_int), value :: nbatch
    end function cublasDgemmBatched

    integer(c_int) function &
        & cublasDgemmStridedBatched(handle, transa, transb, m, n, k, alpha, &
        & dA, ldda, strideA, dB, lddb, strideB, beta, dC, lddc, strideC, nbatch) &
        & bind(c, name="cublasDgemmStridedBatched")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: handle
      integer(c_int), value :: transa
      integer(c_int), value :: transb
      integer(c_int), value :: m
      integer(c_int), value :: n
      integer(c_int), value :: k
      real(c_double) :: alpha
      type(c_ptr), value :: dA
      integer(c_int), value :: ldda
      integer(c_int), value :: strideA
      type(c_ptr), value :: dB
      integer(c_int), value :: lddb
      integer(c_int), value :: strideB
      real(c_double) :: beta
      type(c_ptr), value :: dC
      integer(c_int), value :: strideC
      integer(c_int), value :: lddc
      integer(c_int), value :: nbatch
    end function cublasDgemmStridedBatched

    integer(c_int) function &
        & cublasDtrsv(uplo, trans, diag, n, dA, ldda, dx, incx) &
        & bind(c, name="cublasDtrsv")
      use, intrinsic :: iso_c_binding
      character(c_char), value :: uplo
      character(c_char), value :: trans
      character(c_char), value :: diag
      integer(c_int), value :: n
      type(c_ptr), value :: dA
      integer(c_int), value :: ldda
      type(c_ptr), value :: dx
      integer(c_int), value :: incx
    end function cublasDtrsv

    integer(c_int) function &
        & cublasDtrsv_v2(handle, uplo, trans, diag, n, dA, ldda, dx, incx) &
        & bind(c, name="cublasDtrsv_v2")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: handle
      integer(c_int), value :: uplo
      integer(c_int), value :: trans
      integer(c_int), value :: diag
      integer(c_int), value :: n
      type(c_ptr), value :: dA
      integer(c_int), value :: ldda
      type(c_ptr), value :: dx
      integer(c_int), value :: incx
    end function cublasDtrsv_v2

    integer(c_int) function &
        & cublasDtrsm(side, uplo, trans, diag, m, n, alpha, dA, ldda, dB, lddb) &
        & bind(c, name="cublasDtrsm")
      use, intrinsic :: iso_c_binding
      character(c_char), value :: side
      character(c_char), value :: uplo
      character(c_char), value :: trans
      character(c_char), value :: diag
      integer(c_int), value :: m
      integer(c_int), value :: n
      real(c_double), value :: alpha
      type(c_ptr), value :: dA
      integer(c_int), value :: ldda
      type(c_ptr), value :: dB
      integer(c_int), value :: lddb
    end function cublasDtrsm

    integer(c_int) function &
        & cublasDtrsm_v2(handle, side, uplo, trans, diag, m, n, alpha, dA, ldda, dB, lddb) &
        & bind(c, name="cublasDtrsm_v2")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: handle
      integer(c_int), value :: side
      integer(c_int), value :: uplo
      integer(c_int), value :: trans
      integer(c_int), value :: diag
      integer(c_int), value :: m
      integer(c_int), value :: n
      real(c_double), value :: alpha
      type(c_ptr), value :: dA
      integer(c_int), value :: ldda
      type(c_ptr), value :: dB
      integer(c_int), value :: lddb
    end function cublasDtrsm_v2

    integer(c_int) function &
        & cublasDgemm(transa, transb, m, n, k, alpha, dA, ldda, dB, lddb, beta, dC, lddc) &
        & bind(c, name="cublasDgemm")
      use, intrinsic :: iso_c_binding
      character(c_char), value :: transa
      character(c_char), value :: transb
      integer(c_int), value :: m
      integer(c_int), value :: n
      integer(c_int), value :: k
      real(c_double), value :: alpha
      type(c_ptr), value :: dA
      integer(c_int), value :: ldda
      type(c_ptr), value :: dB
      integer(c_int), value :: lddb
      real(c_double), value :: beta
      type(c_ptr), value :: dC
      integer(c_int), value :: lddc
    end function cublasDgemm

    integer(c_int) function &
        & cublasDgemm_v2(handle, transa, transb, m, n, k, alpha, dA, ldda, dB, lddb, beta, dC, lddc) &
        & bind(c, name="cublasDgemm_v2")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: handle
      integer(c_int), value :: transa
      integer(c_int), value :: transb
      integer(c_int), value :: m
      integer(c_int), value :: n
      integer(c_int), value :: k
      real(c_double) :: alpha
      type(c_ptr), value :: dA
      integer(c_int), value :: ldda
      type(c_ptr), value :: dB
      integer(c_int), value :: lddb
      real(c_double) :: beta
      type(c_ptr), value :: dC
      integer(c_int), value :: lddc
    end function cublasDgemm_v2

    integer(c_int) function &
        & cublasDgetrsBatched(handle, trans, n, nrhs, dA, ldda, dP, dB, lddb, hInfo, nbatch) &
        & bind(c, name="cublasDgetrsBatched")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: handle
      integer(c_int), value :: trans
      integer(c_int), value :: n
      integer(c_int), value :: nrhs
      type(c_ptr), value :: dA
      integer(c_int), value :: ldda
      type(c_ptr), value :: dP
      type(c_ptr), value :: dB
      integer(c_int), value :: lddb
      type(c_ptr), value :: hInfo
      integer(c_int), value :: nbatch
    end function cublasDgetrsBatched

    integer(c_int) function &
        & cublasDgeam(handle, transa, transb, m, n, alpha, dA, ldda, beta, dB, lddb, dC, lddc) &
        & bind(c, name="cublasDgeam")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: handle
      integer(c_int), value :: transa
      integer(c_int), value :: transb
      integer(c_int), value :: m
      integer(c_int), value :: n
      real(c_double) :: alpha
      type(c_ptr), value :: dA
      integer(c_int), value :: ldda
      real(c_double) :: beta
      type(c_ptr), value :: dB
      integer(c_int), value :: lddb
      type(c_ptr), value :: dC
      integer(c_int), value :: lddc
    end function cublasDgeam

    integer(c_int) function &
        & cublasDdgmm(handle, mode, m, n, dA, ldda, dx, incx, dC, lddc) &
        & bind(c, name="cublasDdgmm")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: handle
      integer(c_int), value :: mode
      integer(c_int), value :: m
      integer(c_int), value :: n
      type(c_ptr), value :: dA
      integer(c_int), value :: ldda
      type(c_ptr), value :: dx
      integer(c_int), value :: incx
      type(c_ptr), value :: dC
      integer(c_int), value :: lddc
    end function cublasDdgmm

  end interface

end module CublasModule
