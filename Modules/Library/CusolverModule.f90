!***************************************************************************************************
! CusolverModule.f90 10/18/17
! This file contains the module defining Fortran interfaces for the cuSOLVER library
!***************************************************************************************************

module CusolverModule
  !-------------------------------------------------------------------------------------------------
  ! Interface to cuSOLVER routines
  !-------------------------------------------------------------------------------------------------
  use, intrinsic :: iso_c_binding

  type(c_ptr) :: cusolver_handle
  !!$omp threadprivate(cusolver_handle)

  enum, bind(c) !:: cusolverStatus_t
    enumerator :: CUSOLVER_STATUS_SUCCESS = 0
    enumerator :: CUSOLVER_STATUS_NOT_INITIALIZED = 1
    enumerator :: CUSOLVER_STATUS_ALLOC_FAILED = 2
    enumerator :: CUSOLVER_STATUS_INVALID_VALUE = 3
    enumerator :: CUSOLVER_STATUS_ARCH_MISMATCH = 4
    enumerator :: CUSOLVER_STATUS_MAPPING_ERROR = 5
    enumerator :: CUSOLVER_STATUS_EXECUTION_FAILED = 6
    enumerator :: CUSOLVER_STATUS_INTERNAL_ERROR = 7
    enumerator :: CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED = 8
    enumerator :: CUSOLVER_STATUS_NOT_SUPPORTED  = 9
    enumerator :: CUSOLVER_STATUS_ZERO_PIVOT = 10
    enumerator :: CUSOLVER_STATUS_INVALID_LICENSE = 11
  end enum !cusolverStatus_t

  interface

    integer(c_int) function &
        & cusolverDnCreate(handle) &
        & bind(c, name="cusolverDnCreate")
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: handle
    end function cusolverDnCreate

    integer(c_int) function &
        & cusolverDnDestroy(handle) &
        & bind(c, name="cusolverDnDestroy")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: handle
    end function cusolverDnDestroy

    integer(c_int) function &
        & cusolverDnGetStream(handle, stream) &
        & bind(c, name="cusolverDnGetStream")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: handle
      type(c_ptr) :: stream
    end function cusolverDnGetStream

    integer(c_int) function &
        & cusolverDnSetStream(handle, stream) &
        & bind(c, name="cusolverDnSetStream")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: handle
      type(c_ptr), value :: stream
    end function cusolverDnSetStream

    integer(c_int) function &
        & cusolverDnDgeqrf_bufferSize(handle, m, n, A, lda, Lwork) &
        & bind(c, name="cusolverDnDgeqrf_bufferSize")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: handle
      integer(c_int), value :: m
      integer(c_int), value :: n
      type(c_ptr), value :: A
      integer(c_int), value :: lda
      integer(c_int), target :: Lwork
    end function cusolverDnDgeqrf_bufferSize

    integer(c_int) function &
        & cusolverDnDgeqrf(handle, m, n, A, lda, TAU, Workspace, Lwork, devInfo ) &
        & bind(c, name="cusolverDnDgeqrf")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: handle
      integer(c_int), value :: m
      integer(c_int), value :: n
      type(c_ptr), value :: A
      integer(c_int), value :: lda
      type(c_ptr), value :: TAU
      type(c_ptr), value :: Workspace
      integer(c_int), value :: Lwork
      type(c_ptr), value :: devInfo
    end function cusolverDnDgeqrf

    integer(c_int) function &
        & cusolverDnDormqr(handle, side, trans, m, n, k, A, lda, tau, C, ldc, work, lwork, devInfo ) &
        & bind(c, name="cusolverDnDormqr")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: handle
      integer(c_int), value :: side
      integer(c_int), value :: trans
      integer(c_int), value :: m
      integer(c_int), value :: n
      integer(c_int), value :: k
      type(c_ptr), value :: A
      integer(c_int), value :: lda
      type(c_ptr), value :: tau
      type(c_ptr), value :: C
      integer(c_int), value :: ldc
      type(c_ptr), value :: work
      integer(c_int), value :: lwork
      type(c_ptr), value :: devInfo
    end function cusolverDnDormqr

  end interface

end module CusolverModule
