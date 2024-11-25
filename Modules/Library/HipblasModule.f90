!***************************************************************************************************
! HipblasModule.f90 10/18/17
! This file contains the module defining Fortran interfaces for the hipBLAS library
!***************************************************************************************************

module HipblasModule
  !-------------------------------------------------------------------------------------------------
  ! Interface to hipBLAS routines
  !-------------------------------------------------------------------------------------------------
  use, intrinsic :: iso_c_binding
  use hipfort_hipblas_enums, only: &
    HIPBLAS_OP_N, HIPBLAS_OP_T, &
    HIPBLAS_SIDE_LEFT, &
    HIPBLAS_FILL_MODE_UPPER, &
    HIPBLAS_DIAG_NON_UNIT
  use hipfort_hipblas, only: &
    hipblasCreate, &
    hipblasGetStream, &
    hipblasSetStream, &
    hipblasDnrm2, &
    hipblasDaxpy, &
    hipblasDgemm, &
    hipblasDgemmStridedBatched, &
    hipblasDgetrfBatched, &
    hipblasDgetrsBatched, &
    hipblasDgemv, &
    hipblasDtrsv, &
    hipblasDtrsm, &
    hipblasDgeam, &
    hipblasDdgmm

  type(c_ptr) :: hipblas_handle

end module HipblasModule
