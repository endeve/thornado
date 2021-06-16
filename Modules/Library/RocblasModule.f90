!***************************************************************************************************
! RocblasModule.f90 10/18/17
! This file contains the module defining Fortran interfaces for the rocBLAS library
!***************************************************************************************************

module RocblasModule
  !-------------------------------------------------------------------------------------------------
  ! Interface to rocBLAS routines
  !-------------------------------------------------------------------------------------------------
  use, intrinsic :: iso_c_binding
  use hipfort_rocblas_enums, only: &
    rocblas_operation_none, &
    rocblas_operation_transpose, &
    rocblas_side_left, &
    rocblas_fill_upper, &
    rocblas_diagonal_non_unit
  use hipfort_rocblas, only: &
    rocblas_create_handle, &
    rocblas_get_stream, &
    rocblas_set_stream, &
    rocblas_dnrm2, &
    rocblas_daxpy, &
    rocblas_dgemm, &
    rocblas_dgemm_strided_batched, &
    rocblas_dgetrf_batched, &
    rocblas_dgetrs_batched, &
    rocblas_dgemv, &
    rocblas_dtrsv, &
    rocblas_dtrsm, &
    rocblas_dgeam, &
    rocblas_ddgmm

  type(c_ptr) :: rocblas_handle

end module RocblasModule
