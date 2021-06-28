!***************************************************************************************************
! RocsolverModule.f90 10/18/17
! This file contains the module defining Fortran interfaces for the rocSOLVER library
!***************************************************************************************************

module RocsolverModule
  !-------------------------------------------------------------------------------------------------
  ! Interface to rocSOLVER routines
  !-------------------------------------------------------------------------------------------------
  use, intrinsic :: iso_c_binding
  use hipfort_rocsolver, only: &
    rocsolver_dgeqrf, &
    rocsolver_dormqr, &
    rocsolver_dgetrf_batched, &
    rocsolver_dgetrs_batched

  type(c_ptr) :: rocsolver_handle
  !!$omp threadprivate(rocsolver_handle)

end module RocsolverModule
