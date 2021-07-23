!***************************************************************************************************
! HipsparseModule.f90 10/18/17
! this file contains the module defining fortran interfaces for the hipSPARSE library
!***************************************************************************************************

module HipsparseModule
  !-------------------------------------------------------------------------------------------------
  ! Interface to hipSPARSE routines
  !-------------------------------------------------------------------------------------------------
  use, intrinsic :: iso_c_binding
  use hipfort_hipsparse_enums, only: &
    HIPSPARSE_INDEX_BASE_ONE
  use hipfort_hipsparse, only: &
    hipsparseCreate, &
    hipsparseSetStream, &
    hipsparseDgthr

  type(c_ptr) :: hipsparse_handle
  !!$omp threadprivate(hipsparse_handle)

end module HipsparseModule
