!***************************************************************************************************
! RocsparseModule.f90 10/18/17
! this file contains the module defining fortran interfaces for the rocSPARSE library
!***************************************************************************************************

module RocsparseModule
  !-------------------------------------------------------------------------------------------------
  ! Interface to rocSPARSE routines
  !-------------------------------------------------------------------------------------------------
  use, intrinsic :: iso_c_binding
  use hipfort_rocsparse_enums, only: &
    rocsparse_index_base_one
  use hipfort_rocsparse, only: &
    rocsparse_dgthr

  type(c_ptr) :: rocsparse_handle
  !!$omp threadprivate(rocsparse_handle)

end module RocsparseModule
