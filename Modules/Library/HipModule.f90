!***************************************************************************************************
! HipModule.f90 10/18/17
! This file contains the module defining Fortran interfaces for the HIP Runtime API
!***************************************************************************************************

module HipModule
  !-------------------------------------------------------------------------------------------------
  ! Interface to hip Runtime API
  !-------------------------------------------------------------------------------------------------
  use, intrinsic :: iso_c_binding
  use hipfort, only: &
    hipGetDeviceCount, &
    hipSetDevice, &
    hipStreamCreate, &
    hipStreamSynchronize

  integer :: mystream
  type(c_ptr) :: stream, event

  integer :: nstream, nevent
  type(c_ptr), allocatable :: streamArray(:)
  type(c_ptr), allocatable :: eventArray(:)

end module HipModule
