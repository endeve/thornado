!***************************************************************************************************
! HipModule.f90 10/18/17
! This file contains the module defining Fortran interfaces for the HIP Runtime API
!***************************************************************************************************

module HipModule
  !-------------------------------------------------------------------------------------------------
  ! Interface to hip Runtime API
  !-------------------------------------------------------------------------------------------------
  use, intrinsic :: iso_c_binding
  use hipfort_check, only: &
    hipCheck, &
    hipblasCheck, &
    hipsparseCheck, &
    rocblasCheck, &
    rocsparseCheck, &
    rocsolverCheck
  use hipfort, only: &
    hipGetDeviceCount, &
    hipSetDevice, &
    hipStreamCreate, &
    hipStreamSynchronize, &
    hipDeviceSynchronize

  integer :: mystream
  type(c_ptr), target :: stream
  type(c_ptr) :: event

  integer :: nstream, nevent
  type(c_ptr), allocatable :: streamArray(:)
  type(c_ptr), allocatable :: eventArray(:)

end module HipModule
