!***************************************************************************************************
! OpenMPModule.f90 10/18/17
! This file contains the module defining Fortran interfaces for the OpenACC Runtime Library
!***************************************************************************************************

module OpenMPModule
  !-------------------------------------------------------------------------------------------------
  ! Interface to OpenACC routines
  !-------------------------------------------------------------------------------------------------
  use, intrinsic :: iso_c_binding
  use omp_lib, only: &
    omp_set_default_device, &
    omp_get_default_device, &
    omp_is_initial_device
  implicit none

  interface

    integer(c_int) function omp_target_is_present(hostptr,device) &
        bind(c,name="omp_target_is_present" )
      use, intrinsic :: iso_c_binding
      type(c_ptr),    value :: hostptr
      integer(c_int), value :: device
    end function omp_target_is_present

  end interface

end module OpenMPModule
