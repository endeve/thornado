!***************************************************************************************************
! OpenACCModule.f90 10/18/17
! This file contains the module defining Fortran interfaces for the OpenACC Runtime Library
!***************************************************************************************************

module OpenACCModule
  !-------------------------------------------------------------------------------------------------
  ! Interface to OpenACC routines
  !-------------------------------------------------------------------------------------------------
  use, intrinsic :: iso_c_binding
  use openacc, only: acc_init, acc_set_device_num
  implicit none

  enum, bind(c) !:: acc_device_t
    enumerator :: acc_device_none = 0
    enumerator :: acc_device_default = 1
    enumerator :: acc_device_host = 2
    enumerator :: acc_device_not_host = 3
    enumerator :: acc_device_nvidia = 4
    enumerator :: acc_device_radeon = 5
    enumerator :: acc_device_xeonphi = 6
    enumerator :: acc_device_pgi_opencl = 7
    enumerator :: acc_device_nvidia_opencl = 8
    enumerator :: acc_device_opencl = 9
    enumerator :: acc_device_current = 10
  end enum

  interface

    subroutine acc_map_data(hostptr,devptr,bytes) &
        bind(c,name="acc_map_data")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: hostptr
      type(c_ptr), value :: devptr
      integer(c_size_t), value :: bytes
    end subroutine acc_map_data

    subroutine acc_unmap_data(hostptr) &
        bind(c,name="acc_unmap_data")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: hostptr
    end subroutine acc_unmap_data

    type(c_ptr) function acc_deviceptr(hostptr) &
        bind(c,name="acc_deviceptr")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: hostptr
    end function acc_deviceptr

    type(c_ptr) function acc_hostptr(devptr) &
        bind(c,name="acc_hostptr")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: devptr
    end function acc_hostptr

    integer(c_int) function acc_is_present(hostptr,bytes) &
        bind(c,name="acc_is_present" )
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: hostptr
      integer(c_size_t), value :: bytes
    end function acc_is_present

  end interface

end module OpenACCModule
