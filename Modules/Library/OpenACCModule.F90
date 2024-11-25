!***************************************************************************************************
! OpenACCModule.f90 10/18/17
! This file contains the module defining Fortran interfaces for the OpenACC Runtime Library
!***************************************************************************************************

module OpenACCModule
  !-------------------------------------------------------------------------------------------------
  ! Interface to OpenACC routines
  !-------------------------------------------------------------------------------------------------
  use, intrinsic :: iso_c_binding
  use openacc, only: acc_init, acc_set_device_num, acc_get_device_num, &
    acc_device_nvidia, acc_device_host, acc_device_default, acc_async_sync, acc_async_noval
  implicit none

  !integer(c_int) :: acc_async_default
  !integer(c_int) :: acc_queue
  !!$omp threadprivate(acc_queue,acc_async_default)

  !enum, bind(c) !:: acc_device_t
  !  enumerator :: acc_device_none = 0
  !  enumerator :: acc_device_default = 1
  !  enumerator :: acc_device_host = 2
  !  enumerator :: acc_device_not_host = 3
  !  enumerator :: acc_device_nvidia = 4
  !  enumerator :: acc_device_radeon = 5
  !  enumerator :: acc_device_xeonphi = 6
  !  enumerator :: acc_device_pgi_opencl = 7
  !  enumerator :: acc_device_nvidia_opencl = 8
  !  enumerator :: acc_device_opencl = 9
  !  enumerator :: acc_device_current = 10
  !end enum

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

    subroutine acc_set_cuda_stream(async,stream) &
        bind(c,name="acc_set_cuda_stream" )
      use, intrinsic :: iso_c_binding
      integer(c_long_long), value :: async
      type(c_ptr), value :: stream
    end subroutine acc_set_cuda_stream

    type(c_ptr) function acc_get_cuda_stream(async) &
        bind(c,name="acc_get_cuda_stream")
      use, intrinsic :: iso_c_binding
      integer(c_long_long), value :: async
    end function acc_get_cuda_stream

    type(c_ptr) function acc_get_current_cuda_device() &
        bind(c,name="acc_get_current_cuda_device")
      use, intrinsic :: iso_c_binding
    end function acc_get_current_cuda_device

    type(c_ptr) function acc_get_current_cuda_context() &
        bind(c,name="acc_get_current_cuda_context")
      use, intrinsic :: iso_c_binding
    end function acc_get_current_cuda_context

    integer(c_int) function acc_get_default_async() &
        bind(c,name="acc_get_default_async" )
      use, intrinsic :: iso_c_binding
    end function acc_get_default_async

    subroutine acc_set_default_async(async) &
        bind(c,name="acc_set_default_async" )
      use, intrinsic :: iso_c_binding
      integer(c_long_long), value :: async
    end subroutine acc_set_default_async

    integer(c_int) function acc_on_device_i(devicetype) &
        bind(c,name="acc_on_device" )
#if defined(THORNADO_OMP_OL)
      !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
      !$ACC ROUTINE SEQ
#endif
      use, intrinsic :: iso_c_binding
      integer(c_int), value :: devicetype
    end function acc_on_device_i

  end interface

contains

  logical function acc_on_device(devicetype)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    use, intrinsic :: iso_c_binding
    integer(kind(acc_device_host)) :: devicetype
    acc_on_device = ( .not. acc_on_device_i(int(devicetype,c_int)) == 0 )
    return
  end function acc_on_device

end module OpenACCModule
