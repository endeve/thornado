!***************************************************************************************************
! HipModule.f90 10/18/17
! This file contains the module defining Fortran interfaces for the HIP Runtime API
!***************************************************************************************************

module HipModule
  !-------------------------------------------------------------------------------------------------
  ! Interface to hip Runtime API
  !-------------------------------------------------------------------------------------------------
  use, intrinsic :: iso_c_binding

  integer :: mystream
  type(c_ptr) :: stream, event

  integer :: nstream, nevent
  type(c_ptr), allocatable :: streamArray(:)
  type(c_ptr), allocatable :: eventArray(:)

  enum, bind(c) !:: hipMemcpyKind
    enumerator :: hipMemcpyHostToHost = 0
    enumerator :: hipMemcpyHostToDevice = 1
    enumerator :: hipMemcpyDeviceToHost = 2
    enumerator :: hipMemcpyDeviceToDevice = 3
    enumerator :: hipMemcpyDefault = 4
  end enum !hipMemcpyKind

  !enum, bind(c) !:: hipMallocManaged
  !  enumerator :: hipMemAttachGlobal = int(z'01')
  !  enumerator :: hipMemAttachHost = int(z'02')
  !  enumerator :: hipMemAttachSingle = int(z'04')
  !end enum !hipMallocManaged

  enum, bind(c) !:: hipHostMallocFlags
    enumerator :: hipHostMallocDefault = int(z'00')
    enumerator :: hipHostMallocPortable = int(z'01')
    enumerator :: hipHostMallocMapped = int(z'02')
    enumerator :: hipHostMallocWriteCombined = int(z'04')
  end enum !hipHostMallocFlags

  enum, bind(c) !:: hipHostRegisterFlags
    enumerator :: hipHostRegisterDefault = int(z'00')
    enumerator :: hipHostRegisterPortable = int(z'01')
    enumerator :: hipHostRegisterMapped = int(z'02')
    enumerator :: hipHostRegisterIoMemory = int(z'04')
  end enum !hipHostRegisterFlags

  enum, bind(c) !:: hipDeviceFlags
    enumerator :: hipDeviceScheduleAuto = int(z'00')
    enumerator :: hipDeviceScheduleSpin = int(z'01')
    enumerator :: hipDeviceScheduleYield = int(z'02')
    enumerator :: hipDeviceScheduleBlockingSync = int(z'04')
    enumerator :: hipDeviceBlockingSync = int(z'04')
    enumerator :: hipDeviceScheduleMask = int(z'07')
    enumerator :: hipDeviceMapHost = int(z'08')
    enumerator :: hipDeviceLmemResizeToMax = int(z'16')
    !enumerator :: hipDeviceMask = int(z'1f')
  end enum !hipDeviceFlags

  enum, bind(c) !:: hipStreamFlags
    enumerator :: hipStreamDefault = int(z'00')
    enumerator :: hipStreamNonBlocking = int(z'01')
  end enum

  enum, bind(c) !:: hipEventFlags
    enumerator :: hipEventDefault = int(z'00')
    enumerator :: hipEventBlockingSync = int(z'01')
    enumerator :: hipEventDisableTiming = int(z'02')
    enumerator :: hipEventInterprocess = int(z'04')
  end enum

  enum, bind(c) !:: hipSharedMemConfig
    enumerator :: hipSharedMemBankSizeDefault   = 0
    enumerator :: hipSharedMemBankSizeFourByte  = 1
    enumerator :: hipSharedMemBankSizeEightByte = 2
  end enum

  enum, bind(c) !:: hipFuncCache
    enumerator :: hipFuncCachePreferNone   = 0
    enumerator :: hipFuncCachePreferShared = 1
    enumerator :: hipFuncCachePreferL1     = 2
    enumerator :: hipFuncCachePreferEqual  = 3
  end enum

  enum, bind(c) !:: hipError
    enumerator :: hipSuccess                          = 0
    enumerator :: hipErrorInvalidContext              = 1
    enumerator :: hipErrorInvalidKernelFile           = 2
    enumerator :: hipErrorMemoryAllocation            = 3
    enumerator :: hipErrorInitializationError         = 4
    enumerator :: hipErrorLaunchFailure               = 5
    enumerator :: hipErrorLaunchOutOfResources        = 6
    enumerator :: hipErrorInvalidDevice               = 7
    enumerator :: hipErrorInvalidValue                = 8
    enumerator :: hipErrorInvalidDevicePointer        = 9
    enumerator :: hipErrorInvalidMemcpyDirection      = 10
    enumerator :: hipErrorUnknown                     = 11
    enumerator :: hipErrorInvalidResourceHandle       = 12
    enumerator :: hipErrorNotReady                    = 13
    enumerator :: hipErrorNoDevice                    = 14
    enumerator :: hipErrorPeerAccessAlreadyEnabled    = 15
    enumerator :: hipErrorPeerAccessNotEnabled        = 16
    enumerator :: hipErrorRuntimeMemory               = 17
    enumerator :: hipErrorRuntimeOther                = 18
    enumerator :: hipErrorHostMemoryAlreadyRegistered = 19
    enumerator :: hipErrorHostMemoryNotRegistered     = 20
    enumerator :: hipErrorMapBufferObjectFailed       = 21
    enumerator :: hipErrorTbd                         = 22
  end enum !hipError

  interface

    integer(c_int) function &
        & hipHostMalloc(cPtr, size, flags) &
        & bind(c, name="hipHostMalloc")
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: cPtr
      integer(c_size_t), value :: size
      integer(c_int), value :: flags
    end function hipHostMalloc

    integer(c_int) function &
        & hipMallocHost(cPtr, size) &
        & bind(c, name="hipMallocHost")
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: cPtr
      integer(c_size_t), value :: size
    end function hipMallocHost

    !integer(c_int) function &
    !    & hipMallocManaged(dPtr, size, flags) &
    !    & bind(c, name="hipMallocManaged")
    !  use, intrinsic :: iso_c_binding
    !  type(c_ptr) :: dPtr
    !  integer(c_size_t), value :: size
    !  integer(c_int), value :: flags
    !end function hipMallocManaged

    integer(c_int) function &
        & hipHostFree(cPtr) &
        & bind(c, name="hipHostFree")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: cPtr
    end function hipHostFree

    integer(c_int) function &
        & hipMalloc(dPtr, size) &
        & bind(c, name="hipMalloc")
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: dPtr
      integer(c_size_t), value :: size
    end function hipMalloc

    integer(c_int) function &
        & hipFree(dPtr) &
        & bind(c, name="hipFree")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: dPtr
    end function hipFree

    integer(c_int) function &
        & hipMemcpy(dst, src, memSize, cpyKind) &
        & bind(c, name="hipMemcpy")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: dst
      type(c_ptr), value :: src
      integer(c_size_t), value :: memSize
      integer(c_int), value :: cpyKind
    end function hipMemcpy

    integer(c_int) function &
        & hipMemcpyAsync(dst, src, memSize, cpyKind, stream) &
        & bind(c, name="hipMemcpyAsync")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: dst
      type(c_ptr), value :: src
      integer(c_size_t), value :: memSize
      integer(c_int), value :: cpyKind
      type(c_ptr), value :: stream
    end function hipMemcpyAsync

    integer(c_int) function &
        & hipStreamCreate(stream) &
        & bind(c, name="hipStreamCreate")
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: stream
    end function hipStreamCreate

    integer(c_int) function &
        & hipStreamCreateWithFlags(stream, flags) &
        & bind(c, name="hipStreamCreateWithFlags")
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: stream
      integer(c_int), value :: flags
    end function hipStreamCreateWithFlags

    integer(c_int) function &
        & hipStreamDestroy(stream) &
        & bind(c, name="hipStreamDestroy")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: stream
    end function hipStreamDestroy

    integer(c_int) function &
        & hipStreamSynchronize(stream) &
        & bind(c, name="hipStreamSynchronize")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: stream
    end function hipStreamSynchronize

    integer(c_int) function &
        & hipStreamWaitEvent(stream, event) &
        & bind(c, name="hipStreamWaitEvent")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: stream
      type(c_ptr), value :: event
    end function hipStreamWaitEvent

    integer(c_int) function &
        & hipGetDeviceCount(deviceCount) &
        & bind(c, name="hipGetDeviceCount")
      use, intrinsic :: iso_c_binding
      integer(c_int) :: deviceCount
    end function hipGetDeviceCount

    !integer(c_int) function &
    !    & hipGetDeviceProperties(prop, device) &
    !    & bind(c, name="hipGetDeviceProperties")
    !  use, intrinsic :: iso_c_binding
    !  import hipDeviceProp
    !  type(hipDeviceProp) :: prop
    !  integer(c_int), value :: device
    !end function hipGetDeviceProperties

    integer(c_int) function &
        & hipSetDevice(device) &
        & bind(c, name="hipSetDevice")
      use, intrinsic :: iso_c_binding
      integer(c_int), value :: device
    end function hipSetDevice

    integer(c_int) function &
        & hipSetDeviceFlags(flags) &
        & bind(c, name="hipSetDeviceFlags")
      use, intrinsic :: iso_c_binding
      integer(c_int), value :: flags
    end function hipSetDeviceFlags

    integer(c_int) function &
        & hipHostGetDevicePointer(dPtr, cPtr, flags) &
        & bind(c, name="hipHostGetDevicePointer")
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: dPtr
      type(c_ptr), value :: cPtr
      integer(c_int), value :: flags
    end function hipHostGetDevicePointer

    integer(c_int) function &
        & hipDeviceSynchronize() &
        & bind(c, name="hipDeviceSynchronize")
      use, intrinsic :: iso_c_binding
    end function hipDeviceSynchronize

    integer(c_int) function &
        & hipDeviceReset() &
        & bind(c, name="hipDeviceReset")
      use, intrinsic :: iso_c_binding
    end function hipDeviceReset

    integer(c_int) function &
        & hipEventCreate( event ) &
        & bind(c, name="hipEventCreate")
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: event
    end function hipEventCreate

    integer(c_int) function &
        & hipEventCreateWithFlags( event, flags ) &
        & bind(c, name="hipEventCreateWithFlags")
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: event
      integer(c_int), value :: flags
    end function hipEventCreateWithFlags

    integer(c_int) function &
        & hipEventDestroy( event ) &
        & bind(c, name="hipEventDestroy")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: event
    end function hipEventDestroy

    integer(c_int) function &
        & hipEventElapsedTime( eventTime, eventStart, eventStop ) &
        & bind(c, name="hipEventElapsedTime")
      use, intrinsic :: iso_c_binding
      real(c_float) :: eventTime
      type(c_ptr), value :: eventStart
      type(c_ptr), value :: eventStop
    end function hipEventElapsedTime

    integer(c_int) function &
        & hipEventQuery( event ) &
        & bind(c, name="hipEventQuery")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: event
    end function hipEventQuery

    integer(c_int) function &
        & hipEventRecord( event, stream ) &
        & bind(c, name="hipEventRecord")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: event
      type(c_ptr), value :: stream
    end function hipEventRecord

    integer(c_int) function &
        & hipEventSynchronize( event ) &
        & bind(c, name="hipEventSynchronize")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: event
    end function hipEventSynchronize

    integer(c_int) function &
        & hipDeviceSetSharedMemConfig( config ) &
        & bind(c, name="hipDeviceSetSharedMemConfig")
      use, intrinsic :: iso_c_binding
      integer(c_int), value :: config
    end function hipDeviceSetSharedMemConfig

    integer(c_int) function &
        & hipDeviceSetCacheConfig( cacheConfig ) &
        & bind(c, name="hipDeviceSetCacheConfig")
      use, intrinsic :: iso_c_binding
      integer(c_int), value :: cacheConfig
    end function hipDeviceSetCacheConfig

    integer(c_int) function &
        & hipHostRegister(cPtr, size, flags) &
        & bind(c, name="hipHostRegister")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: cPtr
      integer(c_size_t), value :: size
      integer(c_int), value :: flags
    end function hipHostRegister

    integer(c_int) function &
        & hipHostUnregister(cPtr) &
        & bind(c, name="hipHostUnregister")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: cPtr
    end function hipHostUnregister

  end interface

end module HipModule
