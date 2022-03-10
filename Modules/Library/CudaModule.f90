!***************************************************************************************************
! CudaModule.f90 10/18/17
! This file contains the module defining Fortran interfaces for the CUDA Runtime API
!***************************************************************************************************

module CudaModule
  !-------------------------------------------------------------------------------------------------
  ! Interface to CUDA Runtime API
  !-------------------------------------------------------------------------------------------------
  use, intrinsic :: iso_c_binding

  integer :: mystream
  type(c_ptr), target :: stream
  type(c_ptr) :: event

  integer :: nstream, nevent
  type(c_ptr), allocatable :: streamArray(:)
  type(c_ptr), allocatable :: eventArray(:)

  enum, bind(c) !:: cudaMemcpyKind
    enumerator :: cudaMemcpyHostToHost = 0
    enumerator :: cudaMemcpyHostToDevice = 1
    enumerator :: cudaMemcpyDeviceToHost = 2
    enumerator :: cudaMemcpyDeviceToDevice = 3
    enumerator :: cudaMemcpyDefault = 4
  end enum !cudaMemcpyKind

  enum, bind(c) !:: cudaMallocManaged
    enumerator :: cudaMemAttachGlobal = int(z'01')
    enumerator :: cudaMemAttachHost = int(z'02')
    enumerator :: cudaMemAttachSingle = int(z'04')
  end enum !cudaMallocManaged

  enum, bind(c) !:: cudaHostAllocFlags
    enumerator :: cudaHostAllocDefault = int(z'00')
    enumerator :: cudaHostAllocPortable = int(z'01')
    enumerator :: cudaHostAllocMapped = int(z'02')
    enumerator :: cudaHostAllocWriteCombined = int(z'04')
  end enum !cudaHostAllocFlags

  enum, bind(c) !:: cudaHostRegisterFlags
    enumerator :: cudaHostRegisterDefault = int(z'00')
    enumerator :: cudaHostRegisterPortable = int(z'01')
    enumerator :: cudaHostRegisterMapped = int(z'02')
    enumerator :: cudaHostRegisterIoMemory = int(z'04')
  end enum !cudaHostRegisterFlags

  enum, bind(c) !:: cudaDeviceFlags
    enumerator :: cudaDeviceScheduleAuto = int(z'00')
    enumerator :: cudaDeviceScheduleSpin = int(z'01')
    enumerator :: cudaDeviceScheduleYield = int(z'02')
    enumerator :: cudaDeviceScheduleBlockingSync = int(z'04')
    enumerator :: cudaDeviceBlockingSync = int(z'04')
    enumerator :: cudaDeviceScheduleMask = int(z'07')
    enumerator :: cudaDeviceMapHost = int(z'08')
    enumerator :: cudaDeviceLmemResizeToMax = int(z'10')
    enumerator :: cudaDeviceMask = int(z'1f')
  end enum !cudaDeviceFlags

  enum, bind(c) !:: cudaStreamFlags
    enumerator :: cudaStreamDefault = int(z'00')
    enumerator :: cudaStreamNonBlocking = int(z'01')
  end enum

  enum, bind(c) !:: cudaEventFlags
    enumerator :: cudaEventDefault = int(z'00')
    enumerator :: cudaEventBlockingSync = int(z'01')
    enumerator :: cudaEventDisableTiming = int(z'02')
    enumerator :: cudaEventInterprocess = int(z'04')
  end enum

  enum, bind(c) !:: cudaSharedMemConfig
    enumerator :: cudaSharedMemBankSizeDefault   = 0
    enumerator :: cudaSharedMemBankSizeFourByte  = 1
    enumerator :: cudaSharedMemBankSizeEightByte = 2
  end enum

  enum, bind(c) !:: cudaFuncCache
    enumerator :: cudaFuncCachePreferNone   = 0
    enumerator :: cudaFuncCachePreferShared = 1
    enumerator :: cudaFuncCachePreferL1     = 2
    enumerator :: cudaFuncCachePreferEqual  = 3
  end enum

  enum, bind(c) !:: cudaError
    enumerator :: cudaSuccess                           = 0
    enumerator :: cudaErrorMissingConfiguration         = 1
    enumerator :: cudaErrorMemoryAllocation             = 2
    enumerator :: cudaErrorInitializationError          = 3
    enumerator :: cudaErrorLaunchFailure                = 4
    enumerator :: cudaErrorPriorLaunchFailure           = 5
    enumerator :: cudaErrorLaunchTimeout                = 6
    enumerator :: cudaErrorLaunchOutOfResources         = 7
    enumerator :: cudaErrorInvalidDeviceFunction        = 8
    enumerator :: cudaErrorInvalidConfiguration         = 9
    enumerator :: cudaErrorInvalidDevice                = 10
    enumerator :: cudaErrorInvalidValue                 = 11
    enumerator :: cudaErrorInvalidPitchValue            = 12
    enumerator :: cudaErrorInvalidSymbol                = 13
    enumerator :: cudaErrorMapBufferObjectFailed        = 14
    enumerator :: cudaErrorUnmapBufferObjectFailed      = 15
    enumerator :: cudaErrorInvalidHostPointer           = 16
    enumerator :: cudaErrorInvalidDevicePointer         = 17
    enumerator :: cudaErrorInvalidTexture               = 18
    enumerator :: cudaErrorInvalidTextureBinding        = 19
    enumerator :: cudaErrorInvalidChannelDescriptor     = 20
    enumerator :: cudaErrorInvalidMemcpyDirection       = 21
    enumerator :: cudaErrorAddressOfConstant            = 22
    enumerator :: cudaErrorTextureFetchFailed           = 23
    enumerator :: cudaErrorTextureNotBound              = 24
    enumerator :: cudaErrorSynchronizationError         = 25
    enumerator :: cudaErrorInvalidFilterSetting         = 26
    enumerator :: cudaErrorInvalidNormSetting           = 27
    enumerator :: cudaErrorMixedDeviceExecution         = 28
    enumerator :: cudaErrorCudartUnloading              = 29
    enumerator :: cudaErrorUnknown                      = 30
    enumerator :: cudaErrorNotYetImplemented            = 31
    enumerator :: cudaErrorMemoryValueTooLarge          = 32
    enumerator :: cudaErrorInvalidResourceHandle        = 33
    enumerator :: cudaErrorNotReady                     = 34
    enumerator :: cudaErrorInsufficientDriver           = 35
    enumerator :: cudaErrorSetOnActiveProcess           = 36
    enumerator :: cudaErrorInvalidSurface               = 37
    enumerator :: cudaErrorNoDevice                     = 38
    enumerator :: cudaErrorECCUncorrectable             = 39
    enumerator :: cudaErrorSharedObjectSymbolNotFound   = 40
    enumerator :: cudaErrorSharedObjectInitFailed       = 41
    enumerator :: cudaErrorUnsupportedLimit             = 42
    enumerator :: cudaErrorDuplicateVariableName        = 43
    enumerator :: cudaErrorDuplicateTextureName         = 44
    enumerator :: cudaErrorDuplicateSurfaceName         = 45
    enumerator :: cudaErrorDevicesUnavailable           = 46
    enumerator :: cudaErrorInvalidKernelImage           = 47
    enumerator :: cudaErrorNoKernelImageForDevice       = 48
    enumerator :: cudaErrorIncompatibleDriverContext    = 49
    enumerator :: cudaErrorPeerAccessAlreadyEnabled     = 50
    enumerator :: cudaErrorPeerAccessNotEnabled         = 51
    enumerator :: cudaErrorDeviceAlreadyInUse           = 54
    enumerator :: cudaErrorProfilerDisabled             = 55
    enumerator :: cudaErrorProfilerNotInitialized       = 56
    enumerator :: cudaErrorProfilerAlreadyStarted       = 57
    enumerator :: cudaErrorProfilerAlreadyStopped       = 58
    enumerator :: cudaErrorAssert                       = 59
    enumerator :: cudaErrorTooManyPeers                 = 60
    enumerator :: cudaErrorHostMemoryAlreadyRegistered  = 61
    enumerator :: cudaErrorHostMemoryNotRegistered      = 62
    enumerator :: cudaErrorOperatingSystem              = 63
    enumerator :: cudaErrorPeerAccessUnsupported        = 64
    enumerator :: cudaErrorLaunchMaxDepthExceeded       = 65
    enumerator :: cudaErrorLaunchFileScopedTex          = 66
    enumerator :: cudaErrorLaunchFileScopedSurf         = 67
    enumerator :: cudaErrorSyncDepthExceeded            = 68
    enumerator :: cudaErrorLaunchPendingCountExceeded   = 69
    enumerator :: cudaErrorNotPermitted                 = 70
    enumerator :: cudaErrorNotSupported                 = 71
    enumerator :: cudaErrorStartupFailure               = int(z'7f')
    enumerator :: cudaErrorApiFailureBase               = 10000
  end enum !cudaError

  !include "cudaDeviceProp.fh"

  interface

    integer(c_int) function &
        & cudaHostAlloc(cPtr, size, flags) &
        & bind(c, name="cudaHostAlloc")
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: cPtr
      integer(c_size_t), value :: size
      integer(c_int), value :: flags
    end function cudaHostAlloc

    integer(c_int) function &
        & cudaMallocHost(cPtr, size) &
        & bind(c, name="cudaMallocHost")
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: cPtr
      integer(c_size_t), value :: size
    end function cudaMallocHost

    integer(c_int) function &
        & cudaMallocManaged(dPtr, size, flags) &
        & bind(c, name="cudaMallocManaged")
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: dPtr
      integer(c_size_t), value :: size
      integer(c_int), value :: flags
    end function cudaMallocManaged

    integer(c_int) function &
        & cudaFreeHost(cPtr) &
        & bind(c, name="cudaFreeHost")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: cPtr
    end function cudaFreeHost

    integer(c_int) function &
        & cudaMalloc(dPtr, size) &
        & bind(c, name="cudaMalloc")
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: dPtr
      integer(c_size_t), value :: size
    end function cudaMalloc

    integer(c_int) function &
        & cudaFree(dPtr) &
        & bind(c, name="cudaFree")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: dPtr
    end function cudaFree

    integer(c_int) function &
        & cudaMemcpy(dst, src, memSize, cpyKind) &
        & bind(c, name="cudaMemcpy")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: dst
      type(c_ptr), value :: src
      integer(c_size_t), value :: memSize
      integer(c_int), value :: cpyKind
    end function cudaMemcpy

    integer(c_int) function &
        & cudaMemcpyAsync(dst, src, memSize, cpyKind, stream) &
        & bind(c, name="cudaMemcpyAsync")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: dst
      type(c_ptr), value :: src
      integer(c_size_t), value :: memSize
      integer(c_int), value :: cpyKind
      type(c_ptr), value :: stream
    end function cudaMemcpyAsync

    integer(c_int) function &
        & cudaStreamCreate(stream) &
        & bind(c, name="cudaStreamCreate")
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: stream
    end function cudaStreamCreate

    integer(c_int) function &
        & cudaStreamCreateWithFlags(stream, flags) &
        & bind(c, name="cudaStreamCreateWithFlags")
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: stream
      integer(c_int), value :: flags
    end function cudaStreamCreateWithFlags

    integer(c_int) function &
        & cudaStreamDestroy(stream) &
        & bind(c, name="cudaStreamDestroy")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: stream
    end function cudaStreamDestroy

    integer(c_int) function &
        & cudaStreamSynchronize(stream) &
        & bind(c, name="cudaStreamSynchronize")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: stream
    end function cudaStreamSynchronize

    integer(c_int) function &
        & cudaStreamWaitEvent(stream, event) &
        & bind(c, name="cudaStreamWaitEvent")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: stream
      type(c_ptr), value :: event
    end function cudaStreamWaitEvent

    integer(c_int) function &
        & cudaGetDeviceCount(deviceCount) &
        & bind(c, name="cudaGetDeviceCount")
      use, intrinsic :: iso_c_binding
      integer(c_int) :: deviceCount
    end function cudaGetDeviceCount

    !integer(c_int) function &
    !    & cudaGetDeviceProperties(prop, device) &
    !    & bind(c, name="cudaGetDeviceProperties")
    !  use, intrinsic :: iso_c_binding
    !  import cudaDeviceProp
    !  type(cudaDeviceProp) :: prop
    !  integer(c_int), value :: device
    !end function cudaGetDeviceProperties

    integer(c_int) function &
        & cudaSetDevice(device) &
        & bind(c, name="cudaSetDevice")
      use, intrinsic :: iso_c_binding
      integer(c_int), value :: device
    end function cudaSetDevice

    integer(c_int) function &
        & cudaSetDeviceFlags(flags) &
        & bind(c, name="cudaSetDeviceFlags")
      use, intrinsic :: iso_c_binding
      integer(c_int), value :: flags
    end function cudaSetDeviceFlags

    integer(c_int) function &
        & cudaHostGetDevicePointer(dPtr, cPtr, flags) &
        & bind(c, name="cudaHostGetDevicePointer")
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: dPtr
      type(c_ptr), value :: cPtr
      integer(c_int), value :: flags
    end function cudaHostGetDevicePointer

    integer(c_int) function &
        & cudaDeviceSynchronize() &
        & bind(c, name="cudaDeviceSynchronize")
      use, intrinsic :: iso_c_binding
    end function cudaDeviceSynchronize

    integer(c_int) function &
        & cudaDeviceReset() &
        & bind(c, name="cudaDeviceReset")
      use, intrinsic :: iso_c_binding
    end function cudaDeviceReset

    integer(c_int) function &
        & cudaEventCreate( event ) &
        & bind(c, name="cudaEventCreate")
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: event
    end function cudaEventCreate

    integer(c_int) function &
        & cudaEventCreateWithFlags( event, flags ) &
        & bind(c, name="cudaEventCreateWithFlags")
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: event
      integer(c_int), value :: flags
    end function cudaEventCreateWithFlags

    integer(c_int) function &
        & cudaEventDestroy( event ) &
        & bind(c, name="cudaEventDestroy")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: event
    end function cudaEventDestroy

    integer(c_int) function &
        & cudaEventElapsedTime( eventTime, eventStart, eventStop ) &
        & bind(c, name="cudaEventElapsedTime")
      use, intrinsic :: iso_c_binding
      real(c_float) :: eventTime
      type(c_ptr), value :: eventStart
      type(c_ptr), value :: eventStop
    end function cudaEventElapsedTime

    integer(c_int) function &
        & cudaEventQuery( event ) &
        & bind(c, name="cudaEventQuery")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: event
    end function cudaEventQuery

    integer(c_int) function &
        & cudaEventRecord( event, stream ) &
        & bind(c, name="cudaEventRecord")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: event
      type(c_ptr), value :: stream
    end function cudaEventRecord

    integer(c_int) function &
        & cudaEventSynchronize( event ) &
        & bind(c, name="cudaEventSynchronize")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: event
    end function cudaEventSynchronize

    integer(c_int) function &
        & cudaDeviceSetSharedMemConfig( config ) &
        & bind(c, name="cudaDeviceSetSharedMemConfig")
      use, intrinsic :: iso_c_binding
      integer(c_int), value :: config
    end function cudaDeviceSetSharedMemConfig

    integer(c_int) function &
        & cudaDeviceSetCacheConfig( cacheConfig ) &
        & bind(c, name="cudaDeviceSetCacheConfig")
      use, intrinsic :: iso_c_binding
      integer(c_int), value :: cacheConfig
    end function cudaDeviceSetCacheConfig

    integer(c_int) function &
        & cudaHostRegister(cPtr, size, flags) &
        & bind(c, name="cudaHostRegister")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: cPtr
      integer(c_size_t), value :: size
      integer(c_int), value :: flags
    end function cudaHostRegister

    integer(c_int) function &
        & cudaHostUnregister(cPtr) &
        & bind(c, name="cudaHostUnregister")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: cPtr
    end function cudaHostUnregister

  end interface

end module CudaModule
