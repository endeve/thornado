MODULE DeviceModule

  USE, INTRINSIC :: ISO_C_BINDING
#if defined(THORNADO_GPU)
  USE CudaModule, ONLY: &
    stream, &
    cudaGetDeviceCount, &
    cudaSetDevice, &
    cudaStreamCreate
#endif

#if defined(THORNADO_LA_CUBLAS) || defined(THORNADO_LA_MAGMA)
  USE CublasModule, ONLY: &
    cublas_handle, &
    cublasCreate_v2, &
    cublasGetStream_v2
#endif

#if defined(THORNADO_LA_MAGMA)
  USE MagmaModule, ONLY: &
    magma_device, &
    magma_queue, &
    magma_getdevice, &
    magma_init, &
    magma_queue_create_from_cuda
#endif

#if defined(THORNADO_OMP_OL)
  USE OpenMPModule, ONLY: &
    omp_set_default_device, &
    omp_get_default_device, &
    omp_target_is_present
#endif

#if defined(THORNADO_OACC)
  USE OpenACCModule, ONLY: &
    acc_set_device_num, &
    acc_get_device_num, &
    acc_device_default, &
    acc_is_present
#endif

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  INTEGER, PUBLIC :: mydevice, ndevices

  PUBLIC :: InitializeDevice
  PUBLIC :: FinalizeDevice
  PUBLIC :: device_is_present
  PUBLIC :: get_device_num

CONTAINS


  SUBROUTINE InitializeDevice

    INTEGER :: ierr, myrank, nranks

#if defined(THORNADO_GPU)
    CALL MPI_COMM_RANK( MPI_COMM_WORLD, myrank, ierr )
    CALL MPI_COMM_SIZE( MPI_COMM_WORLD, nranks, ierr )
    ierr = cudaGetDeviceCount( ndevices )
    IF ( ndevices > 0 ) THEN
      mydevice = MOD( myrank, ndevices )
    ELSE
      WRITE(*,*) 'No CUDA capable device found'
      CALL MPI_FINALIZE( ierr )
    END IF
    ierr = cudaSetDevice( mydevice )
#else
    mydevice = -1
    ndevices = 0
#endif

#if defined(THORNADO_LA_CUBLAS) || defined(THORNADO_LA_MAGMA)
    ierr = cublasCreate_v2( cublas_handle )
    !ierr = cudaStreamCreate(stream)
    ierr = cublasGetStream_v2( cublas_handle, stream )
#endif

#if defined(THORNADO_LA_MAGMA)
    CALL magma_getdevice( magma_device )
    CALL magma_init()
    CALL magma_queue_create_from_cuda &
           ( magma_device, stream, cublas_handle, C_NULL_PTR, magma_queue )
#endif

#if defined(THORNADO_OMP_OL)
    CALL omp_set_default_device( mydevice )
#endif

#if defined(THORNADO_OACC)
    CALL acc_set_device_num( mydevice, acc_device_default )
#endif

    RETURN
  END SUBROUTINE InitializeDevice


  SUBROUTINE FinalizeDevice
    RETURN
  END SUBROUTINE FinalizeDevice


  LOGICAL FUNCTION device_is_present( hostptr, device, bytes )
    TYPE(C_PTR), INTENT(in) :: hostptr
    INTEGER, INTENT(in) :: device
    INTEGER(C_SIZE_T), INTENT(in) :: bytes
#if defined(THORNADO_OMP_OL)
    device_is_present = ( omp_target_is_present( hostptr, device ) > 0 )
#elif defined(THORNADO_OACC)
    device_is_present = ( acc_is_present( hostptr, bytes ) > 0 )
#else
    device_is_present = .false.
#endif
    RETURN
  END FUNCTION device_is_present


  INTEGER FUNCTION get_device_num()
#if defined(THORNADO_OMP_OL)
    get_device_num = omp_get_default_device()
#elif defined(THORNADO_OACC)
    get_device_num = acc_get_device_num( acc_device_default )
#else
    get_device_num = -1
#endif
    RETURN
  END FUNCTION get_device_num


END MODULE DeviceModule
