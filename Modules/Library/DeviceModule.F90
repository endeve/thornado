MODULE DeviceModule

  USE, INTRINSIC :: ISO_C_BINDING
  USE KindModule, ONLY: &
    DP

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
    cublasGetStream_v2, &
    cublasSetStream_v2
  USE CusolverModule, ONLY: &
    cusolver_handle, &
    cusolverDnCreate, &
    cusolverDnSetStream
  USE CusparseModule, ONLY: &
    cusparse_handle, &
    cusparseCreate, &
    cusparseSetStream
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
    omp_is_initial_device, &
    omp_target_is_present
#endif

#if defined(THORNADO_OACC)
  USE OpenACCModule, ONLY: &
    acc_get_device_num, &
    acc_on_device, &
    acc_is_present, &
    acc_set_cuda_stream, &
    acc_get_default_async, &
    acc_device_host, &
    acc_device_nvidia, &
    acc_async_default
#endif

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  INTEGER, PUBLIC :: mydevice, ndevices

  INTERFACE dev_ptr
    MODULE PROCEDURE dev_ptr_int
    MODULE PROCEDURE dev_ptr_dp
    MODULE PROCEDURE dev_ptr_cptr
  END INTERFACE

  INTERFACE QueryOnGPU
    MODULE PROCEDURE QueryOnGPU_3D_DP_1
    MODULE PROCEDURE QueryOnGPU_3D_DP_2
    MODULE PROCEDURE QueryOnGPU_3D_DP_3
    MODULE PROCEDURE QueryOnGPU_2D_DP_1
    MODULE PROCEDURE QueryOnGPU_2D_DP_2
    MODULE PROCEDURE QueryOnGPU_2D_DP_3
    MODULE PROCEDURE QueryOnGPU_DP_1
    MODULE PROCEDURE QueryOnGPU_DP_2
    MODULE PROCEDURE QueryOnGPU_DP_3
    MODULE PROCEDURE QueryOnGPU_DP_4
    MODULE PROCEDURE QueryOnGPU_DP_5
    MODULE PROCEDURE QueryOnGPU_DP_6
    MODULE PROCEDURE QueryOnGPU_DP_7
  END INTERFACE

  PUBLIC :: InitializeDevice
  PUBLIC :: FinalizeDevice
  PUBLIC :: device_is_present
  PUBLIC :: get_device_num
  PUBLIC :: on_device
  PUBLIC :: dev_ptr
  PUBLIC :: QueryOnGpu

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
    ierr = cudaStreamCreate( stream )
    ierr = cublasSetStream_v2( cublas_handle, stream )
    !ierr = cublasGetStream_v2( cublas_handle, stream )

    ierr = cusolverDnCreate( cusolver_handle )
    ierr = cusolverDnSetStream( cusolver_handle, stream )

    ierr = cusparseCreate( cusparse_handle )
    ierr = cusparseSetStream( cusparse_handle, stream )
#endif

#if defined(THORNADO_LA_MAGMA)
    CALL magma_getdevice( magma_device )
    CALL magma_init()
    CALL magma_queue_create_from_cuda &
           ( magma_device, stream, cublas_handle, cusparse_handle, magma_queue )
#endif

#if defined(THORNADO_OMP_OL)
    CALL omp_set_default_device( mydevice )
#endif

#if defined(THORNADO_OACC)
    acc_async_default = acc_get_default_async()
    ierr = acc_set_cuda_stream( acc_async_default, stream )
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
    get_device_num = acc_get_device_num( acc_device_nvidia )
#else
    get_device_num = -1
#endif
    RETURN
  END FUNCTION get_device_num


  LOGICAL FUNCTION on_device()
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
    on_device = ( .not. omp_is_initial_device() )
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
    on_device = ( .not. acc_on_device( acc_device_host ) )
#else
    on_device = .false.
#endif
    RETURN
  END FUNCTION on_device

  TYPE(C_PTR) FUNCTION dev_ptr_int( a )
    INTEGER, TARGET, INTENT(IN) :: a
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET DATA USE_DEVICE_PTR( a )
#elif defined(THORNADO_OACC)
    !$ACC HOST_DATA USE_DEVICE( a )
#endif
    dev_ptr_int = C_LOC( a )
#if defined(THORNADO_OMP_OL)
    !$OMP END TARGET DATA
#elif defined(THORNADO_OACC)
    !$ACC END HOST_DATA
#endif
  END FUNCTION dev_ptr_int

  TYPE(C_PTR) FUNCTION dev_ptr_dp( a )
    REAL(DP), TARGET, INTENT(IN) :: a
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET DATA USE_DEVICE_PTR( a )
#elif defined(THORNADO_OACC)
    !$ACC HOST_DATA USE_DEVICE( a )
#endif
    dev_ptr_dp = C_LOC( a )
#if defined(THORNADO_OMP_OL)
    !$OMP END TARGET DATA
#elif defined(THORNADO_OACC)
    !$ACC END HOST_DATA
#endif
  END FUNCTION dev_ptr_dp

  TYPE(C_PTR) FUNCTION dev_ptr_cptr( a )
    TYPE(C_PTR), TARGET, INTENT(IN) :: a
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET DATA USE_DEVICE_PTR( a )
#elif defined(THORNADO_OACC)
    !$ACC HOST_DATA USE_DEVICE( a )
#endif
    dev_ptr_cptr = C_LOC( a )
#if defined(THORNADO_OMP_OL)
    !$OMP END TARGET DATA
#elif defined(THORNADO_OACC)
    !$ACC END HOST_DATA
#endif
  END FUNCTION dev_ptr_cptr


  FUNCTION QueryOnGPU_3D_DP_1( X1 ) RESULT( QueryOnGPU )

    REAL(DP), DIMENSION(:,:,:), INTENT(in), TARGET :: X1
    LOGICAL :: QueryOnGPU

    INTEGER(C_SIZE_T) :: SizeOf_X1

    SizeOf_X1 = SIZE(X1) * C_SIZEOF(0.0_DP)

    QueryOnGPU = device_is_present( C_LOC( X1 ), mydevice, SizeOf_X1 )

  END FUNCTION QueryOnGPU_3D_DP_1


  FUNCTION QueryOnGPU_3D_DP_2( X1, X2 ) RESULT( QueryOnGPU )

    REAL(DP), DIMENSION(:,:,:), INTENT(in), TARGET :: X1, X2
    LOGICAL :: QueryOnGPU

    INTEGER(C_SIZE_T) :: SizeOf_X1
    INTEGER(C_SIZE_T) :: SizeOf_X2

    SizeOf_X1 = SIZE(X1) * C_SIZEOF(0.0_DP)
    SizeOf_X2 = SIZE(X2) * C_SIZEOF(0.0_DP)

    QueryOnGPU = device_is_present( C_LOC( X1 ), mydevice, SizeOf_X1 ) &
           .AND. device_is_present( C_LOC( X2 ), mydevice, SizeOf_X2 )

  END FUNCTION QueryOnGPU_3D_DP_2


  FUNCTION QueryOnGPU_3D_DP_3( X1, X2, X3 ) RESULT( QueryOnGPU )

    REAL(DP), DIMENSION(:,:,:), INTENT(in), TARGET :: X1, X2, X3
    LOGICAL :: QueryOnGPU

    INTEGER(C_SIZE_T) :: SizeOf_X1
    INTEGER(C_SIZE_T) :: SizeOf_X2
    INTEGER(C_SIZE_T) :: SizeOf_X3

    SizeOf_X1 = SIZE(X1) * C_SIZEOF(0.0_DP)
    SizeOf_X2 = SIZE(X2) * C_SIZEOF(0.0_DP)
    SizeOf_X3 = SIZE(X3) * C_SIZEOF(0.0_DP)

    QueryOnGPU = device_is_present( C_LOC( X1 ), mydevice, SizeOf_X1 ) &
           .AND. device_is_present( C_LOC( X2 ), mydevice, SizeOf_X2 ) &
           .AND. device_is_present( C_LOC( X3 ), mydevice, SizeOf_X3 )

  END FUNCTION QueryOnGPU_3D_DP_3


  FUNCTION QueryOnGPU_2D_DP_1( X1 ) RESULT( QueryOnGPU )

    REAL(DP), DIMENSION(:,:), INTENT(in), TARGET :: X1
    LOGICAL :: QueryOnGPU

    INTEGER(C_SIZE_T) :: SizeOf_X1

    SizeOf_X1 = SIZE(X1) * C_SIZEOF(0.0_DP)

    QueryOnGPU = device_is_present( C_LOC( X1 ), mydevice, SizeOf_X1 )

  END FUNCTION QueryOnGPU_2D_DP_1


  FUNCTION QueryOnGPU_2D_DP_2( X1, X2 ) RESULT( QueryOnGPU )

    REAL(DP), DIMENSION(:,:), INTENT(in), TARGET :: X1, X2
    LOGICAL :: QueryOnGPU

    INTEGER(C_SIZE_T) :: SizeOf_X1
    INTEGER(C_SIZE_T) :: SizeOf_X2

    SizeOf_X1 = SIZE(X1) * C_SIZEOF(0.0_DP)
    SizeOf_X2 = SIZE(X2) * C_SIZEOF(0.0_DP)

    QueryOnGPU = device_is_present( C_LOC( X1 ), mydevice, SizeOf_X1 ) &
           .AND. device_is_present( C_LOC( X2 ), mydevice, SizeOf_X2 )

  END FUNCTION QueryOnGPU_2D_DP_2


  FUNCTION QueryOnGPU_2D_DP_3( X1, X2, X3 ) RESULT( QueryOnGPU )

    REAL(DP), DIMENSION(:,:), INTENT(in), TARGET :: X1, X2, X3
    LOGICAL :: QueryOnGPU

    INTEGER(C_SIZE_T) :: SizeOf_X1
    INTEGER(C_SIZE_T) :: SizeOf_X2
    INTEGER(C_SIZE_T) :: SizeOf_X3

    SizeOf_X1 = SIZE(X1) * C_SIZEOF(0.0_DP)
    SizeOf_X2 = SIZE(X2) * C_SIZEOF(0.0_DP)
    SizeOf_X3 = SIZE(X3) * C_SIZEOF(0.0_DP)

    QueryOnGPU = device_is_present( C_LOC( X1 ), mydevice, SizeOf_X1 ) &
           .AND. device_is_present( C_LOC( X2 ), mydevice, SizeOf_X2 ) &
           .AND. device_is_present( C_LOC( X3 ), mydevice, SizeOf_X3 )

  END FUNCTION QueryOnGPU_2D_DP_3


  FUNCTION QueryOnGPU_DP_1( X1 ) RESULT( QueryOnGPU )

    REAL(DP), DIMENSION(:), INTENT(in), TARGET :: X1
    LOGICAL :: QueryOnGPU

    INTEGER(C_SIZE_T) :: SizeOf_X1

    SizeOf_X1 = SIZE(X1) * C_SIZEOF(0.0_DP)

    QueryOnGPU = device_is_present( C_LOC( X1 ), mydevice, SizeOf_X1 )

  END FUNCTION QueryOnGPU_DP_1


  FUNCTION QueryOnGPU_DP_2( X1, X2 ) RESULT( QueryOnGPU )

    REAL(DP), DIMENSION(:), INTENT(in), TARGET :: X1, X2
    LOGICAL :: QueryOnGPU

    INTEGER(C_SIZE_T) :: SizeOf_X1
    INTEGER(C_SIZE_T) :: SizeOf_X2

    SizeOf_X1 = SIZE(X1) * C_SIZEOF(0.0_DP)
    SizeOf_X2 = SIZE(X2) * C_SIZEOF(0.0_DP)

    QueryOnGPU = device_is_present( C_LOC( X1 ), mydevice, SizeOf_X1 ) &
           .AND. device_is_present( C_LOC( X2 ), mydevice, SizeOf_X2 )

  END FUNCTION QueryOnGPU_DP_2


  FUNCTION QueryOnGPU_DP_3( X1, X2, X3 ) RESULT( QueryOnGPU )

    REAL(DP), DIMENSION(:), INTENT(in), TARGET :: X1, X2, X3
    LOGICAL :: QueryOnGPU

    INTEGER(C_SIZE_T) :: SizeOf_X1
    INTEGER(C_SIZE_T) :: SizeOf_X2
    INTEGER(C_SIZE_T) :: SizeOf_X3

    SizeOf_X1 = SIZE(X1) * C_SIZEOF(0.0_DP)
    SizeOf_X2 = SIZE(X2) * C_SIZEOF(0.0_DP)
    SizeOf_X3 = SIZE(X3) * C_SIZEOF(0.0_DP)

    QueryOnGPU = device_is_present( C_LOC( X1 ), mydevice, SizeOf_X1 ) &
           .AND. device_is_present( C_LOC( X2 ), mydevice, SizeOf_X2 ) &
           .AND. device_is_present( C_LOC( X3 ), mydevice, SizeOf_X3 )

  END FUNCTION QueryOnGPU_DP_3


  FUNCTION QueryOnGPU_DP_4( X1, X2, X3, X4 ) RESULT( QueryOnGPU )

    REAL(DP), DIMENSION(:), INTENT(in), TARGET :: X1, X2, X3, X4
    LOGICAL :: QueryOnGPU

    INTEGER(C_SIZE_T) :: SizeOf_X1
    INTEGER(C_SIZE_T) :: SizeOf_X2
    INTEGER(C_SIZE_T) :: SizeOf_X3
    INTEGER(C_SIZE_T) :: SizeOf_X4

    SizeOf_X1 = SIZE(X1) * C_SIZEOF(0.0_DP)
    SizeOf_X2 = SIZE(X2) * C_SIZEOF(0.0_DP)
    SizeOf_X3 = SIZE(X3) * C_SIZEOF(0.0_DP)
    SizeOf_X4 = SIZE(X4) * C_SIZEOF(0.0_DP)

    QueryOnGPU = device_is_present( C_LOC( X1 ), mydevice, SizeOf_X1 ) &
           .AND. device_is_present( C_LOC( X2 ), mydevice, SizeOf_X2 ) &
           .AND. device_is_present( C_LOC( X3 ), mydevice, SizeOf_X3 ) &
           .AND. device_is_present( C_LOC( X4 ), mydevice, SizeOf_X4 )

  END FUNCTION QueryOnGPU_DP_4


  FUNCTION QueryOnGPU_DP_5( X1, X2, X3, X4, X5 ) RESULT( QueryOnGPU )

    REAL(DP), DIMENSION(:), INTENT(in), TARGET :: X1, X2, X3, X4, X5
    LOGICAL :: QueryOnGPU

    INTEGER(C_SIZE_T) :: SizeOf_X1
    INTEGER(C_SIZE_T) :: SizeOf_X2
    INTEGER(C_SIZE_T) :: SizeOf_X3
    INTEGER(C_SIZE_T) :: SizeOf_X4
    INTEGER(C_SIZE_T) :: SizeOf_X5

    SizeOf_X1 = SIZE(X1) * C_SIZEOF(0.0_DP)
    SizeOf_X2 = SIZE(X2) * C_SIZEOF(0.0_DP)
    SizeOf_X3 = SIZE(X3) * C_SIZEOF(0.0_DP)
    SizeOf_X4 = SIZE(X4) * C_SIZEOF(0.0_DP)
    SizeOf_X5 = SIZE(X5) * C_SIZEOF(0.0_DP)

    QueryOnGPU = device_is_present( C_LOC( X1 ), mydevice, SizeOf_X1 ) &
           .AND. device_is_present( C_LOC( X2 ), mydevice, SizeOf_X2 ) &
           .AND. device_is_present( C_LOC( X3 ), mydevice, SizeOf_X3 ) &
           .AND. device_is_present( C_LOC( X4 ), mydevice, SizeOf_X4 ) &
           .AND. device_is_present( C_LOC( X5 ), mydevice, SizeOf_X5 )

  END FUNCTION QueryOnGPU_DP_5


  FUNCTION QueryOnGPU_DP_6( X1, X2, X3, X4, X5, X6 ) RESULT( QueryOnGPU )

    REAL(DP), DIMENSION(:), INTENT(in), TARGET :: X1, X2, X3, X4, X5, X6
    LOGICAL :: QueryOnGPU

    INTEGER(C_SIZE_T) :: SizeOf_X1
    INTEGER(C_SIZE_T) :: SizeOf_X2
    INTEGER(C_SIZE_T) :: SizeOf_X3
    INTEGER(C_SIZE_T) :: SizeOf_X4
    INTEGER(C_SIZE_T) :: SizeOf_X5
    INTEGER(C_SIZE_T) :: SizeOf_X6

    SizeOf_X1 = SIZE(X1) * C_SIZEOF(0.0_DP)
    SizeOf_X2 = SIZE(X2) * C_SIZEOF(0.0_DP)
    SizeOf_X3 = SIZE(X3) * C_SIZEOF(0.0_DP)
    SizeOf_X4 = SIZE(X4) * C_SIZEOF(0.0_DP)
    SizeOf_X5 = SIZE(X5) * C_SIZEOF(0.0_DP)
    SizeOf_X6 = SIZE(X6) * C_SIZEOF(0.0_DP)

    QueryOnGPU = device_is_present( C_LOC( X1 ), mydevice, SizeOf_X1 ) &
           .AND. device_is_present( C_LOC( X2 ), mydevice, SizeOf_X2 ) &
           .AND. device_is_present( C_LOC( X3 ), mydevice, SizeOf_X3 ) &
           .AND. device_is_present( C_LOC( X4 ), mydevice, SizeOf_X4 ) &
           .AND. device_is_present( C_LOC( X5 ), mydevice, SizeOf_X5 ) &
           .AND. device_is_present( C_LOC( X6 ), mydevice, SizeOf_X6 )

  END FUNCTION QueryOnGPU_DP_6


  FUNCTION QueryOnGPU_DP_7( X1, X2, X3, X4, X5, X6, X7 ) RESULT( QueryOnGPU )

    REAL(DP), DIMENSION(:), INTENT(in), TARGET :: X1, X2, X3, X4, X5, X6, X7
    LOGICAL :: QueryOnGPU

    INTEGER(C_SIZE_T) :: SizeOf_X1
    INTEGER(C_SIZE_T) :: SizeOf_X2
    INTEGER(C_SIZE_T) :: SizeOf_X3
    INTEGER(C_SIZE_T) :: SizeOf_X4
    INTEGER(C_SIZE_T) :: SizeOf_X5
    INTEGER(C_SIZE_T) :: SizeOf_X6
    INTEGER(C_SIZE_T) :: SizeOf_X7

    SizeOf_X1 = SIZE(X1) * C_SIZEOF(0.0_DP)
    SizeOf_X2 = SIZE(X2) * C_SIZEOF(0.0_DP)
    SizeOf_X3 = SIZE(X3) * C_SIZEOF(0.0_DP)
    SizeOf_X4 = SIZE(X4) * C_SIZEOF(0.0_DP)
    SizeOf_X5 = SIZE(X5) * C_SIZEOF(0.0_DP)
    SizeOf_X6 = SIZE(X6) * C_SIZEOF(0.0_DP)
    SizeOf_X7 = SIZE(X7) * C_SIZEOF(0.0_DP)

    QueryOnGPU = device_is_present( C_LOC( X1 ), mydevice, SizeOf_X1 ) &
           .AND. device_is_present( C_LOC( X2 ), mydevice, SizeOf_X2 ) &
           .AND. device_is_present( C_LOC( X3 ), mydevice, SizeOf_X3 ) &
           .AND. device_is_present( C_LOC( X4 ), mydevice, SizeOf_X4 ) &
           .AND. device_is_present( C_LOC( X5 ), mydevice, SizeOf_X5 ) &
           .AND. device_is_present( C_LOC( X6 ), mydevice, SizeOf_X6 ) &
           .AND. device_is_present( C_LOC( X7 ), mydevice, SizeOf_X7 )

  END FUNCTION QueryOnGPU_DP_7


END MODULE DeviceModule
