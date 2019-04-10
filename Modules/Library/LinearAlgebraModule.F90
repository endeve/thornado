MODULE LinearAlgebraModule

  USE, INTRINSIC :: ISO_C_BINDING
  USE KindModule, ONLY: &
    DP
  USE DeviceModule, ONLY: &
    mydevice, &
    device_is_present

#if defined(THORNADO_LA_CUBLAS)
  USE CudaModule, ONLY: &
    stream, &
    cudaStreamSynchronize
  USE CublasModule, ONLY: &
    cublas_handle, &
    cublasDgemm_v2, &
    cublasDgemv_v2, &
    CUBLAS_OP_N, CUBLAS_OP_T
#endif

#if defined(THORNADO_LA_MAGMA)
  USE MagmaModule, ONLY: &
    magma_queue, &
    magma_queue_sync, &
    magma_dgemm, &
    magma_dgemv, &
    MagmaNoTrans, MagmaTrans
#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MatrixMatrixMultiply
  PUBLIC :: MatrixVectorMultiply

CONTAINS

  INTEGER FUNCTION itrans_from_char( ctrans )
    CHARACTER, INTENT(in) :: ctrans
#if defined(THORNADO_LA_CUBLAS)
    IF ( ctrans == 'T' ) THEN
      itrans_from_char = CUBLAS_OP_T
    ELSE
      itrans_from_char = CUBLAS_OP_N
    END IF
#elif defined(THORNADO_LA_MAGMA)
    IF ( ctrans == 'T' ) THEN
      itrans_from_char = MagmaTrans
    ELSE
      itrans_from_char = MagmaNoTrans
    END IF
#endif
    RETURN
  END FUNCTION itrans_from_char

  SUBROUTINE MatrixMatrixMultiply( transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc )

    CHARACTER                          :: transa, transb
    INTEGER                            :: m, n, k, lda, ldb, ldc
    REAL(DP)                           :: alpha, beta
    REAL(DP), DIMENSION(lda,*), TARGET :: a
    REAL(DP), DIMENSION(ldb,*), TARGET :: b
    REAL(DP), DIMENSION(ldc,*), TARGET :: c

    INTEGER                            :: ierr
    INTEGER(C_INT)                     :: itransa, itransb
    INTEGER(C_SIZE_T)                  :: sizeof_a, sizeof_b, sizeof_c
    TYPE(C_PTR)                        :: ha, hb, hc
    TYPE(C_PTR)                        :: da, db, dc
    LOGICAL                            :: data_on_device

    data_on_device = .false.
    sizeof_a = m * k * c_sizeof(0.0_DP)
    sizeof_b = k * n * c_sizeof(0.0_DP)
    sizeof_c = m * n * c_sizeof(0.0_DP)

    ha = C_LOC( a )
    hb = C_LOC( b )
    hc = C_LOC( c )

    data_on_device = device_is_present( ha, mydevice, sizeof_a ) &
               .AND. device_is_present( hb, mydevice, sizeof_b ) &
               .AND. device_is_present( hc, mydevice, sizeof_c )

    IF ( data_on_device ) THEN

      itransa = itrans_from_char( transa )
      itransb = itrans_from_char( transb )

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET DATA USE_DEVICE_PTR( a, b, c )
#elif defined(THORNADO_OACC)
      !$ACC HOST_DATA USE_DEVICE( a, b, c )
#endif
      da = C_LOC( a )
      db = C_LOC( b )
      dc = C_LOC( c )
#if defined(THORNADO_OMP_OL)
      !$OMP END TARGET DATA
#elif defined(THORNADO_OACC)
      !$ACC END HOST_DATA
#endif

#if defined(THORNADO_LA_CUBLAS)
      ierr = cublasDgemm_v2 &
             ( cublas_handle, itransa, itransb, m, n, k, alpha, da, lda, db, ldb, beta, dc, ldc )
      ierr = cudaStreamSynchronize( stream )
#elif defined(THORNADO_LA_MAGMA)
      CALL magma_dgemm &
             ( itransa, itransb, m, n, k, alpha, da, lda, db, ldb, beta, dc, ldc, magma_queue )
      CALL magma_queue_sync( magma_queue )
#endif

    ELSE

#if defined(THORNADO_GPU)
      WRITE(*,*) '[MatrixMatrixMultiply] Data not present on device'
#endif

      CALL DGEMM &
             ( transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc )

    END IF

  END SUBROUTINE MatrixMatrixMultiply


  SUBROUTINE MatrixVectorMultiply( trans, m, n, alpha, a, lda, x, incx, beta, y, incy )

    CHARACTER                          :: trans
    INTEGER                            :: m, n, lda, incx, incy
    REAL(DP)                           :: alpha, beta
    REAL(DP), DIMENSION(lda,*), TARGET :: a
    REAL(DP), DIMENSION(*)    , TARGET :: x
    REAL(DP), DIMENSION(*)    , TARGET :: y

    INTEGER                            :: ierr
    INTEGER(C_INT)                     :: itrans
    INTEGER(C_SIZE_T)                  :: sizeof_a, sizeof_x
    TYPE(C_PTR)                        :: ha, hx, hy
    TYPE(C_PTR)                        :: da, dx, dy
    LOGICAL                            :: data_on_device

    data_on_device = .false.
    sizeof_a = m * n * c_sizeof(0.0_DP)
    sizeof_x = m     * c_sizeof(0.0_DP)

    ha = C_LOC( a )
    hx = C_LOC( x )
    hy = C_LOC( y )

    data_on_device = device_is_present( ha, mydevice, sizeof_a ) &
               .AND. device_is_present( hx, mydevice, sizeof_x ) &
               .AND. device_is_present( hy, mydevice, sizeof_x )

    IF ( data_on_device ) THEN

      itrans = itrans_from_char( trans )

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET DATA USE_DEVICE_PTR( a, x, y )
#elif defined(THORNADO_OACC)
      !$ACC HOST_DATA USE_DEVICE( a, x, y )
#endif
      da = C_LOC( a )
      dx = C_LOC( x )
      dy = C_LOC( y )
#if defined(THORNADO_OMP_OL)
      !$OMP END TARGET DATA
#elif defined(THORNADO_OACC)
      !$ACC END HOST_DATA
#endif

#if defined(THORNADO_LA_CUBLAS)
      ierr = cublasDgemv_v2 &
             ( cublas_handle, itrans, m, n, alpha, da, lda, dx, incx, beta, dy, incy )
      ierr = cudaStreamSynchronize( stream )
#elif defined(THORNADO_LA_MAGMA)
      CALL magma_dgemm &
             ( itrans, m, n, alpha, da, lda, dx, incx, beta, dy, incy, magma_queue )
      CALL magma_queue_sync( magma_queue )
#endif

    ELSE

#if defined(THORNADO_GPU)
      WRITE(*,*) '[MatrixMatrixMultiply] Data not present on device'
#endif

      CALL DGEMV &
             ( trans, m, n, alpha, a, lda, x, incx, beta, y, incy )

    END IF

  END SUBROUTINE MatrixVectorMultiply

END MODULE LinearAlgebraModule
