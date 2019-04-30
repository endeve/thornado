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
    REAL(DP), DIMENSION(:,:), POINTER  :: pa, pb, pc
    TYPE(C_PTR)                        :: ha, hb, hc
    TYPE(C_PTR)                        :: da, db, dc
    INTEGER                            :: ka, kb
    LOGICAL                            :: data_on_device

    data_on_device = .false.
    sizeof_a = m * k * c_sizeof(0.0_DP)
    sizeof_b = k * n * c_sizeof(0.0_DP)
    sizeof_c = m * n * c_sizeof(0.0_DP)

    IF ( transa == 'N' ) THEN
      ka = k
    ELSE
      ka = m
    END IF

    IF ( transb == 'N' ) THEN
      kb = n
    ELSE
      kb = k
    END IF

    pa(1:lda,1:ka) => a(:,1:ka)
    pb(1:ldb,1:kb) => b(:,1:kb)
    pc(1:ldc,1:n ) => c(:,1:n )

    ha = C_LOC( pa )
    hb = C_LOC( pb )
    hc = C_LOC( pc )

    data_on_device = device_is_present( ha, mydevice, sizeof_a ) &
               .AND. device_is_present( hb, mydevice, sizeof_b ) &
               .AND. device_is_present( hc, mydevice, sizeof_c )

    IF ( data_on_device ) THEN

      itransa = itrans_from_char( transa )
      itransb = itrans_from_char( transb )

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET DATA USE_DEVICE_PTR( pa, pb, pc )
#elif defined(THORNADO_OACC)
      !$ACC HOST_DATA USE_DEVICE( pa, pb, pc )
#endif
      da = C_LOC( pa )
      db = C_LOC( pb )
      dc = C_LOC( pc )
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
      IF ( .not. device_is_present( ha, mydevice, sizeof_a ) ) &
        WRITE(*,*) '[MatrixMatrixMultiply]   A missing'
      IF ( .not. device_is_present( hb, mydevice, sizeof_b ) ) &
        WRITE(*,*) '[MatrixMatrixMultiply]   B missing'
      IF ( .not. device_is_present( hc, mydevice, sizeof_c ) ) &
        WRITE(*,*) '[MatrixMatrixMultiply]   C missing'
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
    INTEGER(C_SIZE_T)                  :: sizeof_a, sizeof_x, sizeof_y
    REAL(DP), DIMENSION(:,:), POINTER  :: pa
    REAL(DP), DIMENSION(:)  , POINTER  :: px, py
    TYPE(C_PTR)                        :: ha, hx, hy
    TYPE(C_PTR)                        :: da, dx, dy
    LOGICAL                            :: data_on_device

    data_on_device = .false.
    sizeof_a = m * n * c_sizeof(0.0_DP)
    sizeof_x = m     * c_sizeof(0.0_DP)
    sizeof_y =     n * c_sizeof(0.0_DP)

    pa(1:lda,1:n) => a(:,1:n)
    px(1:m) => x(1:m)
    py(1:n) => y(1:n)

    ha = C_LOC( pa )
    hx = C_LOC( px )
    hy = C_LOC( py )

    data_on_device = device_is_present( ha, mydevice, sizeof_a ) &
               .AND. device_is_present( hx, mydevice, sizeof_x ) &
               .AND. device_is_present( hy, mydevice, sizeof_y )

    IF ( data_on_device ) THEN

      itrans = itrans_from_char( trans )

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET DATA USE_DEVICE_PTR( pa, px, py )
#elif defined(THORNADO_OACC)
      !$ACC HOST_DATA USE_DEVICE( pa, px, py )
#endif
      da = C_LOC( pa )
      dx = C_LOC( px )
      dy = C_LOC( py )
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
      CALL magma_dgemv &
             ( itrans, m, n, alpha, da, lda, dx, incx, beta, dy, incy, magma_queue )
      CALL magma_queue_sync( magma_queue )
#endif

    ELSE

#if defined(THORNADO_GPU)
      WRITE(*,*) '[MatrixVectorMultiply] Data not present on device'
      IF ( .not. device_is_present( ha, mydevice, sizeof_a ) ) &
        WRITE(*,*) '[MatrixVectorMultiply]   A missing'
      IF ( .not. device_is_present( hx, mydevice, sizeof_x ) ) &
        WRITE(*,*) '[MatrixVectorMultiply]   x missing'
      IF ( .not. device_is_present( hy, mydevice, sizeof_y ) ) &
        WRITE(*,*) '[MatrixVectorMultiply]   y missing'
#endif

      CALL DGEMV &
             ( trans, m, n, alpha, a, lda, x, incx, beta, y, incy )

    END IF

  END SUBROUTINE MatrixVectorMultiply

END MODULE LinearAlgebraModule
