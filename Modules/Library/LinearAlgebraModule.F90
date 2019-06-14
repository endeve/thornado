MODULE LinearAlgebraModule

  USE, INTRINSIC :: ISO_C_BINDING
  USE KindModule, ONLY: &
    DP, &
    Zero, &
    One
  USE DeviceModule, ONLY: &
    mydevice, &
    device_is_present

#if defined(THORNADO_LA_CUBLAS)
  USE CudaModule, ONLY: &
    stream, &
    cudaStreamSynchronize
  USE CublasModule, ONLY: &
    cublas_handle, &
    cublasDnrm2_v2, &
    cublasDgemm_v2, &
    cublasDgemmStridedBatched, &
    cublasDgemv_v2, &
    cublasDtrsv_v2, &
    cublasDtrsm_v2, &
    CUBLAS_OP_N, CUBLAS_OP_T, &
    CUBLAS_SIDE_LEFT, &
    CUBLAS_FILL_MODE_UPPER, &
    CUBLAS_DIAG_NON_UNIT
  USE CusolverModule, ONLY: &
    cusolver_handle, &
    cusolverDnDgeqrf_bufferSize, &
    cusolverDnDgeqrf, &
    cusolverDnDormqr
#endif

#if defined(THORNADO_LA_MAGMA)
  USE MagmaModule, ONLY: &
    magma_queue, &
    magma_queue_sync, &
    magma_dgemm, &
    magmablas_dgemm_batched_strided, &
    magma_dgemv, &
    magma_dgels_gpu, &
    MagmaNoTrans, MagmaTrans
#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MatrixMatrixMultiply
  PUBLIC :: MatrixMatrixMultiplyBatched
  PUBLIC :: MatrixVectorMultiply
  PUBLIC :: VectorNorm2
  PUBLIC :: VectorNorm2_Kernel
  PUBLIC :: LinearLeastSquares_LWORK
  PUBLIC :: LinearLeastSquares

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


  SUBROUTINE MatrixMatrixMultiplyBatched( transa, transb, m, n, k, alpha, a, lda, stridea, &
                                          b, ldb, strideb, beta, c, ldc, stridec, batchcount )

    CHARACTER                          :: transa, transb
    INTEGER                            :: m, n, k, lda, ldb, ldc, stridea, strideb, stridec, batchcount
    REAL(DP)                           :: alpha, beta
    REAL(DP), DIMENSION(lda,*), TARGET :: a
    REAL(DP), DIMENSION(ldb,*), TARGET :: b
    REAL(DP), DIMENSION(ldc,*), TARGET :: c

    INTEGER                            :: ierr, i
    INTEGER(C_INT)                     :: itransa, itransb
    INTEGER(C_SIZE_T)                  :: sizeof_a, sizeof_b, sizeof_c
    REAL(DP), DIMENSION(:,:), POINTER  :: pa, pb, pc
    TYPE(C_PTR)                        :: ha, hb, hc
    TYPE(C_PTR)                        :: da, db, dc
    INTEGER                            :: ka, kb
    INTEGER                            :: osa, osb, osc
    LOGICAL                            :: data_on_device

    data_on_device = .false.
    sizeof_a = m * k * batchcount * c_sizeof(0.0_DP)
    sizeof_b = k * n * batchcount * c_sizeof(0.0_DP)
    sizeof_c = m * n * batchcount * c_sizeof(0.0_DP)

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

    pa(1:lda,1:ka*batchcount) => a(:,1:ka*batchcount)
    pb(1:ldb,1:kb*batchcount) => b(:,1:kb*batchcount)
    pc(1:ldc,1:n *batchcount) => c(:,1:n *batchcount)

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
      ierr = cublasDgemmStridedBatched &
             ( cublas_handle, itransa, itransb, m, n, k, alpha, da, lda, stridea, &
               db, ldb, strideb, beta, dc, ldc, stridec, batchcount )
      ierr = cudaStreamSynchronize( stream )
#elif defined(THORNADO_LA_MAGMA)
      CALL magmablas_dgemm_batched_strided &
             ( itransa, itransb, m, n, k, alpha, da, lda, stridea, &
               db, ldb, strideb, beta, dc, ldc, stridec, batchcount, magma_queue )
      CALL magma_queue_sync( magma_queue )
#endif

    ELSE

#if defined(THORNADO_GPU)
      WRITE(*,*) '[MatrixMatrixMultiplyBatched] Data not present on device'
      IF ( .not. device_is_present( ha, mydevice, sizeof_a ) ) &
        WRITE(*,*) '[MatrixMatrixMultiplyBatched]   A missing'
      IF ( .not. device_is_present( hb, mydevice, sizeof_b ) ) &
        WRITE(*,*) '[MatrixMatrixMultiplyBatched]   B missing'
      IF ( .not. device_is_present( hc, mydevice, sizeof_c ) ) &
        WRITE(*,*) '[MatrixMatrixMultiplyBatched]   C missing'
#endif

      DO i = 1, batchcount
        osa = (i-1) * ka + 1
        osb = (i-1) * kb + 1
        osc = (i-1) * n  + 1
        CALL DGEMM &
               ( transa, transb, m, n, k, alpha, a(1,osa), lda, b(1,osb), ldb, beta, c(1,osc), ldc )
      END DO

    END IF

  END SUBROUTINE MatrixMatrixMultiplyBatched


  SUBROUTINE MatrixVectorMultiply( trans, m, n, alpha, a, lda, x, incx, beta, y, incy )

    CHARACTER                          :: trans
    INTEGER                            :: m, n, lda, incx, incy
    REAL(DP)                           :: alpha, beta
    REAL(DP), DIMENSION(lda,*), TARGET :: a
    REAL(DP), DIMENSION(*)    , TARGET :: x
    REAL(DP), DIMENSION(*)    , TARGET :: y

    INTEGER                            :: ierr, lenx, leny
    INTEGER(C_INT)                     :: itrans
    INTEGER(C_SIZE_T)                  :: sizeof_a, sizeof_x, sizeof_y
    REAL(DP), DIMENSION(:,:), POINTER  :: pa
    REAL(DP), DIMENSION(:)  , POINTER  :: px, py
    TYPE(C_PTR)                        :: ha, hx, hy
    TYPE(C_PTR)                        :: da, dx, dy
    LOGICAL                            :: data_on_device

    data_on_device = .false.

    IF ( trans == 'T' ) THEN
      lenx = m
      leny = n
    ELSE
      lenx = n
      leny = m
    END IF

    sizeof_a = m * n * c_sizeof(0.0_DP)
    sizeof_x = lenx * c_sizeof(0.0_DP)
    sizeof_y = leny * c_sizeof(0.0_DP)

    pa(1:lda,1:n) => a(:,1:n)
    px(1:lenx) => x(1:lenx)
    py(1:leny) => y(1:leny)

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


  SUBROUTINE LinearLeastSquares_LWORK( trans, m, n, nrhs, a, lda, b, ldb, work, lwork )

    CHARACTER                          :: trans
    INTEGER                            :: m, n, nrhs, lda, ldb, lwork
    REAL(DP), DIMENSION(lda,*), TARGET :: a
    REAL(DP), DIMENSION(ldb,*), TARGET :: b
    REAL(DP), DIMENSION(*)    , TARGET :: work

    INTEGER                            :: ierr, info, max_mn
    INTEGER(C_INT)                     :: itrans
    INTEGER(C_SIZE_T)                  :: sizeof_a, sizeof_b, sizeof_work
    REAL(DP), DIMENSION(:,:), POINTER  :: pa, pb
    REAL(DP), DIMENSION(:)  , POINTER  :: pwork
    TYPE(C_PTR)                        :: ha, hb, hwork
    TYPE(C_PTR)                        :: da, db

    lwork = -1

    max_mn = MAX(m,n)

    sizeof_a = m * n * c_sizeof(0.0_DP)
    sizeof_b = max_mn * nrhs * c_sizeof(0.0_DP)
    sizeof_work = lwork * c_sizeof(0.0_DP)

    pa(1:lda,1:n) => a(:,1:n)
    pb(1:ldb,1:nrhs) => b(:,1:nrhs)
    pwork(1:lwork) => work(1:lwork)

    ha = C_LOC( pa )
    hb = C_LOC( pb )
    hwork = C_LOC( pwork )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET DATA USE_DEVICE_PTR( pa, pb )
#elif defined(THORNADO_OACC)
    !$ACC HOST_DATA USE_DEVICE( pa, pb )
#endif
    da = C_LOC( pa )
    db = C_LOC( pb )
#if defined(THORNADO_OMP_OL)
    !$OMP END TARGET DATA
#elif defined(THORNADO_OACC)
    !$ACC END HOST_DATA
#endif

    itrans = itrans_from_char( trans )

#if defined(THORNADO_LA_CUBLAS)
    ierr = cusolverDnDgeqrf_bufferSize &
           ( cusolver_handle, m, n, da, lda, lwork )
#elif defined(THORNADO_LA_MAGMA)
    CALL magma_dgels_gpu &
           ( itrans, m, n, nrhs, da, lda, db, ldb, hwork, lwork, info, magma_queue )
    lwork = INT( work(1) )
#else
    CALL DGELS &
           ( trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info )
    lwork = INT( work(1) )
#endif

  END SUBROUTINE LinearLeastSquares_LWORK


  SUBROUTINE LinearLeastSquares( trans, m, n, nrhs, a, lda, b, ldb, tau, work, lwork, info )

    CHARACTER                          :: trans
    INTEGER                            :: m, n, nrhs, lda, ldb, lwork
    INTEGER                   , TARGET :: info
    REAL(DP), DIMENSION(lda,*), TARGET :: a
    REAL(DP), DIMENSION(ldb,*), TARGET :: b
    REAL(DP), DIMENSION(*)    , TARGET :: tau
    REAL(DP), DIMENSION(*)    , TARGET :: work

    INTEGER                            :: ierr, max_mn, min_mn
    INTEGER(C_INT)                     :: itrans
    INTEGER(C_SIZE_T)                  :: sizeof_a, sizeof_b, sizeof_tau, sizeof_work, sizeof_info
    REAL(DP), DIMENSION(:,:), POINTER  :: pa, pb
    REAL(DP), DIMENSION(:)  , POINTER  :: ptau, pwork
    INTEGER                 , POINTER  :: pinfo
    TYPE(C_PTR)                        :: ha, hb, htau, hwork, hinfo
    TYPE(C_PTR)                        :: da, db, dtau, dwork, dinfo
    LOGICAL                            :: data_on_device

    data_on_device = .false.
    max_mn = MAX(m,n)
    min_mn = MIN(m,n)

    sizeof_a = m * n * c_sizeof(0.0_DP)
    sizeof_b = max_mn * nrhs * c_sizeof(0.0_DP)
    sizeof_tau = min_mn * c_sizeof(0.0_DP)
    sizeof_work = lwork * c_sizeof(0.0_DP)
    sizeof_info = c_sizeof(info)

    pa(1:lda,1:n) => a(:,1:n)
    pb(1:ldb,1:nrhs) => b(:,1:nrhs)
    ptau(1:min_mn) => tau(1:min_mn)
    pwork(1:lwork) => work(1:lwork)
    pinfo => info

    ha = C_LOC( pa )
    hb = C_LOC( pb )
    htau = C_LOC( ptau )
    hwork = C_LOC( pwork )
    hinfo = C_LOC( pinfo )

    data_on_device = device_is_present( ha,    mydevice, sizeof_a    ) &
               .AND. device_is_present( hb,    mydevice, sizeof_b    ) &
               .AND. device_is_present( htau,  mydevice, sizeof_tau  ) &
               .AND. device_is_present( hwork, mydevice, sizeof_work ) &
               .AND. device_is_present( hinfo, mydevice, sizeof_info )

    IF ( data_on_device ) THEN

      itrans = itrans_from_char( trans )

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET DATA USE_DEVICE_PTR( pa, pb, ptau, pwork, pinfo )
#elif defined(THORNADO_OACC)
      !$ACC HOST_DATA USE_DEVICE( pa, pb, ptau, pwork, pinfo )
#endif
      da = C_LOC( pa )
      db = C_LOC( pb )
      dtau = C_LOC( ptau )
      dwork = C_LOC( pwork )
      dinfo = C_LOC( pinfo )
#if defined(THORNADO_OMP_OL)
      !$OMP END TARGET DATA
#elif defined(THORNADO_OACC)
      !$ACC END HOST_DATA
#endif

#if defined(THORNADO_LA_CUBLAS)
      ierr = cusolverDnDgeqrf &
             ( cusolver_handle, m, n, da, lda, dtau, dwork, lwork, dinfo )
      ierr = cusolverDnDormqr &
             ( cusolver_handle, &
               CUBLAS_SIDE_LEFT, CUBLAS_OP_T, &
               m, nrhs, n, da, lda, dtau, db, ldb, dwork, lwork, dinfo )
      ierr = cudaStreamSynchronize( stream )

      IF ( nrhs == 1 ) THEN

        ierr = cublasDtrsv_v2 &
               ( cublas_handle, &
                 CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, &
                 n, da, lda, db, 1 )

      ELSE

        ierr = cublasDtrsm_v2 &
               ( cublas_handle, &
                 CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, &
                 n, nrhs, One, da, lda, db, ldb )

      END IF
      ierr = cudaStreamSynchronize( stream )
#elif defined(THORNADO_LA_MAGMA)
      CALL magma_dgels_gpu &
             ( itrans, m, n, nrhs, da, lda, db, ldb, hwork, lwork, info, magma_queue )
      CALL magma_queue_sync( magma_queue )
#endif

    ELSE

#if defined(THORNADO_GPU)
      WRITE(*,*) '[LinearLeastSquares] Data not present on device'
      IF ( .not. device_is_present( ha, mydevice, sizeof_a ) ) &
        WRITE(*,*) '[LinearLeastSquares]   A missing'
      IF ( .not. device_is_present( hb, mydevice, sizeof_b ) ) &
        WRITE(*,*) '[LinearLeastSquares]   b missing'
      IF ( .not. device_is_present( htau, mydevice, sizeof_tau ) ) &
        WRITE(*,*) '[LinearLeastSquares]   tau missing'
      IF ( .not. device_is_present( hwork, mydevice, sizeof_work ) ) &
        WRITE(*,*) '[LinearLeastSquares]   work missing'
      IF ( .not. device_is_present( hinfo, mydevice, sizeof_info ) ) &
        WRITE(*,*) '[LinearLeastSquares]   info missing'
#endif

      CALL DGELS &
             ( trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info )

    END IF

  END SUBROUTINE LinearLeastSquares


  SUBROUTINE VectorNorm2( n, x, incx, xnorm )

    INTEGER                         :: n, incx
    REAL(DP), DIMENSION(*), TARGET  :: x
    REAL(DP)                        :: xnorm

    INTEGER                         :: ierr
    INTEGER(C_SIZE_T)               :: sizeof_x
    REAL(DP), DIMENSION(:), POINTER :: px
    TYPE(C_PTR)                     :: hx
    TYPE(C_PTR)                     :: dx
    LOGICAL                         :: data_on_device
    REAL(DP), EXTERNAL              :: DNRM2

    data_on_device = .false.
    sizeof_x = n * c_sizeof(0.0_DP)

    px(1:n) => x(1:n)

    hx = C_LOC( px )

    data_on_device = device_is_present( hx, mydevice, sizeof_x )

    IF ( data_on_device ) THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET DATA USE_DEVICE_PTR( px )
#elif defined(THORNADO_OACC)
      !$ACC HOST_DATA USE_DEVICE( px )
#endif
      dx = C_LOC( px )
#if defined(THORNADO_OMP_OL)
      !$OMP END TARGET DATA
#elif defined(THORNADO_OACC)
      !$ACC END HOST_DATA
#endif

#if defined(THORNADO_LA_CUBLAS)
      ierr = cublasDnrm2_v2( cublas_handle, n, dx, incx, xnorm )
      ierr = cudaStreamSynchronize( stream )
#elif defined(THORNADO_LA_MAGMA)
      xnorm = magma_dnrm2( n, dx, incx, magma_queue )
      CALL magma_queue_sync( magma_queue )
#endif

    ELSE

#if defined(THORNADO_GPU)
      WRITE(*,*) '[VectorNorm2] Data not present on device'
      IF ( .not. device_is_present( hx, mydevice, sizeof_x ) ) &
        WRITE(*,*) '[VectorNorm2]   x missing'
#endif

      xnorm = DNRM2( n, x, incx )

    END IF

  END SUBROUTINE VectorNorm2


  SUBROUTINE VectorNorm2_Kernel( n, x, incx, xnorm )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    INTEGER                         :: n, incx
    REAL(DP), DIMENSION(*), TARGET  :: x
    REAL(DP)                        :: xnorm

    INTEGER                         :: ix
    REAL(DP)                        :: xscale, xssq, absxi

    IF ( n < 1 .OR. incx < 1 ) THEN
      xnorm = 0.0d0
    ELSE IF ( n == 1 ) THEN
      xnorm = ABS( x(1) )
    ELSE
      xscale = 0.0d0
      xssq = 1.0d0
      DO ix = 1, 1 + (n-1)*incx, incx
        IF ( x(ix) /= 0.0d0 ) THEN
          absxi = ABS( x(ix) )
          IF ( xscale < absxi ) THEN
            xssq = 1.0d0 + xssq * (xscale/absxi)**2
            xscale = absxi
          ELSE
            xssq = xssq + (absxi/xscale)**2
          END IF
        END IF
      END DO
      xnorm = xscale * SQRT(xssq)
    END IF

  END SUBROUTINE VectorNorm2_Kernel


END MODULE LinearAlgebraModule
