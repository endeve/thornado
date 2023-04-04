#ifdef THORNADO_DEBUG
#define THORNADO_DEBUG_LA
#endif
#if defined(THORNADO_LA_ONEMKL)
include "mkl_omp_offload.f90"
#endif
MODULE LinearAlgebraModule

  USE, INTRINSIC :: ISO_C_BINDING
  USE KindModule, ONLY: &
    DP, &
    Zero, &
    One, &
    Pi
  USE DeviceModule, ONLY: &
    mydevice, &
    device_is_present, &
    dev_ptr, &
    stream_sync, &
    stream

#if defined(THORNADO_LA_CUBLAS)
  USE CublasModule, ONLY: &
    cublas_handle, &
    cublasDnrm2_v2, &
    cublasDaxpy_v2, &
    cublasDgemm_v2, &
    cublasDgemmStridedBatched, &
    cublasDgetrfBatched, &
    cublasDgetrsBatched, &
    cublasDgemv_v2, &
    cublasDtrsv_v2, &
    cublasDtrsm_v2, &
    cublasDgeam, &
    cublasDdgmm, &
    CUBLAS_OP_N, CUBLAS_OP_T, &
    CUBLAS_SIDE_LEFT, &
    CUBLAS_FILL_MODE_UPPER, &
    CUBLAS_DIAG_NON_UNIT
  USE CusolverModule, ONLY: &
    cusolver_handle, &
    cusolverDnDgeqrf_bufferSize, &
    cusolverDnDgeqrf, &
    cusolverDnDormqr
  USE CusparseModule, ONLY: &
    cusparse_handle, &
    cusparseDgthr, &
    CUSPARSE_INDEX_BASE_ONE
#elif defined(THORNADO_LA_ROCM)
  USE HipModule, ONLY: &
    hipCheck, &
    hipblasCheck, &
    hipsparseCheck, &
    rocblasCheck, &
    rocsparseCheck, &
    rocsolverCheck
  USE RocblasModule, ONLY: &
    rocblas_handle, &
    rocblas_dnrm2, &
    rocblas_daxpy, &
    rocblas_dgemm, &
    rocblas_dgemm_strided_batched, &
    rocblas_dgemv, &
    rocblas_dtrsv, &
    rocblas_dtrsm, &
    rocblas_dgeam, &
    rocblas_ddgmm, &
    rocblas_ddot_strided_batched, &
    rocblas_operation_none, &
    rocblas_operation_transpose, &
    rocblas_side_left, &
    rocblas_fill_upper, &
    rocblas_diagonal_non_unit
  USE RocsolverModule, ONLY: &
    rocsolver_handle, &
    rocsolver_dgeqrf, &
    rocsolver_dormqr, &
    rocsolver_dgetrf_batched, &
    rocsolver_dgetrs_batched
  USE RocsparseModule, ONLY: &
    rocsparse_handle, &
    rocsparse_dgthr, &
    rocsparse_index_base_one
  USE HipblasModule, ONLY: &
    hipblas_handle, &
    hipblasDnrm2, &
    hipblasDaxpy, &
    hipblasDgemm, &
    hipblasDgemmStridedBatched, &
    hipblasDgetrfBatched, &
    hipblasDgetrsBatched, &
    hipblasDgemv, &
    hipblasDtrsv, &
    hipblasDtrsm, &
    hipblasDgeam, &
    hipblasDdgmm, &
    HIPBLAS_OP_N, HIPBLAS_OP_T, &
    HIPBLAS_SIDE_LEFT, &
    HIPBLAS_FILL_MODE_UPPER, &
    HIPBLAS_DIAG_NON_UNIT
  USE HipsparseModule, ONLY: &
    hipsparse_handle, &
    hipsparseDgthr, &
    HIPSPARSE_INDEX_BASE_ONE
#elif defined(THORNADO_LA_ONEMKL)
  USE onemkl_blas_omp_offload_lp64
#elif defined(THORNADO_LA_MAGMA)
  USE MagmaModule, ONLY: &
    magma_queue, &
    magma_dnrm2, &
    magma_daxpy, &
    magma_dgemm, &
    magmablas_dgemm_batched_strided, &
    magma_dgetrf_batched, &
    magma_dgetrs_batched, &
    magma_dgemv, &
    magma_dgels_gpu, &
    magmablas_dtranspose, &
    magmablas_dlacpy, &
    magmablas_dlascl2, &
    magmablas_dgeadd2, &
    magma_dmalloc, &
    MagmaNoTrans, MagmaTrans, &
    MagmaFull
#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MatrixMatrixAdd
  PUBLIC :: MatrixMatrixMultiply
  PUBLIC :: MatrixMatrixMultiplyBatched
  PUBLIC :: MatrixVectorMultiply
  PUBLIC :: MatrixDiagScale
  PUBLIC :: VectorDotProductBatched
  PUBLIC :: VectorNorm2
  PUBLIC :: VectorNorm2_Kernel
  PUBLIC :: VectorVectorAdd
  PUBLIC :: LinearLeastSquares_LWORK
  PUBLIC :: LinearLeastSquares
  PUBLIC :: LinearSolveBatched
  PUBLIC :: EigenvaluesSymmetric3

CONTAINS

  INTEGER FUNCTION itrans_from_char( ctrans )
    CHARACTER, INTENT(in) :: ctrans
    itrans_from_char = 0
#if defined(THORNADO_LA_CUBLAS)
    IF ( ctrans == 'T' ) THEN
      itrans_from_char = CUBLAS_OP_T
    ELSE
      itrans_from_char = CUBLAS_OP_N
    END IF
#elif defined(THORNADO_LA_ROCM)
    IF ( ctrans == 'T' ) THEN
      !itrans_from_char = rocblas_operation_transpose
      itrans_from_char = HIPBLAS_OP_T
    ELSE
      !itrans_from_char = rocblas_operation_none
      itrans_from_char = HIPBLAS_OP_N
    END IF
#elif defined(THORNADO_LA_ONEMKL)
#elif defined(THORNADO_LA_MAGMA)
    IF ( ctrans == 'T' ) THEN
      itrans_from_char = MagmaTrans
    ELSE
      itrans_from_char = MagmaNoTrans
    END IF
#endif
    RETURN
  END FUNCTION itrans_from_char


  SUBROUTINE MatrixMatrixAdd( transa, transb, m, n, alpha, a, lda, beta, b, ldb, c, ldc )

    CHARACTER                          :: transa, transb
    INTEGER                            :: m, n, lda, ldb, ldc
    REAL(DP)                           :: alpha, beta
    REAL(DP), DIMENSION(lda,*), TARGET :: a
    REAL(DP), DIMENSION(ldb,*), TARGET :: b
    REAL(DP), DIMENSION(ldc,*), TARGET :: c

    INTEGER                            :: i, j, ierr, info
    INTEGER(C_INT)                     :: itransa, itransb
    INTEGER(C_SIZE_T)                  :: sizeof_a, sizeof_b, sizeof_c, mn
    REAL(DP), DIMENSION(:,:), POINTER  :: pa, pb, pc
    TYPE(C_PTR)                        :: ha, hb, hc
    TYPE(C_PTR)                        :: da, db, dc
    TYPE(C_PTR)                        :: dat, dbt
    INTEGER                            :: ka, kb
    LOGICAL                            :: data_on_device

    data_on_device = .false.
    sizeof_a = m * n * c_sizeof(0.0_DP)
    sizeof_b = m * n * c_sizeof(0.0_DP)
    sizeof_c = m * n * c_sizeof(0.0_DP)
    mn = m * n

    IF ( transa == 'N' ) THEN
      ka = n
    ELSE
      ka = m
    END IF

    IF ( transb == 'N' ) THEN
      kb = n
    ELSE
      kb = m
    END IF

    pa => a(:,1:ka)
    pb => b(:,1:kb)
    pc => c(:,1:n )

    ha = C_LOC( pa )
    hb = C_LOC( pb )
    hc = C_LOC( pc )

    data_on_device = device_is_present( ha, mydevice, sizeof_a ) &
               .AND. device_is_present( hb, mydevice, sizeof_b ) &
               .AND. device_is_present( hc, mydevice, sizeof_c )

    IF ( data_on_device ) THEN

      itransa = itrans_from_char( transa )
      itransb = itrans_from_char( transb )

      da = dev_ptr( pa(1,1) )
      db = dev_ptr( pb(1,1) )
      dc = dev_ptr( pc(1,1) )

#if defined(THORNADO_LA_CUBLAS)
      ierr = cublasDgeam &
             ( cublas_handle, itransa, itransb, m, n, alpha, da, lda, beta, db, ldb, dc, ldc )
#elif defined(THORNADO_LA_ROCM)
      !CALL rocblasCheck( rocblas_dgeam &
      !       ( rocblas_handle, itransa, itransb, m, n, alpha, da, lda, beta, db, ldb, dc, ldc ) )
      CALL hipblasCheck( hipblasDgeam &
             ( hipblas_handle, itransa, itransb, m, n, alpha, da, lda, beta, db, ldb, dc, ldc ) )
#elif defined(THORNADO_LA_ONEMKL)
      !!$OMP TARGET VARIANT DISPATCH USE_DEVICE_PTR( a, b, c )
      !CALL DGEAM &
      !       ( transa, transb, m, n, alpha, a, lda, beta, b, ldb, c, ldc )
      !!$OMP END TARGET VARIANT DISPATCH
#elif defined(THORNADO_LA_MAGMA)
      IF ( transb  == 'N' ) THEN
        CALL magmablas_dlacpy &
               ( MagmaFull, m, n, db, ldb, dc, ldc, magma_queue )
      ELSE
        CALL magmablas_dtranspose &
               ( n, m, db, ldb, dc, ldc, magma_queue )
      END IF
      IF ( transa == 'N' ) THEN
        CALL magmablas_dgeadd2 &
               ( m, n, alpha, da, lda, beta, dc, ldc, magma_queue )
      ELSE
        ierr = magma_dmalloc( dat, mn )
        CALL magmablas_dtranspose &
               ( n, m, da, lda, dat, m, magma_queue )
        CALL magmablas_dgeadd2 &
               ( m, n, alpha, dat, m, beta, dc, ldc, magma_queue )
      END IF
#endif
#if defined(THORNADO_OMP_OL)
      CALL stream_sync( stream )
#endif

    ELSE

#if defined(THORNADO_DEBUG_LA)
#if defined(THORNADO_GPU)
      WRITE(*,*) '[MatrixMatrixAdd] Data not present on device'
      IF ( .not. device_is_present( ha, mydevice, sizeof_a ) ) &
        WRITE(*,*) '[MatrixMatrixAdd]   A missing'
      IF ( .not. device_is_present( hb, mydevice, sizeof_b ) ) &
        WRITE(*,*) '[MatrixMatrixAdd]   B missing'
      IF ( .not. device_is_present( hc, mydevice, sizeof_c ) ) &
        WRITE(*,*) '[MatrixMatrixAdd]   C missing'
#endif
#endif

      IF ( alpha == 0.0_DP .AND. beta == 0.0_DP ) THEN

        DO j = 1, n
          DO i = 1, m
            c(i,j) = 0.0_DP
          END DO
        END DO

      ELSE IF ( alpha == 0.0_DP ) THEN

        IF ( transb == 'N' ) THEN
          DO j = 1, n
            DO i = 1, m
              c(i,j) = + beta * b(i,j)
            END DO
          END DO
        ELSE
          DO j = 1, n
            DO i = 1, m
              c(i,j) = + beta * b(j,i)
            END DO
          END DO
        END IF

      ELSE IF ( beta == 0.0_DP ) THEN

        IF ( transa == 'N' ) THEN
          DO j = 1, n
            DO i = 1, m
              c(i,j) = alpha * a(i,j)
            END DO
          END DO
        ELSE
          DO j = 1, n
            DO i = 1, m
              c(i,j) = alpha * a(j,i)
            END DO
          END DO
        END IF

      ELSE

        IF ( transa == 'N' ) THEN
          IF ( transb == 'N' ) THEN
            DO j = 1, n
              DO i = 1, m
                c(i,j) = alpha * a(i,j) + beta * b(i,j)
              END DO
            END DO
          ELSE
            DO j = 1, n
              DO i = 1, m
                c(i,j) = alpha * a(i,j) + beta * b(j,i)
              END DO
            END DO
          END IF
        ELSE
          IF ( transb == 'N' ) THEN
            DO j = 1, n
              DO i = 1, m
                c(i,j) = alpha * a(j,i) + beta * b(i,j)
              END DO
            END DO
          ELSE
            DO j = 1, n
              DO i = 1, m
                c(i,j) = alpha * a(j,i) + beta * b(j,i)
              END DO
            END DO
          END IF
        END IF

      END IF

    END IF

  END SUBROUTINE MatrixMatrixAdd


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

    pa => a(:,1:ka)
    pb => b(:,1:kb)
    pc => c(:,1:n )

    ha = C_LOC( pa )
    hb = C_LOC( pb )
    hc = C_LOC( pc )

    data_on_device = device_is_present( ha, mydevice, sizeof_a ) &
               .AND. device_is_present( hb, mydevice, sizeof_b ) &
               .AND. device_is_present( hc, mydevice, sizeof_c )

    IF ( data_on_device ) THEN

      itransa = itrans_from_char( transa )
      itransb = itrans_from_char( transb )

      da = dev_ptr( pa(1,1) )
      db = dev_ptr( pb(1,1) )
      dc = dev_ptr( pc(1,1) )

#if defined(THORNADO_LA_CUBLAS)
      ierr = cublasDgemm_v2 &
             ( cublas_handle, itransa, itransb, m, n, k, alpha, da, lda, db, ldb, beta, dc, ldc )
#elif defined(THORNADO_LA_ROCM)
      !CALL rocblasCheck( rocblas_dgemm &
      !       ( rocblas_handle, itransa, itransb, m, n, k, alpha, da, lda, db, ldb, beta, dc, ldc ) )
      CALL hipblasCheck( hipblasDgemm &
             ( hipblas_handle, itransa, itransb, m, n, k, alpha, da, lda, db, ldb, beta, dc, ldc ) )
#elif defined(THORNADO_LA_ONEMKL)
      !$OMP TARGET VARIANT DISPATCH USE_DEVICE_PTR( a, b, c )
      CALL DGEMM &
             ( transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc )
      !$OMP END TARGET VARIANT DISPATCH
#elif defined(THORNADO_LA_MAGMA)
      CALL magma_dgemm &
             ( itransa, itransb, m, n, k, alpha, da, lda, db, ldb, beta, dc, ldc, magma_queue )
#endif
#if defined(THORNADO_OMP_OL)
      CALL stream_sync( stream )
#endif

    ELSE

#if defined(THORNADO_DEBUG_LA)
#if defined(THORNADO_GPU)
      WRITE(*,*) '[MatrixMatrixMultiply] Data not present on device'
      IF ( .not. device_is_present( ha, mydevice, sizeof_a ) ) &
        WRITE(*,*) '[MatrixMatrixMultiply]   A missing'
      IF ( .not. device_is_present( hb, mydevice, sizeof_b ) ) &
        WRITE(*,*) '[MatrixMatrixMultiply]   B missing'
      IF ( .not. device_is_present( hc, mydevice, sizeof_c ) ) &
        WRITE(*,*) '[MatrixMatrixMultiply]   C missing'
#endif
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
    INTEGER(C_INT64_T)                 :: stridea_64, strideb_64, stridec_64
    INTEGER(C_LONG_LONG)               :: stridea_l, strideb_l, stridec_l
    INTEGER(C_SIZE_T)                  :: sizeof_a, sizeof_b, sizeof_c
    REAL(DP), DIMENSION(:,:), POINTER  :: pa, pb, pc
    TYPE(C_PTR)                        :: ha, hb, hc
    TYPE(C_PTR)                        :: da, db, dc
    INTEGER                            :: ka, kb
    INTEGER                            :: osa, osb, osc
    INTEGER                            :: ia, ib, ic
    INTEGER                            :: ja, jb, jc
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

    pa => a(:,1:ka*batchcount)
    pb => b(:,1:kb*batchcount)
    pc => c(:,1:n *batchcount)

    ha = C_LOC( pa )
    hb = C_LOC( pb )
    hc = C_LOC( pc )

    data_on_device = device_is_present( ha, mydevice, sizeof_a ) &
               .AND. device_is_present( hb, mydevice, sizeof_b ) &
               .AND. device_is_present( hc, mydevice, sizeof_c )

    IF ( data_on_device ) THEN

      itransa = itrans_from_char( transa )
      itransb = itrans_from_char( transb )

      da = dev_ptr( pa(1,1) )
      db = dev_ptr( pb(1,1) )
      dc = dev_ptr( pc(1,1) )

#if defined(THORNADO_LA_CUBLAS)
      ierr = cublasDgemmStridedBatched &
             ( cublas_handle, itransa, itransb, m, n, k, alpha, da, lda, stridea, &
               db, ldb, strideb, beta, dc, ldc, stridec, batchcount )
#elif defined(THORNADO_LA_ROCM)
      !stridea_64 = stridea
      !strideb_64 = strideb
      !stridec_64 = stridec
      !CALL rocblasCheck( rocblas_dgemm_strided_batched &
      !       ( rocblas_handle, itransa, itransb, m, n, k, alpha, da, lda, stridea_64, &
      !         db, ldb, strideb_64, beta, dc, ldc, stridec_64, batchcount ) )
      stridea_l = stridea
      strideb_l = strideb
      stridec_l = stridec
      CALL hipblasCheck( hipblasDgemmStridedBatched &
             ( hipblas_handle, itransa, itransb, m, n, k, alpha, da, lda, stridea_l, &
               db, ldb, strideb_l, beta, dc, ldc, stridec_l, batchcount ) )
#elif defined(THORNADO_LA_ONEMKL)
      !$OMP TARGET VARIANT DISPATCH USE_DEVICE_PTR( a, b, c )
      CALL DGEMM_BATCH_STRIDED &
             ( transa, transb, m, n, k, alpha, a, lda, stridea, b, ldb, strideb, beta, c, ldc, stridec, batchcount )
      !$OMP END TARGET VARIANT DISPATCH
#elif defined(THORNADO_LA_MAGMA)
      CALL magmablas_dgemm_batched_strided &
             ( itransa, itransb, m, n, k, alpha, da, lda, stridea, &
               db, ldb, strideb, beta, dc, ldc, stridec, batchcount, magma_queue )
#endif
#if defined(THORNADO_OMP_OL)
      CALL stream_sync( stream )
#endif

    ELSE

#if defined(THORNADO_DEBUG_LA)
#if defined(THORNADO_GPU)
      WRITE(*,*) '[MatrixMatrixMultiplyBatched] Data not present on device'
      IF ( .not. device_is_present( ha, mydevice, sizeof_a ) ) &
        WRITE(*,*) '[MatrixMatrixMultiplyBatched]   A missing'
      IF ( .not. device_is_present( hb, mydevice, sizeof_b ) ) &
        WRITE(*,*) '[MatrixMatrixMultiplyBatched]   B missing'
      IF ( .not. device_is_present( hc, mydevice, sizeof_c ) ) &
        WRITE(*,*) '[MatrixMatrixMultiplyBatched]   C missing'
#endif
#endif

      DO i = 1, batchcount
        osa = (i-1) * stridea + 1
        osb = (i-1) * strideb + 1
        osc = (i-1) * stridec + 1
        ia = mod( (osa-1), lda ) + 1
        ib = mod( (osb-1), ldb ) + 1
        ic = mod( (osc-1), ldc ) + 1
        ja = 1
        jb = 1
        jc = 1
        IF ( stridea /= 0 ) ja = mod( (osa-1)/lda, stridea ) + 1
        IF ( strideb /= 0 ) jb = mod( (osb-1)/ldb, strideb ) + 1
        IF ( stridec /= 0 ) jc = mod( (osc-1)/ldc, stridec ) + 1
        CALL DGEMM &
               ( transa, transb, m, n, k, alpha, a(ia,ja), lda, b(ib,jb), ldb, beta, c(ic,jc), ldc )
      END DO

    END IF

  END SUBROUTINE MatrixMatrixMultiplyBatched


  SUBROUTINE LinearSolveBatched( trans, n, nrhs, a, lda, ipiv, b, ldb, info, batchcount )

    CHARACTER                          :: trans
    INTEGER                            :: n, nrhs, lda, ldb, batchcount
    REAL(DP), DIMENSION(lda,*), TARGET :: a
    REAL(DP), DIMENSION(ldb,*), TARGET :: b
    INTEGER,  DIMENSION(*),     TARGET :: ipiv
    INTEGER,  DIMENSION(*),     TARGET :: info

    INTEGER                                    :: ierr, i, stridea, strideb, strideipiv
    INTEGER(C_INT)                             :: itrans
    INTEGER(C_INT64_T)                         :: strideP_64
    INTEGER(C_SIZE_T)                          :: sizeof_a, sizeof_b, sizeof_ipiv, sizeof_info
    REAL(DP), DIMENSION(:,:), POINTER          :: pa, pb
    INTEGER,  DIMENSION(:),   POINTER          :: pipiv, pinfo
    TYPE(C_PTR)                                :: ha, hb, hipiv, hinfo
    TYPE(C_PTR), DIMENSION(batchcount), TARGET :: da, db, dipiv
    TYPE(C_PTR)                                :: da_array, db_array, dipiv_array, dinfo
    INTEGER                                    :: osa, osb
    LOGICAL                                    :: data_on_device

    data_on_device = .false.
    sizeof_a    = n * n * batchcount * c_sizeof(0.0_DP)
    sizeof_b    = n * nrhs * batchcount * c_sizeof(0.0_DP)
    sizeof_ipiv = n * batchcount * c_sizeof(0)
    sizeof_info = batchcount * c_sizeof(0)

    stridea    = n * n
    strideb    = n * nrhs
    strideipiv = n

    pa => a(:,1:n*batchcount)
    pb => b(:,1:nrhs*batchcount)
    pipiv => ipiv(1:n*batchcount)
    pinfo => info(1:batchcount)

    ha = C_LOC( pa )
    hb = C_LOC( pb )
    hipiv = C_LOC( pipiv )
    hinfo = C_LOC( pinfo )

    data_on_device = device_is_present( ha,    mydevice, sizeof_a    ) &
               .AND. device_is_present( hb,    mydevice, sizeof_b    ) &
               .AND. device_is_present( hipiv, mydevice, sizeof_ipiv ) &
               .AND. device_is_present( hinfo, mydevice, sizeof_info )

    IF ( data_on_device ) THEN

      itrans = itrans_from_char( trans )

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET ENTER DATA &
      !$OMP MAP( alloc: da, db, dipiv, da_array, db_array, dipiv_array )
#elif defined(THORNADO_OACC)
      !$ACC ENTER DATA &
      !$ACC CREATE( da, db, dipiv, da_array, db_array, dipiv_array )
#endif

      dinfo = dev_ptr( pinfo(1) )
      DO i = 1, batchcount
        osa = (i-1) * n + 1
        osb = (i-1) * nrhs + 1
        da(i) = dev_ptr( pa(1,osa) )
        db(i) = dev_ptr( pb(1,osb) )
        dipiv(i) = dev_ptr( pipiv(osa) )
      END DO
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET UPDATE TO( da, db, dipiv )
#elif defined(THORNADO_OACC)
      !$ACC UPDATE DEVICE( da, db, dipiv )
#endif

      da_array = dev_ptr( da(1) )
      db_array = dev_ptr( db(1) )
      dipiv_array = dev_ptr( dipiv(1) )
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET UPDATE TO( da_array, db_array, dipiv )
#elif defined(THORNADO_OACC)
      !$ACC UPDATE DEVICE( da_array, db_array, dipiv )
#endif

#if defined(THORNADO_LA_CUBLAS)
      ierr = cublasDgetrfBatched &
             ( cublas_handle, n, da_array, lda, dipiv(1), dinfo, batchcount )
      ierr = cublasDgetrsBatched &
             ( cublas_handle, itrans, n, nrhs, da_array, lda, dipiv(1), db_array, ldb, hinfo, batchcount )
#elif defined(THORNADO_LA_ROCM)
      !strideP_64 = n
      !CALL rocsolverCheck( rocsolver_dgetrf_batched &
      !       ( rocsolver_handle, n, n, da_array, lda, dipiv(1), strideP_64, dinfo, batchcount ) )
      !CALL rocsolverCheck( rocsolver_dgetrs_batched &
      !       ( rocsolver_handle, itrans, n, nrhs, da_array, lda, dipiv(1), strideP_64, db_array, ldb, batchcount ) )
      CALL hipblasCheck( hipblasDgetrfBatched &
             ( hipblas_handle, n, da_array, lda, dipiv(1), dinfo, batchcount ) )
      CALL hipblasCheck( hipblasDgetrsBatched &
             ( hipblas_handle, itrans, n, nrhs, da_array, lda, dipiv(1), db_array, ldb, hinfo, batchcount ) )
#elif defined(THORNADO_LA_ONEMKL)
      !$OMP TARGET VARIANT DISPATCH USE_DEVICE_PTR( a, b, ipiv )
      CALL DGETRF_BATCH_STRIDED &
             ( n, n, a, lda, stridea, ipiv, strideipiv, batchcount, info )
      CALL DGETRS_BATCH_STRIDED &
             ( trans, n, nrhs, a, lda, stridea, ipiv, strideipiv, b, ldb, strideb, info )
      !$OMP END TARGET VARIANT DISPATCH
#elif defined(THORNADO_LA_MAGMA)
      CALL magma_dgetrf_batched &
             ( n, n, da_array, lda, dipiv_array, dinfo, batchcount, magma_queue )
      CALL magma_dgetrs_batched &
             ( itrans, n, nrhs, da_array, lda, dipiv_array, db_array, ldb, batchcount, magma_queue )
#endif
#if defined(THORNADO_OMP_OL)
      CALL stream_sync( stream )
#endif

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET EXIT DATA &
      !$OMP MAP( release: da, db, dipiv, da_array, db_array, dipiv_array )
#elif defined(THORNADO_OACC)
      !$ACC EXIT DATA &
      !$ACC DELETE( da, db, dipiv, da_array, db_array, dipiv_array )
#endif

    ELSE

#if defined(THORNADO_DEBUG_LA)
#if defined(THORNADO_GPU)
      WRITE(*,*) '[LinearSolveBatched] Data not present on device'
      IF ( .not. device_is_present( ha, mydevice, sizeof_a ) ) &
        WRITE(*,*) '[LinearSolveBatched]   A missing'
      IF ( .not. device_is_present( hb, mydevice, sizeof_b ) ) &
        WRITE(*,*) '[LinearSolveBatched]   B missing'
      IF ( .not. device_is_present( hipiv, mydevice, sizeof_ipiv ) ) &
        WRITE(*,*) '[LinearSolveBatched]   ipiv missing'
      IF ( .not. device_is_present( hinfo, mydevice, sizeof_info ) ) &
        WRITE(*,*) '[LinearSolveBatched]   info missing'
#endif
#endif

      DO i = 1, batchcount
        osa = (i-1) * n + 1
        osb = (i-1) * nrhs + 1
        CALL DGETRF &
               ( n, n, a(1,osa), lda, ipiv(osa), info(i) )
        CALL DGETRS &
               ( trans, n, nrhs, a(1,osa), lda, ipiv(osa), b(1,osb), ldb, info(i) )
      END DO

    END IF

  END SUBROUTINE LinearSolveBatched


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

    pa => a(:,1:n)
    px => x(1:lenx)
    py => y(1:leny)

    ha = C_LOC( pa )
    hx = C_LOC( px )
    hy = C_LOC( py )

    data_on_device = device_is_present( ha, mydevice, sizeof_a ) &
               .AND. device_is_present( hx, mydevice, sizeof_x ) &
               .AND. device_is_present( hy, mydevice, sizeof_y )

    IF ( data_on_device ) THEN

      itrans = itrans_from_char( trans )

      da = dev_ptr( pa(1,1) )
      dx = dev_ptr( px(1) )
      dy = dev_ptr( py(1) )

#if defined(THORNADO_LA_CUBLAS)
      ierr = cublasDgemv_v2 &
             ( cublas_handle, itrans, m, n, alpha, da, lda, dx, incx, beta, dy, incy )
#elif defined(THORNADO_LA_ROCM)
      !CALL rocblasCheck( rocblas_dgemv &
      !       ( rocblas_handle, itrans, m, n, alpha, da, lda, dx, incx, beta, dy, incy ) )
      CALL hipblasCheck( hipblasDgemv &
             ( hipblas_handle, itrans, m, n, alpha, da, lda, dx, incx, beta, dy, incy ) )
#elif defined(THORNADO_LA_ONEMKL)
      !$OMP TARGET VARIANT DISPATCH USE_DEVICE_PTR( a, x, y )
      CALL DGEMV &
             ( trans, m, n, alpha, a, lda, x, incx, beta, y, incy )
      !$OMP END TARGET VARIANT DISPATCH
#elif defined(THORNADO_LA_MAGMA)
      CALL magma_dgemv &
             ( itrans, m, n, alpha, da, lda, dx, incx, beta, dy, incy, magma_queue )
#endif
#if defined(THORNADO_OMP_OL)
      CALL stream_sync( stream )
#endif

    ELSE

#if defined(THORNADO_DEBUG_LA)
#if defined(THORNADO_GPU)
      WRITE(*,*) '[MatrixVectorMultiply] Data not present on device'
      IF ( .not. device_is_present( ha, mydevice, sizeof_a ) ) &
        WRITE(*,*) '[MatrixVectorMultiply]   A missing'
      IF ( .not. device_is_present( hx, mydevice, sizeof_x ) ) &
        WRITE(*,*) '[MatrixVectorMultiply]   x missing'
      IF ( .not. device_is_present( hy, mydevice, sizeof_y ) ) &
        WRITE(*,*) '[MatrixVectorMultiply]   y missing'
#endif
#endif

      CALL DGEMV &
             ( trans, m, n, alpha, a, lda, x, incx, beta, y, incy )

    END IF

  END SUBROUTINE MatrixVectorMultiply


  SUBROUTINE MatrixDiagScale( m, n, a, lda, x, incx, c, ldc )

    CHARACTER                          :: trans
    INTEGER                            :: m, n, lda, incx, ldc
    REAL(DP), DIMENSION(lda,*), TARGET :: a
    REAL(DP), DIMENSION(ldc,*), TARGET :: c
    REAL(DP), DIMENSION(*)    , TARGET :: x

    INTEGER                            :: ierr, info, i, j, ix
    INTEGER(C_SIZE_T)                  :: sizeof_a, sizeof_c, sizeof_x
    REAL(DP), DIMENSION(:,:), POINTER  :: pa, pc
    REAL(DP), DIMENSION(:)  , POINTER  :: px
    TYPE(C_PTR)                        :: ha, hx, hc
    TYPE(C_PTR)                        :: da, dx, dc
    LOGICAL                            :: data_on_device

    data_on_device = .false.

    sizeof_a = m * n * c_sizeof(0.0_DP)
    sizeof_c = m * n * c_sizeof(0.0_DP)
    sizeof_x = m * c_sizeof(0.0_DP)

    pa => a(:,1:n)
    pc => c(:,1:n)
    px => x(1:m)

    ha = C_LOC( pa )
    hc = C_LOC( pc )
    hx = C_LOC( px )

    data_on_device = device_is_present( ha, mydevice, sizeof_a ) &
               .AND. device_is_present( hc, mydevice, sizeof_c ) &
               .AND. device_is_present( hx, mydevice, sizeof_x )

    IF ( data_on_device ) THEN

      da = dev_ptr( pa(1,1) )
      dc = dev_ptr( pc(1,1) )
      dx = dev_ptr( px(1) )

#if defined(THORNADO_LA_CUBLAS)
      ierr = cublasDdgmm &
             ( cublas_handle, CUBLAS_SIDE_LEFT, m, n, da, lda, dx, incx, dc, ldc )
#elif defined(THORNADO_LA_ROCM)
      !CALL rocblasCheck( rocblas_ddgmm &
      !       ( rocblas_handle, rocblas_side_left, m, n, da, lda, dx, incx, dc, ldc ) )
      CALL hipblasCheck( hipblasDdgmm &
             ( hipblas_handle, hipblas_SIDE_LEFT, m, n, da, lda, dx, incx, dc, ldc ) )
#elif defined(THORNADO_LA_ONEMKL)
      !!$OMP TARGET VARIANT DISPATCH USE_DEVICE_PTR( a, x, c )
      !CALL DDGMM &
      !       ( 'L', m, n, a, lda, x, incx, c, ldc )
      !!$OMP END TARGET VARIANT DISPATCH
#elif defined(THORNADO_LA_MAGMA)
      CALL magmablas_dlacpy &
             ( MagmaFull, m, n, da, lda, dc, ldc, magma_queue )
      CALL magmablas_dlascl2 &
             ( MagmaFull, m, n, dx, dc, ldc, magma_queue, info )
#endif
#if defined(THORNADO_OMP_OL)
      CALL stream_sync( stream )
#endif

    ELSE

#if defined(THORNADO_DEBUG_LA)
#if defined(THORNADO_GPU)
      WRITE(*,*) '[MatrixDiagScale] Data not present on device'
      IF ( .not. device_is_present( ha, mydevice, sizeof_a ) ) &
        WRITE(*,*) '[MatrixDiagScale]   A missing'
      IF ( .not. device_is_present( hc, mydevice, sizeof_c ) ) &
        WRITE(*,*) '[MatrixDiagScale]   C missing'
      IF ( .not. device_is_present( hx, mydevice, sizeof_x ) ) &
        WRITE(*,*) '[MatrixDiagScale]   x missing'
#endif
#endif

      IF ( incx == 1 ) THEN
        DO j = 1, n
          DO i = 1, m
            c(i,j) = a(i,j) * x(i)
          END DO
        END DO
      ELSE
        DO j = 1, n
          ix = 1
          IF ( incx < 0 ) ix = (-m+1)*incx + 1
          DO i = 1, m
            c(i,j) = a(i,j) * x(ix)
            ix = ix + incx
          END DO
        END DO
      END IF

    END IF

  END SUBROUTINE MatrixDiagScale


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

    pa => a(:,1:n)
    pb => b(:,1:nrhs)
    pwork => work(1:lwork)

    ha = C_LOC( pa )
    hb = C_LOC( pb )
    hwork = C_LOC( pwork )

    da = dev_ptr( pa(1,1) )
    db = dev_ptr( pb(1,1) )

    itrans = itrans_from_char( trans )

#if defined(THORNADO_LA_CUBLAS)
    ierr = cusolverDnDgeqrf_bufferSize &
           ( cusolver_handle, m, n, da, lda, lwork )
#elif defined(THORNADO_LA_ROCM)
#elif defined(THORNADO_LA_ONEMKL)
    !$OMP TARGET VARIANT DISPATCH USE_DEVICE_PTR( a, b, work )
    CALL DGELS &
           ( trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info )
    !$OMP END TARGET VARIANT DISPATCH
    !$OMP TARGET UPDATE FROM( work(1) )
    lwork = INT( work(1) )
#elif defined(THORNADO_LA_MAGMA)
    CALL magma_dgels_gpu &
           ( itrans, m, n, nrhs, da, lda, db, ldb, hwork, lwork, info )
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

    pa => a(:,1:n)
    pb => b(:,1:nrhs)
    ptau => tau(1:min_mn)
    pwork => work(1:lwork)
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

      da = dev_ptr( pa(1,1) )
      db = dev_ptr( pb(1,1) )
      dtau = dev_ptr( ptau(1) )
      dwork = dev_ptr( pwork(1) )
      dinfo = dev_ptr( pinfo )

#if defined(THORNADO_LA_CUBLAS)
      ierr = cusolverDnDgeqrf &
             ( cusolver_handle, m, n, da, lda, dtau, dwork, lwork, dinfo )
      ierr = cusolverDnDormqr &
             ( cusolver_handle, &
               CUBLAS_SIDE_LEFT, CUBLAS_OP_T, &
               m, nrhs, n, da, lda, dtau, db, ldb, dwork, lwork, dinfo )

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
#elif defined(THORNADO_LA_ROCM)
      CALL rocsolverCheck( rocsolver_dgeqrf &
             ( rocsolver_handle, m, n, da, lda, dtau ) )
      CALL rocsolverCheck( rocsolver_dormqr &
             ( rocsolver_handle, &
               rocblas_side_left, rocblas_operation_transpose, &
               m, nrhs, n, da, lda, dtau, db, ldb ) )

      IF ( nrhs == 1 ) THEN

        !CALL rocblasCheck( rocblas_dtrsv &
        !       ( rocblas_handle, &
        !         rocblas_fill_upper, rocblas_operation_none, rocblas_diagonal_non_unit, &
        !         n, da, lda, db, 1 ) )
        CALL hipblasCheck( hipblasDtrsv &
               ( hipblas_handle, &
                 HIPBLAS_FILL_MODE_UPPER, HIPBLAS_OP_N, HIPBLAS_DIAG_NON_UNIT, &
                 n, da, lda, db, 1 ) )

      ELSE

        !CALL rocblasCheck( rocblas_dtrsm &
        !       ( rocblas_handle, &
        !         rocblas_side_left, rocblas_fill_upper, &
        !         rocblas_operation_none, rocblas_diagonal_non_unit, &
        !         n, nrhs, One, da, lda, db, ldb ) )
        CALL hipblasCheck( hipblasDtrsm &
               ( hipblas_handle, &
                 HIPBLAS_SIDE_LEFT, HIPBLAS_FILL_MODE_UPPER, HIPBLAS_OP_N, HIPBLAS_DIAG_NON_UNIT, &
                 n, nrhs, One, da, lda, db, ldb ) )

      END IF
#elif defined(THORNADO_LA_ONEMKL)
      !$OMP TARGET VARIANT DISPATCH USE_DEVICE_PTR( a, b, work )
      CALL DGELS &
             ( trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info )
      !$OMP END TARGET VARIANT DISPATCH
#elif defined(THORNADO_LA_MAGMA)
      CALL magma_dgels_gpu &
             ( itrans, m, n, nrhs, da, lda, db, ldb, hwork, lwork, info )
#endif
#if defined(THORNADO_OMP_OL)
      CALL stream_sync( stream )
#endif

    ELSE

#if defined(THORNADO_DEBUG_LA)
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
#endif

      CALL DGELS &
             ( trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info )

    END IF

  END SUBROUTINE LinearLeastSquares


  SUBROUTINE VectorDotProductBatched( n, x, incx, stridex, y, incy, stridey, batchcount, xy )

    INTEGER                         :: n, incx, incy, stridex, stridey, batchcount
    REAL(DP), DIMENSION(*), TARGET  :: x, y, xy

    INTEGER                         :: ierr, i
    INTEGER(C_INT64_T)              :: stridex_64, stridey_64
    INTEGER(C_SIZE_T)               :: sizeof_x, sizeof_y, sizeof_xy
    REAL(DP), DIMENSION(:), POINTER :: px, py, pxy
    TYPE(C_PTR)                     :: hx, hy, hxy
    TYPE(C_PTR)                     :: dx, dy, dxy
    INTEGER                         :: osx, osy
    LOGICAL                         :: data_on_device
    REAL(DP), EXTERNAL              :: DDOT

    data_on_device = .false.
    sizeof_x  = n * batchcount * c_sizeof(0.0_DP)
    sizeof_y  = n * batchcount * c_sizeof(0.0_DP)
    sizeof_xy = batchcount * c_sizeof(0.0_DP)

    px => x(1:n*batchcount)
    py => y(1:n*batchcount)
    pxy => xy(1:batchcount)

    hx = C_LOC( px )
    hy = C_LOC( py )
    hxy = C_LOC( pxy )

    data_on_device = device_is_present( hx,  mydevice, sizeof_x  ) &
               .AND. device_is_present( hy,  mydevice, sizeof_y  ) &
               .AND. device_is_present( hxy, mydevice, sizeof_xy )

    IF ( data_on_device ) THEN

      dx = dev_ptr( px(1) )
      dy = dev_ptr( py(1) )
      dxy = dev_ptr( pxy(1) )

#if defined(THORNADO_LA_CUBLAS)
      ! Currently unavailable
      !ierr = cublasDdot_v2( cublas_handle, n, dx, incx, xnorm )
#elif defined(THORNADO_LA_ROCM)
      ! Currently unavailable
      stridex_64 = stridex
      stridey_64 = stridey
      CALL rocblasCheck( rocblas_ddot_strided_batched &
             ( rocblas_handle, n, dx, incx, stridex_64, dy, incy, stridey_64, batchcount, hxy ) )
#elif defined(THORNADO_LA_ONEMKL)
#elif defined(THORNADO_LA_MAGMA)
      ! Currently unavailable
      !xnorm = magma_ddot( n, dx, incx, magma_queue )
#endif
#if defined(THORNADO_OMP_OL)
      CALL stream_sync( stream )
#endif

    ELSE

#if defined(THORNADO_DEBUG_LA)
#if defined(THORNADO_GPU)
      WRITE(*,*) '[VectorDotProductBatched] Data not present on device'
      IF ( .not. device_is_present( hx, mydevice, sizeof_x ) ) &
        WRITE(*,*) '[VectorDotProductBatched]   x missing'
      IF ( .not. device_is_present( hy, mydevice, sizeof_y ) ) &
        WRITE(*,*) '[VectorDotProductBatched]   y missing'
      IF ( .not. device_is_present( hxy, mydevice, sizeof_xy ) ) &
        WRITE(*,*) '[VectorDotProductBatched]  xy missing'
#endif
#endif

      DO i = 1, batchcount
        osx = (i-1) * n + 1
        osy = (i-1) * n + 1
        xy(i) = DDOT( n, x(osx), incx, y(osy), incy )
      END DO

    END IF

  END SUBROUTINE VectorDotProductBatched


  SUBROUTINE VectorNorm2( n, x, incx, xnorm )

    INTEGER                         :: n, incx
    REAL(DP), DIMENSION(*), TARGET  :: x
    REAL(DP),               TARGET  :: xnorm

    INTEGER                         :: ierr
    INTEGER(C_SIZE_T)               :: sizeof_x
    REAL(DP), DIMENSION(:), POINTER :: px
    REAL(DP),               POINTER :: pxnorm
    TYPE(C_PTR)                     :: hx, hxnorm
    TYPE(C_PTR)                     :: dx, dxnorm
    LOGICAL                         :: data_on_device
    REAL(DP), EXTERNAL              :: DNRM2

    data_on_device = .false.
    sizeof_x = n * c_sizeof(0.0_DP)

    px => x(1:n)
    pxnorm => xnorm

    hx = C_LOC( px )
    hxnorm = C_LOC( pxnorm )

    data_on_device = device_is_present( hx, mydevice, sizeof_x )

    IF ( data_on_device ) THEN

      dx = dev_ptr( px(1) )
      dxnorm = dev_ptr( pxnorm )

#if defined(THORNADO_LA_CUBLAS)
      ierr = cublasDnrm2_v2( cublas_handle, n, dx, incx, xnorm )
#elif defined(THORNADO_LA_ROCM)
      ! Currently unavailable
      !CALL rocblasCheck( rocblas_dnrm2( rocblas_handle, n, dx, incx, xnorm ) )
      CALL hipblasCheck( hipblasDnrm2( hipblas_handle, n, dx, incx, hxnorm ) )
#elif defined(THORNADO_LA_ONEMKL)
      !$OMP TARGET VARIANT DISPATCH USE_DEVICE_PTR( x )
      xnorm = DNRM2( n, x, incx )
      !$OMP END TARGET VARIANT DISPATCH
#elif defined(THORNADO_LA_MAGMA)
      xnorm = magma_dnrm2( n, dx, incx, magma_queue )
#endif
#if defined(THORNADO_OMP_OL)
      CALL stream_sync( stream )
#endif

    ELSE

#if defined(THORNADO_DEBUG_LA)
#if defined(THORNADO_GPU)
      WRITE(*,*) '[VectorNorm2] Data not present on device'
      IF ( .not. device_is_present( hx, mydevice, sizeof_x ) ) &
        WRITE(*,*) '[VectorNorm2]   x missing'
#endif
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


  SUBROUTINE VectorVectorAdd( n, alpha, x, incx, y, incy )

    INTEGER                         :: n, incx, incy
    REAL(DP)                        :: alpha
    REAL(DP), DIMENSION(*), TARGET  :: x, y

    INTEGER                         :: ierr
    INTEGER(C_SIZE_T)               :: sizeof_x, sizeof_y
    REAL(DP), DIMENSION(:), POINTER :: px, py
    TYPE(C_PTR)                     :: hx, hy
    TYPE(C_PTR)                     :: dx, dy
    LOGICAL                         :: data_on_device

    data_on_device = .false.

    sizeof_x = n * c_sizeof(0.0_DP)
    sizeof_y = n * c_sizeof(0.0_DP)

    px => x(1:n)
    py => y(1:n)

    hx = C_LOC( px )
    hy = C_LOC( py )

    data_on_device = device_is_present( hx, mydevice, sizeof_x ) &
               .AND. device_is_present( hy, mydevice, sizeof_y )

    IF ( data_on_device ) THEN

      dx = dev_ptr( px(1) )
      dy = dev_ptr( py(1) )

#if defined(THORNADO_LA_CUBLAS)
      ierr = cublasDaxpy_v2( cublas_handle, n, alpha, dx, incx, dy, incy )
#elif defined(THORNADO_LA_ROCM)
      !CALL rocblasCheck( rocblas_daxpy( rocblas_handle, n, alpha, dx, incx, dy, incy ) )
      CALL hipblasCheck( hipblasDaxpy( hipblas_handle, n, alpha, dx, incx, dy, incy ) )
#elif defined(THORNADO_LA_ONEMKL)
      !$OMP TARGET VARIANT DISPATCH USE_DEVICE_PTR( x, y )
      CALL DAXPY( n, alpha, x, incx, y, incy )
      !$OMP END TARGET VARIANT DISPATCH
#elif defined(THORNADO_LA_MAGMA)
      CALL magma_daxpy( n, alpha, dx, incx, dy, incy, magma_queue )
#endif
#if defined(THORNADO_OMP_OL)
      CALL stream_sync( stream )
#endif

    ELSE

#if defined(THORNADO_DEBUG_LA)
#if defined(THORNADO_GPU)
      WRITE(*,*) '[VectorVectorAdd] Data not present on device'
      IF ( .not. device_is_present( hx, mydevice, sizeof_x ) ) &
        WRITE(*,*) '[VectorVectorAdd]   x missing'
#endif
#endif

      CALL DAXPY( n, alpha, x, incx, y, incy )

    END IF

  END SUBROUTINE VectorVectorAdd


  SUBROUTINE EigenvaluesSymmetric3( A, Lambda )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: A(3,3)
    REAL(DP), INTENT(out) :: Lambda(3)

    REAL(DP) :: B11, B22, B33
    REAL(DP) :: B12, B13, B21, B23, B31, B32
    REAL(DP) :: P1, P2, P, Q, R, PHI, DETB

    P1 = A(1,2)**2 + A(1,3)**2 + A(2,3)**2

    IF ( P1 == Zero ) THEN

      Lambda(1) = A(1,1)
      Lambda(2) = A(2,2)
      Lambda(3) = A(3,3)

    ELSE

      Q = ( A(1,1) + A(2,2) + A(3,3) ) / 3.0_DP
      P2 = 2.0_DP * P1 &
           + ( A(1,1) - Q )**2 &
           + ( A(2,2) - Q )**2 &
           + ( A(3,3) - Q )**2
      P = SQRT( P2 / 6.0_DP )

      B11 = ( A(1,1) - Q ) / P
      B22 = ( A(2,2) - Q ) / P
      B33 = ( A(3,3) - Q ) / P
      B12 = A(1,2) / P ; B21 = B12
      B13 = A(1,3) / P ; B31 = B13
      B23 = A(2,3) / P ; B32 = B23
      DETB =   B11 * B22 * B33  &
             - B11 * B23 * B32  &
             - B12 * B21 * B33  &
             + B12 * B23 * B31  &
             + B13 * B21 * B32  &
             - B13 * B22 * B31
      R = DETB * 0.5_DP
      IF ( R <= - One ) THEN
        PHI = Pi
      ELSE IF ( R >= One ) THEN
        PHI = Zero
      ELSE
        PHI = ACOS( R ) / 3.0_DP
      END IF

      Lambda(1) = Q + 2.0_DP * P * COS( PHI )
      Lambda(3) = Q + 2.0_DP * P * COS( PHI + ( 2.0_DP * Pi / 3.0_DP ) )
      Lambda(2) = 3.0_DP * Q - Lambda(1) - Lambda(3)

    END IF

  END SUBROUTINE EigenvaluesSymmetric3




END MODULE LinearAlgebraModule
