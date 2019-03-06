MODULE LinearAlgebraModule

  USE, INTRINSIC :: ISO_C_BINDING
  USE KindModule, ONLY: &
    DP
#ifdef USE_GPU
  USE cublasf, ONLY: &
    cublas_handle, cublasDgemm_v2, CUBLAS_OP_N, CUBLAS_OP_T
#endif

  !$ USE OMP_LIB

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MatrixMatrixMultiply

  !INTERFACE MatrixMatrixMultiply
  !  MODULE PROCEDURE MatrixMatrixMultiply_DP
  !END INTERFACE MatrixMatrixMultiply

  INTERFACE
    INTEGER(C_INT) FUNCTION omp_target_is_present( hostptr, device ) &
        BIND(C, NAME = "omp_target_is_present" )
      USE, INTRINSIC :: ISO_C_BINDING
      TYPE(C_PTR),    VALUE :: hostptr
      INTEGER(C_INT), VALUE :: device
    END FUNCTION omp_target_is_present
  END INTERFACE

CONTAINS

  SUBROUTINE MatrixMatrixMultiply( transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc )

    CHARACTER                          :: transa, transb
    INTEGER                            :: m, n, k, lda, ldb, ldc
    REAL(DP)                           :: alpha, beta
    REAL(DP), DIMENSION(lda,*), TARGET :: a
    REAL(DP), DIMENSION(ldb,*), TARGET :: b
    REAL(DP), DIMENSION(ldc,*), TARGET :: c

    INTEGER                            :: ierr, mydevice
    INTEGER(C_INT)                     :: itransa, itransb
    TYPE(C_PTR)                        :: ha, hb, hc
    TYPE(C_PTR)                        :: da, db, dc

    LOGICAL                            :: use_cublas

    use_cublas = .false.
#ifdef USE_GPU
    ha = C_LOC( a )
    hb = C_LOC( b )
    hc = C_LOC( c )
    mydevice = omp_get_default_device()
    use_cublas = ( omp_target_is_present( ha, mydevice ) > 0 ) &
           .AND. ( omp_target_is_present( hb, mydevice ) > 0 ) &
           .AND. ( omp_target_is_present( hc, mydevice ) > 0 )
    IF ( .NOT. use_cublas ) WRITE(*,*) '[MatrixMatrixMultiply] Data not present on device'

    IF ( use_cublas ) THEN

      IF ( transa == 'T' ) THEN
        itransa = CUBLAS_OP_T
      ELSE
        itransa = CUBLAS_OP_N
      END IF

      IF ( transb == 'T' ) THEN
        itransb = CUBLAS_OP_T
      ELSE
        itransb = CUBLAS_OP_N
      END IF

      !$OMP TARGET DATA USE_DEVICE_PTR( a, b, c )
      da = C_LOC( a )
      db = C_LOC( b )
      dc = C_LOC( c )
      !$OMP END TARGET DATA
    
      ierr = cublasDgemm_v2 &
             ( cublas_handle, itransa, itransb, m, n, k, alpha, da, lda, db, ldb, beta, dc, ldc )

    ELSE

      CALL DGEMM &
             ( transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc )

    END IF

#else

    CALL DGEMM &
           ( transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc )

#endif

  END SUBROUTINE MatrixMatrixMultiply

END MODULE LinearAlgebraModule
