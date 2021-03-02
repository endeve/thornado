MODULE TwoMoment_UtilitiesModule_OrderV

  USE KindModule, ONLY: &
    DP, Zero, Half, One, Two, Three, Five
  USE ProgramHeaderModule, ONLY: &
    nDOFZ, nDOFX, nDOFE
  USE GeometryFieldsModule, ONLY: &
    nGF, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33
  USE FluidFieldsModule, ONLY: &
    nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne
  USE Euler_UtilitiesModule_NonRelativistic, ONLY: &
    ComputePrimitive_Euler_NonRelativistic
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3, &
    nPR, iPR_D, iPR_I1, iPR_I2, iPR_I3
  USE TwoMoment_ClosureModule, ONLY: &
    FluxFactor, &
    EddingtonFactor, &
    HeatFluxFactor
  USE TwoMoment_TimersModule_OrderV, ONLY: &
    TimersStart, &
    TimersStop, &
    Timer_Streaming_NumericalFlux_InOut, &
    Timer_Streaming_NumericalFlux_RHS, &
    Timer_Streaming_NumericalFlux_LS, &
    Timer_Streaming_NumericalFlux_Update
  USE LinearAlgebraModule, ONLY: &
    MatrixMatrixMultiplyBatched

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeComputePrimitive_TwoMoment
  PUBLIC :: FinalizeComputePrimitive_TwoMoment
  PUBLIC :: ComputePrimitive_TwoMoment
  PUBLIC :: ComputePrimitive_TwoMoment_Scalar
  PUBLIC :: ComputePrimitive_TwoMoment_Vector
  PUBLIC :: ComputeConserved_TwoMoment
  PUBLIC :: ComputeFromConserved_TwoMoment
  PUBLIC :: Flux_E
  PUBLIC :: Flux_X1
  PUBLIC :: Flux_X2
  PUBLIC :: Flux_X3
  PUBLIC :: ComputeEddingtonTensorComponents_uu
  PUBLIC :: ComputeEddingtonTensorComponents_dd
  PUBLIC :: ComputeEddingtonTensorComponents_ud
  PUBLIC :: ComputeHeatFluxTensorComponents_udd
  PUBLIC :: NumericalFlux_LLF

  INTERFACE ComputePrimitive_TwoMoment
    MODULE PROCEDURE ComputePrimitive_TwoMoment_Scalar
    MODULE PROCEDURE ComputePrimitive_TwoMoment_Vector
  END INTERFACE ComputePrimitive_TwoMoment

CONTAINS

  SUBROUTINE InitializeComputePrimitive_TwoMoment
  END SUBROUTINE InitializeComputePrimitive_TwoMoment

  SUBROUTINE FinalizeComputePrimitive_TwoMoment
  END SUBROUTINE FinalizeComputePrimitive_TwoMoment

  SUBROUTINE ComputePrimitive_TwoMoment_Vector &
    ( N, G_d_1, G_d_2, G_d_3, D, I_u_1, I_u_2, I_u_3, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, PositionIndexZ, nIterations_Option )

    REAL(DP), DIMENSION(:), INTENT(in)  :: N, G_d_1, G_d_2, G_d_3
    REAL(DP), DIMENSION(:), INTENT(out) :: D, I_u_1, I_u_2, I_u_3
    REAL(DP), DIMENSION(:), INTENT(in)  :: V_u_1, V_u_2, V_u_3
    REAL(DP), DIMENSION(:), INTENT(in)  :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    INTEGER,  DIMENSION(:), INTENT(in)  :: PositionIndexZ
    INTEGER,  DIMENSION(:), INTENT(out), OPTIONAL :: nIterations_Option

    ! --- Parameters ---

    INTEGER,  PARAMETER :: M = 2
    INTEGER,  PARAMETER :: MaxIterations = 100
    REAL(DP), PARAMETER :: Rtol = 1.0d-08

    ! --- Local Variables ---

    INTEGER  :: nZ
    INTEGER  :: iX, iZ
    INTEGER  :: k, Mk, iM, i, j

    REAL(DP) :: FTMP(4,M), GTMP(4,M)
    REAL(DP) :: k_dd(3,3), SUM1, DET, LMAT(4,4)
    REAL(DP) :: A_d_1, A_d_2, A_d_3
    LOGICAL  :: CONVERGED

    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: FVEC, GVEC
    REAL(DP), DIMENSION(:,:),   ALLOCATABLE :: CVEC, UVEC, FVECm, GVECm, Alpha
    LOGICAL,  DIMENSION(:),     ALLOCATABLE :: ITERATE
    INTEGER,  DIMENSION(:),     ALLOCATABLE :: nIterations

    nZ = SIZE( N, 1 )

    CALL TimersStart( Timer_Streaming_NumericalFlux_InOut )

    ALLOCATE( FVEC(4,M,nZ) )
    ALLOCATE( GVEC(4,M,nZ) )

    ALLOCATE( CVEC (4,nZ) )
    ALLOCATE( UVEC (4,nZ) )
    ALLOCATE( FVECm(4,nZ) )
    ALLOCATE( GVECm(4,nZ) )
    ALLOCATE( Alpha(M,nZ) )

    ALLOCATE( ITERATE(nZ) )
    ALLOCATE( nIterations(nZ) )

    ITERATE = .TRUE.

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: ITERATE ) &
    !$OMP MAP( alloc: FVEC, GVEC, CVEC, UVEC, &
    !$OMP             FVECm, GVECm, Alpha, nIterations )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA ASYNC &
    !$ACC COPYIN( ITERATE ) &
    !$ACC CREATE( FVEC, GVEC, CVEC, UVEC, &
    !$ACC         FVECm, GVECm, Alpha, nIterations )
#endif

    ! --- Initial Guess ---

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR ASYNC &
    !$ACC PRESENT( CVEC, UVEC, N, G_d_1, G_d_2, G_d_3, &
    !$ACC          D, I_u_1, I_u_2, I_u_3, &
    !$ACC          Gm_dd_11, Gm_dd_22, Gm_dd_33 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO
#endif
    DO iZ = 1, nZ
      CVEC(iCR_N ,iZ) = N    (iZ)
      CVEC(iCR_G1,iZ) = G_d_1(iZ)
      CVEC(iCR_G2,iZ) = G_d_2(iZ)
      CVEC(iCR_G3,iZ) = G_d_3(iZ)

      D    (iZ) = N(iZ)
      I_u_1(iZ) = Zero
      I_u_2(iZ) = Zero
      I_u_3(iZ) = Zero
    END DO

    CALL TimersStop( Timer_Streaming_NumericalFlux_InOut )

    k = 0
    DO WHILE( ANY( ITERATE ) .AND. k < MaxIterations )

      k = k + 1
      Mk = MIN( M, k )

      CALL TimersStart( Timer_Streaming_NumericalFlux_RHS )

#if   defined( THORNADO_OMP_OL )
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
      !$OMP PRIVATE( iX, A_d_1, A_d_2, A_d_3, k_dd, DET, LMAT )
#elif defined( THORNADO_OACC   )
      !$ACC PARALLEL LOOP GANG VECTOR ASYNC &
      !$ACC PRIVATE( iX, A_d_1, A_d_2, A_d_3, k_dd, DET, LMAT ) &
      !$ACC PRESENT( ITERATE, UVEC, CVEC, GVEC, FVEC, GVECm, FVECm, &
      !$ACC          PositionIndexZ, D, I_u_1, I_u_2, I_u_3, &
      !$ACC          Gm_dd_11, Gm_dd_22, Gm_dd_33, V_u_1, V_u_2, V_u_3 )
#elif defined( THORNADO_OMP    )
      !$OMP PARALLEL DO &
      !$OMP PRIVATE( iX, A_d_1, A_d_2, A_d_3, k_dd, DET, LMAT )
#endif
      DO iZ = 1, nZ
        IF ( ITERATE(iZ) ) THEN

          iX = PositionIndexZ(iZ)

          UVEC(iPR_D ,iZ) = D    (iZ)
          UVEC(iPR_I1,iZ) = I_u_1(iZ) * Gm_dd_11(iX)
          UVEC(iPR_I2,iZ) = I_u_2(iZ) * Gm_dd_22(iX)
          UVEC(iPR_I3,iZ) = I_u_3(iZ) * Gm_dd_33(iX)

          k_dd = EddingtonTensorComponents_dd &
                   ( D(iZ), I_u_1(iZ), I_u_2(iZ), I_u_3(iZ), &
                     Gm_dd_11(iX), Gm_dd_22(iX), Gm_dd_33(iX) )

          A_d_1 &
            =   V_u_1(iX) * k_dd(1,1) &
              + V_u_2(iX) * k_dd(2,1) &
              + V_u_3(iX) * k_dd(3,1)

          A_d_2 &
            =   V_u_1(iX) * k_dd(1,2) &
              + V_u_2(iX) * k_dd(2,2) &
              + V_u_3(iX) * k_dd(3,2)

          A_d_3 &
            =   V_u_1(iX) * k_dd(1,3) &
              + V_u_2(iX) * k_dd(2,3) &
              + V_u_3(iX) * k_dd(3,3)

          DET &
            = One &
              - ( V_u_1(iX) * A_d_1 &
                + V_u_2(iX) * A_d_2 &
                + V_u_3(iX) * A_d_3 )

          LMAT(1,1) = One
          LMAT(2,1) = - V_u_1(iX)
          LMAT(3,1) = - V_u_2(iX)
          LMAT(4,1) = - V_u_3(iX)

          LMAT(1,2) = - A_d_1
          LMAT(2,2) = One - ( V_u_2(iX) * A_d_2 + V_u_3(iX) * A_d_3 )
          LMAT(3,2) = V_u_2(iX) * A_d_1
          LMAT(4,2) = V_u_3(iX) * A_d_1

          LMAT(1,3) = - A_d_2
          LMAT(2,3) = V_u_1(iX) * A_d_2
          LMAT(3,3) = One - ( V_u_1(iX) * A_d_1 + V_u_3(iX) * A_d_3 )
          LMAT(4,3) = V_u_3(iX) * A_d_2

          LMAT(1,4) = - A_d_3
          LMAT(2,4) = V_u_1(iX) * A_d_3
          LMAT(3,4) = V_u_2(iX) * A_d_3
          LMAT(4,4) = One - ( V_u_1(iX) * A_d_1 + V_u_2(iX) * A_d_2 )

          DO i = 1, 4
            SUM1 = Zero
            DO j = 1, 4
              SUM1 = SUM1 + LMAT(j,i) * CVEC(j,iZ) / DET
            END DO
            GVECm(i,iZ) = SUM1
            FVECm(i,iZ) = GVECm(i,iZ) - UVEC(i,iZ)

            GVEC(i,Mk,iZ) = GVECm(i,iZ)
            FVEC(i,Mk,iZ) = FVECm(i,iZ)
          END DO

        END IF
      END DO

      CALL TimersStop( Timer_Streaming_NumericalFlux_RHS )

      CALL TimersStart( Timer_Streaming_NumericalFlux_LS )

      IF ( Mk > 1 ) THEN
        CALL Alpha_LS_Vector &
               ( ITERATE, nZ, M, Mk, FVECm, FVEC, Alpha )

#if   defined( THORNADO_OMP_OL )
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
        !$OMP PRIVATE( SUM1 )
#elif defined( THORNADO_OACC   )
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) ASYNC &
        !$ACC PRIVATE( SUM1 ) &
        !$ACC PRESENT( ITERATE, GVECm, FVECm, GVEC, UVEC, Alpha )
#elif defined( THORNADO_OMP    )
        !$OMP PARALLEL DO COLLAPSE(2) &
        !$OMP PRIVATE( SUM1 )
#endif
        DO iZ = 1, nZ
          DO i = 1, 4
            IF ( ITERATE(iZ) ) THEN
              SUM1 = Zero
              DO iM = 1, Mk
                SUM1 = SUM1 + GVEC(i,iM,iZ) * Alpha(iM,iZ)
              END DO
              GVECm(i,iZ) = SUM1
              FVECm(i,iZ) = GVECm(i,iZ) - UVEC(i,iZ)
            END IF
          END DO
        END DO
      END IF

      CALL TimersStop( Timer_Streaming_NumericalFlux_LS )

      CALL TimersStart( Timer_Streaming_NumericalFlux_Update )

#if   defined( THORNADO_OMP_OL )
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
      !$OMP PRIVATE( iX, CONVERGED, FTMP, GTMP )
#elif defined( THORNADO_OACC   )
      !$ACC PARALLEL LOOP GANG VECTOR ASYNC &
      !$ACC PRIVATE( iX, CONVERGED, FTMP, GTMP ) &
      !$ACC PRESENT( ITERATE, UVEC, CVEC, GVECm, FVECm, GVEC, FVEC, &
      !$ACC          PositionIndexZ, D, I_u_1, I_u_2, I_u_3, &
      !$ACC          Gm_dd_11, Gm_dd_22, Gm_dd_33, nIterations )
#elif defined( THORNADO_OMP    )
      !$OMP PARALLEL DO &
      !$OMP PRIVATE( iX, CONVERGED, FTMP, GTMP )
#endif
      DO iZ = 1, nZ
        IF ( ITERATE(iZ) ) THEN

          iX = PositionIndexZ(iZ)

          D    (iZ) = GVECm(iPR_D ,iZ)
          I_u_1(iZ) = GVECm(iPR_I1,iZ) / Gm_dd_11(iX)
          I_u_2(iZ) = GVECm(iPR_I2,iZ) / Gm_dd_22(iX)
          I_u_3(iZ) = GVECm(iPR_I3,iZ) / Gm_dd_33(iX)

          CONVERGED = ALL( ABS( FVECm(:,iZ) ) <= Rtol * ABS( CVEC(:,iZ) ) )

          IF ( CONVERGED ) THEN
            ITERATE(iZ) = .FALSE.
            nIterations(iZ) = k
          ELSE IF ( Mk == M ) THEN
            DO j = 1, Mk - 1
              DO i = 1, 4
                FTMP(i,j) = FVEC(i,j+1,iZ)
                GTMP(i,j) = GVEC(i,j+1,iZ)
              END DO
            END DO
            DO j = 1, Mk - 1
              DO i = 1, 4
                FVEC(i,j,iZ) = FTMP(i,j)
                GVEC(i,j,iZ) = GTMP(i,j)
              END DO
            END DO
          END IF
        END IF
      END DO
#if   defined( THORNADO_OMP_OL )
      !$OMP TARGET UPDATE FROM( ITERATE )
#elif defined( THORNADO_OACC   )
      !$ACC UPDATE HOST( ITERATE ) ASYNC
      !$ACC WAIT
#endif

      CALL TimersStop( Timer_Streaming_NumericalFlux_Update )

    END DO

    CALL TimersStart( Timer_Streaming_NumericalFlux_InOut )

    IF( PRESENT( nIterations_Option ) ) THEN
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR ASYNC &
      !$ACC PRESENT( nIterations, nIterations_Option )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO
#endif
      DO iZ = 1, nZ
        nIterations_Option(iZ) = nIterations(iZ)
      END DO
    END IF


#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: FVEC, GVEC, CVEC, UVEC, &
    !$OMP               FVECm, GVECm, Alpha, ITERATE, nIterations )
#elif defined( THORNADO_OACC   )
    !$ACC EXIT DATA WAIT &
    !$ACC DELETE( FVEC, GVEC, CVEC, UVEC, &
    !$ACC         FVECm, GVECm, Alpha, ITERATE, nIterations )
#endif

    DEALLOCATE( FVEC, GVEC )
    DEALLOCATE( CVEC, UVEC, FVECm, GVECm, Alpha )
    DEALLOCATE( ITERATE, nIterations )

    CALL TimersStop( Timer_Streaming_NumericalFlux_InOut )

  END SUBROUTINE ComputePrimitive_TwoMoment_Vector


  SUBROUTINE Alpha_LS_Vector &
    ( MASK, nZ, M, Mk, Fm, F, Alpha )

    LOGICAL,  DIMENSION(:),     INTENT(in)    :: MASK
    INTEGER,                    INTENT(in)    :: nZ, M, Mk
    REAL(DP), DIMENSION(:,:),   INTENT(inout) :: Fm
    REAL(DP), DIMENSION(:,:,:), INTENT(inout) :: F
    REAL(DP), DIMENSION(:,:),   INTENT(inout) :: Alpha

    REAL(DP) :: AA11, AA12, AA21, AA22, AB1, AB2, DET_AA, SUM1
    REAL(DP) :: A1, A2, B
    INTEGER  :: iZ, iP, iM

    IF ( Mk > 1 ) THEN

      IF ( Mk == 2 ) THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
        !$OMP PRIVATE( AA11, AB1, B )
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR ASYNC &
        !$ACC PRIVATE( AA11, AB1, B ) &
        !$ACC PRESENT( MASK, Alpha, F, Fm )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO &
        !$OMP PRIVATE( AA11, AB1, B )
#endif
        DO iZ = 1, nZ
          IF ( MASK(iZ) ) THEN

            AA11 = Zero
            AB1 = Zero

            DO iP = 1, 4

              A1 = F(iP,1,iZ) - Fm(iP,iZ)
              B  = - Fm(iP,iZ)

              AA11 = AA11 + A1 * A1
              AB1  = AB1  + A1 * B

            END DO

            Alpha(1,iZ) = AB1 / AA11
            Alpha(2,iZ) = One - Alpha(1,iZ)

          END IF
        END DO

      ELSE IF ( Mk == 3 ) THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
        !$OMP PRIVATE( AA11, AA12, AA22, AB1, AB2, DET_AA, A1, A2, B )
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR ASYNC &
        !$ACC PRIVATE( AA11, AA12, AA22, AB1, AB2, DET_AA, A1, A2, B ) &
        !$ACC PRESENT( MASK, Alpha, F, Fm )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO &
        !$OMP PRIVATE( AA11, AA12, AA22, AB1, AB2, DET_AA, A1, A2, B )
#endif
        DO iZ = 1, nZ
          IF ( MASK(iZ) ) THEN

            AA11 = Zero
            AA12 = Zero
            AA22 = Zero
            AB1  = Zero
            AB2  = Zero

            DO iP = 1, 4

              A1 = F(iP,1,iZ) - Fm(iP,iZ)
              A2 = F(iP,2,iZ) - Fm(iP,iZ)
              B  = - Fm(iP,iZ)

              AA11 = AA11 + A1 * A1
              AA12 = AA12 + A1 * A2
              AA22 = AA22 + A2 * A2

              AB1  = AB1  + A1 * B
              AB2  = AB2  + A2 * B

            END DO

            DET_AA = AA11*AA22 - AA12*AA12

            Alpha(1,iZ) = ( + AA22 * AB1 - AA12 * AB2 ) / DET_AA
            Alpha(2,iZ) = ( - AA12 * AB1 + AA11 * AB2 ) / DET_AA
            Alpha(3,iZ) = One - Alpha(1,iZ) - Alpha(2,iZ)

          END IF
        END DO

      ELSE IF ( Mk > 3 ) THEN

        ! --- Not Implemented ---

      END IF

    END IF

  END SUBROUTINE Alpha_LS_Vector


  SUBROUTINE ComputePrimitive_TwoMoment_Scalar &
    ( N, G_d_1, G_d_2, G_d_3, D, I_u_1, I_u_2, I_u_3, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, nIterations_Option )

#if   defined( THORNADO_OMP_OL )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC   )
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: N, G_d_1, G_d_2, G_d_3 ! --- Index Down
    REAL(DP), INTENT(out) :: D, I_u_1, I_u_2, I_u_3 ! --- Index Up
    REAL(DP), INTENT(in)  ::    V_u_1, V_u_2, V_u_3 ! --- Index Up
    REAL(DP), INTENT(in)  :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    INTEGER, INTENT(out), OPTIONAL :: nIterations_Option

    ! --- Parameters ---

    INTEGER,  PARAMETER :: M = 2
    INTEGER,  PARAMETER :: MaxIterations = 100
    REAL(DP), PARAMETER :: Rtol = 1.0d-08

    ! --- Local Variables ---

    LOGICAL  :: CONVERGED
    INTEGER  :: i, j, k, mk
    REAL(DP) :: I_d_1, I_d_2, I_d_3, A_d_1, A_d_2, A_d_3, k_dd(3,3)
    REAL(DP) :: UVEC(4), CVEC(4)
    REAL(DP) :: GVEC(4,M), GVECm(4)
    REAL(DP) :: FVEC(4,M), FVECm(4)
    REAL(DP) :: LMAT(4,4), DET, Alpha(M)

    CVEC = [ N, G_d_1, G_d_2, G_d_3 ]

    ! --- Initial Guess ---

    D     = N
    I_u_1 = Zero
    I_u_2 = Zero
    I_u_3 = Zero

    I_d_1 = Gm_dd_11 * I_u_1
    I_d_2 = Gm_dd_22 * I_u_2
    I_d_3 = Gm_dd_33 * I_u_3

    k = 0
    CONVERGED = .FALSE.
    DO WHILE( .NOT. CONVERGED .AND. k < MaxIterations )

      k  = k + 1
      mk = MIN( M, k )

      UVEC = [ D, I_d_1, I_d_2, I_d_3 ]

      k_dd = EddingtonTensorComponents_dd &
               ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

      A_d_1 = V_u_1 * k_dd(1,1) + V_u_2 * k_dd(1,2) + V_u_3 * k_dd(1,3)
      A_d_2 = V_u_1 * k_dd(1,2) + V_u_2 * k_dd(2,2) + V_u_3 * k_dd(2,3)
      A_d_3 = V_u_1 * k_dd(1,3) + V_u_2 * k_dd(2,3) + V_u_3 * k_dd(3,3)

      DET = One - ( V_u_1 * A_d_1 + V_u_2 * A_d_2 + V_u_3 * A_d_3 )

      LMAT(1,1) = One
      LMAT(2,1) = - A_d_1
      LMAT(3,1) = - A_d_2
      LMAT(4,1) = - A_d_3

      LMAT(1,2) = - V_u_1
      LMAT(2,2) = One - ( V_u_2 * A_d_2 + V_u_3 * A_d_3 )
      LMAT(3,2) = V_u_1 * A_d_2
      LMAT(4,2) = V_u_1 * A_d_3

      LMAT(1,3) = - V_u_2
      LMAT(2,3) = V_u_2 * A_d_1
      LMAT(3,3) = One - ( V_u_1 * A_d_1 + V_u_3 * A_d_3 )
      LMAT(4,3) = V_u_2 * A_d_3

      LMAT(1,4) = - V_u_3
      LMAT(2,4) = V_u_3 * A_d_1
      LMAT(3,4) = V_u_3 * A_d_2
      LMAT(4,4) = One - ( V_u_1 * A_d_1 + V_u_2 * A_d_2 )

      LMAT = LMAT / DET

      ! --- Multiply LMAT and CVEC to form GVEC ---

      GVEC(:,mk) = Zero

      DO j = 1, 4
      DO i = 1, 4

        GVEC(i,mk) = GVEC(i,mk) + LMAT(i,j) * CVEC(j)

      END DO
      END DO

      FVEC(:,mk) = GVEC(:,mk) - UVEC

      IF( mk == 1 )THEN

        ! --- Picard Iteration ---

        GVECm = GVEC(:,mk)

      ELSE

        ! --- Anderson Accelerated Fixed-Point ---

        Alpha = Alpha_LS( M, mk, FVEC )

        GVECm = Zero
        DO i = 1, mk

          GVECm = GVECm + Alpha(i) * GVEC(:,i)

        END DO

      END IF

      FVECm = GVECm - UVEC

      IF( ALL( ABS( FVECm ) <= Rtol * ABS( CVEC ) ) )THEN

        CONVERGED = .TRUE.

      END IF

      UVEC = GVECm

      IF( mk == M .AND. .NOT. CONVERGED )THEN

        FVEC = ShiftVec( M, mk, FVEC )
        GVEC = ShiftVec( M, mk, GVEC )

      END IF

      D     = UVEC(1)
      I_d_1 = UVEC(2); I_u_1 = I_d_1 / Gm_dd_11
      I_d_2 = UVEC(3); I_u_2 = I_d_2 / Gm_dd_22
      I_d_3 = UVEC(4); I_u_3 = I_d_3 / Gm_dd_33

    END DO

    IF( PRESENT( nIterations_Option ) )THEN

      nIterations_Option = k

    END IF

!    IF( k == MaxIterations )THEN
!
!      PRINT*
!      PRINT*, "ComputePrimitive_TwoMoment"
!      PRINT*
!      PRINT*, "  N     = ", N
!      PRINT*, "  G_d_1 = ", G_d_1
!      PRINT*, "  G_d_2 = ", G_d_2
!      PRINT*, "  G_d_3 = ", G_d_3
!      PRINT*
!      PRINT*, "  V_u_1 = ", V_u_1
!      PRINT*, "  V_u_2 = ", V_u_2
!      PRINT*, "  V_u_3 = ", V_u_3
!      PRINT*
!
!      PRINT*, "  Converged with k = ", k
!
!      PRINT*
!      PRINT*, "  FVECm = ", FVECm
!      PRINT*
!
!      PRINT*
!      PRINT*, "  D     = ", D
!      PRINT*, "  I_u_1 = ", I_u_1
!      PRINT*, "  I_u_2 = ", I_u_2
!      PRINT*, "  I_u_3 = ", I_u_3
!      PRINT*
!
!    END IF

  END SUBROUTINE ComputePrimitive_TwoMoment_Scalar


  FUNCTION Alpha_LS( M, mk, FVEC )

#if   defined( THORNADO_OMP_OL )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC   )
    !$ACC ROUTINE SEQ
#endif

    INTEGER,  INTENT(in) :: M, mk
    REAL(DP), INTENT(in) :: FVEC(4,M)
    REAL(DP)             :: Alpha_LS(M)

    INTEGER  :: i
    REAL(DP) :: BVEC(4), AMAT(4,M)
    REAL(DP) :: AA11, AA12, AA22, AB1, AB2, DET_AA, SUM1

    BVEC = - FVEC(:,mk)

    DO i = 1, mk - 1

      AMAT(:,i) = FVEC(:,i) - FVEC(:,mk)

    END DO

    IF( mk == 2 )THEN

      AA11 = Zero
      AB1  = Zero

      DO i = 1, 4

        AA11 = AA11 + AMAT(i,1) * AMAT(i,1)
        AB1  = AB1  + AMAT(i,1) * BVEC(i)

      END DO

      BVEC(1) = AB1 / AA11

    ELSEIF( mk == 3 )THEN

      AA11 = Zero
      AA12 = Zero
      AA22 = Zero
      AB1  = Zero
      AB2  = Zero

      DO i = 1, 4

        AA11 = AA11 + AMAT(i,1) * AMAT(i,1)
        AA12 = AA12 + AMAT(i,1) * AMAT(i,2)
        AA22 = AA22 + AMAT(i,2) * AMAT(i,2)
        AB1  = AB1  + AMAT(i,1) * BVEC(i)
        AB2  = AB2  + AMAT(i,2) * BVEC(i)

      END DO

      DET_AA = AA11 * AA22 - AA12 * AA12

      BVEC(1) = ( + AA22 * AB1 - AA12 * AB2 ) / DET_AA
      BVEC(2) = ( - AA12 * AB1 + AA11 * AB2 ) / DET_AA

    ELSEIF( mk > 3 )THEN

      ! --- Not Implemented ---

    END IF

    SUM1 = Zero
    DO i = 1, mk - 1

      Alpha_LS(i) = BVEC(i)

      SUM1 = SUM1 + BVEC(i)

    END DO

    Alpha_LS(mk) = One - SUM1

    RETURN
  END FUNCTION Alpha_LS


  FUNCTION ShiftVec( M, mk, Vec )

#if   defined( THORNADO_OMP_OL )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC   )
    !$ACC ROUTINE SEQ
#endif

    INTEGER,  INTENT(in) :: M, mk
    REAL(DP), INTENT(in) :: Vec(4,M)
    REAL(DP)             :: ShiftVec(4,M)

    INTEGER  :: i, j
    REAL(DP) :: VecTMP(4,M)

    DO j = 1, mk - 1
    DO i = 1, 4

      VecTMP(i,j) = Vec(i,j+1)

    END DO
    END DO

    DO j = 1, mk - 1
    DO i = 1, 4

      ShiftVec(i,j) = VecTMP(i,j)

    END DO
    END DO

    RETURN
  END FUNCTION ShiftVec


  SUBROUTINE ComputeConserved_TwoMoment &
    ( D, I_u_1, I_u_2, I_u_3, N, G_d_1, G_d_2, G_d_3, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33 )

#if   defined( THORNADO_OMP_OL )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC   )
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: D, I_u_1, I_u_2, I_u_3 ! --- Index Up
    REAL(DP), INTENT(out) :: N, G_d_1, G_d_2, G_d_3 ! --- Index Down
    REAL(DP), INTENT(in)  ::    V_u_1, V_u_2, V_u_3 ! --- Index Up
    REAL(DP), INTENT(in)  :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    REAL(DP) :: k_dd(3,3)

    k_dd = EddingtonTensorComponents_dd &
             ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    ! --- Conserved Number Density ---

    N = D + Gm_dd_11 * V_u_1 * I_u_1 &
          + Gm_dd_22 * V_u_2 * I_u_2 &
          + Gm_dd_33 * V_u_3 * I_u_3

    ! --- Conserved Number Flux Density (1) ---

    G_d_1 = Gm_dd_11 * I_u_1 &
              + (   V_u_1 * k_dd(1,1) &
                  + V_u_2 * k_dd(1,2) &
                  + V_u_3 * k_dd(1,3) ) * D

    ! --- Conserved Number Flux Density (2) ---

    G_d_2 = Gm_dd_22 * I_u_2 &
              + (   V_u_1 * k_dd(1,2) &
                  + V_u_2 * k_dd(2,2) &
                  + V_u_3 * k_dd(2,3) ) * D

    ! --- Conserved Number Flux Density (3) ---

    G_d_3 = Gm_dd_33 * I_u_3 &
              + (   V_u_1 * k_dd(1,3) &
                  + V_u_2 * k_dd(2,3) &
                  + V_u_3 * k_dd(3,3) ) * D

  END SUBROUTINE ComputeConserved_TwoMoment


  SUBROUTINE ComputeFromConserved_TwoMoment &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, CF, CR, PR )

    INTEGER,  INTENT(in)  :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)  :: &
      GX(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4), &
         1:nGF)
    REAL(DP), INTENT(in)  :: &
      CF(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4), &
         1:nCF)
    REAL(DP), INTENT(in)  :: &
      CR(1:nDOFZ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      PR(1:nDOFZ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4),1:nPR,1:nSpecies)

    INTEGER  :: &
      iZ1, iZ2, iZ3, iZ4, iS, iNodeZ, iNodeX
    REAL(DP) :: &
      PF(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4), &
         1:nPF)

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX

        CALL ComputePrimitive_Euler_NonRelativistic &
               ( CF(iNodeX,iZ2,iZ3,iZ4,iCF_D ), &
                 CF(iNodeX,iZ2,iZ3,iZ4,iCF_S1), &
                 CF(iNodeX,iZ2,iZ3,iZ4,iCF_S2), &
                 CF(iNodeX,iZ2,iZ3,iZ4,iCF_S3), &
                 CF(iNodeX,iZ2,iZ3,iZ4,iCF_E ), &
                 CF(iNodeX,iZ2,iZ3,iZ4,iCF_Ne), &
                 PF(iNodeX,iZ2,iZ3,iZ4,iPF_D ), &
                 PF(iNodeX,iZ2,iZ3,iZ4,iPF_V1), &
                 PF(iNodeX,iZ2,iZ3,iZ4,iPF_V2), &
                 PF(iNodeX,iZ2,iZ3,iZ4,iPF_V3), &
                 PF(iNodeX,iZ2,iZ3,iZ4,iPF_E ), &
                 PF(iNodeX,iZ2,iZ3,iZ4,iPF_Ne), &
                 GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_11), &
                 GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_22), &
                 GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_33) )

      END DO

    END DO
    END DO
    END DO

    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeZ = 1, nDOFZ

        iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1

        CALL ComputePrimitive_TwoMoment &
               ( CR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_N ,iS), &
                 CR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G1,iS), &
                 CR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G2,iS), &
                 CR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G3,iS), &
                 PR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS), &
                 PR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS), &
                 PR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS), &
                 PR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS), &
                 PF(iNodeX,iZ2,iZ3,iZ4,iPF_V1), &
                 PF(iNodeX,iZ2,iZ3,iZ4,iPF_V2), &
                 PF(iNodeX,iZ2,iZ3,iZ4,iPF_V3), &
                 GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_11), &
                 GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_22), &
                 GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_33) )

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE ComputeFromConserved_TwoMoment


  FUNCTION Flux_E &
    ( D, I_u_1, I_u_2, I_u_3, &
      V_u_1, V_u_2, V_u_3, &
      dV_u_1_dX1, dV_u_2_dX1, dV_u_3_dX1, &
      dV_u_1_dX2, dV_u_2_dX2, dV_u_3_dX2, &
      dV_u_1_dX3, dV_u_2_dX3, dV_u_3_dX3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      dGm_dd_11_dX1, dGm_dd_22_dX1, dGm_dd_33_dX1, &
      dGm_dd_11_dX2, dGm_dd_22_dX2, dGm_dd_33_dX2, &
      dGm_dd_11_dX3, dGm_dd_22_dX3, dGm_dd_33_dX3 )

#if   defined( THORNADO_OMP_OL )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC   )
    !$ACC ROUTINE SEQ
#endif

    REAL(DP)             :: Flux_E(4), Flux_E_Old(4)
    REAL(DP), INTENT(in) :: D, I_u_1, I_u_2, I_u_3
    REAL(DP), INTENT(in) :: V_u_1, V_u_2, V_u_3
    REAL(DP), INTENT(in) :: dV_u_1_dX1, dV_u_2_dX1, dV_u_3_dX1
    REAL(DP), INTENT(in) :: dV_u_1_dX2, dV_u_2_dX2, dV_u_3_dX2
    REAL(DP), INTENT(in) :: dV_u_1_dX3, dV_u_2_dX3, dV_u_3_dX3
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), INTENT(in) :: dGm_dd_11_dX1, dGm_dd_22_dX1, dGm_dd_33_dX1
    REAL(DP), INTENT(in) :: dGm_dd_11_dX2, dGm_dd_22_dX2, dGm_dd_33_dX2
    REAL(DP), INTENT(in) :: dGm_dd_11_dX3, dGm_dd_22_dX3, dGm_dd_33_dX3

    REAL(DP) :: k_uu(3,3), l_uuu(3,3,3)
    REAL(DP) :: VdotGradGm_dd_11
    REAL(DP) :: VdotGradGm_dd_22
    REAL(DP) :: VdotGradGm_dd_33

    VdotGradGm_dd_11 &
      = V_u_1 * dGm_dd_11_dX1 + V_u_2 * dGm_dd_11_dX2 + V_u_3 * dGm_dd_11_dX3

    VdotGradGm_dd_22 &
      = V_u_1 * dGm_dd_22_dX1 + V_u_2 * dGm_dd_22_dX2 + V_u_3 * dGm_dd_22_dX3

    VdotGradGm_dd_33 &
      = V_u_1 * dGm_dd_33_dX1 + V_u_2 * dGm_dd_33_dX2 + V_u_3 * dGm_dd_33_dX3

    ! --- Eddington Tensor Components ---

    CALL ComputeEddingtonTensorComponents_uu &
           ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, k_uu )

    ! --- Heat Flux Tensor Components ---

    CALL ComputeHeatFluxTensorComponents_uuu &
           ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, l_uuu )

    ! --- Number Density ---

    Flux_E(1) &
      = - (   Gm_dd_11 * (   k_uu(1,1) * dV_u_1_dX1 &
                           + k_uu(2,1) * dV_u_1_dX2 &
                           + k_uu(3,1) * dV_u_1_dX3 ) &
            + Gm_dd_22 * (   k_uu(1,2) * dV_u_2_dX1 &
                           + k_uu(2,2) * dV_u_2_dX2 &
                           + k_uu(3,2) * dV_u_2_dX3 ) &
            + Gm_dd_33 * (   k_uu(1,3) * dV_u_3_dX1 &
                           + k_uu(2,3) * dV_u_3_dX2 &
                           + k_uu(3,3) * dV_u_3_dX3 ) ) * D

    Flux_E(1) &
      = Flux_E(1) &
          - Half * (   k_uu(1,1) * VdotGradGm_dd_11 &
                     + k_uu(2,2) * VdotGradGm_dd_22 &
                     + k_uu(3,3) * VdotGradGm_dd_33 ) * D

    ! --- Number Flux Density (Component 1) ---

    Flux_E(2) &
      = - Gm_dd_11 * (   Gm_dd_11 * (   l_uuu(1,1,1) * dV_u_1_dX1 &
                                      + l_uuu(1,2,1) * dV_u_1_dX2 &
                                      + l_uuu(1,3,1) * dV_u_1_dX3 ) &
                       + Gm_dd_22 * (   l_uuu(2,1,1) * dV_u_2_dX1 &
                                      + l_uuu(2,2,1) * dV_u_2_dX2 &
                                      + l_uuu(2,3,1) * dV_u_2_dX3 ) &
                       + Gm_dd_33 * (   l_uuu(3,1,1) * dV_u_3_dX1 &
                                      + l_uuu(3,2,1) * dV_u_3_dX2 &
                                      + l_uuu(3,3,1) * dV_u_3_dX3 ) ) * D

    Flux_E(2) &
      = Flux_E(2) &
          - Half * Gm_dd_11 * (   l_uuu(1,1,1) * VdotGradGm_dd_11 &
                                + l_uuu(2,2,1) * VdotGradGm_dd_22 &
                                + l_uuu(3,3,1) * VdotGradGm_dd_33 ) * D

    ! --- Number Flux Density (Component 2) ---

    Flux_E(3) &
      = - Gm_dd_22 * (   Gm_dd_11 * (   l_uuu(1,1,2) * dV_u_1_dX1 &
                                      + l_uuu(1,2,2) * dV_u_1_dX2 &
                                      + l_uuu(1,3,2) * dV_u_1_dX3 ) &
                       + Gm_dd_22 * (   l_uuu(2,1,2) * dV_u_2_dX1 &
                                      + l_uuu(2,2,2) * dV_u_2_dX2 &
                                      + l_uuu(2,3,2) * dV_u_2_dX3 ) &
                       + Gm_dd_33 * (   l_uuu(3,1,2) * dV_u_3_dX1 &
                                      + l_uuu(3,2,2) * dV_u_3_dX2 &
                                      + l_uuu(3,3,2) * dV_u_3_dX3 ) ) * D

    Flux_E(3) &
      = Flux_E(3) &
          - Half * Gm_dd_22 * (   l_uuu(1,1,2) * VdotGradGm_dd_11 &
                                + l_uuu(2,2,2) * VdotGradGm_dd_22 &
                                + l_uuu(3,3,2) * VdotGradGm_dd_33 ) * D

    ! --- Number Flux Density (Component 3) ---

    Flux_E(4) &
      = - Gm_dd_33 * (   Gm_dd_11 * (   l_uuu(3,1,1) * dV_u_1_dX1 &
                                      + l_uuu(3,1,2) * dV_u_1_dX2 &
                                      + l_uuu(3,1,3) * dV_u_1_dX3 ) &
                       + Gm_dd_22 * (   l_uuu(3,2,1) * dV_u_2_dX1 &
                                      + l_uuu(3,2,2) * dV_u_2_dX2 &
                                      + l_uuu(3,2,3) * dV_u_2_dX3 ) &
                       + Gm_dd_33 * (   l_uuu(3,3,1) * dV_u_3_dX1 &
                                      + l_uuu(3,3,2) * dV_u_3_dX2 &
                                      + l_uuu(3,3,3) * dV_u_3_dX3 ) ) * D

    Flux_E(4) &
      = Flux_E(4) &
          - Half * Gm_dd_33 * (   l_uuu(1,1,3) * VdotGradGm_dd_11 &
                                + l_uuu(2,2,3) * VdotGradGm_dd_22 &
                                + l_uuu(3,3,3) * VdotGradGm_dd_33 ) * D

    RETURN
  END FUNCTION Flux_E


  FUNCTION Flux_X1 &
    ( D, I_u_1, I_u_2, I_u_3, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33 )

#if   defined( THORNADO_OMP_OL )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC   )
    !$ACC ROUTINE SEQ
#endif

    REAL(DP)             :: Flux_X1(4)
    REAL(DP), INTENT(in) :: D, I_u_1, I_u_2, I_u_3
    REAL(DP), INTENT(in) ::    V_u_1, V_u_2, V_u_3
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    REAL(DP) :: FF, EF, a, b
    REAL(DP) :: h_u_1, h_d_1, h_d_2, h_d_3
    REAL(DP) :: K_ud_11, K_ud_12, K_ud_13

    FF = FluxFactor( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    EF = EddingtonFactor( D, FF )

    a = Half * ( One - EF )
    b = Half * ( Three * EF - One )

    h_u_1 = I_u_1 / ( FF * D )

    h_d_1 = Gm_dd_11 * I_u_1 / ( FF * D )
    h_d_2 = Gm_dd_22 * I_u_2 / ( FF * D )
    h_d_3 = Gm_dd_33 * I_u_3 / ( FF * D )

    K_ud_11 = ( a + b * h_u_1 * h_d_1 ) * D
    K_ud_12 = (     b * h_u_1 * h_d_2 ) * D
    K_ud_13 = (     b * h_u_1 * h_d_3 ) * D

    Flux_X1(1) = I_u_1   + V_u_1 * D

    Flux_X1(2) = K_ud_11 + V_u_1 * Gm_dd_11 * I_u_1

    Flux_X1(3) = K_ud_12 + V_u_1 * Gm_dd_22 * I_u_2

    Flux_X1(4) = K_ud_13 + V_u_1 * Gm_dd_33 * I_u_3

    RETURN
  END FUNCTION Flux_X1


  FUNCTION Flux_X2 &
    ( D, I_u_1, I_u_2, I_u_3, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33 )

#if   defined( THORNADO_OMP_OL )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC   )
    !$ACC ROUTINE SEQ
#endif

    REAL(DP) :: Flux_X2(4)
    REAL(DP), INTENT(in) :: D, I_u_1, I_u_2, I_u_3
    REAL(DP), INTENT(in) ::    V_u_1, V_u_2, V_u_3
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    REAL(DP) :: FF, EF, a, b
    REAL(DP) :: h_u_2, h_d_1, h_d_2, h_d_3
    REAL(DP) :: K_ud_21, K_ud_22, K_ud_23

    FF = FluxFactor( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    EF = EddingtonFactor( D, FF )

    a = Half * ( One - EF )
    b = Half * ( Three * EF - One )

    h_u_2 = I_u_2 / ( FF * D )

    h_d_1 = Gm_dd_11 * I_u_1 / ( FF * D )
    h_d_2 = Gm_dd_22 * I_u_2 / ( FF * D )
    h_d_3 = Gm_dd_33 * I_u_3 / ( FF * D )

    K_ud_21 = (     b * h_u_2 * h_d_1 ) * D
    K_ud_22 = ( a + b * h_u_2 * h_d_2 ) * D
    K_ud_23 = (     b * h_u_2 * h_d_3 ) * D

    Flux_X2(1) = I_u_2   + V_u_2 * D

    Flux_X2(2) = K_ud_21 + V_u_2 * Gm_dd_11 * I_u_1

    Flux_X2(3) = K_ud_22 + V_u_2 * Gm_dd_22 * I_u_2

    Flux_X2(4) = K_ud_23 + V_u_2 * Gm_dd_33 * I_u_3

    RETURN
  END FUNCTION Flux_X2


  FUNCTION Flux_X3 &
    ( D, I_u_1, I_u_2, I_u_3, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33 )

#if   defined( THORNADO_OMP_OL )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC   )
    !$ACC ROUTINE SEQ
#endif

    REAL(DP) :: Flux_X3(4)
    REAL(DP), INTENT(in) :: D, I_u_1, I_u_2, I_u_3
    REAL(DP), INTENT(in) ::    V_u_1, V_u_2, V_u_3
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    REAL(DP) :: FF, EF, a, b
    REAL(DP) :: h_u_3, h_d_1, h_d_2, h_d_3
    REAL(DP) :: K_ud_31, K_ud_32, K_ud_33

    FF = FluxFactor( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    EF = EddingtonFactor( D, FF )

    a = Half * ( One - EF )
    b = Half * ( Three * EF - One )

    h_u_3 = I_u_3 / ( FF * D )

    h_d_1 = Gm_dd_11 * I_u_1 / ( FF * D )
    h_d_2 = Gm_dd_22 * I_u_2 / ( FF * D )
    h_d_3 = Gm_dd_33 * I_u_3 / ( FF * D )

    K_ud_31 = (     b * h_u_3 * h_d_1 ) * D
    K_ud_32 = (     b * h_u_3 * h_d_2 ) * D
    K_ud_33 = ( a + b * h_u_3 * h_d_3 ) * D

    Flux_X3(1) = I_u_3   + V_u_3 * D

    Flux_X3(2) = K_ud_31 + V_u_3 * Gm_dd_11 * I_u_1

    Flux_X3(3) = K_ud_32 + V_u_3 * Gm_dd_22 * I_u_2

    Flux_X3(4) = K_ud_33 + V_u_3 * Gm_dd_33 * I_u_3

    RETURN
  END FUNCTION Flux_X3


  SUBROUTINE ComputeEddingtonTensorComponents_uu &
    ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, k_uu )

#if   defined( THORNADO_OMP_OL )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC   )
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: &
      D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), INTENT(out) :: &
      k_uu(3,3)

    INTEGER  :: i, j
    REAL(DP) :: FF, EF, a, b
    REAL(DP) :: h_u(3), Gm_uu(3,3)

    FF = FluxFactor( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    EF = EddingtonFactor( D, FF )

    a = Half * ( One - EF )
    b = Half * ( Three * EF - One )

    h_u = [ I_u_1, I_u_2, I_u_3 ] / ( FF * D )

    Gm_uu      = Zero
    Gm_uu(1,1) = One / Gm_dd_11
    Gm_uu(2,2) = One / Gm_dd_22
    Gm_uu(3,3) = One / Gm_dd_33

    DO j = 1, 3
    DO i = 1, 3

      k_uu(i,j) = a * Gm_uu(i,j) + b * h_u(i) * h_u(j)

    END DO
    END DO

  END SUBROUTINE ComputeEddingtonTensorComponents_uu


  SUBROUTINE ComputeEddingtonTensorComponents_dd &
    ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33 )

#if   defined( THORNADO_OMP_OL )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC   )
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: &
      D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), INTENT(out) :: &
      k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33

    REAL(DP) :: FF, EF, a, b
    REAL(DP) :: h_d_1, h_d_2, h_d_3

    FF = FluxFactor( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    EF = EddingtonFactor( D, FF )

    a = Half * ( One - EF )
    b = Half * ( Three * EF - One )

    h_d_1 = Gm_dd_11 * I_u_1 / ( FF * D )
    h_d_2 = Gm_dd_22 * I_u_2 / ( FF * D )
    h_d_3 = Gm_dd_33 * I_u_3 / ( FF * D )

    ! --- Diagonal Eddington Tensor Components ---

    k_dd_11 = a * Gm_dd_11 + b * h_d_1 * h_d_1
    k_dd_22 = a * Gm_dd_22 + b * h_d_2 * h_d_2
    k_dd_33 = a * Gm_dd_33 + b * h_d_3 * h_d_3

    ! --- Off-Diagonal Eddington Tensor Components ---

    k_dd_12 = b * h_d_1 * h_d_2
    k_dd_13 = b * h_d_1 * h_d_3
    k_dd_23 = b * h_d_2 * h_d_3

  END SUBROUTINE ComputeEddingtonTensorComponents_dd


  FUNCTION EddingtonTensorComponents_dd &
    ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

#if   defined( THORNADO_OMP_OL )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC   )
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: &
      D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP)              :: &
      EddingtonTensorComponents_dd(3,3)

    INTEGER  :: i, j
    REAL(DP) :: FF, EF, a, b
    REAL(DP) :: h_d(3), Gm_dd(3,3)

    FF = FluxFactor( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    EF = EddingtonFactor( D, FF )

    a = Half * ( One - EF )
    b = Half * ( Three * EF - One )

    h_d(1) = Gm_dd_11 * I_u_1 / ( FF * D )
    h_d(2) = Gm_dd_22 * I_u_2 / ( FF * D )
    h_d(3) = Gm_dd_33 * I_u_3 / ( FF * D )

    Gm_dd = Zero
    Gm_dd(1,1) = Gm_dd_11
    Gm_dd(2,2) = Gm_dd_22
    Gm_dd(3,3) = Gm_dd_33

    DO j = 1, 3
    DO i = 1, 3

      EddingtonTensorComponents_dd(i,j) &
        = a * Gm_dd(i,j) + b * h_d(i) * h_d(j)

    END DO
    END DO

    RETURN
  END FUNCTION EddingtonTensorComponents_dd


  SUBROUTINE ComputeEddingtonTensorComponents_ud &
    ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      k_ud_11, k_ud_12, k_ud_13, k_ud_22, k_ud_23, k_ud_33 )

    REAL(DP), INTENT(in)  :: &
      D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), INTENT(out) :: &
      k_ud_11, k_ud_12, k_ud_13, k_ud_22, k_ud_23, k_ud_33

    REAL(DP) :: FF, EF, a, b
    REAL(DP) :: h_u_1, h_u_2, h_u_3
    REAL(DP) :: h_d_1, h_d_2, h_d_3

    FF = FluxFactor( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    EF = EddingtonFactor( D, FF )

    a = Half * ( One - EF )
    b = Half * ( Three * EF - One )

    h_u_1 = I_u_1 / ( FF * D )
    h_u_2 = I_u_2 / ( FF * D )
    h_u_3 = I_u_3 / ( FF * D )

    h_d_1 = Gm_dd_11 * h_u_1
    h_d_2 = Gm_dd_22 * h_u_2
    h_d_3 = Gm_dd_33 * h_u_3

    ! --- Diagonal Eddington Tensor Components ---

    k_ud_11 = a + b * h_u_1 * h_d_1
    k_ud_22 = a + b * h_u_2 * h_d_2
    k_ud_33 = a + b * h_u_3 * h_d_3

    ! --- Off-Diagonal Eddington Tensor Components ---

    k_ud_12 = b * h_u_1 * h_d_2
    k_ud_13 = b * h_u_1 * h_d_3
    k_ud_23 = b * h_u_2 * h_d_3

  END SUBROUTINE ComputeEddingtonTensorComponents_ud


  SUBROUTINE ComputeHeatFluxTensorComponents_uuu &
    ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, l_uuu )

#if   defined( THORNADO_OMP_OL )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC   )
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: &
      D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), INTENT(out) :: &
      l_uuu(3,3,3)

    INTEGER  :: i, j, k
    REAL(DP) :: FF, HF, a, b
    REAL(DP) :: h_u(3), Gm_uu(3,3)

    FF = FluxFactor( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    HF = HeatFluxFactor( D, FF )

    a = Half * ( FF - HF )
    b = Half * ( Five * HF - Three * FF )

    h_u = [ I_u_1, I_u_2, I_u_3 ] / ( FF * D )

    Gm_uu      = Zero
    Gm_uu(1,1) = One / Gm_dd_11
    Gm_uu(2,2) = One / Gm_dd_22
    Gm_uu(3,3) = One / Gm_dd_33

    DO k = 1, 3
    DO j = 1, 3
    DO i = 1, 3

      l_uuu(i,j,k) &
        = a * (   h_u(i) * Gm_uu(j,k) &
                + h_u(j) * Gm_uu(i,k) &
                + h_u(k) * Gm_uu(i,j) ) &
          + b * h_u(i) * h_u(j) * h_u(k)

    END DO
    END DO
    END DO

  END SUBROUTINE ComputeHeatFluxTensorComponents_uuu


  SUBROUTINE ComputeHeatFluxTensorComponents_uud &
    ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      l_uud_111, l_uud_121, l_uud_131, l_uud_221, l_uud_231, l_uud_331, &
      l_uud_222, l_uud_223, l_uud_333, l_uud_332 )

    REAL(DP), INTENT(in)  :: &
      D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), INTENT(out) :: &
      l_uud_111, l_uud_121, l_uud_131, &
                 l_uud_221, l_uud_231, &
                            l_uud_331, &
      l_uud_222, l_uud_223, &
      l_uud_333, l_uud_332

    REAL(DP) :: FF, HF, a, b
    REAL(DP) :: h_u_1, h_u_2, h_u_3
    REAL(DP) :: h_d_1, h_d_2, h_d_3

    FF = FluxFactor( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    HF = HeatFluxFactor( D, FF )

    a = Half * ( FF - HF )
    b = Half * ( Five * HF - Three * FF )

    h_u_1 = I_u_1 / ( FF * D )
    h_u_2 = I_u_2 / ( FF * D )
    h_u_3 = I_u_3 / ( FF * D )

    h_d_1 = Gm_dd_11 * h_u_1
    h_d_2 = Gm_dd_22 * h_u_2
    h_d_3 = Gm_dd_33 * h_u_3

    ! --- Diagonal Heat Flux Tensor Components ---

    l_uud_111 &
      = a * ( h_u_1 + Two * h_u_1 ) + b * h_u_1 * h_u_1 * h_d_1
    l_uud_222 &
      = a * ( h_u_2 + Two * h_u_2 ) + b * h_u_2 * h_u_2 * h_d_2
    l_uud_333 &
      = a * ( h_u_3 + Two * h_u_3 ) + b * h_u_3 * h_u_3 * h_d_3

    ! --- Off-Diagonal Heat Flux Tensor Components ---

    l_uud_121 &
      = a * h_u_2 + b * h_u_1 * h_u_2 * h_d_1
    l_uud_131 &
      = a * h_u_3 + b * h_u_1 * h_u_3 * h_d_1
    l_uud_221 &
      = a * h_d_1 / Gm_dd_22 + b * h_u_2 * h_u_2 * h_d_1
    l_uud_231 &
      =                        b * h_u_2 * h_u_3 * h_d_1
    l_uud_331 &
      = a * h_d_1 / Gm_dd_33 + b * h_u_3 * h_u_3 * h_d_1
    l_uud_223 &
      = a * h_d_3 / Gm_dd_22 + b * h_u_2 * h_u_2 * h_d_3
    l_uud_332 &
      = a * h_d_2 / Gm_dd_33 + b * h_u_3 * h_u_3 * h_d_2

  END SUBROUTINE ComputeHeatFluxTensorComponents_uud


  SUBROUTINE ComputeHeatFluxTensorComponents_udd &
    ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      l_udd_111, l_udd_121, l_udd_131, l_udd_221, l_udd_231, l_udd_331, &
      l_udd_222, l_udd_223, l_udd_333, l_udd_332 )

    REAL(DP), INTENT(in)  :: &
      D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), INTENT(out) :: &
      l_udd_111, l_udd_121, l_udd_131, &
                 l_udd_221, l_udd_231, &
                            l_udd_331, &
      l_udd_222, l_udd_223, &
      l_udd_333, l_udd_332

    REAL(DP) :: FF, HF, a, b
    REAL(DP) :: h_u_1, h_u_2, h_u_3
    REAL(DP) :: h_d_1, h_d_2, h_d_3

    FF = FluxFactor( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    HF = HeatFluxFactor( D, FF )

    a = Half * ( FF - HF )
    b = Half * ( Five * HF - Three * FF )

    h_u_1 = I_u_1 / ( FF * D )
    h_u_2 = I_u_2 / ( FF * D )
    h_u_3 = I_u_3 / ( FF * D )

    h_d_1 = Gm_dd_11 * h_u_1
    h_d_2 = Gm_dd_22 * h_u_2
    h_d_3 = Gm_dd_33 * h_u_3

    ! --- Diagonal Heat Flux Tensor Components ---

    l_udd_111 &
      = a * ( h_d_1 + Two * h_d_1 ) + b * h_u_1 * h_d_1 * h_d_1
    l_udd_222 &
      = a * ( h_d_2 + Two * h_d_2 ) + b * h_u_2 * h_d_2 * h_d_2
    l_udd_333 &
      = a * ( h_d_3 + Two * h_d_3 ) + b * h_u_3 * h_d_3 * h_d_3

    ! --- Off-Diagonal Heat Flux Tensor Components ---

    l_udd_121 &
      = a * h_d_2 + b * h_u_1 * h_d_2 * h_d_1
    l_udd_131 &
      = a * h_d_3 + b * h_u_1 * h_d_3 * h_d_1
    l_udd_221 &
      = a * h_d_1 + b * h_u_2 * h_d_2 * h_d_1
    l_udd_231 &
      =             b * h_u_2 * h_d_3 * h_d_1
    l_udd_331 &
      = a * h_d_1 + b * h_u_3 * h_d_3 * h_d_1
    l_udd_223 &
      = a * h_d_3 + b * h_u_2 * h_d_2 * h_d_3
    l_udd_332 &
      = a * h_d_2 + b * h_u_3 * h_d_3 * h_d_2

  END SUBROUTINE ComputeHeatFluxTensorComponents_udd


  FUNCTION NumericalFlux_LLF( u_L, u_R, Flux_L, Flux_R, alpha )

#if   defined( THORNADO_OMP_OL )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC   )
    !$ACC ROUTINE SEQ
#endif

    REAL(DP)             :: NumericalFlux_LLF
    REAL(DP), INTENT(in) :: u_L, u_R, flux_L, flux_R, alpha

    NumericalFlux_LLF &
      = Half * ( flux_L + flux_R - alpha * ( u_R - u_L ) )

    RETURN
  END FUNCTION NumericalFlux_LLF


END MODULE TwoMoment_UtilitiesModule_OrderV
