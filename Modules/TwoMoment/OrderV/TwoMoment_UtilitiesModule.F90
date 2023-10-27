MODULE TwoMoment_UtilitiesModule

  USE KindModule, ONLY: &
    DP, Zero, Half, One, Two, Three, Five, Third, &
    SqrtTiny, FourPi
  USE UnitsModule, ONLY: &
    UnitsActive, &
    SpeedOfLight, &
    PlanckConstant, &
    MeV
  USE QuadratureModule, ONLY: &
    GetQuadrature
  USE ProgramHeaderModule, ONLY: &
    nDOFZ, nDOFX, nDOFE, &
    nNodes, nDimsX
  USE ReferenceElementModuleE, ONLY: &
    WeightsE
  USE ReferenceElementModuleX, ONLY: &
    nDOFX_X1, nDOFX_X2, nDOFX_X3, &
    WeightsX_q, &
    WeightsX_X1, &
    WeightsX_X2, &
    WeightsX_X3
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    dLXdX1_q, LX_X1_Dn, LX_X1_Up, &
    dLXdX2_q, LX_X2_Dn, LX_X2_Up, &
    dLXdX3_q, LX_X3_Dn, LX_X3_Up
  USE MeshModule, ONLY: &
    MeshX, &
    MeshE, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    nGF, &
    iGF_h_1, &
    iGF_h_2, &
    iGF_h_3, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_SqrtGm
  USE FluidFieldsModule, ONLY: &
    nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne
  USE Euler_UtilitiesModule_NonRelativistic, ONLY: &
    ComputePrimitive_Euler_NonRelativistic
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3, &
    nPR, iPR_D, iPR_I1, iPR_I2, iPR_I3, &
    nAR, iAR_F, iAR_K , iAR_Q, &
    nGR, iGR_N, &
         iGR_D, iGR_I1, iGR_I2, iGR_I3, &
         iGR_J, iGR_H1, iGR_H2, iGR_H3, &
         iGR_RMS, iGR_F, iGR_K, iGR_Q
  USE TwoMoment_ClosureModule, ONLY: &
    FluxFactor, &
    EddingtonFactor, &
    HeatFluxFactor
  USE TwoMoment_TimersModule, ONLY: &
    TimersStart, &
    TimersStop, &
    Timer_Streaming_NumericalFlux_InOut, &
    Timer_Streaming_NumericalFlux_RHS, &
    Timer_Streaming_NumericalFlux_LS, &
    Timer_Streaming_NumericalFlux_Update, &
    Timer_Streaming_LinearAlgebra
  USE LinearAlgebraModule, ONLY: &
    MatrixMatrixMultiply, &
    EigenvaluesSymmetric3

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeComputePrimitive_TwoMoment
  PUBLIC :: FinalizeComputePrimitive_TwoMoment
  PUBLIC :: ComputePrimitive_TwoMoment
  PUBLIC :: ComputePrimitive_TwoMoment_Scalar
  PUBLIC :: ComputeConserved_TwoMoment
  PUBLIC :: ComputeFromConserved_TwoMoment
  PUBLIC :: ComputeTimeStep_TwoMoment
  PUBLIC :: ComputeTimeStep_TwoMoment_Realizability
  PUBLIC :: Flux_E
  PUBLIC :: Flux_X1
  PUBLIC :: Flux_X2
  PUBLIC :: Flux_X3
  PUBLIC :: ComputeEddingtonTensorComponents_uu
  PUBLIC :: ComputeEddingtonTensorComponents_dd
  PUBLIC :: EddingtonTensorComponents_dd
  PUBLIC :: ComputeEddingtonTensorComponents_ud
  PUBLIC :: ComputeHeatFluxTensorComponents_udd
  PUBLIC :: NumericalFlux_LLF
  PUBLIC :: ComputeWeakDerivatives_X1
  PUBLIC :: ComputeWeakDerivatives_X2
  PUBLIC :: ComputeWeakDerivatives_X3
  PUBLIC :: FaceVelocity_X1
  PUBLIC :: FaceVelocity_X2
  PUBLIC :: FaceVelocity_X3

  INTERFACE ComputePrimitive_TwoMoment
    MODULE PROCEDURE ComputePrimitive_TwoMoment_Scalar_Richardson
    MODULE PROCEDURE ComputePrimitive_TwoMoment_Vector_Richardson
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
    !$ACC ENTER DATA &
    !$ACC COPYIN( ITERATE ) &
    !$ACC CREATE( FVEC, GVEC, CVEC, UVEC, &
    !$ACC         FVECm, GVECm, Alpha, nIterations )
#endif

    ! --- Initial Guess ---

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR &
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
      !$OMP PRIVATE( iX, A_d_1, A_d_2, A_d_3, k_dd, DET, LMAT, SUM1 )
#elif defined( THORNADO_OACC   )
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC PRIVATE( iX, A_d_1, A_d_2, A_d_3, k_dd, DET, LMAT, SUM1 ) &
      !$ACC PRESENT( ITERATE, UVEC, CVEC, GVEC, FVEC, GVECm, FVECm, &
      !$ACC          PositionIndexZ, D, I_u_1, I_u_2, I_u_3, &
      !$ACC          Gm_dd_11, Gm_dd_22, Gm_dd_33, V_u_1, V_u_2, V_u_3 )
#elif defined( THORNADO_OMP    )
      !$OMP PARALLEL DO &
      !$OMP PRIVATE( iX, A_d_1, A_d_2, A_d_3, k_dd, DET, LMAT, SUM1 )
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
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
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
      !$ACC PARALLEL LOOP GANG VECTOR &
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

          CONVERGED = SQRT( SUM( FVECm(:,iZ)**2 ) ) &
                        <= Rtol * SQRT( SUM( CVEC(:,iZ)**2 ) )

          nIterations(iZ) = k
          IF ( CONVERGED ) THEN
            ITERATE(iZ) = .FALSE.
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
      !$ACC UPDATE HOST( ITERATE )
      !$ACC WAIT
#endif

      CALL TimersStop( Timer_Streaming_NumericalFlux_Update )

    END DO

    CALL TimersStart( Timer_Streaming_NumericalFlux_InOut )

    IF( PRESENT( nIterations_Option ) ) THEN
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR &
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


  SUBROUTINE ComputePrimitive_TwoMoment_Vector_Richardson &
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
    REAL(DP) :: SUM1, k_dd(3,3)
    REAL(DP) :: vMag, Omega, vI, vK
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
    !$ACC ENTER DATA &
    !$ACC COPYIN( ITERATE ) &
    !$ACC CREATE( FVEC, GVEC, CVEC, UVEC, &
    !$ACC         FVECm, GVECm, Alpha, nIterations )
#endif

    ! --- Initial Guess ---

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( iX )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( iX ) &
    !$ACC PRESENT( CVEC, N, G_d_1, G_d_2, G_d_3, &
    !$ACC          PositionIndexZ, D, I_u_1, I_u_2, I_u_3 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( iX )
#endif
    DO iZ = 1, nZ
      CVEC(iCR_N ,iZ) = N    (iZ)
      CVEC(iCR_G1,iZ) = G_d_1(iZ)
      CVEC(iCR_G2,iZ) = G_d_2(iZ)
      CVEC(iCR_G3,iZ) = G_d_3(iZ)

      iX = PositionIndexZ(iZ)

      D    (iZ) = N(iZ)
      I_u_1(iZ) = G_d_1(iZ) / Gm_dd_11(iX)
      I_u_2(iZ) = G_d_2(iZ) / Gm_dd_22(iX)
      I_u_3(iZ) = G_d_3(iZ) / Gm_dd_33(iX)
    END DO

    CALL TimersStop( Timer_Streaming_NumericalFlux_InOut )

    k = 0
    DO WHILE( ANY( ITERATE ) .AND. k < MaxIterations )

      k = k + 1
      Mk = MIN( M, k )

      CALL TimersStart( Timer_Streaming_NumericalFlux_RHS )

#if   defined( THORNADO_OMP_OL )
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
      !$OMP PRIVATE( iX, k_dd, vMag, Omega, vI, vK )
#elif defined( THORNADO_OACC   )
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC PRIVATE( iX, k_dd, vMag, Omega, vI, vK ) &
      !$ACC PRESENT( ITERATE, UVEC, CVEC, GVEC, FVEC, GVECm, FVECm, &
      !$ACC          PositionIndexZ, D, I_u_1, I_u_2, I_u_3, &
      !$ACC          Gm_dd_11, Gm_dd_22, Gm_dd_33, V_u_1, V_u_2, V_u_3 )
#elif defined( THORNADO_OMP    )
      !$OMP PARALLEL DO &
      !$OMP PRIVATE( iX, k_dd, vMag, Omega, vI, vK )
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

          vMag = SQRT(   V_u_1(iX) * Gm_dd_11(iX) * V_u_1(iX) &
                       + V_u_2(iX) * Gm_dd_22(iX) * V_u_2(iX) &
                       + V_u_3(iX) * Gm_dd_33(iX) * V_u_3(iX) )

          Omega = One / ( One + vMag )

          vI =   V_u_1(iX) * UVEC(iPR_I1,iZ) &
               + V_u_2(iX) * UVEC(iPR_I2,iZ) &
               + V_u_3(iX) * UVEC(iPR_I3,iZ)

          GVECm(1,iZ) = (One - Omega) * UVEC(iPR_D,iZ) &
                        + Omega * ( CVEC(iCR_N,iZ) - vI )

          DO j = 1, 3

            vK =   V_u_1(iX) * k_dd(j,1) &
                 + V_u_2(iX) * k_dd(j,2) &
                 + V_u_3(iX) * k_dd(j,3)

            GVECm(j+1,iZ) = (One - Omega) * UVEC(j+1,iZ) &
                            + Omega * ( CVEC(j+1,iZ) - vK * UVEC(iPR_D,iZ) )

          END DO

          DO i = 1, 4

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
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
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
      !$ACC PARALLEL LOOP GANG VECTOR &
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

          CONVERGED = SQRT( SUM( FVECm(:,iZ)**2 ) ) <= &
                                 Rtol * SQRT( SUM( CVEC(:,iZ)**2 ) )

          nIterations(iZ) = k
          IF ( CONVERGED ) THEN
            ITERATE(iZ) = .FALSE.
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
      !$ACC UPDATE HOST( ITERATE )
      !$ACC WAIT
#endif

      CALL TimersStop( Timer_Streaming_NumericalFlux_Update )

    END DO

    CALL TimersStart( Timer_Streaming_NumericalFlux_InOut )

    IF( PRESENT( nIterations_Option ) ) THEN
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR &
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

  END SUBROUTINE ComputePrimitive_TwoMoment_Vector_Richardson


  SUBROUTINE Alpha_LS_Scalar &
    ( M, Mk, Fm, F, Alpha )

#if   defined( THORNADO_OMP_OL )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC   )
    !$ACC ROUTINE SEQ
#endif

    INTEGER,                  INTENT(in)    :: M, Mk
    REAL(DP), DIMENSION(4),   INTENT(inout) :: Fm
    REAL(DP), DIMENSION(4,M), INTENT(inout) :: F
    REAL(DP), DIMENSION(M),   INTENT(inout) :: Alpha

    REAL(DP) :: AA11, AA12, AA22, AB1, AB2, DET_AA
    REAL(DP) :: A1, A2, B
    INTEGER  :: iP

    IF ( Mk > 1 ) THEN

      IF ( Mk == 2 ) THEN

        AA11 = Zero
        AB1 = Zero

        DO iP = 1, 4

          A1 = F(iP,1) - Fm(iP)
          B  = - Fm(iP)

          AA11 = AA11 + A1 * A1
          AB1  = AB1  + A1 * B

        END DO

        Alpha(1) = AB1 / ( AA11 + SqrtTiny )
        Alpha(2) = One - Alpha(1)

      ELSE IF ( Mk == 3 ) THEN

        AA11 = Zero
        AA12 = Zero
        AA22 = Zero
        AB1  = Zero
        AB2  = Zero

        DO iP = 1, 4

          A1 = F(iP,1) - Fm(iP)
          A2 = F(iP,2) - Fm(iP)
          B  = - Fm(iP)

          AA11 = AA11 + A1 * A1
          AA12 = AA12 + A1 * A2
          AA22 = AA22 + A2 * A2

          AB1  = AB1  + A1 * B
          AB2  = AB2  + A2 * B

        END DO

        DET_AA = AA11*AA22 - AA12*AA12

        Alpha(1) = ( + AA22 * AB1 - AA12 * AB2 ) / DET_AA
        Alpha(2) = ( - AA12 * AB1 + AA11 * AB2 ) / DET_AA
        Alpha(3) = One - Alpha(1) - Alpha(2)

      ELSE IF ( Mk > 3 ) THEN

        ! --- Not Implemented ---

      END IF

    END IF

  END SUBROUTINE Alpha_LS_Scalar


  SUBROUTINE Alpha_LS_Vector &
    ( MASK, nZ, M, Mk, Fm, F, Alpha )

    LOGICAL,  DIMENSION(:),     INTENT(in)    :: MASK
    INTEGER,                    INTENT(in)    :: nZ, M, Mk
    REAL(DP), DIMENSION(:,:),   INTENT(inout) :: Fm
    REAL(DP), DIMENSION(:,:,:), INTENT(inout) :: F
    REAL(DP), DIMENSION(:,:),   INTENT(inout) :: Alpha

    REAL(DP) :: AA11, AA12, AA22, AB1, AB2, DET_AA
    REAL(DP) :: A1, A2, B
    INTEGER  :: iZ, iP

    IF ( Mk > 1 ) THEN

      IF ( Mk == 2 ) THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
        !$OMP PRIVATE( AA11, AB1, A1, B, iP )
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR &
        !$ACC PRIVATE( AA11, AB1, A1, B, iP ) &
        !$ACC PRESENT( MASK, Alpha, F, Fm )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO &
        !$OMP PRIVATE( AA11, AB1, A1, B, iP )
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

            Alpha(1,iZ) = AB1 / ( AA11 + SqrtTiny )
            Alpha(2,iZ) = One - Alpha(1,iZ)

          END IF
        END DO

      ELSE IF ( Mk == 3 ) THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
        !$OMP PRIVATE( AA11, AA12, AA22, AB1, AB2, DET_AA, A1, A2, B, iP )
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR &
        !$ACC PRIVATE( AA11, AA12, AA22, AB1, AB2, DET_AA, A1, A2, B, iP ) &
        !$ACC PRESENT( MASK, Alpha, F, Fm )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO &
        !$OMP PRIVATE( AA11, AA12, AA22, AB1, AB2, DET_AA, A1, A2, B, iP )
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

      IF( SQRT( SUM( FVECm**2 ) ) <= Rtol * SQRT( SUM( CVEC**2 ) ) )THEN

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


  SUBROUTINE ComputePrimitive_TwoMoment_Scalar_Richardson &
    ( N, G_d_1, G_d_2, G_d_3, D, I_u_1, I_u_2, I_u_3, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, nIterations_Option )

#if   defined( THORNADO_OMP_OL )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC   )
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: N, G_d_1, G_d_2, G_d_3
    REAL(DP), INTENT(out) :: D, I_u_1, I_u_2, I_u_3
    REAL(DP), INTENT(in)  :: V_u_1, V_u_2, V_u_3
    REAL(DP), INTENT(in)  :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    INTEGER,  INTENT(out), OPTIONAL :: nIterations_Option

    ! --- Parameters ---

    INTEGER,  PARAMETER :: M = 2
    INTEGER,  PARAMETER :: MaxIterations = 100
    REAL(DP), PARAMETER :: Rtol = 1.0d-08

    ! --- Local Variables ---

    INTEGER  :: k, Mk, iM, i, j
    REAL(DP) :: SUM1, k_dd(3,3)
    REAL(DP) :: vMag, Omega, vI, vK
    REAL(DP) :: I_d_1, I_d_2, I_d_3
    REAL(DP) :: UVEC(4), CVEC(4)
    REAL(DP) :: GVEC(4,M), GVECm(4)
    REAL(DP) :: FVEC(4,M), FVECm(4)
    REAL(DP) :: Alpha(M)

    LOGICAL  :: CONVERGED

    ! --- Initial Guess ---

    CVEC = [ N, G_d_1, G_d_2, G_d_3 ]

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

      k = k + 1
      Mk = MIN( M, k )

      UVEC = [ D, I_d_1, I_d_2, I_d_3 ]

      k_dd = EddingtonTensorComponents_dd &
               ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

      vMag = SQRT(   V_u_1 * Gm_dd_11 * V_u_1 &
                   + V_u_2 * Gm_dd_22 * V_u_2 &
                   + V_u_3 * Gm_dd_33 * V_u_3 )

      Omega = One / ( One + vMag )

      vI =   V_u_1 * UVEC(2) &
           + V_u_2 * UVEC(3) &
           + V_u_3 * UVEC(4)

      GVECm(1) = (One - Omega) *   UVEC(iPR_D) &
                      + Omega  * ( CVEC(iCR_N) - vI )

      DO j = 1, 3

        vK =   V_u_1 * k_dd(j,1) &
             + V_u_2 * k_dd(j,2) &
             + V_u_3 * k_dd(j,3)

        GVECm(j+1) = (One - Omega) *   UVEC(j+1) &
                          + Omega  * ( CVEC(j+1) - vK * UVEC(iPR_D) )

      END DO

      DO i = 1, 4

        FVECm(i) = GVECm(i) - UVEC(i)

        GVEC(i,Mk) = GVECm(i)
        FVEC(i,Mk) = FVECm(i)

      END DO


      IF ( Mk > 1 ) THEN

        CALL Alpha_LS_Scalar( M, Mk, FVECm, FVEC, Alpha )

        DO i = 1, 4
          SUM1 = Zero
          DO iM = 1, Mk
            SUM1 = SUM1 + GVEC(i,iM) * Alpha(iM)
          END DO
          GVECm(i) = SUM1
          FVECm(i) = GVECm(i) - UVEC(i)
        END DO

      END IF

      D     = GVECm(1)
      I_d_1 = GVECm(2); I_u_1 = I_d_1 / Gm_dd_11
      I_d_2 = GVECm(3); I_u_2 = I_d_2 / Gm_dd_22
      I_d_3 = GVECm(4); I_u_3 = I_d_3 / Gm_dd_33

      CONVERGED = SQRT( SUM( FVECm**2 ) ) <= &
                             Rtol * SQRT( SUM( CVEC**2 ) )

      IF ( Mk == M .AND. .NOT. CONVERGED ) THEN

        FVEC = ShiftVec( M, mk, FVEC )
        GVEC = ShiftVec( M, mk, GVEC )

      END IF

    END DO

    IF( PRESENT( nIterations_Option ) ) THEN

      nIterations_Option = k

    END IF

  END SUBROUTINE ComputePrimitive_TwoMoment_Scalar_Richardson


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
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, CF, CR, PR, AR, GR )

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
    REAL(DP), INTENT(out) :: &
      AR(1:nDOFZ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4),1:nAR,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      GR(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4), &
         1:nGR,1:nSpecies)

    INTEGER  :: &
      iZ1, iZ2, iZ3, iZ4, iS, iNodeZ, iNodeX

    REAL(DP), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: PF

    ALLOCATE( PF(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nPF) )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE FROM( GX, CF, CR )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE HOST( GX, CF, CR )
#endif

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

    CALL ComputeAuxiliary_TwoMoment &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, PR, AR )           

    CALL ComputeGray_TwoMoment &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, PF, CR, PR, AR, GR )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE TO( PR, AR, GR )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE DEVICE( PR, AR, GR )
#endif

  END SUBROUTINE ComputeFromConserved_TwoMoment


  SUBROUTINE ComputeAuxiliary_TwoMoment &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, PR, AR )

    INTEGER,  INTENT(in)  :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)  :: &
      GX(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4), &
         1:nGF)
    REAL(DP), INTENT(in)  :: &
      PR(1:nDOFZ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4),1:nPR,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      AR(1:nDOFZ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4),1:nAR,1:nSpecies)

    INTEGER  :: iZ1, iZ2, iZ3, iZ4, iS, iNodeZ, iNodeE, iNodeX
    REAL(DP) :: FF

    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeX = 1, nDOFX
      DO iNodeE = 1, nDOFE

        iNodeZ = (iNodeX-1) * nDOFE + iNodeE

        FF = FluxFactor &
               ( PR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS), &
                 PR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS), &
                 PR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS), &
                 PR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS), &
                 GX(iNodeX    ,iZ2,iZ3,iZ4,iGF_Gm_dd_11), &
                 GX(iNodeX    ,iZ2,iZ3,iZ4,iGF_Gm_dd_22), &
                 GX(iNodeX    ,iZ2,iZ3,iZ4,iGF_Gm_dd_33) )

        AR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iAR_F,iS) &
          = FF

        AR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iAR_K,iS) &
          = EddingtonFactor( PR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D,iS), FF )

        AR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iAR_Q,iS) &
          = HeatFluxFactor ( PR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D,iS), FF )

      END DO
      END DO

    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE ComputeAuxiliary_TwoMoment


  SUBROUTINE ComputeGray_TwoMoment &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, PF, CR, PR, AR, GR )

    INTEGER,  INTENT(in)  :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)  :: &
      GX(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4), &
         1:nGF)
    REAL(DP), INTENT(in)  :: &
      PF(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4), &
         1:nPF)
    REAL(DP), INTENT(in)  :: &
      CR(1:nDOFZ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(in)  :: &
      PR(1:nDOFZ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4),1:nPR,1:nSpecies)
    REAL(DP), INTENT(in)  :: &
      AR(1:nDOFZ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4),1:nAR,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      GR(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4), &
         1:nGR,1:nSpecies)

    INTEGER  :: iZ1, iZ2, iZ3, iZ4, iGR, iS, iNodeZ, iNodeE, iNodeX
    REAL(DP) :: hc3, E_0, E, RMS_Int3, RMS_Int5
    REAL(DP) :: W2(1:nDOFE,iZ_B0(1):iZ_E0(1))
    REAL(DP) :: W3(1:nDOFE,iZ_B0(1):iZ_E0(1))
    REAL(DP) :: W3_RMS(1:nDOFE,iZ_B0(1):iZ_E0(1))
    REAL(DP) :: W5_RMS(1:nDOFE,iZ_B0(1):iZ_E0(1))

    IF( UnitsActive )THEN

      hc3 = ( PlanckConstant * SpeedOfLight )**3

    ELSE

      hc3 = One

    END IF

    ! --- Integration Weights ---

    ASSOCIATE( dZ1 => MeshE % Width )

    E_0 = NodeCoordinate( MeshE, iZ_B0(1), 1 )

    DO iZ1    = iZ_B0(1), iZ_E0(1)
    DO iNodeE = 1, nDOFE

      E = NodeCoordinate( MeshE, iZ1, iNodeE )

      W2(iNodeE,iZ1) = FourPi * WeightsE(iNodeE) * ( dZ1(iZ1) * E**2 / hc3 )

      W3(iNodeE,iZ1) = W2(iNodeE,iZ1) * E

      W3_RMS(iNodeE,iZ1) = W2(iNodeE,iZ1) * ( E / E_0 )

      W5_RMS(iNodeE,iZ1) = W2(iNodeE,iZ1) * ( E / E_0 )**3

    END DO
    END DO

    END ASSOCIATE ! dZ1

    ! --- Initialize Gray Radiation Fields ---

    DO iS  = 1, nSpecies
    DO iGR = 1, nGR
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX

        GR(iNodeX,iZ2,iZ3,iZ4,iGR,iS) = Zero

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO

    ! --- Integrate Over Energy ---

    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX

        RMS_Int3 = Zero
        RMS_Int5 = Zero

        DO iZ1    = iZ_B0(1), iZ_E0(1)
        DO iNodeE = 1, nDOFE

          iNodeZ = (iNodeX-1) * nDOFE + iNodeE

          GR(iNodeX,iZ2,iZ3,iZ4,iGR_N,iS) &
            = GR(iNodeX,iZ2,iZ3,iZ4,iGR_N,iS) &
                + W2(iNodeE,iZ1) * CR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_N,iS)

          GR(iNodeX,iZ2,iZ3,iZ4,iGR_D,iS) &
            = GR(iNodeX,iZ2,iZ3,iZ4,iGR_D,iS) &
                + W2(iNodeE,iZ1) * PR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D,iS)

          GR(iNodeX,iZ2,iZ3,iZ4,iGR_I1,iS) &
            = GR(iNodeX,iZ2,iZ3,iZ4,iGR_I1,iS) &
                + W2(iNodeE,iZ1) * PR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS)

          GR(iNodeX,iZ2,iZ3,iZ4,iGR_I2,iS) &
            = GR(iNodeX,iZ2,iZ3,iZ4,iGR_I2,iS) &
                + W2(iNodeE,iZ1) * PR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS)

          GR(iNodeX,iZ2,iZ3,iZ4,iGR_I3,iS) &
            = GR(iNodeX,iZ2,iZ3,iZ4,iGR_I3,iS) &
                + W2(iNodeE,iZ1) * PR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS)

          GR(iNodeX,iZ2,iZ3,iZ4,iGR_J,iS) &
            = GR(iNodeX,iZ2,iZ3,iZ4,iGR_J,iS) &
                + W3(iNodeE,iZ1) * PR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D,iS)

          GR(iNodeX,iZ2,iZ3,iZ4,iGR_H1,iS) &
            = GR(iNodeX,iZ2,iZ3,iZ4,iGR_H1,iS) &
                + W3(iNodeE,iZ1) * PR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS)

          GR(iNodeX,iZ2,iZ3,iZ4,iGR_H2,iS) &
            = GR(iNodeX,iZ2,iZ3,iZ4,iGR_H2,iS) &
                + W3(iNodeE,iZ1) * PR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS)

          GR(iNodeX,iZ2,iZ3,iZ4,iGR_H3,iS) &
            = GR(iNodeX,iZ2,iZ3,iZ4,iGR_H3,iS) &
                + W3(iNodeE,iZ1) * PR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS)

          RMS_Int3 &
            = RMS_Int3 &
                + W3_RMS(iNodeE,iZ1) * PR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D,iS)

          RMS_Int5 &
            = RMS_Int5 &
                + W5_RMS(iNodeE,iZ1) * PR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D,iS)

          GR(iNodeX,iZ2,iZ3,iZ4,iGR_F,iS) &
            = GR(iNodeX,iZ2,iZ3,iZ4,iGR_F,iS) &
                + W2(iNodeE,iZ1) * PR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D,iS) &
                    * AR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iAR_F,iS)

          GR(iNodeX,iZ2,iZ3,iZ4,iGR_K,iS) &
            = GR(iNodeX,iZ2,iZ3,iZ4,iGR_K,iS) &
                + W2(iNodeE,iZ1) * PR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D,iS) &
                    * AR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iAR_K,iS)

          GR(iNodeX,iZ2,iZ3,iZ4,iGR_Q,iS) &
            = GR(iNodeX,iZ2,iZ3,iZ4,iGR_Q,iS) &
                + W2(iNodeE,iZ1) * PR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D,iS) &
                    * AR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iAR_Q,iS)

        END DO
        END DO

        GR(iNodeX,iZ2,iZ3,iZ4,iGR_RMS,iS) &
          = E_0 * SQRT( RMS_Int5 / RMS_Int3 )

        GR(iNodeX,iZ2,iZ3,iZ4,iGR_F,iS) &
          = GR(iNodeX,iZ2,iZ3,iZ4,iGR_F,iS) &
              / GR(iNodeX,iZ2,iZ3,iZ4,iGR_D,iS)

        GR(iNodeX,iZ2,iZ3,iZ4,iGR_K,iS) &
          = GR(iNodeX,iZ2,iZ3,iZ4,iGR_K,iS) &
              / GR(iNodeX,iZ2,iZ3,iZ4,iGR_D,iS)

        GR(iNodeX,iZ2,iZ3,iZ4,iGR_Q,iS) &
          = GR(iNodeX,iZ2,iZ3,iZ4,iGR_Q,iS) &
              / GR(iNodeX,iZ2,iZ3,iZ4,iGR_D,iS)

      END DO

    END DO
    END DO
    END DO
    END DO
    
  END SUBROUTINE ComputeGray_TwoMoment


  SUBROUTINE ComputeTimeStep_TwoMoment &
    ( iX_B0, iX_E0, iX_B1, iX_E1, GX, CFL, TimeStep )

    INTEGER,  INTENT(in)  :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)  :: &
      GX(1:nDOFX, &
         iX_B1(1):iX_E1(1), &
         iX_B1(2):iX_E1(2), &
         iX_B1(3):iX_E1(3), &
         1:nGF)
    REAL(DP), INTENT(in)  :: &
      CFL
    REAL(DP), INTENT(out) :: &
      TimeStep

    INTEGER  :: iX1, iX2, iX3
    REAL(DP) :: dt(3)

    TimeStep = HUGE( One )
    dt       = HUGE( One )

    ASSOCIATE &
      ( dX1 => MeshX(1) % Width, &
        dX2 => MeshX(2) % Width, &
        dX3 => MeshX(3) % Width )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      dt(1) = dX1(iX1) * MINVAL( GX(:,iX1,iX2,iX3,iGF_h_1) )

      IF( iX_E0(2) .GT. iX_B0(2) )THEN

        dt(2) = dX2(iX2) * MINVAL( GX(:,iX1,iX2,iX3,iGF_h_2) )

      END IF

      IF( iX_E0(3) .GT. iX_B0(3) )THEN

        dt(3) = dX3(iX3) * MINVAL( GX(:,iX1,iX2,iX3,iGF_h_3) )

      END IF

      TimeStep = MIN( TimeStep, MINVAL( dt ) )

    END DO
    END DO
    END DO

    END ASSOCIATE ! dX1, etc.

    TimeStep = MAX( CFL * TimeStep, SqrtTiny )

  END SUBROUTINE ComputeTimeStep_TwoMoment


  SUBROUTINE ComputeTimeStep_TwoMoment_Realizability &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, U_F, CFL, TimeStep, Verbose_Option )

    INTEGER,  INTENT(in)  :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)  :: &
      GX (1:nDOFX, &
          iZ_B1(2):iZ_E1(2), &
          iZ_B1(3):iZ_E1(3), &
          iZ_B1(4):iZ_E1(4), &
          1:nGF)
    REAL(DP), INTENT(in)  :: &
      U_F(1:nDOFX, &
          iZ_B1(2):iZ_E1(2), &
          iZ_B1(3):iZ_E1(3), &
          iZ_B1(4):iZ_E1(4), &
          1:nCF)
    REAL(DP), INTENT(in)  :: &
      CFL
    REAL(DP), INTENT(out) :: &
      TimeStep
    LOGICAL, INTENT(in), OPTIONAL :: &
      Verbose_Option

    REAL(DP), PARAMETER :: Gamma_X = Half
    REAL(DP), PARAMETER :: Gamma_E = Half

    LOGICAL  :: Verbose
    INTEGER  :: iX_B0(3), iX_E0(3)
    INTEGER  :: iX_B1(3), iX_E1(3)
    INTEGER  :: iX1, iX2, iX3, iNodeX, iE
    REAL(DP) :: xQ(nNodes), wQ(nNodes)
    REAL(DP) :: TimeStep_X, dt_X(3)
    REAL(DP) :: TimeStep_E, dt_E
    REAL(DP) :: h_d_1, h_d_2, h_d_3
    REAL(DP) :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP) :: V_d_1, V_d_2, V_d_3, AbsV
    REAL(DP) :: CFL_Eff_X, CFL_Eff_E
    REAL(DP) :: dE_Min, A(3,3), Lambda(3), Alpha_E

    REAL(DP), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: dV_u_dX1, dV_d_dX1, dGm_dd_dX1
    REAL(DP), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: dV_u_dX2, dV_d_dX2, dGm_dd_dX2
    REAL(DP), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: dV_u_dX3, dV_d_dX3, dGm_dd_dX3

    ALLOCATE(   dV_u_dX1(nDOFX,3,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)) )
    ALLOCATE(   dV_d_dX1(nDOFX,3,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)) )
    ALLOCATE( dGm_dd_dX1(nDOFX,3,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)) )
    ALLOCATE(   dV_u_dX2(nDOFX,3,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)) )
    ALLOCATE(   dV_d_dX2(nDOFX,3,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)) )
    ALLOCATE( dGm_dd_dX2(nDOFX,3,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)) )
    ALLOCATE(   dV_u_dX3(nDOFX,3,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)) )
    ALLOCATE(   dV_d_dX3(nDOFX,3,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)) )
    ALLOCATE( dGm_dd_dX3(nDOFX,3,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)) )

    IF( PRESENT( Verbose_Option ) )THEN
      Verbose = Verbose_Option
    ELSE
      Verbose = .FALSE.
    END IF

    ASSOCIATE &
      ( dX1 => MeshX(1) % Width, &
        dX2 => MeshX(2) % Width, &
        dX3 => MeshX(3) % Width, &
        dE  => MeshE    % Width, &
        E_C => MeshE    % Center )

    iX_B0 = iZ_B0(2:4); iX_E0 = iZ_E0(2:4)
    iX_B1 = iZ_B1(2:4); iX_E1 = iZ_E1(2:4)

    TimeStep_X = HUGE( One )
    TimeStep_E = HUGE( One )
    dt_X       = HUGE( One )
    dt_E       = HUGE( One )

    CALL GetQuadrature( nNodes, xQ, wQ, 'Lobatto' )

    CFL_Eff_X = Gamma_X * wQ(nNodes) / DBLE( nDimsX )
    CFL_Eff_E = Gamma_E * wQ(nNodes)

    dE_Min = HUGE( One ) ! --- Min of dE / E_H
    DO iE = iZ_B0(1), iZ_E0(1)
      dE_Min = MIN( dE_Min, dE(iE) / ( E_C(iE) + Half * dE(iE) ) )
    END DO

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: E_C, dE, dX1, dX2, dX3, iZ_B0, iZ_E0, iZ_B1, iZ_E1, &
    !$OMP          iX_B0, iX_E0, iX_B1, iX_E1, GX, U_F, &
    !$OMP          TimeStep_X, TimeStep_E, dt_X, dt_E,  &
    !$OMP          CFL, CFL_Eff_X, CFL_Eff_E, dE_Min ) &
    !$OMP MAP( alloc: dV_u_dX1, dV_d_dX1, dGm_dd_dX1, &
    !$OMP             dV_u_dX2, dV_d_dX2, dGm_dd_dX2, &
    !$OMP             dV_u_dX3, dV_d_dX3, dGm_dd_dX3 )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA &
    !$ACC COPYIN( E_C, dE, dX1, dX2, dX3, iZ_B0, iZ_E0, iZ_B1, iZ_E1, &
    !$ACC         iX_B0, iX_E0, iX_B1, iX_E1, GX, U_F, &
    !$ACC         TimeStep_X, TimeStep_E, dt_X, dt_E,  &
    !$ACC         CFL, CFL_Eff_X, CFL_Eff_E, dE_Min ) &
    !$ACC CREATE( dV_u_dX1, dV_d_dX1, dGm_dd_dX1, &
    !$ACC         dV_u_dX2, dV_d_dX2, dGm_dd_dX2, &
    !$ACC         dV_u_dX3, dV_d_dX3, dGm_dd_dX3 )
#endif

    CALL ComputeWeakDerivatives_X1 &
           ( iX_B0, iX_E0, iX_B1, iX_E1, GX, U_F, &
             dV_u_dX1, dV_d_dX1, dGm_dd_dX1 )

    CALL ComputeWeakDerivatives_X2 &
           ( iX_B0, iX_E0, iX_B1, iX_E1, GX, U_F, &
             dV_u_dX2, dV_d_dX2, dGm_dd_dX2 )

    CALL ComputeWeakDerivatives_X3 &
           ( iX_B0, iX_E0, iX_B1, iX_E1, GX, U_F, &
             dV_u_dX3, dV_d_dX3, dGm_dd_dX3 )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4) &
    !$OMP PRIVATE( A, Lambda, dt_X, AbsV, Alpha_E, dt_E,  &
    !$OMP          h_d_1, h_d_2, h_d_3, &
    !$OMP          Gm_dd_11, Gm_dd_22, Gm_dd_33, &
    !$OMP          V_d_1, V_d_2, V_d_3 ) &
    !$OMP REDUCTION( MIN : TimeStep_X, TimeStep_E )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRIVATE( A, Lambda, dt_X, AbsV, Alpha_E, dt_E,  &
    !$ACC          h_d_1, h_d_2, h_d_3, &
    !$ACC          Gm_dd_11, Gm_dd_22, Gm_dd_33, &
    !$ACC          V_d_1, V_d_2, V_d_3 ) &
    !$ACC REDUCTION( MIN : TimeStep_X, TimeStep_E )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(4) &
    !$OMP PRIVATE( A, Lambda, dt_X, AbsV, Alpha_E, dt_E,  &
    !$OMP          h_d_1, h_d_2, h_d_3, &
    !$OMP          Gm_dd_11, Gm_dd_22, Gm_dd_33, &
    !$OMP          V_d_1, V_d_2, V_d_3 ) &
    !$OMP REDUCTION( MIN : TimeStep_X, TimeStep_E )
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        h_d_1 = GX(iNodeX,iX1,iX2,iX3,iGF_h_1)
        h_d_2 = GX(iNodeX,iX1,iX2,iX3,iGF_h_2)
        h_d_3 = GX(iNodeX,iX1,iX2,iX3,iGF_h_3)

        Gm_dd_11 = MAX( h_d_1**2, SqrtTiny )
        Gm_dd_22 = MAX( h_d_2**2, SqrtTiny )
        Gm_dd_33 = MAX( h_d_3**2, SqrtTiny )

        V_d_1 = U_F(iNodeX,iX1,iX2,iX3,iCF_S1) / U_F(iNodeX,iX1,iX2,iX3,iCF_D)
        V_d_2 = U_F(iNodeX,iX1,iX2,iX3,iCF_S2) / U_F(iNodeX,iX1,iX2,iX3,iCF_D)
        V_d_3 = U_F(iNodeX,iX1,iX2,iX3,iCF_S3) / U_F(iNodeX,iX1,iX2,iX3,iCF_D)

        AbsV = SQRT(   V_d_1**2 / Gm_dd_11 &
                     + V_d_2**2 / Gm_dd_22 &
                     + V_d_3**2 / Gm_dd_33 )

        ! --- Time Step from Spatial Divergence ---

        dt_X(1) = CFL_Eff_X * ( One - AbsV ) * dX1(iX1) * h_d_1
        dt_X(2) = CFL_Eff_X * ( One - AbsV ) * dX2(iX2) * h_d_2
        dt_X(3) = CFL_Eff_X * ( One - AbsV ) * dX3(iX3) * h_d_3

        TimeStep_X = MIN( TimeStep_X, MINVAL( dt_X ) )

        ! --- Eigenvalues of Quadratic Form Matrix ---

        A(:,1) = Half * [ Two * dV_u_dX1(iNodeX,1,iX1,iX2,iX3), &
                                dV_u_dX2(iNodeX,1,iX1,iX2,iX3)  &
                              + dV_u_dX1(iNodeX,2,iX1,iX2,iX3), &
                                dV_u_dX3(iNodeX,1,iX1,iX2,iX3)  &
                              + dV_u_dX1(iNodeX,3,iX1,iX2,iX3) ]
        A(:,2) = Half * [       dV_u_dX1(iNodeX,2,iX1,iX2,iX3)  &
                              + dV_u_dX2(iNodeX,1,iX1,iX2,iX3), &
                          Two * dV_u_dX2(iNodeX,2,iX1,iX2,iX3), &
                                dV_u_dX3(iNodeX,2,iX1,iX2,iX3)  &
                              + dV_u_dX2(iNodeX,3,iX1,iX2,iX3) ]
        A(:,3) = Half * [       dV_u_dX1(iNodeX,3,iX1,iX2,iX3)  &
                              + dV_u_dX3(iNodeX,1,iX1,iX2,iX3), &
                                dV_u_dX2(iNodeX,3,iX1,iX2,iX3)  &
                              + dV_u_dX3(iNodeX,2,iX1,iX2,iX3), &
                          Two * dV_u_dX3(iNodeX,3,iX1,iX2,iX3) ]

        CALL EigenvaluesSymmetric3( A, Lambda )

        Alpha_E = MAX( MAXVAL( ABS( Lambda ) ), SqrtTiny )

        ! --- Time Step from Energy Divergence ---

        dt_E = CFL_Eff_E * ( One - AbsV ) * dE_Min / Alpha_E

        TimeStep_E = MIN( TimeStep_E, dt_E )

      END DO

    END DO
    END DO
    END DO

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from: TimeStep_X, TimeStep_E ) &
    !$OMP MAP( release: E_C, dE, dX1, dX2, dX3, iZ_B0, iZ_E0, iZ_B1, iZ_E1, &
    !$OMP               iX_B0, iX_E0, iX_B1, iX_E1, GX, U_F, &
    !$OMP               dt_X, dt_E,  &
    !$OMP               CFL, CFL_Eff_X, CFL_Eff_E, dE_Min, &
    !$OMP               dV_u_dX1, dV_d_dX1, dGm_dd_dX1, &
    !$OMP               dV_u_dX2, dV_d_dX2, dGm_dd_dX2, &
    !$OMP               dV_u_dX3, dV_d_dX3, dGm_dd_dX3 )
#elif defined( THORNADO_OACC   )
    !$ACC EXIT DATA &
    !$ACC COPYOUT( TimeStep_X, TimeStep_E ) &
    !$ACC DELETE( E_C, dE, dX1, dX2, dX3, iZ_B0, iZ_E0, iZ_B1, iZ_E1, &
    !$ACC         iX_B0, iX_E0, iX_B1, iX_E1, GX, U_F, &
    !$ACC         dt_X, dt_E,  &
    !$ACC         CFL, CFL_Eff_X, CFL_Eff_E, dE_Min, &
    !$ACC         dV_u_dX1, dV_d_dX1, dGm_dd_dX1, &
    !$ACC         dV_u_dX2, dV_d_dX2, dGm_dd_dX2, &
    !$ACC         dV_u_dX3, dV_d_dX3, dGm_dd_dX3 )
#endif

    TimeStep = MAX( CFL * MIN( TimeStep_X, TimeStep_E ), SqrtTiny )

    IF( Verbose )THEN
      WRITE(*,'(A8,A7,ES12.6E2,A8,ES12.6E2)') &
        '', 'dt_X = ', TimeStep_X, ' dt_E = ', TimeStep_E
    END IF

    END ASSOCIATE ! dX1, etc.

  END SUBROUTINE ComputeTimeStep_TwoMoment_Realizability


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

    REAL(DP)             :: Flux_E(4)
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

    IF ( FF <= SqrtTiny ) THEN
      EF = Third
      a = Third
      b = Zero
    ELSE
      EF = EddingtonFactor( D, FF )
      a = Half * ( One - EF )
      b = Half * ( Three * EF - One )
    END IF

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

    IF ( FF <= SqrtTiny ) THEN
      EF = Third
      a = Third
      b = Zero
    ELSE
      EF = EddingtonFactor( D, FF )
      a = Half * ( One - EF )
      b = Half * ( Three * EF - One )
    END IF

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

    IF ( FF <= SqrtTiny ) THEN
      EF = Third
      a = Third
      b = Zero
    ELSE
      EF = EddingtonFactor( D, FF )
      a = Half * ( One - EF )
      b = Half * ( Three * EF - One )
    END IF

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

    IF ( FF <= SqrtTiny ) THEN
      EF = Third
      a = Third
      b = Zero
    ELSE
      EF = EddingtonFactor( D, FF )
      a = Half * ( One - EF )
      b = Half * ( Three * EF - One )
    END IF

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

    IF ( FF <= SqrtTiny ) THEN
      EF = Third
      a = Third
      b = Zero
    ELSE
      EF = EddingtonFactor( D, FF )
      a = Half * ( One - EF )
      b = Half * ( Three * EF - One )
    END IF

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

    IF ( FF <= SqrtTiny ) THEN
      EF = Third
      a = Third
      b = Zero
    ELSE
      EF = EddingtonFactor( D, FF )
      a = Half * ( One - EF )
      b = Half * ( Three * EF - One )
    END IF

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

    IF ( FF <= SqrtTiny ) THEN
      EF = Third
      a = Third
      b = Zero
    ELSE
      EF = EddingtonFactor( D, FF )
      a = Half * ( One - EF )
      b = Half * ( Three * EF - One )
    END IF

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

    IF ( FF <= SqrtTiny ) THEN
      HF = Zero
      a = Zero
      b = Zero
    ELSE
      HF = HeatFluxFactor( D, FF )
      a = Half * ( FF - HF )
      b = Half * ( Five * HF - Three * FF )
    END IF

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

    IF ( FF <= SqrtTiny ) THEN
      HF = Zero
      a = Zero
      b = Zero
    ELSE
      HF = HeatFluxFactor( D, FF )
      a = Half * ( FF - HF )
      b = Half * ( Five * HF - Three * FF )
    END IF

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

    IF ( FF <= SqrtTiny ) THEN
      HF = Zero
      a = Zero
      b = Zero
    ELSE
      HF = HeatFluxFactor( D, FF )
      a = Half * ( FF - HF )
      b = Half * ( Five * HF - Three * FF )
    END IF

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


  SUBROUTINE ComputeWeakDerivatives_X1 &
    ( iX_B0, iX_E0, iX_B1, iX_E1, GX, U_F, dV_u_dX1_Out, dV_d_dX1_Out, &
      dGm_dd_dX1_Out )

    INTEGER,  INTENT(in)  :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)  :: &
      GX (1:nDOFX, &
          iX_B1(1):iX_E1(1), &
          iX_B1(2):iX_E1(2), &
          iX_B1(3):iX_E1(3), &
          1:nGF)
    REAL(DP), INTENT(in)  :: &
      U_F(1:nDOFX, &
          iX_B1(1):iX_E1(1), &
          iX_B1(2):iX_E1(2), &
          iX_B1(3):iX_E1(3), &
          1:nCF)
    REAL(DP), INTENT(out) :: &
      dV_u_dX1_Out &
        (1:nDOFX,1:3, &
         iX_B0(1):iX_E0(1), &
         iX_B0(2):iX_E0(2), &
         iX_B0(3):iX_E0(3)), &
      dV_d_dX1_Out &
        (1:nDOFX,1:3, &
         iX_B0(1):iX_E0(1), &
         iX_B0(2):iX_E0(2), &
         iX_B0(3):iX_E0(3)), &
      dGm_dd_dX1_Out &
        (1:nDOFX,1:3, &
         iX_B0(1):iX_E0(1), &
         iX_B0(2):iX_E0(2), &
         iX_B0(3):iX_E0(3))

    INTEGER  :: iNodeX
    INTEGER  :: iX1, iX2, iX3, i
    INTEGER  :: iCF, iCF_S
    INTEGER  :: iGF, iGF_h, iGF_Gm_dd
    INTEGER  :: nX(3), nX_X1(3), nK_X, nX1_X
    REAL(DP) :: uV_L(3), uV_R(3), uV_F(3), uV_K

    ! --- Geometry Fields ---

    REAL(DP), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: GX_K, GX_F, h_d_F, h_d_K, dh_d_dX1

    ! --- Conserved Fluid Fields ---

    REAL(DP), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: U_F_K, U_F_L, U_F_R

    ! --- Velocities ---

    REAL(DP), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: V_u_X1, V_d_X1, V_u_K, V_d_K, dV_u_dX1, dV_d_dX1

    ALLOCATE( GX_K    (nDOFX   ,nGF,iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),iX_B0(1)-1:iX_E0(1)+1) )
    ALLOCATE( GX_F    (nDOFX_X1,nGF,iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),iX_B0(1)  :iX_E0(1)+1) )
    ALLOCATE( h_d_F   (nDOFX_X1,3  ,iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),iX_B0(1)  :iX_E0(1)+1) )
    ALLOCATE( h_d_K   (nDOFX   ,3  ,iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),iX_B0(1)  :iX_E0(1)  ) )
    ALLOCATE( dh_d_dX1(nDOFX   ,3  ,iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),iX_B0(1)  :iX_E0(1)  ) )
    ALLOCATE( U_F_K   (nDOFX   ,nCF,iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),iX_B0(1)-1:iX_E0(1)+1) )
    ALLOCATE( U_F_L   (nDOFX_X1,nCF,iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),iX_B0(1)  :iX_E0(1)+1) )
    ALLOCATE( U_F_R   (nDOFX_X1,nCF,iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),iX_B0(1)  :iX_E0(1)+1) )
    ALLOCATE( V_u_X1  (nDOFX_X1,3  ,iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),iX_B0(1)  :iX_E0(1)+1) )
    ALLOCATE( V_d_X1  (nDOFX_X1,3  ,iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),iX_B0(1)  :iX_E0(1)+1) )
    ALLOCATE( V_u_K   (nDOFX   ,3  ,iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),iX_B0(1)  :iX_E0(1)  ) )
    ALLOCATE( V_d_K   (nDOFX   ,3  ,iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),iX_B0(1)  :iX_E0(1)  ) )
    ALLOCATE( dV_u_dX1(nDOFX   ,3  ,iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),iX_B0(1)  :iX_E0(1)  ) )
    ALLOCATE( dV_d_dX1(nDOFX   ,3  ,iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),iX_B0(1)  :iX_E0(1)  ) )

    IF( iX_E0(1) .EQ. iX_B0(1) )THEN

#if   defined( THORNADO_OMP_OL )
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5) &
      !$OMP MAP( to: iX_B0, iX_E0 )
#elif defined( THORNADO_OACC   )
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
      !$ACC COPYIN( iX_B0, iX_E0 ) &
      !$ACC PRESENT( dV_u_dX1_Out, dV_d_dX1_Out, dGm_dd_dX1_Out )
#elif defined( THORNADO_OMP    )
      !$OMP PARALLEL DO COLLAPSE(5)
#endif
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO i = 1, 3
      DO iNodeX = 1, nDOFX

          dV_u_dX1_Out  (iNodeX,i,iX1,iX2,iX3) = Zero
          dV_d_dX1_Out  (iNodeX,i,iX1,iX2,iX3) = Zero
          dGm_dd_dX1_Out(iNodeX,i,iX1,iX2,iX3) = Zero

      END DO
      END DO
      END DO
      END DO
      END DO

      RETURN
    END IF

    nX    = iX_E0 - iX_B0 + 1 ! --- Number of Elements per Spatial Dimension
    nX_X1 = nX + [ 1, 0, 0 ]  ! --- Number of X1 Faces per Spatial Dimension
    nK_X  = PRODUCT( nX )     ! --- Number of Elements in Position Space
    nX1_X = PRODUCT( nX_X1 )  ! --- Number of X1 Faces in Position Space

    ASSOCIATE( dX1 => MeshX(1) % Width )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: dX1, iX_B0, iX_E0 ) &
    !$OMP MAP( alloc: GX_K, GX_F, h_d_F, h_d_K, dh_d_dX1, &
    !$OMP             U_F_K, U_F_L, U_F_R, V_u_X1, V_d_X1, V_u_k, V_d_k, &
    !$OMP             dV_u_dX1, dV_d_dX1 )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA &
    !$ACC COPYIN( dX1, iX_B0, iX_E0 ) &
    !$ACC CREATE( GX_K, GX_F, h_d_F, h_d_K, dh_d_dX1, &
    !$ACC         U_F_K, U_F_L, U_F_R, V_u_X1, V_d_X1, V_u_k, V_d_k, &
    !$ACC         dV_u_dX1, dV_d_dX1 )
#endif

    ! --- Permute Geometry Fields ---

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( GX_K, GX, iX_B0, iX_E0 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iX1 = iX_B0(1)-1, iX_E0(1)+1
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)

      DO iGF = 1, nGF
      DO iNodeX = 1, nDOFX

        GX_K(iNodeX,iGF,iX2,iX3,iX1) = GX(iNodeX,iX1,iX2,iX3,iGF)

      END DO
      END DO

    END DO
    END DO
    END DO

    !---------------------
    ! --- Surface Term ---
    !---------------------

    CALL TimersStart( Timer_Streaming_LinearAlgebra )

    ! --- Interpolate Geometry Fields on Shared Face ---

    ! --- Face States (Average of Left and Right States) ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nX1_X*nGF, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
             GX_K(1,1,iX_B0(2),iX_B0(3),iX_B0(1)-1), nDOFX, Zero, &
             GX_F(1,1,iX_B0(2),iX_B0(3),iX_B0(1)  ), nDOFX_X1 )

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nX1_X*nGF, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
             GX_K(1,1,iX_B0(2),iX_B0(3),iX_B0(1)  ), nDOFX, Half, &
             GX_F(1,1,iX_B0(2),iX_B0(3),iX_B0(1)  ), nDOFX_X1 )

    CALL TimersStop( Timer_Streaming_LinearAlgebra )

    ! --- Compute Metric Components from Scale Factors ---

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( GX_F, h_d_F, iX_B0, iX_E0 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(4)
#endif
    DO iX1  = iX_B0(1), iX_E0(1)+1
    DO iX3  = iX_B0(3), iX_E0(3)
    DO iX2  = iX_B0(2), iX_E0(2)

      DO iNodeX = 1, nDOFX_X1

        GX_F(iNodeX,iGF_Gm_dd_11,iX2,iX3,iX1) &
          = MAX( GX_F(iNodeX,iGF_h_1,iX2,iX3,iX1)**2, SqrtTiny )
        GX_F(iNodeX,iGF_Gm_dd_22,iX2,iX3,iX1) &
          = MAX( GX_F(iNodeX,iGF_h_2,iX2,iX3,iX1)**2, SqrtTiny )
        GX_F(iNodeX,iGF_Gm_dd_33,iX2,iX3,iX1) &
          = MAX( GX_F(iNodeX,iGF_h_3,iX2,iX3,iX1)**2, SqrtTiny )
        GX_F(iNodeX,iGF_SqrtGm,iX2,iX3,iX1) &
          = SQRT(   GX_F(iNodeX,iGF_Gm_dd_11,iX2,iX3,iX1) &
                  * GX_F(iNodeX,iGF_Gm_dd_22,iX2,iX3,iX1) &
                  * GX_F(iNodeX,iGF_Gm_dd_33,iX2,iX3,iX1) )

        h_d_F(iNodeX,1,iX2,iX3,iX1) &
          = GX_F(iNodeX,iGF_h_1,iX2,iX3,iX1) * WeightsX_X1(iNodeX)
        h_d_F(iNodeX,2,iX2,iX3,iX1) &
          = GX_F(iNodeX,iGF_h_2,iX2,iX3,iX1) * WeightsX_X1(iNodeX)
        h_d_F(iNodeX,3,iX2,iX3,iX1) &
          = GX_F(iNodeX,iGF_h_3,iX2,iX3,iX1) * WeightsX_X1(iNodeX)

      END DO

    END DO
    END DO
    END DO

    CALL TimersStart( Timer_Streaming_LinearAlgebra )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 3*nK_X, nDOFX_X1, - One, LX_X1_Dn, nDOFX_X1, &
             h_d_F(1,1,iX_B0(2),iX_B0(3),iX_B0(1)  ), nDOFX_X1, Zero, &
             dh_d_dX1, nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 3*nK_X, nDOFX_X1, + One, LX_X1_Up, nDOFX_X1, &
             h_d_F(1,1,iX_B0(2),iX_B0(3),iX_B0(1)+1), nDOFX_X1, One,  &
             dh_d_dX1, nDOFX )

    CALL TimersStop( Timer_Streaming_LinearAlgebra )

    ! --- Permute Fluid Fields ---

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( U_F_K, U_F, iX_B0, iX_E0 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iX1 = iX_B0(1)-1, iX_E0(1)+1
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)

      DO iCF = 1, nCF
      DO iNodeX = 1, nDOFX

        U_F_K(iNodeX,iCF,iX2,iX3,iX1) = U_F(iNodeX,iX1,iX2,iX3,iCF)

      END DO
      END DO

    END DO
    END DO
    END DO

    ! --- Interpolate Fluid Fields ---

    CALL TimersStart( Timer_Streaming_LinearAlgebra )

      ! --- Interpolate Left State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nX1_X*nCF, nDOFX, One, LX_X1_Up, nDOFX_X1, &
             U_F_K(1,1,iX_B0(2),iX_B0(3),iX_B0(1)-1), nDOFX, Zero, &
             U_F_L(1,1,iX_B0(2),iX_B0(3),iX_B0(1)  ), nDOFX_X1 )

    ! --- Interpolate Right State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nX1_X*nCF, nDOFX, One, LX_X1_Dn, nDOFX_X1, &
             U_F_K(1,1,iX_B0(2),iX_B0(3),iX_B0(1)  ), nDOFX, Zero, &
             U_F_R(1,1,iX_B0(2),iX_B0(3),iX_B0(1)  ), nDOFX_X1 )

    CALL TimersStop( Timer_Streaming_LinearAlgebra )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4) &
    !$OMP PRIVATE( uV_L, uV_R, uV_F, iCF_S, iGF_Gm_dd )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRIVATE( uV_L, uV_R, uV_F, iCF_S, iGF_Gm_dd ) &
    !$ACC PRESENT( iX_B0, iX_E0, U_F_L, U_F_R, GX_F, &
    !$ACC          V_u_X1, V_d_X1, WeightsX_X1 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(4) &
    !$OMP PRIVATE( uV_L, uV_R, uV_F, iCF_S, iGF_Gm_dd )
#endif
    DO iX1 = iX_B0(1), iX_E0(1)+1
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)

      DO iNodeX = 1, nDOFX_X1

        DO i = 1, 3

          iCF_S     = iCF_S1       + i - 1
          iGF_Gm_dd = iGF_Gm_dd_11 + i - 1

          ! --- Left States ---

          uV_L(i) = U_F_L(iNodeX,iCF_S,iX2,iX3,iX1) &
                    / ( GX_F (iNodeX,iGF_Gm_dd,iX2,iX3,iX1) &
                      * U_F_L(iNodeX,iCF_D ,iX2,iX3,iX1) )

          ! --- Right States ---

          uV_R(i) = U_F_R(iNodeX,iCF_S,iX2,iX3,iX1) &
                    / ( GX_F (iNodeX,iGF_Gm_dd,iX2,iX3,iX1) &
                      * U_F_R(iNodeX,iCF_D ,iX2,iX3,iX1) )

        END DO

        CALL FaceVelocity_X1 &
               ( uV_L(1), uV_L(2), uV_L(3), &
                 uV_R(1), uV_R(2), uV_R(3), &
                 uV_F(1), uV_F(2), uV_F(3) )

        DO i = 1, 3

          iGF_Gm_dd = iGF_Gm_dd_11 + i - 1

          V_u_X1(iNodeX,i,iX2,iX3,iX1) &
            = uV_F(i) * WeightsX_X1(iNodeX)

          V_d_X1(iNodeX,i,iX2,iX3,iX1) &
            = uV_F(i) * WeightsX_X1(iNodeX) * GX_F(iNodeX,iGF_Gm_dd,iX2,iX3,iX1)

        END DO

      END DO

    END DO
    END DO
    END DO

    CALL TimersStart( Timer_Streaming_LinearAlgebra )

    ! --- Surface Contributions ---

    ! --- Contribution from Left Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 3*nK_X, nDOFX_X1, - One, LX_X1_Dn, nDOFX_X1, &
             V_u_X1(1,1,iX_B0(2),iX_B0(3),iX_B0(1)  ), nDOFX_X1, Zero, &
             dV_u_dX1, nDOFX )

    ! --- Contribution from Right Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 3*nK_X, nDOFX_X1, + One, LX_X1_Up, nDOFX_X1, &
             V_u_X1(1,1,iX_B0(2),iX_B0(3),iX_B0(1)+1), nDOFX_X1, One,  &
             dV_u_dX1, nDOFX )

    ! --- Contribution from Left Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 3*nK_X, nDOFX_X1, - One, LX_X1_Dn, nDOFX_X1, &
             V_d_X1(1,1,iX_B0(2),iX_B0(3),iX_B0(1)  ), nDOFX_X1, Zero, &
             dV_d_dX1, nDOFX )

    ! --- Contribution from Right Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 3*nK_X, nDOFX_X1, + One, LX_X1_Up, nDOFX_X1, &
             V_d_X1(1,1,iX_B0(2),iX_B0(3),iX_B0(1)+1), nDOFX_X1, One,  &
             dV_d_dX1, nDOFX )

    CALL TimersStop( Timer_Streaming_LinearAlgebra )

    ! -------------------
    ! --- Volume Term ---
    ! -------------------

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5) &
    !$OMP PRIVATE( uV_K, iCF, iCF_S, iGF_Gm_dd, iGF_h )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRIVATE( uV_K, iCF, iCF_S, iGF_Gm_dd, iGF_h ) &
    !$ACC PRESENT( iX_B0, iX_E0, U_F_K, GX_K, h_d_K, &
    !$ACC          V_u_K, V_d_K, WeightsX_q )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5) &
    !$OMP PRIVATE( uV_K, iCF, iCF_S, iGF_Gm_dd, iGF_h )
#endif
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)

      DO i = 1, 3
      DO iNodeX = 1, nDOFX

        iCF       = iCF_S1       + i - 1
        iGF_Gm_dd = iGF_Gm_dd_11 + i - 1
        iGF_h     = iGF_h_1      + i - 1

        h_d_K(iNodeX,i,iX2,iX3,iX1) &
          = WeightsX_q(iNodeX) * GX_K(iNodeX,iGF_h,iX2,iX3,iX1)

        uV_K &
          = U_F_K(iNodeX,iCF,iX2,iX3,iX1) &
            / ( GX_K (iNodeX,iGF_Gm_dd,iX2,iX3,iX1) &
              * U_F_K(iNodeX,iCF_D ,iX2,iX3,iX1) )

        V_u_K(iNodeX,i,iX2,iX3,iX1) &
          = uV_K * WeightsX_q(iNodeX)

        V_d_K(iNodeX,i,iX2,iX3,iX1) &
          = uV_K * WeightsX_q(iNodeX) * GX_K(iNodeX,iGF_Gm_dd,iX2,iX3,iX1)

      END DO
      END DO

    END DO
    END DO
    END DO

    CALL TimersStart( Timer_Streaming_LinearAlgebra )

    ! --- Volume Contributions ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 3*nK_X, nDOFX, - One, dLXdX1_q, nDOFX, &
             h_d_K, nDOFX, One, dh_d_dX1, nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 3*nK_X, nDOFX, - One, dLXdX1_q, nDOFX, &
             V_u_K, nDOFX, One, dV_u_dX1, nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 3*nK_X, nDOFX, - One, dLXdX1_q, nDOFX, &
             V_d_K, nDOFX, One, dV_d_dX1, nDOFX )

    CALL TimersStop( Timer_Streaming_LinearAlgebra )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( iX_B0, iX_E0, dX1, &
    !$ACC          dh_d_dX1, dV_u_dX1, dV_d_dX1, WeightsX_q )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)

      DO i      = 1, 3
      DO iNodeX = 1, nDOFX

        dh_d_dX1(iNodeX,i,iX2,iX3,iX1) &
          = dh_d_dX1(iNodeX,i,iX2,iX3,iX1) &
              / ( WeightsX_q(iNodeX) * dX1(iX1) )

        dV_u_dX1(iNodeX,i,iX2,iX3,iX1) &
         = dV_u_dX1(iNodeX,i,iX2,iX3,iX1) &
             / ( WeightsX_q(iNodeX) * dX1(iX1) )

        dV_d_dX1(iNodeX,i,iX2,iX3,iX1) &
         = dV_d_dX1(iNodeX,i,iX2,iX3,iX1) &
             / ( WeightsX_q(iNodeX) * dX1(iX1) )

      END DO
      END DO

    END DO
    END DO
    END DO

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5) &
    !$OMP PRIVATE( iGF_h )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRIVATE( iGF_h ) &
    !$ACC PRESENT( iX_B0, iX_E0, GX_K, dGm_dd_dX1_Out, dh_d_dX1, &
    !$ACC          dV_u_dX1_Out, dV_u_dX1, dV_d_dX1_Out, dV_d_dX1 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5) &
    !$OMP PRIVATE( iGF_h )
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO i = 1, 3
      DO iNodeX = 1, nDOFX

        iGF_h = iGF_h_1 + i - 1

        dGm_dd_dX1_Out(iNodeX,i,iX1,iX2,iX3) &
          = Two * GX_K(iNodeX,iGF_h,iX2,iX3,iX1) &
              * dh_d_dX1(iNodeX,i,iX2,iX3,iX1)

        dV_u_dX1_Out(iNodeX,i,iX1,iX2,iX3) &
          = dV_u_dX1(iNodeX,i,iX2,iX3,iX1)

        dV_d_dX1_Out(iNodeX,i,iX1,iX2,iX3) &
          = dV_d_dX1(iNodeX,i,iX2,iX3,iX1)

      END DO
      END DO

    END DO
    END DO
    END DO

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: dX1, iX_B0, iX_E0, &
    !$OMP               GX_K, GX_F, h_d_F, h_d_K, dh_d_dX1, &
    !$OMP               U_F_K, U_F_L, U_F_R, V_u_X1, V_d_X1, V_u_k, V_d_k, &
    !$OMP               dV_u_dX1, dV_d_dX1 )
#elif defined( THORNADO_OACC   )
    !$ACC EXIT DATA &
    !$ACC DELETE( dX1, iX_B0, iX_E0, &
    !$ACC         GX_K, GX_F, h_d_F, h_d_K, dh_d_dX1, &
    !$ACC         U_F_K, U_F_L, U_F_R, V_u_X1, V_d_X1, V_u_k, V_d_k, &
    !$ACC         dV_u_dX1, dV_d_dX1 )
#endif

    END ASSOCIATE

  END SUBROUTINE ComputeWeakDerivatives_X1


  SUBROUTINE ComputeWeakDerivatives_X2 &
    ( iX_B0, iX_E0, iX_B1, iX_E1, GX, U_F, dV_u_dX2_Out, dV_d_dX2_Out, &
      dGm_dd_dX2_Out )

    INTEGER,  INTENT(in)  :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)  :: &
      GX (1:nDOFX, &
          iX_B1(1):iX_E1(1), &
          iX_B1(2):iX_E1(2), &
          iX_B1(3):iX_E1(3), &
          1:nGF)
    REAL(DP), INTENT(in)  :: &
      U_F(1:nDOFX, &
          iX_B1(1):iX_E1(1), &
          iX_B1(2):iX_E1(2), &
          iX_B1(3):iX_E1(3), &
          1:nCF)
    REAL(DP), INTENT(out) :: &
      dV_u_dX2_Out &
        (1:nDOFX,1:3, &
         iX_B0(1):iX_E0(1), &
         iX_B0(2):iX_E0(2), &
         iX_B0(3):iX_E0(3)), &
      dV_d_dX2_Out &
        (1:nDOFX,1:3, &
         iX_B0(1):iX_E0(1), &
         iX_B0(2):iX_E0(2), &
         iX_B0(3):iX_E0(3)), &
      dGm_dd_dX2_Out &
        (1:nDOFX,1:3, &
         iX_B0(1):iX_E0(1), &
         iX_B0(2):iX_E0(2), &
         iX_B0(3):iX_E0(3))

    INTEGER  :: iNodeX
    INTEGER  :: iX1, iX2, iX3, i
    INTEGER  :: iCF, iCF_S
    INTEGER  :: iGF, iGF_h, iGF_Gm_dd
    INTEGER  :: nX(3), nX_X2(3), nK_X, nX2_X
    REAL(DP) :: uV_L(3), uV_R(3), uV_F(3), uV_K

    ! --- Geometry Fields ---

    REAL(DP), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: GX_K, GX_F, h_d_F, h_d_K, dh_d_dX2

    ! --- Conserved Fluid Fields ---

    REAL(DP), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: U_F_K, U_F_L, U_F_R

    ! --- Velocities ---

    REAL(DP), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: V_u_X2, V_d_X2, V_u_K, V_d_K, dV_u_dX2, dV_d_dX2

    ALLOCATE( GX_K    (nDOFX   ,nGF,iX_B0(1):iX_E0(1),iX_B0(3):iX_E0(3),iX_B0(2)-1:iX_E0(2)+1) )
    ALLOCATE( GX_F    (nDOFX_X2,nGF,iX_B0(1):iX_E0(1),iX_B0(3):iX_E0(3),iX_B0(2)  :iX_E0(2)+1) )
    ALLOCATE( h_d_F   (nDOFX_X2,3  ,iX_B0(1):iX_E0(1),iX_B0(3):iX_E0(3),iX_B0(2)  :iX_E0(2)+1) )
    ALLOCATE( h_d_K   (nDOFX   ,3  ,iX_B0(1):iX_E0(1),iX_B0(3):iX_E0(3),iX_B0(2)  :iX_E0(2)  ) )
    ALLOCATE( dh_d_dX2(nDOFX   ,3  ,iX_B0(1):iX_E0(1),iX_B0(3):iX_E0(3),iX_B0(2)  :iX_E0(2)  ) )
    ALLOCATE( U_F_K   (nDOFX   ,nCF,iX_B0(1):iX_E0(1),iX_B0(3):iX_E0(3),iX_B0(2)-1:iX_E0(2)+1) )
    ALLOCATE( U_F_L   (nDOFX_X2,nCF,iX_B0(1):iX_E0(1),iX_B0(3):iX_E0(3),iX_B0(2)  :iX_E0(2)+1) )
    ALLOCATE( U_F_R   (nDOFX_X2,nCF,iX_B0(1):iX_E0(1),iX_B0(3):iX_E0(3),iX_B0(2)  :iX_E0(2)+1) )
    ALLOCATE( V_u_X2  (nDOFX_X2,3  ,iX_B0(1):iX_E0(1),iX_B0(3):iX_E0(3),iX_B0(2)  :iX_E0(2)+1) )
    ALLOCATE( V_d_X2  (nDOFX_X2,3  ,iX_B0(1):iX_E0(1),iX_B0(3):iX_E0(3),iX_B0(2)  :iX_E0(2)+1) )
    ALLOCATE( V_u_K   (nDOFX   ,3  ,iX_B0(1):iX_E0(1),iX_B0(3):iX_E0(3),iX_B0(2)  :iX_E0(2)  ) )
    ALLOCATE( V_d_K   (nDOFX   ,3  ,iX_B0(1):iX_E0(1),iX_B0(3):iX_E0(3),iX_B0(2)  :iX_E0(2)  ) )
    ALLOCATE( dV_u_dX2(nDOFX   ,3  ,iX_B0(1):iX_E0(1),iX_B0(3):iX_E0(3),iX_B0(2)  :iX_E0(2)  ) )
    ALLOCATE( dV_d_dX2(nDOFX   ,3  ,iX_B0(1):iX_E0(1),iX_B0(3):iX_E0(3),iX_B0(2)  :iX_E0(2)  ) )

    IF( iX_E0(2) .EQ. iX_B0(2) )THEN

#if   defined( THORNADO_OMP_OL )
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5) &
      !$OMP MAP( to: iX_B0, iX_E0 )
#elif defined( THORNADO_OACC   )
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
      !$ACC COPYIN( iX_B0, iX_E0 ) &
      !$ACC PRESENT( dV_u_dX2_Out, dV_d_dX2_Out, dGm_dd_dX2_Out )
#elif defined( THORNADO_OMP    )
      !$OMP PARALLEL DO COLLAPSE(5)
#endif
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO i = 1, 3
      DO iNodeX = 1, nDOFX

          dV_u_dX2_Out  (iNodeX,i,iX1,iX2,iX3) = Zero
          dV_d_dX2_Out  (iNodeX,i,iX1,iX2,iX3) = Zero
          dGm_dd_dX2_Out(iNodeX,i,iX1,iX2,iX3) = Zero

      END DO
      END DO
      END DO
      END DO
      END DO

      RETURN
    END IF

    nX    = iX_E0 - iX_B0 + 1 ! --- Number of Elements per Spatial Dimension
    nX_X2 = nX + [ 0, 1, 0 ]  ! --- Number of X2 Faces per Spatial Dimension
    nK_X  = PRODUCT( nX )     ! --- Number of Elements in Position Space
    nX2_X = PRODUCT( nX_X2 )  ! --- Number of X2 Faces in Position Space

    ASSOCIATE( dX2 => MeshX(2) % Width )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: dX2, iX_B0, iX_E0 ) &
    !$OMP MAP( alloc: GX_K, GX_F, h_d_F, h_d_K, dh_d_dX2, &
    !$OMP             U_F_K, U_F_L, U_F_R, V_u_X2, V_d_X2, V_u_k, V_d_k, &
    !$OMP             dV_u_dX2, dV_d_dX2 )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA &
    !$ACC COPYIN( dX2, iX_B0, iX_E0 ) &
    !$ACC CREATE( GX_K, GX_F, h_d_F, h_d_K, dh_d_dX2, &
    !$ACC         U_F_K, U_F_L, U_F_R, V_u_X2, V_d_X2, V_u_k, V_d_k, &
    !$ACC         dV_u_dX2, dV_d_dX2 )
#endif

    ! --- Permute Geometry Fields ---

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( GX_K, GX, iX_B0, iX_E0 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iX2 = iX_B0(2)-1, iX_E0(2)+1
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iGF = 1, nGF
      DO iNodeX = 1, nDOFX

        GX_K(iNodeX,iGF,iX1,iX3,iX2) = GX(iNodeX,iX1,iX2,iX3,iGF)

      END DO
      END DO

    END DO
    END DO
    END DO

    !---------------------
    ! --- Surface Term ---
    !---------------------

    CALL TimersStart( Timer_Streaming_LinearAlgebra )

    ! --- Interpolate Geometry Fields on Shared Face ---

    ! --- Face States (Average of Left and Right States) ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nX2_X*nGF, nDOFX, One,  LX_X2_Up, nDOFX_X2, &
             GX_K(1,1,iX_B0(1),iX_B0(3),iX_B0(2)-1), nDOFX, Zero, &
             GX_F(1,1,iX_B0(1),iX_B0(3),iX_B0(2)  ), nDOFX_X2 )

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nX2_X*nGF, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
             GX_K(1,1,iX_B0(1),iX_B0(3),iX_B0(2)  ), nDOFX, Half, &
             GX_F(1,1,iX_B0(1),iX_B0(3),iX_B0(2)  ), nDOFX_X2 )

    CALL TimersStop( Timer_Streaming_LinearAlgebra )

    ! --- Compute Metric Components from Scale Factors ---

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( GX_F, h_d_F, iX_B0, iX_E0 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(4)
#endif
    DO iX2  = iX_B0(2), iX_E0(2)+1
    DO iX3  = iX_B0(3), iX_E0(3)
    DO iX1  = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX_X2

        GX_F(iNodeX,iGF_Gm_dd_11,iX1,iX3,iX2) &
          = MAX( GX_F(iNodeX,iGF_h_1,iX1,iX3,iX2)**2, SqrtTiny )
        GX_F(iNodeX,iGF_Gm_dd_22,iX1,iX3,iX2) &
          = MAX( GX_F(iNodeX,iGF_h_2,iX1,iX3,iX2)**2, SqrtTiny )
        GX_F(iNodeX,iGF_Gm_dd_33,iX1,iX3,iX2) &
          = MAX( GX_F(iNodeX,iGF_h_3,iX1,iX3,iX2)**2, SqrtTiny )
        GX_F(iNodeX,iGF_SqrtGm,iX1,iX3,iX2) &
          = SQRT(   GX_F(iNodeX,iGF_Gm_dd_11,iX1,iX3,iX2) &
                  * GX_F(iNodeX,iGF_Gm_dd_22,iX1,iX3,iX2) &
                  * GX_F(iNodeX,iGF_Gm_dd_33,iX1,iX3,iX2) )

        h_d_F(iNodeX,1,iX1,iX3,iX2) &
          = GX_F(iNodeX,iGF_h_1,iX1,iX3,iX2) * WeightsX_X2(iNodeX)
        h_d_F(iNodeX,2,iX1,iX3,iX2) &
          = GX_F(iNodeX,iGF_h_2,iX1,iX3,iX2) * WeightsX_X2(iNodeX)
        h_d_F(iNodeX,3,iX1,iX3,iX2) &
          = GX_F(iNodeX,iGF_h_3,iX1,iX3,iX2) * WeightsX_X2(iNodeX)

      END DO

    END DO
    END DO
    END DO

    CALL TimersStart( Timer_Streaming_LinearAlgebra )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 3*nK_X, nDOFX_X2, - One, LX_X2_Dn, nDOFX_X2, &
             h_d_F(1,1,iX_B0(1),iX_B0(3),iX_B0(2)  ), nDOFX_X2, Zero, &
             dh_d_dX2, nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 3*nK_X, nDOFX_X2, + One, LX_X2_Up, nDOFX_X2, &
             h_d_F(1,1,iX_B0(1),iX_B0(3),iX_B0(2)+1), nDOFX_X2, One,  &
             dh_d_dX2, nDOFX )

    CALL TimersStop( Timer_Streaming_LinearAlgebra )

    ! --- Permute Fluid Fields ---

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( U_F_K, U_F, iX_B0, iX_E0 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iX2 = iX_B0(2)-1, iX_E0(2)+1
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iCF = 1, nCF
      DO iNodeX = 1, nDOFX

        U_F_K(iNodeX,iCF,iX1,iX3,iX2) = U_F(iNodeX,iX1,iX2,iX3,iCF)

      END DO
      END DO

    END DO
    END DO
    END DO

    ! --- Interpolate Fluid Fields ---

    CALL TimersStart( Timer_Streaming_LinearAlgebra )

      ! --- Interpolate Left State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nX2_X*nCF, nDOFX, One, LX_X2_Up, nDOFX_X2, &
             U_F_K(1,1,iX_B0(1),iX_B0(3),iX_B0(2)-1), nDOFX, Zero, &
             U_F_L(1,1,iX_B0(1),iX_B0(3),iX_B0(2)  ), nDOFX_X2 )

    ! --- Interpolate Right State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nX2_X*nCF, nDOFX, One, LX_X2_Dn, nDOFX_X2, &
             U_F_K(1,1,iX_B0(1),iX_B0(3),iX_B0(2)  ), nDOFX, Zero, &
             U_F_R(1,1,iX_B0(1),iX_B0(3),iX_B0(2)  ), nDOFX_X2 )

    CALL TimersStop( Timer_Streaming_LinearAlgebra )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4) &
    !$OMP PRIVATE( uV_L, uV_R, uV_F, iCF_S, iGF_Gm_dd )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRIVATE( uV_L, uV_R, uV_F, iCF_S, iGF_Gm_dd ) &
    !$ACC PRESENT( iX_B0, iX_E0, U_F_L, U_F_R, GX_F, &
    !$ACC          V_u_X2, V_d_X2, WeightsX_X2 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(4) &
    !$OMP PRIVATE( uV_L, uV_R, uV_F, iCF_S, iGF_Gm_dd )
#endif
    DO iX2 = iX_B0(2), iX_E0(2)+1
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX_X2

        DO i = 1, 3

          iCF_S     = iCF_S1       + i - 1
          iGF_Gm_dd = iGF_Gm_dd_11 + i - 1

          ! --- Left States ---

          uV_L(i) = U_F_L(iNodeX,iCF_S,iX1,iX3,iX2) &
                    / ( GX_F (iNodeX,iGF_Gm_dd,iX1,iX3,iX2) &
                      * U_F_L(iNodeX,iCF_D ,iX1,iX3,iX2) )

          ! --- Right States ---

          uV_R(i) = U_F_R(iNodeX,iCF_S,iX1,iX3,iX2) &
                    / ( GX_F (iNodeX,iGF_Gm_dd,iX1,iX3,iX2) &
                      * U_F_R(iNodeX,iCF_D ,iX1,iX3,iX2) )

        END DO

        CALL FaceVelocity_X2 &
               ( uV_L(1), uV_L(2), uV_L(3), &
                 uV_R(1), uV_R(2), uV_R(3), &
                 uV_F(1), uV_F(2), uV_F(3) )

        DO i = 1, 3

          iGF_Gm_dd = iGF_Gm_dd_11 + i - 1

          V_u_X2(iNodeX,i,iX1,iX3,iX2) &
            = uV_F(i) * WeightsX_X2(iNodeX)

          V_d_X2(iNodeX,i,iX1,iX3,iX2) &
            = uV_F(i) * WeightsX_X2(iNodeX) * GX_F(iNodeX,iGF_Gm_dd,iX1,iX3,iX2)

        END DO

      END DO

    END DO
    END DO
    END DO

    CALL TimersStart( Timer_Streaming_LinearAlgebra )

    ! --- Surface Contributions ---

    ! --- Contribution from Left Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 3*nK_X, nDOFX_X2, - One, LX_X2_Dn, nDOFX_X2, &
             V_u_X2(1,1,iX_B0(1),iX_B0(3),iX_B0(2)  ), nDOFX_X2, Zero, &
             dV_u_dX2, nDOFX )

    ! --- Contribution from Right Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 3*nK_X, nDOFX_X2, + One, LX_X2_Up, nDOFX_X2, &
             V_u_X2(1,1,iX_B0(1),iX_B0(3),iX_B0(2)+1), nDOFX_X2, One,  &
             dV_u_dX2, nDOFX )

    ! --- Contribution from Left Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 3*nK_X, nDOFX_X2, - One, LX_X2_Dn, nDOFX_X2, &
             V_d_X2(1,1,iX_B0(1),iX_B0(3),iX_B0(2)  ), nDOFX_X2, Zero, &
             dV_d_dX2, nDOFX )

    ! --- Contribution from Right Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 3*nK_X, nDOFX_X2, + One, LX_X2_Up, nDOFX_X2, &
             V_d_X2(1,1,iX_B0(1),iX_B0(3),iX_B0(2)+1), nDOFX_X2, One,  &
             dV_d_dX2, nDOFX )

    CALL TimersStop( Timer_Streaming_LinearAlgebra )

    ! -------------------
    ! --- Volume Term ---
    ! -------------------

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5) &
    !$OMP PRIVATE( uV_K, iCF, iCF_S, iGF_Gm_dd, iGF_h )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRIVATE( uV_K, iCF, iCF_S, iGF_Gm_dd, iGF_h ) &
    !$ACC PRESENT( iX_B0, iX_E0, U_F_K, GX_K, h_d_K, &
    !$ACC          V_u_K, V_d_K, WeightsX_q )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5) &
    !$OMP PRIVATE( uV_K, iCF, iCF_S, iGF_Gm_dd, iGF_h )
#endif
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO i = 1, 3
      DO iNodeX = 1, nDOFX

        iCF       = iCF_S1       + i - 1
        iGF_Gm_dd = iGF_Gm_dd_11 + i - 1
        iGF_h     = iGF_h_1      + i - 1

        h_d_K(iNodeX,i,iX1,iX3,iX2) &
          = WeightsX_q(iNodeX) * GX_K(iNodeX,iGF_h,iX1,iX3,iX2)

        uV_K &
          = U_F_K(iNodeX,iCF,iX1,iX3,iX2) &
            / ( GX_K (iNodeX,iGF_Gm_dd,iX1,iX3,iX2) &
              * U_F_K(iNodeX,iCF_D ,iX1,iX3,iX2) )

        V_u_K(iNodeX,i,iX1,iX3,iX2) &
          = uV_K * WeightsX_q(iNodeX)

        V_d_K(iNodeX,i,iX1,iX3,iX2) &
          = uV_K * WeightsX_q(iNodeX) * GX_K(iNodeX,iGF_Gm_dd,iX1,iX3,iX2)

      END DO
      END DO

    END DO
    END DO
    END DO

    CALL TimersStart( Timer_Streaming_LinearAlgebra )

    ! --- Volume Contributions ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 3*nK_X, nDOFX, - One, dLXdX2_q, nDOFX, &
             h_d_K, nDOFX, One, dh_d_dX2, nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 3*nK_X, nDOFX, - One, dLXdX2_q, nDOFX, &
             V_u_K, nDOFX, One, dV_u_dX2, nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 3*nK_X, nDOFX, - One, dLXdX2_q, nDOFX, &
             V_d_K, nDOFX, One, dV_d_dX2, nDOFX )

    CALL TimersStop( Timer_Streaming_LinearAlgebra )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( iX_B0, iX_E0, dX2, &
    !$ACC          dh_d_dX2, dV_u_dX2, dV_d_dX2, WeightsX_q )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO i      = 1, 3
      DO iNodeX = 1, nDOFX

        dh_d_dX2(iNodeX,i,iX1,iX3,iX2) &
          = dh_d_dX2(iNodeX,i,iX1,iX3,iX2) &
              / ( WeightsX_q(iNodeX) * dX2(iX2) )

        dV_u_dX2(iNodeX,i,iX1,iX3,iX2) &
         = dV_u_dX2(iNodeX,i,iX1,iX3,iX2) &
             / ( WeightsX_q(iNodeX) * dX2(iX2) )

        dV_d_dX2(iNodeX,i,iX1,iX3,iX2) &
         = dV_d_dX2(iNodeX,i,iX1,iX3,iX2) &
             / ( WeightsX_q(iNodeX) * dX2(iX2) )

      END DO
      END DO

    END DO
    END DO
    END DO

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5) &
    !$OMP PRIVATE( iGF_h )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRIVATE( iGF_h ) &
    !$ACC PRESENT( iX_B0, iX_E0, GX_K, dGm_dd_dX2_Out, dh_d_dX2, &
    !$ACC          dV_u_dX2_Out, dV_u_dX2, dV_d_dX2_Out, dV_d_dX2 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5) &
    !$OMP PRIVATE( iGF_h )
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO i = 1, 3
      DO iNodeX = 1, nDOFX

        iGF_h = iGF_h_1 + i - 1

        dGm_dd_dX2_Out(iNodeX,i,iX1,iX2,iX3) &
          = Two * GX_K(iNodeX,iGF_h,iX1,iX3,iX2) &
              * dh_d_dX2(iNodeX,i,iX1,iX3,iX2)

        dV_u_dX2_Out(iNodeX,i,iX1,iX2,iX3) &
          = dV_u_dX2(iNodeX,i,iX1,iX3,iX2)

        dV_d_dX2_Out(iNodeX,i,iX1,iX2,iX3) &
          = dV_d_dX2(iNodeX,i,iX1,iX3,iX2)


      END DO
      END DO

    END DO
    END DO
    END DO

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: dX2, iX_B0, iX_E0, &
    !$OMP               GX_K, GX_F, h_d_F, h_d_K, dh_d_dX2, &
    !$OMP               U_F_K, U_F_L, U_F_R, V_u_X2, V_d_X2, V_u_k, V_d_k, &
    !$OMP               dV_u_dX2, dV_d_dX2 )
#elif defined( THORNADO_OACC   )
    !$ACC EXIT DATA &
    !$ACC DELETE( dX2, iX_B0, iX_E0, &
    !$ACC         GX_K, GX_F, h_d_F, h_d_K, dh_d_dX2, &
    !$ACC         U_F_K, U_F_L, U_F_R, V_u_X2, V_d_X2, V_u_k, V_d_k, &
    !$ACC         dV_u_dX2, dV_d_dX2 )
#endif

    END ASSOCIATE

  END SUBROUTINE ComputeWeakDerivatives_X2


  SUBROUTINE ComputeWeakDerivatives_X3 &
    ( iX_B0, iX_E0, iX_B1, iX_E1, GX, U_F, dV_u_dX3_Out, dV_d_dX3_Out, &
      dGm_dd_dX3_Out )

    INTEGER,  INTENT(in)  :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)  :: &
      GX (1:nDOFX, &
          iX_B1(1):iX_E1(1), &
          iX_B1(2):iX_E1(2), &
          iX_B1(3):iX_E1(3), &
          1:nGF)
    REAL(DP), INTENT(in)  :: &
      U_F(1:nDOFX, &
          iX_B1(1):iX_E1(1), &
          iX_B1(2):iX_E1(2), &
          iX_B1(3):iX_E1(3), &
          1:nCF)
    REAL(DP), INTENT(out) :: &
      dV_u_dX3_Out &
        (1:nDOFX,1:3, &
         iX_B0(1):iX_E0(1), &
         iX_B0(2):iX_E0(2), &
         iX_B0(3):iX_E0(3)), &
      dV_d_dX3_Out &
        (1:nDOFX,1:3, &
         iX_B0(1):iX_E0(1), &
         iX_B0(2):iX_E0(2), &
         iX_B0(3):iX_E0(3)), &
      dGm_dd_dX3_Out &
        (1:nDOFX,1:3, &
         iX_B0(1):iX_E0(1), &
         iX_B0(2):iX_E0(2), &
         iX_B0(3):iX_E0(3))

    INTEGER  :: iNodeX
    INTEGER  :: iX2, iX3, iX1, i
    INTEGER  :: iCF, iCF_S
    INTEGER  :: iGF, iGF_h, iGF_Gm_dd
    INTEGER  :: nX(3), nX_X3(3), nK_X, nX3_X
    REAL(DP) :: uV_L(3), uV_R(3), uV_F(3), uV_K

    ! --- Geometry Fields ---

    REAL(DP), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: GX_K, GX_F, h_d_F, h_d_K, dh_d_dX3

    ! --- Conserved Fluid Fields ---

    REAL(DP), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: U_F_K, U_F_L, U_F_R

    ! --- Velocities ---

    REAL(DP), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: V_u_X3, V_d_X3, V_u_K, V_d_K, dV_u_dX3, dV_d_dX3

    ALLOCATE( GX_K    (nDOFX   ,nGF,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3)-1:iX_E0(3)+1) )
    ALLOCATE( GX_F    (nDOFX_X3,nGF,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3)  :iX_E0(3)+1) )
    ALLOCATE( h_d_F   (nDOFX_X3,3  ,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3)  :iX_E0(3)+1) )
    ALLOCATE( h_d_K   (nDOFX   ,3  ,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3)  :iX_E0(3)  ) )
    ALLOCATE( dh_d_dX3(nDOFX   ,3  ,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3)  :iX_E0(3)  ) )
    ALLOCATE( U_F_K   (nDOFX   ,nCF,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3)-1:iX_E0(3)+1) )
    ALLOCATE( U_F_L   (nDOFX_X3,nCF,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3)  :iX_E0(3)+1) )
    ALLOCATE( U_F_R   (nDOFX_X3,nCF,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3)  :iX_E0(3)+1) )
    ALLOCATE( V_u_X3  (nDOFX_X3,3  ,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3)  :iX_E0(3)+1) )
    ALLOCATE( V_d_X3  (nDOFX_X3,3  ,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3)  :iX_E0(3)+1) )
    ALLOCATE( V_u_K   (nDOFX   ,3  ,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3)  :iX_E0(3)  ) )
    ALLOCATE( V_d_K   (nDOFX   ,3  ,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3)  :iX_E0(3)  ) )
    ALLOCATE( dV_u_dX3(nDOFX   ,3  ,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3)  :iX_E0(3)  ) )
    ALLOCATE( dV_d_dX3(nDOFX   ,3  ,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3)  :iX_E0(3)  ) )

    IF( iX_E0(3) .EQ. iX_B0(3) )THEN

#if   defined( THORNADO_OMP_OL )
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5) &
      !$OMP MAP( to: iX_B0, iX_E0 )
#elif defined( THORNADO_OACC   )
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
      !$ACC COPYIN( iX_B0, iX_E0 ) &
      !$ACC PRESENT( dV_u_dX3_Out, dV_d_dX3_Out, dGm_dd_dX3_Out )
#elif defined( THORNADO_OMP    )
      !$OMP PARALLEL DO COLLAPSE(5)
#endif
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO i = 1, 3
      DO iNodeX = 1, nDOFX

          dV_u_dX3_Out  (iNodeX,i,iX1,iX2,iX3) = Zero
          dV_d_dX3_Out  (iNodeX,i,iX1,iX2,iX3) = Zero
          dGm_dd_dX3_Out(iNodeX,i,iX1,iX2,iX3) = Zero

      END DO
      END DO
      END DO
      END DO
      END DO

      RETURN
    END IF

    nX    = iX_E0 - iX_B0 + 1 ! --- Number of Elements per Spatial Dimension
    nX_X3 = nX + [ 0, 0, 1 ]  ! --- Number of X3 Faces per Spatial Dimension
    nK_X  = PRODUCT( nX )     ! --- Number of Elements in Position Space
    nX3_X = PRODUCT( nX_X3 )  ! --- Number of X3 Faces in Position Space

    ASSOCIATE( dX3 => MeshX(3) % Width )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: dX3, iX_B0, iX_E0 ) &
    !$OMP MAP( alloc: GX_K, GX_F, h_d_F, h_d_K, dh_d_dX3, &
    !$OMP             U_F_K, U_F_L, U_F_R, V_u_X3, V_d_X3, V_u_k, V_d_k, &
    !$OMP             dV_u_dX3, dV_d_dX3 )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA &
    !$ACC COPYIN( dX3, iX_B0, iX_E0 ) &
    !$ACC CREATE( GX_K, GX_F, h_d_F, h_d_K, dh_d_dX3, &
    !$ACC         U_F_K, U_F_L, U_F_R, V_u_X3, V_d_X3, V_u_k, V_d_k, &
    !$ACC         dV_u_dX3, dV_d_dX3 )
#endif

    ! --- Permute Geometry Fields ---

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( GX_K, GX, iX_B0, iX_E0 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iX3 = iX_B0(3)-1, iX_E0(3)+1
    DO iX2 = iX_B0(2)  , iX_E0(2)
    DO iX1 = iX_B0(1)  , iX_E0(1)

      DO iGF    = 1, nGF
      DO iNodeX = 1, nDOFX

        GX_K(iNodeX,iGF,iX1,iX2,iX3) = GX(iNodeX,iX1,iX2,iX3,iGF)

      END DO
      END DO

    END DO
    END DO
    END DO

    !---------------------
    ! --- Surface Term ---
    !---------------------

    CALL TimersStart( Timer_Streaming_LinearAlgebra )

    ! --- Interpolate Geometry Fields on Shared Face ---

    ! --- Face States (Average of Left and Right States) ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X3, nX3_X*nGF, nDOFX, One,  LX_X3_Up, nDOFX_X3, &
             GX_K(1,1,iX_B0(1),iX_B0(2),iX_B0(3)-1), nDOFX, Zero, &
             GX_F(1,1,iX_B0(1),iX_B0(2),iX_B0(3)  ), nDOFX_X3 )

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X3, nX3_X*nGF, nDOFX, Half, LX_X3_Dn, nDOFX_X3, &
             GX_K(1,1,iX_B0(1),iX_B0(2),iX_B0(3)  ), nDOFX, Half, &
             GX_F(1,1,iX_B0(1),iX_B0(2),iX_B0(3)  ), nDOFX_X3 )

    CALL TimersStop( Timer_Streaming_LinearAlgebra )

    ! --- Compute Metric Components from Scale Factors ---

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( GX_F, h_d_F, iX_B0, iX_E0 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(4)
#endif
    DO iX3  = iX_B0(3), iX_E0(3)+1
    DO iX2  = iX_B0(2), iX_E0(2)
    DO iX1  = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX_X3

        GX_F(iNodeX,iGF_Gm_dd_11,iX1,iX2,iX3) &
          = MAX( GX_F(iNodeX,iGF_h_1,iX1,iX2,iX3)**2, SqrtTiny )
        GX_F(iNodeX,iGF_Gm_dd_22,iX1,iX2,iX3) &
          = MAX( GX_F(iNodeX,iGF_h_2,iX1,iX2,iX3)**2, SqrtTiny )
        GX_F(iNodeX,iGF_Gm_dd_33,iX1,iX2,iX3) &
          = MAX( GX_F(iNodeX,iGF_h_3,iX1,iX2,iX3)**2, SqrtTiny )
        GX_F(iNodeX,iGF_SqrtGm,iX1,iX2,iX3) &
          = SQRT(   GX_F(iNodeX,iGF_Gm_dd_11,iX1,iX2,iX3) &
                  * GX_F(iNodeX,iGF_Gm_dd_22,iX1,iX2,iX3) &
                  * GX_F(iNodeX,iGF_Gm_dd_33,iX1,iX2,iX3) )

        h_d_F(iNodeX,1,iX1,iX2,iX3) &
          = GX_F(iNodeX,iGF_h_1,iX1,iX2,iX3) * WeightsX_X3(iNodeX)
        h_d_F(iNodeX,2,iX1,iX2,iX3) &
          = GX_F(iNodeX,iGF_h_2,iX1,iX2,iX3) * WeightsX_X3(iNodeX)
        h_d_F(iNodeX,3,iX1,iX2,iX3) &
          = GX_F(iNodeX,iGF_h_3,iX1,iX2,iX3) * WeightsX_X3(iNodeX)

      END DO

    END DO
    END DO
    END DO

    CALL TimersStart( Timer_Streaming_LinearAlgebra )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 3*nK_X, nDOFX_X3, - One, LX_X3_Dn, nDOFX_X3, &
             h_d_F(1,1,iX_B0(1),iX_B0(2),iX_B0(3)  ), nDOFX_X3, Zero, &
             dh_d_dX3, nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 3*nK_X, nDOFX_X3, + One, LX_X3_Up, nDOFX_X3, &
             h_d_F(1,1,iX_B0(1),iX_B0(2),iX_B0(3)+1), nDOFX_X3, One,  &
             dh_d_dX3, nDOFX )

    CALL TimersStop( Timer_Streaming_LinearAlgebra )

    ! --- Permute Fluid Fields ---

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( U_F_K, U_F, iX_B0, iX_E0 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iX3 = iX_B0(3)-1, iX_E0(3)+1
    DO iX2 = iX_B0(2)  , iX_E0(2)
    DO iX1 = iX_B0(1)  , iX_E0(1)

      DO iCF    = 1, nCF
      DO iNodeX = 1, nDOFX

        U_F_K(iNodeX,iCF,iX1,iX2,iX3) = U_F(iNodeX,iX1,iX2,iX3,iCF)

      END DO
      END DO

    END DO
    END DO
    END DO

    ! --- Interpolate Fluid Fields ---

    CALL TimersStart( Timer_Streaming_LinearAlgebra )

      ! --- Interpolate Left State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X3, nX3_X*nCF, nDOFX, One, LX_X3_Up, nDOFX_X3, &
             U_F_K(1,1,iX_B0(1),iX_B0(2),iX_B0(3)-1), nDOFX, Zero, &
             U_F_L(1,1,iX_B0(1),iX_B0(2),iX_B0(3)  ), nDOFX_X3 )

    ! --- Interpolate Right State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X3, nX3_X*nCF, nDOFX, One, LX_X3_Dn, nDOFX_X3, &
             U_F_K(1,1,iX_B0(1),iX_B0(2),iX_B0(3)  ), nDOFX, Zero, &
             U_F_R(1,1,iX_B0(1),iX_B0(2),iX_B0(3)  ), nDOFX_X3 )

    CALL TimersStop( Timer_Streaming_LinearAlgebra )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4) &
    !$OMP PRIVATE( uV_L, uV_R, uV_F, iCF_S, iGF_Gm_dd )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRIVATE( uV_L, uV_R, uV_F, iCF_S, iGF_Gm_dd ) &
    !$ACC PRESENT( iX_B0, iX_E0, U_F_L, U_F_R, GX_F, &
    !$ACC          V_u_X3, V_d_X3, WeightsX_X3 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(4) &
    !$OMP PRIVATE( uV_L, uV_R, uV_F, iCF_S, iGF_Gm_dd )
#endif
    DO iX3 = iX_B0(3), iX_E0(3)+1
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX_X3

        DO i = 1, 3

          iCF_S     = iCF_S1       + i - 1
          iGF_Gm_dd = iGF_Gm_dd_11 + i - 1

          ! --- Left States ---

          uV_L(i) = U_F_L(iNodeX,iCF_S,iX1,iX2,iX3) &
                    / ( GX_F (iNodeX,iGF_Gm_dd,iX1,iX2,iX3) &
                         * U_F_L(iNodeX,iCF_D ,iX1,iX2,iX3) )

          ! --- Right States ---

          uV_R(i) = U_F_R(iNodeX,iCF_S,iX1,iX2,iX3) &
                    / ( GX_F (iNodeX,iGF_Gm_dd,iX1,iX2,iX3) &
                         * U_F_R(iNodeX,iCF_D ,iX1,iX2,iX3) )

        END DO

        CALL FaceVelocity_X3 &
               ( uV_L(1), uV_L(2), uV_L(3), &
                 uV_R(1), uV_R(2), uV_R(3), &
                 uV_F(1), uV_F(2), uV_F(3) )

        DO i = 1, 3

          iGF_Gm_dd = iGF_Gm_dd_11 + i - 1

          V_u_X3(iNodeX,i,iX1,iX2,iX3) &
            = uV_F(i) * WeightsX_X3(iNodeX)

          V_d_X3(iNodeX,i,iX1,iX2,iX3) &
            = uV_F(i) * WeightsX_X3(iNodeX) &
                * GX_F(iNodeX,iGF_Gm_dd,iX1,iX2,iX3)

        END DO

      END DO

    END DO
    END DO
    END DO

    CALL TimersStart( Timer_Streaming_LinearAlgebra )

    ! --- Surface Contributions ---

    ! --- Contribution from Left Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 3*nK_X, nDOFX_X3, - One, LX_X3_Dn, nDOFX_X3, &
             V_u_X3(1,1,iX_B0(1),iX_B0(2),iX_B0(3)  ), nDOFX_X3, Zero, &
             dV_u_dX3, nDOFX )

    ! --- Contribution from Right Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 3*nK_X, nDOFX_X3, + One, LX_X3_Up, nDOFX_X3, &
             V_u_X3(1,1,iX_B0(1),iX_B0(2),iX_B0(3)+1), nDOFX_X3, One,  &
             dV_u_dX3, nDOFX )

    ! --- Contribution from Left Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 3*nK_X, nDOFX_X3, - One, LX_X3_Dn, nDOFX_X3, &
             V_d_X3(1,1,iX_B0(1),iX_B0(2),iX_B0(3)  ), nDOFX_X3, Zero, &
             dV_d_dX3, nDOFX )

    ! --- Contribution from Right Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 3*nK_X, nDOFX_X3, + One, LX_X3_Up, nDOFX_X3, &
             V_d_X3(1,1,iX_B0(1),iX_B0(2),iX_B0(3)+1), nDOFX_X3, One,  &
             dV_d_dX3, nDOFX )

    CALL TimersStop( Timer_Streaming_LinearAlgebra )

    ! -------------------
    ! --- Volume Term ---
    ! -------------------

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5) &
    !$OMP PRIVATE( uV_K, iCF, iCF_S, iGF_Gm_dd, iGF_h )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRIVATE( uV_K, iCF, iCF_S, iGF_Gm_dd, iGF_h ) &
    !$ACC PRESENT( iX_B0, iX_E0, U_F_K, GX_K, h_d_K, &
    !$ACC          V_u_K, V_d_K, WeightsX_q )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5) &
    !$OMP PRIVATE( uV_K, iCF, iCF_S, iGF_Gm_dd, iGF_h )
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO i = 1, 3
      DO iNodeX = 1, nDOFX

        iCF       = iCF_S1       + i - 1
        iGF_Gm_dd = iGF_Gm_dd_11 + i - 1
        iGF_h     = iGF_h_1      + i - 1

        h_d_K(iNodeX,i,iX1,iX2,iX3) &
          = WeightsX_q(iNodeX) * GX_K(iNodeX,iGF_h,iX1,iX2,iX3)

        uV_K &
          = U_F_K(iNodeX,iCF,iX1,iX2,iX3) &
            / ( GX_K (iNodeX,iGF_Gm_dd,iX1,iX2,iX3) &
                 * U_F_K(iNodeX,iCF_D ,iX1,iX2,iX3) )

        V_u_K(iNodeX,i,iX1,iX2,iX3) &
          = uV_K * WeightsX_q(iNodeX)

        V_d_K(iNodeX,i,iX1,iX2,iX3) &
          = uV_K * WeightsX_q(iNodeX) * GX_K(iNodeX,iGF_Gm_dd,iX1,iX2,iX3)

      END DO
      END DO

    END DO
    END DO
    END DO

    CALL TimersStart( Timer_Streaming_LinearAlgebra )

    ! --- Volume Contributions ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 3*nK_X, nDOFX, - One, dLXdX3_q, nDOFX, &
             h_d_K, nDOFX, One, dh_d_dX3, nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 3*nK_X, nDOFX, - One, dLXdX3_q, nDOFX, &
             V_u_K, nDOFX, One, dV_u_dX3, nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 3*nK_X, nDOFX, - One, dLXdX3_q, nDOFX, &
             V_d_K, nDOFX, One, dV_d_dX3, nDOFX )

    CALL TimersStop( Timer_Streaming_LinearAlgebra )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( iX_B0, iX_E0, dX3, &
    !$ACC          dh_d_dX3, dV_u_dX3, dV_d_dX3, WeightsX_q )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO i      = 1, 3
      DO iNodeX = 1, nDOFX

        dh_d_dX3(iNodeX,i,iX1,iX2,iX3) &
          = dh_d_dX3(iNodeX,i,iX1,iX2,iX3) &
              / ( WeightsX_q(iNodeX) * dX3(iX3) )

        dV_u_dX3(iNodeX,i,iX1,iX2,iX3) &
         = dV_u_dX3(iNodeX,i,iX1,iX2,iX3) &
             / ( WeightsX_q(iNodeX) * dX3(iX3) )

        dV_d_dX3(iNodeX,i,iX1,iX2,iX3) &
         = dV_d_dX3(iNodeX,i,iX1,iX2,iX3) &
             / ( WeightsX_q(iNodeX) * dX3(iX3) )

      END DO
      END DO

    END DO
    END DO
    END DO

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5) &
    !$OMP PRIVATE( iGF_h )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRIVATE( iGF_h ) &
    !$ACC PRESENT( iX_B0, iX_E0, GX_K, dGm_dd_dX3_Out, dh_d_dX3, &
    !$ACC          dV_u_dX3_Out, dV_u_dX3, dV_d_dX3_Out, dV_d_dX3 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5) &
    !$OMP PRIVATE( iGF_h )
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO i      = 1, 3
      DO iNodeX = 1, nDOFX

        iGF_h = iGF_h_1 + i - 1

        dGm_dd_dX3_Out(iNodeX,i,iX1,iX2,iX3) &
          = Two * GX_K(iNodeX,iGF_h,iX1,iX2,iX3) &
              * dh_d_dX3(iNodeX,i,iX1,iX2,iX3)

        dV_u_dX3_Out(iNodeX,i,iX1,iX2,iX3) &
          = dV_u_dX3(iNodeX,i,iX1,iX2,iX3)

        dV_d_dX3_Out(iNodeX,i,iX1,iX2,iX3) &
          = dV_d_dX3(iNodeX,i,iX1,iX2,iX3)

      END DO
      END DO

    END DO
    END DO
    END DO

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: dX3, iX_B0, iX_E0, &
    !$OMP               GX_K, GX_F, h_d_F, h_d_K, dh_d_dX3, &
    !$OMP               U_F_K, U_F_L, U_F_R, V_u_X3, V_d_X3, V_u_k, V_d_k, &
    !$OMP               dV_u_dX3, dV_d_dX3 )
#elif defined( THORNADO_OACC   )
    !$ACC EXIT DATA &
    !$ACC DELETE( dX3, iX_B0, iX_E0, &
    !$ACC         GX_K, GX_F, h_d_F, h_d_K, dh_d_dX3, &
    !$ACC         U_F_K, U_F_L, U_F_R, V_u_X3, V_d_X3, V_u_k, V_d_k, &
    !$ACC         dV_u_dX3, dV_d_dX3 )
#endif

    END ASSOCIATE

  END SUBROUTINE ComputeWeakDerivatives_X3


  SUBROUTINE FaceVelocity_X1 &
    ( V1_L, V2_L, V3_L, V1_R, V2_R, V3_R, V1_F, V2_F, V3_F )

#if   defined( THORNADO_OMP_OL )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC   )
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: V1_L, V2_L, V3_L
    REAL(DP), INTENT(in)  :: V1_R, V2_R, V3_R
    REAL(DP), INTENT(out) :: V1_F, V2_F, V3_F

    ! --- Average Left and Right States ---

    V1_F = Half * ( V1_L + V1_R )
    V2_F = Half * ( V2_L + V2_R )
    V3_F = Half * ( V3_L + V3_R )

    RETURN
  END SUBROUTINE FaceVelocity_X1


  SUBROUTINE FaceVelocity_X2 &
    ( V1_L, V2_L, V3_L, V1_R, V2_R, V3_R, V1_F, V2_F, V3_F )

#if   defined( THORNADO_OMP_OL )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC   )
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: V1_L, V2_L, V3_L
    REAL(DP), INTENT(in)  :: V1_R, V2_R, V3_R
    REAL(DP), INTENT(out) :: V1_F, V2_F, V3_F

    ! --- Average Left and Right States ---

    V1_F = Half * ( V1_L + V1_R )
    V2_F = Half * ( V2_L + V2_R )
    V3_F = Half * ( V3_L + V3_R )

    RETURN
  END SUBROUTINE FaceVelocity_X2


  SUBROUTINE FaceVelocity_X3 &
    ( V1_L, V2_L, V3_L, V1_R, V2_R, V3_R, V1_F, V2_F, V3_F )

#if   defined( THORNADO_OMP_OL )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC   )
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: V1_L, V2_L, V3_L
    REAL(DP), INTENT(in)  :: V1_R, V2_R, V3_R
    REAL(DP), INTENT(out) :: V1_F, V2_F, V3_F

    ! --- Average Left and Right States ---

    V1_F = Half * ( V1_L + V1_R )
    V2_F = Half * ( V2_L + V2_R )
    V3_F = Half * ( V3_L + V3_R )

    RETURN
  END SUBROUTINE FaceVelocity_X3


END MODULE TwoMoment_UtilitiesModule
