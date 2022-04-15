MODULE TwoMoment_UtilitiesModule_Relativistic

  USE KindModule, ONLY: &
    DP, Zero, Half, One, Two, Three, Five, SqrtTiny
  USE ProgramHeaderModule, ONLY: &
    nDOFZ, nDOFX, nDOFE
  USE GeometryFieldsModule, ONLY: &
    nGF, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_Alpha, &
    iGF_Beta_1, &
    iGF_Beta_2, &
    iGF_Beta_3
  USE FluidFieldsModule, ONLY: &
    nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne
    USE Euler_UtilitiesModule_Relativistic, ONLY: &
    ComputePrimitive_Euler_Relativistic
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3, &
    nPR, iPR_D, iPR_I1, iPR_I2, iPR_I3
  USE TwoMoment_ClosureModule, ONLY: &
    FluxFactor_Relativistic, &
    EddingtonFactor, &
    HeatFluxFactor

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputePrimitive_TwoMoment
  PUBLIC :: ComputePrimitive_TwoMoment_Vector_Richardson
  PUBLIC :: ComputeConserved_TwoMoment
  PUBLIC :: ComputeFromConserved_TwoMoment
  PUBLIC :: ComputeEddingtonTensorComponents_dd
  PUBLIC :: ComputeEddingtonTensorComponents_uu
  PUBLIC :: ComputeEddingtonTensorComponents_ud
  PUBLIC :: ComputeHeatFluxTensorComponents_ddd_Lagrangian 
  PUBLIC :: ComputeHeatFluxTensorComponents_uud_Lagrangian
  PUBLIC :: ComputeHeatFluxTensorComponents_udd_Lagrangian
  PUBLIC :: Flux_X1
  PUBLIC :: Flux_X2
  PUBLIC :: Flux_X3
  PUBLIC :: Flux_E
  PUBLIC :: Source_E
  PUBLIC :: NumericalFlux_LLF

CONTAINS

  SUBROUTINE ComputePrimitive_TwoMoment &
    ( N, G_d_1, G_d_2, G_d_3, D, I_u_1, I_u_2, I_u_3, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, Gm_dd_12, Gm_dd_13, Gm_dd_23, alp, B_u_1, B_u_2, B_u_3, nIterations_Option )

    REAL(DP), INTENT(in)  :: N, G_d_1, G_d_2, G_d_3 ! --- Index Down
    REAL(DP), INTENT(out) :: D, I_u_1, I_u_2, I_u_3 ! --- Index Up
    REAL(DP), INTENT(in)  ::    V_u_1, V_u_2, V_u_3 ! --- Index Up
    REAL(DP), INTENT(in)  :: Gm_dd_11, Gm_dd_22, Gm_dd_33, Gm_dd_12, Gm_dd_13, Gm_dd_23
    REAL(DP), INTENT(in)  :: B_u_1, B_u_2, B_u_3, alp
    INTEGER, INTENT(out), OPTIONAL :: nIterations_Option

    ! --- Parameters ---

    INTEGER,  PARAMETER :: M = 2
    INTEGER,  PARAMETER :: LWORK = 2 * M
    INTEGER,  PARAMETER :: MaxIterations = 100
    REAL(DP), PARAMETER :: Rtol = 1.0d-08

    ! --- Local Variables ---

    LOGICAL  :: CONVERGED
    INTEGER  :: i, k, mk, INFO
    REAL(DP) :: I_d_1, I_d_2, I_d_3, A_d_1, A_d_2, A_d_3
    REAL(DP) :: k_dd_ij(1:3,1:3)
    REAL(DP) :: B_d_1, B_d_2, B_d_3
    REAL(DP) :: UVEC(4), CVEC(4)
    REAL(DP) :: GVEC(4,M), GVECm(4)
    REAL(DP) :: FVEC(4,M), FVECm(4)
    REAL(DP) :: LMAT(4,4), DET, Alpha(M)
    REAL(DP) :: BVEC(4), AMAT(4,M), WORK(LWORK)
    REAL(DP) :: W, DT
    
    B_d_1 = Gm_dd_11 * B_u_1
    B_d_2 = Gm_dd_22 * B_u_2
    B_d_3 = Gm_dd_33 * B_u_3


    W = 1.0_DP / SQRT( 1.0_DP - (Gm_dd_11 * V_u_1 * V_u_1 &
               + Gm_dd_22 * V_u_2 * V_u_2 &  
               + Gm_dd_33 * V_u_3 * V_u_3) )

    DT = 1.0_DP / ( B_d_1 * V_u_1 + B_d_2 * V_u_2 + B_d_3 * V_u_3 - alp )

    CVEC = [ N, G_d_1, G_d_2, G_d_3 ]

    ! --- Initial Guess ---
    D     = N
    I_u_1 = G_d_1
    I_u_2 = G_d_2
    I_u_3 = G_d_3


    I_d_1 = DT * ( B_d_2 * V_u_2 + B_d_3 * V_u_3 - alp ) * Gm_dd_11 * I_u_1 &
          - DT * ( B_d_1 * V_u_2 *Gm_dd_22 ) * I_u_2 - DT * ( B_d_1 * V_u_3 * Gm_dd_33 ) * I_u_3 
    I_d_2 = DT * ( B_d_1 * V_u_1 + B_d_3 * V_u_3 - alp ) * Gm_dd_22 * I_u_2 &
          - DT * ( B_d_2 * V_u_1 * Gm_dd_11 ) * I_u_1 - DT * ( Gm_dd_33 * I_u_3 * B_d_2 * V_u_3 ) 
    I_d_3 = DT * ( B_d_1 * V_u_1 + B_d_2 * V_u_2 - alp ) * Gm_dd_33 * I_u_3 &
          - DT * ( Gm_dd_11 * I_u_1 * B_d_3 * V_u_1 ) - DT * ( Gm_dd_22 * I_u_2 * B_d_3 * V_u_2 )
    k = 0
    CONVERGED = .FALSE.
    DO WHILE( .NOT. CONVERGED .AND. k < MaxIterations )

      k  = k + 1
      mk = MIN( M, k )

      UVEC = [ D, I_d_1, I_d_2, I_d_3 ]

      CALL ComputeEddingtonTensorComponents_dd &
             ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
               alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3, k_dd_ij   )

      A_d_1 = V_u_1 * k_dd_ij(1,1) + V_u_2 * k_dd_ij(1,2) + V_u_3 * k_dd_ij(1,3)
      A_d_2 = V_u_1 * k_dd_ij(1,2) + V_u_2 * k_dd_ij(2,2) + V_u_3 * k_dd_ij(2,3)
      A_d_3 = V_u_1 * k_dd_ij(1,3) + V_u_2 * k_dd_ij(2,3) + V_u_3 * k_dd_ij(3,3)

      DET = W**2 - ( V_u_1 * A_d_1 + V_u_2 * A_d_2 + V_u_3 * A_d_3 )

      LMAT(1,1) = W
      LMAT(2,1) = - A_d_1
      LMAT(3,1) = - A_d_2
      LMAT(4,1) = - A_d_3

      LMAT(1,2) = - V_u_1
      LMAT(2,2) = W - ( V_u_2 * A_d_2 + V_u_3 * A_d_3 ) / W
      LMAT(3,2) = V_u_1 * A_d_2 / W
      LMAT(4,2) = V_u_1 * A_d_3 / W

      LMAT(1,3) = - V_u_2
      LMAT(2,3) = V_u_2 * A_d_1 / W
      LMAT(3,3) = W - ( V_u_1 * A_d_1 + V_u_3 * A_d_3 ) / W
      LMAT(4,3) = V_u_2 * A_d_3 / W

      LMAT(1,4) = - V_u_3
      LMAT(2,4) = V_u_3 * A_d_1 / W
      LMAT(3,4) = V_u_3 * A_d_2 / W
      LMAT(4,4) = W - ( V_u_1 * A_d_1 + V_u_2 * A_d_2 ) / W

      LMAT = LMAT / DET

      CALL DGEMV( 'N', 4, 4, One, LMAT, 4, CVEC, 1, Zero, GVEC(:,mk), 1 )

      FVEC(:,mk) = GVEC(:,mk) - UVEC

      IF( mk == 1 )THEN

        ! --- Picard Iteration ---

        GVECm = GVEC(:,mk)

      ELSE

        ! --- Anderson Accelerated Fixed-Point ---

        BVEC = - FVEC(:,mk)

        AMAT(:,1:mk-1) &
          = FVEC(:,1:mk-1) - SPREAD( FVEC(:,mk), DIM = 2, NCOPIES = mk-1 )

        CALL DGELS( 'N', 4, mk-1, 1, AMAT(:,1:mk-1), 4, BVEC, 4, &
                    WORK, LWORK, INFO )

        Alpha(1:mk-1) = BVEC(1:mk-1)
        Alpha(mk)     = One - SUM( Alpha(1:mk-1) )

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

        GVEC = CSHIFT( GVEC, SHIFT = + 1, DIM = 2 )
        FVEC = CSHIFT( FVEC, SHIFT = + 1, DIM = 2 )

      END IF

     
      D     = UVEC(1)
      I_d_1 = UVEC(2); I_u_1 = ( 1.0_DP - B_d_1 * V_u_1 / alp ) * I_d_1 / Gm_dd_11  &
            - I_d_2 * B_d_1 * V_u_2 / ( alp *Gm_dd_11 ) - I_d_3 * B_d_1 * V_u_3 / ( Gm_dd_11 * alp )
      I_d_2 = UVEC(3); I_u_2 = ( 1.0_DP - B_d_2 * V_u_2 / alp ) * I_d_2 / Gm_dd_22  &
            - I_d_1 * B_d_2 * V_u_1 / ( alp *Gm_dd_22 ) - I_d_3 * B_d_2 * V_u_3 / ( Gm_dd_22 * alp )
      I_d_3 = UVEC(4); I_u_3 = ( 1.0_DP - B_d_3 * V_u_3 / alp ) * I_d_3 / Gm_dd_33  &
            - I_d_1 * B_d_3 * V_u_1 / ( alp *Gm_dd_33 ) - I_d_2 * B_d_3 * V_u_2 / ( Gm_dd_33 * alp )
    END DO

    IF( PRESENT( nIterations_Option ) )THEN

      nIterations_Option = k

    END IF

    IF( k == MaxIterations )THEN

      PRINT*
      PRINT*, "ComputePrimitive_TwoMoment"
      PRINT*
      PRINT*, "  N     = ", N
      PRINT*, "  G_d_1 = ", G_d_1
      PRINT*, "  G_d_2 = ", G_d_2
      PRINT*, "  G_d_3 = ", G_d_3
      PRINT*
      PRINT*, "  V_u_1 = ", V_u_1
      PRINT*, "  V_u_2 = ", V_u_2
      PRINT*, "  V_u_3 = ", V_u_3
      PRINT*

      PRINT*, "  Converged with k = ", k

      PRINT*
      PRINT*, "  FVECm = ", FVECm
      PRINT*

      PRINT*
      PRINT*, "  D     = ", D
      PRINT*, "  I_u_1 = ", I_u_1
      PRINT*, "  I_u_2 = ", I_u_2
      PRINT*, "  I_u_3 = ", I_u_3
      PRINT*
    END IF

  END SUBROUTINE ComputePrimitive_TwoMoment

  SUBROUTINE ComputePrimitive_TwoMoment_Vector_Richardson &
    ( N, G_d_1, G_d_2, G_d_3, D, I_u_1, I_u_2, I_u_3, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      alp, B_u_1, B_u_2, B_u_3, PositionIndexZ, nIterations_Option )

    REAL(DP), DIMENSION(:), INTENT(in)  :: N, G_d_1, G_d_2, G_d_3
    REAL(DP), DIMENSION(:), INTENT(out) :: D, I_u_1, I_u_2, I_u_3
    REAL(DP), DIMENSION(:), INTENT(in)  :: V_u_1, V_u_2, V_u_3
    REAL(DP), DIMENSION(:), INTENT(in)  :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), DIMENSION(:), INTENT(in)  :: B_u_1, B_u_2, B_u_3, alp
    INTEGER,  DIMENSION(:), INTENT(in)  :: PositionIndexZ
    INTEGER,  DIMENSION(:), INTENT(out), OPTIONAL :: nIterations_Option


  ! --- Parameters ---

    INTEGER,  PARAMETER :: M = 2
    INTEGER,  PARAMETER :: MaxIterations = 1000
    REAL(DP), PARAMETER :: Rtol = 1.0d-08

    ! --- Local Variables ---

    INTEGER  :: nZ
    INTEGER  :: iX, iZ
    INTEGER  :: k, Mk, iM, i, j

    REAL(DP) :: FTMP(4,M), GTMP(4,M)
    REAL(DP) :: SUM1, k_dd(3,3)
    REAL(DP) :: B_d_1, B_d_2, B_d_3
    REAL(DP) :: I_d_1, I_d_2, I_d_3
    REAL(DP) :: vK, vI
    REAL(DP) :: W, DT, absV
    REAL(DP) :: Omega
    LOGICAL  :: CONVERGED

    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: FVEC, GVEC
    REAL(DP), DIMENSION(:,:),   ALLOCATABLE :: CVEC, UVEC, FVECm, GVECm, Alpha
    LOGICAL,  DIMENSION(:),     ALLOCATABLE :: ITERATE
    INTEGER,  DIMENSION(:),     ALLOCATABLE :: nIterations

    nZ = SIZE( N, 1 )


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
    !$ACC PRESENT( CVEC, N, G_d_1, G_d_2, G_d_3, &
    !$ACC          D, I_u_1, I_u_2, I_u_3 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO
#endif
    DO iZ = 1, nZ
      CVEC(iCR_N ,iZ) = N    (iZ)
      CVEC(iCR_G1,iZ) = G_d_1(iZ)
      CVEC(iCR_G2,iZ) = G_d_2(iZ)
      CVEC(iCR_G3,iZ) = G_d_3(iZ)

      D    (iZ) = N(iZ)
    !  I_u_1(iZ) = G_d_1(iZ)
    !  I_u_2(iZ) = G_d_2(iZ)
    !  I_u_3(iZ) = G_d_3(iZ)

      I_u_1(iZ) = 0.0_DP
      I_u_2(iZ) = 0.0_DP
      I_u_3(iZ) = 0.0_DP

    END DO


    k = 0
    DO WHILE( ANY( ITERATE ) .AND. k < MaxIterations )

      k = k + 1
      Mk = MIN( M, k )


#if   defined( THORNADO_OMP_OL )
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
      !$OMP PRIVATE( iX, k_dd, absV, Omega, vI, vK )
#elif defined( THORNADO_OACC   )
      !$ACC PARALLEL LOOP GANG VECTOR ASYNC &
      !$ACC PRIVATE( iX, k_dd, absV, Omega, vI, vK, I_d_1, I_d_2, I_d_3, &
      !$ACC          B_d_1, B_d_2, B_d_3, W, DT ) &
      !$ACC PRESENT( ITERATE, UVEC, CVEC, GVEC, FVEC, GVECm, FVECm, &
      !$ACC          PositionIndexZ, D, I_u_1, I_u_2, I_u_3, &
      !$ACC          Gm_dd_11, Gm_dd_22, Gm_dd_33, V_u_1, V_u_2, V_u_3, &
      !$ACC          B_u_1, B_u_2, B_u_3, alp )
#elif defined( THORNADO_OMP    )
      !$OMP PARALLEL DO &
      !$OMP PRIVATE( iX, k_dd, absV, Omega, vI, vK, I_d_1, I_d_2, I_d_3, &
      !$OMP          B_d_1, B_d_2, B_d_3, W, DT ) 
#endif
      DO iZ = 1, nZ
        IF ( ITERATE(iZ) ) THEN

          iX = PositionIndexZ(iZ)

          B_d_1 = Gm_dd_11(iX) * B_u_1(iX)
          B_d_2 = Gm_dd_22(iX) * B_u_2(iX)
          B_d_3 = Gm_dd_33(iX) * B_u_3(iX)

          DT = 1.0_DP / ( B_d_1 * V_u_1(iX) + B_d_2 * V_u_2(iX) + B_d_3 * V_u_3(iX) - alp(iX) )



          I_d_1 = DT * ( B_d_2 * V_u_2(iX) + B_d_3 * V_u_3(iX) - alp(iX) ) * Gm_dd_11(iX) * I_u_1(iZ) &
          - DT * ( B_d_1 * V_u_2(iX) *Gm_dd_22(iX) ) * I_u_2(iZ) - DT * ( B_d_1 * V_u_3(iX) * Gm_dd_33(iX) ) * I_u_3(iZ) 
          I_d_2 = DT * ( B_d_1 * V_u_1(iX) + B_d_3 * V_u_3(iX) - alp(iX) ) * Gm_dd_22(iX) * I_u_2(iZ) &
          - DT * ( B_d_2 * V_u_1(iX) * Gm_dd_11(iX) ) * I_u_1(iZ) - DT * ( Gm_dd_33(iX) * I_u_3(iZ) * B_d_2 * V_u_3(iX) ) 
          I_d_3 = DT * ( B_d_1 * V_u_1(iX) + B_d_2 * V_u_2(iX) - alp(iX) ) * Gm_dd_33(iX) * I_u_3(iZ) &
          - DT * ( Gm_dd_11(iX) * I_u_1(iZ) * B_d_3 * V_u_1(iX) ) - DT * ( Gm_dd_22(iX) * I_u_2(iZ) * B_d_3 * V_u_2(iX) )
    
          UVEC(iPR_D ,iZ) = D    (iZ)
          UVEC(iPR_I1,iZ) = I_d_1
          UVEC(iPR_I2,iZ) = I_d_2
          UVEC(iPR_I3,iZ) = I_d_3

          CALL ComputeEddingtonTensorComponents_dd &
             ( D(iZ), I_u_1(iZ), I_u_2(iZ), I_u_3(iZ), Gm_dd_11(iX), Gm_dd_22(iX), Gm_dd_33(iX), &
               alp(iX), B_u_1(iX), B_u_2(iX), B_u_3(iX), V_u_1(iX), V_u_2(iX), V_u_3(iX), k_dd   )


          absV = SQRT(   Gm_dd_11(iX) * V_u_1(iX) * V_u_1(iX) &
                       + Gm_dd_22(iX) * V_u_2(iX) * V_u_2(iX) &  
                       + Gm_dd_33(iX) * V_u_3(iX) * V_u_3(iX) ) 

          W = 1.0_DP / SQRT( 1.0_DP - absV**2)




          Omega = 1.0_DP / ( absV + W )

 
          vI =   V_u_1(iX) * UVEC(iPR_I1,iZ) &
               + V_u_2(iX) * UVEC(iPR_I2,iZ) &
               + V_u_3(iX) * UVEC(iPR_I3,iZ)


          GVECm(1,iZ) = (One - W * Omega) * UVEC(iPR_D,iZ) &
                        + Omega * ( CVEC(iCR_N,iZ) - vI )

          DO j = 1, 3

            vK =   V_u_1(iX) * k_dd(j,1) &
                 + V_u_2(iX) * k_dd(j,2) &
                 + V_u_3(iX) * k_dd(j,3)

            GVECm(j+1,iZ) = (One - W * Omega) * UVEC(j+1,iZ) &
                            + Omega * ( CVEC(j+1,iZ) - vK * UVEC(iPR_D,iZ) )

          END DO

          DO i = 1, 4

            FVECm(i,iZ) = GVECm(i,iZ) - UVEC(i,iZ)

            GVEC(i,Mk,iZ) = GVECm(i,iZ)
            FVEC(i,Mk,iZ) = FVECm(i,iZ)

          END DO


        END IF
      END DO



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



#if   defined( THORNADO_OMP_OL )
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
      !$OMP PRIVATE( iX, CONVERGED, FTMP, GTMP )
#elif defined( THORNADO_OACC   )
      !$ACC PARALLEL LOOP GANG VECTOR ASYNC &
      !$ACC PRIVATE( iX, CONVERGED, FTMP, GTMP, B_d_1, B_d_2, B_d_3 ) &
      !$ACC PRESENT( ITERATE, UVEC, CVEC, GVECm, FVECm, GVEC, FVEC, &
      !$ACC          PositionIndexZ, D, I_u_1, I_u_2, I_u_3, &
      !$ACC          Gm_dd_11, Gm_dd_22, Gm_dd_33, B_u_1, B_u_2, B_u_3, alp, nIterations )
#elif defined( THORNADO_OMP    )
      !$OMP PARALLEL DO &
      !$OMP PRIVATE( iX, CONVERGED, FTMP, GTMP, B_d_1, B_d_2, B_d_3 )
#endif
      DO iZ = 1, nZ
        IF ( ITERATE(iZ) ) THEN

          iX = PositionIndexZ(iZ)

          B_d_1 = Gm_dd_11(iX) * B_u_1(iX)
          B_d_2 = Gm_dd_22(iX) * B_u_2(iX)
          B_d_3 = Gm_dd_33(iX) * B_u_3(iX)

          D    (iZ) = GVECm(iPR_D ,iZ)
          I_u_1(iZ) = ( 1.0_DP - B_d_1 * V_u_1(iX) / alp(iX) ) * GVECm(iPR_I1,iZ) / Gm_dd_11(iX)  &
                    - GVECm(iPR_I2,iZ) * B_d_1 * V_u_2(iX) / ( alp(iX) *Gm_dd_11(iX) ) &
                    - GVECm(iPR_I3,iZ) * B_d_1 * V_u_3(iX) / ( Gm_dd_11(iX) * alp(iX) )

          I_u_2(iZ) =( 1.0_DP - B_d_2 * V_u_2(iX) / alp(iX) ) * GVECm(iPR_I2,iZ) / Gm_dd_22(iX)  &
                    - GVECm(iPR_I1,iZ) * B_d_2 * V_u_1(iX) / ( alp(iX) *Gm_dd_22(iX) ) &
                    - GVECm(iPR_I3,iZ) * B_d_2 * V_u_3(iX) / ( Gm_dd_22(iX) * alp(iX) )

          I_u_3(iZ) =( 1.0_DP - B_d_3 * V_u_3(iX) / alp(iX) ) * GVECm(iPR_I3,iZ) / Gm_dd_33(iX)  &
                    - GVECm(iPR_I1,iZ) * B_d_3 * V_u_1(iX) / ( alp(iX) *Gm_dd_33(iX) ) &
                    - GVECm(iPR_I2,iZ) * B_d_3 * V_u_2(iX) / ( Gm_dd_33(iX) * alp(iX) )


          CONVERGED = SQRT( SUM( FVECm(:,iZ)**2 ) ) <= &
                                 Rtol * SQRT( SUM( CVEC(:,iZ)**2 ) )



          !CONVERGED = (k==MaxIterations)

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


    END DO


    IF( PRESENT( nIterations_Option ) ) THEN
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
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




  END SUBROUTINE ComputePrimitive_TwoMoment_Vector_Richardson


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

            Alpha(1,iZ) = AB1 / ( AA11 + SqrtTiny )
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

  SUBROUTINE ComputeConserved_TwoMoment &
    ( D, I_u_1, I_u_2, I_u_3, N, G_d_1, G_d_2, G_d_3, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, Gm_dd_12, Gm_dd_13, Gm_dd_23, alp, B_u_1, B_u_2, B_u_3 )

    REAL(DP), INTENT(in)  :: D, I_u_1, I_u_2, I_u_3 ! --- Index Up
    REAL(DP), INTENT(out) :: N, G_d_1, G_d_2, G_d_3 ! --- Index Down
    REAL(DP), INTENT(in)  ::    V_u_1, V_u_2, V_u_3 ! --- Index Up
    REAL(DP), INTENT(in)  :: Gm_dd_11, Gm_dd_22, Gm_dd_33, Gm_dd_12, Gm_dd_13, Gm_dd_23
    REAL(DP), INTENT(in)  :: B_u_1, B_u_2, B_u_3, alp
    REAL(DP) :: k_dd_ij(1:3,1:3)
    REAL(DP) :: W  
    REAL(DP) :: DT 
    REAL(DP) :: B_d_1, B_d_2, B_d_3
    REAL(DP) :: I_d_1, I_d_2, I_d_3 ! --- Index Up

    CALL ComputeEddingtonTensorComponents_dd &
           ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
             alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3, k_dd_ij  )



    B_d_1 = Gm_dd_11 * B_u_1
    B_d_2 = Gm_dd_22 * B_u_2
    B_d_3 = Gm_dd_33 * B_u_3

    W = 1.0_DP / SQRT( 1.0_DP - (Gm_dd_11 * V_u_1 * V_u_1 &
               + Gm_dd_22 * V_u_2 * V_u_2 &  
               + Gm_dd_33 * V_u_3 * V_u_3) )
   
    DT = 1.0_DP / ( B_d_1 * V_u_1 + B_d_2 * V_u_2 + B_d_3 * V_u_3 - alp )

 
    I_d_1 = DT * ( B_d_2 * V_u_2 + B_d_3 * V_u_3 - alp ) * Gm_dd_11 * I_u_1 &
          - DT * ( B_d_1 * V_u_2 *Gm_dd_22 ) * I_u_2 - DT * ( B_d_1 * V_u_3 * Gm_dd_33 ) * I_u_3 
    I_d_2 = DT * ( B_d_1 * V_u_1 + B_d_3 * V_u_3 - alp ) * Gm_dd_22 * I_u_2 &
          - DT * ( B_d_2 * V_u_1 * Gm_dd_11 ) * I_u_1 - DT * ( Gm_dd_33 * I_u_3 * B_d_2 * V_u_3 ) 
    I_d_3 = DT * ( B_d_1 * V_u_1 + B_d_2 * V_u_2 - alp ) * Gm_dd_33 * I_u_3 &
          - DT * ( Gm_dd_11 * I_u_1 * B_d_3 * V_u_1 ) - DT * ( Gm_dd_22 * I_u_2 * B_d_3 * V_u_2 )

    ! --- Conserved Number Density ---

    N = W * D + V_u_1 * I_d_1 &
              + V_u_2 * I_d_2 &
              + V_u_3 * I_d_3

    ! --- Conserved Number Flux Density (1) ---

    G_d_1 = W * I_d_1 &
                 + (   V_u_1 * k_dd_ij(1,1) &
                     + V_u_2 * k_dd_ij(1,2) &
                     + V_u_3 * k_dd_ij(1,3) ) * D

    ! --- Conserved Number Flux Density (2) ---

    G_d_2 = W * I_d_2 &
                 + (   V_u_1 * k_dd_ij(1,2) &
                     + V_u_2 * k_dd_ij(2,2) &
                     + V_u_3 * k_dd_ij(2,3) ) * D

    ! --- Conserved Number Flux Density (3) ---

    G_d_3 = W * I_d_3 &
                 + (   V_u_1 * k_dd_ij(1,3) &
                     + V_u_2 * k_dd_ij(2,3) &
                     + V_u_3 * k_dd_ij(3,3) ) * D


  END SUBROUTINE ComputeConserved_TwoMoment


  SUBROUTINE ComputeFromConserved_TwoMoment &
     ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, CF, CR, PR )

    INTEGER,  INTENT(in)  :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)  :: &
      GX(1:nDOFX, &
         iZ_B1(2):iZ_E1(2), &
         iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4), &
         1:nGF)
    REAL(DP), INTENT(in)  :: &
      CF(1:nDOFX, &
         iZ_B1(2):iZ_E1(2), &
         iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4), &
         1:nCF)
    REAL(DP), INTENT(in)  :: &
      CR(1:nDOFZ, &
         iZ_B1(1):iZ_E1(1), &
         iZ_B1(2):iZ_E1(2), &
         iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4), &
         1:nCR,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      PR(1:nDOFZ, &
         iZ_B1(1):iZ_E1(1), &
         iZ_B1(2):iZ_E1(2), &
         iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4), &
         1:nPR,1:nSpecies)

    INTEGER  :: &
      iZ1, iZ2, iZ3, iZ4, iS, iNodeZ, iNodeX
    REAL(DP) :: &
      PF(1:nDOFX, &
         iZ_B1(2):iZ_E1(2), &
         iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4), &
         1:nPF)

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX

        CALL ComputePrimitive_Euler_Relativistic &
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
                 GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_33), &
                 Zero, Zero, Zero, &
                 GX(iNodeX,iZ2,iZ3,iZ4,iGF_Alpha ), &
                 GX(iNodeX,iZ2,iZ3,iZ4,iGF_Beta_1), &
                 GX(iNodeX,iZ2,iZ3,iZ4,iGF_Beta_2), &
                 GX(iNodeX,iZ2,iZ3,iZ4,iGF_Beta_3) )

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE ComputeFromConserved_TwoMoment


  SUBROUTINE ComputeEddingtonTensorComponents_dd &
    ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3, k_dd_ij  )

#if   defined( THORNADO_OMP_OL )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC   )
    !$ACC ROUTINE SEQ
#endif


    REAL(DP), INTENT(in)  :: &
      D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3
    REAL(DP), INTENT(out) :: &
      k_dd_ij(1:3,1:3)

    REAL(DP) :: FF, EF, a, b, DT
    REAL(DP) :: h_d_1, h_d_2, h_d_3, I_d_1, I_d_2, I_d_3
    REAL(DP) :: B_d_1, B_d_2, B_d_3, V_d_1, V_d_2, V_d_3
    REAL(DP) :: u_d_1, u_d_2, u_d_3, W
    REAL(DP) :: k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33

    FF = FluxFactor_Relativistic( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                                  alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3 )
    EF = EddingtonFactor( D, FF )
    W = 1.0_DP / SQRT( 1.0_DP - (Gm_dd_11 * V_u_1 * V_u_1 &
               + Gm_dd_22 * V_u_2 * V_u_2 &  
               + Gm_dd_33 * V_u_3 * V_u_3) )
   
    B_d_1 = Gm_dd_11 * B_u_1
    B_d_2 = Gm_dd_22 * B_u_2
    B_d_3 = Gm_dd_33 * B_u_3

    V_d_1 = Gm_dd_11 * V_u_1
    V_d_2 = Gm_dd_22 * V_u_2
    V_d_3 = Gm_dd_33 * V_u_3


    u_d_1 = W * V_d_1
    u_d_2 = W * V_d_2
    u_d_3 = W * V_d_3
 
    DT = 1.0_DP / ( B_d_1 * V_u_1 + B_d_2 * V_u_2 + B_d_3 * V_u_3 - alp )

    I_d_1 = DT * ( B_d_2 * V_u_2 + B_d_3 * V_u_3 - alp ) * Gm_dd_11 * I_u_1 &
          - DT * ( B_d_1 * V_u_2 *Gm_dd_22 ) * I_u_2 - DT * ( B_d_1 * V_u_3 * Gm_dd_33 ) * I_u_3 
    I_d_2 = DT * ( B_d_1 * V_u_1 + B_d_3 * V_u_3 - alp ) * Gm_dd_22 * I_u_2 &
          - DT * ( B_d_2 * V_u_1 * Gm_dd_11 ) * I_u_1 - DT * ( Gm_dd_33 * I_u_3 * B_d_2 * V_u_3 ) 
    I_d_3 = DT * ( B_d_1 * V_u_1 + B_d_2 * V_u_2 - alp ) * Gm_dd_33 * I_u_3 &
          - DT * ( Gm_dd_11 * I_u_1 * B_d_3 * V_u_1 ) - DT * ( Gm_dd_22 * I_u_2 * B_d_3 * V_u_2 )
    h_d_1 = I_d_1 / ( FF * D )
    h_d_2 = I_d_2 / ( FF * D )
    h_d_3 = I_d_3 / ( FF * D )

    a = Half * ( One - EF )
    b = Half * ( Three * EF - One )
   
    ! --- Diagonal Eddington Tensor Components ---

    k_dd_11 = a * ( Gm_dd_11 + u_d_1 * u_d_1 ) + b * h_d_1 * h_d_1
    k_dd_22 = a * ( Gm_dd_22 + u_d_2 * u_d_2 ) + b * h_d_2 * h_d_2
    k_dd_33 = a * ( Gm_dd_33 + u_d_3 * u_d_3 ) + b * h_d_3 * h_d_3
    ! --- Off-Diagonal Eddington Tensor Components ---
    k_dd_12 = a * u_d_1 * u_d_2 + b * h_d_1 * h_d_2
    k_dd_13 = a * u_d_1 * u_d_3 + b * h_d_1 * h_d_3
    k_dd_23 = a * u_d_2 * u_d_3 + b * h_d_2 * h_d_3
   
    k_dd_ij(1,1) = k_dd_11
 
    k_dd_ij(1,2) = k_dd_12
 
    k_dd_ij(1,3) = k_dd_13
 
    k_dd_ij(2,1) = k_dd_12
 
    k_dd_ij(2,2) = k_dd_22
 
    k_dd_ij(2,3) = k_dd_23
 
    k_dd_ij(3,1) = k_dd_13
 
    k_dd_ij(3,2) = k_dd_23
 
    k_dd_ij(3,3) = k_dd_33

  END SUBROUTINE ComputeEddingtonTensorComponents_dd

  SUBROUTINE ComputeEddingtonTensorComponents_uu &
    ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3, k_uu_munu  )


	    REAL(DP), INTENT(in)  :: &
      D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3
	    REAL(DP), INTENT(out) :: &
                k_uu_munu(0:3,0:3)


	    REAL(DP) :: FF, EF, a, b, Gm_uu_11, Gm_uu_22, Gm_uu_33
	    REAL(DP) :: h_u_1, h_u_2, h_u_3
	    REAL(DP) :: B_d_1, B_d_2, B_d_3, V_d_1, V_d_2, V_d_3, V_0
	    REAL(DP) :: u_u_1, u_u_2, u_u_3, W
            REAL(DP) :: k_uu_11, k_uu_12, k_uu_13, k_uu_22, k_uu_23, k_uu_33, k_uu_10, k_uu_20, k_uu_30, k_uu_00

    FF = FluxFactor_Relativistic( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                                  alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3 )
    EF = EddingtonFactor( D, FF )
    W = 1.0_DP / SQRT( 1.0_DP - (Gm_dd_11 * V_u_1 * V_u_1 &
               + Gm_dd_22 * V_u_2 * V_u_2 &  
               + Gm_dd_33 * V_u_3 * V_u_3) )
   

    Gm_uu_11 = 1.0_DP / Gm_dd_11 - B_u_1**2 / alp**2 
    Gm_uu_22 = 1.0_DP / Gm_dd_22 - B_u_2**2 / alp**2 
    Gm_uu_33 = 1.0_DP / Gm_dd_33 - B_u_3**2 / alp**2 

    V_d_1 = Gm_dd_11 * V_u_1
    V_d_2 = Gm_dd_22 * V_u_2
    V_d_3 = Gm_dd_33 * V_u_3

    B_d_1 = Gm_dd_11 * B_u_1
    B_d_2 = Gm_dd_22 * B_u_2
    B_d_3 = Gm_dd_33 * B_u_3

    V_0 = B_d_1 * V_u_1 + B_d_2 * V_u_2 + B_d_3 * V_u_3

    u_u_1 = W * ( V_u_1 - B_u_1 / alp ) 
    u_u_2 = W * ( V_u_2 - B_u_2 / alp ) 
    u_u_3 = W * ( V_u_2 - B_u_2 / alp ) 


    h_u_1 = I_u_1 / ( FF * D )
    h_u_2 = I_u_2 / ( FF * D )
    h_u_3 = I_u_3 / ( FF * D )


    a = Half * ( One - EF )
    b = Half * ( Three * EF - One )
   
    ! --- Diagonal Euuington Tensor Components ---

    k_uu_11 = a * ( Gm_uu_11 + u_u_1 * u_u_1 ) + b * h_u_1 * h_u_1
    k_uu_22 = a * ( Gm_uu_22 + u_u_2 * u_u_2 ) + b * h_u_2 * h_u_2
    k_uu_33 = a * ( Gm_uu_33 + u_u_3 * u_u_3 ) + b * h_u_3 * h_u_3
    ! --- Off-Diagonal Euuington Tensor Components ---
    k_uu_12 = a * u_u_1 * u_u_2 + b * h_u_1 * h_u_2
    k_uu_13 = a * u_u_1 * u_u_3 + b * h_u_1 * h_u_3
    k_uu_23 = a * u_u_2 * u_u_3 + b * h_u_2 * h_u_3

    k_uu_10 = ( 1.0_DP/( 1.0_DP - V_0 / alp ) ) * ( V_d_1 * k_uu_11 + V_d_2 * k_uu_12 + V_d_3 * k_uu_13 ) / alp

    k_uu_20 = ( 1.0_DP/( 1.0_DP - V_0 / alp ) ) * ( V_d_1 * k_uu_12 + V_d_2 * k_uu_22 + V_d_3 * k_uu_23 ) / alp

    k_uu_30 = ( 1.0_DP/( 1.0_DP - V_0 / alp ) ) * ( V_d_1 * k_uu_13 + V_d_2 * k_uu_23 + V_d_3 * k_uu_33 ) / alp
    
    k_uu_00 =  ( 1.0_DP/( 1.0_DP - V_0 / alp )**2 ) *( ( V_d_1**2 * k_uu_11 + V_d_2**2 * k_uu_22 + V_d_3**2 * k_uu_33 ) &
               + 2.0_DP * ( V_d_1 * V_d_2 * k_uu_12 + V_d_1 * V_d_3 * k_uu_13 + V_d_2 * V_d_3 * k_uu_23 ))/alp**2

    k_uu_munu(0,0) = k_uu_00
 
    k_uu_munu(0,1) = k_uu_10
 
    k_uu_munu(0,2) = k_uu_20
 
    k_uu_munu(0,3) = k_uu_30
 
    k_uu_munu(1,0) = k_uu_10
 
    k_uu_munu(1,1) = k_uu_11
 
    k_uu_munu(1,2) = k_uu_12
 
    k_uu_munu(1,3) = k_uu_13
 
    k_uu_munu(2,0) = k_uu_20
 
    k_uu_munu(2,1) = k_uu_12
 
    k_uu_munu(2,2) = k_uu_22
 
    k_uu_munu(2,3) = k_uu_23
 
    k_uu_munu(3,0) = k_uu_30
 
    k_uu_munu(3,1) = k_uu_13
 
    k_uu_munu(3,2) = k_uu_23
 
    k_uu_munu(3,3) = k_uu_33


  END SUBROUTINE ComputeEddingtonTensorComponents_uu



  SUBROUTINE ComputeEddingtonTensorComponents_ud &
    ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3, k_ud_munu )
   REAL(DP), INTENT(in)  :: &
      D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3
    REAL(DP), INTENT(out) :: &
               k_ud_munu(0:3,0:3)
    REAL(DP) :: FF, EF, a, b, DT
    REAL(DP) :: h_u_1, h_u_2, h_u_3, h_d_1, h_d_2, h_d_3, I_d_1, I_d_2, I_d_3
    REAL(DP) :: B_d_1, B_d_2, B_d_3, V_d_1, V_d_2, V_d_3, V_0
    REAL(DP) :: u_d_1, u_d_2, u_d_3, u_u_1, u_u_2, u_u_3, W, x, y
    REAL(DP) :: k_ud_11, k_ud_12, k_ud_13, k_ud_22, k_ud_23, k_ud_33, k_ud_21, k_ud_31, k_ud_32, & 
                k_ud_10, k_ud_20, k_ud_30, k_ud_01, k_ud_02, k_ud_03,  k_ud_00 

    FF = FluxFactor_Relativistic( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                                  alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3 )

    EF = EddingtonFactor( D, FF )


    W = 1.0_DP / SQRT( 1.0_DP - (Gm_dd_11 * V_u_1 * V_u_1 &
               + Gm_dd_22 * V_u_2 * V_u_2 &  
               + Gm_dd_33 * V_u_3 * V_u_3) )
   
    B_d_1 = Gm_dd_11 * B_u_1
    B_d_2 = Gm_dd_22 * B_u_2
    B_d_3 = Gm_dd_33 * B_u_3

    V_d_1 = Gm_dd_11 * V_u_1
    V_d_2 = Gm_dd_22 * V_u_2
    V_d_3 = Gm_dd_33 * V_u_3

    V_0 = B_d_1 * V_u_1 + B_d_2 * V_u_2 + B_d_3 * V_u_3

    u_d_1 = W * V_d_1
    u_d_2 = W * V_d_2
    u_d_3 = W * V_d_3

    u_u_1 = W * ( V_u_1 - B_u_1 / alp ) 
    u_u_2 = W * ( V_u_2 - B_u_2 / alp ) 
    u_u_3 = W * ( V_u_3 - B_u_3 / alp ) 
 
    DT = 1.0_DP / ( B_d_1 * V_u_1 + B_d_2 * V_u_2 + B_d_3 * V_u_3 - alp )


    I_d_1 = DT * ( B_d_2 * V_u_2 + B_d_3 * V_u_3 - alp ) * Gm_dd_11 * I_u_1 &
          - DT * ( B_d_1 * V_u_2 *Gm_dd_22 ) * I_u_2 - DT * ( B_d_1 * V_u_3 * Gm_dd_33 ) * I_u_3 
    I_d_2 = DT * ( B_d_1 * V_u_1 + B_d_3 * V_u_3 - alp ) * Gm_dd_22 * I_u_2 &
          - DT * ( B_d_2 * V_u_1 * Gm_dd_11 ) * I_u_1 - DT * ( Gm_dd_33 * I_u_3 * B_d_2 * V_u_3 ) 
    I_d_3 = DT * ( B_d_1 * V_u_1 + B_d_2 * V_u_2 - alp ) * Gm_dd_33 * I_u_3 &
          - DT * ( Gm_dd_11 * I_u_1 * B_d_3 * V_u_1 ) - DT * ( Gm_dd_22 * I_u_2 * B_d_3 * V_u_2 )

    h_u_1 = I_u_1 / ( FF * D )
    h_u_2 = I_u_2 / ( FF * D )
    h_u_3 = I_u_3 / ( FF * D )

    h_d_1 = I_d_1 / ( FF * D )
    h_d_2 = I_d_2 / ( FF * D )
    h_d_3 = I_d_3 / ( FF * D )

    a = Half * ( One - EF )
    b = Half * ( Three * EF - One )

    ! --- Diagonal Eddington Tensor Components ---

    k_ud_11 = a * ( 1.0_DP + u_u_1 * u_d_1 ) + b * h_u_1 * h_d_1
    k_ud_22 = a * ( 1.0_DP + u_u_2 * u_d_2 ) + b * h_u_2 * h_d_2
    k_ud_33 = a * ( 1.0_DP + u_u_3 * u_d_3 ) + b * h_u_3 * h_d_3

    ! --- Off-Diagonal Eddington Tensor Components ---

    k_ud_12 = a * u_u_1 * u_d_2 + b * h_u_1 * h_d_2
    k_ud_13 = a * u_u_1 * u_d_3 + b * h_u_1 * h_d_3
    k_ud_23 = a * u_u_2 * u_d_3 + b * h_u_2 * h_d_3

    k_ud_21 = a * u_u_2 * u_d_1 + b * h_u_2 * h_d_1
    k_ud_31 = a * u_u_3 * u_d_1 + b * h_u_3 * h_d_1
    k_ud_32 = a * u_u_3 * u_d_2 + b * h_u_3 * h_d_2

    x = ( 1.0_DP / alp ) * ( 1.0_DP / ( 1.0_DP - V_0 / alp ) ) 
    y = ( 1.0_DP / alp ) 
 
    k_ud_01 = x * ( V_d_1 * k_ud_11 + V_d_2 * k_ud_21 + V_d_3 * k_ud_31 )
    k_ud_02 = x * ( V_d_1 * k_ud_12 + V_d_2 * k_ud_22 + V_d_3 * k_ud_32 )
    k_ud_03 = x * ( V_d_1 * k_ud_13 + V_d_2 * k_ud_23 + V_d_3 * k_ud_33 )

    k_ud_10 = y * ( V_u_1 * k_ud_11 + V_u_2 * k_ud_12 + V_u_3 * k_ud_13 )
    k_ud_20 = y * ( V_u_1 * k_ud_21 + V_u_2 * k_ud_22 + V_u_3 * k_ud_23 )
    k_ud_30 = y * ( V_u_1 * k_ud_31 + V_u_2 * k_ud_32 + V_u_3 * k_ud_33 )

    k_ud_00 = x * y * ( V_u_1 * V_d_1 * k_ud_11 + V_u_1 * V_d_2 * k_ud_21 + V_u_1 * V_d_3 * k_ud_31 &
            + V_u_2 * V_d_1 * k_ud_12 + V_u_2 * V_d_2 * k_ud_22 + V_u_2 * V_d_3 * k_ud_32 &
            + V_u_3 * V_d_1 * k_ud_13 + V_u_3 * V_d_2 * k_ud_23 + V_u_3 * V_d_3 * k_ud_33 )
            
    k_ud_munu(0,0) = k_ud_00
 
    k_ud_munu(0,1) = k_ud_01
 
    k_ud_munu(0,2) = k_ud_02
 
    k_ud_munu(0,3) = k_ud_03
 
    k_ud_munu(1,0) = k_ud_10
 
    k_ud_munu(1,1) = k_ud_11
 
    k_ud_munu(1,2) = k_ud_12
 
    k_ud_munu(1,3) = k_ud_13
 
    k_ud_munu(2,0) = k_ud_20
 
    k_ud_munu(2,1) = k_ud_21
 
    k_ud_munu(2,2) = k_ud_22
 
    k_ud_munu(2,3) = k_ud_23
 
    k_ud_munu(3,0) = k_ud_30
 
    k_ud_munu(3,1) = k_ud_31
 
    k_ud_munu(3,2) = k_ud_32
 
    k_ud_munu(3,3) = k_ud_33

  END SUBROUTINE ComputeEddingtonTensorComponents_ud

  SUBROUTINE ComputeHeatFluxTensorComponents_ddd_Lagrangian &
    ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3, l_ddd_ijk )
 
    REAL(DP), INTENT(in)  :: &
      D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3

    REAL(DP), INTENT(out) :: &
      l_ddd_ijk(1:3,1:3,1:3)

    REAL(DP) :: FF, HF, a, b, DT
    REAL(DP) :: h_d_1, h_d_2, h_d_3
    REAL(DP) :: I_d_1, I_d_2, I_d_3
    REAL(DP) :: B_d_1, B_d_2, B_d_3, V_d_1, V_d_2, V_d_3
    REAL(DP) :: u_d_1, u_d_2, u_d_3, W
    REAL(DP) :: h_dd_11, h_dd_22, h_dd_33, h_dd_12, h_dd_13, h_dd_23, h_dd_21, h_dd_31, h_dd_32
    REAL(DP) :: &
      l_ddd_111, l_ddd_112, l_ddd_113, l_ddd_121, l_ddd_122, & 
      l_ddd_123, l_ddd_131, l_ddd_132, l_ddd_133, l_ddd_211, & 
      l_ddd_212, l_ddd_213, l_ddd_221, l_ddd_222, l_ddd_223, & 
      l_ddd_231, l_ddd_232, l_ddd_233, l_ddd_311, l_ddd_312, & 
      l_ddd_313, l_ddd_321, l_ddd_322, l_ddd_323, l_ddd_331, & 
      l_ddd_332, l_ddd_333
   

    FF = FluxFactor_Relativistic( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                                  alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3 )
    HF = HeatFluxFactor( D, FF )


    W = 1.0_DP / SQRT( 1.0_DP - (Gm_dd_11 * V_u_1 * V_u_1 &
               + Gm_dd_22 * V_u_2 * V_u_2 &  
               + Gm_dd_33 * V_u_3 * V_u_3) )
   
    B_d_1 = Gm_dd_11 * B_u_1
    B_d_2 = Gm_dd_22 * B_u_2
    B_d_3 = Gm_dd_33 * B_u_3

    V_d_1 = Gm_dd_11 * V_u_1
    V_d_2 = Gm_dd_22 * V_u_2
    V_d_3 = Gm_dd_33 * V_u_3

    u_d_1 = W * V_d_1
    u_d_2 = W * V_d_2
    u_d_3 = W * V_d_3


    DT = 1.0_DP / ( B_d_1 * V_u_1 + B_d_2 * V_u_2 + B_d_3 * V_u_3 - alp )


    I_d_1 = DT * ( B_d_2 * V_u_2 + B_d_3 * V_u_3 - alp ) * Gm_dd_11 * I_u_1 &
          - DT * ( B_d_1 * V_u_2 *Gm_dd_22 ) * I_u_2 - DT * ( B_d_1 * V_u_3 * Gm_dd_33 ) * I_u_3 
    I_d_2 = DT * ( B_d_1 * V_u_1 + B_d_3 * V_u_3 - alp ) * Gm_dd_22 * I_u_2 &
          - DT * ( B_d_2 * V_u_1 * Gm_dd_11 ) * I_u_1 - DT * ( Gm_dd_33 * I_u_3 * B_d_2 * V_u_3 ) 
    I_d_3 = DT * ( B_d_1 * V_u_1 + B_d_2 * V_u_2 - alp ) * Gm_dd_33 * I_u_3 &
          - DT * ( Gm_dd_11 * I_u_1 * B_d_3 * V_u_1 ) - DT * ( Gm_dd_22 * I_u_2 * B_d_3 * V_u_2 )


    h_dd_11 = Gm_dd_11 + u_d_1 * u_d_1
    h_dd_22 = Gm_dd_22 + u_d_2 * u_d_2
    h_dd_33 = Gm_dd_33 + u_d_3 * u_d_3
    h_dd_12 = u_d_1 * u_d_2
    h_dd_13 = u_d_1 * u_d_3
    h_dd_23 = u_d_2 * u_d_3
    h_dd_21 = u_d_1 * u_d_2
    h_dd_31 = u_d_1 * u_d_3
    h_dd_32 = u_d_2 * u_d_3
 
    a = Half * ( FF - HF )
    b = Half * ( Five * HF - Three * FF )

    h_d_1 = I_d_1 / ( FF * D )
    h_d_2 = I_d_2 / ( FF * D )
    h_d_3 = I_d_3 / ( FF * D )


     l_ddd_111 & 
      = a * ( h_d_1 * h_dd_11 + h_d_1 * h_dd_11 + h_d_1 * h_dd_11 ) + b * h_d_1 * h_d_1 * h_d_1

    l_ddd_112 & 
      = a * ( h_d_1 * h_dd_12 + h_d_1 * h_dd_12 + h_d_2 * h_dd_11 ) + b * h_d_1 * h_d_1 * h_d_2

    l_ddd_113 & 
      = a * ( h_d_1 * h_dd_13 + h_d_1 * h_dd_13 + h_d_3 * h_dd_11 ) + b * h_d_1 * h_d_1 * h_d_3

    l_ddd_121 & 
      = a * ( h_d_1 * h_dd_21 + h_d_2 * h_dd_11 + h_d_1 * h_dd_12 ) + b * h_d_1 * h_d_2 * h_d_1

    l_ddd_122 & 
      = a * ( h_d_1 * h_dd_22 + h_d_2 * h_dd_12 + h_d_2 * h_dd_12 ) + b * h_d_1 * h_d_2 * h_d_2

    l_ddd_123 & 
      = a * ( h_d_1 * h_dd_23 + h_d_2 * h_dd_13 + h_d_3 * h_dd_12 ) + b * h_d_1 * h_d_2 * h_d_3

    l_ddd_131 & 
      = a * ( h_d_1 * h_dd_31 + h_d_3 * h_dd_11 + h_d_1 * h_dd_13 ) + b * h_d_1 * h_d_3 * h_d_1

    l_ddd_132 & 
      = a * ( h_d_1 * h_dd_32 + h_d_3 * h_dd_12 + h_d_2 * h_dd_13 ) + b * h_d_1 * h_d_3 * h_d_2

    l_ddd_133 & 
      = a * ( h_d_1 * h_dd_33 + h_d_3 * h_dd_13 + h_d_3 * h_dd_13 ) + b * h_d_1 * h_d_3 * h_d_3

    l_ddd_211 & 
      = a * ( h_d_2 * h_dd_11 + h_d_1 * h_dd_21 + h_d_1 * h_dd_21 ) + b * h_d_2 * h_d_1 * h_d_1

    l_ddd_212 & 
      = a * ( h_d_2 * h_dd_12 + h_d_1 * h_dd_22 + h_d_2 * h_dd_21 ) + b * h_d_2 * h_d_1 * h_d_2

    l_ddd_213 & 
      = a * ( h_d_2 * h_dd_13 + h_d_1 * h_dd_23 + h_d_3 * h_dd_21 ) + b * h_d_2 * h_d_1 * h_d_3

    l_ddd_221 & 
      = a * ( h_d_2 * h_dd_21 + h_d_2 * h_dd_21 + h_d_1 * h_dd_22 ) + b * h_d_2 * h_d_2 * h_d_1

    l_ddd_222 & 
      = a * ( h_d_2 * h_dd_22 + h_d_2 * h_dd_22 + h_d_2 * h_dd_22 ) + b * h_d_2 * h_d_2 * h_d_2

    l_ddd_223 & 
      = a * ( h_d_2 * h_dd_23 + h_d_2 * h_dd_23 + h_d_3 * h_dd_22 ) + b * h_d_2 * h_d_2 * h_d_3

    l_ddd_231 & 
      = a * ( h_d_2 * h_dd_31 + h_d_3 * h_dd_21 + h_d_1 * h_dd_23 ) + b * h_d_2 * h_d_3 * h_d_1

    l_ddd_232 & 
      = a * ( h_d_2 * h_dd_32 + h_d_3 * h_dd_22 + h_d_2 * h_dd_23 ) + b * h_d_2 * h_d_3 * h_d_2

    l_ddd_233 & 
      = a * ( h_d_2 * h_dd_33 + h_d_3 * h_dd_23 + h_d_3 * h_dd_23 ) + b * h_d_2 * h_d_3 * h_d_3

    l_ddd_311 & 
      = a * ( h_d_3 * h_dd_11 + h_d_1 * h_dd_31 + h_d_1 * h_dd_31 ) + b * h_d_3 * h_d_1 * h_d_1

    l_ddd_312 & 
      = a * ( h_d_3 * h_dd_12 + h_d_1 * h_dd_32 + h_d_2 * h_dd_31 ) + b * h_d_3 * h_d_1 * h_d_2

    l_ddd_313 & 
      = a * ( h_d_3 * h_dd_13 + h_d_1 * h_dd_33 + h_d_3 * h_dd_31 ) + b * h_d_3 * h_d_1 * h_d_3

    l_ddd_321 & 
      = a * ( h_d_3 * h_dd_21 + h_d_2 * h_dd_31 + h_d_1 * h_dd_32 ) + b * h_d_3 * h_d_2 * h_d_1

    l_ddd_322 & 
      = a * ( h_d_3 * h_dd_22 + h_d_2 * h_dd_32 + h_d_2 * h_dd_32 ) + b * h_d_3 * h_d_2 * h_d_2

    l_ddd_323 & 
      = a * ( h_d_3 * h_dd_23 + h_d_2 * h_dd_33 + h_d_3 * h_dd_32 ) + b * h_d_3 * h_d_2 * h_d_3

    l_ddd_331 & 
      = a * ( h_d_3 * h_dd_31 + h_d_3 * h_dd_31 + h_d_1 * h_dd_33 ) + b * h_d_3 * h_d_3 * h_d_1

    l_ddd_332 & 
      = a * ( h_d_3 * h_dd_32 + h_d_3 * h_dd_32 + h_d_2 * h_dd_33 ) + b * h_d_3 * h_d_3 * h_d_2

    l_ddd_333 & 
      = a * ( h_d_3 * h_dd_33 + h_d_3 * h_dd_33 + h_d_3 * h_dd_33 ) + b * h_d_3 * h_d_3 * h_d_3

    l_ddd_ijk(1,1,1) = l_ddd_111
 
    l_ddd_ijk(1,1,2) = l_ddd_112
 
    l_ddd_ijk(1,1,3) = l_ddd_113
 
    l_ddd_ijk(1,2,1) = l_ddd_121
 
    l_ddd_ijk(1,2,2) = l_ddd_122
 
    l_ddd_ijk(1,2,3) = l_ddd_123
 
    l_ddd_ijk(1,3,1) = l_ddd_131
 
    l_ddd_ijk(1,3,2) = l_ddd_132
 
    l_ddd_ijk(1,3,3) = l_ddd_133
 
    l_ddd_ijk(2,1,1) = l_ddd_211
 
    l_ddd_ijk(2,1,2) = l_ddd_212
 
    l_ddd_ijk(2,1,3) = l_ddd_213
 
    l_ddd_ijk(2,2,1) = l_ddd_221
 
    l_ddd_ijk(2,2,2) = l_ddd_222
 
    l_ddd_ijk(2,2,3) = l_ddd_223
 
    l_ddd_ijk(2,3,1) = l_ddd_231
 
    l_ddd_ijk(2,3,2) = l_ddd_232
 
    l_ddd_ijk(2,3,3) = l_ddd_233
 
    l_ddd_ijk(3,1,1) = l_ddd_311
 
    l_ddd_ijk(3,1,2) = l_ddd_312
 
    l_ddd_ijk(3,1,3) = l_ddd_313
 
    l_ddd_ijk(3,2,1) = l_ddd_321
 
    l_ddd_ijk(3,2,2) = l_ddd_322
 
    l_ddd_ijk(3,2,3) = l_ddd_323
 
    l_ddd_ijk(3,3,1) = l_ddd_331
 
    l_ddd_ijk(3,3,2) = l_ddd_332
 
    l_ddd_ijk(3,3,3) = l_ddd_333

  END SUBROUTINE ComputeHeatFluxTensorComponents_ddd_Lagrangian

  SUBROUTINE ComputeHeatFluxTensorComponents_uud_Lagrangian &
    ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3, l_uud_munurho )
    REAL(DP), INTENT(in)  :: &
      D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3
     REAL(DP), INTENT(out)  :: &  
      l_uud_munurho(0:3,0:3,0:3)

    REAL(DP) :: FF, HF, a, b, DT, V_0, x, y
    REAL(DP) :: h_u_1, h_u_2, h_u_3
    REAL(DP) :: h_d_1, h_d_2, h_d_3
    REAL(DP) :: I_d_1, I_d_2, I_d_3
    REAL(DP) :: B_d_1, B_d_2, B_d_3
    REAL(DP) :: V_d_1, V_d_2, V_d_3
    REAL(DP) :: u_d_1, u_d_2, u_d_3, W, u_u_1, u_u_2, u_u_3
    REAL(DP) :: h_uu_11, h_uu_22, h_uu_33, h_uu_12, h_uu_13, h_uu_23, h_uu_21, h_uu_31, h_uu_32
    REAL(DP) :: h_ud_11, h_ud_22, h_ud_33, h_ud_12, h_ud_21, h_ud_13, h_ud_31, h_ud_23, h_ud_32
    REAL(DP) :: Gm_uu_11, Gm_uu_22, Gm_uu_33  
 
    REAL(DP) :: &
      l_uud_000, l_uud_001, l_uud_002, l_uud_003, & 
      l_uud_010, l_uud_011, l_uud_012, l_uud_013, & 
      l_uud_020, l_uud_021, l_uud_022, l_uud_023, & 
      l_uud_030, l_uud_031, l_uud_032, l_uud_033, & 
      l_uud_100, l_uud_101, l_uud_102, l_uud_103, & 
      l_uud_110, l_uud_111, l_uud_112, l_uud_113, & 
      l_uud_120, l_uud_121, l_uud_122, l_uud_123, & 
      l_uud_130, l_uud_131, l_uud_132, l_uud_133, & 
      l_uud_200, l_uud_201, l_uud_202, l_uud_203, & 
      l_uud_210, l_uud_211, l_uud_212, l_uud_213, & 
      l_uud_220, l_uud_221, l_uud_222, l_uud_223, & 
      l_uud_230, l_uud_231, l_uud_232, l_uud_233, & 
      l_uud_300, l_uud_301, l_uud_302, l_uud_303, & 
      l_uud_310, l_uud_311, l_uud_312, l_uud_313, & 
      l_uud_320, l_uud_321, l_uud_322, l_uud_323, & 
      l_uud_330, l_uud_331, l_uud_332, l_uud_333

    FF = FluxFactor_Relativistic( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                                  alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3 )
    HF = HeatFluxFactor( D, FF )


    W = 1.0_DP / SQRT( 1.0_DP - (Gm_dd_11 * V_u_1 * V_u_1 &
               + Gm_dd_22 * V_u_2 * V_u_2 &  
               + Gm_dd_33 * V_u_3 * V_u_3) )

    Gm_uu_11 = 1.0_DP /  Gm_dd_11 - B_u_1**2 / alp**2 
    Gm_uu_22 = 1.0_DP /  Gm_dd_22 - B_u_2**2 / alp**2 
    Gm_uu_33 = 1.0_DP /  Gm_dd_33 - B_u_3**2 / alp**2 

    B_d_1 = Gm_dd_11 * B_u_1
    B_d_2 = Gm_dd_22 * B_u_2
    B_d_3 = Gm_dd_33 * B_u_3

    V_d_1 = Gm_dd_11 * V_u_1
    V_d_2 = Gm_dd_22 * V_u_2
    V_d_3 = Gm_dd_33 * V_u_3

    u_d_1 = W * V_d_1
    u_d_2 = W * V_d_2
    u_d_3 = W * V_d_3

    u_u_1 = W * ( V_u_1 - B_u_1 / alp ) 
    u_u_2 = W * ( V_u_2 - B_u_2 / alp ) 
    u_u_3 = W * ( V_u_3 - B_u_3 / alp ) 

    DT = 1.0_DP / ( B_d_1 * V_u_1 + B_d_2 * V_u_2 + B_d_3 * V_u_3 - alp )


    I_d_1 = DT * ( B_d_2 * V_u_2 + B_d_3 * V_u_3 - alp ) * Gm_dd_11 * I_u_1 &
          - DT * ( B_d_1 * V_u_2 *Gm_dd_22 ) * I_u_2 - DT * ( B_d_1 * V_u_3 * Gm_dd_33 ) * I_u_3 
    I_d_2 = DT * ( B_d_1 * V_u_1 + B_d_3 * V_u_3 - alp ) * Gm_dd_22 * I_u_2 &
          - DT * ( B_d_2 * V_u_1 * Gm_dd_11 ) * I_u_1 - DT * ( Gm_dd_33 * I_u_3 * B_d_2 * V_u_3 ) 
    I_d_3 = DT * ( B_d_1 * V_u_1 + B_d_2 * V_u_2 - alp ) * Gm_dd_33 * I_u_3 &
          - DT * ( Gm_dd_11 * I_u_1 * B_d_3 * V_u_1 ) - DT * ( Gm_dd_22 * I_u_2 * B_d_3 * V_u_2 )

    h_uu_11 = Gm_uu_11 + u_u_1 * u_u_1
    h_uu_22 = Gm_uu_22 + u_u_2 * u_u_2
    h_uu_33 = Gm_uu_33 + u_u_3 * u_u_3
    h_uu_12 = u_u_1 * u_u_2
    h_uu_13 = u_u_1 * u_u_3
    h_uu_23 = u_u_2 * u_u_3
    h_uu_21 = u_u_1 * u_u_2
    h_uu_31 = u_u_1 * u_u_3
    h_uu_32 = u_u_2 * u_u_3


    h_ud_11 = 1.0_DP + u_u_1 * u_d_1
    h_ud_22 = 1.0_DP + u_u_2 * u_d_2
    h_ud_33 = 1.0_DP + u_u_3 * u_d_3
    h_ud_12 = u_u_1 * u_d_2
    h_ud_13 = u_u_1 * u_d_3
    h_ud_23 = u_u_2 * u_d_3
    h_ud_21 = u_u_2 * u_d_1
    h_ud_31 = u_u_3 * u_d_1
    h_ud_32 = u_u_3 * u_d_2
 




    a = Half * ( FF - HF )
    b = Half * ( Five * HF - Three * FF )

    h_u_1 = I_u_1 / ( FF * D )
    h_u_2 = I_u_2 / ( FF * D )
    h_u_3 = I_u_3 / ( FF * D )

    h_d_1 = I_d_1 / ( FF * D )
    h_d_2 = I_d_2 / ( FF * D )
    h_d_3 = I_d_3 / ( FF * D )

    ! --- Diagonal Heat Flux Tensor Components ---
    l_uud_111 & 
      = a * ( h_u_1 * h_ud_11 + h_u_1 * h_ud_11 + h_d_1 * h_uu_11 ) + b * h_u_1 * h_u_1 * h_d_1

    l_uud_112 & 
      = a * ( h_u_1 * h_ud_12 + h_u_1 * h_ud_12 + h_d_2 * h_uu_11 ) + b * h_u_1 * h_u_1 * h_d_2

    l_uud_113 & 
      = a * ( h_u_1 * h_ud_13 + h_u_1 * h_ud_13 + h_d_3 * h_uu_11 ) + b * h_u_1 * h_u_1 * h_d_3

    l_uud_121 & 
      = a * ( h_u_1 * h_ud_21 + h_u_2 * h_ud_11 + h_d_1 * h_uu_12 ) + b * h_u_1 * h_u_2 * h_d_1

    l_uud_122 & 
      = a * ( h_u_1 * h_ud_22 + h_u_2 * h_ud_12 + h_d_2 * h_uu_12 ) + b * h_u_1 * h_u_2 * h_d_2

    l_uud_123 & 
      = a * ( h_u_1 * h_ud_23 + h_u_2 * h_ud_13 + h_d_3 * h_uu_12 ) + b * h_u_1 * h_u_2 * h_d_3

    l_uud_131 & 
      = a * ( h_u_1 * h_ud_31 + h_u_3 * h_ud_11 + h_d_1 * h_uu_13 ) + b * h_u_1 * h_u_3 * h_d_1

    l_uud_132 & 
      = a * ( h_u_1 * h_ud_32 + h_u_3 * h_ud_12 + h_d_2 * h_uu_13 ) + b * h_u_1 * h_u_3 * h_d_2

    l_uud_133 & 
      = a * ( h_u_1 * h_ud_33 + h_u_3 * h_ud_13 + h_d_3 * h_uu_13 ) + b * h_u_1 * h_u_3 * h_d_3

    l_uud_211 & 
      = a * ( h_u_2 * h_ud_11 + h_u_1 * h_ud_21 + h_d_1 * h_uu_21 ) + b * h_u_2 * h_u_1 * h_d_1

    l_uud_212 & 
      = a * ( h_u_2 * h_ud_12 + h_u_1 * h_ud_22 + h_d_2 * h_uu_21 ) + b * h_u_2 * h_u_1 * h_d_2

    l_uud_213 & 
      = a * ( h_u_2 * h_ud_13 + h_u_1 * h_ud_23 + h_d_3 * h_uu_21 ) + b * h_u_2 * h_u_1 * h_d_3

    l_uud_221 & 
      = a * ( h_u_2 * h_ud_21 + h_u_2 * h_ud_21 + h_d_1 * h_uu_22 ) + b * h_u_2 * h_u_2 * h_d_1

    l_uud_222 & 
      = a * ( h_u_2 * h_ud_22 + h_u_2 * h_ud_22 + h_d_2 * h_uu_22 ) + b * h_u_2 * h_u_2 * h_d_2

    l_uud_223 & 
      = a * ( h_u_2 * h_ud_23 + h_u_2 * h_ud_23 + h_d_3 * h_uu_22 ) + b * h_u_2 * h_u_2 * h_d_3

    l_uud_231 & 
      = a * ( h_u_2 * h_ud_31 + h_u_3 * h_ud_21 + h_d_1 * h_uu_23 ) + b * h_u_2 * h_u_3 * h_d_1

    l_uud_232 & 
      = a * ( h_u_2 * h_ud_32 + h_u_3 * h_ud_22 + h_d_2 * h_uu_23 ) + b * h_u_2 * h_u_3 * h_d_2

    l_uud_233 & 
      = a * ( h_u_2 * h_ud_33 + h_u_3 * h_ud_23 + h_d_3 * h_uu_23 ) + b * h_u_2 * h_u_3 * h_d_3

    l_uud_311 & 
      = a * ( h_u_3 * h_ud_11 + h_u_1 * h_ud_31 + h_d_1 * h_uu_31 ) + b * h_u_3 * h_u_1 * h_d_1

    l_uud_312 & 
      = a * ( h_u_3 * h_ud_12 + h_u_1 * h_ud_32 + h_d_2 * h_uu_31 ) + b * h_u_3 * h_u_1 * h_d_2

    l_uud_313 & 
      = a * ( h_u_3 * h_ud_13 + h_u_1 * h_ud_33 + h_d_3 * h_uu_31 ) + b * h_u_3 * h_u_1 * h_d_3

    l_uud_321 & 
      = a * ( h_u_3 * h_ud_21 + h_u_2 * h_ud_31 + h_d_1 * h_uu_32 ) + b * h_u_3 * h_u_2 * h_d_1

    l_uud_322 & 
      = a * ( h_u_3 * h_ud_22 + h_u_2 * h_ud_32 + h_d_2 * h_uu_32 ) + b * h_u_3 * h_u_2 * h_d_2

    l_uud_323 & 
      = a * ( h_u_3 * h_ud_23 + h_u_2 * h_ud_33 + h_d_3 * h_uu_32 ) + b * h_u_3 * h_u_2 * h_d_3

    l_uud_331 & 
      = a * ( h_u_3 * h_ud_31 + h_u_3 * h_ud_31 + h_d_1 * h_uu_33 ) + b * h_u_3 * h_u_3 * h_d_1

    l_uud_332 & 
      = a * ( h_u_3 * h_ud_32 + h_u_3 * h_ud_32 + h_d_2 * h_uu_33 ) + b * h_u_3 * h_u_3 * h_d_2

    l_uud_333 & 
      = a * ( h_u_3 * h_ud_33 + h_u_3 * h_ud_33 + h_d_3 * h_uu_33 ) + b * h_u_3 * h_u_3 * h_d_3

    V_0 = B_d_1 * V_u_1 + B_d_2 * V_u_2 + B_d_3 * V_u_3

    x = ( 1.0_DP / alp ) * ( 1.0_DP / ( 1.0_DP - V_0 / alp ) ) 
    y = ( 1.0_DP / alp ) 

    l_uud_011 = x * (V_d_1 * l_uud_111 + V_d_2 * l_uud_211 + V_d_3 * l_uud_311)
 
    l_uud_012 = x * (V_d_1 * l_uud_112 + V_d_2 * l_uud_212 + V_d_3 * l_uud_312) 
 
    l_uud_013 = x * (V_d_1 * l_uud_113 + V_d_2 * l_uud_213 + V_d_3 * l_uud_313) 
 
    l_uud_021 = x * (V_d_1 * l_uud_121 + V_d_2 * l_uud_221 + V_d_3 * l_uud_321) 
 
    l_uud_022 = x * (V_d_1 * l_uud_122 + V_d_2 * l_uud_222 + V_d_3 * l_uud_322) 
 
    l_uud_023 = x * (V_d_1 * l_uud_123 + V_d_2 * l_uud_223 + V_d_3 * l_uud_323) 
 
    l_uud_031 = x * (V_d_1 * l_uud_131 + V_d_2 * l_uud_231 + V_d_3 * l_uud_331) 
 
    l_uud_032 = x * (V_d_1 * l_uud_132 + V_d_2 * l_uud_232 + V_d_3 * l_uud_332) 

    l_uud_033 = x * (V_d_1 * l_uud_133 + V_d_2 * l_uud_233 + V_d_3 * l_uud_333)
 
    l_uud_101 = x * (V_d_1 * l_uud_111 + V_d_2 * l_uud_121 + V_d_3 * l_uud_131) 
 
    l_uud_102 = x * (V_d_1 * l_uud_112 + V_d_2 * l_uud_122 + V_d_3 * l_uud_132) 
 
    l_uud_103 = x * (V_d_1 * l_uud_113 + V_d_2 * l_uud_123 + V_d_3 * l_uud_133) 
 
    l_uud_201 = x * (V_d_1 * l_uud_211 + V_d_2 * l_uud_221 + V_d_3 * l_uud_231) 
 
    l_uud_202 = x * (V_d_1 * l_uud_212 + V_d_2 * l_uud_222 + V_d_3 * l_uud_232) 
 
    l_uud_203 = x * (V_d_1 * l_uud_213 + V_d_2 * l_uud_223 + V_d_3 * l_uud_233) 
 
    l_uud_301 = x * (V_d_1 * l_uud_311 + V_d_2 * l_uud_321 + V_d_3 * l_uud_331) 
 
    l_uud_302 = x * (V_d_1 * l_uud_312 + V_d_2 * l_uud_322 + V_d_3 * l_uud_332) 
 
    l_uud_303 = x * (V_d_1 * l_uud_313 + V_d_2 * l_uud_323 + V_d_3 * l_uud_333)
 
    l_uud_110 = y * (V_u_1 * l_uud_111 + V_u_2 * l_uud_112 + V_u_3 * l_uud_113) 
 
    l_uud_120 = y * (V_u_1 * l_uud_121 + V_u_2 * l_uud_122 + V_u_3 * l_uud_123) 
 
    l_uud_130 = y * (V_u_1 * l_uud_131 + V_u_2 * l_uud_132 + V_u_3 * l_uud_133)
 
    l_uud_210 = y * (V_u_1 * l_uud_211 + V_u_2 * l_uud_212 + V_u_3 * l_uud_213)
 
    l_uud_220 = y * (V_u_1 * l_uud_221 + V_u_2 * l_uud_222 + V_u_3 * l_uud_223)
 
    l_uud_230 = y * (V_u_1 * l_uud_231 + V_u_2 * l_uud_232 + V_u_3 * l_uud_233)
 
    l_uud_310 = y * (V_u_1 * l_uud_311 + V_u_2 * l_uud_312 + V_u_3 * l_uud_313)
 
    l_uud_320 = y * (V_u_1 * l_uud_321 + V_u_2 * l_uud_322 + V_u_3 * l_uud_323)
 
    l_uud_330 = y * (V_u_1 * l_uud_331 + V_u_2 * l_uud_332 + V_u_3 * l_uud_333)

    l_uud_001 = x**2 * (  V_d_1 * V_d_1 * l_uud_111 + V_d_1 * V_d_2 * l_uud_121 + V_d_1 * V_d_3 * l_uud_131 & 
              + V_d_2 * V_d_1 * l_uud_211 + V_d_2 * V_d_2 * l_uud_221 + V_d_2 * V_d_3 * l_uud_231 & 
              + V_d_3 * V_d_1 * l_uud_311 + V_d_3 * V_d_2 * l_uud_321 + V_d_3 * V_d_3 * l_uud_331 ) 
             
 
    l_uud_002 = x**2 * (  V_d_1 * V_d_1 * l_uud_112 + V_d_1 * V_d_2 * l_uud_122 + V_d_1 * V_d_3 * l_uud_132 & 
              + V_d_2 * V_d_1 * l_uud_212 + V_d_2 * V_d_2 * l_uud_222 + V_d_2 * V_d_3 * l_uud_232 & 
              + V_d_3 * V_d_1 * l_uud_312 + V_d_3 * V_d_2 * l_uud_322 + V_d_3 * V_d_3 * l_uud_332 ) 
             
 
    l_uud_003 = x**2 * (  V_d_1 * V_d_1 * l_uud_113 + V_d_1 * V_d_2 * l_uud_123 + V_d_1 * V_d_3 * l_uud_133 & 
              + V_d_2 * V_d_1 * l_uud_213 + V_d_2 * V_d_2 * l_uud_223 + V_d_2 * V_d_3 * l_uud_233 & 
              + V_d_3 * V_d_1 * l_uud_313 + V_d_3 * V_d_2 * l_uud_323 + V_d_3 * V_d_3 * l_uud_333 ) 

    l_uud_100 = x * y * (  V_d_1 * V_u_1 * l_uud_111 + V_d_1 * V_u_2 * l_uud_112 + V_d_1 * V_u_3 * l_uud_113 & 
              + V_d_2 * V_u_1 * l_uud_121 + V_d_2 * V_u_2 * l_uud_122 + V_d_2 * V_u_3 * l_uud_123 & 
              + V_d_3 * V_u_1 * l_uud_131 + V_d_3 * V_u_2 * l_uud_132 + V_d_3 * V_u_3 * l_uud_133 ) 
             
 
    l_uud_200 = x * y * (  V_d_1 * V_u_1 * l_uud_211 + V_d_1 * V_u_2 * l_uud_212 + V_d_1 * V_u_3 * l_uud_213 & 
              + V_d_2 * V_u_1 * l_uud_221 + V_d_2 * V_u_2 * l_uud_222 + V_d_2 * V_u_3 * l_uud_223 & 
              + V_d_3 * V_u_1 * l_uud_231 + V_d_3 * V_u_2 * l_uud_232 + V_d_3 * V_u_3 * l_uud_233 ) 
             
    l_uud_300 = x * y * (  V_d_1 * V_u_1 * l_uud_311 + V_d_1 * V_u_2 * l_uud_312 + V_d_1 * V_u_3 * l_uud_313 & 
              + V_d_2 * V_u_1 * l_uud_321 + V_d_2 * V_u_2 * l_uud_322 + V_d_2 * V_u_3 * l_uud_323 & 
              + V_d_3 * V_u_1 * l_uud_331 + V_d_3 * V_u_2 * l_uud_332 + V_d_3 * V_u_3 * l_uud_333 ) 


    l_uud_010 = x * y * (  V_d_1 * V_u_1 * l_uud_111 + V_d_1 * V_u_2 * l_uud_112 + V_d_1 * V_u_3 * l_uud_113 & 
              + V_d_2 * V_u_1 * l_uud_211 + V_d_2 * V_u_2 * l_uud_212 + V_d_2 * V_u_3 * l_uud_213 & 
              + V_d_3 * V_u_1 * l_uud_311 + V_d_3 * V_u_2 * l_uud_312 + V_d_3 * V_u_3 * l_uud_313 ) 
             
 
    l_uud_020 = x * y * (  V_d_1 * V_u_1 * l_uud_121 + V_d_1 * V_u_2 * l_uud_122 + V_d_1 * V_u_3 * l_uud_123 & 
              + V_d_2 * V_u_1 * l_uud_221 + V_d_2 * V_u_2 * l_uud_222 + V_d_2 * V_u_3 * l_uud_223 & 
              + V_d_3 * V_u_1 * l_uud_321 + V_d_3 * V_u_2 * l_uud_322 + V_d_3 * V_u_3 * l_uud_323 ) 
             
 
    l_uud_030 = x * y * (  V_d_1 * V_u_1 * l_uud_131 + V_d_1 * V_u_2 * l_uud_132 + V_d_1 * V_u_3 * l_uud_133 & 
              + V_d_2 * V_u_1 * l_uud_231 + V_d_2 * V_u_2 * l_uud_232 + V_d_2 * V_u_3 * l_uud_233 & 
              + V_d_3 * V_u_1 * l_uud_331 + V_d_3 * V_u_2 * l_uud_332 + V_d_3 * V_u_3 * l_uud_333 ) 

    l_uud_000 &
      = x**2 * y * (   V_d_1 * V_d_1 * V_u_1 * l_uud_111 &
                     + V_d_1 * V_d_1 * V_u_2 * l_uud_112 &
                     + V_d_1 * V_d_1 * V_u_3 * l_uud_113 & 
                     + V_d_1 * V_d_2 * V_u_1 * l_uud_121 &
                     + V_d_1 * V_d_2 * V_u_2 * l_uud_122 &
                     + V_d_1 * V_d_2 * V_u_3 * l_uud_123 & 
                     + V_d_1 * V_d_3 * V_u_1 * l_uud_131 &
                     + V_d_1 * V_d_3 * V_u_2 * l_uud_132 &
                     + V_d_1 * V_d_3 * V_u_3 * l_uud_133 & 
                     + V_d_2 * V_d_1 * V_u_1 * l_uud_211 &
                     + V_d_2 * V_d_1 * V_u_2 * l_uud_212 &
                     + V_d_2 * V_d_1 * V_u_3 * l_uud_213 & 
                     + V_d_2 * V_d_2 * V_u_1 * l_uud_221 &
                     + V_d_2 * V_d_2 * V_u_2 * l_uud_222 &
                     + V_d_2 * V_d_2 * V_u_3 * l_uud_223 & 
                     + V_d_2 * V_d_3 * V_u_1 * l_uud_231 &
                     + V_d_2 * V_d_3 * V_u_2 * l_uud_232 &
                     + V_d_2 * V_d_3 * V_u_3 * l_uud_233 & 
                     + V_d_3 * V_d_1 * V_u_1 * l_uud_311 &
                     + V_d_3 * V_d_1 * V_u_2 * l_uud_312 &
                     + V_d_3 * V_d_1 * V_u_3 * l_uud_313 & 
                     + V_d_3 * V_d_2 * V_u_1 * l_uud_321 &
                     + V_d_3 * V_d_2 * V_u_2 * l_uud_322 &
                     + V_d_3 * V_d_2 * V_u_3 * l_uud_323 & 
                     + V_d_3 * V_d_3 * V_u_1 * l_uud_331 &
                     + V_d_3 * V_d_3 * V_u_2 * l_uud_332 &
                     + V_d_3 * V_d_3 * V_u_3 * l_uud_333 )

    l_uud_munurho(0,0,0) = l_uud_000
 
    l_uud_munurho(0,0,1) = l_uud_001
 
    l_uud_munurho(0,0,2) = l_uud_002
 
    l_uud_munurho(0,0,3) = l_uud_003
 
    l_uud_munurho(0,1,0) = l_uud_010
 
    l_uud_munurho(0,1,1) = l_uud_011
 
    l_uud_munurho(0,1,2) = l_uud_012
 
    l_uud_munurho(0,1,3) = l_uud_013
 
    l_uud_munurho(0,2,0) = l_uud_020
 
    l_uud_munurho(0,2,1) = l_uud_021
 
    l_uud_munurho(0,2,2) = l_uud_022
 
    l_uud_munurho(0,2,3) = l_uud_023
 
    l_uud_munurho(0,3,0) = l_uud_030
 
    l_uud_munurho(0,3,1) = l_uud_031
 
    l_uud_munurho(0,3,2) = l_uud_032
 
    l_uud_munurho(0,3,3) = l_uud_033
 
    l_uud_munurho(1,0,0) = l_uud_100
 
    l_uud_munurho(1,0,1) = l_uud_101
 
    l_uud_munurho(1,0,2) = l_uud_102
 
    l_uud_munurho(1,0,3) = l_uud_103
 
    l_uud_munurho(1,1,0) = l_uud_110
 
    l_uud_munurho(1,1,1) = l_uud_111
 
    l_uud_munurho(1,1,2) = l_uud_112
 
    l_uud_munurho(1,1,3) = l_uud_113
 
    l_uud_munurho(1,2,0) = l_uud_120
 
    l_uud_munurho(1,2,1) = l_uud_121
 
    l_uud_munurho(1,2,2) = l_uud_122
 
    l_uud_munurho(1,2,3) = l_uud_123
 
    l_uud_munurho(1,3,0) = l_uud_130
 
    l_uud_munurho(1,3,1) = l_uud_131
 
    l_uud_munurho(1,3,2) = l_uud_132
 
    l_uud_munurho(1,3,3) = l_uud_133
 
    l_uud_munurho(2,0,0) = l_uud_200
 
    l_uud_munurho(2,0,1) = l_uud_201
 
    l_uud_munurho(2,0,2) = l_uud_202
 
    l_uud_munurho(2,0,3) = l_uud_203
 
    l_uud_munurho(2,1,0) = l_uud_210
 
    l_uud_munurho(2,1,1) = l_uud_211
 
    l_uud_munurho(2,1,2) = l_uud_212
 
    l_uud_munurho(2,1,3) = l_uud_213
 
    l_uud_munurho(2,2,0) = l_uud_220
 
    l_uud_munurho(2,2,1) = l_uud_221
 
    l_uud_munurho(2,2,2) = l_uud_222
 
    l_uud_munurho(2,2,3) = l_uud_223
 
    l_uud_munurho(2,3,0) = l_uud_230
 
    l_uud_munurho(2,3,1) = l_uud_231
 
    l_uud_munurho(2,3,2) = l_uud_232
 
    l_uud_munurho(2,3,3) = l_uud_233
 
    l_uud_munurho(3,0,0) = l_uud_300
 
    l_uud_munurho(3,0,1) = l_uud_301
 
    l_uud_munurho(3,0,2) = l_uud_302
 
    l_uud_munurho(3,0,3) = l_uud_303
 
    l_uud_munurho(3,1,0) = l_uud_310
 
    l_uud_munurho(3,1,1) = l_uud_311
 
    l_uud_munurho(3,1,2) = l_uud_312
 
    l_uud_munurho(3,1,3) = l_uud_313
 
    l_uud_munurho(3,2,0) = l_uud_320
 
    l_uud_munurho(3,2,1) = l_uud_321
 
    l_uud_munurho(3,2,2) = l_uud_322
 
    l_uud_munurho(3,2,3) = l_uud_323
 
    l_uud_munurho(3,3,0) = l_uud_330
 
    l_uud_munurho(3,3,1) = l_uud_331
 
    l_uud_munurho(3,3,2) = l_uud_332
 
    l_uud_munurho(3,3,3) = l_uud_333


  END SUBROUTINE ComputeHeatFluxTensorComponents_uud_Lagrangian





  SUBROUTINE ComputeHeatFluxTensorComponents_udd_Lagrangian &
    ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3, l_udd_ijk )

    REAL(DP), INTENT(in)  :: &
      D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3
    REAL(DP), INTENT(out) :: &
      l_udd_ijk(1:3,1:3,1:3)      

    REAL(DP) :: &
      l_udd_111, l_udd_112, l_udd_113, l_udd_121, l_udd_122, & 
      l_udd_123, l_udd_131, l_udd_132, l_udd_133, l_udd_211, & 
      l_udd_212, l_udd_213, l_udd_221, l_udd_222, l_udd_223, & 
      l_udd_231, l_udd_232, l_udd_233, l_udd_311, l_udd_312, & 
      l_udd_313, l_udd_321, l_udd_322, l_udd_323, l_udd_331, & 
      l_udd_332, l_udd_333 
    REAL(DP) :: FF, HF, a, b, DT
    REAL(DP) :: h_u_1, h_u_2, h_u_3
    REAL(DP) :: h_d_1, h_d_2, h_d_3
    REAL(DP) :: I_d_1, I_d_2, I_d_3
    REAL(DP) :: B_d_1, B_d_2, B_d_3, V_d_1, V_d_2, V_d_3
    REAL(DP) :: u_d_1, u_d_2, u_d_3, W, u_u_1, u_u_2, u_u_3
    REAL(DP) :: h_dd_11, h_dd_22, h_dd_33, h_dd_12, h_dd_13, h_dd_23, h_dd_21, h_dd_31, h_dd_32
    REAL(DP) :: h_ud_11, h_ud_22, h_ud_33, h_ud_12, h_ud_21, h_ud_13, h_ud_31, h_ud_23, h_ud_32
   

    FF = FluxFactor_Relativistic( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                                  alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3 )
    HF = HeatFluxFactor( D, FF )


    W = 1.0_DP / SQRT( 1.0_DP - (Gm_dd_11 * V_u_1 * V_u_1 &
               + Gm_dd_22 * V_u_2 * V_u_2 &  
               + Gm_dd_33 * V_u_3 * V_u_3) )
   
    B_d_1 = Gm_dd_11 * B_u_1
    B_d_2 = Gm_dd_22 * B_u_2
    B_d_3 = Gm_dd_33 * B_u_3

    V_d_1 = Gm_dd_11 * V_u_1
    V_d_2 = Gm_dd_22 * V_u_2
    V_d_3 = Gm_dd_33 * V_u_3

    u_d_1 = W * V_d_1
    u_d_2 = W * V_d_2
    u_d_3 = W * V_d_3

    u_u_1 = W * ( V_u_1 - B_u_1 / alp ) 
    u_u_2 = W * ( V_u_2 - B_u_2 / alp ) 
    u_u_3 = W * ( V_u_3 - B_u_3 / alp ) 

    DT = 1.0_DP / ( B_d_1 * V_u_1 + B_d_2 * V_u_2 + B_d_3 * V_u_3 - alp )


    I_d_1 = DT * ( B_d_2 * V_u_2 + B_d_3 * V_u_3 - alp ) * Gm_dd_11 * I_u_1 &
          - DT * ( B_d_1 * V_u_2 *Gm_dd_22 ) * I_u_2 - DT * ( B_d_1 * V_u_3 * Gm_dd_33 ) * I_u_3 
    I_d_2 = DT * ( B_d_1 * V_u_1 + B_d_3 * V_u_3 - alp ) * Gm_dd_22 * I_u_2 &
          - DT * ( B_d_2 * V_u_1 * Gm_dd_11 ) * I_u_1 - DT * ( Gm_dd_33 * I_u_3 * B_d_2 * V_u_3 ) 
    I_d_3 = DT * ( B_d_1 * V_u_1 + B_d_2 * V_u_2 - alp ) * Gm_dd_33 * I_u_3 &
          - DT * ( Gm_dd_11 * I_u_1 * B_d_3 * V_u_1 ) - DT * ( Gm_dd_22 * I_u_2 * B_d_3 * V_u_2 )

    h_dd_11 = Gm_dd_11 + u_d_1 * u_d_1
    h_dd_22 = Gm_dd_22 + u_d_2 * u_d_2
    h_dd_33 = Gm_dd_33 + u_d_3 * u_d_3
    h_dd_12 = u_d_1 * u_d_2
    h_dd_13 = u_d_1 * u_d_3
    h_dd_23 = u_d_2 * u_d_3
    h_dd_21 = u_d_1 * u_d_2
    h_dd_31 = u_d_1 * u_d_3
    h_dd_32 = u_d_2 * u_d_3


    h_ud_11 = 1.0_DP + u_u_1 * u_d_1
    h_ud_22 = 1.0_DP + u_u_2 * u_d_2
    h_ud_33 = 1.0_DP + u_u_3 * u_d_3
    h_ud_12 = u_u_1 * u_d_2
    h_ud_13 = u_u_1 * u_d_3
    h_ud_23 = u_u_2 * u_d_3
    h_ud_21 = u_u_2 * u_d_1
    h_ud_31 = u_u_3 * u_d_1
    h_ud_32 = u_u_3 * u_d_2
 




    a = Half * ( FF - HF )
    b = Half * ( Five * HF - Three * FF )

    h_u_1 = I_u_1 / ( FF * D )
    h_u_2 = I_u_2 / ( FF * D )
    h_u_3 = I_u_3 / ( FF * D )

    h_d_1 = I_d_1 / ( FF * D )
    h_d_2 = I_d_2 / ( FF * D )
    h_d_3 = I_d_3 / ( FF * D )

    ! --- Diagonal Heat Flux Tensor Components ---

    l_udd_111 & 
      = a * ( h_u_1 * h_dd_11 + h_d_1 * h_ud_11 + h_d_1 * h_ud_11 ) + b * h_u_1 * h_d_1 * h_d_1

    l_udd_112 & 
      = a * ( h_u_1 * h_dd_12 + h_d_1 * h_ud_12 + h_d_2 * h_ud_11 ) + b * h_u_1 * h_d_1 * h_d_2

    l_udd_113 & 
      = a * ( h_u_1 * h_dd_13 + h_d_1 * h_ud_13 + h_d_3 * h_ud_11 ) + b * h_u_1 * h_d_1 * h_d_3

    l_udd_121 & 
      = a * ( h_u_1 * h_dd_21 + h_d_2 * h_ud_11 + h_d_1 * h_ud_12 ) + b * h_u_1 * h_d_2 * h_d_1

    l_udd_122 & 
      = a * ( h_u_1 * h_dd_22 + h_d_2 * h_ud_12 + h_d_2 * h_ud_12 ) + b * h_u_1 * h_d_2 * h_d_2

    l_udd_123 & 
      = a * ( h_u_1 * h_dd_23 + h_d_2 * h_ud_13 + h_d_3 * h_ud_12 ) + b * h_u_1 * h_d_2 * h_d_3

    l_udd_131 & 
      = a * ( h_u_1 * h_dd_31 + h_d_3 * h_ud_11 + h_d_1 * h_ud_13 ) + b * h_u_1 * h_d_3 * h_d_1

    l_udd_132 & 
      = a * ( h_u_1 * h_dd_32 + h_d_3 * h_ud_12 + h_d_2 * h_ud_13 ) + b * h_u_1 * h_d_3 * h_d_2

    l_udd_133 & 
      = a * ( h_u_1 * h_dd_33 + h_d_3 * h_ud_13 + h_d_3 * h_ud_13 ) + b * h_u_1 * h_d_3 * h_d_3

    l_udd_211 & 
      = a * ( h_u_2 * h_dd_11 + h_d_1 * h_ud_21 + h_d_1 * h_ud_21 ) + b * h_u_2 * h_d_1 * h_d_1

    l_udd_212 & 
      = a * ( h_u_2 * h_dd_12 + h_d_1 * h_ud_22 + h_d_2 * h_ud_21 ) + b * h_u_2 * h_d_1 * h_d_2

    l_udd_213 & 
      = a * ( h_u_2 * h_dd_13 + h_d_1 * h_ud_23 + h_d_3 * h_ud_21 ) + b * h_u_2 * h_d_1 * h_d_3

    l_udd_221 & 
      = a * ( h_u_2 * h_dd_21 + h_d_2 * h_ud_21 + h_d_1 * h_ud_22 ) + b * h_u_2 * h_d_2 * h_d_1

    l_udd_222 & 
      = a * ( h_u_2 * h_dd_22 + h_d_2 * h_ud_22 + h_d_2 * h_ud_22 ) + b * h_u_2 * h_d_2 * h_d_2

    l_udd_223 & 
      = a * ( h_u_2 * h_dd_23 + h_d_2 * h_ud_23 + h_d_3 * h_ud_22 ) + b * h_u_2 * h_d_2 * h_d_3

    l_udd_231 & 
      = a * ( h_u_2 * h_dd_31 + h_d_3 * h_ud_21 + h_d_1 * h_ud_23 ) + b * h_u_2 * h_d_3 * h_d_1

    l_udd_232 & 
      = a * ( h_u_2 * h_dd_32 + h_d_3 * h_ud_22 + h_d_2 * h_ud_23 ) + b * h_u_2 * h_d_3 * h_d_2

    l_udd_233 & 
      = a * ( h_u_2 * h_dd_33 + h_d_3 * h_ud_23 + h_d_3 * h_ud_23 ) + b * h_u_2 * h_d_3 * h_d_3

    l_udd_311 & 
      = a * ( h_u_3 * h_dd_11 + h_d_1 * h_ud_31 + h_d_1 * h_ud_31 ) + b * h_u_3 * h_d_1 * h_d_1

    l_udd_312 & 
      = a * ( h_u_3 * h_dd_12 + h_d_1 * h_ud_32 + h_d_2 * h_ud_31 ) + b * h_u_3 * h_d_1 * h_d_2

    l_udd_313 & 
      = a * ( h_u_3 * h_dd_13 + h_d_1 * h_ud_33 + h_d_3 * h_ud_31 ) + b * h_u_3 * h_d_1 * h_d_3

    l_udd_321 & 
      = a * ( h_u_3 * h_dd_21 + h_d_2 * h_ud_31 + h_d_1 * h_ud_32 ) + b * h_u_3 * h_d_2 * h_d_1

    l_udd_322 & 
      = a * ( h_u_3 * h_dd_22 + h_d_2 * h_ud_32 + h_d_2 * h_ud_32 ) + b * h_u_3 * h_d_2 * h_d_2

    l_udd_323 & 
      = a * ( h_u_3 * h_dd_23 + h_d_2 * h_ud_33 + h_d_3 * h_ud_32 ) + b * h_u_3 * h_d_2 * h_d_3

    l_udd_331 & 
      = a * ( h_u_3 * h_dd_31 + h_d_3 * h_ud_31 + h_d_1 * h_ud_33 ) + b * h_u_3 * h_d_3 * h_d_1

    l_udd_332 & 
      = a * ( h_u_3 * h_dd_32 + h_d_3 * h_ud_32 + h_d_2 * h_ud_33 ) + b * h_u_3 * h_d_3 * h_d_2

    l_udd_333 & 
      = a * ( h_u_3 * h_dd_33 + h_d_3 * h_ud_33 + h_d_3 * h_ud_33 ) + b * h_u_3 * h_d_3 * h_d_3

    l_udd_ijk(1,1,1) = l_udd_111
 
    l_udd_ijk(1,1,2) = l_udd_112
 
    l_udd_ijk(1,1,3) = l_udd_113
 
    l_udd_ijk(1,2,1) = l_udd_121
 
    l_udd_ijk(1,2,2) = l_udd_122
 
    l_udd_ijk(1,2,3) = l_udd_123
 
    l_udd_ijk(1,3,1) = l_udd_131
 
    l_udd_ijk(1,3,2) = l_udd_132
 
    l_udd_ijk(1,3,3) = l_udd_133
 
    l_udd_ijk(2,1,1) = l_udd_211
 
    l_udd_ijk(2,1,2) = l_udd_212
 
    l_udd_ijk(2,1,3) = l_udd_213
 
    l_udd_ijk(2,2,1) = l_udd_221
 
    l_udd_ijk(2,2,2) = l_udd_222
 
    l_udd_ijk(2,2,3) = l_udd_223
 
    l_udd_ijk(2,3,1) = l_udd_231
 
    l_udd_ijk(2,3,2) = l_udd_232
 
    l_udd_ijk(2,3,3) = l_udd_233
 
    l_udd_ijk(3,1,1) = l_udd_311
 
    l_udd_ijk(3,1,2) = l_udd_312
 
    l_udd_ijk(3,1,3) = l_udd_313
 
    l_udd_ijk(3,2,1) = l_udd_321
 
    l_udd_ijk(3,2,2) = l_udd_322
 
    l_udd_ijk(3,2,3) = l_udd_323
 
    l_udd_ijk(3,3,1) = l_udd_331
 
    l_udd_ijk(3,3,2) = l_udd_332
 
    l_udd_ijk(3,3,3) = l_udd_333
 

  END SUBROUTINE ComputeHeatFluxTensorComponents_udd_Lagrangian







  FUNCTION Flux_X1( D, I_u_1, I_u_2, I_u_3, V_u_1, V_u_2, V_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                      alp, B_u_1, B_u_2, B_u_3 )

    REAL(DP)             :: Flux_X1(4)
    REAL(DP), INTENT(in) :: D, I_u_1, I_u_2, I_u_3
    REAL(DP), INTENT(in) ::    V_u_1, V_u_2, V_u_3
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), INTENT(in) :: alp, B_u_1, B_u_2, B_u_3

    REAL(DP) :: k_ud_ij(0:3,0:3)
    REAL(DP) :: W, DT, I_d_1, I_d_2, I_d_3, B_d_1, B_d_2, B_d_3


    W = 1.0_DP / SQRT( 1.0_DP - (Gm_dd_11 * V_u_1 * V_u_1 &
               + Gm_dd_22 * V_u_2 * V_u_2 &  
               + Gm_dd_33 * V_u_3 * V_u_3) )

    B_d_1 = Gm_dd_11 * B_u_1
    B_d_2 = Gm_dd_22 * B_u_2
    B_d_3 = Gm_dd_33 * B_u_3

    DT = 1.0_DP / ( B_d_1 * V_u_1 + B_d_2 * V_u_2 + B_d_3 * V_u_3 - alp )

    I_d_1 = DT * ( B_d_2 * V_u_2 + B_d_3 * V_u_3 - alp ) * Gm_dd_11 * I_u_1 &
          - DT * ( B_d_1 * V_u_2 *Gm_dd_22 ) * I_u_2 - DT * ( B_d_1 * V_u_3 * Gm_dd_33 ) * I_u_3 
    I_d_2 = DT * ( B_d_1 * V_u_1 + B_d_3 * V_u_3 - alp ) * Gm_dd_22 * I_u_2 &
          - DT * ( B_d_2 * V_u_1 * Gm_dd_11 ) * I_u_1 - DT * ( Gm_dd_33 * I_u_3 * B_d_2 * V_u_3 ) 
    I_d_3 = DT * ( B_d_1 * V_u_1 + B_d_2 * V_u_2 - alp ) * Gm_dd_33 * I_u_3 &
          - DT * ( Gm_dd_11 * I_u_1 * B_d_3 * V_u_1 ) - DT * ( Gm_dd_22 * I_u_2 * B_d_3 * V_u_2 )


    CALL ComputeEddingtonTensorComponents_ud &
    ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3, k_ud_ij )

    Flux_X1(1) = alp * I_u_1 + ( alp * V_u_1 - B_u_1 ) * W * D
    Flux_X1(2) = alp * k_ud_ij(1,1) * D + (alp * V_u_1 - B_u_1 ) * W * I_d_1
    Flux_X1(3) = alp * k_ud_ij(1,2) * D + (alp * V_u_1 - B_u_1 ) * W * I_d_2
    Flux_X1(4) = alp * k_ud_ij(1,3) * D + (alp * V_u_1 - B_u_1 ) * W * I_d_3
    RETURN

  END FUNCTION FLUX_X1

  FUNCTION Flux_X2( D, I_u_1, I_u_2, I_u_3, V_u_1, V_u_2, V_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                      alp, B_u_1, B_u_2, B_u_3 )

    REAL(DP)             :: Flux_X2(4)
    REAL(DP), INTENT(in) :: D, I_u_1, I_u_2, I_u_3
    REAL(DP), INTENT(in) ::    V_u_1, V_u_2, V_u_3
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), INTENT(in) :: alp, B_u_1, B_u_2, B_u_3

    REAL(DP) :: k_ud_ij(0:3,0:3)
    REAL(DP) :: W, DT, I_d_1, I_d_2, I_d_3, B_d_1, B_d_2, B_d_3


    W = 1.0_DP / SQRT( 1.0_DP - (Gm_dd_11 * V_u_1 * V_u_1 &
               + Gm_dd_22 * V_u_2 * V_u_2 &  
               + Gm_dd_33 * V_u_3 * V_u_3) )

    B_d_1 = Gm_dd_11 * B_u_1
    B_d_2 = Gm_dd_22 * B_u_2
    B_d_3 = Gm_dd_33 * B_u_3

    DT = 1.0_DP / ( B_d_1 * V_u_1 + B_d_2 * V_u_2 + B_d_3 * V_u_3 - alp )

    I_d_1 = DT * ( B_d_2 * V_u_2 + B_d_3 * V_u_3 - alp ) * Gm_dd_11 * I_u_1 &
          - DT * ( B_d_1 * V_u_2 *Gm_dd_22 ) * I_u_2 - DT * ( B_d_1 * V_u_3 * Gm_dd_33 ) * I_u_3 
    I_d_2 = DT * ( B_d_1 * V_u_1 + B_d_3 * V_u_3 - alp ) * Gm_dd_22 * I_u_2 &
          - DT * ( B_d_2 * V_u_1 * Gm_dd_11 ) * I_u_1 - DT * ( Gm_dd_33 * I_u_3 * B_d_2 * V_u_3 ) 
    I_d_3 = DT * ( B_d_1 * V_u_1 + B_d_2 * V_u_2 - alp ) * Gm_dd_33 * I_u_3 &
          - DT * ( Gm_dd_11 * I_u_1 * B_d_3 * V_u_1 ) - DT * ( Gm_dd_22 * I_u_2 * B_d_3 * V_u_2 )

    CALL ComputeEddingtonTensorComponents_ud &
    ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3, k_ud_ij )

    Flux_X2(1) = alp * I_u_2 + ( alp * V_u_2 - B_u_2 ) * W * D
    Flux_X2(2) = alp * k_ud_ij(1,2) * D + (alp * V_u_2 - B_u_2 ) * W * I_d_1
    Flux_X2(3) = alp * k_ud_ij(2,2) * D + (alp * V_u_2 - B_u_2 ) * W * I_d_2
    Flux_X2(4) = alp * k_ud_ij(2,3) * D + (alp * V_u_2 - B_u_2 ) * W * I_d_3




    RETURN

  END FUNCTION FLUX_X2

  FUNCTION Flux_X3( D, I_u_1, I_u_2, I_u_3, V_u_1, V_u_2, V_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                      alp, B_u_1, B_u_2, B_u_3 )

    REAL(DP)             :: Flux_X3(4)
    REAL(DP), INTENT(in) :: D, I_u_1, I_u_2, I_u_3
    REAL(DP), INTENT(in) ::    V_u_1, V_u_2, V_u_3
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), INTENT(in) :: alp, B_u_1, B_u_2, B_u_3

    REAL(DP) :: k_ud_ij(0:3,0:3)
    REAL(DP) :: W, DT, I_d_1, I_d_2, I_d_3, B_d_1, B_d_2, B_d_3


    W = 1.0_DP / SQRT( 1.0_DP - (Gm_dd_11 * V_u_1 * V_u_1 &
               + Gm_dd_22 * V_u_2 * V_u_2 &  
               + Gm_dd_33 * V_u_3 * V_u_3) )

    B_d_1 = Gm_dd_11 * B_u_1
    B_d_2 = Gm_dd_22 * B_u_2
    B_d_3 = Gm_dd_33 * B_u_3

    DT = 1.0_DP / ( B_d_1 * V_u_1 + B_d_2 * V_u_2 + B_d_3 * V_u_3 - alp )

    I_d_1 = DT * ( B_d_2 * V_u_2 + B_d_3 * V_u_3 - alp ) * Gm_dd_11 * I_u_1 &
          - DT * ( B_d_1 * V_u_2 *Gm_dd_22 ) * I_u_2 - DT * ( B_d_1 * V_u_3 * Gm_dd_33 ) * I_u_3 
    I_d_2 = DT * ( B_d_1 * V_u_1 + B_d_3 * V_u_3 - alp ) * Gm_dd_22 * I_u_2 &
          - DT * ( B_d_2 * V_u_1 * Gm_dd_11 ) * I_u_1 - DT * ( Gm_dd_33 * I_u_3 * B_d_2 * V_u_3 ) 
    I_d_3 = DT * ( B_d_1 * V_u_1 + B_d_2 * V_u_2 - alp ) * Gm_dd_33 * I_u_3 &
          - DT * ( Gm_dd_11 * I_u_1 * B_d_3 * V_u_1 ) - DT * ( Gm_dd_22 * I_u_2 * B_d_3 * V_u_2 )

    CALL ComputeEddingtonTensorComponents_ud &
    ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3, k_ud_ij )


    Flux_X3(1) = alp * I_u_3 + ( alp * V_u_3 - B_u_3 ) * W * D
    Flux_X3(2) = alp * k_ud_ij(1,3) * D + (alp * V_u_3 - B_u_3 ) * W * I_d_1
    Flux_X3(3) = alp * k_ud_ij(2,3) * D + (alp * V_u_3 - B_u_3 ) * W * I_d_2
    Flux_X3(4) = alp * k_ud_ij(3,3) * D + (alp * V_u_3 - B_u_3 ) * W * I_d_3




    RETURN

  END FUNCTION FLUX_X3

  FUNCTION Flux_E( D, I_u_1, I_u_2, I_u_3, V_u_1, V_u_2, V_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                      alp, B_u_1, B_u_2, B_u_3, U, dU_dX0, dU_dX1, dU_dX2, dU_dX3)

! negative sign has been incorperated into Flux_E

    REAL(DP)             :: Flux_E(4)
    REAL(DP), INTENT(in) :: D, I_u_1, I_u_2, I_u_3
    REAL(DP), INTENT(in) ::    V_u_1, V_u_2, V_u_3
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), INTENT(in) :: alp, B_u_1, B_u_2, B_u_3
    REAL(DP), INTENT(in) :: U(0:3), dU_dX0(0:3), dU_dX1(0:3), dU_dX2(0:3), dU_dX3(0:3)

    REAL(DP) :: V_0, B_d_1, B_d_2, B_d_3, V_d_1, V_d_2, V_d_3
    REAL(DP) :: I(0:3), k_uu_munu(0:3,0:3), k_ud_munu(0:3,0:3), l_uud_munurho(0:3,0:3,0:3), dU_dX(0:3,0:3)
    INTEGER :: mu, nu

    CALL ComputeEddingtonTensorComponents_uu &
    ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3, k_uu_munu  )


    CALL ComputeEddingtonTensorComponents_ud &
    ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3, k_ud_munu )

    CALL ComputeHeatFluxTensorComponents_uud_Lagrangian &
    ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3, l_uud_munurho )

    B_d_1 = Gm_dd_11 * B_u_1 
    B_d_2 = Gm_dd_22 * B_u_2 
    B_d_3 = Gm_dd_33 * B_u_3 

    V_d_1 = Gm_dd_11 * V_u_1 
    V_d_2 = Gm_dd_22 * V_u_2 
    V_d_3 = Gm_dd_33 * V_u_3
 
    V_0 = B_d_1 * V_u_1 + B_d_2 * V_u_2 + B_d_3 * V_u_3   

    I(1) = I_u_1
    I(2) = I_u_2
    I(3) = I_u_3
    I(0) = 1.0_DP / ( 1.0_DP - V_0 / alp ) * ( V_d_1 * I(1) + V_d_2 * I(2) + V_d_3 * I(3) ) / alp



    dU_dX(0,0) = dU_dX0(0) 
    dU_dX(0,1) = dU_dX0(1) 
    dU_dX(0,2) = dU_dX0(2) 
    dU_dX(0,3) = dU_dX0(3) 

    dU_dX(1,0) = dU_dX1(0) 
    dU_dX(1,1) = dU_dX1(1) 
    dU_dX(1,2) = dU_dX1(2) 
    dU_dX(1,3) = dU_dX1(3) 

    dU_dX(2,0) = dU_dX2(0) 
    dU_dX(2,1) = dU_dX2(1) 
    dU_dX(2,2) = dU_dX2(2) 
    dU_dX(2,3) = dU_dX2(3)
 
    dU_dX(3,0) = dU_dX3(0) 
    dU_dX(3,1) = dU_dX3(1) 
    dU_dX(3,2) = dU_dX3(2) 
    dU_dX(3,3) = dU_dX3(3) 

    Flux_E = 0.0_DP
   
    DO mu = 0,3
    DO nu = 0,3

      Flux_E(1) = Flux_E(1) + ( I(nu) * U(mu) + k_uu_munu(mu,nu) * D ) * dU_dX(mu,nu) 
      Flux_E(2) = Flux_E(2) + ( k_ud_munu(nu,1) * D * U(mu) + l_uud_munurho(mu,nu,1) * D ) * dU_dX(mu,nu)
      Flux_E(3) = Flux_E(3) + ( k_ud_munu(nu,2) * D * U(mu) + l_uud_munurho(mu,nu,2) * D ) * dU_dX(mu,nu)
      Flux_E(4) = Flux_E(4) + ( k_ud_munu(nu,3) * D * U(mu) + l_uud_munurho(mu,nu,3) * D ) * dU_dX(mu,nu)

    END DO
    END DO


    Flux_E = - Flux_E

    RETURN

  END FUNCTION Flux_E

  FUNCTION Source_E( D, I_u_1, I_u_2, I_u_3, V_u_1, V_u_2, V_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                      alp, B_u_1, B_u_2, B_u_3, U_u, U_d, dU_dX0, dU_dX1, dU_dX2, dU_dX3, &
                      dU_dX0_COV, dU_dX1_COV, dU_dX2_COV, dU_dX3_COV, &
                      dG_dd_dX1, dG_dd_dX2, dG_dd_dX3, E, R  )

! negative sign has been incorperated into Flux_E

    REAL(DP)             :: Source_E(3)
    REAL(DP), INTENT(in) :: D, I_u_1, I_u_2, I_u_3
    REAL(DP), INTENT(in) ::    V_u_1, V_u_2, V_u_3
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), INTENT(in) :: alp, B_u_1, B_u_2, B_u_3
    REAL(DP), INTENT(in) :: U_u(0:3), U_d(0:3) 
    REAL(DP), INTENT(in) :: dU_dX0(0:3), dU_dX1(0:3), dU_dX2(0:3), dU_dX3(0:3)
    REAL(DP), INTENT(in) :: dU_dX0_COV(0:3), dU_dX1_COV(0:3), dU_dX2_COV(0:3), dU_dX3_COV(0:3)
    REAL(DP), INTENT(in) :: dG_dd_dX1(0:3,0:3), dG_dd_dX2(0:3,0:3), dG_dd_dX3(0:3,0:3), R
    REAL(DP) :: dG_dd_dX(1:3,0:3,0:3), E
    REAL(DP) :: V_0, B_d_1, B_d_2, B_d_3, V_d_1, V_d_2, V_d_3
    REAL(DP) :: I(0:3), k_uu_munu(0:3,0:3), k_ud_munu(0:3,0:3), l_uud_munurho(0:3,0:3,0:3), dU_dX(0:3,0:3), dU_dX_COV(0:3,0:3)
    INTEGER :: mu, nu, j

    CALL ComputeEddingtonTensorComponents_uu &
    ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3, k_uu_munu  )


    CALL ComputeEddingtonTensorComponents_ud &
    ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3, k_ud_munu )

    CALL ComputeHeatFluxTensorComponents_uud_Lagrangian &
    ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3, l_uud_munurho )

    B_d_1 = Gm_dd_11 * B_u_1 
    B_d_2 = Gm_dd_22 * B_u_2 
    B_d_3 = Gm_dd_33 * B_u_3 

    V_d_1 = Gm_dd_11 * V_u_1 
    V_d_2 = Gm_dd_22 * V_u_2 
    V_d_3 = Gm_dd_33 * V_u_3
 
    V_0 = B_d_1 * V_u_1 + B_d_2 * V_u_2 + B_d_3 * V_u_3   

    I(1) = I_u_1
    I(2) = I_u_2
    I(3) = I_u_3
    I(0) = 1.0_DP / ( 1.0_DP - V_0 / alp ) * ( V_d_1 * I(1) + V_d_2 * I(2) + V_d_3 * I(3) ) / alp

    dU_dX(0,0) = dU_dX0(0) 
    dU_dX(0,1) = dU_dX0(1) 
    dU_dX(0,2) = dU_dX0(2) 
    dU_dX(0,3) = dU_dX0(3) 

    dU_dX(1,0) = dU_dX1(0) 
    dU_dX(1,1) = dU_dX1(1) 
    dU_dX(1,2) = dU_dX1(2) 
    dU_dX(1,3) = dU_dX1(3) 

    dU_dX(2,0) = dU_dX2(0) 
    dU_dX(2,1) = dU_dX2(1) 
    dU_dX(2,2) = dU_dX2(2) 
    dU_dX(2,3) = dU_dX2(3)
 
    dU_dX(3,0) = dU_dX3(0) 
    dU_dX(3,1) = dU_dX3(1) 
    dU_dX(3,2) = dU_dX3(2) 
    dU_dX(3,3) = dU_dX3(3) 

    dU_dX_COV(0,0) = dU_dX0_COV(0) 
    dU_dX_COV(0,1) = dU_dX0_COV(1) 
    dU_dX_COV(0,2) = dU_dX0_COV(2) 
    dU_dX_COV(0,3) = dU_dX0_COV(3) 

    dU_dX_COV(1,0) = dU_dX1_COV(0) 
    dU_dX_COV(1,1) = dU_dX1_COV(1) 
    dU_dX_COV(1,2) = dU_dX1_COV(2) 
    dU_dX_COV(1,3) = dU_dX1_COV(3) 

    dU_dX_COV(2,0) = dU_dX2_COV(0) 
    dU_dX_COV(2,1) = dU_dX2_COV(1) 
    dU_dX_COV(2,2) = dU_dX2_COV(2) 
    dU_dX_COV(2,3) = dU_dX2_COV(3)
 
    dU_dX_COV(3,0) = dU_dX3_COV(0) 
    dU_dX_COV(3,1) = dU_dX3_COV(1) 
    dU_dX_COV(3,2) = dU_dX3_COV(2) 
    dU_dX_COV(3,3) = dU_dX3_COV(3) 


    dG_dd_dX(1,:,:) = dG_dd_dX1(:,:)
    dG_dd_dX(2,:,:) = dG_dd_dX2(:,:)
    dG_dd_dX(3,:,:) = dG_dd_dX3(:,:)

    Source_E = 0.0_DP



    DO j = 1,3
    DO mu = 0,3
    DO nu = 0,3

      Source_E(j) = Source_E(j) - ( I(nu) * U_d(j) * U_u(mu) + l_uud_munurho(mu,nu,j) * D &
                    + k_ud_munu(nu,j) * D * U_u(mu) + k_uu_munu(nu,mu) * D * U_d(j) ) * dU_dX_COV(mu,nu)

    END DO
    END DO
    END DO

    DO j = 1,3
    DO mu = 0,3
    DO nu = 0,3

      Source_E(j) = Source_E(j) - 0.5_DP * ( D * U_u(mu) * U_u(nu) + I(mu) * U_u(nu) &
                                + I(nu) * U_u(mu) + k_uu_munu(mu,nu) * D ) * dG_dd_dX(j,mu,nu)
    END DO
    END DO
    END DO

    DO j = 1,3
    DO mu = 0,3

      Source_E(j) = Source_E(j) + ( D * U_u(mu) + I(mu) ) * dU_dX(mu,j)

    END DO
    END DO
    RETURN
  END FUNCTION Source_E


  FUNCTION NumericalFlux_LLF( u_L, u_R, Flux_L, Flux_R, alpha )

    REAL(DP)             :: NumericalFlux_LLF
    REAL(DP), INTENT(in) :: u_L, u_R, flux_L, flux_R, alpha

    NumericalFlux_LLF &
      = Half * ( flux_L + flux_R - alpha * ( u_R - u_L ) )

    RETURN
  END FUNCTION NumericalFlux_LLF


END MODULE TwoMoment_UtilitiesModule_Relativistic!
