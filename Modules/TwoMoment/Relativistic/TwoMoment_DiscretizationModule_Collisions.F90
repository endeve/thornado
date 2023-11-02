MODULE TwoMoment_DiscretizationModule_Collisions

  USE KindModule, ONLY: &
    DP, Zero, One
  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    nDOFE, &
    nDOFZ
  USE GeometryFieldsModuleE, ONLY: &
    nGE
  USE GeometryFieldsModule, ONLY: &
    nGF, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33, &
    iGF_Alpha, iGF_Beta_1, iGF_Beta_2, iGF_Beta_3
  USE FluidFieldsModule, ONLY: &
    nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne
  USE Euler_UtilitiesModule_Relativistic, ONLY: &
    ComputePrimitive_Euler_Relativistic
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3
  USE TwoMoment_UtilitiesModule, ONLY: &
    ComputeEddingtonTensorComponents_dd
  USE TwoMoment_OpacityModule, ONLY: &
    uOP, iOP_D0, iOP_Chi, iOP_Sigma, nOP

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeIncrement_TwoMoment_Implicit

  INTEGER :: iE_B0,    iE_E0
  INTEGER :: iX_B0(3), iX_E0(3)
  INTEGER :: nZ(4), nE, nX(3), nE_G, nX_G

  REAL(DP)              :: PF_N(nPF)
  REAL(DP), ALLOCATABLE :: GX_N(:,:)
  REAL(DP), ALLOCATABLE :: CF_N(:,:)
  REAL(DP), ALLOCATABLE :: CR_N(:,:,:,:)
  REAL(DP), ALLOCATABLE :: dCR_N(:,:,:,:)
  REAL(DP), ALLOCATABLE :: OP_N(:,:,:,:)

CONTAINS


  SUBROUTINE ComputeIncrement_TwoMoment_Implicit &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, dt, GE, GX, U_F, U_R, dU_R, Verbose_Option )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)  :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)  :: &
      dt
    REAL(DP), INTENT(in)  :: &
      GE  (1:nDOFE,iZ_B1(1):iZ_E1(1),1:nGE)
    REAL(DP), INTENT(in)  :: &
      GX  (1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
                   iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(in)  :: &
      U_F (1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
                   iZ_B1(4):iZ_E1(4),1:nCF)
    REAL(DP), INTENT(in)  :: &
      U_R (1:nDOFZ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
                   iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      dU_R(1:nDOFZ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
                   iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    LOGICAL,          INTENT(in), OPTIONAL :: Verbose_Option

    INTEGER :: iZ1, iZ2, iZ3, iZ4, iCR, iS, iGF, iCF, iOP
    INTEGER :: iX1, iX2, iX3, iE
    INTEGER :: iNodeZ, iNodeX, iNodeE, iN_X, iN_E
    LOGICAL :: Verbose

    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    IF (Verbose) THEN
      PRINT*, "      ComputeIncrement_TwoMoment_Implicit (Start)"
    END IF


    CALL InitializeCollisions( iZ_B0, iZ_E0, iZ_B1, iZ_E1 )
    DO iS  = 1, nSpecies
    DO iCR = 1, nCR
    DO iZ4 = iZ_B1(4), iZ_E1(4)
    DO iZ3 = iZ_B1(3), iZ_E1(3)
    DO iZ2 = iZ_B1(2), iZ_E1(2)
    DO iZ1 = iZ_B1(1), iZ_E1(1)

      DO iNodeZ = 1, nDOFZ

        dU_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS) = Zero

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    ! --- Arrange Geometry Fields ---

    DO iN_X = 1, nX_G
    DO iGF  = 1, nGF

      iX3    = MOD( (iN_X-1) / ( nDOFX * nX(1) * nX(2) ), nX(3) ) + iX_B0(3)
      iX2    = MOD( (iN_X-1) / ( nDOFX * nX(1)         ), nX(2) ) + iX_B0(2)
      iX1    = MOD( (iN_X-1) / ( nDOFX                 ), nX(1) ) + iX_B0(1)
      iNodeX = MOD( (iN_X-1)                            , nDOFX ) + 1

      GX_N(iGF,iN_X) = GX(iNodeX,iX1,iX2,iX3,iGF)

    END DO
    END DO
    ! --- Arrange Fluid Fields ---

    DO iN_X = 1, nX_G
    DO iCF  = 1, nCF

      iX3    = MOD( (iN_X-1) / ( nDOFX * nX(1) * nX(2) ), nX(3) ) + iX_B0(3)
      iX2    = MOD( (iN_X-1) / ( nDOFX * nX(1)         ), nX(2) ) + iX_B0(2)
      iX1    = MOD( (iN_X-1) / ( nDOFX                 ), nX(1) ) + iX_B0(1)
      iNodeX = MOD( (iN_X-1)                            , nDOFX ) + 1

      CF_N(iCF,iN_X) = U_F(iNodeX,iX1,iX2,iX3,iCF)

    END DO
    END DO

    ! --- Arrange Radiation Fields ---

    DO iS   = 1, nSpecies
    DO iCR  = 1, nCR
    DO iN_X = 1, nX_G
    DO iN_E = 1, nE_G

      iE     = MOD( (iN_E-1) / nDOFE, nZ(1) ) + iE_B0
      iNodeE = MOD( (iN_E-1)        , nDOFE ) + 1

      iX3    = MOD( (iN_X-1) / ( nDOFX * nX(1) * nX(2) ), nX(3) ) + iX_B0(3)
      iX2    = MOD( (iN_X-1) / ( nDOFX * nX(1)         ), nX(2) ) + iX_B0(2)
      iX1    = MOD( (iN_X-1) / ( nDOFX                 ), nX(1) ) + iX_B0(1)
      iNodeX = MOD( (iN_X-1)                            , nDOFX ) + 1

      iNodeZ = iNodeE + ( iNodeX - 1 ) * nDOFE

      CR_N(iCR,iS,iN_E,iN_X) = U_R(iNodeZ,iE,iX1,iX2,iX3,iCR,iS)

    END DO
    END DO
    END DO
    END DO

    ! --- Arrange Opacities ---

    DO iS   = 1, nSpecies
    DO iOP  = 1, nOP
    DO iN_X = 1, nX_G
    DO iN_E = 1, nE_G

      iE     = MOD( (iN_E-1) / nDOFE, nZ(1) ) + iE_B0
      iNodeE = MOD( (iN_E-1)        , nDOFE ) + 1

      iX3    = MOD( (iN_X-1) / ( nDOFX * nX(1) * nX(2) ), nX(3) ) + iX_B0(3)
      iX2    = MOD( (iN_X-1) / ( nDOFX * nX(1)         ), nX(2) ) + iX_B0(2)
      iX1    = MOD( (iN_X-1) / ( nDOFX                 ), nX(1) ) + iX_B0(1)
      iNodeX = MOD( (iN_X-1)                            , nDOFX ) + 1

      iNodeZ = iNodeE + ( iNodeX - 1 ) * nDOFE

      OP_N(iOP,iS,iN_E,iN_X) = uOP(iNodeZ,iE,iX1,iX2,iX3,iOP,iS)
    END DO
    END DO
    END DO
    END DO


    DO iN_X = 1, nX_G

      CALL ComputePrimitive_Euler_Relativistic &
             ( CF_N(iCF_D ,iN_X), &
               CF_N(iCF_S1,iN_X), &
               CF_N(iCF_S2,iN_X), &
               CF_N(iCF_S3,iN_X), &
               CF_N(iCF_E ,iN_X), &
               CF_N(iCF_Ne,iN_X), &
               PF_N(iPF_D ), &
               PF_N(iPF_V1), &
               PF_N(iPF_V2), &
               PF_N(iPF_V3), &
               PF_N(iPF_E ), &
               PF_N(iPF_Ne), &
               GX_N(iGF_Gm_dd_11,iN_X), &
               GX_N(iGF_Gm_dd_22,iN_X), &
               GX_N(iGF_Gm_dd_33,iN_X) )

      DO iN_E = 1, nE_G
      DO iS   = 1, nSpecies
        CALL ComputeIncrement_FixedPoint &
               ( dt, &
                 CR_N(iCR_N ,iS,iN_E,iN_X), &
                 CR_N(iCR_G1,iS,iN_E,iN_X), &
                 CR_N(iCR_G2,iS,iN_E,iN_X), &
                 CR_N(iCR_G3,iS,iN_E,iN_X), &
                 PF_N(iPF_V1), &
                 PF_N(iPF_V2), &
                 PF_N(iPF_V3), &
                 GX_N(iGF_Gm_dd_11,iN_X), &
                 GX_N(iGF_Gm_dd_22,iN_X), &
                 GX_N(iGF_Gm_dd_33,iN_X), &
                 OP_N(iOP_D0   ,iS,iN_E,iN_X), &
                 OP_N(iOP_Chi  ,iS,iN_E,iN_X), &
                 OP_N(iOP_Sigma,iS,iN_E,iN_X), &
                 GX_N(iGF_Alpha   ,iN_X), &
                 GX_N(iGF_Beta_1  ,iN_X), &
                 GX_N(iGF_Beta_2  ,iN_X), &
                 GX_N(iGF_Beta_3  ,iN_X), &
                 dCR_N(iCR_N ,iS,iN_E,iN_X), &
                 dCR_N(iCR_G1,iS,iN_E,iN_X), &
                 dCR_N(iCR_G2,iS,iN_E,iN_X), &
                 dCR_N(iCR_G3,iS,iN_E,iN_X) )
                 

      END DO
      END DO

    END DO

    ! --- Revert Radiation Increment ---

    DO iS  = 1, nSpecies
    DO iCR = 1, nCR
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iE  = iE_B0, iE_E0

      DO iNodeZ = 1, nDOFZ

        iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1
        iNodeE = MOD( (iNodeZ-1)        , nDOFE ) + 1

        iN_X = iNodeX &
                 + (iX1-iX_B0(1)) * nDOFX &
                 + (iX2-iX_B0(2)) * nDOFX * nX(1) &
                 + (iX3-iX_B0(3)) * nDOFX * nX(1) * nX(2)

        iN_E = iNodeE &
                 + (iE-iE_B0) * nDOFE

        dU_R(iNodeZ,iE,iX1,iX2,iX3,iCR,iS) = dCR_N(iCR,iS,iN_E,iN_X)

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    CALL FinalizeCollisions
    
    IF (Verbose) THEN
      PRINT*, "      ComputeIncrement_TwoMoment_Implicit (End)"
    END IF

  END SUBROUTINE ComputeIncrement_TwoMoment_Implicit

!Will need to change this module
  SUBROUTINE ComputeIncrement_FixedPoint &
    ( dt, N, G_d_1, G_d_2, G_d_3, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, D_0, Chi, Sigma, &
      alp, B_u_1, B_u_2, B_u_3,  &
      dN, dG_d_1, dG_d_2, dG_d_3 )

    REAL(DP), INTENT(in)  :: dt
    REAL(DP), INTENT(in)  :: N, G_d_1, G_d_2, G_d_3
    REAL(DP), INTENT(in)  ::    V_u_1, V_u_2, V_u_3
    REAL(DP), INTENT(in)  :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), INTENT(in)  :: D_0, Chi, Sigma
    REAL(DP), INTENT(in)  :: B_u_1, B_u_2, B_u_3, alp
    REAL(DP), INTENT(out) :: dN, dG_d_1, dG_d_2, dG_d_3

    ! --- Parameters ---

    INTEGER,  PARAMETER :: M = 2
    INTEGER,  PARAMETER :: LWORK = 2 * M
    INTEGER,  PARAMETER :: MaxIterations = 100
    REAL(DP), PARAMETER :: Rtol = 1.0d-08

    ! --- Local Variables ---

    LOGICAL  :: CONVERGED
    INTEGER  :: i, k, mk, INFO
    REAL(DP) :: D, I_d_1, I_d_2, I_d_3, Kappa
    REAL(DP) ::    I_u_1, I_u_2, I_u_3
    REAL(DP) :: A_d_1, A_d_2, A_d_3
    REAL(DP) :: k_dd_ij(1:3,1:3)
    REAL(DP) :: D_00, D_ii
    REAL(DP) :: UVEC(4), CVEC(4)
    REAL(DP) :: GVEC(4,M), GVECm(4)
    REAL(DP) :: FVEC(4,M), FVECm(4)
    REAL(DP) :: LMAT(4,4), DET, Alpha(M)
    REAL(DP) :: BVEC(4), AMAT(4,M), WORK(LWORK)
    REAL(DP) :: W, DTE, B_d_1, B_d_2, B_d_3


    B_d_1 = Gm_dd_11 * B_u_1
    B_d_2 = Gm_dd_22 * B_u_2
    B_d_3 = Gm_dd_33 * B_u_3

    W = 1.0_DP / SQRT( 1.0_DP - (Gm_dd_11 * V_u_1 * V_u_1 &
               + Gm_dd_22 * V_u_2 * V_u_2 &  
               + Gm_dd_33 * V_u_3 * V_u_3) )

    DTE = 1.0_DP / ( B_d_1 * V_u_1 + B_d_2 * V_u_2 + B_d_3 * V_u_3 - alp )

    Kappa = Chi + Sigma

    D_00 = W + alp * dt * Chi
    D_ii = W + alp * dt * Kappa
    ! --- Constant Vector ---

    CVEC = [ N + alp * dt * Chi * D_0, G_d_1, G_d_2, G_d_3 ]

    ! --- Initial Guess ---

    D     = N
    I_u_1 = Zero
    I_u_2 = Zero
    I_u_3 = Zero

    I_d_1 = DTE * ( B_d_2 * V_u_2 + B_d_3 * V_u_3 - alp ) * Gm_dd_11 * I_u_1 &
          - DTE * ( B_d_1 * V_u_2 *Gm_dd_22 ) * I_u_2 - DTE * ( B_d_1 * V_u_3 * Gm_dd_33 ) * I_u_3 
    I_d_2 = DTE * ( B_d_1 * V_u_1 + B_d_3 * V_u_3 - alp ) * Gm_dd_22 * I_u_2 &
          - DTE * ( B_d_2 * V_u_1 * Gm_dd_11 ) * I_u_1 - DTE * ( Gm_dd_33 * I_u_3 * B_d_2 * V_u_3 ) 
    I_d_3 = DTE * ( B_d_1 * V_u_1 + B_d_2 * V_u_2 - alp ) * Gm_dd_33 * I_u_3 &
          - DTE * ( Gm_dd_11 * I_u_1 * B_d_3 * V_u_1 ) - DTE * ( Gm_dd_22 * I_u_2 * B_d_3 * V_u_2 )

    k = 0
    CONVERGED = .FALSE.
    DO WHILE( .NOT. CONVERGED .AND. k < MaxIterations )

      k  = k + 1
      mk = MIN( M, k )

      UVEC = [ D, I_d_1, I_d_2, I_d_3 ]

      CALL ComputeEddingtonTensorComponents_dd &
             ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
               alp, B_u_1, B_u_2, B_u_3, V_u_1, V_u_2, V_u_3, k_dd_ij  )

      A_d_1 = V_u_1 * k_dd_ij(1,1) + V_u_2 * k_dd_ij(1,2) + V_u_3 * k_dd_ij(1,3)
      A_d_2 = V_u_1 * k_dd_ij(1,2) + V_u_2 * k_dd_ij(2,2) + V_u_3 * k_dd_ij(2,3)
      A_d_3 = V_u_1 * k_dd_ij(1,3) + V_u_2 * k_dd_ij(2,3) + V_u_3 * k_dd_ij(3,3)

      DET = ( D_00 * D_ii &
              - ( V_u_1 * A_d_1 + V_u_2 * A_d_2 + V_u_3 * A_d_3 ) ) * D_ii**2

      LMAT(1,1) = D_ii**3
      LMAT(2,1) = - A_d_1 * D_ii**2
      LMAT(3,1) = - A_d_2 * D_ii**2
      LMAT(4,1) = - A_d_3 * D_ii**2

      LMAT(1,2) = - V_u_1 * D_ii**2
      LMAT(2,2) = D_00 * D_ii**2 &
                    - ( V_u_2 * A_d_2 + V_u_3 * A_d_3 ) * D_ii
      LMAT(3,2) = V_u_1 * A_d_2 * D_ii
      LMAT(4,2) = V_u_1 * A_d_3 * D_ii

      LMAT(1,3) = - V_u_2 * D_ii**2
      LMAT(2,3) = V_u_2 * A_d_1 * D_ii
      LMAT(3,3) = D_00 * D_ii**2 &
                    - ( V_u_1 * A_d_1 + V_u_3 * A_d_3 ) * D_ii
      LMAT(4,3) = V_u_2 * A_d_3 * D_ii

      LMAT(1,4) = - V_u_3 * D_ii**2
      LMAT(2,4) = V_u_3 * A_d_1 * D_ii
      LMAT(3,4) = V_u_3 * A_d_2 * D_ii
      LMAT(4,4) = D_00 * D_ii**2 &
                    - ( V_u_1 * A_d_1 + V_u_2 * A_d_2 ) * D_ii

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
    dN     = alp * Chi * ( D_0 - D )
    dG_d_1 = - alp * Kappa * I_d_1
    dG_d_2 = - alp * Kappa * I_d_2
    dG_d_3 = - alp * Kappa * I_d_3
    IF( k == MaxIterations )THEN

      PRINT*
      PRINT*, "ComputeIncrement_FixedPoint"
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

  END SUBROUTINE ComputeIncrement_FixedPoint


  SUBROUTINE InitializeCollisions( iZ_B0, iZ_E0, iZ_B1, iZ_E1 )

    INTEGER, INTENT(in) :: iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)

    iE_B0 = iZ_B0(1);   iE_E0 = iZ_E0(1)
    iX_B0 = iZ_B0(2:4); iX_E0 = iZ_E0(2:4)

    nZ = iZ_E0 - iZ_B0 + 1
    nE = iE_E0 - iE_B0 + 1
    nX = iX_E0 - iX_B0 + 1

    nE_G = nDOFE * nE
    nX_G = nDOFX * PRODUCT( nX )

    ALLOCATE( GX_N(nGF,nX_G) )
    ALLOCATE( CF_N(nCF,nX_G) )
    ALLOCATE( CR_N (nCR,nSpecies,nE_G,nX_G) )
    ALLOCATE( dCR_N(nCR,nSpecies,nE_G,nX_G) )
    ALLOCATE( OP_N (nOP,nSpecies,nE_G,nX_G) )

  END SUBROUTINE InitializeCollisions


  SUBROUTINE FinalizeCollisions

    DEALLOCATE( GX_N, CF_N, CR_N, dCR_N, OP_N )

  END SUBROUTINE FinalizeCollisions


END MODULE TwoMoment_DiscretizationModule_Collisions
