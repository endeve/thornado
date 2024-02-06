MODULE TwoMoment_DiscretizationModule_Collisions_FMC

  USE KindModule, ONLY: &
    DP, Zero, One, Half, Three
  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    nDOFE, &
    nDOFZ
  USE FluidFieldsModule, ONLY: &
    nPF, iPF_V1, iPF_V2, iPF_V3
  USE GeometryFieldsModuleE, ONLY: &
    nGE
  USE GeometryFieldsModule, ONLY: &
    nGF, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33
  USE TwoMoment_FieldsModule_FMC, ONLY: &
    nSpecies, &
    nCM, iCM_E, iCM_F1, iCM_F2, iCM_F3
  USE TwoMoment_UtilitiesModule_FMC, ONLY: &
    EddingtonTensorComponents_dd, ComputeHatMomentsFromConserved
  USE TwoMoment_OpacityModule_FMC, ONLY: &
    uOP, iOP_J0, iOP_Chi, iOP_Sigma, nOP
  USE TwoMoment_ClosureModule, ONLY: &
    FluxFactor, &
    EddingtonFactor, &
    HeatFluxFactor
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeIncrement_TwoMoment_Implicit

  INTEGER :: iE_B0,    iE_E0
  INTEGER :: iX_B0(3), iX_E0(3)
  INTEGER :: nZ(4), nE, nX(3), nE_G, nX_G

  INTEGER , ALLOCATABLE :: nIterations(:,:,:)
  REAL(DP), ALLOCATABLE :: GX_N(:,:)
  REAL(DP), ALLOCATABLE :: PF_N(:,:)
  REAL(DP), ALLOCATABLE :: CM_N(:,:,:,:)
  REAL(DP), ALLOCATABLE :: dCM_N(:,:,:,:)
  REAL(DP), ALLOCATABLE :: OP_N(:,:,:,:)

  CONTAINS

  SUBROUTINE ComputeIncrement_TwoMoment_Implicit &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, dt, GE, GX, U_F, U_M, dU_M)

    ! --- Input/Output variables ---
    INTEGER, INTENT(in) :: iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in) :: dt
    REAL(DP), INTENT(in) :: GE(1:nDOFE, iZ_B1(1):iZ_E1(1), 1:nGE)
    REAL(DP), INTENT(in) :: GX(1:nDOFX, &
                               iZ_B1(2):iZ_E1(2), &
                               iZ_B1(3):iZ_E1(3), &
                               iZ_B1(4):iZ_E1(4), &
                               1:nGF)
    REAL(DP), INTENT(in) :: U_F(1:nDOFX, &
                                iZ_B1(2):iZ_E1(2), &
                                iZ_B1(3):iZ_E1(3), &
                                iZ_B1(4):iZ_E1(4), &
                                1:nPF)
    REAL(DP), INTENT(in) :: U_M(1:nDOFZ, &
                                iZ_B1(1):iZ_E1(1), &
                                iZ_B1(2):iZ_E1(2), &
                                iZ_B1(3):iZ_E1(3), &
                                iZ_B1(4):iZ_E1(4), &
                                1:nCM, &
                                1:nSpecies)
    REAL(DP), INTENT(out) :: dU_M(1:nDOFZ, &
                                  iZ_B1(1):iZ_E1(1), &
                                  iZ_B1(2):iZ_E1(2), &
                                  iZ_B1(3):iZ_E1(3), &
                                  iZ_B1(4):iZ_E1(4), &
                                  1:nCM, &
                                  1:nSpecies)

    INTEGER :: iZ1, iZ2, iZ3, iZ4, iCM, iS, iGF, iPF, iOP
    INTEGER :: iX1, iX2, iX3, iE
    INTEGER :: iNodeZ, iNodeX, iNodeE, iN_X, iN_E

    Write(*,*)
    print *,'ComputeIncrement_TwoMoment_Implicit'

    CALL InitializeCollisions( iZ_B0, iZ_E0, iZ_B1, iZ_E1 )

    DO iS  = 1, nSpecies
    DO iCM = 1, nCM
    DO iZ4 = iZ_B1(4), iZ_E1(4)
    DO iZ3 = iZ_B1(3), iZ_E1(3)
    DO iZ2 = iZ_B1(2), iZ_E1(2)
    DO iZ1 = iZ_B1(1), iZ_E1(1)

      DO iNodeZ = 1, nDOFZ

        dU_M(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM,iS) = Zero

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
    DO iPF  = 1, nPF

      iX3    = MOD( (iN_X-1) / ( nDOFX * nX(1) * nX(2) ), nX(3) ) + iX_B0(3)
      iX2    = MOD( (iN_X-1) / ( nDOFX * nX(1)         ), nX(2) ) + iX_B0(2)
      iX1    = MOD( (iN_X-1) / ( nDOFX                 ), nX(1) ) + iX_B0(1)
      iNodeX = MOD( (iN_X-1)                            , nDOFX ) + 1

      PF_N(iPF,iN_X) = U_F(iNodeX,iX1,iX2,iX3,iPF)

    END DO
    END DO

    ! --- Arrange Two-Moment Fields ---

    DO iS   = 1, nSpecies
    DO iCM  = 1, nCM
    DO iN_X = 1, nX_G
    DO iN_E = 1, nE_G

      iE     = MOD( (iN_E-1) / nDOFE, nZ(1) ) + iE_B0
      iNodeE = MOD( (iN_E-1)        , nDOFE ) + 1

      iX3    = MOD( (iN_X-1) / ( nDOFX * nX(1) * nX(2) ), nX(3) ) + iX_B0(3)
      iX2    = MOD( (iN_X-1) / ( nDOFX * nX(1)         ), nX(2) ) + iX_B0(2)
      iX1    = MOD( (iN_X-1) / ( nDOFX                 ), nX(1) ) + iX_B0(1)
      iNodeX = MOD( (iN_X-1)                            , nDOFX ) + 1

      iNodeZ = iNodeE + ( iNodeX - 1 ) * nDOFE

      CM_N(iCM,iS,iN_E,iN_X) = U_M(iNodeZ,iE,iX1,iX2,iX3,iCM,iS)

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
    DO iN_E = 1, nE_G
    DO iS   = 1, nSpecies

      CALL ComputeIncrement_FixedPoint_Richardson &
             ( dt, &
               CM_N (iCM_E       ,iS,iN_E,iN_X), &
               CM_N (iCM_F1      ,iS,iN_E,iN_X), &
               CM_N (iCM_F2      ,iS,iN_E,iN_X), &
               CM_N (iCM_F3      ,iS,iN_E,iN_X), &
               PF_N (iPF_V1              ,iN_X), &
               PF_N (iPF_V2              ,iN_X), &
               PF_N (iPF_V3              ,iN_X), &
               GX_N (iGF_Gm_dd_11        ,iN_X), &
               GX_N (iGF_Gm_dd_22        ,iN_X), &
               GX_N (iGF_Gm_dd_33        ,iN_X), &
               OP_N (iOP_J0      ,iS,iN_E,iN_X), &
               OP_N (iOP_Chi     ,iS,iN_E,iN_X), &
               OP_N (iOP_Sigma   ,iS,iN_E,iN_X), &
               dCM_N(iCM_E       ,iS,iN_E,iN_X), &
               dCM_N(iCM_F1      ,iS,iN_E,iN_X), &
               dCM_N(iCM_F2      ,iS,iN_E,iN_X), &
               dCM_N(iCM_F3      ,iS,iN_E,iN_X), &
               nIterations(       iS,iN_E,iN_X) )

    END DO
    END DO
    END DO

    ! --- Revert Two-Moment Increment ---

    DO iS  = 1, nSpecies
    DO iCM = 1, nCM
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

        dU_M(iNodeZ,iE,iX1,iX2,iX3,iCM,iS) = dCM_N(iCM,iS,iN_E,iN_X)

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    CALL FinalizeCollisions

  END SUBROUTINE ComputeIncrement_TwoMoment_Implicit
  
  SUBROUTINE ComputeIncrement_FixedPoint_Richardson &
    ( dt, E, F_d_1, F_d_2, F_d_3, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, J_0, Chi, Sigma, &
      dE, dF_d_1, dF_d_2, dF_d_3, nIterations )

    REAL(DP), INTENT(in)  :: dt
    REAL(DP), INTENT(in)  :: E, F_d_1, F_d_2, F_d_3
    REAL(DP), INTENT(in)  ::    V_u_1, V_u_2, V_u_3
    REAL(DP), INTENT(in)  :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), INTENT(in)  :: J_0, Chi, Sigma
    REAL(DP), INTENT(out) :: dE, dF_d_1, dF_d_2, dF_d_3
    INTEGER , INTENT(out) :: nIterations

    ! --- Parameters ---

    INTEGER,  PARAMETER :: MaxIterations = 100
    REAL(DP), PARAMETER :: TOL = 1.0d-08

    ! --- Local Variables ---

    LOGICAL  :: CONVERGED
    INTEGER :: iteration
    REAL(DP) :: J, H_d_1, H_d_2, H_d_3
    REAL(DP) :: Kappa, Mu_Chi, Mu_Kappa
    ! REAL(DP) ::    H_u_1, H_u_2, H_u_3
    REAL(DP) :: k_dd(3,3)
    REAL(DP) :: vMagSq, vDotH, vDotk1, vDotk2, vDotk3, lambda, W
    REAL(DP) :: E_hat, F_hat_d_1, F_hat_d_2, F_hat_d_3
    REAL(DP) :: fvec_0, fvec_1, fvec_2, fvec_3

    Kappa = Chi + Sigma

    ! --- Initial guess ---

    CALL ComputeHatMomentsFromConserved &
      ( E, F_d_1, F_d_2, F_d_3, E_hat, F_hat_d_1, F_hat_d_2, F_hat_d_3, &
        V_u_1, V_u_2, V_u_3, &
        Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    ! --- Richardson update ---

    vMagSq = V_u_1 * Gm_dd_11 * V_u_1 &
           + V_u_2 * Gm_dd_22 * V_u_2 &
           + V_u_3 * Gm_dd_33 * V_u_3
    W = One / SQRT ( One - vMagSq )
    lambda = One / ( One + SQRT( vMagSq ) )

    Mu_Chi = One / (W + lambda * dt * Chi)
    Mu_Kappa = One / (W + lambda * dt * Kappa)
    
    J = E_hat
    H_d_1 = F_hat_d_1
    H_d_2 = F_hat_d_2
    H_d_3 = F_hat_d_3

    iteration = 0
    CONVERGED = .FALSE.
    DO WHILE( .NOT. CONVERGED .AND. iteration < MaxIterations )

      iteration = iteration + 1

      k_dd = EddingtonTensorComponents_dd &
        ( J, H_d_1, H_d_2, H_d_3, V_u_1, V_u_2, V_u_3, &
          Gm_dd_11, Gm_dd_22, Gm_dd_33 )
      
      vDotH = V_u_1 * H_d_1 + V_u_2 * H_d_2 + V_u_3 * H_d_3

      vDotk1 = V_u_1 * k_dd(1,1) + V_u_2 * k_dd(1,2) + V_u_3 * k_dd(1,3)
      vDotk2 = V_u_1 * k_dd(2,1) + V_u_2 * k_dd(2,2) + V_u_3 * k_dd(2,3)
      vDotk3 = V_u_1 * k_dd(3,1) + V_u_2 * k_dd(3,2) + V_u_3 * k_dd(3,3)

      ! --- Compute components of vector function f ---
      ! If converged, f = 0

      fvec_0 = E_hat + dt * Chi * J_0 - (W + dt * Chi) * J - vDotH
      fvec_1 = F_hat_d_1 - (W + dt * Kappa) * H_d_1 - vDotk1 * J
      fvec_2 = F_hat_d_2 - (W + dt * Kappa) * H_d_2 - vDotk2 * J
      fvec_3 = F_hat_d_3 - (W + dt * Kappa) * H_d_3 - vDotk3 * J

      ! --- Compute next iterate ---

      J = J + lambda * Mu_Chi * fvec_0
      H_d_1 = H_d_1 + lambda * Mu_Kappa * fvec_1
      H_d_2 = H_d_2 + lambda * Mu_Kappa * fvec_2
      H_d_3 = H_d_3 + lambda * Mu_Kappa * fvec_3
      
      CONVERGED = SQRT( fvec_0**2 + fvec_1**2 + fvec_2**2 +fvec_3**2 ) < TOL

    END DO

    dE = W * Chi * (J_0 - J) - Kappa * (V_u_1 * H_d_1 + V_u_2 * H_d_2 + V_u_3 * H_d_3)
    dF_d_1 = W * V_u_1 * Chi * (J_0 - J) - Kappa * H_d_1
    dF_d_2 = W * V_u_2 * Chi * (J_0 - J) - Kappa * H_d_2
    dF_d_3 = W * V_u_3 * Chi * (J_0 - J) - Kappa * H_d_3
    ! dE = 1000000.0_DP
    ! dF_d_1 = 1000000.0_DP
    ! dF_d_2 = 1000000.0_DP
    ! dF_d_3 = 1000000.0_DP

    ! print *, dE
    ! print *, dF_d_1
    ! Write(*,*)

    nIterations = iteration

    IF( iteration >= MaxIterations )THEN

      print *, 'Warning! Max iterations reached.'

    END IF

  END SUBROUTINE ComputeIncrement_FixedPoint_Richardson

  SUBROUTINE InitializeCollisions( iZ_B0, iZ_E0, iZ_B1, iZ_E1 )

    INTEGER, INTENT(in) :: iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)

    iE_B0 = iZ_B0(1);   iE_E0 = iZ_E0(1)
    iX_B0 = iZ_B0(2:4); iX_E0 = iZ_E0(2:4)

    nZ = iZ_E0 - iZ_B0 + 1
    nE = iE_E0 - iE_B0 + 1
    nX = iX_E0 - iX_B0 + 1

    nE_G = nDOFE * nE
    nX_G = nDOFX * PRODUCT( nX )

    ALLOCATE( nIterations(nSpecies,nE_G,nX_G) )
    ALLOCATE( GX_N(nGF,nX_G) )
    ALLOCATE( PF_N(nPF,nX_G) )
    ALLOCATE( CM_N (nCM,nSpecies,nE_G,nX_G) )
    ALLOCATE( dCM_N(nCM,nSpecies,nE_G,nX_G) )
    ALLOCATE( OP_N (nOP,nSpecies,nE_G,nX_G) )

  END SUBROUTINE InitializeCollisions


  SUBROUTINE FinalizeCollisions

    DEALLOCATE( nIterations, GX_N, PF_N, CM_N, dCM_N, OP_N )

  END SUBROUTINE FinalizeCollisions

END MODULE TwoMoment_DiscretizationModule_Collisions_FMC