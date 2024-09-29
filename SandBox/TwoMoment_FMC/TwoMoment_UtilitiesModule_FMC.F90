MODULE TwoMoment_UtilitiesModule_FMC

  USE KindModule, ONLY: &
    DP, Third, Half, Zero, One, Two, Three, Four, Five, SqrtTiny
  USE QuadratureModule, ONLY: &
    GetQuadrature
  USE UtilitiesModule, ONLY: &
    WriteVector
  USE ProgramHeaderModule, ONLY: &
    nDOFZ, nDOFX, nDOFE, &
    nNOdes, nDimsX, nDims
  USE ReferenceElementModuleX, ONLY: &
    nDOFX_X1, nDOFX_X2, nDOFX_X3, &
    WeightsX_q, &
    WeightsX1, &
    WeightsX_X1, &
    WeightsX_X2, &
    WeightsX_X3, &
    NodeNumberTableX
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    dLXdX1_q, LX_X1_Dn, LX_X1_Up, &
    dLXdX2_q, LX_X2_Dn, LX_X2_Up, &
    dLXdX3_q, LX_X3_Dn, LX_X3_Up
  USE MeshModule, ONLY: &
    MeshX, &
    MeshE, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    nGF, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33, &
    iGF_h_1, iGF_h_2, iGF_h_3, iGF_SqrtGm
  USE FluidFieldsModule, ONLY: &
    nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3
  USE Euler_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_Euler
  USE TwoMoment_FieldsModule_FMC, ONLY: &
    nSpecies, &
    nCM, iCM_E, iCM_F1, iCM_F2, iCM_F3, &
    nPM, iPM_J, iPM_H1, iPM_H2, iPM_H3, &
    nAM, &
    nGM
  USE TwoMoment_ClosureModule, ONLY: &
    FluxFactor_Relativistic, EddingtonFactor, HeatFluxFactor
  USE LinearAlgebraModule, ONLY: &
    MatrixMatrixMultiply, &
    EigenvaluesSymmetric3

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Initialize_MomentConversion
  PUBLIC :: ComputeConserved_TwoMoment_FMC
  PUBLIC :: ComputeFromConserved_TwoMoment_FMC
  PUBLIC :: ComputePrimitive_TwoMoment_FMC
  PUBLIC :: ComputePrimitive_TwoMoment_Richardson_FMC
  PUBLIC :: ComputePrimitive_TwoMoment_Newton_FMC
  PUBLIC :: MomentConversion_Jacobian_Inverse
  PUBLIC :: ComputeHatMomentsFromConserved
  PUBLIC :: ComputeHatMomentsFromPrimitive
  PUBLIC :: EddingtonTensorComponents_dd
  PUBLIC :: HeatFluxTensorComponents_uuu
  PUBLIC :: ComputeHeatFluxTensorComponents_ddd_Lagrangian
  PUBLIC :: ComputeHeatFluxTensorComponents_uud_Lagrangian
  PUBLIC :: Flux_X1
  PUBLIC :: Flux_X2
  PUBLIC :: Flux_X3
  PUBLIC :: Flux_E
  PUBLIC :: NumericalFlux_LLF
  PUBLIC :: ComputeTimeStep_TwoMoment
  PUBLIC :: ComputeTimeStep_TwoMoment_Realizable
  PUBLIC :: ComputeWeakDerivatives_X1
  PUBLIC :: ComputeWeakDerivatives_X2
  PUBLIC :: ComputeWeakDerivatives_X3
  PUBLIC :: FaceVelocity_X1
  PUBLIC :: FaceFourVelocity_X1

  LOGICAL :: Newton_Solver

CONTAINS

  SUBROUTINE Initialize_MomentConversion ( Newtons_Option )

    LOGICAL, INTENT(in), OPTIONAL :: Newtons_Option

    Newton_Solver = .FALSE.
    IF( PRESENT(Newtons_Option) )THEN

      Newton_Solver = Newtons_Option

    END IF

    WRITE(*,*)
    WRITE(*,'(A)') '  INFO: Initialize_MomentConversion:'
    WRITE(*,'(A)') '  --------------------------------------------'
    WRITE(*,*)
    WRITE(*,'(A4,A32,L1)') &
      '', 'Use Newton''s Method: ', Newton_Solver

  END SUBROUTINE Initialize_MomentConversion

  SUBROUTINE ComputeConserved_TwoMoment_FMC &
    ( J, H_d_1, H_d_2, H_d_3, E, F_d_1, F_d_2, F_d_3, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    ! --- Input/Output variables ---
    REAL(DP), INTENT(in)  :: J, H_d_1, H_d_2, H_d_3 ! --- Index down change to up
    REAL(DP), INTENT(out) :: E, F_d_1, F_d_2, F_d_3 ! --- Index down
    REAL(DP), INTENT(in)  ::    V_u_1, V_u_2, V_u_3 ! --- Index up 
    REAL(DP), INTENT(in)  :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    ! --- Local variables ---
    REAL(DP) :: k_dd(3,3), vMagSq, W, vDotFhat
    REAL(DP) :: E_hat, F_hat_d_1, F_hat_d_2, F_hat_d_3 ! --- Index down

    k_dd = EddingtonTensorComponents_dd &
      ( J, H_d_1, H_d_2, H_d_3, V_u_1, V_u_2, V_u_3, &
        Gm_dd_11, Gm_dd_22, Gm_dd_33)

    vMagSq = V_u_1 * Gm_dd_11 * V_u_1 &
           + V_u_2 * Gm_dd_22 * V_u_2 &
           + V_u_3 * Gm_dd_33 * V_u_3
    W = One / SQRT ( One - vMagSq )

    CALL ComputeHatMomentsFromPrimitive &
           ( J, H_d_1, H_d_2, H_d_3, E_hat, F_hat_d_1, F_hat_d_2, F_hat_d_3, &
             V_u_1, V_u_2, V_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33 ) !change H indices to up

    vDotFhat = V_u_1 * F_hat_d_1 + V_u_2 * F_hat_d_2 + V_u_3 * F_hat_d_3

    ! --- Eulerian Energy Density ---
    
    E = W * E_hat + vDotFhat

    ! --- Eulerian Momentum Density (1) ---
    
    F_d_1 = F_hat_d_1 + W * Gm_dd_11 * V_u_1 * E_hat

    ! --- Eulerian Momentum Density (2) ---
        
    F_d_2 = F_hat_d_2 + W * Gm_dd_22 * V_u_2 * E_hat

    ! --- Eulerian Momentum Density (3) ---
        
    F_d_3 = F_hat_d_3 + W * Gm_dd_33 * V_u_3 * E_hat

  END SUBROUTINE ComputeConserved_TwoMoment_FMC

  SUBROUTINE ComputeFromConserved_TwoMoment_FMC &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, PF, CM, PM, AM, GM )
    
    ! --- Input/Output variables ---
    INTEGER,  INTENT(in)  :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)  :: &
      GX(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4), &
         1:nGF)
    REAL(DP), INTENT(in)  :: &
      PF(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4), &
         1:nPF)
    REAL(DP), INTENT(in)  :: &
      CM(1:nDOFZ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4),1:nCM,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      PM(1:nDOFZ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4),1:nPM,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      AM(1:nDOFZ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4),1:nAM,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      GM(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4), &
         1:nGM,1:nSpecies)

    ! --- Local variables ---
    INTEGER  :: &
      iZ1, iZ2, iZ3, iZ4, iS, iNodeZ, iNodeX

    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeZ = 1, nDOFZ

        iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1

        ! IF ( CM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_E ,iS) < &
        !      SQRT ( CM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_F1 ,iS)**2 + &
        !             CM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_F2 ,iS)**2 + &
        !             CM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_F3 ,iS)**2 ) ) THEN
          ! print *, 'E  = ', CM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_E ,iS)
          ! print *, 'F1 = ', CM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_F1 ,iS)
          ! STOP
        ! END IF

        CALL ComputePrimitive_TwoMoment_FMC &
               ( CM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_E ,iS), &
                 CM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_F1,iS), &
                 CM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_F2,iS), &
                 CM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_F3,iS), &
                 PM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_J ,iS), &
                 PM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H1,iS), &
                 PM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H2,iS), &
                 PM(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPM_H3,iS), &
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

  END SUBROUTINE ComputeFromConserved_TwoMoment_FMC

  SUBROUTINE ComputePrimitive_TwoMoment_FMC &
    ( E, F_d_1, F_d_2, F_d_3, J, H_d_1, H_d_2, H_d_3, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, nIterations_Option )

    ! --- Input/Output variables ---
    REAL(DP), INTENT(in) :: E, F_d_1, F_d_2, F_d_3
    REAL(DP), INTENT(out) :: J, H_d_1, H_d_2, H_d_3
    REAL(DP), INTENT(in) :: V_u_1, V_u_2, V_u_3
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    ! LOGICAL, INTENT(in) :: UseNewtons
    INTEGER, INTENT(out), OPTIONAL :: nIterations_Option

    IF( Newton_Solver )THEN

      CALL ComputePrimitive_TwoMoment_Newton_FMC &
        ( E, F_d_1, F_d_2, F_d_3, J, H_d_1, H_d_2, H_d_3, V_u_1, V_u_2, V_u_3, &
          Gm_dd_11, Gm_dd_22, Gm_dd_33, nIterations_Option )

    ELSE

      CALL ComputePrimitive_TwoMoment_Richardson_FMC &
        ( E, F_d_1, F_d_2, F_d_3, J, H_d_1, H_d_2, H_d_3, V_u_1, V_u_2, V_u_3, &
          Gm_dd_11, Gm_dd_22, Gm_dd_33, nIterations_Option )

    END IF

  END SUBROUTINE ComputePrimitive_TwoMoment_FMC

  SUBROUTINE ComputePrimitive_TwoMoment_Richardson_FMC &
    ( E, F_d_1, F_d_2, F_d_3, J, H_d_1, H_d_2, H_d_3, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, nIterations_Option )

    ! --- Input/Output variables ---
    REAL(DP), INTENT(in) :: E, F_d_1, F_d_2, F_d_3
    REAL(DP), INTENT(out) :: J, H_d_1, H_d_2, H_d_3
    REAL(DP), INTENT(in) :: V_u_1, V_u_2, V_u_3
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    INTEGER, INTENT(out), OPTIONAL :: nIterations_Option

    ! --- Parameters ---
    INTEGER, PARAMETER :: MaxIterations = 200
    REAL(DP), PARAMETER :: TOL = 1.0d-08

    ! --- Local variables ---
    INTEGER :: iteration
    REAL(DP) :: vMagSq, vDotH, vDotk1, vDotk2, vDotk3, lambda, W
    REAL(DP) :: E_hat, F_hat_d_1, F_hat_d_2, F_hat_d_3
    REAL(DP) :: fvec_0, fvec_1, fvec_2, fvec_3
    REAL(DP) :: k_dd(3,3)

    LOGICAL :: CONVERGED

    ! --- Initial guess ---

    CALL ComputeHatMomentsFromConserved &
      ( E, F_d_1, F_d_2, F_d_3, E_hat, F_hat_d_1, F_hat_d_2, F_hat_d_3, &
        V_u_1, V_u_2, V_u_3, &
        Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    ! print *, J, H_d_1

    ! J = SQRT( F_d_1**2 + F_d_2**2 + F_d_3**2 )
    ! H_d_1 = F_d_1
    ! H_d_2 = F_d_2
    ! H_d_3 = F_d_3

    ! --- Richardson update ---

    vMagSq = V_u_1 * Gm_dd_11 * V_u_1 &
           + V_u_2 * Gm_dd_22 * V_u_2 &
           + V_u_3 * Gm_dd_33 * V_u_3
    W = One / SQRT ( One - vMagSq )
    lambda = One / ( W + W * SQRT( vMagSq ) )
    
    J = E_hat / W
    H_d_1 = F_hat_d_1 / W
    H_d_2 = F_hat_d_2 / W
    H_d_3 = F_hat_d_3 / W

    iteration = 0
    CONVERGED = .FALSE.
    DO WHILE( .NOT. CONVERGED .AND. iteration < MaxIterations )

      iteration = iteration + 1

      ! ! print *, "V_u_1 = ", V_u_1
      ! print *, "J     = ", J
      ! print *, "H_d_1 = ", H_d_1
      ! Write(*,*)

      k_dd = EddingtonTensorComponents_dd &
        ( J, H_d_1, H_d_2, H_d_3, V_u_1, V_u_2, V_u_3, &
          Gm_dd_11, Gm_dd_22, Gm_dd_33 )
      
      vDotH = V_u_1 * H_d_1 + V_u_2 * H_d_2 + V_u_3 * H_d_3

      vDotk1 = V_u_1 * k_dd(1,1) + V_u_2 * k_dd(1,2) + V_u_3 * k_dd(1,3)
      vDotk2 = V_u_1 * k_dd(2,1) + V_u_2 * k_dd(2,2) + V_u_3 * k_dd(2,3)
      vDotk3 = V_u_1 * k_dd(3,1) + V_u_2 * k_dd(3,2) + V_u_3 * k_dd(3,3)

      ! --- Compute components of vector function f ---
      ! If converged, f = 0

      fvec_0 = W * J + vDotH - E_hat
      fvec_1 = W * H_d_1 + vDotk1 * J - F_hat_d_1
      fvec_2 = W * H_d_2 + vDotk2 * J - F_hat_d_2
      fvec_3 = W * H_d_3 + vDotk3 * J - F_hat_d_3

      ! --- Compute next iterate ---

      J = J - lambda * fvec_0
      H_d_1 = H_d_1 - lambda * fvec_1
      H_d_2 = H_d_2 - lambda * fvec_2
      H_d_3 = H_d_3 - lambda * fvec_3
      
      CONVERGED = SQRT( fvec_0**2 + fvec_1**2 + fvec_2**2 +fvec_3**2 ) < TOL

    END DO

    IF( PRESENT( nIterations_Option ) ) THEN

      nIterations_Option = iteration

    END IF

    IF( iteration >= MaxIterations )THEN

      print *, 'Warning! Max iterations reached.'

    END IF

  END SUBROUTINE ComputePrimitive_TwoMoment_Richardson_FMC

  SUBROUTINE ComputePrimitive_TwoMoment_Newton_FMC &
    ( E, F_d_1, F_d_2, F_d_3, J, H_d_1, H_d_2, H_d_3, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, nIterations_Option )

    ! --- Input/Output variables ---
    REAL(DP), INTENT(in) :: E, F_d_1, F_d_2, F_d_3
    REAL(DP), INTENT(out) :: J, H_d_1, H_d_2, H_d_3
    REAL(DP), INTENT(in) :: V_u_1, V_u_2, V_u_3
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    INTEGER, INTENT(out), OPTIONAL :: nIterations_Option

    ! --- Parameters ---
    INTEGER, PARAMETER :: MaxIterations = 200
    REAL(DP), PARAMETER :: TOL = 1.0d-08

    ! --- Local variables ---
    INTEGER  :: iteration
    REAL(DP) :: vMagSq, vDotH, vDotk1, vDotk2, vDotk3, W
    REAL(DP) :: E_hat, F_hat_d_1, F_hat_d_2, F_hat_d_3
    REAL(DP) :: fvec_0, fvec_1, fvec_2, fvec_3
    REAL(DP) :: k_dd(3,3), Jacobian_Inv(4,4)
    LOGICAL  :: CONVERGED

    ! --- Initial guess ---

    CALL ComputeHatMomentsFromConserved &
      ( E, F_d_1, F_d_2, F_d_3, E_hat, F_hat_d_1, F_hat_d_2, F_hat_d_3, &
        V_u_1, V_u_2, V_u_3, &
        Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    ! --- Newton update ---

    vMagSq = V_u_1 * Gm_dd_11 * V_u_1 &
           + V_u_2 * Gm_dd_22 * V_u_2 &
           + V_u_3 * Gm_dd_33 * V_u_3
    W = One / SQRT ( One - vMagSq )
    
    J = E_hat / W
    H_d_1 = F_hat_d_1 / W
    H_d_2 = F_hat_d_2 / W
    H_d_3 = F_hat_d_3 / W

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

      fvec_0 = W * J + vDotH - E_hat
      fvec_1 = W * H_d_1 + vDotk1 * J - F_hat_d_1
      fvec_2 = W * H_d_2 + vDotk2 * J - F_hat_d_2
      fvec_3 = W * H_d_3 + vDotk3 * J - F_hat_d_3

      ! --- Compute next iterate ---

      Jacobian_Inv = MomentConversion_Jacobian_Inverse &
        ( J, H_d_1, H_d_2, H_d_3, k_dd, V_u_1, V_u_2, V_u_3, &
          Gm_dd_11, Gm_dd_22, Gm_dd_33 )

      ! print *, 'Inverse Matrix', Jacobian_Inv
      ! STOP

      J     = - (Jacobian_Inv(1,1) * fvec_0 + Jacobian_Inv(1,2) * fvec_1 + &
                 Jacobian_Inv(1,3) * fvec_2 + Jacobian_Inv(1,4) * fvec_3 ) + J
      H_d_1 = - (Jacobian_Inv(2,1) * fvec_0 + Jacobian_Inv(2,2) * fvec_1 + &
                 Jacobian_Inv(2,3) * fvec_2 + Jacobian_Inv(2,4) * fvec_3 ) + H_d_1
      H_d_2 = - (Jacobian_Inv(3,1) * fvec_0 + Jacobian_Inv(3,2) * fvec_1 + &
                 Jacobian_Inv(3,3) * fvec_2 + Jacobian_Inv(3,4) * fvec_3 ) + H_d_2
      H_d_3 = - (Jacobian_Inv(4,1) * fvec_0 + Jacobian_Inv(4,2) * fvec_1 + &
                 Jacobian_Inv(4,3) * fvec_2 + Jacobian_Inv(4,4) * fvec_3 ) + H_d_3
      
      CONVERGED = SQRT( fvec_0**2 + fvec_1**2 + fvec_2**2 +fvec_3**2 ) < TOL

    END DO

    IF( PRESENT( nIterations_Option ) ) THEN

      nIterations_Option = iteration

    END IF

    IF( iteration >= MaxIterations )THEN

      print *, 'Warning! Max iterations reached.'

    END IF

    IF ( J - SQRT(-(V_u_1*H_d_1+V_u_2*H_d_2+V_u_3*H_d_3)**2 + H_d_1**2+H_d_2**2+H_d_3**2)<0 )THEN

      print *, 'Warning: non-realizable iterate'

    END IF

  END SUBROUTINE ComputePrimitive_TwoMoment_Newton_FMC

  FUNCTION MomentConversion_Jacobian_Inverse &
    ( J, H_d_1, H_d_2, H_d_3, k_dd, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    ! --- Input/Output variables
    REAL(DP), INTENT(in) :: J, H_d_1, H_d_2, H_d_3
    REAL(DP), INTENT(in) :: k_dd(3,3)
    REAL(DP), INTENT(in) :: V_u_1, V_u_2, V_u_3
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP) :: MomentConversion_Jacobian_Inverse(4,4)

    ! --- Local variables ---
    INTEGER  :: i, k
    REAL(DP) :: FF, EF, W, vMagSq, HSq
    REAL(DP) :: h_hat_d(0:3), h_hat_d_H1(1:3), h_hat_d_H2(1:3), h_hat_d_H3(1:3)
    REAL(DP) :: u_d(3), h_dd(3,3), Gm_dd(3,3)
    REAL(DP) :: EF_J, EF_H1, EF_H2, EF_H3
    REAL(DP) :: K_dd_J(3,3), K_dd_H1(3,3), K_dd_H2(3,3), K_dd_H3(3,3)
    REAL(DP) :: vDotK1_J, vDotK2_J, vDotK3_J
    REAL(DP) :: vDotK1_H1, vDotK2_H1, vDotK3_H1
    REAL(DP) :: vDotK1_H2, vDotK2_H2, vDotK3_H2
    REAL(DP) :: vDotK1_H3, vDotK2_H3, vDotK3_H3
    REAL(DP) :: Jacobian(4,4), detinv

    HSq = - ( V_u_1 * H_d_1 + V_u_2 * H_d_2 + V_u_3 * H_d_3 )**2 + &
          H_d_1**2 + H_d_2**2 + H_d_3**2

    FF = MIN( MAX( SQRT( HSq ) / MAX( ABS( J ), SqrtTiny ), &
                   SqrtTiny ), &
              One )

    IF ( FF <= SqrtTiny ) THEN
      EF = Third
    ELSE
      EF = EddingtonFactor( J, FF )
    END IF

    ! print *, 'FF ', FF
    ! print *, 'EF ', EF
    ! print *, 'HSQ', HSQ

    h_hat_d(1) = H_d_1 / ( MAX( FF * ABS(J), SqrtTiny ) ) ! Was running in to division by zero.
    h_hat_d(2) = H_d_2 / ( MAX( FF * ABS(J), SqrtTiny ) ) ! Is this a good fix?
    h_hat_d(3) = H_d_3 / ( MAX( FF * ABS(J), SqrtTiny ) )
    h_hat_d(0) = - ( V_u_1 * h_hat_d(1) + V_u_2 * h_hat_d(2) + V_u_3 * h_hat_d(3) )

    ! --- Partial derivatives of \hat{h}_i with respect of \mathcal{H}_i ---
    h_hat_d_H1(1) = ( SQRT(HSq) - H_d_1 * (h_hat_d(1) + h_hat_d(0) * V_u_1) ) / MAX(HSq,SqrtTiny)
    h_hat_d_H1(2) = (           - H_d_2 * (h_hat_d(1) + h_hat_d(0) * V_u_1) ) / MAX(HSq,SqrtTiny)
    h_hat_d_H1(3) = (           - H_d_3 * (h_hat_d(1) + h_hat_d(0) * V_u_1) ) / MAX(HSq,SqrtTiny)

    h_hat_d_H2(1) = (           - H_d_1 * (h_hat_d(2) + h_hat_d(0) * V_u_2) ) / MAX(HSq,SqrtTiny)
    h_hat_d_H2(2) = ( SQRT(HSq) - H_d_2 * (h_hat_d(2) + h_hat_d(0) * V_u_2) ) / MAX(HSq,SqrtTiny)
    h_hat_d_H2(3) = (           - H_d_3 * (h_hat_d(2) + h_hat_d(0) * V_u_2) ) / MAX(HSq,SqrtTiny)

    h_hat_d_H3(1) = (           - H_d_1 * (h_hat_d(3) + h_hat_d(0) * V_u_3) ) / MAX(HSq,SqrtTiny)
    h_hat_d_H3(2) = (           - H_d_2 * (h_hat_d(3) + h_hat_d(0) * V_u_3) ) / MAX(HSq,SqrtTiny)
    h_hat_d_H3(3) = ( SQRT(HSq) - H_d_3 * (h_hat_d(3) + h_hat_d(0) * V_u_3) ) / MAX(HSq,SqrtTiny)

    vMagSq = V_u_1 * Gm_dd_11 * V_u_1 &
           + V_u_2 * Gm_dd_22 * V_u_2 &
           + V_u_3 * Gm_dd_33 * V_u_3
    W = One / SQRT ( One - vMagSq )

    ! --- Fluid four-velocity i=1,2,3 components ===
    u_d(1) = W * Gm_dd_11 * V_u_1
    u_d(2) = W * Gm_dd_22 * V_u_2
    u_d(3) = W * Gm_dd_33 * V_u_3

    Gm_dd = Zero
    Gm_dd(1,1) = Gm_dd_11
    Gm_dd(2,2) = Gm_dd_22
    Gm_dd(3,3) = Gm_dd_33

    ! --- Projection operator for fluid four-velocity
    DO k = 1, 3
      DO i = 1, 3

        h_dd(i,k) = Gm_dd(i,k) + u_d(i) * u_d(k)

      END DO
    END DO

    ! --- Partial derivatives of the Eddington Factor ---
    EF_J  = 2.0_DP / 15.0_DP * &
            ( - 6.0_DP * FF + 3.0_DP * FF**2 - 12.0_DP * FF**3 ) * &
            FF / J
    EF_H1 = 2.0_DP / 15.0_DP * &
            (   6.0_DP * FF - 3.0_DP * FF**2 + 12.0_DP * FF**3 ) * &
            ( h_hat_d(1) + h_hat_d(0) * V_u_1 ) / J
    EF_H2 = 2.0_DP / 15.0_DP * &
            (   6.0_DP * FF - 3.0_DP * FF**2 + 12.0_DP * FF**3 ) * &
            ( h_hat_d(2) + h_hat_d(0) * V_u_2 ) / J
    EF_H3 = 2.0_DP / 15.0_DP * &
            (   6.0_DP * FF - 3.0_DP * FF**2 + 12.0_DP * FF**3 ) * &
            ( h_hat_d(3) + h_hat_d(0) * V_u_3 ) / J

    ! print *, 'EF_J ', EF_J
    ! print *, 'EF_H1', EF_H1
    ! print *, 'EF_H2', EF_H2
    ! print *, 'EF_H3', EF_H3

    ! --- Partial derivative of \mathcal{K}_{ij} wrt to \mathcal{J} ---
    DO k = 1, 3
    DO i = 1, 3

      K_dd_J(i,k) = k_dd(i,k) + &
                    J / Two * ( 3.0_DP * h_hat_d(i) * h_hat_d(k) - h_dd(i,k) ) * EF_J

    END DO
    END DO

    ! print *, 'K_dd_J', K_dd_J

    ! --- Partial derivative of \mathcal{K}_{ij} wrt to \mathcal{H}_1 ---
    DO k = 1, 3
    DO i = 1, 3

      K_dd_H1(i,k) = J / Two * ( EF_H1 * (3.0_DP * h_hat_d(i) * h_hat_d(k) - h_dd(i,k)) + &
                                 (3.0_DP * EF - One) * &
                                 (h_hat_d_H1(i) * h_hat_d(k) + h_hat_d(i) * h_hat_d_H1(k)) )

    END DO
    END DO

    ! print *, 'K_dd_H1', K_dd_H1

    ! --- Partial derivative of \mathcal{K}_{ij} wrt to \mathcal{H}_2 ---
    DO k = 1, 3
    DO i = 1, 3

      K_dd_H2(i,k) = J / Two * ( EF_H2 * (3.0_DP * h_hat_d(i) * h_hat_d(k) - h_dd(i,k)) + &
                                 (3.0_DP * EF - One) * &
                                 (h_hat_d_H2(i) * h_hat_d(k) + h_hat_d(i) * h_hat_d_H2(k)) )

    END DO
    END DO

    ! --- Partial derivative of \mathcal{K}_{ij} wrt to \mathcal{H}_3 ---
    DO k = 1, 3
    DO i = 1, 3

      K_dd_H3(i,k) = J / Two * ( EF_H3 * (3.0_DP * h_hat_d(i) * h_hat_d(k) - h_dd(i,k)) + &
                                 (3.0_DP * EF - One) * &
                                 (h_hat_d_H3(i) * h_hat_d(k) + h_hat_d(i) * h_hat_d_H3(k)) )

    END DO
    END DO

    ! --- Set components of Jacobian matrix for moment conversion ---

    vDotK1_J = V_u_1 * K_dd_J(1,1) + V_u_2 * K_dd_J(1,2) + V_u_3 * K_dd_J(1,3)
    vDotK2_J = V_u_1 * K_dd_J(2,1) + V_u_2 * K_dd_J(2,2) + V_u_3 * K_dd_J(2,3)
    vDotK3_J = V_u_1 * K_dd_J(3,1) + V_u_2 * K_dd_J(3,2) + V_u_3 * K_dd_J(3,3)

    vDotK1_H1 = V_u_1 * K_dd_H1(1,1) + V_u_2 * K_dd_H1(1,2) + V_u_3 * K_dd_H1(1,3)
    vDotK2_H1 = V_u_1 * K_dd_H1(2,1) + V_u_2 * K_dd_H1(2,2) + V_u_3 * K_dd_H1(2,3)
    vDotK3_H1 = V_u_1 * K_dd_H1(3,1) + V_u_2 * K_dd_H1(3,2) + V_u_3 * K_dd_H1(3,3)

    vDotK1_H2 = V_u_1 * K_dd_H2(1,1) + V_u_2 * K_dd_H2(1,2) + V_u_3 * K_dd_H2(1,3)
    vDotK2_H2 = V_u_1 * K_dd_H2(2,1) + V_u_2 * K_dd_H2(2,2) + V_u_3 * K_dd_H2(2,3)
    vDotK3_H2 = V_u_1 * K_dd_H2(3,1) + V_u_2 * K_dd_H2(3,2) + V_u_3 * K_dd_H2(3,3)

    vDotK1_H3 = V_u_1 * K_dd_H3(1,1) + V_u_2 * K_dd_H3(1,2) + V_u_3 * K_dd_H3(1,3)
    vDotK2_H3 = V_u_1 * K_dd_H3(2,1) + V_u_2 * K_dd_H3(2,2) + V_u_3 * K_dd_H3(2,3)
    vDotK3_H3 = V_u_1 * K_dd_H3(3,1) + V_u_2 * K_dd_H3(3,2) + V_u_3 * K_dd_H3(3,3)

    Jacobian(1,1) = W
    Jacobian(1,2) = V_u_1
    Jacobian(1,3) = V_u_2
    Jacobian(1,4) = V_u_3
    Jacobian(2,1) = vDotK1_J
    Jacobian(2,2) = vDotK1_H1 + W
    Jacobian(2,3) = vDotK1_H2
    Jacobian(2,4) = vDotK1_H3
    Jacobian(3,1) = vDotK2_J
    Jacobian(3,2) = vDotK2_H1
    Jacobian(3,3) = vDotK2_H2 + W
    Jacobian(3,4) = vDotK2_H3
    Jacobian(4,1) = vDotK3_J
    Jacobian(4,2) = vDotK3_H1
    Jacobian(4,3) = vDotK3_H2
    Jacobian(4,4) = vDotK3_H3 + W

    ! print *, 'Matrix', Jacobian

    ! --- Set components of inverse Jacobian matrix for moment conversion ---

    detinv = &
      1/(Jacobian(1,1)*(Jacobian(2,2)*(Jacobian(3,3)*Jacobian(4,4)-Jacobian(3,4)*Jacobian(4,3))+&
                        Jacobian(2,3)*(Jacobian(3,4)*Jacobian(4,2)-Jacobian(3,2)*Jacobian(4,4))+&
                        Jacobian(2,4)*(Jacobian(3,2)*Jacobian(4,3)-Jacobian(3,3)*Jacobian(4,2)))&
       - Jacobian(1,2)*(Jacobian(2,1)*(Jacobian(3,3)*Jacobian(4,4)-Jacobian(3,4)*Jacobian(4,3))+&
                        Jacobian(2,3)*(Jacobian(3,4)*Jacobian(4,1)-Jacobian(3,1)*Jacobian(4,4))+&
                        Jacobian(2,4)*(Jacobian(3,1)*Jacobian(4,3)-Jacobian(3,3)*Jacobian(4,1)))&
       + Jacobian(1,3)*(Jacobian(2,1)*(Jacobian(3,2)*Jacobian(4,4)-Jacobian(3,4)*Jacobian(4,2))+&
                        Jacobian(2,2)*(Jacobian(3,4)*Jacobian(4,1)-Jacobian(3,1)*Jacobian(4,4))+&
                        Jacobian(2,4)*(Jacobian(3,1)*Jacobian(4,2)-Jacobian(3,2)*Jacobian(4,1)))&
       - Jacobian(1,4)*(Jacobian(2,1)*(Jacobian(3,2)*Jacobian(4,3)-Jacobian(3,3)*Jacobian(4,2))+&
                        Jacobian(2,2)*(Jacobian(3,3)*Jacobian(4,1)-Jacobian(3,1)*Jacobian(4,3))+&
                        Jacobian(2,3)*(Jacobian(3,1)*Jacobian(4,2)-Jacobian(3,2)*Jacobian(4,1))))

    MomentConversion_Jacobian_Inverse(1,1) = &
      detinv*(Jacobian(2,2)*(Jacobian(3,3)*Jacobian(4,4)-Jacobian(3,4)*Jacobian(4,3))+&
              Jacobian(2,3)*(Jacobian(3,4)*Jacobian(4,2)-Jacobian(3,2)*Jacobian(4,4))+&
              Jacobian(2,4)*(Jacobian(3,2)*Jacobian(4,3)-Jacobian(3,3)*Jacobian(4,2)))
    MomentConversion_Jacobian_Inverse(2,1) = &
      detinv*(Jacobian(2,1)*(Jacobian(3,4)*Jacobian(4,3)-Jacobian(3,3)*Jacobian(4,4))+&
              Jacobian(2,3)*(Jacobian(3,1)*Jacobian(4,4)-Jacobian(3,4)*Jacobian(4,1))+&
              Jacobian(2,4)*(Jacobian(3,3)*Jacobian(4,1)-Jacobian(3,1)*Jacobian(4,3)))
    MomentConversion_Jacobian_Inverse(3,1) = &
      detinv*(Jacobian(2,1)*(Jacobian(3,2)*Jacobian(4,4)-Jacobian(3,4)*Jacobian(4,2))+&
              Jacobian(2,2)*(Jacobian(3,4)*Jacobian(4,1)-Jacobian(3,1)*Jacobian(4,4))+&
              Jacobian(2,4)*(Jacobian(3,1)*Jacobian(4,2)-Jacobian(3,2)*Jacobian(4,1)))
    MomentConversion_Jacobian_Inverse(4,1) = &
      detinv*(Jacobian(2,1)*(Jacobian(3,3)*Jacobian(4,2)-Jacobian(3,2)*Jacobian(4,3))+&
              Jacobian(2,2)*(Jacobian(3,1)*Jacobian(4,3)-Jacobian(3,3)*Jacobian(4,1))+&
              Jacobian(2,3)*(Jacobian(3,2)*Jacobian(4,1)-Jacobian(3,1)*Jacobian(4,2)))
    MomentConversion_Jacobian_Inverse(1,2) = &
      detinv*(Jacobian(1,2)*(Jacobian(3,4)*Jacobian(4,3)-Jacobian(3,3)*Jacobian(4,4))+&
              Jacobian(1,3)*(Jacobian(3,2)*Jacobian(4,4)-Jacobian(3,4)*Jacobian(4,2))+&
              Jacobian(1,4)*(Jacobian(3,3)*Jacobian(4,2)-Jacobian(3,2)*Jacobian(4,3)))
    MomentConversion_Jacobian_Inverse(2,2) = &
      detinv*(Jacobian(1,1)*(Jacobian(3,3)*Jacobian(4,4)-Jacobian(3,4)*Jacobian(4,3))+&
              Jacobian(1,3)*(Jacobian(3,4)*Jacobian(4,1)-Jacobian(3,1)*Jacobian(4,4))+&
              Jacobian(1,4)*(Jacobian(3,1)*Jacobian(4,3)-Jacobian(3,3)*Jacobian(4,1)))
    MomentConversion_Jacobian_Inverse(3,2) = &
      detinv*(Jacobian(1,1)*(Jacobian(3,4)*Jacobian(4,2)-Jacobian(3,2)*Jacobian(4,4))+&
              Jacobian(1,2)*(Jacobian(3,1)*Jacobian(4,4)-Jacobian(3,4)*Jacobian(4,1))+&
              Jacobian(1,4)*(Jacobian(3,2)*Jacobian(4,1)-Jacobian(3,1)*Jacobian(4,2)))
    MomentConversion_Jacobian_Inverse(4,2) = &
      detinv*(Jacobian(1,1)*(Jacobian(3,2)*Jacobian(4,3)-Jacobian(3,3)*Jacobian(4,2))+&
              Jacobian(1,2)*(Jacobian(3,3)*Jacobian(4,1)-Jacobian(3,1)*Jacobian(4,3))+&
              Jacobian(1,3)*(Jacobian(3,1)*Jacobian(4,2)-Jacobian(3,2)*Jacobian(4,1)))
    MomentConversion_Jacobian_Inverse(1,3) = &
      detinv*(Jacobian(1,2)*(Jacobian(2,3)*Jacobian(4,4)-Jacobian(2,4)*Jacobian(4,3))+&
              Jacobian(1,3)*(Jacobian(2,4)*Jacobian(4,2)-Jacobian(2,2)*Jacobian(4,4))+&
              Jacobian(1,4)*(Jacobian(2,2)*Jacobian(4,3)-Jacobian(2,3)*Jacobian(4,2)))
    MomentConversion_Jacobian_Inverse(2,3) = &
      detinv*(Jacobian(1,1)*(Jacobian(2,4)*Jacobian(4,3)-Jacobian(2,3)*Jacobian(4,4))+&
              Jacobian(1,3)*(Jacobian(2,1)*Jacobian(4,4)-Jacobian(2,4)*Jacobian(4,1))+&
              Jacobian(1,4)*(Jacobian(2,3)*Jacobian(4,1)-Jacobian(2,1)*Jacobian(4,3)))
    MomentConversion_Jacobian_Inverse(3,3) = &
      detinv*(Jacobian(1,1)*(Jacobian(2,2)*Jacobian(4,4)-Jacobian(2,4)*Jacobian(4,2))+&
              Jacobian(1,2)*(Jacobian(2,4)*Jacobian(4,1)-Jacobian(2,1)*Jacobian(4,4))+&
              Jacobian(1,4)*(Jacobian(2,1)*Jacobian(4,2)-Jacobian(2,2)*Jacobian(4,1)))
    MomentConversion_Jacobian_Inverse(4,3) = &
      detinv*(Jacobian(1,1)*(Jacobian(2,3)*Jacobian(4,2)-Jacobian(2,2)*Jacobian(4,3))+&
              Jacobian(1,2)*(Jacobian(2,1)*Jacobian(4,3)-Jacobian(2,3)*Jacobian(4,1))+&
              Jacobian(1,3)*(Jacobian(2,2)*Jacobian(4,1)-Jacobian(2,1)*Jacobian(4,2)))
    MomentConversion_Jacobian_Inverse(1,4) = &
      detinv*(Jacobian(1,2)*(Jacobian(2,4)*Jacobian(3,3)-Jacobian(2,3)*Jacobian(3,4))+&
              Jacobian(1,3)*(Jacobian(2,2)*Jacobian(3,4)-Jacobian(2,4)*Jacobian(3,2))+&
              Jacobian(1,4)*(Jacobian(2,3)*Jacobian(3,2)-Jacobian(2,2)*Jacobian(3,3)))
    MomentConversion_Jacobian_Inverse(2,4) = &
      detinv*(Jacobian(1,1)*(Jacobian(2,3)*Jacobian(3,4)-Jacobian(2,4)*Jacobian(3,3))+&
              Jacobian(1,3)*(Jacobian(2,4)*Jacobian(3,1)-Jacobian(2,1)*Jacobian(3,4))+&
              Jacobian(1,4)*(Jacobian(2,1)*Jacobian(3,3)-Jacobian(2,3)*Jacobian(3,1)))
    MomentConversion_Jacobian_Inverse(3,4) = &
      detinv*(Jacobian(1,1)*(Jacobian(2,4)*Jacobian(3,2)-Jacobian(2,2)*Jacobian(3,4))+&
              Jacobian(1,2)*(Jacobian(2,1)*Jacobian(3,4)-Jacobian(2,4)*Jacobian(3,1))+&
              Jacobian(1,4)*(Jacobian(2,2)*Jacobian(3,1)-Jacobian(2,1)*Jacobian(3,2)))
    MomentConversion_Jacobian_Inverse(4,4) = &
      detinv*(Jacobian(1,1)*(Jacobian(2,2)*Jacobian(3,3)-Jacobian(2,3)*Jacobian(3,2))+&
              Jacobian(1,2)*(Jacobian(2,3)*Jacobian(3,1)-Jacobian(2,1)*Jacobian(3,3))+&
              Jacobian(1,3)*(Jacobian(2,1)*Jacobian(3,2)-Jacobian(2,2)*Jacobian(3,1)))
  
    RETURN
  END FUNCTION MomentConversion_Jacobian_Inverse

  SUBROUTINE ComputeHatMomentsFromConserved &
    ( E, F_d_1, F_d_2, F_d_3, E_hat, F_hat_d_1, F_hat_d_2, F_hat_d_3, &
      V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    ! --- Input/Output variables ---
    REAL(DP), INTENT(in) :: E, F_d_1, F_d_2, F_d_3
    REAL(DP), INTENT(out) :: E_hat, F_hat_d_1, F_hat_d_2, F_hat_d_3
    REAL(DP), INTENT(in) :: V_u_1, V_u_2, V_u_3
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    ! --- Local variables ---
    REAL(DP) :: W, vMagSq, vDotF
    
    vMagSq = V_u_1 * Gm_dd_11 * V_u_1 &
           + V_u_2 * Gm_dd_22 * V_u_2 &
           + V_u_3 * Gm_dd_33 * V_u_3
    W = One / SQRT ( One - vMagSq )
    vDotF = V_u_1 * F_d_1 + V_u_2 * F_d_2 + V_u_3 * F_d_3

    E_hat = W * ( E - vDotF )
    F_hat_d_1 = F_d_1 - W * Gm_dd_11 * V_u_1 * E_hat
    F_hat_d_2 = F_d_2 - W * Gm_dd_22 * V_u_2 * E_hat
    F_hat_d_3 = F_d_3 - W * Gm_dd_33 * V_u_3 * E_hat

  END SUBROUTINE ComputeHatMomentsFromConserved

  SUBROUTINE ComputeHatMomentsFromPrimitive &
    ( J, H_d_1, H_d_2, H_d_3, E_hat, F_hat_d_1, F_hat_d_2, F_hat_d_3, &
      V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33)
    
    ! --- Input/Output variables ---
    REAL(DP), INTENT(in) :: J, H_d_1, H_d_2, H_d_3
    REAL(DP), INTENT(out) :: E_hat, F_hat_d_1, F_hat_d_2, F_hat_d_3
    REAL(DP), INTENT(in) :: V_u_1, V_u_2, V_u_3
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    ! --- Local variables ---
    REAL(DP) :: k_dd(3,3)
    REAL(DP) :: vMagSq, W, vDotH, vDotk1, vDotk2, vDotk3

    k_dd = EddingtonTensorComponents_dd &
      ( J, H_d_1, H_d_2, H_d_3, V_u_1, V_u_2, V_u_3, &
        Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    vMagSq = V_u_1 * Gm_dd_11 * V_u_1 &
           + V_u_2 * Gm_dd_22 * V_u_2 &
           + V_u_3 * Gm_dd_33 * V_u_3
    W = One / SQRT ( One - vMagSq )

    vDotH = V_u_1 * H_d_1 + V_u_2 * H_d_2 + V_u_3 * H_d_3

    vDotk1 = V_u_1 * k_dd(1,1) + V_u_2 * k_dd(1,2) + V_u_3 * k_dd(1,3)
    vDotk2 = V_u_1 * k_dd(2,1) + V_u_2 * k_dd(2,2) + V_u_3 * k_dd(2,3)
    vDotk3 = V_u_1 * k_dd(3,1) + V_u_2 * k_dd(3,2) + V_u_3 * k_dd(3,3)

    E_hat = W * J + vDotH
    F_hat_d_1 = W * H_d_1 + vDotk1 * J
    F_hat_d_2 = W * H_d_2 + vDotk2 * J
    F_hat_d_3 = W * H_d_3 + vDotk3 * J

  END SUBROUTINE ComputeHatMomentsFromPrimitive

  FUNCTION EddingtonTensorComponents_dd &
    ( J, H_d_1, H_d_2, H_d_3, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    ! --- Input/Output variables ---
    REAL(DP), INTENT(in) :: J, H_d_1, H_d_2, H_d_3
    REAL(DP), INTENT(in) :: V_u_1, V_u_2, V_u_3
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP) :: EddingtonTensorComponents_dd(3,3)

    ! --- Local variables ---
    INTEGER :: i, k
    REAL(DP) :: FF, EF, a, b, W, vMagSq
    REAL(DP) :: h_hat_d(3), u_d(3), Gm_dd(3,3) 

    ! FF = FluxFactor_Relativistic &
    !   ( J, H_d_1, H_d_2, H_d_3, &
    !     Gm_dd_11, Gm_dd_22, Gm_dd_33, &
    !     One, Zero, Zero, Zero, &
    !     V_u_1, V_u_2, V_u_3 )

    ! IF (J < SQRT( -(V_u_1 * H_d_1 + V_u_2 * H_d_2 + V_u_3 * H_d_3)**2 + H_d_1**2 + H_d_2**2 + H_d_3**2 ) ) THEN
    !   ! print *, "V_u_1 = ", V_u_1
    !   ! print *, "V_u_2 = ", V_u_2
    !   ! print *, "V_u_3 = ", V_u_3
    !   print *, "J     = ", J
    !   print *, "H_d_1 = ", H_d_1
    !   ! print *, "H_d_2 = ", H_d_2
    !   ! print *, "H_d_3 = ", H_d_3
    !   Write(*,*)
      
    ! END IF

    FF = MIN( MAX( SQRT( -(V_u_1 * H_d_1 + V_u_2 * H_d_2 + V_u_3 * H_d_3)**2 &
                         + H_d_1**2 + H_d_2**2 + H_d_3**2 ) &
                         / MAX( ABS( J ), SqrtTiny ), &
                   SqrtTiny ), &
              One )
    ! FF = One

    ! print *, "FF = ", FF

    IF ( FF <= SqrtTiny ) THEN
      EF = Third
      a = Third
      b = Zero
    ELSE
      EF = EddingtonFactor( J, FF )
      a = Half * ( One - EF )
      b = Half * ( Three * EF - One )
    END IF

    ! EF = EddingtonFactor( J, FF )

    ! a = Half * ( One - EF )
    ! b= Half * ( Three * EF - One )

    h_hat_d(1) = H_d_1 / ( MAX( FF * ABS(J), SqrtTiny ) ) ! Was running in to division by zero.
    h_hat_d(2) = H_d_2 / ( MAX( FF * ABS(J), SqrtTiny ) ) ! Is this a good fix?
    h_hat_d(3) = H_d_3 / ( MAX( FF * ABS(J), SqrtTiny ) )

    vMagSq = V_u_1 * Gm_dd_11 * V_u_1 &
           + V_u_2 * Gm_dd_22 * V_u_2 &
           + V_u_3 * Gm_dd_33 * V_u_3
    W = One / SQRT ( One - vMagSq )

    u_d(1) = W * Gm_dd_11 * V_u_1
    u_d(2) = W * Gm_dd_22 * V_u_2
    u_d(3) = W * Gm_dd_33 * V_u_3

    Gm_dd = Zero
    Gm_dd(1,1) = Gm_dd_11
    Gm_dd(2,2) = Gm_dd_22
    Gm_dd(3,3) = Gm_dd_33

    DO k = 1, 3
      DO i = 1, 3

        EddingtonTensorComponents_dd(i,k) &
          = a * ( Gm_dd(i,k) + u_d(i) * u_d(k) ) + b * h_hat_d(i) * h_hat_d(k)

      END DO
    END DO

    RETURN
  END FUNCTION EddingtonTensorComponents_dd

  SUBROUTINE HeatFluxTensorComponents_uuu &
    ( J, H_u_1, H_u_2, H_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      V_u_1, V_u_2, V_u_3, l_uuu_munurho )

    ! Loops to shorten code and get rid of symmetric variables?

    ! --- Input/Output variables ---
    REAL(DP), INTENT(in) :: J, H_u_1, H_u_2, H_u_3
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), INTENT(in) :: V_u_1, V_u_2, V_u_3
    REAL(DP), INTENT(out) :: l_uuu_munurho(0:3,0:3,0:3)

    ! --- Local variables ---
    REAL(DP) :: FF, HF
    REAL(DP) :: W, vMagSq, H_u_0
    REAL(DP) :: h_hat_u_1, h_hat_u_2, h_hat_u_3
    REAL(DP) :: h_uu_11, h_uu_12, h_uu_13, h_uu_22, h_uu_23, h_uu_33
    REAL(DP) :: a, b
    REAL(DP) :: l_uuu_111, l_uuu_112, l_uuu_113, &
                l_uuu_121, l_uuu_122, l_uuu_123, &
                l_uuu_131, l_uuu_132, l_uuu_133, &
                l_uuu_211, l_uuu_212, l_uuu_213, &
                l_uuu_221, l_uuu_222, l_uuu_223, &
                l_uuu_231, l_uuu_232, l_uuu_233, &
                l_uuu_311, l_uuu_312, l_uuu_313, &
                l_uuu_321, l_uuu_322, l_uuu_323, &
                l_uuu_331, l_uuu_332, l_uuu_333
    REAL(DP) :: l_uuu_011, l_uuu_012, l_uuu_013, &
                l_uuu_021, l_uuu_022, l_uuu_023, &
                l_uuu_031, l_uuu_032, l_uuu_033, &
                l_uuu_101, l_uuu_102, l_uuu_103, &
                l_uuu_201, l_uuu_202, l_uuu_203, &
                l_uuu_301, l_uuu_302, l_uuu_303, &
                l_uuu_110, l_uuu_120, l_uuu_130, &
                l_uuu_210, l_uuu_220, l_uuu_230, &
                l_uuu_310, l_uuu_320, l_uuu_330, &
                l_uuu_001, l_uuu_002, l_uuu_003, &
                l_uuu_100, l_uuu_200, l_uuu_300, &
                l_uuu_010, l_uuu_020, l_uuu_030, &
                l_uuu_000


    ! FF = FluxFactor_Relativistic &
    !   ( J, H_u_1, H_u_2, H_u_3, &
    !       Gm_dd_11, Gm_dd_22, Gm_dd_33, &
    !       One, Zero, Zero, Zero, &
    !       V_u_1, V_u_2, V_u_3 )
    ! FF = SQRT( -(V_u_1 * H_u_1 + V_u_2 * H_u_2 + V_u_3 * H_u_3)**2 &
    !            + H_u_1**2 + H_u_2**2 + H_u_3**2 ) / J
    FF = MIN( MAX( SQRT( -(V_u_1 * H_u_1 + V_u_2 * H_u_2 + V_u_3 * H_u_3)**2 &
                         + H_u_1**2 + H_u_2**2 + H_u_3**2 ) &
                         / MAX( J, SqrtTiny ), &
                   SqrtTiny ), &
              One )
    HF = HeatFluxFactor(J, FF)

    vMagSq = V_u_1 * Gm_dd_11 * V_u_1 &
           + V_u_2 * Gm_dd_22 * V_u_2 &
           + V_u_3 * Gm_dd_33 * V_u_3
    W = One / SQRT ( One - vMagSq )

    h_hat_u_1 = H_u_1 / ( MAX( FF*J, SqrtTiny ) )!( FF * MAX( ABS(J), SqrtTiny ) )
    h_hat_u_2 = H_u_2 / ( MAX( FF*J, SqrtTiny ) )!( FF * MAX( ABS(J), SqrtTiny ) )
    h_hat_u_3 = H_u_3 / ( MAX( FF*J, SqrtTiny ) )!( FF * MAX( ABS(J), SqrtTiny ) )

    h_uu_11 = Gm_dd_11 + W**2 * V_u_1 * V_u_1
    h_uu_12 = W**2 * V_u_1 * V_u_2
    h_uu_13 = W**2 * V_u_1 * V_u_3
    h_uu_22 = Gm_dd_22 + W**2 * V_u_2 * V_u_2
    h_uu_23 = W**2 * V_u_2 * V_u_3
    h_uu_33 = Gm_dd_33 + W**2 * V_u_3 * V_u_3

    a = Half * (FF - HF)
    b = Half * (Five * HF - Three * FF)

    ! --- ijk components ---

    l_uuu_111 = a * (h_hat_u_1 * h_uu_11 + h_hat_u_1 * h_uu_11 + h_hat_u_1 * h_uu_11) &
                + b * (h_hat_u_1 * h_hat_u_1 * h_hat_u_1)
    l_uuu_112 = a * (h_hat_u_1 * h_uu_12 + h_hat_u_1 * h_uu_12 + h_hat_u_2 * h_uu_11) &
                + b * (h_hat_u_1 * h_hat_u_1 * h_hat_u_2)
    l_uuu_113 = a * (h_hat_u_1 * h_uu_13 + h_hat_u_1 * h_uu_13 + h_hat_u_3 * h_uu_11) &
                + b * (h_hat_u_1 * h_hat_u_1 * h_hat_u_3)
    l_uuu_121 = l_uuu_112
    l_uuu_122 = a * (h_hat_u_1 * h_uu_22 + h_hat_u_2 * h_uu_12 + h_hat_u_2 * h_uu_12) &
                + b * (h_hat_u_1 * h_hat_u_2 * h_hat_u_2)
    l_uuu_123 = a * (h_hat_u_1 * h_uu_23 + h_hat_u_2 * h_uu_13 + h_hat_u_3 * h_uu_12) &
                + b * (h_hat_u_1 * h_hat_u_2 * h_hat_u_3)
    l_uuu_131 = l_uuu_113
    l_uuu_132 = l_uuu_123
    l_uuu_133 = a * (h_hat_u_1 * h_uu_33 + h_hat_u_3 * h_uu_13 + h_hat_u_3 * h_uu_13) &
                + b * (h_hat_u_1 * h_hat_u_3 * h_hat_u_3)
    l_uuu_211 = l_uuu_112
    l_uuu_212 = l_uuu_122
    l_uuu_213 = l_uuu_123
    l_uuu_221 = l_uuu_122
    l_uuu_222 = a * (h_hat_u_2 * h_uu_22 + h_hat_u_2 * h_uu_22 + h_hat_u_2 * h_uu_22) &
                + b * (h_hat_u_2 * h_hat_u_2 * h_hat_u_2)
    l_uuu_223 = a * (h_hat_u_2 * h_uu_23 + h_hat_u_2 * h_uu_23 + h_hat_u_3 * h_uu_22) &
                + b * (h_hat_u_2 * h_hat_u_2 * h_hat_u_3)
    l_uuu_231 = l_uuu_123
    l_uuu_232 = l_uuu_223
    l_uuu_233 = a * (h_hat_u_2 * h_uu_33 + h_hat_u_3 * h_uu_23 + h_hat_u_3 * h_uu_23) &
                + b * (h_hat_u_2 * h_hat_u_3 * h_hat_u_3)
    l_uuu_311 = l_uuu_113
    l_uuu_312 = l_uuu_123
    l_uuu_313 = l_uuu_133
    l_uuu_321 = l_uuu_123
    l_uuu_322 = l_uuu_223
    l_uuu_323 = l_uuu_233
    l_uuu_331 = l_uuu_133
    l_uuu_332 = l_uuu_233
    l_uuu_333 = a * (h_hat_u_3 * h_uu_33 + h_hat_u_3 * h_uu_33 + h_hat_u_3 * h_uu_33) &
                + b * (h_hat_u_3 * h_hat_u_3 * h_hat_u_3)

    ! --- 0th components ---

    l_uuu_011 = V_u_1 * l_uuu_111 + V_u_2 * l_uuu_211 + V_u_3 * l_uuu_311
    l_uuu_012 = V_u_1 * l_uuu_112 + V_u_2 * l_uuu_212 + V_u_3 * l_uuu_312
    l_uuu_013 = V_u_1 * l_uuu_113 + V_u_2 * l_uuu_213 + V_u_3 * l_uuu_313
    l_uuu_021 = l_uuu_012
    l_uuu_022 = V_u_1 * l_uuu_122 + V_u_2 * l_uuu_222 + V_u_3 * l_uuu_322
    l_uuu_023 = V_u_1 * l_uuu_123 + V_u_2 * l_uuu_223 + V_u_3 * l_uuu_323
    l_uuu_031 = l_uuu_013
    l_uuu_032 = l_uuu_023
    l_uuu_033 = V_u_1 * l_uuu_133 + V_u_2 * l_uuu_233 + V_u_3 * l_uuu_333
    l_uuu_101 = l_uuu_011
    l_uuu_102 = l_uuu_012
    l_uuu_103 = l_uuu_013
    l_uuu_201 = l_uuu_012
    l_uuu_202 = l_uuu_022
    l_uuu_203 = l_uuu_023
    l_uuu_301 = l_uuu_013
    l_uuu_302 = l_uuu_023
    l_uuu_303 = l_uuu_033
    l_uuu_110 = l_uuu_011
    l_uuu_120 = l_uuu_012
    l_uuu_130 = l_uuu_013
    l_uuu_210 = l_uuu_012
    l_uuu_220 = l_uuu_022
    l_uuu_230 = l_uuu_023
    l_uuu_310 = l_uuu_013
    l_uuu_320 = l_uuu_023
    l_uuu_330 = l_uuu_033
    l_uuu_001 = V_u_1 * l_uuu_011 + V_u_2 * l_uuu_021 + V_u_3 * l_uuu_031
    l_uuu_002 = V_u_1 * l_uuu_012 + V_u_2 * l_uuu_022 + V_u_3 * l_uuu_032
    l_uuu_003 = V_u_1 * l_uuu_013 + V_u_2 * l_uuu_023 + V_u_3 * l_uuu_033
    l_uuu_100 = l_uuu_001
    l_uuu_200 = l_uuu_002
    l_uuu_300 = l_uuu_003
    l_uuu_010 = l_uuu_001
    l_uuu_020 = l_uuu_002
    l_uuu_030 = l_uuu_003
    l_uuu_000 = V_u_1 * l_uuu_001 + V_u_2 * l_uuu_002 + V_u_3 * l_uuu_003

    ! --- Array assignment ---

    l_uuu_munurho(0,0,0) = l_uuu_000
    l_uuu_munurho(0,0,1) = l_uuu_001
    l_uuu_munurho(0,0,2) = l_uuu_002
    l_uuu_munurho(0,0,3) = l_uuu_003
    l_uuu_munurho(0,1,0) = l_uuu_010
    l_uuu_munurho(0,1,1) = l_uuu_011
    l_uuu_munurho(0,1,2) = l_uuu_012
    l_uuu_munurho(0,1,3) = l_uuu_013
    l_uuu_munurho(0,2,0) = l_uuu_020
    l_uuu_munurho(0,2,1) = l_uuu_021
    l_uuu_munurho(0,2,2) = l_uuu_022
    l_uuu_munurho(0,2,3) = l_uuu_023
    l_uuu_munurho(0,3,0) = l_uuu_030
    l_uuu_munurho(0,3,1) = l_uuu_031
    l_uuu_munurho(0,3,2) = l_uuu_032
    l_uuu_munurho(0,3,3) = l_uuu_033
    l_uuu_munurho(1,0,0) = l_uuu_100
    l_uuu_munurho(1,0,1) = l_uuu_101
    l_uuu_munurho(1,0,2) = l_uuu_102
    l_uuu_munurho(1,0,3) = l_uuu_103
    l_uuu_munurho(1,1,0) = l_uuu_110
    l_uuu_munurho(1,1,1) = l_uuu_111
    l_uuu_munurho(1,1,2) = l_uuu_112
    l_uuu_munurho(1,1,3) = l_uuu_113
    l_uuu_munurho(1,2,0) = l_uuu_120
    l_uuu_munurho(1,2,1) = l_uuu_121
    l_uuu_munurho(1,2,2) = l_uuu_122
    l_uuu_munurho(1,2,3) = l_uuu_123
    l_uuu_munurho(1,3,0) = l_uuu_130
    l_uuu_munurho(1,3,1) = l_uuu_131
    l_uuu_munurho(1,3,2) = l_uuu_132
    l_uuu_munurho(1,3,3) = l_uuu_133
    l_uuu_munurho(2,0,0) = l_uuu_200
    l_uuu_munurho(2,0,1) = l_uuu_201
    l_uuu_munurho(2,0,2) = l_uuu_202
    l_uuu_munurho(2,0,3) = l_uuu_203
    l_uuu_munurho(2,1,0) = l_uuu_210
    l_uuu_munurho(2,1,1) = l_uuu_211
    l_uuu_munurho(2,1,2) = l_uuu_212
    l_uuu_munurho(2,1,3) = l_uuu_213
    l_uuu_munurho(2,2,0) = l_uuu_220
    l_uuu_munurho(2,2,1) = l_uuu_221
    l_uuu_munurho(2,2,2) = l_uuu_222
    l_uuu_munurho(2,2,3) = l_uuu_223
    l_uuu_munurho(2,3,0) = l_uuu_230
    l_uuu_munurho(2,3,1) = l_uuu_231
    l_uuu_munurho(2,3,2) = l_uuu_232
    l_uuu_munurho(2,3,3) = l_uuu_233
    l_uuu_munurho(3,0,0) = l_uuu_300
    l_uuu_munurho(3,0,1) = l_uuu_301
    l_uuu_munurho(3,0,2) = l_uuu_302
    l_uuu_munurho(3,0,3) = l_uuu_303
    l_uuu_munurho(3,1,0) = l_uuu_310
    l_uuu_munurho(3,1,1) = l_uuu_311
    l_uuu_munurho(3,1,2) = l_uuu_312
    l_uuu_munurho(3,1,3) = l_uuu_313
    l_uuu_munurho(3,2,0) = l_uuu_320
    l_uuu_munurho(3,2,1) = l_uuu_321
    l_uuu_munurho(3,2,2) = l_uuu_322
    l_uuu_munurho(3,2,3) = l_uuu_323
    l_uuu_munurho(3,3,0) = l_uuu_330
    l_uuu_munurho(3,3,1) = l_uuu_331
    l_uuu_munurho(3,3,2) = l_uuu_332
    l_uuu_munurho(3,3,3) = l_uuu_333
  
  END SUBROUTINE HeatFluxTensorComponents_uuu

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

  FUNCTION Flux_X1 &
    ( J, H_u_1, H_u_2, H_u_3, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33)

    ! --- Input/Output variables ---
    REAL(DP) :: Flux_X1(4)
    REAL(DP), INTENT(in) :: J, H_u_1, H_u_2, H_u_3
    REAL(DP), INTENT(in) :: V_u_1, V_u_2, V_u_3
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    ! --- Local variables ---
    REAL(DP) :: vMagSq, W, vDotH
    REAL(DP) :: k_ud(3,3)
    REAL(DP) :: vK_u_1

    vMagSq = V_u_1 * Gm_dd_11 * V_u_1 &
           + V_u_2 * Gm_dd_22 * V_u_2 &
           + V_u_3 * Gm_dd_33 * V_u_3
    W = One / SQRT ( One - vMagSq )
    vDotH = V_u_1 * H_u_1 + V_u_2 * H_u_2 + V_u_3 * H_u_3

    k_ud = EddingtonTensorComponents_dd ( J, H_u_1, H_u_2, H_u_3, &
      V_u_1, V_u_2, V_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33)
    vK_u_1 = V_u_1 * k_ud(1,1) + V_u_2 * k_ud(1,2) + V_u_3 * k_ud(1,3)

    Flux_X1(1) = W * H_u_1 + W * V_u_1 * (W * J + vDotH) + vK_u_1 * J
    Flux_X1(2) = k_ud(1,1) * J + W * (H_u_1 * Gm_dd_11 * V_u_1 + V_u_1 * H_u_1) &
                 + W**2 * V_u_1 * Gm_dd_11 * V_u_1 * J
    Flux_X1(3) = k_ud(1,2) * J + W * (H_u_1 * Gm_dd_22 * V_u_2 + V_u_1 * H_u_2) &
                 + W**2 * V_u_1 * Gm_dd_22 * V_u_2 * J
    Flux_X1(4) = k_ud(1,3) * J + W * (H_u_1 * Gm_dd_33 * V_u_3 + V_u_1 * H_u_3) &
                 + W**2 * V_u_1 * Gm_dd_33 * V_u_3 * J
    
    RETURN
  END FUNCTION Flux_X1

  FUNCTION Flux_X2 ( J, H_u_1, H_u_2, H_u_3, V_u_1, V_u_2, V_u_3, &
    Gm_dd_11, Gm_dd_22, Gm_dd_33)

    ! --- Input/Output variables ---
    REAL(DP) :: Flux_X2(4)
    REAL(DP), INTENT(in) :: J, H_u_1, H_u_2, H_u_3
    REAL(DP), INTENT(in) :: V_u_1, V_u_2, V_u_3
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    ! --- Local variables ---
    REAL(DP) :: vMagSq, W, vDotH
    REAL(DP) :: k_ud(3,3)
    REAL(DP) :: vK_u_2

    vMagSq = V_u_1 * Gm_dd_11 * V_u_1 &
           + V_u_2 * Gm_dd_22 * V_u_2 &
           + V_u_3 * Gm_dd_33 * V_u_3
    W = One / SQRT ( One - vMagSq )
    vDotH = V_u_1 * H_u_1 + V_u_2 * H_u_2 + V_u_3 * H_u_3

    k_ud = EddingtonTensorComponents_dd ( J, H_u_1, H_u_2, H_u_3, &
      V_u_1, V_u_2, V_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33)
    vK_u_2 = V_u_1 * k_ud(2,1) + V_u_2 * k_ud(2,2) + V_u_3 * k_ud(2,3)

    Flux_X2(1) = W * H_u_2 + W * V_u_2 * (W * J + vDotH) + vK_u_2 * J
    Flux_X2(2) = k_ud(2,1) * J + W * (H_u_2 * Gm_dd_11 * V_u_1 + V_u_2 * H_u_1) &
                 + W**2 * V_u_2 * Gm_dd_11 * V_u_1 * J
    Flux_X2(3) = k_ud(2,2) * J + W * (H_u_2 * Gm_dd_22 * V_u_2 + V_u_2 * H_u_2) &
                 + W**2 * V_u_2 * Gm_dd_22 * V_u_2 * J
    Flux_X2(4) = k_ud(2,3) * J + W * (H_u_2 * Gm_dd_33 * V_u_3 + V_u_2 * H_u_3) &
                 + W**2 * V_u_2 * Gm_dd_33 * V_u_3 * J
    
    RETURN
  END FUNCTION Flux_X2

  FUNCTION Flux_X3 ( J, H_u_1, H_u_2, H_u_3, V_u_1, V_u_2, V_u_3, &
    Gm_dd_11, Gm_dd_22, Gm_dd_33)

    ! --- Input/Output variables ---
    REAL(DP) :: Flux_X3(4)
    REAL(DP), INTENT(in) :: J, H_u_1, H_u_2, H_u_3
    REAL(DP), INTENT(in) :: V_u_1, V_u_2, V_u_3
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    ! --- Local variables ---
    REAL(DP) :: vMagSq, W, vDotH
    REAL(DP) :: k_ud(3,3)
    REAL(DP) :: vK_u_3

    vMagSq = V_u_1 * Gm_dd_11 * V_u_1 &
           + V_u_2 * Gm_dd_22 * V_u_2 &
           + V_u_3 * Gm_dd_33 * V_u_3
    W = One / SQRT ( One - vMagSq )
    vDotH = V_u_1 * H_u_1 + V_u_2 * H_u_2 + V_u_3 * H_u_3

    k_ud = EddingtonTensorComponents_dd ( J, H_u_1, H_u_2, H_u_3, &
      V_u_1, V_u_2, V_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33)
    vK_u_3 = V_u_1 * k_ud(3,1) + V_u_2 * k_ud(3,2) + V_u_3 * k_ud(3,3)

    Flux_X3(1) = W * H_u_3 + W * V_u_3 * (W * J + vDotH) + vK_u_3 * J
    Flux_X3(2) = k_ud(3,1) * J + W * (H_u_3 * Gm_dd_11 * V_u_1 + V_u_3 * H_u_1) &
                 + W**2 * V_u_3 * Gm_dd_11 * V_u_1 * J
    Flux_X3(3) = k_ud(3,2) * J + W * (H_u_3 * Gm_dd_22 * V_u_2 + V_u_3 * H_u_2) &
                 + W**2 * V_u_3 * Gm_dd_22 * V_u_2 * J
    Flux_X3(4) = k_ud(3,3) * J + W * (H_u_3 * Gm_dd_33 * V_u_3 + V_u_3 * H_u_3) &
                 + W**2 * V_u_3 * Gm_dd_33 * V_u_3 * J
    
    RETURN
  END FUNCTION Flux_X3

  FUNCTION Flux_E( J, H_d_1, H_d_2, H_d_3, &
    V_u_1, V_u_2, V_u_3, &
    dV_d_0_dX1, dV_d_1_dX1, dV_d_2_dX1, dV_d_3_dX1, &
    dV_d_0_dX2, dV_d_1_dX2, dV_d_2_dX2, dV_d_3_dX2, &
    dV_d_0_dX3, dV_d_1_dX3, dV_d_2_dX3, dV_d_3_dX3, &
    Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    ! --- Input/Output variables ---
    REAL(DP) :: Flux_E(4)
    REAL(DP), INTENT(in) :: J, H_d_1, H_d_2, H_d_3
    REAL(DP), INTENT(in) :: V_u_1, V_u_2, V_u_3
    REAL(DP), INTENT(in) :: dV_d_0_dX1, dV_d_1_dX1, dV_d_2_dX1, dV_d_3_dX1
    REAL(DP), INTENT(in) :: dV_d_0_dX2, dV_d_1_dX2, dV_d_2_dX2, dV_d_3_dX2
    REAL(DP), INTENT(in) :: dV_d_0_dX3, dV_d_1_dX3, dV_d_2_dX3, dV_d_3_dX3
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    ! --- Local variables ---
    REAL(DP) :: vMagSq, W, u_u(0:3)
    REAL(DP) :: vH, H_u(0:3)
    REAL(DP) :: k_uu(0:3,0:3), l_uuu(0:3,0:3,0:3), Q_uuu(0:3,0:3,0:3)
    REAL(DP) :: Jacobian_U(0:3,0:3)
    INTEGER :: mu, nu, rho

    vMagSq = V_u_1 * Gm_dd_11 * V_u_1 &
           + V_u_2 * Gm_dd_22 * V_u_2 &
           + V_u_3 * Gm_dd_33 * V_u_3
    W = One / SQRT ( One - vMagSq )
    u_u(0) = W
    u_u(1) = W * V_u_1
    u_u(2) = W * V_u_2
    u_u(3) = W * V_u_3

    vH = V_u_1 * H_d_1 + V_u_2 * H_d_2 + V_u_3 * H_d_3
    H_u(0) = vH
    H_u(1) = H_d_1
    H_u(2) = H_d_2
    H_u(3) = H_d_3

    k_uu(1:3,1:3) = EddingtonTensorComponents_dd &
      ( J, H_d_1, H_d_2, H_d_3, V_u_1, V_u_2, V_u_3, &
        Gm_dd_11, Gm_dd_22, Gm_dd_33 )
    DO nu = 1,3

      k_uu(0,nu) = V_u_1 * k_uu(1,nu) + V_u_2 * k_uu(2,nu) + V_u_3 * k_uu(3,nu)
      k_uu(nu,0) = k_uu(0,nu)

    END DO
    k_uu(0,0) = V_u_1 * k_uu(1,0) + V_u_2 * k_uu(2,0) + V_u_3 * k_uu(3,0)

    CALL HeatFluxTensorComponents_uuu &
      ( J, H_d_1, H_d_2, H_d_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
        V_u_1, V_u_2, V_u_3, l_uuu )

    Jacobian_U(0,:) = Zero ! Assuming velocity is constan wrt to time
    Jacobian_U(1,:) = [ dV_d_0_dX1, dV_d_1_dX1, dV_d_2_dX1, dV_d_3_dX1 ]
    Jacobian_U(2,:) = [ dV_d_0_dX2, dV_d_1_dX2, dV_d_2_dX2, dV_d_3_dX2 ]
    Jacobian_U(3,:) = [ dV_d_0_dX3, dV_d_1_dX3, dV_d_2_dX3, dV_d_3_dX3 ]

    Flux_E = Zero
    DO rho = 0,3
      DO nu = 0,3
        DO mu = 0,3

          Q_uuu(mu,nu,rho) = J*u_u(mu)*u_u(nu)*u_u(rho)+ &
                             H_u(mu)*u_u(nu)*u_u(rho)+ &
                             H_u(nu)*u_u(mu)*u_u(rho)+ &
                             H_u(rho)*u_u(mu)*u_u(nu)+ &
                             k_uu(mu,nu)*J*u_u(rho)+ &
                             k_uu(mu,rho)*J*u_u(nu)+ &
                             k_uu(nu,rho)*J*u_u(mu)+ &
                             l_uuu(mu,nu,rho)*J

        END DO
        Flux_E(1) = Flux_E(1) + Q_uuu(0,nu,rho) * Jacobian_U(nu,rho)
        Flux_E(2) = Flux_E(2) + Q_uuu(1,nu,rho) * Jacobian_U(nu,rho)
        Flux_E(3) = Flux_E(3) + Q_uuu(2,nu,rho) * Jacobian_U(nu,rho)
        Flux_E(4) = Flux_E(4) + Q_uuu(3,nu,rho) * Jacobian_U(nu,rho)
      END DO
    END DO
    Flux_E = -Flux_E
    RETURN

  END FUNCTION Flux_E

  FUNCTION NumericalFlux_LLF( u_L, u_R, Flux_L, Flux_R, alpha )

    REAL(DP)             :: NumericalFlux_LLF
    REAL(DP), INTENT(in) :: u_L, u_R, flux_L, flux_R, alpha

    NumericalFlux_LLF &
      = Half * ( flux_L + flux_R - alpha * ( u_R - u_L ) )

    RETURN
  END FUNCTION NumericalFlux_LLF

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

  SUBROUTINE ComputeTimeStep_TwoMoment_Realizable &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, U_F, CFL, TimeStep, Verbose_Option )

    INTEGER,  INTENT(in)  :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)  :: &
      GX (1:nDOFX, &
          iZ_B1(2):iZ_E1(2), &
          iZ_B1(3):iZ_E1(3), &
          iZ_B1(4):iZ_E1(4), &
          1:nGF)
    REAL(DP), INTENT(inout)  :: &
      U_F(1:nDOFX, &
          iZ_B1(2):iZ_E1(2), &
          iZ_B1(3):iZ_E1(3), &
          iZ_B1(4):iZ_E1(4), &
          1:nPF)
    REAL(DP), INTENT(in)  :: &
      CFL
    REAL(DP), INTENT(out) :: &
      TimeStep
    LOGICAL, INTENT(in), OPTIONAL :: &
      Verbose_Option

    LOGICAL  :: Verbose
    INTEGER  :: iX_B0(3), iX_E0(3)
    INTEGER  :: iX_B1(3), iX_E1(3)
    INTEGER  :: iX1, iX2, iX3, iNodeX, iE
    INTEGER  :: nQ
    REAL(DP), ALLOCATABLE :: xQ(:), wQ(:)
    REAL(DP) :: TimeStep_X, dt_X(3)
    REAL(DP) :: TimeStep_E, dt_E
    REAL(DP) :: V_d_1, V_d_2, V_d_3, vMag, W
    REAL(DP) :: dE_Min, A, B_i(1:3), C_ij(1:3,1:3), Lambda(3), Alpha_E
    REAL(DP), DIMENSION(nDOFX,0:3, &
                        iZ_B0(2):iZ_E0(2), &
                        iZ_B0(3):iZ_E0(3), &
                        iZ_B0(4):iZ_E0(4)) :: &
      dV_d_dX1, &
      dV_d_dX2, &
      dV_d_dX3

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

      nQ = (nNodes+4)/2+MOD(nNodes+4,2)

      ALLOCATE( xQ(nQ) )
      ALLOCATE( wQ(nQ) )
  
      CALL GetQuadrature( nQ, xQ, wQ, 'Lobatto' )

      dE_Min = HUGE( One ) ! --- Min of dE / E_H
      DO iE = iZ_B0(1), iZ_E0(1)

        dE_Min = MIN( dE_Min, dE(iE) / ( E_C(iE) + Half * dE(iE) ) )

      END DO

      CALL ApplyBoundaryConditions_Euler &
           ( iX_B0, iX_E0, iX_B1, iX_E1, U_F )

      CALL ComputeWeakDerivatives_X1 &
           ( iX_B0, iX_E0, iX_B1, iX_E1, GX, U_F, &
             dV_d_dX1 )

      CALL ComputeWeakDerivatives_X2 &
            ( iX_B0, iX_E0, iX_B1, iX_E1, GX, U_F, &
              dV_d_dX2 )

      CALL ComputeWeakDerivatives_X3 &
            ( iX_B0, iX_E0, iX_B1, iX_E1, GX, U_F, &
              dV_d_dX3 )

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
    
        DO iNodeX = 1, nDOFX
   
          V_d_1 = U_F(iNodeX,iX1,iX2,iX3,iPF_V1)
          V_d_2 = U_F(iNodeX,iX1,iX2,iX3,iPF_V2)
          V_d_3 = U_F(iNodeX,iX1,iX2,iX3,iPF_V3)
  
          vMag = SQRT( V_d_1**2 + V_d_2**2 + V_d_3**2 )
          W = One / SQRT ( One - vMag**2 )
  
          ! --- Time Step from Spatial Divergence ---
  
          dt_X(1) = wQ(nQ) * dX1(iX1) / DBLE(nDims) ! change to division by nDims and convert to REAL
          dt_X(2) = wQ(nQ) * dX2(iX2) / DBLE(nDims)
          dt_X(3) = wQ(nQ) * dX3(iX3) / DBLE(nDims)
  
          TimeStep_X = MIN( TimeStep_X, MINVAL( dt_X ) )
  
          ! --- Quadratic Form Matrix ---
          ! Assuming velocity is constant with respect to time

          A = Zero

          B_i = Half * [ dV_d_dX1(iNodeX,0,iX1,iX2,iX3), &
                         dV_d_dX2(iNodeX,0,iX1,iX2,iX3), &
                         dV_d_dX3(iNodeX,0,iX1,iX2,iX3) ]

          C_ij(:,1) = Half * [ Two * dV_d_dX1(iNodeX,1,iX1,iX2,iX3), &
                                     dV_d_dX2(iNodeX,1,iX1,iX2,iX3)  &
                                   + dV_d_dX1(iNodeX,2,iX1,iX2,iX3), &
                                     dV_d_dX3(iNodeX,1,iX1,iX2,iX3)  &
                                   + dV_d_dX1(iNodeX,3,iX1,iX2,iX3) ]

          C_ij(:,2) = Half * [       dV_d_dX1(iNodeX,2,iX1,iX2,iX3)  &
                                   + dV_d_dX2(iNodeX,1,iX1,iX2,iX3), &
                               Two * dV_d_dX2(iNodeX,2,iX1,iX2,iX3), &
                                     dV_d_dX3(iNodeX,2,iX1,iX2,iX3)  &
                                   + dV_d_dX2(iNodeX,3,iX1,iX2,iX3) ]

          C_ij(:,3) = Half * [       dV_d_dX1(iNodeX,3,iX1,iX2,iX3)  &
                                   + dV_d_dX3(iNodeX,1,iX1,iX2,iX3), &
                                     dV_d_dX2(iNodeX,3,iX1,iX2,iX3)  &
                                   + dV_d_dX3(iNodeX,2,iX1,iX2,iX3), &
                               Two * dV_d_dX3(iNodeX,3,iX1,iX2,iX3) ]

  
          CALL EigenvaluesSymmetric3( C_ij, Lambda )
  
          Alpha_E = MAX( MAXVAL( ABS( Lambda ) ), SqrtTiny )

          Alpha_E = MAXVAL( ABS( Lambda ) ) &
                      + Two * SQRT ( B_i(1)**2 + B_i(2)**2 + B_i(3)**2 ) &
                      + ABS( A )
          Alpha_E = Alpha_E * W**2 * ( One + vMag )**2
  
          ! --- Time Step from Energy Divergence ---
  
          dt_E = W * ( One - vMag ) * dE_Min * wQ(nQ) &
                 / ( MAX(Alpha_E,SqrtTiny) * DBLE(nDims) )
  
          TimeStep_E = MIN( TimeStep_E, dt_E )
  
        END DO
    
      END DO
      END DO
      END DO

      TimeStep = MAX( CFL * MIN( TimeStep_X, TimeStep_E ), SqrtTiny )

      IF( Verbose )THEN
        WRITE(*,'(A8,A7,ES12.6E2,A8,ES12.6E2)') &
          '', 'dt_X = ', TimeStep_X, ' dt_E = ', TimeStep_E
      END IF

      END ASSOCIATE ! dX1, etc.

      DEALLOCATE( xQ )
      DEALLOCATE( wQ )

  END SUBROUTINE ComputeTimeStep_TwoMoment_Realizable

  SUBROUTINE ComputeWeakDerivatives_X1 &
    ( iX_B0, iX_E0, iX_B1, iX_E1, GX, U_F, dU_d_dX1_Out )

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
          1:nPF)
    REAL(DP), INTENT(out) :: &
      dU_d_dX1_Out &
        (1:nDOFX,0:3, &
         iX_B0(1):iX_E0(1), &
         iX_B0(2):iX_E0(2), &
         iX_B0(3):iX_E0(3))

    INTEGER  :: iNodeX
    INTEGER  :: iX1, iX2, iX3, i
    INTEGER  :: iPF, iPF_V
    INTEGER  :: iGF, iGF_h, iGF_Gm_dd
    INTEGER  :: nX(3), nX_X1(3), nK_X, nX1_X
    REAL(DP) :: uV_L(0:3), uV_R(0:3), uV_F(0:3), uV_K
    REAL(DP) :: W_K, W_L, W_R, vMagSq_K, vMagSq_L, vMagSq_R

    ! --- Geometry Fields ---

    REAL(DP) :: &
      GX_K   (nDOFX,nGF, &
              iX_B0(2)  :iX_E0(2)  , &
              iX_B0(3)  :iX_E0(3)  , &
              iX_B0(1)-1:iX_E0(1)+1)
    REAL(DP) :: &
      GX_F   (nDOFX_X1,nGF, &
              iX_B0(2)  :iX_E0(2)  , &
              iX_B0(3)  :iX_E0(3)  , &
              iX_B0(1)  :iX_E0(1)+1)

    ! --- Primitive Fluid Fields ---

    REAL(DP) :: &
      U_F_K(nDOFX,0:3, &
            iX_B0(2)  :iX_E0(2)  , &
            iX_B0(3)  :iX_E0(3)  , &
            iX_B0(1)-1:iX_E0(1)+1)
    REAL(DP) :: &
      U_F_L(nDOFX_X1,0:3, &
            iX_B0(2)  :iX_E0(2)  , &
            iX_B0(3)  :iX_E0(3)  , &
            iX_B0(1)  :iX_E0(1)+1)
    REAL(DP) :: &
      U_F_R(nDOFX_X1,0:3, &
            iX_B0(2)  :iX_E0(2)  , &
            iX_B0(3)  :iX_E0(3)  , &
            iX_B0(1)  :iX_E0(1)+1)

    ! --- Velocities ---

    REAL(DP) :: &
      U_d_X1  (nDOFX_X1,0:3, &
               iX_B0(2)  :iX_E0(2)  , &
               iX_B0(3)  :iX_E0(3)  , &
               iX_B0(1)  :iX_E0(1)+1)
    REAL(DP) :: &
      U_d_K   (nDOFX,0:3, &
               iX_B0(2)  :iX_E0(2)  , &
               iX_B0(3)  :iX_E0(3)  , &
               iX_B0(1)  :iX_E0(1)  )
    REAL(DP) :: &
      dU_d_dX1(nDOFX,0:3, &
               iX_B0(2)  :iX_E0(2)  , &
               iX_B0(3)  :iX_E0(3)  , &
               iX_B0(1)  :iX_E0(1)  )

    ! --- Structures for computing derivative exactly in Transparent Shock ---
    
    REAL(DP) :: &
      dV_d_dX1(nDOFX, &
               iX_B0(2)  :iX_E0(2)  , &
               iX_B0(3)  :iX_E0(3)  , &
               iX_B0(1)  :iX_E0(1)  )

    REAL(DP) :: &
      dW_dX1(nDOFX, &
               iX_B0(2)  :iX_E0(2)  , &
               iX_B0(3)  :iX_E0(3)  , &
               iX_B0(1)  :iX_E0(1)  )

    REAL(DP) :: &
      U_F_P(nDOFX,0:3, &
               iX_B0(2)  :iX_E0(2)  , &
               iX_B0(3)  :iX_E0(3)  , &
               iX_B0(1)  :iX_E0(1)  )

    REAL(DP) :: &
      U_F_cell(nDOFX,0:3, &
            iX_B0(2)  :iX_E0(2)  , &
            iX_B0(3)  :iX_E0(3)  , &
            iX_B0(1)-1:iX_E0(1)+1)

    INTEGER  :: iNodeX1
    REAL(DP) :: X1

    IF( iX_E0(1) .EQ. iX_B0(1) )THEN

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO i = 0, 3
      DO iNodeX = 1, nDOFX

          dU_d_dX1_Out  (iNodeX,i,iX1,iX2,iX3) = Zero

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

      ! --- Permute Geometry Fields ---

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

      ! --- Compute Metric Components from Scale Factors ---

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

        END DO

      END DO
      END DO
      END DO

      ! --- Permute Four-velocity ---

      DO iX1 = iX_B0(1)-1, iX_E0(1)+1
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
  
        ! DO iPF = 1, nPF
        DO iNodeX = 1, nDOFX
  
          U_F_K(iNodeX,0,iX2,iX3,iX1) = One / &
                                        SQRT( One &
                                            - U_F(iNodeX,iX1,iX2,iX3,iPF_V1)**2 &
                                            - U_F(iNodeX,iX1,iX2,iX3,iPF_V2)**2 & 
                                            - U_F(iNodeX,iX1,iX2,iX3,iPF_V3)**2 )
          U_F_K(iNodeX,1,iX2,iX3,iX1) = U_F_K(iNodeX,0,iX2,iX3,iX1) * U_F(iNodeX,iX1,iX2,iX3,iPF_V1)
          U_F_K(iNodeX,2,iX2,iX3,iX1) = U_F_K(iNodeX,0,iX2,iX3,iX1) * U_F(iNodeX,iX1,iX2,iX3,iPF_V2)
          U_F_K(iNodeX,3,iX2,iX3,iX1) = U_F_K(iNodeX,0,iX2,iX3,iX1) * U_F(iNodeX,iX1,iX2,iX3,iPF_V3)
          U_F_K(iNodeX,0,iX2,iX3,iX1) = - U_F_K(iNodeX,0,iX2,iX3,iX1)
  
        END DO
        ! END DO
  
      END DO
      END DO
      END DO

    ! --- Interpolate Left State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nX1_X*4, nDOFX, One, LX_X1_Up, nDOFX_X1, &
             U_F_K(1,0,iX_B0(2),iX_B0(3),iX_B0(1)-1), nDOFX, Zero, &
             U_F_L(1,0,iX_B0(2),iX_B0(3),iX_B0(1)  ), nDOFX_X1 )

    ! --- Interpolate Right State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X1, nX1_X*4, nDOFX, One, LX_X1_Dn, nDOFX_X1, &
             U_F_K(1,0,iX_B0(2),iX_B0(3),iX_B0(1)  ), nDOFX, Zero, &
             U_F_R(1,0,iX_B0(2),iX_B0(3),iX_B0(1)  ), nDOFX_X1 )

    DO iX1 = iX_B0(1), iX_E0(1)+1
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)

      DO iNodeX = 1, nDOFX_X1

          ! --- Left States ---

          uV_L(0) = U_F_L(iNodeX,0,iX2,iX3,iX1)
          uV_L(1) = U_F_L(iNodeX,1,iX2,iX3,iX1)
          uV_L(2) = U_F_L(iNodeX,2,iX2,iX3,iX1)
          uV_L(3) = U_F_L(iNodeX,3,iX2,iX3,iX1)

          ! --- Right States ---

          uV_R(0) = U_F_R(iNodeX,0,iX2,iX3,iX1)
          uV_R(1) = U_F_R(iNodeX,1,iX2,iX3,iX1)
          uV_R(2) = U_F_R(iNodeX,2,iX2,iX3,iX1)
          uV_R(3) = U_F_R(iNodeX,3,iX2,iX3,iX1)

        CALL FaceFourVelocity_X1 &
               ( uV_L(0), uV_L(1), uV_L(2), uV_L(3), &
                 uV_R(0), uV_R(1), uV_R(2), uV_R(3), &
                 uV_F(0), uV_F(1), uV_F(2), uV_F(3) )

        U_d_X1(iNodeX,0,iX2,iX3,iX1) &
          = uV_F(0) * WeightsX_X1(iNodeX)
        U_d_X1(iNodeX,1,iX2,iX3,iX1) &
          = uV_F(1) * WeightsX_X1(iNodeX)
        U_d_X1(iNodeX,2,iX2,iX3,iX1) &
          = uV_F(2) * WeightsX_X1(iNodeX)
        U_d_X1(iNodeX,3,iX2,iX3,iX1) &
          = uV_F(3) * WeightsX_X1(iNodeX)

      END DO

    END DO
    END DO
    END DO

    ! --- Surface Contributions ---

    ! --- Contribution from Left Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 4*nK_X, nDOFX_X1, - One, LX_X1_Dn, nDOFX_X1, &
             U_d_X1(1,0,iX_B0(2),iX_B0(3),iX_B0(1)  ), nDOFX_X1, Zero, &
             dU_d_dX1(:,0:3,:,:,:), nDOFX )

    ! --- Contribution from Right Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 4*nK_X, nDOFX_X1, + One, LX_X1_Up, nDOFX_X1, &
             U_d_X1(1,0,iX_B0(2),iX_B0(3),iX_B0(1)+1), nDOFX_X1, One,  &
             dU_d_dX1(:,0:3,:,:,:), nDOFX )

    ! -------------------
    ! --- Volume Term ---
    ! -------------------

    DO iX1 = iX_B0(1), iX_E0(1)
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)

      DO iNodeX = 1, nDOFX

        U_d_K(iNodeX,0,iX2,iX3,iX1) &
          = U_F_K(iNodeX,0,iX2,iX3,iX1) * WeightsX_q(iNodeX)
        U_d_K(iNodeX,1,iX2,iX3,iX1) &
          = U_F_K(iNodeX,1,iX2,iX3,iX1) * WeightsX_q(iNodeX)
        U_d_K(iNodeX,2,iX2,iX3,iX1) &
          = U_F_K(iNodeX,2,iX2,iX3,iX1) * WeightsX_q(iNodeX)
        U_d_K(iNodeX,3,iX2,iX3,iX1) &
          = U_F_K(iNodeX,3,iX2,iX3,iX1) * WeightsX_q(iNodeX)

      END DO

    END DO
    END DO
    END DO

    ! --- Volume Contributions ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 4*nK_X, nDOFX, - One, dLXdX1_q, nDOFX, &
             U_d_K(:,0:3,:,:,:), nDOFX, One, dU_d_dX1(:,0:3,:,:,:), nDOFX )

    DO iX1 = iX_B0(1), iX_E0(1)
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)

      DO i      = 0, 3
      DO iNodeX = 1, nDOFX

        dU_d_dX1(iNodeX,i,iX2,iX3,iX1) &
         = dU_d_dX1(iNodeX,i,iX2,iX3,iX1) &
             / ( WeightsX_q(iNodeX) * dX1(iX1) )

      END DO
      END DO

    END DO
    END DO
    END DO

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        dU_d_dX1_Out(iNodeX,0,iX1,iX2,iX3) &
          = dU_d_dX1(iNodeX,0,iX2,iX3,iX1)
        dU_d_dX1_Out(iNodeX,1,iX1,iX2,iX3) &
          = dU_d_dX1(iNodeX,1,iX2,iX3,iX1)
        dU_d_dX1_Out(iNodeX,2,iX1,iX2,iX3) &
          = dU_d_dX1(iNodeX,2,iX2,iX3,iX1)
        dU_d_dX1_Out(iNodeX,3,iX1,iX2,iX3) &
          = dU_d_dX1(iNodeX,3,iX2,iX3,iX1)

      END DO

    END DO
    END DO
    END DO

    END ASSOCIATE

    ! ! --- Compute Derivative exactly for Transparent Shock ---

    ! DO iX1 = iX_B0(1), iX_E0(1)
    ! DO iX3 = iX_B0(3), iX_E0(3)
    ! DO iX2 = iX_B0(2), iX_E0(2)

    !   DO iNodeX = 1, nDOFX

    !     ! --- Compute four-velocity ---
    !     U_F_P(iNodeX,0,iX2,iX3,iX1) = One / SQRT( One &
    !                                               - U_F(iNodeX,iX1,iX2,iX3,iPF_V1)**2 &
    !                                               - U_F(iNodeX,iX1,iX2,iX3,iPF_V2)**2 &
    !                                               - U_F(iNodeX,iX1,iX2,iX3,iPF_V3)**2 )
    !     U_F_P(iNodeX,1,iX2,iX3,iX1) = U_F(iNodeX,iX1,iX2,iX3,iPF_V1)
    !     U_F_P(iNodeX,2,iX2,iX3,iX1) = U_F(iNodeX,iX1,iX2,iX3,iPF_V2)
    !     U_F_P(iNodeX,3,iX2,iX3,iX1) = U_F(iNodeX,iX1,iX2,iX3,iPF_V3)

    !     ! --- Compute derivative of three-velocity ---
    !     iNodeX1 = NodeNumberTableX(1,iNodeX)
    !     X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

    !     dV_d_dX1(iNodeX,iX2,iX3,iX1) = -0.1_DP / ( COSH( (X1 - One) / 1.0d-3 )**2 * Two * 1.0d-3)
    !     dW_dX1  (iNodeX,iX2,iX3,iX1) = U_F_P(iNodeX,1,iX2,iX3,iX1) * &
    !                                    dV_d_dX1(iNodeX,iX2,iX3,iX1) * &
    !                                    U_F_P(iNodeX,0,iX2,iX3,iX1)**3

    !   END DO

    ! END DO
    ! END DO
    ! END DO

    ! DO iX1 = iX_B0(1), iX_E0(1)
    ! DO iX3 = iX_B0(3), iX_E0(3)
    ! DO iX2 = iX_B0(2), iX_E0(2)

    !   DO iNodeX = 1, nDOFX

    !     dU_d_dX1(iNodeX,0,iX2,iX3,iX1) &
    !       = - dW_dX1(iNodeX,iX2,iX3,iX1)
    !     dU_d_dX1(iNodeX,1,iX2,iX3,iX1) &
    !       = U_F_P(iNodeX,1,iX2,iX3,iX1) * dW_dX1(iNodeX,iX2,iX3,iX1) + &
    !         U_F_P(iNodeX,0,iX2,iX3,iX1) * dV_d_dX1(iNodeX,iX2,iX3,iX1)

    !   END DO

    ! END DO
    ! END DO
    ! END DO

    ! DO iX3 = iX_B0(3), iX_E0(3)
    ! DO iX2 = iX_B0(2), iX_E0(2)
    ! DO iX1 = iX_B0(1), iX_E0(1)

    !   DO iNodeX = 1, nDOFX

    !     dU_d_dX1_Out(iNodeX,0,iX1,iX2,iX3) &
    !       = dU_d_dX1(iNodeX,0,iX2,iX3,iX1)
    !     dU_d_dX1_Out(iNodeX,1,iX1,iX2,iX3) &
    !       = dU_d_dX1(iNodeX,1,iX2,iX3,iX1)
    !     dU_d_dX1_Out(iNodeX,2,iX1,iX2,iX3) = Zero
    !     dU_d_dX1_Out(iNodeX,3,iX1,iX2,iX3) = Zero

    !   END DO

    ! END DO
    ! END DO
    ! END DO

    ! ! --- Compute Derivative using cell average ---

    ! DO iX1 = iX_B0(1)-1, iX_E0(1)+1
    !   DO iX3 = iX_B0(3), iX_E0(3)
    !   DO iX2 = iX_B0(2), iX_E0(2)
  
    !     DO iNodeX = 1, nDOFX
  
    !       U_F_cell(iNodeX,0,iX2,iX3,iX1) = Zero
    !       U_F_cell(iNodeX,1,iX2,iX3,iX1) = Zero
    !       U_F_cell(iNodeX,2,iX2,iX3,iX1) = Zero
    !       U_F_cell(iNodeX,3,iX2,iX3,iX1) = Zero
  
    !     END DO
  
    !   END DO
    !   END DO
    !   END DO

    ! DO iX1 = iX_B0(1)-1, iX_E0(1)+1
    !   DO iX3 = iX_B0(3), iX_E0(3)
    !   DO iX2 = iX_B0(2), iX_E0(2)
  
    !     DO iNodeX = 1, nDOFX
  
    !       ! --- Compute four-velocity ---
    !       U_F_cell(iNodeX,0,iX2,iX3,iX1) = One / SQRT( One &
    !                                                 - U_F(iNodeX,iX1,iX2,iX3,iPF_V1)**2 &
    !                                                 - U_F(iNodeX,iX1,iX2,iX3,iPF_V2)**2 &
    !                                                 - U_F(iNodeX,iX1,iX2,iX3,iPF_V3)**2 )
    !       U_F_cell(iNodeX,1,iX2,iX3,iX1) = U_F_cell(iNodeX,0,iX2,iX3,iX1) * U_F(iNodeX,iX1,iX2,iX3,iPF_V1)
    !       U_F_cell(iNodeX,2,iX2,iX3,iX1) = U_F_cell(iNodeX,0,iX2,iX3,iX1) * U_F(iNodeX,iX1,iX2,iX3,iPF_V2)
    !       U_F_cell(iNodeX,3,iX2,iX3,iX1) = U_F_cell(iNodeX,0,iX2,iX3,iX1) * U_F(iNodeX,iX1,iX2,iX3,iPF_V3)
  
    !       ! ! --- Compute cellaverage of four-velocity ---
    !       ! iNodeX1 = NodeNumberTableX(1,iNodeX)
    !       ! X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
  
    !       ! U_F_cell(iNodeX,0,iX2,iX3,iX1) = - U_F_cell(iNodeX,0,iX2,iX3,iX1) * WeightsX1(iNodeX1)
    !       ! U_F_cell(iNodeX,1,iX2,iX3,iX1) =   U_F_cell(iNodeX,1,iX2,iX3,iX1) * WeightsX1(iNodeX1)
    !       ! U_F_cell(iNodeX,2,iX2,iX3,iX1) =   U_F_cell(iNodeX,2,iX2,iX3,iX1) * WeightsX1(iNodeX1)
    !       ! U_F_cell(iNodeX,3,iX2,iX3,iX1) =   U_F_cell(iNodeX,3,iX2,iX3,iX1) * WeightsX1(iNodeX1)
  
    !     END DO
  
    !   END DO
    !   END DO
    !   END DO
  
    !   DO iX1 = iX_B0(1), iX_E0(1)
    !   DO iX3 = iX_B0(3), iX_E0(3)
    !   DO iX2 = iX_B0(2), iX_E0(2)
  
    !     DO iNodeX = 1, nDOFX
  
    !       dU_d_dX1(iNodeX,0,iX2,iX3,iX1) &
    !         = (-U_F_cell(iNodeX,0,iX2,iX3,iX1+1)+U_F_cell(iNodeX,0,iX2,iX3,iX1-1))/dX1(iX1)
    !       dU_d_dX1(iNodeX,1,iX2,iX3,iX1) &
    !         = (U_F_cell(iNodeX,1,iX2,iX3,iX1+1)-U_F_cell(iNodeX,1,iX2,iX3,iX1-1))/dX1(iX1)
  
    !     END DO
  
    !   END DO
    !   END DO
    !   END DO
  
    !   DO iX3 = iX_B0(3), iX_E0(3)
    !   DO iX2 = iX_B0(2), iX_E0(2)
    !   DO iX1 = iX_B0(1), iX_E0(1)
  
    !     DO iNodeX = 1, nDOFX
  
    !       dU_d_dX1_Out(iNodeX,0,iX1,iX2,iX3) &
    !         = dU_d_dX1(iNodeX,0,iX2,iX3,iX1)
    !       dU_d_dX1_Out(iNodeX,1,iX1,iX2,iX3) &
    !         = dU_d_dX1(iNodeX,1,iX2,iX3,iX1)
    !       dU_d_dX1_Out(iNodeX,2,iX1,iX2,iX3) = Zero
    !       dU_d_dX1_Out(iNodeX,3,iX1,iX2,iX3) = Zero
  
    !     END DO
  
    !   END DO
    !   END DO
    !   END DO
    
    ! END ASSOCIATE

  END SUBROUTINE ComputeWeakDerivatives_X1

  SUBROUTINE ComputeWeakDerivatives_X2 &
    ( iX_B0, iX_E0, iX_B1, iX_E1, GX, U_F, dU_d_dX2_Out )

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
          1:nPF)
    REAL(DP), INTENT(out) :: &
      dU_d_dX2_Out &
        (1:nDOFX,0:3, &
         iX_B0(1):iX_E0(1), &
         iX_B0(2):iX_E0(2), &
         iX_B0(3):iX_E0(3))

    INTEGER  :: iNodeX
    INTEGER  :: iX1, iX2, iX3, i
    INTEGER  :: iPF, iPF_V
    INTEGER  :: iGF, iGF_h, iGF_Gm_dd
    INTEGER  :: nX(3), nX_X2(3), nK_X, nX2_X
    REAL(DP) :: uV_L(0:3), uV_R(0:3), uV_F(0:3), uV_K
    REAL(DP) :: W_K, W_L, W_R, vMagSq_K, vMagSq_L, vMagSq_R

    ! --- Geometry Fields ---

    REAL(DP) :: &
      GX_K   (nDOFX,nGF, &
              iX_B0(1)  :iX_E0(1)  , &
              iX_B0(3)  :iX_E0(3)  , &
              iX_B0(2)-1:iX_E0(2)+1)
    REAL(DP) :: &
      GX_F   (nDOFX_X2,nGF, &
              iX_B0(1)  :iX_E0(1)  , &
              iX_B0(3)  :iX_E0(3)  , &
              iX_B0(2)  :iX_E0(2)+1)

    ! --- Primitive Fluid Fields ---

    REAL(DP) :: &
      U_F_K(nDOFX,0:3, &
            iX_B0(1)  :iX_E0(1)  , &
            iX_B0(3)  :iX_E0(3)  , &
            iX_B0(2)-1:iX_E0(2)+1)
    REAL(DP) :: &
      U_F_L(nDOFX_X2,0:3, &
            iX_B0(1)  :iX_E0(1)  , &
            iX_B0(3)  :iX_E0(3)  , &
            iX_B0(2)  :iX_E0(2)+1)
    REAL(DP) :: &
      U_F_R(nDOFX_X2,0:3, &
            iX_B0(1)  :iX_E0(1)  , &
            iX_B0(3)  :iX_E0(3)  , &
            iX_B0(2)  :iX_E0(2)+1)

    ! --- Velocities ---

    REAL(DP) :: &
      U_d_X2  (nDOFX_X2,0:3, &
               iX_B0(1)  :iX_E0(1)  , &
               iX_B0(3)  :iX_E0(3)  , &
               iX_B0(2)  :iX_E0(2)+1)
    REAL(DP) :: &
      U_d_K   (nDOFX,0:3, &
               iX_B0(1)  :iX_E0(1)  , &
               iX_B0(3)  :iX_E0(3)  , &
               iX_B0(2)  :iX_E0(2)  )
    REAL(DP) :: &
      dU_d_dX2(nDOFX,0:3, &
               iX_B0(1)  :iX_E0(1)  , &
               iX_B0(3)  :iX_E0(3)  , &
               iX_B0(2)  :iX_E0(2)  )

    IF( iX_E0(2) .EQ. iX_B0(2) )THEN

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO i = 0, 3
      DO iNodeX = 1, nDOFX

          dU_d_dX2_Out  (iNodeX,i,iX1,iX2,iX3) = Zero

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

      ! --- Permute Geometry Fields ---

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

      ! --- Compute Metric Components from Scale Factors ---

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

        END DO

      END DO
      END DO
      END DO

      ! --- Permute Four-velocity ---

      DO iX2 = iX_B0(2)-1, iX_E0(2)+1
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX1 = iX_B0(1), iX_E0(1)
  
        ! DO iPF = 1, nPF
        DO iNodeX = 1, nDOFX
  
          U_F_K(iNodeX,0,iX1,iX3,iX2) = One / &
                                        SQRT( One &
                                            - U_F(iNodeX,iX1,iX2,iX3,iPF_V1)**2 &
                                            - U_F(iNodeX,iX1,iX2,iX3,iPF_V2)**2 &
                                            - U_F(iNodeX,iX1,iX2,iX3,iPF_V1)**2 )
          U_F_K(iNodeX,1,iX1,iX3,iX2) = U_F_K(iNodeX,0,iX1,iX3,iX2) * U_F(iNodeX,iX1,iX2,iX3,iPF_V1)
          U_F_K(iNodeX,2,iX1,iX3,iX2) = U_F_K(iNodeX,0,iX1,iX3,iX2) * U_F(iNodeX,iX1,iX2,iX3,iPF_V2)
          U_F_K(iNodeX,3,iX1,iX3,iX2) = U_F_K(iNodeX,0,iX1,iX3,iX2) * U_F(iNodeX,iX1,iX2,iX3,iPF_V3)
          U_F_K(iNodeX,0,iX1,iX3,iX2) = - U_F_K(iNodeX,0,iX1,iX3,iX2)
  
        END DO
        ! END DO
  
      END DO
      END DO
      END DO

    ! --- Interpolate Left State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nX2_X*4, nDOFX, One, LX_X2_Up, nDOFX_X2, &
             U_F_K(1,0,iX_B0(1),iX_B0(3),iX_B0(2)-1), nDOFX, Zero, &
             U_F_L(1,0,iX_B0(1),iX_B0(3),iX_B0(2)  ), nDOFX_X2 )

    ! --- Interpolate Right State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X2, nX2_X*4, nDOFX, One, LX_X2_Dn, nDOFX_X2, &
             U_F_K(1,0,iX_B0(1),iX_B0(3),iX_B0(2)  ), nDOFX, Zero, &
             U_F_R(1,0,iX_B0(1),iX_B0(3),iX_B0(2)  ), nDOFX_X2 )

    DO iX2 = iX_B0(2), iX_E0(2)+1
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX_X2

          ! --- Left States ---

          uV_L(0) = U_F_L(iNodeX,0,iX1,iX3,iX2)
          uV_L(1) = U_F_L(iNodeX,1,iX1,iX3,iX2)
          uV_L(2) = U_F_L(iNodeX,2,iX1,iX3,iX2)
          uV_L(3) = U_F_L(iNodeX,3,iX1,iX3,iX2)

          ! --- Right States ---

          uV_R(0) = U_F_R(iNodeX,0,iX1,iX3,iX2)
          uV_R(1) = U_F_R(iNodeX,1,iX1,iX3,iX2)
          uV_R(2) = U_F_R(iNodeX,2,iX1,iX3,iX2)
          uV_R(3) = U_F_R(iNodeX,3,iX1,iX3,iX2)

        CALL FaceFourVelocity_X1 &
               ( uV_L(0), uV_L(1), uV_L(2), uV_L(3), &
                 uV_R(0), uV_R(1), uV_R(2), uV_R(3), &
                 uV_F(0), uV_F(1), uV_F(2), uV_F(3) )

        U_d_X2(iNodeX,0,iX1,iX3,iX2) &
          = uV_F(0) * WeightsX_X2(iNodeX)
        U_d_X2(iNodeX,1,iX1,iX3,iX2) &
          = uV_F(1) * WeightsX_X2(iNodeX)
        U_d_X2(iNodeX,2,iX1,iX3,iX2) &
          = uV_F(2) * WeightsX_X2(iNodeX)
        U_d_X2(iNodeX,3,iX1,iX3,iX2) &
          = uV_F(3) * WeightsX_X2(iNodeX)

      END DO

    END DO
    END DO
    END DO

    ! --- Surface Contributions ---

    ! --- Contribution from Left Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 4*nK_X, nDOFX_X2, - One, LX_X2_Dn, nDOFX_X2, &
             U_d_X2(1,0,iX_B0(1),iX_B0(3),iX_B0(2)  ), nDOFX_X2, Zero, &
             dU_d_dX2(:,0:3,:,:,:), nDOFX )

    ! --- Contribution from Right Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 4*nK_X, nDOFX_X2, + One, LX_X2_Up, nDOFX_X2, &
             U_d_X2(1,0,iX_B0(1),iX_B0(3),iX_B0(2)+1), nDOFX_X2, One,  &
             dU_d_dX2(:,0:3,:,:,:), nDOFX )

    ! -------------------
    ! --- Volume Term ---
    ! -------------------

    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        U_d_K(iNodeX,0,iX1,iX3,iX2) &
          = U_F_K(iNodeX,0,iX1,iX3,iX2) * WeightsX_q(iNodeX)
        U_d_K(iNodeX,1,iX1,iX3,iX2) &
          = U_F_K(iNodeX,1,iX1,iX3,iX2) * WeightsX_q(iNodeX)
        U_d_K(iNodeX,2,iX1,iX3,iX2) &
          = U_F_K(iNodeX,2,iX1,iX3,iX2) * WeightsX_q(iNodeX)
        U_d_K(iNodeX,3,iX1,iX3,iX2) &
          = U_F_K(iNodeX,3,iX1,iX3,iX2) * WeightsX_q(iNodeX)

      END DO

    END DO
    END DO
    END DO

    ! --- Volume Contributions ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 4*nK_X, nDOFX, - One, dLXdX2_q, nDOFX, &
             U_d_K(:,0:3,:,:,:), nDOFX, One, dU_d_dX2(:,0:3,:,:,:), nDOFX )

    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO i      = 0, 3
      DO iNodeX = 1, nDOFX

        dU_d_dX2(iNodeX,i,iX1,iX3,iX2) &
         = dU_d_dX2(iNodeX,i,iX1,iX3,iX2) &
             / ( WeightsX_q(iNodeX) * dX2(iX2) )

      END DO
      END DO

    END DO
    END DO
    END DO

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        dU_d_dX2_Out(iNodeX,0,iX1,iX2,iX3) &
          = dU_d_dX2(iNodeX,0,iX1,iX3,iX2)
        dU_d_dX2_Out(iNodeX,1,iX1,iX2,iX3) &
          = dU_d_dX2(iNodeX,1,iX1,iX3,iX2)
        dU_d_dX2_Out(iNodeX,2,iX1,iX2,iX3) &
          = dU_d_dX2(iNodeX,2,iX1,iX3,iX2)
        dU_d_dX2_Out(iNodeX,3,iX1,iX2,iX3) &
          = dU_d_dX2(iNodeX,3,iX1,iX3,iX2)

      END DO

    END DO
    END DO
    END DO

    END ASSOCIATE

  END SUBROUTINE ComputeWeakDerivatives_X2

  SUBROUTINE ComputeWeakDerivatives_X3 &
    ( iX_B0, iX_E0, iX_B1, iX_E1, GX, U_F, dU_d_dX3_Out )

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
          1:nPF)
    REAL(DP), INTENT(out) :: &
      dU_d_dX3_Out &
        (1:nDOFX,0:3, &
         iX_B0(1):iX_E0(1), &
         iX_B0(2):iX_E0(2), &
         iX_B0(3):iX_E0(3))

    INTEGER  :: iNodeX
    INTEGER  :: iX1, iX2, iX3, i
    INTEGER  :: iPF, iPF_V
    INTEGER  :: iGF, iGF_h, iGF_Gm_dd
    INTEGER  :: nX(3), nX_X3(3), nK_X, nX3_X
    REAL(DP) :: uV_L(0:3), uV_R(0:3), uV_F(0:3), uV_K
    REAL(DP) :: W_K, W_L, W_R, vMagSq_K, vMagSq_L, vMagSq_R

    ! --- Geometry Fields ---

    REAL(DP) :: &
      GX_K   (nDOFX,nGF, &
              iX_B0(1)  :iX_E0(1)  , &
              iX_B0(2)  :iX_E0(2)  , &
              iX_B0(3)-1:iX_E0(3)+1)
    REAL(DP) :: &
      GX_F   (nDOFX_X3,nGF, &
              iX_B0(1)  :iX_E0(1)  , &
              iX_B0(2)  :iX_E0(2)  , &
              iX_B0(3)  :iX_E0(3)+1)

    ! --- Primitive Fluid Fields ---

    REAL(DP) :: &
      U_F_K(nDOFX,0:3, &
            iX_B0(1)  :iX_E0(1)  , &
            iX_B0(2)  :iX_E0(2)  , &
            iX_B0(3)-1:iX_E0(3)+1)
    REAL(DP) :: &
      U_F_L(nDOFX_X3,0:3, &
            iX_B0(1)  :iX_E0(1)  , &
            iX_B0(2)  :iX_E0(2)  , &
            iX_B0(3)  :iX_E0(3)+1)
    REAL(DP) :: &
      U_F_R(nDOFX_X3,0:3, &
            iX_B0(1)  :iX_E0(1)  , &
            iX_B0(2)  :iX_E0(2)  , &
            iX_B0(3)  :iX_E0(3)+1)

    ! --- Velocities ---

    REAL(DP) :: &
      U_d_X3  (nDOFX_X3,0:3, &
               iX_B0(1)  :iX_E0(1)  , &
               iX_B0(2)  :iX_E0(2)  , &
               iX_B0(3)  :iX_E0(3)+1)
    REAL(DP) :: &
      U_d_K   (nDOFX,0:3, &
               iX_B0(1)  :iX_E0(1)  , &
               iX_B0(2)  :iX_E0(2)  , &
               iX_B0(3)  :iX_E0(3)  )
    REAL(DP) :: &
      dU_d_dX3(nDOFX,0:3, &
               iX_B0(1)  :iX_E0(1)  , &
               iX_B0(2)  :iX_E0(2)  , &
               iX_B0(3)  :iX_E0(3)  )

    IF( iX_E0(3) .EQ. iX_B0(3) )THEN

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO i = 0, 3
      DO iNodeX = 1, nDOFX

          dU_d_dX3_Out  (iNodeX,i,iX1,iX2,iX3) = Zero

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

      ! --- Permute Geometry Fields ---

      DO iX3 = iX_B0(3)-1, iX_E0(3)+1
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
  
        DO iGF = 1, nGF
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

      ! --- Compute Metric Components from Scale Factors ---

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

        END DO

      END DO
      END DO
      END DO

      ! --- Permute Four-velocity ---

      DO iX3 = iX_B0(3)-1, iX_E0(3)+1
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
  
        ! DO iPF = 1, nPF
        DO iNodeX = 1, nDOFX
  
          U_F_K(iNodeX,iPF,iX1,iX2,iX3) = U_F(iNodeX,iX1,iX2,iX3,iPF)
          U_F_K(iNodeX,0,iX1,iX2,iX3) = One / &
                                        SQRT( One &
                                            - U_F(iNodeX,iX1,iX2,iX3,iPF_V1)**2 &
                                            - U_F(iNodeX,iX1,iX2,iX3,iPF_V2)**2 &
                                            - U_F(iNodeX,iX1,iX2,iX3,iPF_V3)**2 )
          U_F_K(iNodeX,1,iX1,iX2,iX3) = U_F_K(iNodeX,0,iX1,iX2,iX3) * U_F(iNodeX,iX1,iX2,iX3,iPF_V1)
          U_F_K(iNodeX,2,iX1,iX2,iX3) = U_F_K(iNodeX,0,iX1,iX2,iX3) * U_F(iNodeX,iX1,iX2,iX3,iPF_V2)
          U_F_K(iNodeX,3,iX1,iX2,iX3) = U_F_K(iNodeX,0,iX1,iX2,iX3) * U_F(iNodeX,iX1,iX2,iX3,iPF_V3)
          U_F_K(iNodeX,0,iX1,iX2,iX3) = - U_F_K(iNodeX,0,iX1,iX2,iX3)
  
        END DO
        ! END DO
  
      END DO
      END DO
      END DO

    ! --- Interpolate Left State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X3, nX3_X*4, nDOFX, One, LX_X3_Up, nDOFX_X3, &
             U_F_K(1,0,iX_B0(1),iX_B0(2),iX_B0(3)-1), nDOFX, Zero, &
             U_F_L(1,0,iX_B0(1),iX_B0(2),iX_B0(3)  ), nDOFX_X3 )

    ! --- Interpolate Right State ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX_X3, nX3_X*4, nDOFX, One, LX_X3_Dn, nDOFX_X3, &
             U_F_K(1,0,iX_B0(1),iX_B0(2),iX_B0(3)  ), nDOFX, Zero, &
             U_F_R(1,0,iX_B0(1),iX_B0(2),iX_B0(3)  ), nDOFX_X3 )

    DO iX3 = iX_B0(3), iX_E0(3)+1
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX_X3

          ! --- Left States ---

          uV_L(0) = U_F_L(iNodeX,0,iX1,iX2,iX3)
          uV_L(1) = U_F_L(iNodeX,1,iX1,iX2,iX3)
          uV_L(2) = U_F_L(iNodeX,2,iX1,iX2,iX3)
          uV_L(3) = U_F_L(iNodeX,3,iX1,iX2,iX3)

          ! --- Right States ---

          uV_R(0) = U_F_R(iNodeX,0,iX1,iX2,iX3)
          uV_R(1) = U_F_R(iNodeX,1,iX1,iX2,iX3)
          uV_R(2) = U_F_R(iNodeX,2,iX1,iX2,iX3)
          uV_R(3) = U_F_R(iNodeX,3,iX1,iX2,iX3)

        CALL FaceFourVelocity_X1 &
               ( uV_L(0), uV_L(1), uV_L(2), uV_L(3), &
                 uV_R(0), uV_R(1), uV_R(2), uV_R(3), &
                 uV_F(0), uV_F(1), uV_F(2), uV_F(3) )

        U_d_X3(iNodeX,0,iX1,iX2,iX3) &
          = uV_F(0) * WeightsX_X3(iNodeX)
        U_d_X3(iNodeX,1,iX1,iX2,iX3) &
          = uV_F(1) * WeightsX_X3(iNodeX)
        U_d_X3(iNodeX,2,iX1,iX2,iX3) &
          = uV_F(2) * WeightsX_X3(iNodeX)
        U_d_X3(iNodeX,3,iX1,iX2,iX3) &
          = uV_F(3) * WeightsX_X3(iNodeX)

      END DO

    END DO
    END DO
    END DO

    ! --- Surface Contributions ---

    ! --- Contribution from Left Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 4*nK_X, nDOFX_X3, - One, LX_X3_Dn, nDOFX_X3, &
             U_d_X3(1,0,iX_B0(1),iX_B0(2),iX_B0(3)  ), nDOFX_X3, Zero, &
             dU_d_dX3(:,0:3,:,:,:), nDOFX )

    ! --- Contribution from Right Face ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 4*nK_X, nDOFX_X3, + One, LX_X3_Up, nDOFX_X3, &
             U_d_X3(1,0,iX_B0(1),iX_B0(2),iX_B0(3)+1), nDOFX_X3, One,  &
             dU_d_dX3(:,0:3,:,:,:), nDOFX )

    ! -------------------
    ! --- Volume Term ---
    ! -------------------

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        U_d_K(iNodeX,0,iX1,iX2,iX3) &
          = U_F_K(iNodeX,0,iX1,iX2,iX3) * WeightsX_q(iNodeX)
        U_d_K(iNodeX,1,iX1,iX2,iX3) &
          = U_F_K(iNodeX,1,iX1,iX2,iX3) * WeightsX_q(iNodeX)
        U_d_K(iNodeX,2,iX1,iX2,iX3) &
          = U_F_K(iNodeX,2,iX1,iX2,iX3) * WeightsX_q(iNodeX)
        U_d_K(iNodeX,3,iX1,iX2,iX3) &
          = U_F_K(iNodeX,3,iX1,iX2,iX3) * WeightsX_q(iNodeX)

      END DO

    END DO
    END DO
    END DO

    ! --- Volume Contributions ---

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 4*nK_X, nDOFX, - One, dLXdX3_q, nDOFX, &
             U_d_K(:,0:3,:,:,:), nDOFX, One, dU_d_dX3(:,0:3,:,:,:), nDOFX )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO i      = 0, 3
      DO iNodeX = 1, nDOFX

        dU_d_dX3(iNodeX,i,iX1,iX2,iX3) &
         = dU_d_dX3(iNodeX,i,iX1,iX2,iX3) &
             / ( WeightsX_q(iNodeX) * dX3(iX3) )

      END DO
      END DO

    END DO
    END DO
    END DO

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        dU_d_dX3_Out(iNodeX,0,iX1,iX2,iX3) &
          = dU_d_dX3(iNodeX,0,iX1,iX2,iX3)
        dU_d_dX3_Out(iNodeX,1,iX1,iX2,iX3) &
          = dU_d_dX3(iNodeX,1,iX1,iX2,iX3)
        dU_d_dX3_Out(iNodeX,2,iX1,iX2,iX3) &
          = dU_d_dX3(iNodeX,2,iX1,iX2,iX3)
        dU_d_dX3_Out(iNodeX,3,iX1,iX2,iX3) &
          = dU_d_dX3(iNodeX,3,iX1,iX2,iX3)

      END DO

    END DO
    END DO
    END DO

    END ASSOCIATE

  END SUBROUTINE ComputeWeakDerivatives_X3

  SUBROUTINE FaceVelocity_X1 &
    ( V1_L, V2_L, V3_L, V1_R, V2_R, V3_R, V1_F, V2_F, V3_F )

    REAL(DP), INTENT(in)  :: V1_L, V2_L, V3_L
    REAL(DP), INTENT(in)  :: V1_R, V2_R, V3_R
    REAL(DP), INTENT(out) :: V1_F, V2_F, V3_F

    ! --- Average Left and Right States ---

    V1_F = Half * ( V1_L + V1_R )
    V2_F = Half * ( V2_L + V2_R )
    V3_F = Half * ( V3_L + V3_R )

    RETURN

  END SUBROUTINE FaceVelocity_X1

  SUBROUTINE FaceFourVelocity_X1 &
    ( U0_L, U1_L, U2_L, U3_L, U0_R, U1_R, U2_R, U3_R, U0_F, U1_F, U2_F, U3_F )

    REAL(DP), INTENT(in)  :: U0_L, U1_L, U2_L, U3_L
    REAL(DP), INTENT(in)  :: U0_R, U1_R, U2_R, U3_R
    REAL(DP), INTENT(out) :: U0_F, U1_F, U2_F, U3_F

    ! --- Average Left and Right States ---

    U0_F = Half * ( U0_L + U0_R )
    U1_F = Half * ( U1_L + U1_R )
    U2_F = Half * ( U2_L + U2_R )
    U3_F = Half * ( U3_L + U3_R )

    RETURN

  END SUBROUTINE FaceFourVelocity_X1


END MODULE TwoMoment_UtilitiesModule_FMC
