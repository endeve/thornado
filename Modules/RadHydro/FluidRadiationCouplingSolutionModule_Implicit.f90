MODULE FluidRadiationCouplingSolutionModule_Implicit

  USE KindModule, ONLY: &
    DP, Pi, FourPi
  USE UnitsModule, ONLY: &
    SpeedOfLight, &
    BoltzmannConstant, &
    PlanckConstant, &
    Gram, &
    Centimeter, &
    Kelvin
  USE ProgramHeaderModule, ONLY: &
    nX, nNodesX, nDOFX, &
    nE, nNodesE, nDOFE
  USE UtilitiesModule, ONLY: &
    WriteVector, &
    WriteMatrix
  USE FluidFieldsModule, ONLY: &
    iPF_D, iPF_E, iPF_Ne, nPF, &
    iAF_T, iAF_E, iAF_Ye, iAF_Me, iAF_Mp, iAF_Mn, nAF
  USE RadiationFieldsModule, ONLY: &
    iPR_D, iPR_I1, iPR_I2, iPR_I3, nPR
  USE MomentEquationsUtilitiesModule, ONLY: &
    ComputeConserved, &
    ComputePrimitive
  USE FluidRadiationCouplingUtilitiesModule, ONLY: &
    InitializeNodes, &
    InitializeWeights, &
    InitializeFluidFields, &
    FinalizeFluidFields, &
    InitializeRadiationFields, &
    FinalizeRadiationFields
  USE EquationOfStateModule, ONLY: &
    BaryonMass, &
    ComputeThermodynamicStates_Auxiliary, &
    ComputeSpecificInternalEnergy, &
    ComputeElectronChemicalPotential, &
    ComputeProtonChemicalPotential, &
    ComputeNeutronChemicalPotential, &
    ApplyEquationOfState
  USE OpacityModule, ONLY: &
    ComputeAbsorptionCoefficients

  IMPLICIT NONE
  PRIVATE

  INTEGER :: nNodesX_G, nNodesE_G, iFRC_Ne, iFRC_E
  INTEGER, PARAMETER :: iOld = 0, iNew = 1
  REAL(DP), DIMENSION(:),     ALLOCATABLE :: E_N, W2_N, W3_N
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE :: uPF_N, uAF_N
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE :: duAFdT_N, duAFdYe_N
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE :: Chi, dChidT, dChidYe
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE :: f_FD, df_FDdT, df_FDdYe
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: uPR_N

  REAL(DP), DIMENSION(:),     ALLOCATABLE :: Floor_FRC
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE :: FVEC_FRC, dU_FRC
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: U_FRC, FJAC_FRC

  PUBLIC :: CoupleFluidRadiation_Implicit_EmissionAbsorption

CONTAINS


  SUBROUTINE CoupleFluidRadiation_Implicit_EmissionAbsorption &
               ( dt, iX_Begin, iX_End )

    REAL(DP),              INTENT(in) :: dt
    INTEGER, DIMENSION(3), INTENT(in) :: iX_Begin, iX_End

    CALL ComputePrimitive( iX_Begin = iX_Begin, iX_End = iX_End )

    CALL InitializeFluidRadiationCoupling

    CALL CoupleFluidRadiation_EmissionAbsorption( dt )

    CALL FinalizeFluidRadiationCoupling

    CALL ApplyEquationOfState

    CALL ComputeConserved( iX_Begin = iX_Begin, iX_End = iX_End )

  END SUBROUTINE CoupleFluidRadiation_Implicit_EmissionAbsorption


  SUBROUTINE InitializeFluidRadiationCoupling

    nNodesX_G = PRODUCT(nX) * nDOFX
    nNodesE_G =         nE  * nDOFE

    iFRC_Ne = nNodesE_G + 1
    iFRC_E  = nNodesE_G + 2

    ALLOCATE( E_N(nNodesE_G) )
    CALL InitializeNodes( E_N )

    ALLOCATE( W2_N(nNodesE_G), W3_N(nNodesE_G) )
    CALL InitializeWeights( W2_N, W3_N )

    ALLOCATE( uPF_N(nPF, nNodesX_G) )
    ALLOCATE( uAF_N(nAF, nNodesX_G) )
    CALL InitializeFluidFields( uPF_N, uAF_N )

    ALLOCATE( duAFdT_N (nAF, nNodesX_G) )
    ALLOCATE( duAFdYe_N(nAF, nNodesX_G) )

    ALLOCATE( uPR_N(nNodesE_G, nPR, nNodesX_G) )
    CALL InitializeRadiationFields( uPR_N )

    ALLOCATE &
      ( Chi    (nNodesE_G, nNodesX_G), &
        dChidT (nNodesE_G, nNodesX_G), &
        dChidYe(nNodesE_G, nNodesX_G) )

    ALLOCATE &
      ( f_FD    (nNodesE_G, nNodesX_G), &
        df_FDdT (nNodesE_G, nNodesX_G), &
        df_FDdYe(nNodesE_G, nNodesX_G) )

    ALLOCATE( Floor_FRC(nNodesE_G+2) )
    ALLOCATE( U_FRC    (nNodesE_G+2, nNodesX_G, 0:1) )
    ALLOCATE( FVEC_FRC (nNodesE_G+2, nNodesX_G) )
    ALLOCATE( dU_FRC   (nNodesE_G+2, nNodesX_G) )
    ALLOCATE( FJAC_FRC (nNodesE_G+2, nNodesE_G+2, nNodesX_G) )

    Floor_FRC = 0.0_DP
    Floor_FRC(1:nNodesE_G) = 1.0d-16

  END SUBROUTINE InitializeFluidRadiationCoupling


  SUBROUTINE FinalizeFluidRadiationCoupling

    CALL FinalizeFluidFields( uPF_N, uAF_N )

    CALL FinalizeRadiationFields( uPR_N )

    DEALLOCATE( E_N, W2_N, W3_N )
    DEALLOCATE( uPF_N, uPR_N )
    DEALLOCATE( uAF_N, duAFdT_N, duAFdYe_N )
    DEALLOCATE( Chi, dChidT, dChidYe )
    DEALLOCATE( f_FD, df_FDdT, df_FDdYe )

    DEALLOCATE( Floor_FRC, U_FRC, FVEC_FRC, dU_FRC, FJAC_FRC )

  END SUBROUTINE FinalizeFluidRadiationCoupling


  SUBROUTINE CoupleFluidRadiation_EmissionAbsorption( dt )

    REAL(DP), INTENT(in) :: dt

    LOGICAL, DIMENSION(1:nNodesX_G) :: Converged
    INTEGER :: iter, iX, iX_MAX
    REAL(DP) :: Norm, MaxNorm
    REAL(DP), PARAMETER :: NewtonTol = 1.0d-8

    CALL SetStates_FRC &
           ( uPR_N(:,iPR_D,:), uPF_N(iPF_Ne,:), uPF_N(iPF_E,:), iOld )

    Converged  = .FALSE.
    iter       = 0

    DO WHILE( .NOT. ALL( Converged ) )

      iter = iter + 1

      CALL ComputeThermodynamicStates_Auxiliary &
             ( uPF_N(iPF_D,:), uPF_N(iPF_E,:), uPF_N(iPF_Ne,:), &
               uAF_N(iAF_T,:), uAF_N(iAF_E,:), uAF_N(iAF_Ye,:) )

      CALL SetStates_FRC &
             ( uPR_N(:,iPR_D,:), uPF_N(iPF_Ne,:), uPF_N(iPF_E,:), iNew )

      CALL SetRates_EmissionAbsorption ! --- Computes Absoption Coefficients
                                       !       and Derivatives

      CALL SetEquilibrium_EmissionAbsorption ! --- Computes FD Distribution
                                             !       and Derivatives

      CALL ComputeSpecificInternalEnergy & ! --- Computes Derivatives
             ( uPF_N(iPF_D, :), uAF_N(iAF_T, :), uAF_N(iAF_Ye,:), &
               uAF_N(iAF_E,:), duAFdT_N(iAF_E,:), duAFdYe_N(iAF_E,:) )

      CALL SetFVEC_EmissionAbsorption( dt ) ! --- Equations Vector

      CALL SetFJAC_EmissionAbsorption( dt ) ! --- Jacobian Matrix

      CALL SolveLinearSystems_FRC( Converged )

      MaxNorm = 0.0_DP
      DO iX = 1, nNodesX_G

        U_FRC(:,iX,iNew) = U_FRC(:,iX,iNew) + dU_FRC(:,iX)

        Norm = ENORM( dU_FRC(:,iX) / ( U_FRC(:,iX,iNew) + Floor_FRC(:) ) )
        IF( Norm <= NewtonTol ) Converged(iX) = .TRUE.

        IF( Norm >= MaxNorm )THEN

          MaxNorm = Norm
          iX_MAX  = iX

        END IF

      END DO

      CALL GetStates_FRC &
             ( uPR_N(:,iPR_D,:), uPF_N(iPF_Ne,:), uPF_N(iPF_E,:), iNew )

      IF( ALL( Converged ) )THEN

        CALL ComputeThermodynamicStates_Auxiliary &
               ( uPF_N(iPF_D,:), uPF_N(iPF_E,:), uPF_N(iPF_Ne,:), &
                 uAF_N(iAF_T,:), uAF_N(iAF_E,:), uAF_N(iAF_Ye,:) )

      END IF

      IF( MOD( iter, 10 ) == 0 )THEN

        WRITE(*,*)
        WRITE(*,'(A8,A)') ' ', 'Emission/Absorption'
        WRITE(*,*)
        WRITE(*,'(A10,A12,I6.6,A2,A11,ES10.4E2,A2,A9,I6.6)') &
          ' ', 'Iteration = ', iter, &
          ', ', '||dU/U|| = ', MaxNorm, &
          ', ', 'iX_MAX = ', iX_MAX
        WRITE(*,*)
        WRITE(*,'(A12,A4,ES10.4E2,A2,A4,ES10.4E2,A2,A4,ES10.4E2)') &
          '', 'D = ', uPF_N(iPF_D, iX_MAX) / ( Gram / Centimeter**3 ), &
          '', 'T = ', uAF_N(iAF_T, iX_MAX) / Kelvin, &
          '', 'Y = ', uAF_N(iAF_Ye,iX_MAX)
        WRITE(*,*)

      END IF

    END DO

    CALL UpdateNumberFlux_EmissionAbsorption( dt )

  END SUBROUTINE CoupleFluidRadiation_EmissionAbsorption


  SUBROUTINE SetRates_EmissionAbsorption

    INTEGER :: iX

    ASSOCIATE &
      ( D_N  => uPF_N(iPF_D, 1:nNodesX_G), &
        T_N  => uAF_N(iAF_T, 1:nNodesX_G), &
        Ye_N => uAF_N(iAF_Ye,1:nNodesX_G) )

    DO iX = 1, nNodesX_G

      CALL ComputeAbsorptionCoefficients &
             ( E_N, [ D_N(iX) ], [ T_N(iX) ], [ Ye_N(iX) ], &
               Chi(:,iX), dChidT(:,iX), dChidYe(:,iX) )

    END DO

    END ASSOCIATE ! D_N, etc.

  END SUBROUTINE SetRates_EmissionAbsorption


  SUBROUTINE SetEquilibrium_EmissionAbsorption

    INTEGER  :: iX
    REAL(DP) :: Mnu, dMnudT, dMnudYe

    CALL ComputeElectronChemicalPotential &
           ( uPF_N(iPF_D, :), uAF_N(iAF_T, :), uAF_N(iAF_Ye,:), &
             uAF_N(iAF_Me,:), duAFdT_N(iAF_Me,:), duAFdYe_N(iAF_Me,:) )

    CALL ComputeProtonChemicalPotential &
           ( uPF_N(iPF_D, :), uAF_N(iAF_T, :), uAF_N(iAF_Ye,:), &
             uAF_N(iAF_Mp,:), duAFdT_N(iAF_Mp,:), duAFdYe_N(iAF_Mp,:) )

    CALL ComputeNeutronChemicalPotential &
           ( uPF_N(iPF_D, :), uAF_N(iAF_T, :), uAF_N(iAF_Ye,:), &
             uAF_N(iAF_Mn,:), duAFdT_N(iAF_Mn,:), duAFdYe_N(iAF_Mn,:) )

    DO iX = 1, nNodesX_G

      Mnu     = uAF_N(iAF_Me,iX) + uAF_N(iAF_Mp,iX) &
                  - uAF_N(iAF_Mn,iX)
      dMnudT  = duAFdT_N(iAF_Me,iX) + duAFdT_N(iAF_Mp,iX) &
                  - duAFdT_N(iAF_Mn,iX)
      dMnudYe = duAFdYe_N(iAF_Me,iX) + duAFdYe_N(iAF_Mp,iX) &
                  - duAFdYe_N(iAF_Mn,iX)

      f_FD(:,iX) &
        = FermiDirac &
            ( E_N, Mnu, BoltzmannConstant * uAF_N(iAF_T,iX) )

      df_FDdT(:,iX) &
        = dFermiDiracdT &
            ( E_N, Mnu, BoltzmannConstant * uAF_N(iAF_T,iX), dMnudT, &
              uAF_N(iAF_T,iX) )

      df_FDdYe(:,iX) &
        = dFermiDiracdYe &
            ( E_N, Mnu, BoltzmannConstant * uAF_N(iAF_T,iX), dMnudYe )

    END DO

  END SUBROUTINE SetEquilibrium_EmissionAbsorption


  SUBROUTINE SetFVEC_EmissionAbsorption( dt )

    REAL(DP), INTENT(in) :: dt

    INTEGER :: iX, iE

    ASSOCIATE &
      ( hc3 => ( PlanckConstant * SpeedOfLight )**3 )

    FVEC_FRC = 0.0_DP

    DO iX = 1, nNodesX_G

      DO iE = 1, nNodesE_G

        FVEC_FRC(iE,iX) &
          = ( 1.0_DP + dt * Chi(iE,iX) ) * U_FRC(iE,iX,iNew) &
              - ( U_FRC(iE,iX,iOld) &
                    + dt * Chi(iE,iX) * FourPi * f_FD(iE,iX) )

      END DO

      FVEC_FRC(iFRC_Ne,iX) &
        = ( U_FRC(iFRC_Ne,iX,iNew) - U_FRC(iFRC_Ne,iX,iOld) )

      FVEC_FRC(iFRC_Ne,iX) &
        = FVEC_FRC(iFRC_Ne,iX) &
          + dt / hc3 &
            * DOT_PRODUCT &
              ( W2_N, Chi(:,iX) * ( FourPi * f_FD(:,iX) &
                                    - U_FRC(1:nNodesE_G,iX,iNew) ) )

      FVEC_FRC(iFRC_E,iX) &
        = ( U_FRC(iFRC_E, iX,iNew) - U_FRC(iFRC_E, iX,iOld) )

      FVEC_FRC(iFRC_E,iX) &
        = FVEC_FRC(iFRC_E,iX) &
          + dt / hc3 &
            * DOT_PRODUCT &
              ( W3_N, Chi(:,iX) * ( FourPi * f_FD(:,iX) &
                                    - U_FRC(1:nNodesE_G,iX,iNew) ) )

    END DO

    END ASSOCIATE ! hc3

  END SUBROUTINE SetFVEC_EmissionAbsorption


  SUBROUTINE SetFJAC_EmissionAbsorption( dt )

    REAL(DP), INTENT(in) :: dt

    INTEGER :: iX, iE

    ASSOCIATE &
      ( hc3 => ( PlanckConstant * SpeedOfLight )**3 )

    ASSOCIATE &
      ( mB   => BaryonMass, &
        D    => uPF_N(iPF_D,:), &
        dEdT => duAFdT_N(iAF_E,:) )

    FJAC_FRC = 0.0_DP

    DO iX = 1, nNodesX_G

      ! --- Neutrino Number Equation ---

      DO iE = 1, nNodesE_G

        ! --- Diagonal Components:
        FJAC_FRC(iE,iE,iX) &
          = 1.0_DP + dt * Chi(iE,iX)

        ! --- Derivative wrt Electron Density:
        FJAC_FRC(iE,iFRC_Ne,iX) &
          = dt * ( ( U_FRC(iE,iX,iNew) &
                     - FourPi * f_FD(iE,iX) ) * dChidYe(iE,iX) &
                   - Chi(iE,iX) * FourPi * df_FDdYe(iE,iX) ) &
               * ( mB / D(iX) )

        ! --- Derivative wrt Internal Energy Density:
        FJAC_FRC(iE,iFRC_E, iX) &
          = dt * ( ( U_FRC(iE,iX,iNew) &
                     - FourPi * f_FD(iE,iX) ) * dChidT(iE,iX) &
                   - Chi(iE,iX) * FourPi * df_FDdT(iE,iX) ) &
               / ( dEdT(iX) * D(iX) )

      END DO

      ! --- Electron Density Equation ---

      FJAC_FRC(iFRC_Ne,1:nNodesE_G,iX) &
        = - dt / hc3 * W2_N(:) * Chi(:,iX)

      FJAC_FRC(iFRC_Ne,iFRC_Ne,iX) &
        = 1.0_DP

      FJAC_FRC(iFRC_Ne,iFRC_Ne,iX) &
        = FJAC_FRC(iFRC_Ne,iFRC_Ne,iX) &
            + dt / hc3 &
              * DOT_PRODUCT &
                ( W2_N, dChidYe(:,iX) &
                        * ( FourPi * f_FD(:,iX) &
                            - U_FRC(1:nNodesE_G,iX,iNew) ) &
                        + Chi(:,iX) * FourPi * df_FDdYe(:,iX) ) &
              * ( mB / D(iX) )

      FJAC_FRC(iFRC_Ne,iFRC_E, iX) &
        = dt / hc3 &
          * DOT_PRODUCT &
            ( W2_N, dChidT(:,iX) &
                    * ( FourPi * f_FD(:,iX) &
                        - U_FRC(1:nNodesE_G,iX,iNew) ) &
                    + Chi(:,iX) * FourPi * df_FDdT(:,iX) ) &
          / ( dEdT(iX) * D(ix) )

      ! --- Internal Energy Equation ---

      FJAC_FRC(iFRC_E,1:nNodesE_G,iX) &
        = - dt / hc3 * W3_N(:) * Chi(:,iX)

      FJAC_FRC(iFRC_E, iFRC_Ne,iX) &
        = dt / hc3 &
          * DOT_PRODUCT &
            ( W3_N, dChidYe(:,iX) &
                    * ( FourPi * f_FD(:,iX) &
                        - U_FRC(1:nNodesE_G,iX,iNew) ) &
                    + Chi(:,iX) * FourPi * df_FDdYe(:,iX) ) &
          * ( mB / D(iX) )

      FJAC_FRC(iFRC_E, iFRC_E, iX) &
        = 1.0_DP

      FJAC_FRC(iFRC_E, iFRC_E, iX) &
        = FJAC_FRC(iFRC_E, iFRC_E, iX) &
            + dt / hc3 &
              * DOT_PRODUCT &
                ( W3_N, dChidT(:,iX) &
                        * ( FourPi * f_FD(:,iX) &
                            - U_FRC(1:nNodesE_G,iX,iNew) ) &
                        + Chi(:,iX) * FourPi * df_FDdT(:,iX) ) &
              / ( dEdT(iX) * D(ix) )

    END DO

    END ASSOCIATE ! mB, etc.

    END ASSOCIATE ! hc3

  END SUBROUTINE SetFJAC_EmissionAbsorption


  SUBROUTINE SolveLinearSystems_FRC( Converged )

    LOGICAL, DIMENSION(:), INTENT(in) :: Converged

    INTEGER :: &
      iX, LDA, INFO, iE
    INTEGER,  DIMENSION(nNodesE_G+2) :: &
      IPIV
    REAL(DP), DIMENSION(nNodesE_G+2) :: &
      b
    REAL(DP), DIMENSION(nNodesE_G+2,nNodesE_G+2) :: &
      A

    LDA = nNodesE_G + 2

    dU_FRC = 0.0_DP

    DO iX = 1, nNodesX_G

      IF( Converged(iX) ) CYCLE

      b(:)   = - FVEC_FRC(:,iX)
      A(:,:) =   FJAC_FRC(:,:,iX)

      CALL DGESV( LDA, 1, A, LDA, IPIV, b, LDA, INFO )

      dU_FRC(:,iX) = b(:)

    END DO

  END SUBROUTINE SolveLinearSystems_FRC


  PURE ELEMENTAL REAL(DP) FUNCTION FermiDirac( E, Mu, kT )

    REAL(DP), INTENT(in) :: E, Mu, kT

    FermiDirac = 1.0_DP / ( EXP( ( E - Mu ) / kT ) + 1.0_DP )

    RETURN
  END FUNCTION FermiDirac


  PURE ELEMENTAL REAL(DP) FUNCTION dFermiDiracdT( E, Mu, kT, dMudT, T )

    REAL(DP), INTENT(in) :: E, Mu, kT, dMudT, T

    dFermiDiracdT &
      = FermiDirac( E, Mu, kT )**2 * EXP( ( E - Mu ) / kT ) &
          * ( dMudT + ( E - Mu ) / T ) / kT

    RETURN
  END FUNCTION dFermiDiracdT


  PURE ELEMENTAL REAL(DP) FUNCTION dFermiDiracdYe( E, Mu, kT, dMudYe )

    REAL(DP), INTENT(in) :: E, Mu, kT, dMudYe

    dFermiDiracdYe &
      = FermiDirac( E, Mu, kT )**2 * EXP( ( E - Mu ) / kT ) * dMudYe / kT

    RETURN
  END FUNCTION dFermiDiracdYe


  SUBROUTINE SetStates_FRC( N, Ne, E, iState )

    REAL(DP), DIMENSION(nNodesE_G,nNodesX_G), INTENT(in) :: N
    REAL(DP), DIMENSION          (nNodesX_G), INTENT(in) :: Ne, E
    INTEGER,                                  INTENT(in) :: iState

    INTEGER :: iX

    DO iX = 1, nNodesX_G

      U_FRC(:,iX,iState) = State_FRC( N(:,iX), Ne(iX), E(iX) )

    END DO

  END SUBROUTINE SetStates_FRC


  PURE FUNCTION State_FRC( N, Ne, E )

    REAL(DP), DIMENSION(nNodesE_G+2)           :: State_FRC
    REAL(DP), DIMENSION(nNodesE_G), INTENT(in) :: N
    REAL(DP),                       INTENT(in) :: Ne, E

    State_FRC = [ N(:), Ne, E ]

    RETURN
  END FUNCTION State_FRC


  SUBROUTINE GetStates_FRC( N, Ne, E, iState )

    REAL(DP), DIMENSION(nNodesE_G,nNodesX_G), INTENT(out) :: N
    REAL(DP), DIMENSION          (nNodesX_G), INTENT(out) :: Ne, E
    INTEGER,                                  INTENT(in)  :: iState

    INTEGER :: iX, iE

    DO iX = 1, nNodesX_G
      DO iE = 1, nNodesE_G
        N(iE,iX) = U_FRC(iE,iX,iState)
      END DO
      Ne(iX) = U_FRC(iFRC_Ne,iX,iState)
      E (iX) = U_FRC(iFRC_E, iX,iState)
    END DO

  END SUBROUTINE GetStates_FRC


  SUBROUTINE UpdateNumberFlux_EmissionAbsorption( dt )

    REAL(DP), INTENT(in) :: dt

    INTEGER :: iX, iE

    DO iX = 1, nNodesX_G
      DO iE = 1, nNodesE_G

        uPR_N(iE,iPR_I1,iX) &
          = uPR_N(iE,iPR_I1,iX) / ( 1.0_DP + dt * Chi(iE,iX) )

        uPR_N(iE,iPR_I2,iX) &
          = uPR_N(iE,iPR_I2,iX) / ( 1.0_DP + dt * Chi(iE,iX) )

        uPR_N(iE,iPR_I3,iX) &
          = uPR_N(iE,iPR_I3,iX) / ( 1.0_DP + dt * Chi(iE,iX) )

      END DO
    END DO

  END SUBROUTINE UpdateNumberFlux_EmissionAbsorption


  PURE REAL(DP) FUNCTION ENORM( X )

    REAL(DP), DIMENSION(:), INTENT(in) :: X

    ENORM = SQRT( DOT_PRODUCT( X, X ) )

    RETURN
  END FUNCTION ENORM


END MODULE FluidRadiationCouplingSolutionModule_Implicit
