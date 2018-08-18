MODULE TwoMoment_DiscretizationModule_Collisions_Neutrinos

  USE KindModule, ONLY: &
    DP, Half, One
  USE UnitsModule, ONLY: &
    SpeedOfLight, &
    PlanckConstant, &
    BoltzmannConstant, &
    AtomicMassUnit
  USE ProgramHeaderModule, ONLY: &
    nNodesE, &
    nNodesX, &
    nNodesZ
  USE ReferenceElementModuleE, ONLY: &
    WeightsE
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX3D
  USE ReferenceElementModule, ONLY: &
    NodeNumberTable4D
  USE MeshModule, ONLY: &
    MeshE, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    nGF, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33
  USE FluidFieldsModule, ONLY: &
    nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne, &
    nAF, iAF_T, iAF_Ye, iAF_E
  USE RadiationFieldsModule, ONLY: &
    nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3, &
    nSpecies
  USE EquationOfStateModule_TABLE, ONLY: &
    ComputeThermodynamicStates_Auxiliary_TABLE, &
    ComputeThermodynamicStates_Primitive_TABLE, &
    ComputeTemperatureFromSpecificInternalEnergy_TABLE, &
    ComputeElectronChemicalPotential_TABLE, &
    ComputeProtonChemicalPotential_TABLE, &
    ComputeNeutronChemicalPotential_TABLE, &
    ComputeSpecificInternalEnergy_TABLE
  USE NeutrinoOpacitiesComputationModule, ONLY: &
    ComputeNeutrinoOpacities_EC_Point, &
    ComputeNeutrinoOpacities_ES_Point, &
    FermiDirac, &
    dFermiDiracdT, &
    dFermiDiracdY

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  PUBLIC :: ComputeIncrement_TwoMoment_Implicit

  LOGICAL, PARAMETER :: ReportConvergenceData = .TRUE.
  INTEGER  :: Iterations_Min
  INTEGER  :: Iterations_Max
  INTEGER  :: Iterations_Ave

  INTEGER  :: nE_G, nX_G
  INTEGER  :: iE_B0,    iE_E0
  INTEGER  :: iE_B1,    iE_E1
  INTEGER  :: iX_B0(3), iX_E0(3)
  INTEGER  :: iX_B1(3), iX_E1(3)
  REAL(DP) :: wTime
  REAL(DP), ALLOCATABLE :: E_N(:)        ! --- Energy Grid
  REAL(DP), ALLOCATABLE :: W2_N(:)       ! --- Ingegration Weights (E^2)
  REAL(DP), ALLOCATABLE :: W3_N(:)       ! --- Integration Weights (E^3)
  REAL(DP), ALLOCATABLE :: GX_N(:,:)     ! --- Spatial Geometry
  REAL(DP), ALLOCATABLE :: CF_N(:,:)     ! --- Conserved Fluid
  REAL(DP), ALLOCATABLE :: dF_N(:,:)     ! --- Conserved Fluid Increment
  REAL(DP), ALLOCATABLE :: CR_N(:,:,:,:) ! --- Conserved Radiation
  REAL(DP), ALLOCATABLE :: dR_N(:,:,:,:) ! --- Conserved Radiation Increment

CONTAINS


  SUBROUTINE ComputeIncrement_TwoMoment_Implicit &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, dt, GE, GX, U_F, dU_F, U_R, dU_R )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)    :: &
      dt
    REAL(DP), INTENT(in)    :: &
      GE(1:,iZ_B1(1):,1:)
    REAL(DP), INTENT(in)    :: &
      GX(1:,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:)
    REAL(DP), INTENT(in)    :: &
      U_F (1:,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:)
    REAL(DP), INTENT(inout) :: &
      dU_F(1:,iZ_B0(2):,iZ_B0(3):,iZ_B0(4):,1:)
    REAL(DP), INTENT(in)    :: &
      U_R (1:,iZ_B1(1):,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:,1:)
    REAL(DP), INTENT(inout) :: &
      dU_R(1:,iZ_B0(1):,iZ_B0(2):,iZ_B0(3):,iZ_B0(4):,1:,1:)

    INTEGER               :: iCR, iS, iCF, iGF, iX_G, iE_G
    REAL(DP)              :: PF_N(1:nPF)
    REAL(DP)              :: AF_N(1:nAF)
    REAL(DP), ALLOCATABLE :: Chi(:,:) ! --- Absorption Opacity
    REAL(DP), ALLOCATABLE :: Sig(:,:) ! --- Scattering Opacity
    REAL(DP), ALLOCATABLE :: fEQ(:,:) ! --- Equilibrium Distribution

    iE_B0 = iZ_B0(1);   iE_E0 = iZ_E0(1)
    iE_B1 = iZ_B1(1);   iE_E1 = iZ_E1(1)
    iX_B0 = iZ_B0(2:4); iX_E0 = iZ_E0(2:4)
    iX_B1 = iZ_B1(2:4); iX_E1 = iZ_E1(2:4)

    wTime = MPI_WTIME( )

    CALL InitializeCollisions( iZ_B0, iZ_E0 )

    wTime = MPI_WTIME( ) - wTime

!    WRITE(*,'(A4,A32,ES10.4E2)') '', 'InitializeCollisions: ', wTime

    ! --- Energy and Integration Weights ---

    wTime = MPI_WTIME( )

    CALL ComputePointsAndWeightsE( E_N, W2_N, W3_N )

    wTime = MPI_WTIME( ) - wTime

!    WRITE(*,'(A4,A32,ES10.4E2)') '', 'ComputePointsAndWeightsE: ', wTime

    ! --- Map Spatial Geometry Data ---

    DO iGF = iGF_Gm_dd_11, iGF_Gm_dd_33

      CALL MapForward_GX &
             ( iX_B0, iX_E0, &
               GX(:,iX_B0(1):iX_E0(1), &
                    iX_B0(2):iX_E0(2), &
                    iX_B0(3):iX_E0(3),iGF), &
               GX_N(iGF,1:nX_G) )

    END DO

    ! --- Map Fluid Data for Collision Update ---

    DO iCF = 1, nCF

      CALL MapForward_F &
             ( iX_B0, iX_E0, &
               U_F(:,iX_B0(1):iX_E0(1), &
                     iX_B0(2):iX_E0(2), &
                     iX_B0(3):iX_E0(3),iCF), &
               CF_N(iCF,1:nX_G) )

    END DO

    ! --- Map Radiation Data for Collision Update ---

    wTime = MPI_WTIME( )

    DO iS = 1, nSpecies
      DO iCR = 1, nCR

        CALL MapForward_R &
               ( iZ_B0, iZ_E0, &
                 U_R(:,iZ_B0(1):iZ_E0(1), &
                       iZ_B0(2):iZ_E0(2), &
                       iZ_B0(3):iZ_E0(3), &
                       iZ_B0(4):iZ_E0(4),iCR,iS), &
                 CR_N(1:nE_G,iCR,iS,1:nX_G) )

      END DO
    END DO

    wTime = MPI_WTIME( ) - wTime

!    WRITE(*,'(A4,A32,ES10.4E2)') '', 'MapForward_R: ', wTime

    ! --- Allocate Local Opacities ---

    ALLOCATE( Chi(nE_G,nSpecies), Sig(nE_G,nSpecies), fEQ(nE_G,nSpecies) )

    ! --- Implicit Update ---

    Iterations_Min = + HUGE( 1 )
    Iterations_Max = - HUGE( 1 )
    Iterations_Ave = 0

    DO iX_G = 1, nX_G

      dF_N(1:nCF,iX_G) = CF_N(1:nCF,iX_G)

      wTime = MPI_WTIME( )

      CALL ComputePrimitive_Euler &
             ( CF_N(iCF_D ,iX_G), CF_N(iCF_S1,iX_G), CF_N(iCF_S2,iX_G), &
               CF_N(iCF_S3,iX_G), CF_N(iCF_E ,iX_G), CF_N(iCF_Ne,iX_G), &
               PF_N(iPF_D ),      PF_N(iPF_V1),      PF_N(iPF_V2),      &
               PF_N(iPF_V3),      PF_N(iPF_E ),      PF_N(iPF_Ne),      &
               GX_N(iGF_Gm_dd_11,iX_G), &
               GX_N(iGF_Gm_dd_22,iX_G), &
               GX_N(iGF_Gm_dd_33,iX_G) )

      wTime = MPI_WTIME( ) - wTime

!      WRITE(*,'(A4,A32,ES10.4E2)') '', 'ComputePrimitive_Euler: ', wTime

      wTime = MPI_WTIME( )

      CALL ComputeThermodynamicStates_Auxiliary_TABLE &
             ( PF_N(iPF_D:iPF_D), PF_N(iPF_E:iPF_E), PF_N(iPF_Ne:iPF_Ne), &
               AF_N(iAF_T:iAF_T), AF_N(iAF_E:iAF_E), AF_N(iAF_Ye:iAF_Ye) )

      wTime = MPI_WTIME( ) - wTime

!      WRITE(*,'(A4,A32,ES10.4E2)') '', 'ComputeTS_Auxiliary: ', wTime

      wTime = MPI_WTIME( )

      ! --- Electron Capture Opacities ---

      CALL ComputeNeutrinoOpacities_EC_Point &
             ( 1, nE_G, E_N(1:nE_G), PF_N(iPF_D), AF_N(iAF_T), AF_N(iAF_Ye), &
               Chi(1:nE_G,1), iSpecies = 1 )

      ! --- Elastic Scattering Opacities ---

      CALL ComputeNeutrinoOpacities_ES_Point &
             ( 1, nE_G, E_N(1:nE_G), PF_N(iPF_D), AF_N(iAF_T), AF_N(iAF_Ye), &
               Sig(1:nE_G,1), iSpecies = 1 )

      wTime = MPI_WTIME( ) - wTime

!      WRITE(*,'(A4,A32,ES10.4E2)') '', 'ComputeNuOp_Point: ', wTime

      wTime = MPI_WTIME( )

      CALL SolveMatterEquations_EmAb &
             ( CR_N(:,iCR_N,1,iX_G), dt * Chi(:,1), fEQ(:,1), &
               PF_N(iPF_D), AF_N(iAF_T), AF_N(iAF_Ye), AF_N(iAF_E) )

      wTime = MPI_WTIME( ) - wTime

!      WRITE(*,'(A4,A32,ES10.4E2)') '', 'SolveMatterEquations_EmAb: ', wTime

      CALL ComputeThermodynamicStates_Primitive_TABLE &
             ( PF_N(iPF_D:iPF_D), AF_N(iAF_T:iAF_T), AF_N(iAF_Ye:iAF_Ye), &
               PF_N(iPF_E:iPF_E), AF_N(iAF_E:iAF_E), PF_N(iPF_Ne:iPF_Ne) )

      CALL ComputeConserved_Euler &
             ( PF_N(iPF_D ),      PF_N(iPF_V1),      PF_N(iPF_V2),      &
               PF_N(iPF_V3),      PF_N(iPF_E ),      PF_N(iPF_Ne),      &
               CF_N(iCF_D ,iX_G), CF_N(iCF_S1,iX_G), CF_N(iCF_S2,iX_G), &
               CF_N(iCF_S3,iX_G), CF_N(iCF_E ,iX_G), CF_N(iCF_Ne,iX_G), &
               GX_N(iGF_Gm_dd_11,iX_G), &
               GX_N(iGF_Gm_dd_22,iX_G), &
               GX_N(iGF_Gm_dd_33,iX_G) )

      ! --- Conserved Fluid Increment ---

      dF_N(1:nCF,iX_G) = ( CF_N(1:nCF,iX_G) - dF_N(1:nCF,iX_G) ) / dt

      ! --- Update Radiation Fields ---

      DO iS = 1, nSpecies

        ! --- Number Density ---

        CR_N(:,iCR_N,iS,iX_G) &
          = ( dt * Chi(:,iS) * fEQ(:,iS) &
              + CR_N(:,iCR_N,iS,iX_G) ) / ( One + dt * Chi(:,iS) )

        ! --- Number Flux (1) ---

        CR_N(:,iCR_G1,iS,iX_G) &
          = CR_N(:,iCR_G1,iS,iX_G) &
              / ( One + dt * ( Chi(:,iS) + Sig(:,iS) ) )

        ! --- Number Flux (2) ---

        CR_N(:,iCR_G2,iS,iX_G) &
          = CR_N(:,iCR_G2,iS,iX_G) &
              / ( One + dt * ( Chi(:,iS) + Sig(:,iS) ) )

        ! --- Number Flux (3) ---

        CR_N(:,iCR_G3,iS,iX_G) &
          = CR_N(:,iCR_G3,iS,iX_G) &
              / ( One + dt * ( Chi(:,iS) + Sig(:,iS) ) )

        ! --- Increments ---

        dR_N(:,iCR_N,iS,iX_G) &
          = Chi(:,iS) * ( fEQ(:,iS) - CR_N(:,iCR_N,iS,iX_G) )

        dR_N(:,iCR_G1,iS,iX_G) &
          = - ( Chi(:,iS) + Sig(:,iS) ) * CR_N(:,iCR_G1,iS,iX_G)

        dR_N(:,iCR_G2,iS,iX_G) &
          = - ( Chi(:,iS) + Sig(:,iS) ) * CR_N(:,iCR_G2,iS,iX_G)

        dR_N(:,iCR_G3,iS,iX_G) &
          = - ( Chi(:,iS) + Sig(:,iS) ) * CR_N(:,iCR_G3,iS,iX_G)

      END DO ! iS

    END DO ! iX_G

    IF( ReportConvergenceData )THEN

      WRITE(*,*)
      WRITE(*,'(A4,A)') &
        '', 'Convergence Data:'
      WRITE(*,*)
      WRITE(*,'(A6,A18,I4.4)') &
        '', 'Iterations (Min): ', Iterations_Min
      WRITE(*,'(A6,A18,I4.4)') &
        '', 'Iterations (Max): ', Iterations_Max
      WRITE(*,'(A6,A18,I4.4)') &
        '', 'Iterations (Ave): ', Iterations_Ave / nX_G
      WRITE(*,*)

    END IF

    ! --- Deallocate Local Opacities

    DEALLOCATE( Chi, Sig, fEQ )

    ! --- Map Increments Back ---

    ! --- Fluid ---

    DO iCF = 1, nCF

      CALL MapBackward_F &
             ( iX_B0, iX_E0, &
               dU_F(:,iX_B0(1):iX_E0(1), &
                      iX_B0(2):iX_E0(2), &
                      iX_B0(3):iX_E0(3),iCF), &
               dF_N(iCF,1:nX_G) )

    END DO

    ! --- Radiation ---

    DO iS = 1, nSpecies
      DO iCR = 1, nCR

        CALL MapBackward_R &
               ( iZ_B0, iZ_E0, &
                 dU_R(:,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
                        iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),iCR,iS), &
                 dR_N(1:nE_G,iCR,iS,1:nX_G) )

      END DO
    END DO

    CALL FinalizeCollisions

  END SUBROUTINE ComputeIncrement_TwoMoment_Implicit


  SUBROUTINE SolveMatterEquations_EmAb( J, Chi, J0, D, T, Y, E )

    REAL(DP), INTENT(in)    :: J  (1:nE_G)
    REAL(DP), INTENT(in)    :: Chi(1:nE_G)
    REAL(DP), INTENT(inout) :: J0 (1:nE_G)
    REAL(DP), INTENT(inout) :: D, T, Y, E

    INTEGER,  PARAMETER :: MaxIter = 20
    REAL(DP), PARAMETER :: Rtol = 1.0d-08
    REAL(DP), PARAMETER :: Utol = 1.0d-14

    LOGICAL  :: CONVERGED
    INTEGER  :: k
    REAL(DP) :: Yold, Eold, N_B
    REAL(DP) :: Theta2_N(1:nE_G), Theta3_N(1:nE_G)
    REAL(DP) :: Mnu, dMnudT, dMnudY
    REAL(DP) :: dEdT, dEdY
    REAL(DP) :: TMP(1), dTMPdT(1), dTMPdY(1)
    REAL(DP) :: dJ0dT_Y(1:nE_G), dJ0dY_T(1:nE_G)
    REAL(DP) :: dJ0dE_Y(1:nE_G), dJ0dY_E(1:nE_G)
    REAL(DP) :: U(2), dU(2), C(2), FVEC(2), FVEC0(2)
    REAL(DP) :: FJAC(2,2), IJAC(2,2), DJAC

    Yold = Y; Eold = E

    ! --- Auxiliary Variables ---

    N_B = D / AtomicMassUnit
    Theta2_N = W2_N * Chi / ( One + Chi )
    Theta3_N = W3_N * Chi / ( One + Chi )

    ! --- Neutrino Chemical Potential and Derivatives ---

    CALL ComputeNeutrinoChemicalPotentials( D, T, Y, Mnu, dMnudT, dMnudY )

    ! --- Equilibrium Distribution ---

    J0 = FermiDirac( E_N, Mnu, BoltzmannConstant * T )

    ! --- Initial Guess ---

    U(1:2) = [ Yold, Eold ]

    ! --- Old States (Constant) ---

    C(1) = DOT_PRODUCT( Theta2_N(:), J(:) ) + N_B * U(1)
    C(2) = DOT_PRODUCT( Theta3_N(:), J(:) ) + N_B * U(2)

    ! --- Electron Fraction Equation ---

    FVEC(1) = DOT_PRODUCT( Theta2_N(:), J0(:) ) + N_B * U(1) - C(1)

    ! --- Internal Energy Equation ---

    FVEC(2) = DOT_PRODUCT( Theta3_N(:), J0(:) ) + N_B * U(2) - C(2)

    ! --- Scale Equations and Save Initial Evaluation ---

    FVEC(:) = FVEC(:) / C(:); FVEC0(:) = FVEC(:);

    k = 0
    CONVERGED = .FALSE.
    DO WHILE( .NOT. CONVERGED )

      k = k + 1

      CALL ComputeSpecificInternalEnergy_TABLE &
             ( [ D ], [ T ], [ Y ], E = TMP, &
               dEdT_Option = dTMPdT, dEdY_Option = dTMPdY )

      dEdT = dTMPdT(1); dEdY = dTMPdY(1)

      ! --- Derivative of J0 wrt. T (Constant Y) ---

      dJ0dT_Y = dFermiDiracdT( E_N, Mnu, BoltzmannConstant * T, dMnudT, T )

      ! --- Derivative of J0 wrt. T (Constant T) ---

      dJ0dY_T = dFermiDiracdY( E_N, Mnu, BoltzmannConstant * T, dMnudY, T )

      ! --- Derivative of J0 wrt. E (Constant Y) ---

      dJ0dE_Y = dJ0dT_Y / dEdT

      ! --- Derivative of J0 wrt. Y (Constant E) ---

      dJ0dY_E = dJ0dY_T - dJ0dT_Y * dEdY / dEdT

      ! --- Jacobian ---

      FJAC(1,1) = DOT_PRODUCT( Theta2_N(:), dJ0dY_E(:) ) + N_B

      FJAC(1,2) = DOT_PRODUCT( Theta2_N(:), dJ0dE_Y(:) )

      FJAC(2,1) = DOT_PRODUCT( Theta3_N(:), dJ0dY_E(:) )

      FJAC(2,2) = DOT_PRODUCT( Theta3_N(:), dJ0dE_Y(:) ) + N_B

      ! --- Scale Jacobian ---

      FJAC(:,1) = FJAC(:,1) / C(:)

      FJAC(:,2) = FJAC(:,2) / C(:)

      ! --- Determinant of Jacobian ---

      DJAC = FJAC(1,1) * FJAC(2,2) - FJAC(2,1) * FJAC(1,2)

      ! --- Invert Jacobian ---

      IJAC(1,1) =   FJAC(2,2) / DJAC
      IJAC(2,1) = - FJAC(2,1) / DJAC
      IJAC(1,2) = - FJAC(1,2) / DJAC
      IJAC(2,2) =   FJAC(1,1) / DJAC

      ! --- Correction ---

      dU = - MATMUL( IJAC, FVEC )

      ! --- Apply Correction ---

      U = U + dU

      Y = U(1); E = U(2)

      CALL ComputeTemperatureFromSpecificInternalEnergy_TABLE &
             ( [ D ], [ E ], [ Y ], TMP ); T = TMP(1)

      ! --- Neutrino Chemical Potential and Derivatives ---

      CALL ComputeNeutrinoChemicalPotentials( D, T, Y, Mnu, dMnudT, dMnudY )

      ! --- Equilibrium Distribution ---

      J0 = FermiDirac( E_N, Mnu, BoltzmannConstant * T )

      ! --- Electron Fraction Equation ---

      FVEC(1) = DOT_PRODUCT( Theta2_N(:), J0(:) ) + N_B * U(1) - C(1)

      ! --- Internal Energy Equation ---

      FVEC(2) = DOT_PRODUCT( Theta3_N(:), J0(:) ) + N_B * U(2) - C(2)

      ! --- Scale Equations ---

      FVEC(:) = FVEC(:) / C(:)

      ! --- Check for Convergence ---

      IF( ENORM( FVEC ) < Rtol * ENORM( FVEC0 ) .OR. &
          ENORM( dU/U ) < Utol )THEN

        CONVERGED = .TRUE.

        Iterations_Min = MIN( Iterations_Min, k )
        Iterations_Max = MAX( Iterations_Max, k )
        Iterations_Ave = Iterations_Ave + k

      END IF

      IF( ( k .EQ. MaxIter ) .AND. ( .NOT. CONVERGED ) )THEN

        WRITE(*,*)
        WRITE(*,'(A4,A)') &
          '', 'SolveMatterEquations_EmAb:'
        WRITE(*,'(A6,A20,I4.4,A11)') &
          '', 'Did not converge in ', k, ' iterations'
        WRITE(*,'(A6,A)') &
          '', 'Exiting with unconverged result'
        WRITE(*,*)

      END IF

    END DO

  END SUBROUTINE SolveMatterEquations_EmAb


  SUBROUTINE ComputeNeutrinoChemicalPotentials( D, T, Y, M, dMdT, dMdY )

    REAL(DP), INTENT(in)  :: D, T, Y
    REAL(DP), INTENT(out) :: M, dMdT, dMdY

    REAL(DP) :: Me, dMedT, dMedY
    REAL(DP) :: Mp, dMpdT, dMpdY
    REAL(DP) :: Mn, dMndT, dMndY
    REAL(DP) :: TMP(1), dTMPdT(1), dTMPdY(1)

    ! --- Matter Chemical Potentials and Derivatives ---

    CALL ComputeElectronChemicalPotential_TABLE &
           ( [ D ], [ T ], [ Y ], M = TMP, &
             dMdT_Option = dTMPdT, dMdY_Option = dTMPdY )

    Me = TMP(1); dMedT = dTMPdT(1); dMedY = dTMPdY(1)

    CALL ComputeProtonChemicalPotential_TABLE &
           ( [ D ], [ T ], [ Y ], M = TMP, &
             dMdT_Option = dTMPdT, dMdY_Option = dTMPdY )

    Mp = TMP(1); dMpdT = dTMPdT(1); dMpdY = dTMPdY(1)

    CALL ComputeNeutronChemicalPotential_TABLE &
           ( [ D ], [ T ], [ Y ], M = TMP, &
             dMdT_Option = dTMPdT, dMdY_Option = dTMPdY )

    Mn = TMP(1); dMndT = dTMPdT(1); dMndY = dTMPdY(1)

    ! --- Neutrino Chemical Potential and Derivatives ---

    M = ( Me + Mp ) - Mn
    dMdT = ( dMedT + dMpdT ) - dMndT
    dMdY = ( dMedY + dMpdY ) - dMndY

  END SUBROUTINE ComputeNeutrinoChemicalPotentials


  PURE REAL(DP) FUNCTION ENORM( X )

    REAL(DP), DIMENSION(:), INTENT(in) :: X

    ENORM = SQRT( DOT_PRODUCT( X, X ) )

    RETURN
  END FUNCTION ENORM


  SUBROUTINE InitializeCollisions( iZ_B, iZ_E )

    INTEGER, INTENT(in) :: iZ_B(4), iZ_E(4)

    INTEGER :: nZ(4)

    nZ = iZ_E - iZ_B + 1

    nE_G = nZ(1) * nNodesZ(1)
    nX_G = PRODUCT( nZ(2:4) * nNodesZ(2:4) )

    ALLOCATE( E_N (nE_G) )
    ALLOCATE( W2_N(nE_G) )
    ALLOCATE( W3_N(nE_G) )
    ALLOCATE( GX_N(nGF,nX_G) )
    ALLOCATE( CF_N(nCF,nX_G) )
    ALLOCATE( dF_N(nCF,nX_G) )
    ALLOCATE( CR_N(nE_G,nCR,nSpecies,nX_G) )
    ALLOCATE( dR_N(nE_G,nCR,nSpecies,nX_G) )

  END SUBROUTINE InitializeCollisions


  SUBROUTINE FinalizeCollisions

    DEALLOCATE( E_N, W2_N, W3_N )
    DEALLOCATE( GX_N )
    DEALLOCATE( CF_N, dF_N )
    DEALLOCATE( CR_N, dR_N )

  END SUBROUTINE FinalizeCollisions


  SUBROUTINE ComputePointsAndWeightsE( E, W2, W3 )

    REAL(DP), INTENT(out) :: &
      E(:), W2(:), W3(:)

    INTEGER  :: iE_G, iE, iN
    REAL(DP) :: hc

    hc = PlanckConstant * SpeedOfLight

    ASSOCIATE( dE => MeshE % Width(iE_B0:iE_E0) )

    iE_G = 0
    DO iE = iE_B0, iE_E0
    DO iN = 1, nNodesE

      iE_G = iE_G + 1

      E (iE_G) = NodeCoordinate( MeshE, iE, iN )
      W2(iE_G) = WeightsE(iN) * ( dE(iE) / hc ) * ( E(iE_G) / hc )**2
      W3(iE_G) = W2(iE_G) * ( E(iE_G) / AtomicMassUnit )

    END DO
    END DO

    END ASSOCIATE ! -- dE

  END SUBROUTINE ComputePointsAndWeightsE


  SUBROUTINE MapForward_GX( iX_B, iX_E, GX, GX_N )

    INTEGER,  INTENT(in)  :: &
      iX_B(3), iX_E(3)
    REAL(DP), INTENT(in)  :: &
      GX(:,iX_B(1):,iX_B(2):,iX_B(3):)
    REAL(DP), INTENT(out) :: &
      GX_N(:)

    INTEGER :: iX1, iX2, iX3, iX
    INTEGER :: iN1, iN2, iN3, iN

    iX = 0
    DO iX3 = iX_B(3), iX_E(3)
    DO iX2 = iX_B(2), iX_E(2)
    DO iX1 = iX_B(1), iX_E(1)
    DO iN3 = 1, nNodesX(3)
    DO iN2 = 1, nNodesX(2)
    DO iN1 = 1, nNodesX(1)

      iX = iX + 1

      iN = NodeNumberTableX3D(iN1,iN2,iN3)

      GX_N(iX) = GX(iN,iX1,iX2,iX3)

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE MapForward_GX


  SUBROUTINE MapForward_F( iX_B, iX_E, FF, FF_N )

    INTEGER,  INTENT(in)  :: &
      iX_B(3), iX_E(3)
    REAL(DP), INTENT(in)  :: &
      FF(:,iX_B(1):,iX_B(2):,iX_B(3):)
    REAL(DP), INTENT(out) :: &
      FF_N(:)

    INTEGER :: iX1, iX2, iX3, iX
    INTEGER :: iN1, iN2, iN3, iN

    iX = 0
    DO iX3 = iX_B(3), iX_E(3)
    DO iX2 = iX_B(2), iX_E(2)
    DO iX1 = iX_B(1), iX_E(1)
    DO iN3 = 1, nNodesX(3)
    DO iN2 = 1, nNodesX(2)
    DO iN1 = 1, nNodesX(1)

      iX = iX + 1

      iN = NodeNumberTableX3D(iN1,iN2,iN3)

      FF_N(iX) = FF(iN,iX1,iX2,iX3)

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE MapForward_F


  SUBROUTINE MapBackward_F( iX_B, iX_E, FF, FF_N )

    INTEGER,  INTENT(in)  :: &
      iX_B(3), iX_E(3)
    REAL(DP), INTENT(out) :: &
      FF(:,iX_B(1):,iX_B(2):,iX_B(3):)
    REAL(DP), INTENT(in)  :: &
      FF_N(:)

    INTEGER :: iX1, iX2, iX3, iX
    INTEGER :: iN1, iN2, iN3, iN

    iX = 0
    DO iX3 = iX_B(3), iX_E(3)
    DO iX2 = iX_B(2), iX_E(2)
    DO iX1 = iX_B(1), iX_E(1)
    DO iN3 = 1, nNodesX(3)
    DO iN2 = 1, nNodesX(2)
    DO iN1 = 1, nNodesX(1)

      iX = iX + 1

      iN = NodeNumberTableX3D(iN1,iN2,iN3)

      FF(iN,iX1,iX2,iX3) = FF_N(iX)

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE MapBackward_F


  SUBROUTINE MapForward_R( iZ_B, iZ_E, RF, RF_N )

    INTEGER,  INTENT(in)  :: &
      iZ_B(4), iZ_E(4)
    REAL(DP), INTENT(in)  :: &
      RF(:,iZ_B(1):,iZ_B(2):,iZ_B(3):,iZ_B(4):)
    REAL(DP), INTENT(out) :: &
      RF_N(:,:)

    INTEGER :: iZ1, iZ2, iZ3, iZ4, iX, iE
    INTEGER :: iN1, iN2, iN3, iN4, iN

    iX = 0
    DO iZ4 = iZ_B(4), iZ_E(4)
    DO iZ3 = iZ_B(3), iZ_E(3)
    DO iZ2 = iZ_B(2), iZ_E(2)
    DO iN4 = 1, nNodesZ(4)
    DO iN3 = 1, nNodesZ(3)
    DO iN2 = 1, nNodesZ(2)

      iX = iX + 1

      iE = 0
      DO iZ1 = iZ_B(1), iZ_E(1)
      DO iN1 = 1, nNodesZ(1)

        iE = iE + 1

        iN = NodeNumberTable4D(iN1,iN2,iN3,iN4)

        RF_N(iE,iX) = RF(iN,iZ1,iZ2,iZ3,iZ4)

      END DO
      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE MapForward_R


  SUBROUTINE MapBackward_R( iZ_B, iZ_E, RF, RF_N )

    INTEGER,  INTENT(in)  :: &
      iZ_B(4), iZ_E(4)
    REAL(DP), INTENT(out) :: &
      RF(:,iZ_B(1):,iZ_B(2):,iZ_B(3):,iZ_B(4):)
    REAL(DP), INTENT(in)  :: &
      RF_N(:,:)

    INTEGER :: iZ1, iZ2, iZ3, iZ4, iX, iE
    INTEGER :: iN1, iN2, iN3, iN4, iN

    iX = 0
    DO iZ4 = iZ_B(4), iZ_E(4)
    DO iZ3 = iZ_B(3), iZ_E(3)
    DO iZ2 = iZ_B(2), iZ_E(2)
    DO iN4 = 1, nNodesZ(4)
    DO iN3 = 1, nNodesZ(3)
    DO iN2 = 1, nNodesZ(2)

      iX = iX + 1

      iE = 0
      DO iZ1 = iZ_B(1), iZ_E(1)
      DO iN1 = 1, nNodesZ(1)

        iE = iE + 1

        iN = NodeNumberTable4D(iN1,iN2,iN3,iN4)

        RF(iN,iZ1,iZ2,iZ3,iZ4) = RF_N(iE,iX)

      END DO
      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE MapBackward_R


  SUBROUTINE ComputePrimitive_Euler &
               ( N, S_1, S_2, S_3, G, Ne, D, V_1, V_2, V_3, E, De, &
                 Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    REAL(DP), INTENT(in)  :: N, S_1, S_2, S_3, G, Ne
    REAL(DP), INTENT(out) :: D, V_1, V_2, V_3, E, De
    REAL(DP), INTENT(in)  :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    ! --- Three-Velocity: Index Up   ---
    ! --- Three-Momentum: Index Down ---

    D   = N
    V_1 = S_1 / ( Gm_dd_11 * N )
    V_2 = S_2 / ( Gm_dd_22 * N )
    V_3 = S_3 / ( Gm_dd_33 * N )
    E   = G - Half * ( S_1**2 / Gm_dd_11 &
                       + S_2**2 / Gm_dd_22 &
                       + S_3**2 / Gm_dd_33 ) / N
    De  = Ne

  END SUBROUTINE ComputePrimitive_Euler


  SUBROUTINE ComputeConserved_Euler &
               ( D, V_1, V_2, V_3, E, De, N, S_1, S_2, S_3, G, Ne, &
                 Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    REAL(DP), INTENT(in)  :: D, V_1, V_2, V_3, E, De
    REAL(DP), INTENT(out) :: N, S_1, S_2, S_3, G, Ne
    REAL(DP), INTENT(in)  :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    ! --- Three-Velocity: Index Up   ---
    ! --- Three-Momentum: Index Down ---

    N   = D
    S_1 = D * V_1
    S_2 = D * V_2
    S_3 = D * V_3
    G   = E + Half * D * ( Gm_dd_11 * V_1**2 &
                           + Gm_dd_22 * V_2**2 &
                           + Gm_dd_33 * V_3**2 )
    Ne  = De

  END SUBROUTINE ComputeConserved_Euler


END MODULE TwoMoment_DiscretizationModule_Collisions_Neutrinos
