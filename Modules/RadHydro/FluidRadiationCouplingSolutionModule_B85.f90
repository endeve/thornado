MODULE FluidRadiationCouplingSolutionModule_B85

  USE KindModule, ONLY: &
    DP, FourPi
  USE UnitsModule, ONLY: &
    BoltzmannConstant, &
    PlanckConstant, &
    SpeedOfLight
  USE ProgramHeaderModule, ONLY: &
    nX, nNodesX, nDOFX, &
    nE, nNodesE, nDOFE
  USE UtilitiesModule, ONLY: &
    WriteVector, &
    WriteMatrix
  USE FluidFieldsModule, ONLY: &
    uCF, nCF, &
    uPF, nPF, iPF_D, iPF_E, iPF_Ne, &
    uAF, nAF, iAF_P, iAF_T, iAF_E, iAF_Ye, iAF_S, &
    iAF_Me, iAF_Mp, iAF_Mn, iAF_Gm, iAF_Cs
  USE EulerEquationsUtilitiesModule, ONLY: &
    ComputeConserved, &
    ComputePrimitive
  USE RadiationFieldsModule, ONLY: &
    nPR, iPR_D, iPR_I1, iPR_I2, iPR_I3
  USE MomentEquationsUtilitiesModule, ONLY: &
    ComputeConservedMoments, &
    ComputePrimitiveMoments
  USE FluidRadiationCouplingUtilitiesModule, ONLY: &
    InitializeNodes, &
    InitializeNodesX, &
    InitializeWeights, &
    InitializeFluidFields, &
    InitializeRadiationFields, &
    FinalizeFluidFields, &
    FinalizeRadiationFields, &
    FermiDirac, &
    dFermiDiracdT, &
    dFermiDiracdY, &
    ENORM
  USE EquationOfStateModule, ONLY: &
    BaryonMass, &
    ComputeAuxiliary_Fluid, &
    ComputeSpecificInternalEnergy, &
    ComputeTemperatureFromSpecificInternalEnergy, &
    ComputeThermodynamicStates_Primitive, &
    ComputeElectronChemicalPotential, &
    ComputeProtonChemicalPotential, &
    ComputeNeutronChemicalPotential
  USE OpacityModule, ONLY: &
    ComputeAbsorptionOpacity, &
    ComputeScatteringOpacity_ES, &
    ComputeScatteringOpacity_NES

  IMPLICIT NONE
  PRIVATE

  LOGICAL            :: EvolveFluid
  INTEGER            :: nNodesX_G, nNodesE_G
  INTEGER            :: i_Y, i_E
  INTEGER, PARAMETER :: iOld = 1, iNew = 2
  REAL(DP), DIMENSION(:),     ALLOCATABLE :: E_N, W2_N, W3_N
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE :: X_N, uPF_N, uAF_N
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE :: Chi, FD
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE :: dFDdY_E
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE :: dFDdE_Y
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE :: Sigma
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: uPR_N
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: R0_in, R0_Out

  REAL(DP), DIMENSION(:,:),   ALLOCATABLE :: FVEC
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE :: DVEC
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: UVEC
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: FJAC

  PUBLIC :: CoupleFluidRadiation_B85

CONTAINS


  SUBROUTINE CoupleFluidRadiation_B85 &
               ( dt, iX_B0, iX_E0, iX_B1, iX_E1, U_F, dU_F, U_R, dU_R, &
                 EvolveFluid_Option )

    REAL(DP), INTENT(in)  :: &
      dt
    INTEGER,  INTENT(in)  :: &
      iX_B0(3), iX_B1(3), iX_E0(3), iX_E1(3)
    REAL(DP), INTENT(in)  :: &
      U_F (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(out) :: &
      dU_F(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)
    REAL(DP), INTENT(in)  :: &
      U_R (1:,1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:,1:)
    REAL(DP), INTENT(out) :: &
      dU_R(1:,1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:,1:)
    LOGICAL,  INTENT(in), OPTIONAL :: &
      EvolveFluid_Option

    INTEGER :: iX1, iX2, iX3

    EvolveFluid = .TRUE.
    IF( PRESENT( EvolveFluid_Option ) ) &
      EvolveFluid = EvolveFluid_Option

    CALL ComputePrimitiveMoments &
           ( iX_Begin = iX_B0, iX_End = iX_E0 )

    DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)

          CALL ComputePrimitive &
                 ( uCF(:,iX1,iX2,iX3,1:nCF), uPF(:,iX1,iX2,iX3,1:nPF) )

          CALL ComputeAuxiliary_Fluid &
                 ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_E ), &
                   uPF(:,iX1,iX2,iX3,iPF_Ne), uAF(:,iX1,iX2,iX3,iAF_P ), &
                   uAF(:,iX1,iX2,iX3,iAF_T ), uAF(:,iX1,iX2,iX3,iAF_Ye), &
                   uAF(:,iX1,iX2,iX3,iAF_S ), uAF(:,iX1,iX2,iX3,iAF_E ), &
                   uAF(:,iX1,iX2,iX3,iAF_Gm), uAF(:,iX1,iX2,iX3,iAF_Cs) )

        END DO
      END DO
    END DO

    CALL InitializeFluidRadiationCoupling

    CALL CoupleFluidRadiation( dt )

    CALL FinalizeFluidRadiationCoupling

    DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)

          CALL ComputeThermodynamicStates_Primitive &
                 ( uPF(1:nDOFX,iX1,iX2,iX3,iPF_D ), &
                   uAF(1:nDOFX,iX1,iX2,iX3,iAF_T ), &
                   uAF(1:nDOFX,iX1,iX2,iX3,iAF_Ye), &
                   uPF(1:nDOFX,iX1,iX2,iX3,iPF_E ), &
                   uAF(1:nDOFX,iX1,iX2,iX3,iAF_E ), &
                   uPF(1:nDOFX,iX1,iX2,iX3,iPF_Ne) )

          CALL ComputeConserved &
                 ( uPF(:,iX1,iX2,iX3,1:nPF), uCF(:,iX1,iX2,iX3,1:nCF) )

        END DO
      END DO
    END DO

    CALL ComputeConservedMoments &
           ( iX_Begin = iX_B0, iX_End = iX_E0 )

  END SUBROUTINE CoupleFluidRadiation_B85


  SUBROUTINE InitializeFluidRadiationCoupling

    nNodesX_G = PRODUCT(nX) * nDOFX
    nNodesE_G =         nE  * nDOFE

    ALLOCATE( E_N(nNodesE_G) )
    CALL InitializeNodes( E_N )

    ALLOCATE( W2_N(nNodesE_G), W3_N(nNodesE_G) )
    CALL InitializeWeights( W2_N, W3_N )

    ALLOCATE( X_N(nNodesX_G,3) )
    CALL InitializeNodesX( X_N )

    ALLOCATE( uPF_N(nPF, nNodesX_G) )
    ALLOCATE( uAF_N(nAF, nNodesX_G) )
    CALL InitializeFluidFields( uPF_N, uAF_N )

    ALLOCATE( uPR_N(nNodesE_G, nPR, nNodesX_G) )
    CALL InitializeRadiationFields( uPR_N )

    ALLOCATE &
      ( FD     (1:nNodesE_G,1:nNodesX_G), &
        dFDdY_E(1:nNodesE_G,1:nNodesX_G), &
        dFDdE_Y(1:nNodesE_G,1:nNodesX_G) )

    ALLOCATE &
      ( Chi   (1:nNodesE_G,1:nNodesX_G), &
        Sigma (1:nNodesE_G,1:nNodesX_G), &
        R0_In (1:nNodesE_G,1:nNodesE_G,1:nNodesX_G), &
        R0_Out(1:nNodesE_G,1:nNodesE_G,1:nNodesX_G) )

    i_Y = nNodesE_G + 1
    i_E = nNodesE_G + 2

    ALLOCATE &
      ( FVEC(nNodesE_G+2,nNodesX_G),   &
        DVEC(nNodesE_G+2,nNodesX_G),   &
        UVEC(nNodesE_G+2,nNodesX_G,2), &
        FJAC(nNodesE_G+2,nNodesE_G+2,nNodesX_G))

  END SUBROUTINE InitializeFluidRadiationCoupling


  SUBROUTINE FinalizeFluidRadiationCoupling

    CALL FinalizeFluidFields( uPF_N, uAF_N )

    CALL FinalizeRadiationFields( uPR_N )

    DEALLOCATE( E_N, W2_N, W3_N )
    DEALLOCATE( X_N, uPF_N, uAF_N, uPR_N )

    DEALLOCATE( FD, dFDdY_E, dFDdE_Y )

    DEALLOCATE( Chi, Sigma, R0_In, R0_Out )

    DEALLOCATE( FVEC, DVEC, UVEC, FJAC )

  END SUBROUTINE FinalizeFluidRadiationCoupling


  SUBROUTINE CoupleFluidRadiation( dt )

    REAL(DP), INTENT(in) :: dt

    LOGICAL,  DIMENSION(nNodesX_G) :: Converged
    INTEGER                        :: iter, iX, iE, iX_Max
    REAL(DP)                       :: Norm, MaxNorm, RES
    REAL(DP),            PARAMETER :: NewtonTol = 1.0d-10
    REAL(DP),            PARAMETER :: ResTol    = 1.0d-8
    REAL(DP), DIMENSION(nNodesX_G) :: RES_0
    REAL(DP), DIMENSION(nNodesE_G,nNodesX_G) :: absI

    absI = SQRT( uPR_N(1:nNodesE_G,iPR_I1,1:nNodesX_G)**2 &
                 + uPR_N(1:nNodesE_G,iPR_I2,1:nNodesX_G)**2 &
                 + uPR_N(1:nNodesE_G,iPR_I3,1:nNodesX_G)**2 )

    CALL SetUVEC( uPR_N(:,iPR_D,:), uAF_N(iAF_Ye,:), uAF_N(iAF_E,:), iOld )

    CALL SetRates

    Converged = .FALSE.
    iter      = 0

    DO WHILE ( .NOT. ALL( Converged ) )

      iter = iter + 1

      CALL SetUVEC( uPR_N(:,iPR_D,:), uAF_N(iAF_Ye,:), uAF_N(iAF_E,:), iNew )

      ! --- Compute FD Distribution and Derivatives ---

      CALL SetEquilibrium( Converged )

      ! --- Set Equation Vectors ---

      CALL SetFVEC( dt, Converged )

      ! --- Set Jacobian Matrix ---

      CALL SetFJAC( dt, Converged )

      ! --- Invert for Correction ---

      CALL SolveLinearSystems( Converged )

      MaxNorm = 0.0_DP
      DO iX = 1, nNodesX_G

        IF( iter == 1 ) &
          RES_0(iX) = ENORM( FVEC(:,iX) )

        DO iE = 1, nNodesE_G
          UVEC(iE,iX,iNew) &
            = MAX( UVEC(iE,iX,iNew) + DVEC(iE,iX), absI(iE,iX) )
        END DO
        UVEC(i_Y,iX,iNew) = UVEC(i_Y,iX,iNew) + DVEC(i_Y,iX)
        UVEC(i_E,iX,iNew) = UVEC(i_E,iX,iNew) + DVEC(i_E,iX)

        RES  = ENORM( FVEC(:,iX) )
        Norm = ENORM( DVEC(:,iX) / ( UVEC(:,iX,iNew) + EPSILON(1.0_DP) ) )
        IF( Norm <= NewtonTol .OR. RES < ResTol * RES_0(iX) ) &
          Converged(iX) = .TRUE.

        IF( Norm >= MaxNorm )THEN

          MaxNorm = Norm
          iX_Max  = iX

        END IF

      END DO

      CALL GetUVEC( uPR_N(:,iPR_D,:), uAF_N(iAF_Ye,:), uAF_N(iAF_E,:), iNew )

      CALL ComputeTemperatureFromSpecificInternalEnergy &
             ( uPF_N(iPF_D,:), uAF_N(iAF_E,:), uAF_N(iAF_Ye,:), uAF_N(iAF_T,:) )

      IF( iter > 9 )THEN

        WRITE(*,*)
        WRITE(*,'(A8,A)') ' ', 'CoupleFluidRadiation_B85'
        WRITE(*,*)
        WRITE(*,'(A10,A12,I6.6,A2,A11,ES10.4E2,A2,A9,I6.6)') &
          ' ', 'Iteration = ', iter, &
          ', ', '||dU/U|| = ', MaxNorm, &
          ', ', 'iX_MAX = ', iX_MAX
        WRITE(*,*)
        WRITE(*,'(A12,I4.4,A19,I4.4)') &
          '', COUNT( Converged .EQV. .TRUE. ), &
          ' converged, out of ', SIZE( Converged )
        WRITE(*,*)
        WRITE(*,'(A12,A22,ES10.4E2,A3,ES10.4E2)') &
          '', 'MIN/MAX    ||dN/N|| = ', &
          MINVAL &
            ( ABS( DVEC(1:nNodesE_G,iX_MAX) &
              / (UVEC(1:nNodesE_G,iX_MAX,iNew)+EPSILON(1.0_DP)) ) ),' / ', &
          MAXVAL &
            ( ABS( DVEC(1:nNodesE_G,iX_MAX) &
              / (UVEC(1:nNodesE_G,iX_MAX,iNew)+EPSILON(1.0_DP)) ) )

        WRITE(*,'(A12,A22,ES10.4E2,A3,ES10.4E2)') &
          '', 'MIN/MAX      ||dN|| = ', &
          MINVAL( ABS( DVEC(1:nNodesE_G,iX_MAX) ) ),' / ', &
          MAXVAL( ABS( DVEC(1:nNodesE_G,iX_MAX) ) )

        WRITE(*,'(A12,A22,ES10.4E2,A3,ES10.4E2)') &
          '', 'MIN/MAX       ||N|| = ', &
          MINVAL( ABS( UVEC(1:nNodesE_G,iX_MAX,iNew) ) ),' / ', &
          MAXVAL( ABS( UVEC(1:nNodesE_G,iX_MAX,iNew) ) )

        WRITE(*,'(A12,A22,ES10.4E2,A3,ES10.4E2)') &
          '', 'MIN/MAX ||FVEC(N)|| = ', &
          MINVAL( ABS( FVEC(1:nNodesE_G,iX_MAX) ) ), ' / ', &
          MAXVAL( ABS( FVEC(1:nNodesE_G,iX_MAX) ) )
        WRITE(*,*)
        WRITE(*,'(A12,A22,ES10.4E2)') &
          '', '           ||dY/Y|| = ', &
          ABS( DVEC(i_Y,iX_MAX) / UVEC(i_Y,iX_MAX,iNew) )
        WRITE(*,'(A12,A22,ES10.4E2)') &
          '', '        ||FVEC(Y)|| = ', ABS( FVEC(i_Y,iX_MAX) )
        WRITE(*,*)
        WRITE(*,'(A12,A22,ES10.4E2)') &
          '', '           ||dE/E|| = ', &
          ABS( DVEC(i_E,iX_MAX) / UVEC(i_E,iX_MAX,iNew) )
        WRITE(*,'(A12,A22,ES10.4E2)') &
          '', '        ||FVEC(E)|| = ', ABS( FVEC(i_E,iX_MAX) )
        WRITE(*,*)

        IF( iter > 99 )THEN

          CALL WriteVector &
                 ( nNodesE_G, E_N,    'EVEC.dat' )
          CALL WriteVector &
                 ( nNodesE_G, UVEC(1:nNodesE_G,iX_MAX,iNew), 'UVEC.dat' )
          PRINT*, 'absI = ', absI(:,iX_MAX)

          STOP
        END IF

      END IF

    END DO

    CALL UpdateNumberFlux( dt )

!!$    CALL WriteVector( nNodesE_G+2, FVEC(:,1), 'FVEC.dat' )
!!$
!!$    CALL WriteMatrix( nNodesE_G+2, nNodesE_G+2, FJAC(:,:,1), 'FJAC.dat' )

  END SUBROUTINE CoupleFluidRadiation


  SUBROUTINE SetUVEC( D, Y, E, iState )

    REAL(DP), DIMENSION(nNodesE_G,nNodesX_G), INTENT(in) :: D
    REAL(DP), DIMENSION(nNodesX_G),           INTENT(in) :: Y, E
    INTEGER,                                  INTENT(in) :: iState

    INTEGER :: iX, iE

    DO iX = 1, nNodesX_G
      DO iE = 1, nNodesE_G
        UVEC(iE,iX,iState) = D(iE,iX)
      END DO
      UVEC(i_Y,iX,iState) = Y(iX)
      UVEC(i_E,iX,iState) = E(iX)
    END DO

  END SUBROUTINE SetUVEC


  SUBROUTINE GetUVEC( D, Y, E, iState )

    REAL(DP), DIMENSION(nNodesE_G,nNodesX_G), INTENT(out) :: D
    REAL(DP), DIMENSION(nNodesX_G),           INTENT(out) :: Y, E
    INTEGER,                                  INTENT(in)  :: iState

    INTEGER :: iX, iE

    DO iX = 1, nNodesX_G
      DO iE = 1, nNodesE_G
        D(iE,iX) = UVEC(iE,iX,iState)
      END DO
      Y(iX) = UVEC(i_Y,iX,iState)
      E(iX) = UVEC(i_E,iX,iState)
    END DO

  END SUBROUTINE GetUVEC


  SUBROUTINE SetRates

    ASSOCIATE &
      ( kT => BoltzmannConstant * uAF_N(iAF_T,1:nNodesX_G) )
    ASSOCIATE &
      ( D_N   => uPF_N(iPF_D, 1:nNodesX_G), &
        T_N   => uAF_N(iAF_T, 1:nNodesX_G), &
        Y_N   => uAF_N(iAF_Ye,1:nNodesX_G), &
        Eta_N => uAF_N(iAF_Me,1:nNodesX_G) / kT )

    CALL ComputeAbsorptionOpacity &
           ( E_N, D_N, T_N, Y_N, X_N(:,1), X_N(:,2), X_N(:,3), Chi )

    CALL ComputeScatteringOpacity_ES &
           ( E_N, D_N, T_N, Y_N, X_N(:,1), X_N(:,2), X_N(:,3), Sigma )

!    CALL ComputeScatteringOpacity_NES &
!           ( E_N, T_N, Eta_N, R0_In, R0_Out )

    END ASSOCIATE ! D_N, etc.
    END ASSOCIATE ! kT

  END SUBROUTINE SetRates


  SUBROUTINE SetEquilibrium( Converged )

    LOGICAL, DIMENSION(nNodesX_G), INTENT(in) :: Converged

    INTEGER  :: iX
    REAL(DP) :: Mnu, dMnudT, dMnudY
    REAL(DP), DIMENSION(1)         :: TMP, dEdT, dEdY
    REAL(DP), DIMENSION(nNodesX_G) :: dMedT, dMedY
    REAL(DP), DIMENSION(nNodesX_G) :: dMpdT, dMpdY
    REAL(DP), DIMENSION(nNodesX_G) :: dMndT, dMndY
    REAL(DP), DIMENSION(nNodesE_G) :: dFDdT_Y
    REAL(DP), DIMENSION(nNodesE_G) :: dFDdY_T

    CALL ComputeElectronChemicalPotential &
           ( uPF_N(iPF_D, :), uAF_N(iAF_T, :), &
             uAF_N(iAF_Ye,:), uAF_N(iAF_Me,:), &
             dMdT_Option = dMedT(:), &
             dMdY_Option = dMedY(:) )

    CALL ComputeProtonChemicalPotential &
           ( uPF_N(iPF_D, :), uAF_N(iAF_T, :), &
             uAF_N(iAF_Ye,:), uAF_N(iAF_Mp,:), &
             dMdT_Option = dMpdT(:), &
             dMdY_Option = dMpdY(:) )

    CALL ComputeNeutronChemicalPotential &
           ( uPF_N(iPF_D, :), uAF_N(iAF_T, :), &
             uAF_N(iAF_Ye,:), uAF_N(iAF_Mn,:), &
             dMdT_Option = dMndT(:), &
             dMdY_Option = dMndY(:) )

    ASSOCIATE( kT => BoltzmannConstant * uAF_N(iAF_T,:) )

    DO iX = 1, nNodesX_G

      IF( Converged(iX) ) CYCLE

      Mnu    = uAF_N(iAF_Me,iX) + uAF_N(iAF_Mp,iX) - uAF_N(iAF_Mn,iX)
      dMnudT = dMedT(iX)        + dMpdT(iX)        - dMndT(iX)
      dMnudY = dMedY(iX)        + dMpdY(iX)        - dMndY(iX)

      CALL ComputeSpecificInternalEnergy &
             ( [ uPF_N(iPF_D, iX) ], [ uAF_N(iAF_T,iX) ], &
               [ uAF_N(iAF_Ye,iX) ], TMP, &
               dEdT_Option = dEdT, &
               dEdY_Option = dEdY )

      FD(:,iX) = FermiDirac( E_N, Mnu, kT(iX) )
      dFDdT_Y  = dFermiDiracdT( E_N, Mnu, kT(iX), dMnudT, uAF_N(iAF_T,iX) )
      dFDdY_T  = dFermiDiracdY( E_N, Mnu, kT(iX), dMnudY, uAF_N(iAF_T,iX) )

      dFDdY_E(:,iX) &
        = dFDdY_T(:) - dFDdT_Y(:) * dEdY(1) / dEdT(1)

      dFDdE_Y(:,iX) &
        = dFDdT_Y(:) / dEdT(1)

    END DO

    END ASSOCIATE ! kT

  END SUBROUTINE SetEquilibrium


  SUBROUTINE SetFVEC( dt, Converged )

    REAL(DP),                      INTENT(in) :: dt
    LOGICAL, DIMENSION(nNodesX_G), INTENT(in) :: Converged

    INTEGER :: iX

    ASSOCIATE &
      ( hc3 => ( PlanckConstant * SpeedOfLight )**3, &
        mB  => BaryonMass, &
        D_N => uPF_N(iPF_D,:) )

    DO iX = 1, nNodesX_G

      IF( Converged(iX) ) CYCLE

      ! --- Number Equation ---

      FVEC(1:nNodesE_G,iX) &
        = ( UVEC(1:nNodesE_G,iX,iNew) - UVEC(1:nNodesE_G,iX,iOld) ) &
            - dt * RHSVEC &
                     ( nNodesE_G, W2_N(:), UVEC(1:nNodesE_G,iX,iNew), &
                       FD(:,iX), Chi(:,iX), R0_In(:,:,iX), R0_Out(:,:,iX) )

      ! --- Electron Fraction Equation ---

      FVEC(i_Y,iX) &
        = DOT_PRODUCT &
            ( W2_N , ( UVEC(1:nNodesE_G,iX,iNew) &
                       - UVEC(1:nNodesE_G,iX,iOld) ) ) &
            * ( mB / D_N(iX) ) / hc3 &
          + ( UVEC(i_Y,iX,iNew) - UVEC(i_Y,iX,iOld) )

      ! --- Energy Equation ---

      FVEC(i_E,iX) &
        = DOT_PRODUCT &
            ( W3_N , ( UVEC(1:nNodesE_G,iX,iNew) &
                       - UVEC(1:nNodesE_G,iX,iOld) ) ) &
            * ( 1.0_DP / D_N(iX) ) / hc3 &
          + ( UVEC(i_E,iX,iNew) - UVEC(i_E,iX,iOld) )

    END DO

    END ASSOCIATE ! mB, D_N

  END SUBROUTINE SetFVEC


  SUBROUTINE SetFJAC( dt, Converged )

    REAL(DP),              INTENT(in) :: dt
    LOGICAL, DIMENSION(:), INTENT(in) :: Converged

    INTEGER :: iX, iE

    ASSOCIATE &
      ( hc3 => ( PlanckConstant * SpeedOfLight )**3, &
        mB  => BaryonMass, &
        D_N => uPF_N(iPF_D,:) )

    DO iX = 1, nNodesX_G

      FJAC(:,:,iX) = 0.0_DP

      IF( Converged(iX) ) CYCLE

      ! --- Number Equation ---

      ! ------ Diagonal Elements ------

      DO iE = 1, nNodesE_G

        FJAC(iE,iE,iX) &
          = 1.0_DP + dt * Chi(iE,iX) &
            + dt * SUM( W2_N(1:nNodesE_G) * R0_In(1:nNodesE_G,iE,iX) &
                        * UVEC(1:nNodesE_G,iX,iNew) ) &
            + dt * SUM( W2_N(1:nNodesE_G) * R0_Out(1:nNodesE_G,iE,iX) &
                        * ( FourPi - UVEC(1:nNodesE_G,iX,iNew) ) )

      END DO

      ! ------ Off-Diagonal Elements ------

      DO iE = 1, nNodesE_G

        FJAC(1:nNodesE_G,iE,iX) &
          = FJAC(1:nNodesE_G,iE,iX) &
            - dt * W2_N(iE)  &
              * ( UVEC(1:nNodesE_G,iX,iNew) *  R0_In(:,iE,iX) &
                  + ( FourPi - UVEC(1:nNodesE_G,iX,iNew) ) * R0_Out(:,iE,iX) )

      END DO

      FJAC(1:nNodesE_G,i_Y,iX) &
        = - dt * Chi(1:nNodesE_G,iX) * FourPi * dFDdY_E(1:nNodesE_G,iX)

      FJAC(1:nNodesE_G,i_E,iX) &
        = - dt * Chi(1:nNodesE_G,iX) * FourPi * dFDdE_Y(1:nNodesE_G,iX)

      ! --- Electron Fraction Equation ---

      FJAC(i_Y,1:nNodesE_G,iX) &
        = W2_N(:) * ( mB / D_N(iX) ) / hc3

      FJAC(i_Y,i_Y,iX) &
        = 1.0_DP

      FJAC(i_Y,i_E,iX) &
        = 0.0_DP

      ! --- Energy Equation ---

      FJAC(i_E,1:nNodesE_G,iX) &
        = W3_N(:) * ( 1.0_DP / D_N(iX) ) / hc3

      FJAC(i_E,i_Y,iX) &
        = 0.0_DP

      FJAC(i_E,i_E,iX) &
        = 1.0_DP

    END DO

    END ASSOCIATE ! hc3, etc.

  END SUBROUTINE SetFJAC


  SUBROUTINE SolveLinearSystems( Converged )

    LOGICAL, DIMENSION(nNodesX_G), INTENT(in) :: Converged

    INTEGER                         :: iX, INFO
    INTEGER, DIMENSION(nNodesE_G+2) :: IPIV

    DO iX = 1, nNodesX_G

      DVEC(:,iX) = 0.0_DP

      IF( Converged(iX) ) CYCLE

      DVEC(:,iX) = - FVEC(:,iX)

      CALL DGESV &
             ( nNodesE_G+2, 1, FJAC(:,:,iX), nNodesE_G+2, &
               IPIV, DVEC(:,iX), nNodesE_G+2, INFO )

    END DO

  END SUBROUTINE SolveLinearSystems


  SUBROUTINE UpdateNumberFlux( dt )

    REAL(DP), INTENT(in) :: dt

    INTEGER  :: iX, iE
    REAL(DP) :: Sigma_NES
    REAL(DP) :: DampFactor

    DO iX = 1, nNodesX_G
      DO iE = 1, nNodesE_G

        Sigma_NES &
          = SUM( W2_N(:) * R0_In(:,iE,iX) * uPR_N(:,iPR_D,iX) ) &
            + SUM( W2_N(:) * R0_Out(:,iE,iX) * ( FourPi - uPR_N(:,iPR_D,iX) ) )

        DampFactor &
          = 1.0 / ( 1.0_DP + dt * ( Chi(iE,iX) + FourPi * Sigma(iE,iX) &
                                    + Sigma_NES ) )

        uPR_N(iE,iPR_I1,iX) &
          = DampFactor * uPR_N(iE,iPR_I1,iX)

        uPR_N(iE,iPR_I2,iX) &
          = DampFactor * uPR_N(iE,iPR_I2,iX)

        uPR_N(iE,iPR_I3,iX) &
          = DampFactor * uPR_N(iE,iPR_I3,iX)

      END DO
    END DO

  END SUBROUTINE UpdateNumberFlux


  PURE FUNCTION RHSVEC( N, W, D, D_0, Chi, R_In, R_Out )

    INTEGER,                  INTENT(in) :: N
    REAL(DP), DIMENSION(N),   INTENT(in) :: W, D, D_0, Chi
    REAL(DP), DIMENSION(N,N), INTENT(in) :: R_In, R_Out
    REAL(DP), DIMENSION(N)               :: RHSVEC

    INTEGER :: i

    DO i = 1, N

      ! --- Emission/Absorption ---

      RHSVEC(i) &
        = Chi(i) * ( FourPi * D_0(i) - D(i) )

      ! --- Neutrino-Electron Scattering ---

      RHSVEC(i) &
        = RHSVEC(i) &
          ! --- In-Scattering Term ---
          + ( FourPi - D(i) ) * SUM( W(:) * R_In(:,i) * D(:) ) &
          ! --- Out-Scattering Term ---
          - D(i) * SUM( W(:) * R_Out(:,i) * ( FourPi - D(:) ) )

    END DO

    RETURN
  END FUNCTION RHSVEC


END MODULE FluidRadiationCouplingSolutionModule_B85
