MODULE InitializationModule_Relativistic

  USE KindModule, ONLY: &
    DP, &
    Zero, &
    Half, &
    One, &
    Two, &
    Three, &
    Four, &
    FourPi
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    nNodesX, &
    nDOFX, &
    iX_B0, &
    iX_B1, &
    iX_E0, &
    iX_E1
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    uGF, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_Alpha, &
    iGF_Psi
  USE FluidFieldsModule, ONLY: &
    uPF, &
    iPF_D, &
    iPF_V1, &
    iPF_V2, &
    iPF_V3, &
    iPF_E, &
    iPF_Ne, &
    uCF, &
    iCF_D, &
    iCF_S1, &
    iCF_S2, &
    iCF_S3, &
    iCF_E, &
    iCF_Ne, &
    uAF, &
    iAF_P
  USE Euler_BoundaryConditionsModule, ONLY: &
    ExpD, &
    ExpE
  USE EquationOfStateModule_IDEAL, ONLY: &
    Gamma_IDEAL, &
    ComputePressureFromPrimitive_IDEAL
  USE Euler_UtilitiesModule_Relativistic, ONLY: &
    ComputeConserved_Euler_Relativistic
  USE UnitsModule, ONLY: &
    GravitationalConstant, &
    SpeedOfLight, &
    Kilometer, &
    SolarMass, &
    Gram, &
    Centimeter, &
    Erg, &
    Second
  USE UtilitiesModule, ONLY: &
    NodeNumberX

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFields_Relativistic


CONTAINS


  SUBROUTINE InitializeFields_Relativistic &
               ( MassPNS_Option, &
                 ShockRadius_Option, &
                 AccretionRate_Option, &
                 PolytropicConstant_Option, &
                 ApplyPerturbation_Option, &
                 PerturbationOrder_Option, &
                 PerturbationAmplitude_Option, &
                 rPerturbationInner_Option, &
                 rPerturbationOuter_Option, &
                 InitializeFluidFromFile_Option, &
                 InitialConditionsFileName_Option )

    REAL(DP),           INTENT(in), OPTIONAL :: MassPNS_Option
    REAL(DP),           INTENT(in), OPTIONAL :: ShockRadius_Option
    REAL(DP),           INTENT(in), OPTIONAL :: AccretionRate_Option
    REAL(DP),           INTENT(in), OPTIONAL :: PolytropicConstant_Option
    LOGICAL,            INTENT(in), OPTIONAL :: ApplyPerturbation_Option
    INTEGER,            INTENT(in), OPTIONAL :: PerturbationOrder_Option
    REAL(DP),           INTENT(in), OPTIONAL :: PerturbationAmplitude_Option
    REAL(DP),           INTENT(in), OPTIONAL :: rPerturbationInner_Option
    REAL(DP),           INTENT(in), OPTIONAL :: rPerturbationOuter_Option
    LOGICAL,            INTENT(in), OPTIONAL :: InitializeFluidFromFile_Option
    CHARACTER(LEN=128), INTENT(in), OPTIONAL :: InitialConditionsFileName_Option

    ! --- Standing Accretion Shock (Defaults) ---

    REAL(DP)           :: &
      MassPNS = 1.4_DP * SolarMass
    REAL(DP)           :: &
      ShockRadius = 180.0_DP * Kilometer
    REAL(DP)           :: &
      AccretionRate = 0.3_DP * ( SolarMass / Second )
    REAL(DP),PARAMETER :: &
      PolytropicConstant2 = 2.0e14_DP &
                              * ( ( Erg / Centimeter**3 ) &
                              / ( Gram / Centimeter**3 )**( Four / Three ) )
    REAL(DP)       :: PolytropicConstant        = PolytropicConstant2
    LOGICAL        :: ApplyPerturbation         = .FALSE.
    INTEGER        :: PerturbationOrder         = 0
    REAL(DP)       :: PerturbationAmplitude     = Zero
    REAL(DP)       :: rPerturbationInner        = Zero
    REAL(DP)       :: rPerturbationOuter        = Zero
    LOGICAL        :: InitializeFluidFromFile   = .FALSE.
    CHARACTER(128) :: InitialConditionsFileName = 'InitialConditions.dat'

    uPF(:,:,:,:,iPF_Ne) = Zero

    IF( PRESENT( MassPNS_Option ) ) &
      MassPNS = MassPNS_Option
    IF( PRESENT( ShockRadius_Option ) ) &
      ShockRadius = ShockRadius_Option
    IF( PRESENT( AccretionRate_Option ) ) &
      AccretionRate = AccretionRate_Option
    IF( PRESENT( PolytropicConstant_Option ) ) &
      PolytropicConstant = PolytropicConstant_Option
    IF( PRESENT( ApplyPerturbation_Option ) ) &
      ApplyPerturbation = ApplyPerturbation_Option
    IF( PRESENT( PerturbationOrder_Option ) ) &
      PerturbationOrder = PerturbationOrder_Option
    IF( PRESENT( PerturbationAmplitude_Option ) ) &
      PerturbationAmplitude = PerturbationAmplitude_Option
    IF( PRESENT( rPerturbationInner_Option ) ) &
      rPerturbationInner = rPerturbationInner_Option
    IF( PRESENT( rPerturbationOuter_Option ) ) &
      rPerturbationOuter = rPerturbationOuter_Option
    IF( PRESENT( InitializeFluidFromFile_Option ) ) &
      InitializeFluidFromFile = InitializeFluidFromFile_Option
    IF( PRESENT( InitialConditionsFileName_Option ) ) &
      InitialConditionsFileName = InitialConditionsFileName_Option

    WRITE(*,*)
    WRITE(*,'(A,A)') '    INFO: ', TRIM( ProgramName )

    SELECT CASE ( TRIM( ProgramName ) )

      CASE( 'StandingAccretionShock' )

        CALL InitializeFields_StandingAccretionShock &
               ( MassPNS, &
                 ShockRadius, &
                 AccretionRate, &
                 PolytropicConstant, &
                 ApplyPerturbation, &
                 PerturbationOrder, &
                 PerturbationAmplitude, &
                 rPerturbationInner, &
                 rPerturbationOuter, &
                 InitializeFluidFromFile, &
                 InitialConditionsFileName )

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,'(A21,A)') 'Invalid ProgramName: ', ProgramName
        WRITE(*,'(A)')     'Stopping...'
        STOP

    END SELECT

#if defined(THORNADO_OMP_OL)
  !$OMP TARGET UPDATE TO( uCF )
#elif defined(THORNADO_OACC)
  !$ACC UPDATE DEVICE   ( uCF )
#endif

  END SUBROUTINE InitializeFields_Relativistic


  SUBROUTINE InitializeFields_StandingAccretionShock &
    ( MassPNS, &
      ShockRadius, &
      AccretionRate, &
      PolytropicConstant, &
      ApplyPerturbation, &
      PerturbationOrder, &
      PerturbationAmplitude, &
      rPerturbationInner, &
      rPerturbationOuter, &
      InitializeFluidFromFile, &
      InitialConditionsFileName )

    REAL(DP),           INTENT(in) :: MassPNS
    REAL(DP),           INTENT(in) :: ShockRadius
    REAL(DP),           INTENT(in) :: AccretionRate
    REAL(DP),           INTENT(in) :: PolytropicConstant
    LOGICAL,            INTENT(in) :: ApplyPerturbation
    INTEGER,            INTENT(in) :: PerturbationOrder
    REAL(DP),           INTENT(in) :: PerturbationAmplitude
    REAL(DP),           INTENT(in) :: rPerturbationInner
    REAL(DP),           INTENT(in) :: rPerturbationOuter
    LOGICAL,            INTENT(in) :: InitializeFluidFromFile
    CHARACTER(LEN=128), INTENT(in) :: InitialConditionsFileName

    INTEGER  :: iX1, iX2, iX3, iNodeX1, iNodeX2, iNodeX3, iNodeX
    INTEGER  :: iX1_1, iX1_2, iNodeX1_1, iNodeX1_2
    REAL(DP) :: X1_1, X1_2, D_1, D_2, V_1, V_2, P_1, P_2, C1, C2, C3
    REAL(DP) :: D0, V0, P0
    REAL(DP) :: X1, X2, dX1, Ka, Kb, Mdot, W
    REAL(DP) :: D    (1:nNodesX(1),iX_B1(1):iX_E1(1))
    REAL(DP) :: V    (1:nNodesX(1),iX_B1(1):iX_E1(1))
    REAL(DP) :: P    (1:nNodesX(1),iX_B1(1):iX_E1(1))
    REAL(DP) :: Alpha(1:nNodesX(1),iX_B1(1):iX_E1(1))
    REAL(DP) :: Psi  (1:nNodesX(1),iX_B1(1):iX_E1(1))
    LOGICAL  :: FirstPreShockElement = .FALSE.

    INTEGER, PARAMETER :: nX_LeastSquares = 5
    REAL(DP)           :: lnR(nX_LeastSquares,nNodesX(1)), &
                          lnD(nX_LeastSquares,nNodesX(1)), &
                          lnE(nX_LeastSquares,nNodesX(1))

    ! --- Quantities with (1) are pre-shock, those with (2) are post-shock ---

    Mdot = AccretionRate
    Ka   = PolytropicConstant

    WRITE(*,*)
    WRITE(*,'(6x,A,L)') &
      'Initialize From File:            ', &
      InitializeFluidFromFile
    WRITE(*,'(6x,A,ES9.2E3,A)') &
      'Shock radius:                    ', &
      ShockRadius / Kilometer, ' km'
    WRITE(*,'(6x,A,ES9.2E3,A)') &
      'PNS Mass:                        ', &
      MassPNS / SolarMass, ' Msun'
    WRITE(*,'(6x,A,ES9.2E3,A)') &
      'Accretion Rate:                  ', &
      Mdot / ( SolarMass / Second ), ' Msun/s'
    WRITE(*,'(6x,A,ES9.2E3,A)') &
      'Polytropic Constant (pre-shock): ', &
        Ka / ( ( Erg / Centimeter**3 ) &
               / ( Gram / Centimeter**3 )**( Gamma_IDEAL ) ), &
        ' erg/cm^3 / (g/cm^3)^( Gamma )'
    WRITE(*,*)
    WRITE(*,'(6x,A,L)') &
      'Apply Perturbation:           ', ApplyPerturbation
    WRITE(*,'(6x,A,I1)') &
      'Perturbation order:           ', PerturbationOrder
    WRITE(*,'(6x,A,ES9.2E3)') &
      'Perturbation amplitude:       ', PerturbationAmplitude
    WRITE(*,'(6x,A,ES9.2E3,A)') &
      'Inner radius of perturbation: ', rPerturbationInner / Kilometer, ' km'
    WRITE(*,'(6x,A,ES9.2E3,A)') &
      'Outer radius of perturbation: ', rPerturbationOuter / Kilometer, ' km'

    ! --- Make local copies of Lapse and Conformal Factor ---

    DO iX3 = iX_B1(3), iX_E1(3)
    DO iX2 = iX_B1(2), iX_E1(2)
    DO iX1 = iX_B1(1), iX_E1(1)

      DO iNodeX3 = 1, nNodesX(3)
      DO iNodeX2 = 1, nNodesX(2)
      DO iNodeX1 = 1, nNodesX(1)

        iNodeX = NodeNumberX(iNodeX1,iNodeX2,iNodeX3)

        Alpha(iNodeX1,iX1) = uGF(iNodeX,iX1,iX2,iX3,iGF_Alpha)
        Psi  (iNodeX1,iX1) = uGF(iNodeX,iX1,iX2,iX3,iGF_Psi)

      END DO
      END DO
      END DO

    END DO
    END DO
    END DO

    IF( InitializeFluidFromFile )THEN

      CALL ReadFluidFieldsFromFile &
             ( iX_B1, iX_E1, D, V, P, InitialConditionsFileName )

    ELSE

      ! --- Locate first element/node containing un-shocked fluid ---

      X1 = Zero

      DO iX1 = iX_B1(1), iX_E1(1)

        DO iNodeX1 = 1, nNodesX(1)

          dX1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 ) - X1
          X1  = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

          IF( X1 .LE. ShockRadius ) CYCLE

          IF( X1 .GT. ShockRadius .AND. .NOT. FirstPreShockElement )THEN

            iX1_1     = iX1
            iNodeX1_1 = iNodeX1
            X1_1      = X1
            X1_2      = X1 - dX1

            IF( iNodeX1_1 .EQ. 1 )THEN

              iX1_2     = iX1_1 - 1
              iNodeX1_2 = nNodesX(1)

            ELSE

              iX1_2     = iX1_1
              iNodeX1_2 = iNodeX1_1 - 1

            END IF

            FirstPreShockElement = .TRUE.

          END IF

        END DO

      END DO

      ! --- Pre-shock Fields ---

      X1 = NodeCoordinate( MeshX(1), iX_E1(1), nNodesX(1) )

      ! --- Use Newtonian values as initial guesses ---

      V0 = -SQRT( Two * GravitationalConstant * MassPNS / X1 )
      D0 = -Mdot / ( FourPi * X1**2 * V0 )
      P0 = Ka * D0**( Gamma_IDEAL )

      DO iX1 = iX_E1(1), iX1_1, -1

        DO iNodeX1 = nNodesX(1), 1, -1

          X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

          IF( X1 .LE. ShockRadius ) CYCLE

          CALL NewtonRaphson_SAS &
                 ( X1, MassPNS, Ka, Mdot, &
                   Alpha(iNodeX1,iX1), Psi(iNodeX1,iX1), D0, V0, P0, &
                   D(iNodeX1,iX1), V(iNodeX1,iX1), P(iNodeX1,iX1) )

          D0 = D(iNodeX1,iX1)
          V0 = V(iNodeX1,iX1)
          P0 = P(iNodeX1,iX1)

        END DO

      END DO

      ! --- Apply Jump Conditions ---

      D_1 = D(iNodeX1_1,iX1_1)
      V_1 = V(iNodeX1_1,iX1_1)
      P_1 = P(iNodeX1_1,iX1_1)

      W = LorentzFactor( Psi(iNodeX1_1,iX1), V_1 )

      C1 = D_1 * W * V_1
      C2 = D_1 &
             * ( SpeedOfLight**2 + Gamma_IDEAL / ( Gamma_IDEAL - One ) &
                   * P_1 / D_1  ) * W**2 * V_1**2 / SpeedOfLight**2 &
             + Psi(iNodeX1_1,iX1)**( -4 ) * P_1
      C3 = D_1 &
             * ( SpeedOfLight**2 + Gamma_IDEAL / ( Gamma_IDEAL - One ) &
                   * P_1 / D_1  ) * W**2 * V_1

      CALL ApplyJumpConditions_SAS &
             ( Psi(iNodeX1_1,iX1), V_1, C1, C2, C3, D_2, V_2, P_2 )

      Kb = P_2 / D_2**( Gamma_IDEAL )

      WRITE(*,*)
      WRITE(*,'(6x,A)') 'Jump Conditions'
      WRITE(*,'(6x,A)') '---------------'
      WRITE(*,*)
      WRITE(*,'(8x,A)') 'Pre-shock:'
      WRITE(*,'(10x,A,I4.4)')       'iX1      = ', iX1_1
      WRITE(*,'(10x,A,I2.2)')       'iNodeX1  = ', iNodeX1_1
      WRITE(*,'(10x,A,ES13.6E3,A)') 'X1       = ', X1_1 / Kilometer, '  km'
      WRITE(*,'(10x,A,ES13.6E3,A)') 'Density  = ', &
        D_1 / ( Gram / Centimeter**3 ), '  g/cm^3'
      WRITE(*,'(10x,A,ES14.6E3,A)') 'Velocity = ', &
        V_1 / ( Kilometer / Second ), ' km/s'
      WRITE(*,'(10x,A,ES13.6E3,A)') 'Pressure = ', &
        P_1 / ( Erg / Centimeter**3 ), '  erg/cm^3'
      WRITE(*,*)
      WRITE(*,'(8x,A)') 'Post-shock:'
      WRITE(*,'(10x,A,I4.4)')       'iX1      = ', iX1_2
      WRITE(*,'(10x,A,I2.2)')       'iNodeX1  = ', iNodeX1_2
      WRITE(*,'(10x,A,ES13.6E3,A)') 'X1       = ', X1_2 / Kilometer, '  km'
      WRITE(*,'(10x,A,ES13.6E3,A)') 'Density  = ', &
        D_2 / ( Gram / Centimeter**3 ), '  g/cm^3'
      WRITE(*,'(10x,A,ES14.6E3,A)') 'Velocity = ', &
        V_2 / ( Kilometer / Second ), ' km/s'
      WRITE(*,'(10x,A,ES13.6E3,A)') 'Pressure = ', &
        P_2 / ( Erg / Centimeter**3 ), '  erg/cm^3'
      WRITE(*,*)

      ! --- Post-shock Fields ---

      D0 = D_2
      V0 = V_2
      P0 = P_2

      DO iX1 = iX1_2, iX_B1(1), -1

        DO iNodeX1 = nNodesX(1), 1, -1

          X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

          IF( X1 .GT. ShockRadius ) CYCLE

          CALL NewtonRaphson_SAS &
                 ( X1, MassPNS, Kb, Mdot, &
                   Alpha(iNodeX1,iX1), Psi(iNodeX1,iX1), D0, V0, P0, &
                   D(iNodeX1,iX1), V(iNodeX1,iX1), P(iNodeX1,iX1) )

          D0 = D(iNodeX1,iX1)
          V0 = V(iNodeX1,iX1)
          P0 = P(iNodeX1,iX1)

        END DO

      END DO

    END IF

    ! --- Map to 3D Domain ---

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B1(1), iX_E1(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)
        iNodeX2 = NodeNumberTableX(2,iNodeX)

        IF( ApplyPerturbation )THEN

          X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
          X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

          IF( X1 .GE. rPerturbationInner .AND. X1 .LE. rPerturbationOuter )THEN

            IF( PerturbationOrder .EQ. 0 ) &
              uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
                = D(iNodeX1,iX1) * ( One + PerturbationAmplitude )

            IF( PerturbationOrder .EQ. 1 ) &
              uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
                = D(iNodeX1,iX1) * ( One + PerturbationAmplitude * COS( X2 ) )

          ELSE

            uPF(iNodeX,iX1,iX2,iX3,iPF_D) = D(iNodeX1,iX1)

          END IF

        ELSE

          uPF(iNodeX,iX1,iX2,iX3,iPF_D) = D(iNodeX1,iX1)

        END IF

        uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = V(iNodeX1,iX1)
        uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = Zero
        uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = Zero
        uPF(iNodeX,iX1,iX2,iX3,iPF_E ) = P(iNodeX1,iX1) / ( Gamma_IDEAL - One )

      END DO

      CALL ComputePressureFromPrimitive_IDEAL &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_E), &
               uPF(:,iX1,iX2,iX3,iPF_Ne), uAF(:,iX1,iX2,iX3,iAF_P) )

      CALL ComputeConserved_Euler_Relativistic &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_V1), &
               uPF(:,iX1,iX2,iX3,iPF_V2), uPF(:,iX1,iX2,iX3,iPF_V3), &
               uPF(:,iX1,iX2,iX3,iPF_E ), uPF(:,iX1,iX2,iX3,iPF_Ne), &
               uCF(:,iX1,iX2,iX3,iCF_D ), uCF(:,iX1,iX2,iX3,iCF_S1), &
               uCF(:,iX1,iX2,iX3,iCF_S2), uCF(:,iX1,iX2,iX3,iCF_S3), &
               uCF(:,iX1,iX2,iX3,iCF_E ), uCF(:,iX1,iX2,iX3,iCF_Ne), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33), &
               uAF(:,iX1,iX2,iX3,iAF_P) )

      IF( iX1 .LT. nX_LeastSquares+1 .AND. iX1 .GT. 0 )THEN

        DO iNodeX1 = 1, nNodesX(1)

          lnD(iX1,iNodeX1) = LOG( uCF(iNodeX1,iX1,iX2,iX3,iCF_D) )
          lnE(iX1,iNodeX1) = LOG( uCF(iNodeX1,iX1,iX2,iX3,iCF_E) )

        END DO

      END IF

    END DO
    END DO
    END DO

    DO iX1 = 1, nX_LeastSquares

      DO iNodeX1 = 1, nNodesX(1)

        lnR(iX1,iNodeX1) = LOG( NodeCoordinate( MeshX(1), iX1, iNodeX1 ) )

      END DO

    END DO

    ! --- Expression for exponents from:
    !     https://mathworld.wolfram.com/LeastSquaresFittingPowerLaw.html ---

    ExpD = -( nX_LeastSquares * nNodesX(1) * SUM( lnR * lnD ) &
               - SUM( lnR ) * SUM( lnD ) ) &
             / ( nX_LeastSquares * nNodesX(1) * SUM( lnR**2 ) &
               - SUM( lnR )**2 )
    ExpE = -( nX_LeastSquares * nNodesX(1) * SUM( lnR * lnE ) &
               - SUM( lnR ) * SUM( lnE ) ) &
             / ( nX_LeastSquares * nNodesX(1) * SUM( lnR**2 ) &
               - SUM( lnR )**2 )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE TO( ExpD, ExpE )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE DEVICE   ( ExpD, ExpE )
#endif

    WRITE(*,'(6x,A,I2.2)') &
      'nX_LeastSquares: ', &
      nX_LeastSquares
    WRITE(*,'(6x,A,F8.6)') &
      'ExpD:            ', &
      ExpD
    WRITE(*,'(6x,A,F8.6)') &
      'ExpE:            ', &
      ExpE

  END SUBROUTINE InitializeFields_StandingAccretionShock


  ! --- Auxiliary functions for standing accretion shock problem ---


  SUBROUTINE NewtonRaphson_SAS &
    ( X1, MassPNS, K, Mdot, Alpha, Psi, D0, V0, P0, D, V, P )

    REAL(DP), INTENT(in)  :: X1, MassPNS, K, &
                             Mdot, Alpha, Psi, D0, V0, P0
    REAL(DP), INTENT(out) :: D ,V ,P

    REAL(DP) :: W
    REAL(DP) :: Jac(3,3), invJac(3,3)
    REAL(DP) :: f(3), uO(3), uN(3), du(3)

    LOGICAL             :: CONVERGED
    INTEGER             :: ITER
    REAL(DP), PARAMETER :: Tolu = 1.0e-16_DP
    REAL(DP), PARAMETER :: Tolf = 1.0e-16_DP
    INTEGER,  PARAMETER :: MAX_ITER = 4 - INT( LOG( Tolu ) /  LOG( Two ) )

    uO(1) = One
    uO(2) = One
    uO(3) = One

    CONVERGED = .FALSE.
    ITER      = 0
    DO WHILE( .NOT. CONVERGED .AND. ITER .LT. MAX_ITER )

      ITER = ITER + 1

      W = LorentzFactor( Psi, V0 * uO(2) )

      f(1) &
        = FourPi * Alpha * Psi**6 * X1**2 / Mdot * D0 * V0 &
            * uO(1) * W * uO(2) + One
      f(2) &
        = Alpha * ( One + Gamma_IDEAL / ( Gamma_IDEAL - One ) &
                      * P0 / ( D0 * SpeedOfLight**2 ) * uO(3) / uO(1) ) * W &
            - One
      f(3) &
        = P0 * D0**( -Gamma_IDEAL ) / K * uO(3) * uO(1)**( -Gamma_IDEAL ) - One

      Jac(1,1) = FourPi * Alpha * Psi**6 * X1**2 / Mdot * D0 * V0 &
                   * W * uO(2)
      Jac(1,2) = FourPi * Alpha * Psi**6 * X1**2 / Mdot * D0 * V0 &
                   * uO(1) * W**3
      Jac(1,3) = Zero
      Jac(2,1) = -Alpha * Gamma_IDEAL / ( Gamma_IDEAL - One ) &
                   * P0 / ( D0 * SpeedOfLight**2 ) * uO(3) / uO(1)**2 * W
      Jac(2,2) = Alpha * ( One + Gamma_IDEAL / ( Gamma_IDEAL - One ) &
                   * P0 / ( D0 * SpeedOfLight**2 ) * uO(3) / uO(1) ) &
                   * W**3 * Psi**4 * V0**2 / SpeedOfLight**2 * uO(2)
      Jac(2,3) = Alpha * Gamma_IDEAL / ( Gamma_IDEAL - One ) &
                   * P0 / ( D0 * SpeedOfLight**2 ) / uO(1) * W
      Jac(3,1) = -P0 * D0**( -Gamma_IDEAL ) / K &
                   * Gamma_IDEAL * uO(3) * uO(1)**( -Gamma_IDEAL - One )

      Jac(3,2) = Zero
      Jac(3,3) = P0 * D0**( -Gamma_IDEAL ) / K &
                   * uO(1)**( -Gamma_IDEAL )

      InvJac = Inv3x3( Jac )

      uN = uO - MATMUL( InvJac, f )

      du = uN - uO

      IF( MAXVAL( ABS( du / uO ) ) .LT. Tolu ) CONVERGED = .TRUE.

      uO = uN

    END DO

    D = uN(1) * D0
    V = uN(2) * V0
    P = uN(3) * P0

  END SUBROUTINE NewtonRaphson_SAS


  SUBROUTINE ApplyJumpConditions_SAS( Psi, V_1, C1, C2, C3, D_2, V_2, P_2 )

    REAL(DP), INTENT(in)  :: Psi, V_1, C1, C2, C3
    REAL(DP), INTENT(out) :: D_2, V_2, P_2

    REAL(DP) :: A, B, C, D, E
    REAL(DP) :: dx, xa, xb, xc, fa, fb, fc, W

    INTEGER             :: ITER
    INTEGER,  PARAMETER :: MAX_ITER = 1000
    REAL(DP), PARAMETER :: TolChi = 1.0e-16_DP

    LOGICAL :: CONVERGED

    A = SpeedOfLight**( -4 ) * ( C3 / C1 )**2 - One
    B = -Two * SpeedOfLight**( 3 ) * C2 * C3 / C1**2 &
          * Gamma_IDEAL / ( Gamma_IDEAL - One ) * Psi**4
    C = SpeedOfLight**( -2 ) * ( C2 / C1 )**2 &
          * ( Gamma_IDEAL / ( Gamma_IDEAL - One ) )**2 * Psi**8 &
          + Two * SpeedOfLight**( -4 ) * ( C3 / C1 )**2 &
          / ( Gamma_IDEAL - One ) * Psi**4 + Psi**4
    D = -Two * SpeedOfLight**( -3 ) * C2 * C3 / C1**2 &
          * Gamma_IDEAL / ( Gamma_IDEAL - One )**2 * Psi**8
    E = SpeedOfLight**( -4 ) * ( C3 / C1 )**2 &
          / ( Gamma_IDEAL - One )**2 * Psi**8

    ! --- Solve with bisection ---

    ! Add 1 km/s to exclude smooth solution
    xa = ( V_1 + One * Kilometer / Second ) / SpeedOfLight

    xb = 1.0e-10_DP * xa

    fa = A * xa**( -2 ) + B * xa**( -1 ) &
             + C + D * xa + E * xa**2
    fb = A * xb**( -2 ) + B * xb**( -1 ) &
             + C + D * xb + E * xb**2

    dx = xb - xa

    ITER      = 0
    CONVERGED = .FALSE.
    DO WHILE( .NOT. CONVERGED .AND. ITER .LT. MAX_ITER )

      ITER = ITER + 1

      dx = Half * dx

      xc = xa + dx

      fc = A * xc**( -2 ) + B * xc**( -1 ) + C + D * xc + E * xc**2

      IF( fa * fc .LT. Zero )THEN

        xb = xc
        fb = fc

      ELSE IF( fa * fc .GT. Zero )THEN

        xa = xc
        fa = fc

      ELSE

        CONVERGED = .TRUE.

      END IF

      IF( ABS( dx / xa ) .LT. TolChi ) CONVERGED = .TRUE.

    END DO

    V_2 = xc * SpeedOfLight

    W = LorentzFactor( Psi, V_2 )

    D_2 = SpeedOfLight**( -1 ) * ABS( C1 ) * SQRT( xc**( -2 ) - Psi**4 )
    P_2 = ( C3 - D_2 * SpeedOfLight**2 * W**2 * V_2 ) &
            / ( Gamma_IDEAL / ( Gamma_IDEAL - One ) * W**2 * V_2 )

  END SUBROUTINE ApplyJumpConditions_SAS


  ! --- From: http://fortranwiki.org/fortran/show/Matrix+inversion ---
  FUNCTION Inv3x3( A ) RESULT( invA )

    ! --- Performs a direct calculation of the inverse of a 3×3 matrix ---

    REAL(DP), INTENT(in) :: A   (3,3)
    REAL(DP)             :: invA(3,3)
    REAL(DP)             :: InvDet

    ! --- Calculate the inverse of the determinant of the matrix ---

    InvDet = One / ( A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)     &
                       - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1) &
                       + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1) )

    ! --- Calculate the inverse of the matrix ---

    invA(1,1) = +InvDet * ( A(2,2)*A(3,3) - A(2,3)*A(3,2) )
    invA(2,1) = -InvDet * ( A(2,1)*A(3,3) - A(2,3)*A(3,1) )
    invA(3,1) = +InvDet * ( A(2,1)*A(3,2) - A(2,2)*A(3,1) )
    invA(1,2) = -InvDet * ( A(1,2)*A(3,3) - A(1,3)*A(3,2) )
    invA(2,2) = +InvDet * ( A(1,1)*A(3,3) - A(1,3)*A(3,1) )
    invA(3,2) = -InvDet * ( A(1,1)*A(3,2) - A(1,2)*A(3,1) )
    invA(1,3) = +InvDet * ( A(1,2)*A(2,3) - A(1,3)*A(2,2) )
    invA(2,3) = -InvDet * ( A(1,1)*A(2,3) - A(1,3)*A(2,1) )
    invA(3,3) = +InvDet * ( A(1,1)*A(2,2) - A(1,2)*A(2,1) )

    RETURN
  END FUNCTION Inv3x3


  REAL(DP) FUNCTION LorentzFactor( Psi, V )

    REAL(DP), INTENT(in) :: Psi, V

    LorentzFactor = One / SQRT( One - Psi**4 * ( V / SpeedOfLight )**2 )

    RETURN
  END FUNCTION LorentzFactor


  SUBROUTINE ReadFluidFieldsFromFile &
    ( iX_B1, iX_E1, D, V, P, InitialConditionsFileName )

    INTEGER,  INTENT(in)           :: &
      iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(out)          :: &
      D(1:,iX_B1(1):), V(1:,iX_B1(1):), P(1:,iX_B1(1):)
    CHARACTER(LEN=128), INTENT(in) :: &
      InitialConditionsFileName

    CHARACTER(LEN=16) :: FMT
    INTEGER           :: iX1

    D = Zero
    V = Zero
    P = Zero

    WRITE(*,*)
    WRITE(*,'(6x,A,A)') &
      'Initial Conditions File: ', TRIM( InitialConditionsFileName )
    WRITE(*,*)

    OPEN( UNIT = 101, FILE = TRIM( InitialConditionsFileName ) // '_D.dat' )
    OPEN( UNIT = 102, FILE = TRIM( InitialConditionsFileName ) // '_V.dat' )
    OPEN( UNIT = 103, FILE = TRIM( InitialConditionsFileName ) // '_P.dat' )

    READ(101,*) FMT
    READ(102,*) FMT
    READ(103,*) FMT

    DO iX1 = iX_B1(1), iX_E1(1)

      READ(101,TRIM(FMT)) D(:,iX1)
      READ(102,TRIM(FMT)) V(:,iX1)
      READ(103,TRIM(FMT)) P(:,iX1)

    END DO

    CLOSE( 103 )
    CLOSE( 102 )
    CLOSE( 101 )

  END SUBROUTINE ReadFluidFieldsFromFile


END MODULE InitializationModule_Relativistic
