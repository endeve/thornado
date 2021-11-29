MODULE InitializationModule

  USE KindModule, ONLY: &
    DP, Zero, Half, One, Pi, TwoPi
  USE UnitsModule, ONLY: &
    Gram, Centimeter, &
    Kilometer, Erg, Second, Kelvin, &
    SpeedOfLight
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    iX_B0, iX_E0, iX_B1, iX_E1, &
    nNodesX, nDOFX
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    uGF, iGF_Phi_N, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33
  USE FluidFieldsModule, ONLY: &
    uPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne, &
    uCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    uAF, iAF_P, iAF_Ye, iAF_T, iAF_E, iAF_S, iAF_Me, &
    iAF_Mp, iAF_Mn, iAF_Xp, iAF_Xn, iAF_Xa, iAF_Xh, iAF_Gm
  USE Euler_UtilitiesModule_NonRelativistic, ONLY: &
    ComputeConserved_Euler_NonRelativistic
  USE EquationOfStateModule, ONLY: &
    ComputeTemperatureFromPressure, &
    ComputeThermodynamicStates_Primitive, &
    ApplyEquationOfState
  USE ProgenitorModule, ONLY: &
    ProgenitorType1D, &
    ReadProgenitor1D
  USE UtilitiesModule, ONLY: &
    Locate, &
    NodeNumberX, &
    Interpolate1D_Linear

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFields

CONTAINS


  SUBROUTINE InitializeFields &
    ( AdvectionProfile_Option, RiemannProblemName_Option, &
      ProgenitorFileName_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: &
      AdvectionProfile_Option, &
      RiemannProblemName_Option, &
      ProgenitorFileName_Option

    WRITE(*,*)
    WRITE(*,'(A2,A6,A)') '', 'INFO: ', TRIM( ProgramName )

    SELECT CASE ( TRIM( ProgramName ) )

      CASE ( 'Advection' )

        CALL InitializeFields_Advection &
               ( AdvectionProfile_Option &
                   = AdvectionProfile_Option )

      CASE ( 'RiemannProblem' )

        CALL InitializeFields_RiemannProblem &
               ( RiemannProblemName_Option &
                   = RiemannProblemName_Option )

      CASE ( 'RiemannProblemSpherical' )

        CALL InitializeFields_RiemannProblemSpherical &
               ( RiemannProblemName_Option &
                   = RiemannProblemName_Option )

      CASE( 'RiemannProblemCylindrical' )

        CALL InitializeFields_RiemannProblemCylindrical &
               ( RiemannProblemName_Option &
                   = RiemannProblemName_Option )

      CASE( 'Implosion' )

        CALL InitializeFields_Implosion

      CASE( 'Jet' )

        CALL InitializeFields_Jet

      CASE( 'GravitationalCollapse' )

        CALL InitializeFields_GravitationalCollapse &
               ( ProgenitorFileName_Option &
                   = ProgenitorFileName_Option )

      CASE( 'ShockEntropyWave' )

        CALL InitiaizeFields_InitializeShockEntropyWaveInteraction1D

    END SELECT

  END SUBROUTINE InitializeFields


  SUBROUTINE InitializeFields_Advection &
    ( AdvectionProfile_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: &
      AdvectionProfile_Option

    REAL(DP), PARAMETER :: D_0    = 1.0d12 * Gram / Centimeter**3
    REAL(DP), PARAMETER :: Amp    = 1.0d11 * Gram / Centimeter**3
    REAL(DP), PARAMETER :: L      = 1.0d02 * Kilometer
    REAL(DP), PARAMETER :: Ye_0   = 0.0125_DP
    REAL(DP), PARAMETER :: Amp_Ye = 0.030_DP

    CHARACTER(32) :: AdvectionProfile
    INTEGER       :: iX1, iX2, iX3
    INTEGER       :: iNodeX, iNodeX1
    REAL(DP)      :: X1

    IF( PRESENT( AdvectionProfile_Option ) )THEN
      AdvectionProfile = TRIM( AdvectionProfile_Option )
    ELSE
      AdvectionProfile = 'SineWave'
    END IF

    WRITE(*,*)
    WRITE(*,'(A4,A,A)') &
      '', 'Advection Profile: ', TRIM( AdvectionProfile )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

        SELECT CASE ( TRIM( AdvectionProfile ) )

          CASE( 'SineWave' )

            uPF(iNodeX,iX1,iX2,iX3,iPF_D ) &
              = D_0 + Amp * SIN( TwoPi * X1 / L )
            uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
              = 0.1_DP * SpeedOfLight
            uPF(iNodeX,iX1,iX2,iX3,iPF_V2) &
              = 0.0_DP * Kilometer / Second
            uPF(iNodeX,iX1,iX2,iX3,iPF_V3) &
              = 0.0_DP * Kilometer / Second
            uAF(iNodeX,iX1,iX2,iX3,iAF_P ) &
              = 1.0d-2 * D_0 * SpeedOfLight**2
            uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) &
              = 0.3_DP

          CASE( 'QuarticSineWave' )

            ! SIN^4 profile modified from Suresh and Huynh (1997),
            ! JCP 136, 83-99.

            uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
              = D_0 + Amp * SIN( Pi * X1 / L )**4
            uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
              = 0.1_DP * SpeedOfLight
            uPF(iNodeX,iX1,iX2,iX3,iPF_V2) &
              = 0.0_DP * Kilometer / Second
            uPF(iNodeX,iX1,iX2,iX3,iPF_V3) &
              = 0.0_DP * Kilometer / Second
            uAF(iNodeX,iX1,iX2,iX3,iAF_P) &
              = 1.0d-2 * D_0 * SpeedOfLight**2
            uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) &
              = 0.3_DP

          CASE( 'DiscontinuousMultiWave' )

            ! Discontinuous profile with multiple types of waves
            ! from Suresh and Huynh (1997), JCP 136, 83-99.

            IF( ( X1 .GE. -0.8_DP * L ) .AND. ( X1 .LE. -0.6_DP * L ) )THEN

              uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
                = D_0 + Amp &
                          * EXP( ( -LOG(2.0_DP) / ( 9.0d-4 ) ) &
                                 * ( ( X1 / L ) + 7.0d-1 )**2 )
              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
                = 0.1_DP * SpeedOfLight
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) &
                = 0.0_DP * Kilometer / Second
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) &
                = 0.0_DP * Kilometer / Second
              uAF(iNodeX,iX1,iX2,iX3,iAF_P) &
                = 1.0d-2 * D_0 * SpeedOfLight**2
              uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) &
                = 0.3_DP

            ELSE IF( ( X1 .GE. -0.4_DP * L ) .AND. ( X1 .LE. -0.2_DP * L ) )THEN

              uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
                = D_0 + Amp
              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
                = 0.1_DP * SpeedOfLight
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) &
                = 0.0_DP * Kilometer / Second
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) &
                = 0.0_DP * Kilometer / Second
              uAF(iNodeX,iX1,iX2,iX3,iAF_P) &
                = 1.0d-2 * D_0 * SpeedOfLight**2
              uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) &
                = 0.3_DP

            ELSE IF( ( X1 .GE. 0.0_DP * L ) .AND. ( X1 .LE. 0.2_DP * L ) )THEN

              uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
                = D_0 + Amp &
                          * ( 1.0_DP &
                                - ABS( 10.0_DP * ( ( X1 / L) - 0.1_DP ) ) )
              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
                = 0.1_DP * SpeedOfLight
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) &
                = 0.0_DP * Kilometer / Second
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) &
                = 0.0_DP * Kilometer / Second
              uAF(iNodeX,iX1,iX2,iX3,iAF_P) &
                = 1.0d-2 * D_0 * SpeedOfLight**2
              uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) &
                = 0.3_DP

            ELSE IF( ( X1 .GE. 0.4_DP * L ) .AND. ( X1 .LE. 0.6_DP * L ) )THEN

              uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
                = D_0 + Amp &
                        * SQRT( 1.0_DP - 1.0d2 * ( ( X1 / L ) - 0.5_DP )**2 )
              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
                = 0.1_DP * SpeedOfLight
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) &
                = 0.0_DP * Kilometer / Second
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) &
                = 0.0_DP * Kilometer / Second
              uAF(iNodeX,iX1,iX2,iX3,iAF_P) &
                = 1.0d-2 * D_0 * SpeedOfLight**2
              uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) &
                = 0.3_DP

            ELSE

              uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
                = D_0
              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
                = 0.1_DP * SpeedOfLight
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) &
                = 0.0_DP * Kilometer / Second
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) &
                = 0.0_DP * Kilometer / Second
              uAF(iNodeX,iX1,iX2,iX3,iAF_P) &
                = 1.0d-2 * D_0 * SpeedOfLight**2
              uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) &
                = 0.3_DP

            END IF

          CASE( 'TopHat' )

            IF( ( X1 .GE. -0.4_DP * L ) .AND. ( X1 .LE. 0.4_DP * L ) )THEN

              uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
                = D_0 + Amp
              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
                = 0.1_DP * SpeedOfLight
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) &
                = 0.0_DP * Kilometer / Second
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) &
                = 0.0_DP * Kilometer / Second
              uAF(iNodeX,iX1,iX2,iX3,iAF_P) &
                = 1.5d-2 * D_0 * SpeedOfLight**2
              uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) &
                = Ye_0 + Amp_Ye

            ELSE

              uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
                = D_0
              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
                = 0.1_DP * SpeedOfLight
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) &
                = 0.0_DP * Kilometer / Second
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) &
                = 0.0_DP * Kilometer / Second
              uAF(iNodeX,iX1,iX2,iX3,iAF_P) &
                = 1.5d-2 * D_0 * SpeedOfLight**2
              uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) &
                = Ye_0

            END IF

        END SELECT

        CALL ComputeTemperatureFromPressure &
               ( uPF(iNodeX,iX1,iX2,iX3,iPF_D ), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_P ), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_Ye), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_T ) )

        CALL ComputeThermodynamicStates_Primitive &
               ( uPF(iNodeX,iX1,iX2,iX3,iPF_D ), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_T ), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_Ye), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_E ), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_E ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_Ne) )

        CALL ComputeConserved_Euler_NonRelativistic &
               ( uPF(iNodeX,iX1,iX2,iX3,iPF_D ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V1), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V2), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V3), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_E ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_Ne), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_D ), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S1), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S2), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S3), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_E ), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_Ne), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33) )

      END DO

    END DO
    END DO
    END DO

  END SUBROUTINE InitializeFields_Advection


  SUBROUTINE InitializeFields_RiemannProblem &
    ( RiemannProblemName_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: &
      RiemannProblemName_Option

    CHARACTER(32) :: RiemannProblemName
    INTEGER       :: iX1, iX2, iX3
    INTEGER       :: iNodeX, iNodeX1
    REAL(DP)      :: X1

    RiemannProblemName = 'Sod'
    IF( PRESENT( RiemannProblemName_Option ) ) &
       RiemannProblemName = TRIM( RiemannProblemName_Option )

    WRITE(*,*)
    WRITE(*,'(A4,A,A)') &
      '', 'Riemann Problem Name: ', TRIM( RiemannProblemName )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

        SELECT CASE ( TRIM( RiemannProblemName ) )

          CASE( 'Sod' )

            IF( X1 <= Zero )THEN

              uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = 1.0d12 * Gram / Centimeter**3
              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP * Kilometer / Second
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP * Kilometer / Second
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP * Kilometer / Second
              uAF(iNodeX,iX1,iX2,iX3,iAF_P ) = 1.00d32 * Erg / Centimeter**3
              uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) = 0.4_DP

            ELSE

              uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = 1.25d11 * Gram / Centimeter**3
              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP * Kilometer / Second
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP * Kilometer / Second
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP * Kilometer / Second
              uAF(iNodeX,iX1,iX2,iX3,iAF_P ) = 1.00d31 * Erg / Centimeter**3
              uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) = 0.3_DP

            END IF

          CASE( 'Pochik' )

            IF( X1 <= Zero )THEN

              uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = 1.0d13 * Gram / Centimeter**3
              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP * Kilometer / Second
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP * Kilometer / Second
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP * Kilometer / Second
              uAF(iNodeX,iX1,iX2,iX3,iAF_P ) = 1.07d31 * Erg / Centimeter**3
              uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) = 0.04_DP

            ELSE

              uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = 1.25d12 * Gram / Centimeter**3
              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP * Kilometer / Second
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP * Kilometer / Second
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP * Kilometer / Second
              uAF(iNodeX,iX1,iX2,iX3,iAF_P ) = 1.023d30 * Erg / Centimeter**3
              uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) = 0.1_DP

            END IF

        END SELECT

        CALL ComputeTemperatureFromPressure &
               ( uPF(iNodeX,iX1,iX2,iX3,iPF_D ), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_P ), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_Ye), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_T ) )

        CALL ComputeThermodynamicStates_Primitive &
               ( uPF(iNodeX,iX1,iX2,iX3,iPF_D ), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_T ), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_Ye), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_E ), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_E ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_Ne) )

        CALL ComputeConserved_Euler_NonRelativistic &
               ( uPF(iNodeX,iX1,iX2,iX3,iPF_D ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V1), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V2), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V3), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_E ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_Ne), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_D ), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S1), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S2), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S3), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_E ), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_Ne), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33) )

      END DO

    END DO
    END DO
    END DO

  END SUBROUTINE InitializeFields_RiemannProblem


  SUBROUTINE InitializeFields_RiemannProblemSpherical &
    ( RiemannProblemName_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: &
      RiemannProblemName_Option

    CHARACTER(32) :: RiemannProblemName
    INTEGER       :: iX1, iX2, iX3
    INTEGER       :: iNodeX, iNodeX1
    REAL(DP)      :: X1

    RiemannProblemName = 'SphericalSod'
    IF( PRESENT( RiemannProblemName_Option ) ) &
       RiemannProblemName = TRIM( RiemannProblemName_Option )

    WRITE(*,*)
    WRITE(*,'(A4,A,A)') &
      '', 'Riemann Problem Name: ', TRIM( RiemannProblemName )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

        SELECT CASE ( TRIM( RiemannProblemName ) )

          CASE( 'SphericalSod' )

            IF( X1 <= 5.0_DP * Kilometer )THEN

              uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = 1.00d12 * Gram / Centimeter**3
              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP * Kilometer / Second
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP * Kilometer / Second
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP * Kilometer / Second
              uAF(iNodeX,iX1,iX2,iX3,iAF_P ) = 1.0d32 * Erg / Centimeter**3
              uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) = 0.4_DP

            ELSE

              uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = 1.25d11 * Gram / Centimeter**3
              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP * Kilometer / Second
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP * Kilometer / Second
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP * Kilometer / Second
              uAF(iNodeX,iX1,iX2,iX3,iAF_P ) = 1.0d31 * Erg / Centimeter**3
              uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) = 0.4_DP

            END IF

        END SELECT

      END DO

      CALL ComputeTemperatureFromPressure &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uAF(:,iX1,iX2,iX3,iAF_P), &
               uAF(:,iX1,iX2,iX3,iAF_Ye), uAF(:,iX1,iX2,iX3,iAF_T) )

      CALL ComputeThermodynamicStates_Primitive &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uAF(:,iX1,iX2,iX3,iAF_T ), &
               uAF(:,iX1,iX2,iX3,iAF_Ye), uPF(:,iX1,iX2,iX3,iPF_E ), &
               uAF(:,iX1,iX2,iX3,iAF_E ), uPF(:,iX1,iX2,iX3,iPF_Ne) )

      CALL ComputeConserved_Euler_NonRelativistic &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_V1), &
               uPF(:,iX1,iX2,iX3,iPF_V2), uPF(:,iX1,iX2,iX3,iPF_V3), &
               uPF(:,iX1,iX2,iX3,iPF_E ), uPF(:,iX1,iX2,iX3,iPF_Ne), &
               uCF(:,iX1,iX2,iX3,iCF_D ), uCF(:,iX1,iX2,iX3,iCF_S1), &
               uCF(:,iX1,iX2,iX3,iCF_S2), uCF(:,iX1,iX2,iX3,iCF_S3), &
               uCF(:,iX1,iX2,iX3,iCF_E ), uCF(:,iX1,iX2,iX3,iCF_Ne), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33) )

    END DO
    END DO
    END DO

  END SUBROUTINE InitializeFields_RiemannProblemSpherical


  SUBROUTINE InitiaizeFields_InitializeShockEntropyWaveInteraction1D

    INTEGER       :: iX1, iX2, iX3
    INTEGER       :: iNodeX, iNodeX1
    REAL(DP)      :: X1
    REAL(DP)      :: X_D, Amplitude, Wavenumber

    X_D        = - 4.0_DP * Kilometer
    Amplitude  = + 0.2d12
    Wavenumber = + 5.0_DP / Kilometer

    WRITE(*,*)
    WRITE(*,'(A2,A6,A)') &
      '', 'INFO: ', TRIM( ProgramName )
    WRITE(*,*)
    WRITE(*,'(A7,A6,ES10.3E2)') &
      '', 'X_D = ', X_D
    WRITE(*,*)
    WRITE(*,'(A7,A13,ES10.3E2)') '', 'Amplitude  = ', Amplitude
    WRITE(*,'(A7,A13,ES10.3E2)') '', 'Wavenumber = ', Wavenumber
    WRITE(*,*)

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
  
      DO iNodeX = 1, nDOFX
  
        iNodeX1 = NodeNumberTableX(1,iNodeX)
  
        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
      
        IF( X1 <= X_D )THEN

          uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = 3.60632d12 * Gram / Centimeter**3
          uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 7.425d4 * Kilometer / Second !2.629369d4
          uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP * Kilometer / Second
          uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP * Kilometer / Second
          uAF(iNodeX,iX1,iX2,iX3,iAF_P ) = 10.333333d31 * Erg / Centimeter**3
          uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) = 0.5_DP

        ELSE

          uPF(iNodeX,iX1,iX2,iX3,iPF_D ) &
            = ( 1.00d12 + Amplitude * SIN( Wavenumber * X1 ) ) * Gram / Centimeter**3
          uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP * Kilometer / Second
          uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP * Kilometer / Second
          uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP * Kilometer / Second
          uAF(iNodeX,iX1,iX2,iX3,iAF_P ) = 1.00d31 * Erg / Centimeter**3
          uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) = 0.5_DP

        END IF

        CALL ComputeTemperatureFromPressure &
               ( uPF(iNodeX,iX1,iX2,iX3,iPF_D ), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_P ), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_Ye), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_T ) )
  
        CALL ComputeThermodynamicStates_Primitive &
               ( uPF(iNodeX,iX1,iX2,iX3,iPF_D ), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_T ), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_Ye), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_E ), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_E ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_Ne) )
  
        CALL ComputeConserved_Euler_NonRelativistic &
               ( uPF(iNodeX,iX1,iX2,iX3,iPF_D ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V1), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V2), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V3), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_E ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_Ne), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_D ), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S1), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S2), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S3), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_E ), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_Ne), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33) )

        END DO

      END DO
      END DO
      END DO

  END SUBROUTINE InitiaizeFields_InitializeShockEntropyWaveInteraction1D


  SUBROUTINE InitializeFields_RiemannProblemCylindrical &
    ( RiemannProblemName_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: &
      RiemannProblemName_Option

    CHARACTER(32) :: RiemannProblemName
    INTEGER       :: iX1, iX2, iX3
    INTEGER       :: iNodeX, iNodeX1, iNodeX2
    REAL(DP)      :: X1, X2

    RiemannProblemName = 'CylindricalSod'
    IF( PRESENT( RiemannProblemName_Option ) ) &
       RiemannProblemName = TRIM( RiemannProblemName_Option )

    WRITE(*,*)
    WRITE(*,'(A4,A,A)') &
      '', 'Riemann Problem Name: ', TRIM( RiemannProblemName )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)
        iNodeX2 = NodeNumberTableX(2,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
        X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

        SELECT CASE ( TRIM( RiemannProblemName ) )

          CASE( 'CylindricalSod' )

            IF( SQRT( X1**2 + X2**2 ) <= 5.0_DP * Kilometer )THEN

              uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = 1.00d12 * Gram / Centimeter**3
              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP * Kilometer / Second
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP * Kilometer / Second
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP * Kilometer / Second
              uAF(iNodeX,iX1,iX2,iX3,iAF_P ) = 1.0d32 * Erg / Centimeter**3
              uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) = 0.4_DP

            ELSE

              uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = 1.25d11 * Gram / Centimeter**3
              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP * Kilometer / Second
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP * Kilometer / Second
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP * Kilometer / Second
              uAF(iNodeX,iX1,iX2,iX3,iAF_P ) = 1.0d31 * Erg / Centimeter**3
              uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) = 0.4_DP

            END IF

        END SELECT

      END DO

      CALL ComputeTemperatureFromPressure &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uAF(:,iX1,iX2,iX3,iAF_P), &
               uAF(:,iX1,iX2,iX3,iAF_Ye), uAF(:,iX1,iX2,iX3,iAF_T) )

      CALL ComputeThermodynamicStates_Primitive &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uAF(:,iX1,iX2,iX3,iAF_T ), &
               uAF(:,iX1,iX2,iX3,iAF_Ye), uPF(:,iX1,iX2,iX3,iPF_E ), &
               uAF(:,iX1,iX2,iX3,iAF_E ), uPF(:,iX1,iX2,iX3,iPF_Ne) )

      CALL ComputeConserved_Euler_NonRelativistic &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_V1), &
               uPF(:,iX1,iX2,iX3,iPF_V2), uPF(:,iX1,iX2,iX3,iPF_V3), &
               uPF(:,iX1,iX2,iX3,iPF_E ), uPF(:,iX1,iX2,iX3,iPF_Ne), &
               uCF(:,iX1,iX2,iX3,iCF_D ), uCF(:,iX1,iX2,iX3,iCF_S1), &
               uCF(:,iX1,iX2,iX3,iCF_S2), uCF(:,iX1,iX2,iX3,iCF_S3), &
               uCF(:,iX1,iX2,iX3,iCF_E ), uCF(:,iX1,iX2,iX3,iCF_Ne), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33) )

    END DO
    END DO
    END DO

  END SUBROUTINE InitializeFields_RiemannProblemCylindrical


  SUBROUTINE InitializeFields_Implosion

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX, iNodeX1, iNodeX2
    REAL(DP) :: X1, X2
    REAL(DP), PARAMETER :: D_0  = 1.25d11 * ( Gram / Centimeter**3 )
    REAL(DP), PARAMETER :: P_0  = 1.0d31  * ( Erg / Centimeter**3 )
    REAL(DP), PARAMETER :: Ye_0 = 0.3_DP
    REAL(DP), PARAMETER :: D_1  = 1.0d12  * ( Gram / Centimeter**3 )
    REAL(DP), PARAMETER :: P_1  = 1.0d32  * ( Erg / Centimeter**3 )
    REAL(DP), PARAMETER :: Ye_1 = 0.3_DP


!    REAL(DP), PARAMETER :: D_0  = 1.25d13 * ( Gram / Centimeter**3 )
!    REAL(DP), PARAMETER :: P_0  = 1.0d32  * ( Erg / Centimeter**3 )
!    REAL(DP), PARAMETER :: Ye_0 = 1.35d-1
!    REAL(DP), PARAMETER :: D_1  = 1.0d14  * ( Gram / Centimeter**3 )
!    REAL(DP), PARAMETER :: P_1  = 1.0d33  * ( Erg / Centimeter**3 )
!    REAL(DP), PARAMETER :: Ye_1 = 1.5d-1

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)
        iNodeX2 = NodeNumberTableX(2,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
        X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

        IF( X1 + X2 .LT. 0.15_DP * Kilometer )THEN

           uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
             = D_0
           uAF(iNodeX,iX1,iX2,iX3,iAF_P) &
             = P_0
           uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) &
             = Ye_0

         ELSE

           uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
             = D_1
           uAF(iNodeX,iX1,iX2,iX3,iAF_P) &
             = P_1
           uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) &
             = Ye_1

         ENDIF

         uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
           = Zero
         uPF(iNodeX,iX1,iX2,iX3,iPF_V2) &
           = Zero
         uPF(iNodeX,iX1,iX2,iX3,iPF_V3) &
           = Zero

        CALL ComputeTemperatureFromPressure &
               ( uPF(iNodeX,iX1,iX2,iX3,iPF_D ), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_P ), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_Ye), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_T ) )

        CALL ComputeThermodynamicStates_Primitive &
               ( uPF(iNodeX,iX1,iX2,iX3,iPF_D ), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_T ), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_Ye), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_E ), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_E ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_Ne) )

        CALL ComputeConserved_Euler_NonRelativistic &
               ( uPF(iNodeX,iX1,iX2,iX3,iPF_D ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V1), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V2), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V3), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_E ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_Ne), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_D ), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S1), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S2), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S3), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_E ), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_Ne), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33) )

      END DO

    END DO
    END DO
    END DO

  END SUBROUTINE InitializeFields_Implosion


  SUBROUTINE InitializeFields_Jet

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX, iNodeX1, iNodeX2
    REAL(DP) :: X1, X2

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)
        iNodeX2 = NodeNumberTableX(2,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 ) !* Kilometer
        X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 ) !* Kilometer

        IF( X1 .LE. Half * Kilometer .AND. X2 .LE. Half * Kilometer )THEN

          !SW
          uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = 0.80d12 * Gram / Centimeter**3
          uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP * Kilometer / Second
          uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP * Kilometer / Second
          uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP * Kilometer / Second
          uAF(iNodeX,iX1,iX2,iX3,iAF_P ) = 1.0d32 * Erg / Centimeter**3
          uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) = 0.3_DP

        ELSE IF ( X1 .LE. Half * Kilometer .AND. X2 .GT. Half * Kilometer )THEN

          !NW
          uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = 1.0d12 * Gram / Centimeter**3
          uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 7.275d4 * Kilometer / Second
          uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP * Kilometer / Second
          uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP * Kilometer / Second
          uAF(iNodeX,iX1,iX2,iX3,iAF_P ) = 1.0d32 * Erg / Centimeter**3
          uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) = 0.3_DP

        ELSE IF( X1 .GT. Half * Kilometer .AND. X2 .GT. Half * Kilometer )THEN

          !NE
          uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = 0.5313d12 * Gram / Centimeter**3
          uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP * Kilometer / Second
          uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP * Kilometer / Second
          uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP * Kilometer / Second
          uAF(iNodeX,iX1,iX2,iX3,iAF_P ) = 0.4d32 * Erg / Centimeter**3
          uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) = 0.3_DP

        ELSE

          !SE
          uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = 1.0d12 * Gram / Centimeter**3
          uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP * Kilometer / Second
          uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 7.275d4 * Kilometer / Second
          uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP * Kilometer / Second
          uAF(iNodeX,iX1,iX2,iX3,iAF_P ) = 1.0d32 * Erg / Centimeter**3
          uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) = 0.3_DP

        END IF

      END DO

      CALL ComputeTemperatureFromPressure &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uAF(:,iX1,iX2,iX3,iAF_P), &
               uAF(:,iX1,iX2,iX3,iAF_Ye), uAF(:,iX1,iX2,iX3,iAF_T) )

      CALL ComputeThermodynamicStates_Primitive &
             ( uPF(:,iX1,iX2,iX3,iPF_D),  uAF(:,iX1,iX2,iX3,iAF_T), &
               uAF(:,iX1,iX2,iX3,iAF_Ye), uPF(:,iX1,iX2,iX3,iPF_E ),&
               uAF(:,iX1,iX2,iX3,iAF_E),  uPF(:,iX1,iX2,iX3,iPF_Ne) )

      CALL ComputeConserved_Euler_NonRelativistic &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_V1), &
               uPF(:,iX1,iX2,iX3,iPF_V2), uPF(:,iX1,iX2,iX3,iPF_V3), &
               uPF(:,iX1,iX2,iX3,iPF_E ), uPF(:,iX1,iX2,iX3,iPF_Ne), &
               uCF(:,iX1,iX2,iX3,iCF_D ), uCF(:,iX1,iX2,iX3,iCF_S1), &
               uCF(:,iX1,iX2,iX3,iCF_S2), uCF(:,iX1,iX2,iX3,iCF_S3), &
               uCF(:,iX1,iX2,iX3,iCF_E ), uCF(:,iX1,iX2,iX3,iCF_Ne), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33) )

    END DO
    END DO
    END DO

  END SUBROUTINE InitializeFields_Jet


  SUBROUTINE InitializeFields_GravitationalCollapse( ProgenitorFileName_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: ProgenitorFileName_Option

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX1, iNodeX2, iNodeX3, iNodeX
    REAL(DP) :: X1
    CHARACTER(LEN=32) :: ProgenitorFile
    TYPE(ProgenitorType1D) :: P1D

    WRITE(*,*)
    WRITE(*,'(A2,A6,A)') '', 'INFO: ', TRIM( ProgramName )
    WRITE(*,*)

    ProgenitorFile = 'WH07_15M_Sun.h5'
    IF ( PRESENT( ProgenitorFileName_Option) ) &
      ProgenitorFile = ProgenitorFileName_Option

    CALL ReadProgenitor1D( TRIM( ProgenitorFile ), P1D )

    ! --- Initialize Fluid Fields ---

    ASSOCIATE &
      ( R1D => P1D % Radius, &
        D1D => P1D % MassDensity, &
        V1D => P1D % RadialVelocity, &
        T1D => P1D % Temperature, &
        Y1D => P1D % ElectronFraction )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E1(1)

      DO iNodeX3 = 1, nNodesX(3)
      DO iNodeX2 = 1, nNodesX(2)
      DO iNodeX1 = 1, nNodesX(1)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

        iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

        uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
          = Interpolate1D( R1D, D1D, SIZE( R1D ), X1 )

        uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
          = Interpolate1D( R1D, V1D, SIZE( R1D ), X1 )

        uPF(iNodeX,iX1,iX2,iX3,iPF_V2) &
          = Zero

        uPF(iNodeX,iX1,iX2,iX3,iPF_V3) &
          = Zero

        uAF(iNodeX,iX1,iX2,iX3,iAF_T) &
          = Interpolate1D( R1D, T1D, SIZE( R1D ), X1 )

        uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) &
          = Interpolate1D( R1D, Y1D, SIZE( R1D ), X1 )

        CALL ComputeThermodynamicStates_Primitive &
               ( uPF(iNodeX,iX1,iX2,iX3,iPF_D ), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_T ), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_Ye), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_E ), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_E ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_Ne) )

        CALL ApplyEquationOfState &
                ( uPF(iNodeX,iX1,iX2,iX3,iPF_D ), &
                  uAF(iNodeX,iX1,iX2,iX3,iAF_T ), &
                  uAF(iNodeX,iX1,iX2,iX3,iAF_Ye), &
                  uAF(iNodeX,iX1,iX2,iX3,iAF_P ), &
                  uAF(iNodeX,iX1,iX2,iX3,iAF_S ), &
                  uAF(iNodeX,iX1,iX2,iX3,iAF_E ), &
                  uAF(iNodeX,iX1,iX2,iX3,iAF_Me), &
                  uAF(iNodeX,iX1,iX2,iX3,iAF_Mp), &
                  uAF(iNodeX,iX1,iX2,iX3,iAF_Mn), &
                  uAF(iNodeX,iX1,iX2,iX3,iAF_Xp), &
                  uAF(iNodeX,iX1,iX2,iX3,iAF_Xn), &
                  uAF(iNodeX,iX1,iX2,iX3,iAF_Xa), &
                  uAF(iNodeX,iX1,iX2,iX3,iAF_Xh), &
                  uAF(iNodeX,iX1,iX2,iX3,iAF_Gm) )

        CALL ComputeConserved_Euler_NonRelativistic &
               ( uPF(iNodeX,iX1,iX2,iX3,iPF_D ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V1), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V2), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V3), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_E ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_Ne), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_D ), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S1), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S2), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S3), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_E ), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_Ne), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33) )

      END DO
      END DO
      END DO

    END DO
    END DO
    END DO

    END ASSOCIATE ! R1D, etc

  END SUBROUTINE InitializeFields_GravitationalCollapse


  REAL(DP) FUNCTION Interpolate1D( x, y, n, xq )

    INTEGER,                INTENT(in) :: n
    REAL(DP), DIMENSION(n), INTENT(in) :: x, y
    REAL(DP),               INTENT(in) :: xq

    INTEGER :: i

    i = Locate( xq, x, n )

    IF( i == 0 )THEN

      ! --- Extrapolate Left ---

      Interpolate1D &
        = Interpolate1D_Linear( xq, x(1), x(2), y(1), y(2) )

    ELSE IF( i == n )THEN

      ! --- Extrapolate Right ---

      Interpolate1D &
        = Interpolate1D_Linear( xq, x(n-1), x(n), y(n-1), y(n) )

    ELSE

      Interpolate1D &
        = Interpolate1D_Linear( xq, x(i), x(i+1), y(i), y(i+1) )

    END IF

    RETURN

  END FUNCTION Interpolate1D

END MODULE InitializationModule
