MODULE InitializationModule

  USE KindModule, ONLY: &
    DP, Zero, Half, One, Three, &
    Pi, TwoPi, FourPi
  USE UnitsModule, ONLY: &
    Gram, Centimeter, &
    Kilometer, Erg, Second, Kelvin
  USE UtilitiesModule, ONLY: &
    Locate
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    nX, nNodesX, &
    nDOFX
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
    uAF, iAF_P, iAF_Ye, iAF_T, iAF_E
  USE Euler_UtilitiesModule, ONLY: &
    ComputeConserved
  USE EquationOfStateModule, ONLY: &
    ComputeTemperatureFromPressure, &
    ComputeThermodynamicStates_Primitive

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFields

CONTAINS


  SUBROUTINE InitializeFields &
    ( AdvectionProfile_Option, Direction_Option, RiemannProblemName_Option, &
      nDetCells_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: &
      AdvectionProfile_Option
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: &
      Direction_Option
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: &
      RiemannProblemName_Option

    INTEGER, INTENT(in), OPTIONAL :: nDetCells_Option

    WRITE(*,*)
    WRITE(*,'(A2,A6,A)') '', 'INFO: ', TRIM( ProgramName )

    SELECT CASE ( TRIM( ProgramName ) )

      CASE( 'Advection' )

        CALL InitializeFields_Advection &
               ( AdvectionProfile_Option &
                   = AdvectionProfile_Option, &
                 Direction_Option &
                   = Direction_Option )

      CASE ( 'RiemannProblem' )

        CALL InitializeFields_RiemannProblem &
               ( RiemannProblemName_Option &
                   = RiemannProblemName_Option )

      CASE ( 'Jet' )

        CALL InitializeFields_Jet(  )

    END SELECT


  END SUBROUTINE InitializeFields

  SUBROUTINE InitializeFields_Advection &
    ( AdvectionProfile_Option, Direction_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: &
      AdvectionProfile_Option
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: &
      Direction_Option

    CHARACTER(32) :: AdvectionProfile
    CHARACTER(32) :: Direction
    INTEGER       :: iX1, iX2, iX3
    INTEGER       :: iNodeX, iNodeX1, iNodeX2
    REAL(DP)      :: X1, X2

    AdvectionProfile = 'SineWave'
    IF( PRESENT( AdvectionProfile_Option ) ) &
      AdvectionProfile = TRIM( AdvectionProfile_Option )

    Direction = 'X'
    IF( PRESENT( Direction_Option ) ) &
      Direction = Direction_Option

    WRITE(*,*)
    WRITE(*,'(A4,A20,A)') &
      '', 'Advection Profile: ', TRIM( AdvectionProfile )
    WRITE(*,'(A4,A20,A)') &
      '', 'Direction: ', TRIM( Direction )

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          DO iNodeX = 1, nDOFX

            iNodeX1 = NodeNumberTableX(1,iNodeX)
            iNodeX2 = NodeNumberTableX(2,iNodeX)

            X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
            X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

            SELECT CASE ( TRIM( AdvectionProfile ) )

              CASE( 'SineWave' )

                IF( TRIM( Direction ) .EQ. 'X' )THEN

                  uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
                    = One + Half * SIN( TwoPi * X1 )

                ELSEIF( TRIM( Direction ) .EQ. 'Y' )THEN

                  uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
                    = One + Half * SIN( TwoPi * X2 )

                END IF

              CASE( 'TopHat' )

                IF( TRIM( Direction ) .EQ. 'X' )THEN

                  IF( ABS( X1 - 0.5_DP ) .LE. 0.25_DP )THEN

                    uPF(iNodeX,iX1,iX2,iX3,iPF_D) = 2.0_DP

                  ELSE

                    uPF(iNodeX,iX1,iX2,iX3,iPF_D) = 1.0_DP

                  END IF

                ELSEIF( TRIM( Direction ) .EQ. 'Y' )THEN

                  IF( ABS( X2 - 0.5_DP ) .LE. 0.25_DP )THEN

                    uPF(iNodeX,iX1,iX2,iX3,iPF_D) = 2.0_DP

                  ELSE

                    uPF(iNodeX,iX1,iX2,iX3,iPF_D) = 1.0_DP

                  END IF

                ELSEIF( TRIM( Direction ) .EQ. 'XY' )THEN

                  IF( ( ABS( X1 - 0.5_DP ) .LE. 0.25_DP ) &
                      .AND. ( ABS( X2 - 0.5_DP ) .LE. 0.25_DP ) )THEN

                    uPF(iNodeX,iX1,iX2,iX3,iPF_D) = 2.0_DP

                  ELSE

                    uPF(iNodeX,iX1,iX2,iX3,iPF_D) = 1.0_DP

                  END IF

                END IF

            END SELECT

            IF( TRIM( Direction ) .EQ. 'X' )THEN

              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 1.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP

            ELSEIF( TRIM( Direction ) .EQ. 'Y' )THEN

              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 1.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP

            ELSEIF( TRIM( Direction ) .EQ. 'XY' )THEN

              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 1.0_DP / SQRT( 2.0_DP )
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.5_DP / SQRT( 2.0_DP )
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP

            END IF

            uPF(iNodeX,iX1,iX2,iX3,iPF_E) = 1.0d-2

          END DO

          CALL ComputeConserved &
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

  END SUBROUTINE InitializeFields_Advection


  SUBROUTINE InitializeFields_RiemannProblem &
    ( RiemannProblemName_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: &
      RiemannProblemName_Option

    CHARACTER(32) :: RiemannProblemName
    INTEGER       :: iX1, iX2, iX3
    INTEGER       :: iNodeX, iNodeX1, iNodeX2
    REAL(DP)      :: X1, X2

    RiemannProblemName = 'Sod'
    IF( PRESENT( RiemannProblemName_Option ) ) &
       RiemannProblemName = TRIM( RiemannProblemName_Option )

    WRITE(*,*)
    WRITE(*,'(A4,A,A)') &
      '', 'Riemann Problem Name: ', TRIM( RiemannProblemName )

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          DO iNodeX = 1, nDOFX

            iNodeX1 = NodeNumberTableX(1,iNodeX)
            iNodeX2 = NodeNumberTableX(2,iNodeX)

            X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
            X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

            SELECT CASE ( TRIM( RiemannProblemName ) )

              CASE( 'Sod' )

                IF( X1 <= 0.0_DP )THEN
                 uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 1.0d14 * Gram / Centimeter**3
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP * Kilometer / Second
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP * Kilometer / Second
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP * Kilometer / Second
                 uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = 1.0d33 * Erg / Centimeter**3
                 uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) = 1.5d-1
               ELSE
                 uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 1.25d13 * Gram / Centimeter**3
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP * Kilometer / Second
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP * Kilometer / Second
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP * Kilometer / Second
                 uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = 1.0d32 * Erg / Centimeter**3
                 uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) = 1.35d-1

                ! IF( X1 <= Half )THEN
                !   uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 1.0d15 * Gram / Centimeter**3
                !   uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP * Kilometer / Second
                !   uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP * Kilometer / Second
                !   uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP * Kilometer / Second
                !   uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = 4.18d35 * Erg / Centimeter**3
                !   uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) = 3.0d-1
  
                ! ELSE
  
                !   uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 1.25d14 * Gram / Centimeter**3
                !   uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP * Kilometer / Second
                !   uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP * Kilometer / Second
                !   uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP * Kilometer / Second
                !   uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = 4.5d34 * Erg / Centimeter**3
                !   uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) = 3.00d-1  

                END IF

              CASE( 'ChimeraProfile' )

                IF( X1 <= 0.0_DP )THEN

                  uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 1.0d14 * Gram / Centimeter**3
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP * Kilometer / Second
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP * Kilometer / Second
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP * Kilometer / Second
                  uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = 1.0d33 * Erg / Centimeter**3!1.5d32 * Erg / Centimeter**3
                  !uAF(iNodeX,iX1,iX2,iX3,iAF_T)  = 1.26d11 * Kelvin
                  uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) = 1.5d-1

                ELSE

                  uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 1.25d13 * Gram / Centimeter**3
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP * Kilometer / Second
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP * Kilometer / Second
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP * Kilometer / Second
                  uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = 1.0d32 * Erg / Centimeter**3
                  !uAF(iNodeX,iX1,iX2,iX3,iAF_T)  = 6.9d9 * Kelvin
                  uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) = 1.35d-1

                END IF

              CASE( 'SteadyContact' )

                IF( X1 <= Half )THEN

                  uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 1.4_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_E)  = 1.0_DP / 0.4_DP

                ELSE

                  uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 1.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_E)  = 1.0_DP / 0.4_DP

                END IF

              CASE( 'MovingContact' )

                IF( X1 <= Half )THEN

                  uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 1.4_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.1_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_E)  = 1.0_DP / 0.4_DP

                ELSE

                  uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 1.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.1_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_E)  = 1.0_DP / 0.4_DP

                END IF

              CASE( 'Toro' )

              CASE( 'InteractingBlasts' )

              CASE( 'ShockEntropyWave' )

                IF( X1 < -4.0_DP )THEN

                  uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 3.857143_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 2.629369_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_E)  = 10.333333_DP / 0.4_DP

                ELSE

                  uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 1.0_DP + 0.2 * SIN(5*X1)
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_E)  = 1.0_DP / 0.4_DP

                END IF

            END SELECT

          END DO

          ! CALL TempFromPRessure
          CALL ComputeTemperatureFromPressure &
                 ( uPF(:,iX1,iX2,iX3,iPF_D ), uAF(:,iX1,iX2,iX3,iAF_P), &
                   uAF(:,iX1,iX2,iX3,iAF_Ye), uAF(:,iX1,iX2,iX3,iAF_T) )

          CALL ComputeThermodynamicStates_Primitive &
                 ( uPF(:,iX1,iX2,iX3,iPF_D),  uAF(:,iX1,iX2,iX3,iAF_T), &
                   uAF(:,iX1,iX2,iX3,iAF_Ye), uPF(:,iX1,iX2,iX3,iPF_E ),&
                   uAF(:,iX1,iX2,iX3,iAF_E),  uPF(:,iX1,iX2,iX3,iPF_Ne) )
                   ! iAF_E = Em. Ev = iPF_E

          CALL ComputeConserved &
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

  END SUBROUTINE InitializeFields_RiemannProblem

  SUBROUTINE InitializeFields_Jet

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX, iNodeX1, iNodeX2
    REAL(DP) :: X1, X2

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          DO iNodeX = 1, nDOFX

           iNodeX1 = NodeNumberTableX(1,iNodeX)
           iNodeX2 = NodeNumberTableX(2,iNodeX)

           X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 ) !* Kilometer
           X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 ) !* Kilometer

           IF( X1 .LE. Half*Kilometer .AND. X2 .LE. Half*Kilometer )THEN

            !SW
             uPF(iNodeX,iX1,iX2,iX3,iPF_D)   = 0.80d12 * Gram / Centimeter**3
             uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP * Kilometer / Second
             uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP * Kilometer / Second
             uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP * Kilometer / Second
             uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = 1.0d32 * Erg / Centimeter**3
             uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) = 0.3_DP

           !NW
           ELSE IF ( X1 .LE. Half*Kilometer .AND. X2 .GT. Half*Kilometer )THEN

             uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 1.0d12 * Gram / Centimeter**3
             uPF(iNodeX,iX1,iX2,iX3,iPF_V1)  = 7.275d4 * Kilometer / Second
             uPF(iNodeX,iX1,iX2,iX3,iPF_V2)  = 0.0_DP * Kilometer / Second
             uPF(iNodeX,iX1,iX2,iX3,iPF_V3)  = 0.0_DP * Kilometer / Second
             uAF(iNodeX,iX1,iX2,iX3,iAF_P)   = 1.0d32 * Erg / Centimeter**3
             uAF(iNodeX,iX1,iX2,iX3,iAF_Ye)  = 0.3_DP

           !NE
           ELSE IF( X1 .GT. Half*Kilometer .AND. X2 .GT. Half*Kilometer )THEN

             uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 0.5313d12 * Gram / Centimeter**3
             uPF(iNodeX,iX1,iX2,iX3,iPF_V1)  = 0.0_DP * Kilometer / Second
             uPF(iNodeX,iX1,iX2,iX3,iPF_V2)  = 0.0_DP * Kilometer / Second
             uPF(iNodeX,iX1,iX2,iX3,iPF_V3)  = 0.0_DP * Kilometer / Second
             uAF(iNodeX,iX1,iX2,iX3,iAF_P)   = 0.4d32 * Erg / Centimeter**3
             uAF(iNodeX,iX1,iX2,iX3,iAF_Ye)  = 0.3_DP

           !SE 
           ELSE

             uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 1.0d12 * Gram / Centimeter**3
             uPF(iNodeX,iX1,iX2,iX3,iPF_V1)  = 0.0_DP * Kilometer / Second
             uPF(iNodeX,iX1,iX2,iX3,iPF_V2)  = 7.275d4 * Kilometer / Second
             uPF(iNodeX,iX1,iX2,iX3,iPF_V3)  = 0.0_DP * Kilometer / Second
             uAF(iNodeX,iX1,iX2,iX3,iAF_P)   = 1.0d32 * Erg / Centimeter**3
             uAF(iNodeX,iX1,iX2,iX3,iAF_Ye)  = 0.3_DP

          END IF

          END DO

         ! CALL TempFromPRessure
         CALL ComputeTemperatureFromPressure &
                ( uPF(:,iX1,iX2,iX3,iPF_D ), uAF(:,iX1,iX2,iX3,iAF_P), &
                  uAF(:,iX1,iX2,iX3,iAF_Ye), uAF(:,iX1,iX2,iX3,iAF_T) )

         CALL ComputeThermodynamicStates_Primitive &
                ( uPF(:,iX1,iX2,iX3,iPF_D),  uAF(:,iX1,iX2,iX3,iAF_T), &
                  uAF(:,iX1,iX2,iX3,iAF_Ye), uPF(:,iX1,iX2,iX3,iPF_E ),&
                  uAF(:,iX1,iX2,iX3,iAF_E),  uPF(:,iX1,iX2,iX3,iPF_Ne) )
                  ! iAF_E = Em. Ev = iPF_E

         CALL ComputeConserved &
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

END MODULE InitializationModule
