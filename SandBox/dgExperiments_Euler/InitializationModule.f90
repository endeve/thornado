MODULE InitializationModule

  USE KindModule, ONLY: &
    DP, Zero, Half, One, TwoPi, FourPi
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
    uGF, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33
  USE FluidFieldsModule, ONLY: &
    uPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne, &
    uCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne
  USE EulerEquationsUtilitiesModule_Beta, ONLY: &
    ComputeConserved

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFields

CONTAINS


  SUBROUTINE InitializeFields &
    ( AdvectionProfile_Option, Direction_Option, RiemannProblemName_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: &
      AdvectionProfile_Option
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: &
      Direction_Option
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: &
      RiemannProblemName_Option

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

      CASE( 'KelvinHelmholtz' )

        CALL InitializeFields_KelvinHelmholtz

      CASE( 'Implosion' )

        CALL InitializeFields_Implosion   

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
    INTEGER       :: iNodeX, iNodeX1
    REAL(DP)      :: X1

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

            X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

            SELECT CASE ( TRIM( RiemannProblemName ) )

              CASE( 'Sod' )

                IF( X1 <= Half )THEN

                  uPF(iNodeX1,iX1,iX2,iX3,iPF_D)  = 1.0_DP
                  uPF(iNodeX1,iX1,iX2,iX3,iPF_V1) = 0.0_DP
                  uPF(iNodeX1,iX1,iX2,iX3,iPF_V2) = 0.0_DP
                  uPF(iNodeX1,iX1,iX2,iX3,iPF_V3) = 0.0_DP
                  uPF(iNodeX1,iX1,iX2,iX3,iPF_E)  = 1.0_DP / 0.4_DP

                ELSE

                  uPF(iNodeX1,iX1,iX2,iX3,iPF_D)  = 0.125_DP
                  uPF(iNodeX1,iX1,iX2,iX3,iPF_V1) = 0.0_DP
                  uPF(iNodeX1,iX1,iX2,iX3,iPF_V2) = 0.0_DP
                  uPF(iNodeX1,iX1,iX2,iX3,iPF_V3) = 0.0_DP
                  uPF(iNodeX1,iX1,iX2,iX3,iPF_E)  = 0.1_DP / 0.4_DP

                END IF

              CASE( 'SteadyContact' )

                IF( X1 <= Half )THEN

                  uPF(iNodeX1,iX1,iX2,iX3,iPF_D)  = 1.4_DP
                  uPF(iNodeX1,iX1,iX2,iX3,iPF_V1) = 0.0_DP
                  uPF(iNodeX1,iX1,iX2,iX3,iPF_V2) = 0.0_DP
                  uPF(iNodeX1,iX1,iX2,iX3,iPF_V3) = 0.0_DP
                  uPF(iNodeX1,iX1,iX2,iX3,iPF_E)  = 1.0_DP / 0.4_DP

                ELSE

                  uPF(iNodeX1,iX1,iX2,iX3,iPF_D)  = 1.0_DP
                  uPF(iNodeX1,iX1,iX2,iX3,iPF_V1) = 0.0_DP
                  uPF(iNodeX1,iX1,iX2,iX3,iPF_V2) = 0.0_DP
                  uPF(iNodeX1,iX1,iX2,iX3,iPF_V3) = 0.0_DP
                  uPF(iNodeX1,iX1,iX2,iX3,iPF_E)  = 1.0_DP / 0.4_DP

                END IF

              CASE( 'MovingContact' )

                IF( X1 <= Half )THEN

                  uPF(iNodeX1,iX1,iX2,iX3,iPF_D)  = 1.4_DP
                  uPF(iNodeX1,iX1,iX2,iX3,iPF_V1) = 0.1_DP
                  uPF(iNodeX1,iX1,iX2,iX3,iPF_V2) = 0.0_DP
                  uPF(iNodeX1,iX1,iX2,iX3,iPF_V3) = 0.0_DP
                  uPF(iNodeX1,iX1,iX2,iX3,iPF_E)  = 1.0_DP / 0.4_DP

                ELSE

                  uPF(iNodeX1,iX1,iX2,iX3,iPF_D)  = 1.0_DP
                  uPF(iNodeX1,iX1,iX2,iX3,iPF_V1) = 0.1_DP
                  uPF(iNodeX1,iX1,iX2,iX3,iPF_V2) = 0.0_DP
                  uPF(iNodeX1,iX1,iX2,iX3,iPF_V3) = 0.0_DP
                  uPF(iNodeX1,iX1,iX2,iX3,iPF_E)  = 1.0_DP / 0.4_DP

                END IF

              CASE( 'Toro' )

              CASE( 'InteractingBlasts' )

              CASE( 'ShockEntropyWave' )

            END SELECT

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

  END SUBROUTINE InitializeFields_RiemannProblem


  SUBROUTINE InitializeFields_KelvinHelmholtz

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX, iNodeX1, iNodeX2
    REAL(DP) :: X1, X2, D_M, V_M
    REAL(DP), PARAMETER :: D_1 = 1.0_DP
    REAL(DP), PARAMETER :: D_2 = 2.0_DP
    REAL(DP), PARAMETER :: L = 0.025_DP
    REAL(DP), PARAMETER :: V_1 = + 0.5_DP
    REAL(DP), PARAMETER :: V_2 = - 0.5_DP

    D_M = Half * ( D_1 - D_2 )
    V_M = Half * ( V_1 - V_2 )

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          DO iNodeX = 1, nDOFX

            iNodeX1 = NodeNumberTableX(1,iNodeX)
            iNodeX2 = NodeNumberTableX(2,iNodeX)

            X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
            X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

            IF(     ( X2 .GE. 0.00_DP ) .AND. ( X2 .LT. 0.25_DP ) )THEN

              uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
                = D_1 - D_M * EXP( ( X2 - 0.25_DP ) / L )
              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
                = V_1 - V_M * EXP( ( X2 - 0.25_DP ) / L )

            ELSEIF( ( X2 .GE. 0.25_DP ) .AND. ( X2 .LT. 0.50_DP ) )THEN

              uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
                = D_2 + D_M * EXP( ( 0.25_DP - X2 ) / L )
              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
                = V_2 + V_M * EXP( ( 0.25_DP - X2 ) / L )

            ELSEIF( ( X2 .GE. 0.50_DP ) .AND. ( X2 .LT. 0.75_DP ) )THEN

              uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
                = D_2 + D_M * EXP( ( X2 - 0.75_DP ) / L )
              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
                = V_2 + V_M * EXP( ( X2 - 0.75_DP ) / L )

            ELSEIF( ( X2 .GE. 0.75_DP ) .AND. ( X2 .LT. 1.00_DP ) )THEN

              uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
                = D_1 - D_M * EXP( ( 0.75_DP - X2 ) / L )
              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
                = V_1 - V_M * EXP( ( 0.75_DP - X2 ) / L )

            END IF

            uPF(iNodeX,iX1,iX2,iX3,iPF_V2) &
              = 0.01_DP * SIN( FourPi * X1 )
            uPF(iNodeX,iX1,iX2,iX3,iPF_V3) &
              = Zero
            uPF(iNodeX,iX1,iX2,iX3,iPF_E) &
              = 3.75_DP

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

  END SUBROUTINE InitializeFields_KelvinHelmholtz


  SUBROUTINE InitializeFields_Implosion

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX, iNodeX1, iNodeX2
    REAL(DP) :: X1, X2, D_M, V_M
    REAL(DP), PARAMETER :: D_0 = 0.125_DP
    REAL(DP), PARAMETER :: D_1 = 1.0_DP
    REAL(DP), PARAMETER :: E_0 = + 2.5_DP
    REAL(DP), PARAMETER :: E_1 = - 0.35_DP


    DO iX2 = 1, nX(2)
      DO iX1 = 1, nX(1)

        DO iNodeX = 1, nDOFX

         iNodeX1 = NodeNumberTableX(1,iNodeX)
         iNodeX2 = NodeNumberTableX(2,iNodeX)

         X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
         X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

         IF( X1 + X2 .LE. 0.15 )THEN
   
           uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
             = D_0
           uPF(iNodeX,iX1,iX2,iX3,iPF_E) &
             = E_0
           
         ELSEIF( X1 + X2 .GE. 0.15 )THEN

           uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
             = D_1
           uPF(iNodeX,iX1,iX2,iX3,iPF_E) &
             = E_1

         ENDIF

         uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
           = Zero
         uPF(iNodeX,iX1,iX2,iX3,iPF_V2) &
           = Zero
         uPF(iNodeX,iX1,iX2,iX3,iPF_V3) &
           = Zero

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

  END SUBROUTINE InitializeFields_Implosion

END MODULE InitializationModule
