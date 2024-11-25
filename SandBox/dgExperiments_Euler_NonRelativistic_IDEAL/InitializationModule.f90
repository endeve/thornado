MODULE InitializationModule

  USE KindModule, ONLY: &
    DP, Zero, Half, One, Three, &
    Pi, TwoPi, FourPi
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
    uAF
  USE Euler_UtilitiesModule_NonRelativistic, ONLY: &
    ComputeConserved_Euler_NonRelativistic

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFields


CONTAINS


  SUBROUTINE InitializeFields &
    ( AdvectionProfile_Option, Direction_Option, RiemannProblemName_Option, &
      SedovEnergy_Option, nDetCells_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: &
      AdvectionProfile_Option
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: &
      Direction_Option
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: &
      RiemannProblemName_Option

    REAL(DP), INTENT(in), OPTIONAL :: SedovEnergy_Option

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

      CASE ( 'RiemannProblemSpherical' )

        CALL InitializeFields_RiemannProblemSpherical

      CASE ( 'SphericalSedov' )

        CALL InitializeFields_SphericalSedov &
               ( SedovEnergy_Option &
                   = SedovEnergy_Option, &
                 nDetCells_Option &
                   = nDetCells_Option )

      CASE( 'IsentropicVortex' )

        CALL InitializeFields_IsentropicVortex

      CASE( 'KelvinHelmholtz' )

        CALL InitializeFields_KelvinHelmholtz

      CASE( 'RayleighTaylor' )

        CALL InitializeFields_RayleighTaylor

      CASE( 'Implosion' )

        CALL InitializeFields_Implosion

      CASE( 'Explosion' )

        CALL InitializeFields_Explosion

    END SELECT


  END SUBROUTINE InitializeFields


  SUBROUTINE InitializeFields_Advection &
    ( AdvectionProfile_Option, Direction_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: &
      AdvectionProfile_Option
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: &
      Direction_Option

    REAL(DP), PARAMETER :: Amp = 0.1_DP

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
                    = One + Amp * SIN( TwoPi * X1 )

                ELSEIF( TRIM( Direction ) .EQ. 'Y' )THEN

                  uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
                    = One + Amp * SIN( TwoPi * X2 )

                END IF

              CASE( 'QuarticSineWave' )

                IF( TRIM( Direction ) .EQ. 'X' )THEN

                  uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
                    = One + Amp * SIN( Pi * X1 )**4

                ELSEIF( TRIM( Direction ) .EQ. 'Y' )THEN

                  uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
                    = One + Amp * SIN( Pi * X2 )**4

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

                  uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 1.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_E)  = 1.0_DP / 0.4_DP

                ELSE

                  uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 0.125_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_E)  = 0.1_DP / 0.4_DP

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

  END SUBROUTINE InitializeFields_RiemannProblem


  SUBROUTINE InitializeFields_RiemannProblemSpherical

    INTEGER       :: iX1, iX2, iX3
    INTEGER       :: iNodeX, iNodeX1
    REAL(DP)      :: X1

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          DO iNodeX = 1, nDOFX

            iNodeX1 = NodeNumberTableX(1,iNodeX)

            X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

            IF( X1 <= One )THEN

              uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 1.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_E)  = 1.0_DP / 0.4_DP

            ELSE

              uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 0.125_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_E)  = 0.1_DP / 0.4_DP

            END IF

          END DO

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


  SUBROUTINE InitializeFields_SphericalSedov &
    ( SedovEnergy_Option, nDetCells_Option )

    REAL(DP), INTENT(in), OPTIONAL :: SedovEnergy_Option
    INTEGER,  INTENT(in), OPTIONAL :: nDetCells_Option

    INTEGER  :: iX1, iX2, iX3, nDetCells
    INTEGER  :: iNodeX, iNodeX1
    REAL(DP) :: X1, R_0, E_0

    nDetCells = 1
    IF( PRESENT( nDetCells_Option ) ) &
      nDetCells = nDetCells_Option

    E_0 = 1.0_DP
    IF( PRESENT( SedovEnergy_Option ) ) &
      E_0 = SedovEnergy_Option

    R_0 = REAL( nDetCells ) * MeshX(1) % Width(1)

    WRITE(*,*)
    WRITE(*,'(A7,A,ES10.3E2,A2,I2.2,A2,ES10.3E2)') &
      '', 'E_0, # of Detonation Cells, R_0 = ', E_0, ', ', nDetCells, ',', R_0

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          DO iNodeX = 1, nDOFX

            iNodeX1 = NodeNumberTableX(1,iNodeX)

            X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

            uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 1.0_DP
            uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
            uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
            uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP

            IF( X1 < R_0 )THEN

              uPF(iNodeX,iX1,iX2,iX3,iPF_E) = 3.0_DP * E_0 &
                                            / ( 4.0_DP * Pi * R_0**3 )

            ELSE

              uPF(iNodeX,iX1,iX2,iX3,iPF_E) = 1.0d-5

            END IF

          END DO

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

  END SUBROUTINE


  SUBROUTINE InitializeFields_IsentropicVortex

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX, iNodeX1, iNodeX2
    REAL(DP) :: X1, X2, R

    REAL(DP), PARAMETER :: Beta  = 5.0_DP
    REAL(DP), PARAMETER :: Gamma = 1.4_DP

    DO iX3 = 1, nX(3)
    DO iX2 = 1, nX(2)
    DO iX1 = 1, nX(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)
        iNodeX2 = NodeNumberTableX(2,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
        X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

        R  = SQRT( X1**2 + X2**2 )

        uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
          = ( One - ( Gamma - One ) * Beta**2 &
                    / ( 8.0_DP * Gamma * Pi**2 ) * EXP( One - R**2 ) &
            )**( One / ( Gamma - One ) )

        uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
          = One &
            - X2 * ( Beta / TwoPi ) * EXP( Half * ( One - R**2 ) )

        uPF(iNodeX,iX1,iX2,iX3,iPF_V2) &
          = One &
            + X1 * ( Beta / TwoPi ) * EXP( Half * ( One - R**2 ) )

        uPF(iNodeX,iX1,iX2,iX3,iPF_V3) &
          = 0.0_DP

        uPF(iNodeX,iX1,iX2,iX3,iPF_E) &
          = uPF(iNodeX,iX1,iX2,iX3,iPF_D)**Gamma / ( Gamma - One )

      END DO

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

  END SUBROUTINE InitializeFields_IsentropicVortex


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

  END SUBROUTINE InitializeFields_KelvinHelmholtz


  SUBROUTINE InitializeFields_RayleighTaylor

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX, iNodeX1, iNodeX2
    REAL(DP) :: X1, X2

    REAL(DP), PARAMETER :: D_1  = 1.0_DP
    REAL(DP), PARAMETER :: D_2  = 2.0_DP
    REAL(DP), PARAMETER :: X1_b = - 0.75_DP
    REAL(DP), PARAMETER :: E_b  = 2.5_DP / 0.4_DP
    REAL(DP), PARAMETER :: g    = 0.1_DP

    DO iX3 = 1, nX(3)
       DO iX2 = 1, nX(2)
          DO iX1 = 0, nX(1) + 1

             DO iNodeX = 1, nDOFX

                iNodeX1 = NodeNumberTableX(1,iNodeX)
                iNodeX2 = NodeNumberTableX(2,iNodeX)

                X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
                X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

                IF( X1 .LE. 0.0_DP )THEN

                  uPF(iNodeX,iX1,iX2,iX3,iPF_D) = D_1

                ELSE

                  uPF(iNodeX,iX1,iX2,iX3,iPF_D) = D_2

                END IF

                uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
                  = 0.0025_DP * ( One + COS( FourPi * X2 ) ) &
                      * ( One + COS( Three * Pi * X1 ) )
                uPF(iNodeX,iX1,iX2,iX3,iPF_V2) &
                  = 0.0_DP
                uPF(iNodeX,iX1,iX2,iX3,iPF_V3) &
                  = 0.0_DP

                uPF(iNodeX,iX1,iX2,iX3,iPF_E) &
                  = E_b - g * D_1 * ( MIN( X1, 0.0_DP ) - X1_b ) &
                        - g * D_2 * MAX( X1, 0.0_DP )

                uGF(iNodeX,iX1,iX2,iX3,iGF_Phi_N) &
                  = g * X1

             END DO

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


  END SUBROUTINE InitializeFields_RayleighTaylor


  SUBROUTINE InitializeFields_Implosion

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX, iNodeX1, iNodeX2
    REAL(DP) :: X1, X2
    REAL(DP), PARAMETER :: D_0 = 0.125_DP
    REAL(DP), PARAMETER :: E_0 = 0.350_DP
    REAL(DP), PARAMETER :: D_1 = 1.000_DP
    REAL(DP), PARAMETER :: E_1 = 2.500_DP

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          DO iNodeX = 1, nDOFX

           iNodeX1 = NodeNumberTableX(1,iNodeX)
           iNodeX2 = NodeNumberTableX(2,iNodeX)

           X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
           X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

           IF( X1 + X2 .LT. 0.15_DP )THEN

             uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
               = D_0
             uPF(iNodeX,iX1,iX2,iX3,iPF_E) &
               = E_0

           ELSE

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

  END SUBROUTINE InitializeFields_Implosion


  SUBROUTINE InitializeFields_Explosion

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX, iNodeX1, iNodeX2
    REAL(DP) :: X1, X2, R
    INTEGER,  PARAMETER :: p = 150
    REAL(DP), PARAMETER :: R_0 = 0.4_DP
    REAL(DP), PARAMETER :: D_O = 0.125_DP
    REAL(DP), PARAMETER :: E_O = 0.250_DP
    REAL(DP), PARAMETER :: D_I = 1.000_DP
    REAL(DP), PARAMETER :: E_I = 2.500_DP

    DO iX3 = 1, nX(3)
    DO iX2 = 1, nX(2)
    DO iX1 = 1, nX(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)
        iNodeX2 = NodeNumberTableX(2,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
        X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

        R = SQRT( X1**2 + X2**2 )

        uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
          = ( D_I + (R/R_0)**p * D_O ) / ( One + (R/R_0)**p )

        uPF(iNodeX,iX1,iX2,iX3,iPF_E) &
          = ( E_I + (R/R_0)**p * E_O ) / ( One + (R/R_0)**p )

!!$        IF( R <= R_0 )THEN
!!$
!!$          uPF(iNodeX,iX1,iX2,iX3,iPF_D) = D_I
!!$          uPF(iNodeX,iX1,iX2,iX3,iPF_E) = E_I
!!$
!!$        ELSE
!!$
!!$          uPF(iNodeX,iX1,iX2,iX3,iPF_D) = D_O
!!$          uPF(iNodeX,iX1,iX2,iX3,iPF_E) = E_O
!!$
!!$        ENDIF

        uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = Zero
        uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = Zero
        uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = Zero

      END DO

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

  END SUBROUTINE InitializeFields_Explosion


END MODULE InitializationModule
