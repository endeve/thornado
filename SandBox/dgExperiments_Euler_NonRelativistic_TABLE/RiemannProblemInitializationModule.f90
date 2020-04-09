MODULE RiemannProblemInitializationModule

  USE KindModule, ONLY: &
    DP
  USE UnitsModule, ONLY: &
    UnitsDisplay, &
    Gram, &
    Centimeter, &
    Kilometer, &
    Second, &
    Erg
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    nX, nNodesX
  USE UtilitiesModule, ONLY: &
    NodeNumberX
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE FluidFieldsModule, ONLY: &
    uCF, nCF, &
    uPF, nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne, &
    uAF, iAF_P, iAF_T, iAF_Ye, iAF_S, iAF_E, iAF_Me, iAF_Mp, &
    iAF_Mn, iAF_Xp, iAF_Xn, iAF_Xa, iAF_Xh, iAF_Gm, iAF_Cs
  USE InputOutputModule, ONLY: &
    ReadFluidFieldsRestart1D
  USE EquationOfStateModule, ONLY: &
    ComputeTemperatureFromPressure, &
    ComputeThermodynamicStates_Primitive, &
    ApplyEquationOfState
  USE EulerEquationsUtilitiesModule, ONLY: &
    Conserved, &
    ComputeConserved

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeRiemannProblem1D
  PUBLIC :: InitializeShockEntropyWaveInteraction1D

CONTAINS


  SUBROUTINE InitializeRiemannProblem1D &
               ( D_L, V_L, P_L, Ye_L, D_R, V_R, P_R, Ye_R, StartTime, &
                 X_D_Option, RestartNumber_Option )

    REAL(DP),               INTENT(in)  :: D_L, P_L, Ye_L, D_R, P_R, Ye_R
    REAL(DP), DIMENSION(3), INTENT(in)  :: V_L, V_R
    REAL(DP),               INTENT(out) :: StartTime
    REAL(DP),               INTENT(in), OPTIONAL :: X_D_Option
    INTEGER,                INTENT(in), OPTIONAL :: RestartNumber_Option

    LOGICAL  :: Restart
    INTEGER  :: iX1, iX2, iX3, RestartNumber
    INTEGER  :: iNodeX1, iNodeX2, iNodeX3, iNode
    REAL(DP) :: X_D, X1

    X_D = 0.5_DP * Kilometer
    IF( PRESENT( X_D_Option ) ) &
      X_D = X_D_Option

    Restart       = .FALSE.
    RestartNumber = 0
    IF( PRESENT( RestartNumber_Option ) )THEN
      Restart       = .TRUE.
      RestartNumber = RestartNumber_Option
    END IF

    StartTime = 0.0_DP

    ASSOCIATE( U => UnitsDisplay )

    WRITE(*,*)
    WRITE(*,'(A2,A6,A)') &
      '', 'INFO: ', TRIM( ProgramName )
    WRITE(*,*)
    WRITE(*,'(A7,A7,ES10.3E2,A1,A,A22,A7,ES10.3E2,A1,A)') &
      '', 'D_L = ', D_L / U % MassDensityUnit, &
      '', TRIM( U % MassDensityLabel ), &
      '', 'D_R = ', D_R / U % MassDensityUnit, &
      '', TRIM( U % MassDensityLabel )
    WRITE(*,'(A7,A7,3ES10.3E2,A1,A,A4,A7,3ES10.3E2,A1,A)') &
      '', 'V_L = ', V_L / U % VelocityUnitX1, &
      '', TRIM( U % VelocityLabelX1 ), &
      '', 'V_R = ', V_R / U % VelocityUnitX1, &
      '', TRIM( U % VelocityLabelX1 )
    WRITE(*,'(A7,A7,ES10.3E2,A1,A,A20,A7,ES10.3E2,A1,A)') &
      '', 'P_L = ', P_L / U % PressureUnit, &
      '', TRIM( U % PressureLabel ), &
      '', 'P_R = ', P_R / U % PressureUnit, &
      '', TRIM( U % PressureLabel )
    WRITE(*,'(A7,A7,ES10.3E2,A29,A7,ES10.3E2)') &
      '', 'Ye_L = ', Ye_L, '', 'Ye_R = ', Ye_R

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          DO iNodeX3 = 1, nNodesX(3)
            DO iNodeX2 = 1, nNodesX(2)
              DO iNodeX1 = 1, nNodesX(1)

                IF( Restart ) CYCLE

                X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

                iNode = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

                IF( X1 <= X_D )THEN

                  ! -- Left State --

                  uPF(iNode,iX1,iX2,iX3,iPF_D)  = D_L
                  uPF(iNode,iX1,iX2,iX3,iPF_V1) = V_L(1)
                  uPF(iNode,iX1,iX2,iX3,iPF_V2) = V_L(2)
                  uPF(iNode,iX1,iX2,iX3,iPF_V3) = V_L(3)
                  uAF(iNode,iX1,iX2,iX3,iAF_P)  = P_L
                  uAF(iNode,iX1,iX2,iX3,iAF_Ye) = Ye_L

                ELSE

                  ! -- Right State --

                  uPF(iNode,iX1,iX2,iX3,iPF_D)  = D_R
                  uPF(iNode,iX1,iX2,iX3,iPF_V1) = V_R(1)
                  uPF(iNode,iX1,iX2,iX3,iPF_V2) = V_R(2)
                  uPF(iNode,iX1,iX2,iX3,iPF_V3) = V_R(3)
                  uAF(iNode,iX1,iX2,iX3,iAF_P)  = P_R
                  uAF(iNode,iX1,iX2,iX3,iAF_Ye) = Ye_R

                END IF

              END DO
            END DO
          END DO

        END DO
      END DO
    END DO

    IF( Restart )THEN

      CALL ReadFluidFieldsRestart1D &
             ( RestartNumber, StartTime )

    END IF

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          IF( .NOT. Restart ) &
          CALL ComputeTemperatureFromPressure &
                 ( uPF(:,iX1,iX2,iX3,iPF_D ), uAF(:,iX1,iX2,iX3,iAF_P), &
                   uAF(:,iX1,iX2,iX3,iAF_Ye), uAF(:,iX1,iX2,iX3,iAF_T) )

          CALL ComputeThermodynamicStates_Primitive &
                 ( uPF(:,iX1,iX2,iX3,iPF_D ), uAF(:,iX1,iX2,iX3,iAF_T ), &
                   uAF(:,iX1,iX2,iX3,iAF_Ye), uPF(:,iX1,iX2,iX3,iPF_E ), &
                   uAF(:,iX1,iX2,iX3,iAF_E ), uPF(:,iX1,iX2,iX3,iPF_Ne) )

          CALL ApplyEquationOfState &
                 ( uPF(:,iX1,iX2,iX3,iPF_D ), uAF(:,iX1,iX2,iX3,iAF_T ), &
                   uAF(:,iX1,iX2,iX3,iAF_Ye), uAF(:,iX1,iX2,iX3,iAF_P ), &
                   uAF(:,iX1,iX2,iX3,iAF_S ), uAF(:,iX1,iX2,iX3,iAF_E ), &
                   uAF(:,iX1,iX2,iX3,iAF_Me), uAF(:,iX1,iX2,iX3,iAF_Mp), &
                   uAF(:,iX1,iX2,iX3,iAF_Mn), uAF(:,iX1,iX2,iX3,iAF_Xp), &
                   uAF(:,iX1,iX2,iX3,iAF_Xn), uAF(:,iX1,iX2,iX3,iAF_Xa), &
                   uAF(:,iX1,iX2,iX3,iAF_Xh), uAF(:,iX1,iX2,iX3,iAF_Gm) )

          CALL ComputeConserved &
                 ( uPF(:,iX1,iX2,iX3,1:nPF), uCF(:,iX1,iX2,iX3,1:nCF) )

        END DO
      END DO
    END DO

    END ASSOCIATE ! U

  END SUBROUTINE InitializeRiemannProblem1D


  SUBROUTINE InitializeShockEntropyWaveInteraction1D &
               ( StartTime, RestartNumber_Option )

    REAL(DP), INTENT(out)          :: StartTime
    INTEGER,  INTENT(in), OPTIONAL :: RestartNumber_Option

    LOGICAL  :: Restart
    INTEGER  :: iX1, iX2, iX3, RestartNumber
    INTEGER  :: iNodeX1, iNodeX2, iNodeX3, iNode
    REAL(DP) :: X1, Amplitude, Wavenumber

    Restart       = .FALSE.
    RestartNumber = 0
    IF( PRESENT( RestartNumber_Option ) )THEN
      Restart       = .TRUE.
      RestartNumber = RestartNumber_Option
    END IF

    StartTime = 0.0_DP

    Amplitude  = 2.0d11 * Gram / Centimeter**3
    Wavenumber = 2.5d00 / Kilometer

    WRITE(*,*)
    WRITE(*,'(A2,A6,A)') &
      '', 'INFO: ', TRIM( ProgramName )
    WRITE(*,*)

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          DO iNodeX3 = 1, nNodesX(3)
            DO iNodeX2 = 1, nNodesX(2)
              DO iNodeX1 = 1, nNodesX(1)

                IF( Restart ) CYCLE

                X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

                iNode = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

                IF( X1 <= - 8.0_DP * Kilometer )THEN

                  ! -- Left State --

                  uPF(iNode,iX1,iX2,iX3,iPF_D)  &
                    = 3.5740379540d12 * Gram / Centimeter**3
                  uPF(iNode,iX1,iX2,iX3,iPF_V1) &
                    = 8.3695270908d04 * Kilometer / Second
                  uPF(iNode,iX1,iX2,iX3,iPF_V2) &
                    = 0.0_DP
                  uPF(iNode,iX1,iX2,iX3,iPF_V3) &
                    = 0.0_DP
                  uAF(iNode,iX1,iX2,iX3,iAF_P)  &
                    = 1.0726263984d32 * Erg / Centimeter**3
                  uAF(iNode,iX1,iX2,iX3,iAF_Ye) &
                    = 0.3_DP

                ELSE

                  ! -- Right State --

                  uPF(iNode,iX1,iX2,iX3,iPF_D)  &
                    = 1.0d12 * Gram / Centimeter**3 &
                      + Amplitude * SIN( WaveNumber * X1 )
                  uPF(iNode,iX1,iX2,iX3,iPF_V1) = 0.0_DP
                  uPF(iNode,iX1,iX2,iX3,iPF_V2) = 0.0_DP
                  uPF(iNode,iX1,iX2,iX3,iPF_V3) = 0.0_DP
                  uAF(iNode,iX1,iX2,iX3,iAF_P)  = 1.0d31 * Erg / Centimeter**3
                  uAF(iNode,iX1,iX2,iX3,iAF_Ye) = 0.3_DP

                END IF

              END DO
            END DO
          END DO

        END DO
      END DO
    END DO

    IF( Restart )THEN

      CALL ReadFluidFieldsRestart1D &
             ( RestartNumber, StartTime )

    END IF

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          IF( .NOT. Restart ) &
          CALL ComputeTemperatureFromPressure &
                 ( uPF(:,iX1,iX2,iX3,iPF_D ), uAF(:,iX1,iX2,iX3,iAF_P), &
                   uAF(:,iX1,iX2,iX3,iAF_Ye), uAF(:,iX1,iX2,iX3,iAF_T) )

          CALL ComputeThermodynamicStates_Primitive &
                 ( uPF(:,iX1,iX2,iX3,iPF_D ), uAF(:,iX1,iX2,iX3,iAF_T ), &
                   uAF(:,iX1,iX2,iX3,iAF_Ye), uPF(:,iX1,iX2,iX3,iPF_E ), &
                   uAF(:,iX1,iX2,iX3,iAF_E ), uPF(:,iX1,iX2,iX3,iPF_Ne) )

          CALL ApplyEquationOfState &
                 ( uPF(:,iX1,iX2,iX3,iPF_D ), uAF(:,iX1,iX2,iX3,iAF_T ), &
                   uAF(:,iX1,iX2,iX3,iAF_Ye), uAF(:,iX1,iX2,iX3,iAF_P ), &
                   uAF(:,iX1,iX2,iX3,iAF_S ), uAF(:,iX1,iX2,iX3,iAF_E ), &
                   uAF(:,iX1,iX2,iX3,iAF_Me), uAF(:,iX1,iX2,iX3,iAF_Mp), &
                   uAF(:,iX1,iX2,iX3,iAF_Mn), uAF(:,iX1,iX2,iX3,iAF_Xp), &
                   uAF(:,iX1,iX2,iX3,iAF_Xn), uAF(:,iX1,iX2,iX3,iAF_Xa), &
                   uAF(:,iX1,iX2,iX3,iAF_Xh), uAF(:,iX1,iX2,iX3,iAF_Gm) )

          CALL ComputeConserved &
                 ( uPF(:,iX1,iX2,iX3,1:nPF), uCF(:,iX1,iX2,iX3,1:nCF) )

        END DO
      END DO
    END DO

  END SUBROUTINE InitializeShockEntropyWaveInteraction1D


END MODULE RiemannProblemInitializationModule
