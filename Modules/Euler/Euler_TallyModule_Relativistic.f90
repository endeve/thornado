MODULE Euler_TallyModule_Relativistic

  USE KindModule, ONLY: &
    DP, Zero
  USE UnitsModule, ONLY: &
    UnitsDisplay
  USE ProgramHeaderModule, ONLY: &
    nDOFX
  USE ReferenceElementModuleX, ONLY: &
    WeightsX_q
  USE MeshModule, ONLY: &
    MeshX
  USE GeometryFieldsModule, ONLY: &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_SqrtGm
  USE FluidFieldsModule, ONLY: &
    nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne
  USE Euler_UtilitiesModule_Relativistic, ONLY: &
    Euler_ComputePrimitive_Relativistic
  USE EquationOfStateModule, ONLY: &
    ComputePressureFromPrimitive_Relativistic

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Euler_InitializeTally_Relativistic
  PUBLIC :: Euler_FinalizeTally_Relativistic
  PUBLIC :: Euler_ComputeTally_Relativistic

  CHARACTER(256)        :: TallyFileName
  INTEGER               :: nTallies
  INTEGER               :: iTally_Pres
  REAL(DP), ALLOCATABLE :: EulerTally(:,:)


CONTAINS


  SUBROUTINE Euler_InitializeTally_Relativistic( iX_B0, iX_E0, G, U )

    INTEGER,  INTENT(in) :: &
      iX_B0(3), iX_E0(3)
    REAL(DP), INTENT(in) :: &
      G(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)
    REAL(DP), INTENT(in) :: &
      U(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)

    INTEGER :: FileUnit

    nTallies    = nCF + 1
    iTally_Pres = nCF + 1

    ALLOCATE( EulerTally(1:nTallies,0:1) )

    TallyFileName = '../Output/.EulerTally.dat'

    OPEN( NEWUNIT = FileUnit, FILE = TRIM( TallyFileName ) )

    WRITE( FileUnit, '(7(A20,x))' ) &
      'Time', 'D', 'S1', 'S2', 'S3', 'tau', 'P'

    CLOSE( FileUnit )

    CALL Euler_ComputeTally_Relativistic &
           ( iX_B0, iX_E0, G, U, Time = Zero, &
             iState_Option = 0, DisplayTally_Option = .FALSE. )

    CALL Euler_ComputeTally_Relativistic &
           ( iX_B0, iX_E0, G, U, Time = Zero, &
             iState_Option = 1, DisplayTally_Option = .TRUE. )

  END SUBROUTINE Euler_InitializeTally_Relativistic


  SUBROUTINE Euler_FinalizeTally_Relativistic

    DEALLOCATE( EulerTally )

  END SUBROUTINE Euler_FinalizeTally_Relativistic


  SUBROUTINE Euler_ComputeTally_Relativistic &
    ( iX_B0, iX_E0, G, U, Time, iState_Option, DisplayTally_Option )

    INTEGER,  INTENT(in) :: &
      iX_B0(3), iX_E0(3)
    REAL(DP), INTENT(in) :: &
      G(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)
    REAL(DP), INTENT(in) :: &
      U(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)
    REAL(DP), INTENT(in) :: &
      Time
    INTEGER,  INTENT(in), OPTIONAL :: &
      iState_Option
    LOGICAL,  INTENT(in), OPTIONAL :: &
      DisplayTally_Option

    LOGICAL  :: DisplayTally
    INTEGER  :: iState, iCF
    INTEGER  :: iX1, iX2, iX3
    REAL(DP) :: P(nDOFX,nPF), Pressure(nDOFX)

    IF( PRESENT( iState_Option ) )THEN
      iState = iState_Option
    ELSE
      iState = 1
    END IF

    IF( PRESENT( DisplayTally_Option ) )THEN
      DisplayTally = DisplayTally_Option
    ELSE
      DisplayTally = .FALSE.
    END IF

    ASSOCIATE &
      ( dX1 => MeshX(1) % Width(iX_B0(1):iX_E0(1)), &
        dX2 => MeshX(2) % Width(iX_B0(2):iX_E0(2)), &
        dX3 => MeshX(3) % Width(iX_B0(3):iX_E0(3)) )

    EulerTally(:,iState) = Zero

    ! --- Conserved Fluid ---

    DO iCF = 1, nCF

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)

        EulerTally(iCF,iState) &
          = EulerTally(iCF,iState) &
              + SUM( WeightsX_q(:) &
                       * G(:,iX1,iX2,iX3,iGF_SqrtGm) &
                       * U(:,iX1,iX2,iX3,iCF) ) &
                  * dX1(iX1) * dX2(iX2) * dX3(iX3)

      END DO
      END DO
      END DO

    END DO

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      CALL Euler_ComputePrimitive_Relativistic &
             ( U(:,iX1,iX2,iX3,iCF_D) , U(:,iX1,iX2,iX3,iCF_S1), &
               U(:,iX1,iX2,iX3,iCF_S2), U(:,iX1,iX2,iX3,iCF_S3), &
               U(:,iX1,iX2,iX3,iCF_E) , U(:,iX1,iX2,iX3,iCF_Ne), &
               P(:,iPF_D) , P(:,iPF_V1), P(:,iPF_V2), &
               P(:,iPF_V3), P(:,iPF_E) , P(:,iPF_Ne), &
               G(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               G(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               G(:,iX1,iX2,iX3,iGF_Gm_dd_33) )

      CALL ComputePressureFromPrimitive_Relativistic &
             ( P(:,iPF_D), P(:,iPF_E), P(:,iPF_Ne), Pressure )

      ! --- Pressure ---

      EulerTally(iTally_Pres,iState) &
        = EulerTally(iTally_Pres,iState) &
            + SUM( WeightsX_q(:) &
                     * G(:,iX1,iX2,iX3,iGF_SqrtGm) &
                     * Pressure(:) ) &
                * dX1(iX1) * dX2(iX2) * dX3(iX3)

    END DO
    END DO
    END DO

    END ASSOCIATE ! dX1, etc.

    IF( DisplayTally )THEN

      CALL WriteTally_Euler( Time )

    END IF

  END SUBROUTINE Euler_ComputeTally_Relativistic


  SUBROUTINE WriteTally_Euler( Time )

    REAL(DP), INTENT(in) :: Time

    INTEGER :: FileUnit

    ASSOCIATE( U => UnitsDisplay )

    OPEN( NEWUNIT=FileUnit, FILE=TRIM( TallyFileName ), ACCESS='APPEND' )

    WRITE( FileUnit, '(7(ES20.12,x))' ) &
      Time                      / U % TimeUnit, &
      EulerTally(iCF_D,      1) / U % MassUnit, &
      EulerTally(iCF_S1,     1) / U % MomentumUnit, &
      EulerTally(iCF_S2,     1) / U % MomentumUnit, &
      EulerTally(iCF_S3,     1) / U % MomentumUnit, &
      EulerTally(iCF_E,      1) / U % EnergyGlobalUnit, &
      EulerTally(iTally_Pres,1) / U % PressureUnit

    CLOSE( FileUnit )

    END ASSOCIATE ! U

  END SUBROUTINE WriteTally_Euler


END MODULE Euler_TallyModule_Relativistic
