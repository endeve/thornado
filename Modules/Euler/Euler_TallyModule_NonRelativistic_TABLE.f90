MODULE Euler_TallyModule_NonRelativistic_TABLE

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
    iGF_SqrtGm, &
    iGF_Phi_N
  USE FluidFieldsModule, ONLY: &
    nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne
  USE Euler_UtilitiesModule_NonRelativistic, ONLY: &
    Euler_ComputePrimitive_NonRelativistic

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Euler_InitializeTally_NonRelativistic_TABLE
  PUBLIC :: Euler_FinalizeTally_NonRelativistic_TABLE
  PUBLIC :: Euler_ComputeTally_NonRelativistic_TABLE

  CHARACTER(256)        :: TallyFileName
  INTEGER               :: nTallies
  INTEGER               :: iTally_E_i
  INTEGER               :: iTally_E_k
  INTEGER               :: iTally_E_g
  REAL(DP), ALLOCATABLE :: EulerTally(:,:)

CONTAINS


  SUBROUTINE Euler_InitializeTally_NonRelativistic_TABLE( iX_B0, iX_E0, G, U )

    INTEGER,  INTENT(in) :: &
      iX_B0(3), iX_E0(3)
    REAL(DP), INTENT(in) :: &
      G(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)
    REAL(DP), INTENT(in) :: &
      U(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)

    INTEGER :: FileUnit

    nTallies   = nCF + 3
    iTally_E_i = nCF + 1
    iTally_E_k = nCF + 2
    iTally_E_g = nCF + 3

    ALLOCATE( EulerTally(1:nTallies,0:1) )

    TallyFileName = '../Output/EulerTally.dat'

    OPEN( NEWUNIT = FileUnit, FILE = TRIM( TallyFileName ) )

    WRITE( FileUnit, '(9(A20,x))' ) &
      'Time', 'D', 'S1', 'S2', 'S3', 'E', 'E_i', 'E_k', 'E_g'

    CLOSE( FileUnit )

    CALL Euler_ComputeTally_NonRelativistic_TABLE &
           ( iX_B0, iX_E0, G, U, Time = Zero, &
             iState_Option = 0, DisplayTally_Option = .FALSE. )

    CALL Euler_ComputeTally_NonRelativistic_TABLE &
           ( iX_B0, iX_E0, G, U, Time = Zero, &
             iState_Option = 1, DisplayTally_Option = .TRUE. )

  END SUBROUTINE Euler_InitializeTally_NonRelativistic_TABLE


  SUBROUTINE Euler_FinalizeTally_NonRelativistic_TABLE

    DEALLOCATE( EulerTally )

  END SUBROUTINE Euler_FinalizeTally_NonRelativistic_TABLE


  SUBROUTINE Euler_ComputeTally_NonRelativistic_TABLE &
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
    REAL(DP) :: P(nDOFX,nPF)

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

          CALL Euler_ComputePrimitive_NonRelativistic &
                 ( U(:,iX1,iX2,iX3,iCF_D) , U(:,iX1,iX2,iX3,iCF_S1), &
                   U(:,iX1,iX2,iX3,iCF_S2), U(:,iX1,iX2,iX3,iCF_S3), &
                   U(:,iX1,iX2,iX3,iCF_E) , U(:,iX1,iX2,iX3,iCF_Ne), &
                   P(:,iPF_D) , P(:,iPF_V1), P(:,iPF_V2), &
                   P(:,iPF_V3), P(:,iPF_E) , P(:,iPF_Ne), &
                   G(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
                   G(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
                   G(:,iX1,iX2,iX3,iGF_Gm_dd_33) )

          ! --- Internal Energy ---

          EulerTally(iTally_E_i,iState) &
            = EulerTally(iTally_E_i,iState) &
                + SUM( WeightsX_q(:) &
                         * G(:,iX1,iX2,iX3,iGF_SqrtGm) &
                         * P(:,iPF_E) ) &
                    * dX1(iX1) * dX2(iX2) * dX3(iX3)

          ! --- Kinetic Energy ---

          EulerTally(iTally_E_k,iState) &
            = EulerTally(iTally_E_k,iState) &
                + SUM( WeightsX_q(:) &
                         * G(:,iX1,iX2,iX3,iGF_SqrtGm) &
                         * 0.5_DP &
                         * (   U(:,iX1,iX2,iX3,iCF_S1) * P(:,iPF_V1) &
                             + U(:,iX1,iX2,iX3,iCF_S2) * P(:,iPF_V2) &
                             + U(:,iX1,iX2,iX3,iCF_S3) * P(:,iPF_V3) ) ) &
                    * dX1(iX1) * dX2(iX2) * dX3(iX3)

          ! --- Gravitational Energy ---

          EulerTally(iTally_E_g,iState) &
            = EulerTally(iTally_E_g,iState) &
                + SUM( WeightsX_q(:) &
                         * G(:,iX1,iX2,iX3,iGF_SqrtGm) &
                         * 0.5_DP &
                         * U(:,iX1,iX2,iX3,iCF_D) &
                         * G(:,iX1,iX2,iX3,iGF_Phi_N) ) &
                    * dX1(iX1) * dX2(iX2) * dX3(iX3)

        END DO
      END DO
    END DO

    END ASSOCIATE ! dX1, etc.

    IF( DisplayTally )THEN

      CALL WriteTally_Euler( Time )

    END IF

  END SUBROUTINE Euler_ComputeTally_NonRelativistic_TABLE


  SUBROUTINE WriteTally_Euler( Time )

    REAL(DP), INTENT(in) :: Time

    INTEGER :: FileUnit

    ASSOCIATE( U => UnitsDisplay )

    OPEN( NEWUNIT=FileUnit, FILE=TRIM( TallyFileName ), POSITION='APPEND', ACTION='WRITE' )

    WRITE( FileUnit, '(9(ES20.12,x))' ) &
      Time / U % TimeUnit, &
      EulerTally(iCF_D,     1) / U % MassUnit, &
      EulerTally(iCF_S1,    1) / U % MomentumUnit, &
      EulerTally(iCF_S2,    1) / U % MomentumUnit, &
      EulerTally(iCF_S3,    1) / U % MomentumUnit, &
      EulerTally(iCF_E,     1) / U % EnergyGlobalUnit, &
      EulerTally(iTally_E_i,1) / U % EnergyGlobalUnit, &
      EulerTally(iTally_E_k,1) / U % EnergyGlobalUnit, &
      EulerTally(iTally_E_g,1) / U % EnergyGlobalUnit

    CLOSE( FileUnit )

    END ASSOCIATE ! U

  END SUBROUTINE WriteTally_Euler


END MODULE Euler_TallyModule_NonRelativistic_TABLE
