MODULE Euler_TallyModule_Relativistic

  USE KindModule, ONLY: &
    DP, &
    Zero, &
    Two,  &
    Pi
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
    iGF_Psi
  USE FluidFieldsModule, ONLY: &
    nCF, &
    iCF_D, &
    iCF_S1, &
    iCF_S2, &
    iCF_S3, &
    iCF_E,  &
    iCF_Ne, &
    nPF, &
    iPF_D, &
    iPF_V1, &
    iPF_V2, &
    iPF_V3, &
    iPF_E, &
    iPF_Ne
  USE Euler_UtilitiesModule_Relativistic, ONLY: &
    ComputePrimitive_Euler_Relativistic

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeTally_Euler_Relativistic
  PUBLIC :: FinalizeTally_Euler_Relativistic
  PUBLIC :: ComputeTally_Euler_Relativistic

  LOGICAL               :: SuppressTally
  CHARACTER(256)        :: TallyFileName, FMT
  INTEGER               :: nTallies
  INTEGER               :: iTally_Eb ! baryonic mass
  INTEGER               :: iTally_Eg ! gravitational mass
  REAL(DP), ALLOCATABLE :: EulerTally(:)


CONTAINS


  SUBROUTINE InitializeTally_Euler_Relativistic &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, SuppressTally_Option )

    INTEGER,  INTENT(in)           :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)           :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(in)           :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    LOGICAL,  INTENT(in), OPTIONAL :: &
      SuppressTally_Option

    INTEGER :: FileUnit

    SuppressTally = .FALSE.
    IF( PRESENT( SuppressTally_Option ) ) &
      SuppressTally = SuppressTally_Option

    IF( SuppressTally ) RETURN

    iTally_Eb = nCF + 1
    iTally_Eg = nCF + 2
    nTallies  = nCF + 2

    ALLOCATE( EulerTally(1:nTallies) )

    WRITE(FMT,'(A,I2.2,A)') '(',  nTallies, 'ES25.16E3,1x)'

    TallyFileName = '../Output/EulerTally.dat'

    OPEN( NEWUNIT = FileUnit, FILE = TRIM( TallyFileName ) )

    WRITE(FileUnit,'(A)') &
      'Time, D, S1, S2, S3, tau, BaryonicMass, GravitationalMass'

    CLOSE( FileUnit )

  END SUBROUTINE InitializeTally_Euler_Relativistic


  SUBROUTINE FinalizeTally_Euler_Relativistic

    IF( .NOT. SuppressTally ) &
      DEALLOCATE( EulerTally )

  END SUBROUTINE FinalizeTally_Euler_Relativistic


  SUBROUTINE ComputeTally_Euler_Relativistic &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, Time )

    INTEGER,  INTENT(in) :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in) :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(in) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(in) :: &
      Time

    INTEGER  :: iCF, iX1, iX2, iX3, iErr(nDOFX)
    REAL(DP) :: P(nDOFX,nPF), d3X

    IF( SuppressTally ) RETURN

    ASSOCIATE &
      ( dX1 => MeshX(1) % Width(iX_B0(1):iX_E0(1)), &
        dX2 => MeshX(2) % Width(iX_B0(2):iX_E0(2)), &
        dX3 => MeshX(3) % Width(iX_B0(3):iX_E0(3)) )

    EulerTally = Zero

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      d3X = Two / Pi * dX1(iX1) * dX2(iX2) * dX3(iX3) ! Hack for 1D spherical

      DO iCF = 1, nCF

        EulerTally(iCF) &
          = EulerTally(iCF) &
              + SUM( WeightsX_q(:) &
                       * G(:,iX1,iX2,iX3,iGF_SqrtGm) &
                       * U(:,iX1,iX2,iX3,iCF) ) &
                  * d3X
      END DO

      CALL ComputePrimitive_Euler_Relativistic &
             ( U(:,iX1,iX2,iX3,iCF_D) , U(:,iX1,iX2,iX3,iCF_S1), &
               U(:,iX1,iX2,iX3,iCF_S2), U(:,iX1,iX2,iX3,iCF_S3), &
               U(:,iX1,iX2,iX3,iCF_E) , U(:,iX1,iX2,iX3,iCF_Ne), &
               P(:,iPF_D) , P(:,iPF_V1), P(:,iPF_V2), &
               P(:,iPF_V3), P(:,iPF_E) , P(:,iPF_Ne), &
               G(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               G(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               G(:,iX1,iX2,iX3,iGF_Gm_dd_33), &
               iErr )

      ! --- Baryonic Mass ---

      EulerTally(iTally_Eb) &
        = EulerTally(iTally_Eb) &
            + d3X * SUM( WeightsX_q(:) &
                     * U(:,iX1,iX2,iX3,iCF_D) &
                     * G(:,iX1,iX2,iX3,iGF_SqrtGm) )

      ! --- Gravitational Mass ---

      EulerTally(iTally_Eg) &
        = EulerTally(iTally_Eg) &
            + d3X * SUM( WeightsX_q(:) &
                     * ( P(:,iPF_D) + P(:,iPF_E) ) &
                     * G(:,iX1,iX2,iX3,iGF_SqrtGm) / G(:,iX1,iX2,iX3,iGF_Psi) )

    END DO
    END DO
    END DO

    END ASSOCIATE ! dX1, etc.

    CALL WriteTally_Euler( Time )

  END SUBROUTINE ComputeTally_Euler_Relativistic


  SUBROUTINE WriteTally_Euler( Time )

    REAL(DP), INTENT(in) :: Time

    INTEGER :: FileUnit

    ASSOCIATE( U => UnitsDisplay )

    OPEN( NEWUNIT = FileUnit, FILE = TRIM( TallyFileName ), &
          POSITION = 'APPEND', ACTION = 'WRITE' )

    WRITE(FileUnit,FMT) &
      Time                  / U % TimeUnit,         &
      EulerTally(iCF_D    ) / U % MassUnit,         &
      EulerTally(iCF_S1   ) / U % MomentumUnit,     &
      EulerTally(iCF_S2   ) / U % MomentumUnit,     &
      EulerTally(iCF_S3   ) / U % MomentumUnit,     &
      EulerTally(iCF_E    ) / U % EnergyGlobalUnit, &
      EulerTally(iTally_Eb) / U % MassUnit,         &
      EulerTally(iTally_Eg) / U % MassUnit

    CLOSE( FileUnit )

    END ASSOCIATE ! U

  END SUBROUTINE WriteTally_Euler


END MODULE Euler_TallyModule_Relativistic
