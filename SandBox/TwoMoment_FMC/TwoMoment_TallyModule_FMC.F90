MODULE TwoMoment_TallyModule_FMC

  USE KindModule, ONLY: &
    DP, Zero, FourPi
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    nDOFE, &
    nDOFX, &
    nDOFZ
  USE ReferenceElementModule, ONLY: &
    Weights_q
  USE MeshModule, ONLY: &
    MeshE, &
    MeshX
  USE GeometryFieldsModuleE, ONLY: &
    nGE, iGE_Ep2
  USE GeometryFieldsModule, ONLY: &
    nGF, iGF_SqrtGm
  USE TwoMoment_FieldsModule_FMC, ONLY: &
    nCM, iCM_E, iCM_F1, iCM_F2, iCM_F3, &
    nSpecies

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeTally
  PUBLIC :: ComputeTally
  PUBLIC :: IncrementOffGridTally_TwoMoment
  PUBLIC :: FinalizeTally

  CHARACTER(256) :: EulerianEnergy_FileName
  REAL(DP)       :: EulerianEnergy_Interior
  REAL(DP)       :: EulerianEnergy_Initial
  REAL(DP)       :: EulerianEnergy_OffGrid
  REAL(DP)       :: EulerianEnergy_Change

  CHARACTER(256) :: EulerianMomentum1_FileName
  REAL(DP)       :: EulerianMomentum1_Interior
  REAL(DP)       :: EulerianMomentum1_Initial
  REAL(DP)       :: EulerianMomentum1_OffGrid
  REAL(DP)       :: EulerianMomentum1_Change

  CHARACTER(256) :: EulerianMomentum2_FileName
  REAL(DP)       :: EulerianMomentum2_Interior
  REAL(DP)       :: EulerianMomentum2_Initial
  REAL(DP)       :: EulerianMomentum2_OffGrid
  REAL(DP)       :: EulerianMomentum2_Change

  CHARACTER(256) :: EulerianMomentum3_FileName
  REAL(DP)       :: EulerianMomentum3_Interior
  REAL(DP)       :: EulerianMomentum3_Initial
  REAL(DP)       :: EulerianMomentum3_OffGrid
  REAL(DP)       :: EulerianMomentum3_Change

CONTAINS


  SUBROUTINE InitializeTally

    CHARACTER(256) :: BaseFileName

    BaseFileName = '../Output/' // TRIM( ProgramName )

    ! --- Eulerian Energy ---

    EulerianEnergy_FileName &
      = TRIM( BaseFileName ) // '_Tally_EulerianEnergy.dat'

    CALL WriteTally_Header( EulerianEnergy_FileName )

    EulerianEnergy_Interior = Zero
    EulerianEnergy_Initial  = Zero
    EulerianEnergy_OffGrid  = Zero
    EulerianEnergy_Change   = Zero

    ! --- Eulerian Momentum 1 ---

    EulerianMomentum1_FileName &
      = TRIM( BaseFileName ) // '_Tally_EulerianMomentum1.dat'

    CALL WriteTally_Header( EulerianMomentum1_FileName )

    EulerianMomentum1_Interior = Zero
    EulerianMomentum1_Initial  = Zero
    EulerianMomentum1_OffGrid  = Zero
    EulerianMomentum1_Change   = Zero

    ! --- Eulerian Momentum 2 ---

    EulerianMomentum2_FileName &
      = TRIM( BaseFileName ) // '_Tally_EulerianMomentum2.dat'

    CALL WriteTally_Header( EulerianMomentum2_FileName )

    EulerianMomentum2_Interior = Zero
    EulerianMomentum2_Initial  = Zero
    EulerianMomentum2_OffGrid  = Zero
    EulerianMomentum2_Change   = Zero

    ! --- Eulerian Momentum 3 ---

    EulerianMomentum3_FileName &
      = TRIM( BaseFileName ) // '_Tally_EulerianMomentum3.dat'

    CALL WriteTally_Header( EulerianMomentum3_FileName )

    EulerianMomentum3_Interior = Zero
    EulerianMomentum3_Initial  = Zero
    EulerianMomentum3_OffGrid  = Zero
    EulerianMomentum3_Change   = Zero

  END SUBROUTINE InitializeTally


  SUBROUTINE ComputeTally &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, Time, GE, GX, M, &
      SetInitialValues_Option, Verbose_Option )

    INTEGER,  INTENT(in) :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in) :: &
      Time
    REAL(DP), INTENT(in) :: &
      GE(1:nDOFE, &
         iZ_B1(1):iZ_E1(1), &
         1:nGE)
    REAL(DP), INTENT(in) :: &
      GX(1:nDOFX, &
         iZ_B1(2):iZ_E1(2), &
         iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4), &
         1:nGF)
    REAL(DP), INTENT(in) :: &
      M (1:nDOFZ, &
         iZ_B1(1):iZ_E1(1), &
         iZ_B1(2):iZ_E1(2), &
         iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4), &
         1:nCM,1:nSpecies)
    LOGICAL, INTENT(in), OPTIONAL :: &
      SetInitialValues_Option
    LOGICAL, INTENT(in), OPTIONAL :: &
      Verbose_Option

    LOGICAL :: SetInitialValues
    LOGICAL :: Verbose

    SetInitialValues = .FALSE.
    IF( PRESENT( SetInitialValues_Option ) )THEN
      SetInitialValues = SetInitialValues_Option
    END IF

    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) )THEN
      Verbose = Verbose_Option
    END IF

    CALL ComputeTally_TwoMoment &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, M )

    IF( SetInitialValues )THEN

      EulerianEnergy_Initial    = EulerianEnergy_Interior
      EulerianMomentum1_Initial = EulerianMomentum1_Interior
      EulerianMomentum2_Initial = EulerianMomentum2_Interior
      EulerianMomentum3_Initial = EulerianMomentum3_Interior

    END IF

    EulerianEnergy_Change &
      = EulerianEnergy_Interior &
          - ( EulerianEnergy_Initial    + EulerianEnergy_OffGrid    )

    EulerianMomentum1_Change &
      = EulerianMomentum1_Interior &
          - ( EulerianMomentum1_Initial + EulerianMomentum1_OffGrid )

    EulerianMomentum2_Change &
      = EulerianMomentum2_Interior &
          - ( EulerianMomentum2_Initial + EulerianMomentum2_OffGrid )

    EulerianMomentum3_Change &
      = EulerianMomentum3_Interior &
          - ( EulerianMomentum3_Initial + EulerianMomentum3_OffGrid )

    CALL WriteTally( Time )

    IF( Verbose )THEN

      CALL DisplayTally( Time )

    END IF

  END SUBROUTINE ComputeTally


  SUBROUTINE ComputeTally_TwoMoment( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, M )

    INTEGER,  INTENT(in) :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in) :: &
      GE(1:nDOFE, &
         iZ_B1(1):iZ_E1(1), &
         1:nGE)
    REAL(DP), INTENT(in) :: &
      GX(1:nDOFX, &
         iZ_B1(2):iZ_E1(2), &
         iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4), &
         1:nGF)
    REAL(DP), INTENT(in) :: &
      M (1:nDOFZ, &
         iZ_B1(1):iZ_E1(1), &
         iZ_B1(2):iZ_E1(2), &
         iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4), &
         1:nCM,1:nSpecies)

    INTEGER :: iZ1, iZ2, iZ3, iZ4, iS
    INTEGER :: iNodeE, iNodeX, iNodeZ

    ASSOCIATE &
      ( dZ1 => MeshE    % Width, dZ2 => MeshX(1) % Width, &
        dZ3 => MeshX(2) % Width, dZ4 => MeshX(3) % Width )

    EulerianEnergy_Interior    = Zero
    EulerianMomentum1_Interior = Zero
    EulerianMomentum2_Interior = Zero
    EulerianMomentum3_Interior = Zero

    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeX = 1, nDOFX
      DO iNodeE = 1, nDOFE

        iNodeZ = (iNodeX-1) * nDOFE + iNodeE

        EulerianEnergy_Interior                           &
          = EulerianEnergy_Interior                       &
              + dZ1(iZ1) * dZ2(iZ2) * dZ3(iZ3) * dZ4(iZ4) &
                * Weights_q(iNodeZ)                       &
                * GE(iNodeE,iZ1,iGE_Ep2)                  &
                * GX(iNodeX,iZ2,iZ3,iZ4,iGF_SqrtGm)       &
                * M(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_E,iS)

        EulerianMomentum1_Interior                        &
          = EulerianMomentum1_Interior                    &
              + dZ1(iZ1) * dZ2(iZ2) * dZ3(iZ3) * dZ4(iZ4) &
                * Weights_q(iNodeZ)                       &
                * GE(iNodeE,iZ1,iGE_Ep2)                  &
                * GX(iNodeX,iZ2,iZ3,iZ4,iGF_SqrtGm)       &
                * M(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_F1,iS)

        EulerianMomentum2_Interior                        &
          = EulerianMomentum2_Interior                    &
              + dZ1(iZ1) * dZ2(iZ2) * dZ3(iZ3) * dZ4(iZ4) &
                * Weights_q(iNodeZ)                       &
                * GE(iNodeE,iZ1,iGE_Ep2)                  &
                * GX(iNodeX,iZ2,iZ3,iZ4,iGF_SqrtGm)       &
                * M(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_F2,iS)

        EulerianMomentum3_Interior                        &
          = EulerianMomentum3_Interior                    &
              + dZ1(iZ1) * dZ2(iZ2) * dZ3(iZ3) * dZ4(iZ4) &
                * Weights_q(iNodeZ)                       &
                * GE(iNodeE,iZ1,iGE_Ep2)                  &
                * GX(iNodeX,iZ2,iZ3,iZ4,iGF_SqrtGm)       &
                * M(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_F3,iS)

      END DO
      END DO

    END DO
    END DO
    END DO
    END DO
    END DO

    EulerianEnergy_Interior    = FourPi * EulerianEnergy_Interior
    EulerianMomentum1_Interior = FourPi * EulerianMomentum1_Interior
    EulerianMomentum2_Interior = FourPi * EulerianMomentum2_Interior
    EulerianMomentum3_Interior = FourPi * EulerianMomentum3_Interior

    END ASSOCIATE ! dZ1, etc.

  END SUBROUTINE ComputeTally_TwoMoment


  SUBROUTINE IncrementOffGridTally_TwoMoment( dM )

    REAL(DP), INTENT(in) :: dM(nCM)

    EulerianEnergy_OffGrid    = EulerianEnergy_OffGrid    + FourPi * dM(iCM_E )
    EulerianMomentum1_OffGrid = EulerianMomentum1_OffGrid + FourPi * dM(iCM_F1)
    EulerianMomentum2_OffGrid = EulerianMomentum2_OffGrid + FourPi * dM(iCM_F2)
    EulerianMomentum3_OffGrid = EulerianMomentum3_OffGrid + FourPi * dM(iCM_F3)

  END SUBROUTINE IncrementOffGridTally_TwoMoment


  SUBROUTINE WriteTally( Time )

    REAL(DP), INTENT(in) :: Time

    CALL WriteTally_Variable &
           ( Time, &
             EulerianEnergy_Interior, &
             EulerianEnergy_OffGrid , &
             EulerianEnergy_Initial , &
             EulerianEnergy_Change  , &
             EulerianEnergy_FileName )

    CALL WriteTally_Variable &
           ( Time, &
             EulerianMomentum1_Interior, &
             EulerianMomentum1_OffGrid , &
             EulerianMomentum1_Initial , &
             EulerianMomentum1_Change  , &
             EulerianMomentum1_FileName )

    CALL WriteTally_Variable &
           ( Time, &
             EulerianMomentum2_Interior, &
             EulerianMomentum2_OffGrid , &
             EulerianMomentum2_Initial , &
             EulerianMomentum2_Change  , &
             EulerianMomentum2_FileName )

    CALL WriteTally_Variable &
           ( Time, &
             EulerianMomentum3_Interior, &
             EulerianMomentum3_OffGrid , &
             EulerianMomentum3_Initial , &
             EulerianMomentum3_Change  , &
             EulerianMomentum3_FileName )

  END SUBROUTINE WriteTally


  SUBROUTINE WriteTally_Header( FileName )

    CHARACTER(*), INTENT(in) :: FileName

    INTEGER :: FileUnit

    OPEN( NEWUNIT = FileUnit, FILE = TRIM( FileName ) )

    WRITE( FileUnit, '(5(A20,x))' ) &
      'Time', 'Interior', 'Off Grid', 'Initial', 'Change'

    CLOSE( FileUnit )

  END SUBROUTINE WriteTally_Header


  SUBROUTINE WriteTally_Variable &
    ( Time, Interior, OffGrid, Initial, Change, FileName )

    REAL(DP)    , INTENT(in) :: Time, Interior, OffGrid, Initial, Change
    CHARACTER(*), INTENT(in) :: FileName

    INTEGER :: FileUnit

    OPEN( NEWUNIT = FileUnit, FILE = TRIM( FileName ), &
          POSITION = 'APPEND', ACTION = 'WRITE' )

    WRITE( FileUnit, '(5(ES20.12,x))' ) &
      Time, Interior, OffGrid , Initial , Change

    CLOSE( FileUnit )

  END SUBROUTINE WriteTally_Variable


  SUBROUTINE DisplayTally( Time )

    REAL(DP), INTENT(in) :: Time

    WRITE(*,*)
    WRITE(*,'(A8,A,ES8.2E2)') '', 'Two-Moment Tally FMC. t = ', &
      Time
    WRITE(*,*)
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Eulerian Energy Interior........: ', &
      EulerianEnergy_Interior
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Eulerian Energy Initial.........: ', &
      EulerianEnergy_Initial
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Eulerian Energy Off Grid........: ', &
      EulerianEnergy_OffGrid
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Eulerian Energy Change..........: ', &
      EulerianEnergy_Change
    WRITE(*,*)
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Eulerian Momentum 1 Interior....: ', &
      EulerianMomentum1_Interior
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Eulerian Momentum 1 Initial.....: ', &
      EulerianMomentum1_Initial
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Eulerian Momentum 1 Off Grid....: ', &
      EulerianMomentum1_OffGrid
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Eulerian Momentum 1 Change......: ', &
      EulerianMomentum1_Change
    WRITE(*,*)
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Eulerian Momentum 2 Interior....: ', &
      EulerianMomentum2_Interior
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Eulerian Momentum 2 Initial.....: ', &
      EulerianMomentum2_Initial
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Eulerian Momentum 2 Off Grid....: ', &
      EulerianMomentum2_OffGrid
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Eulerian Momentum 2 Change......: ', &
      EulerianMomentum2_Change
    WRITE(*,*)
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Eulerian Momentum 3 Interior....: ', &
      EulerianMomentum3_Interior
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Eulerian Momentum 3 Initial.....: ', &
      EulerianMomentum3_Initial
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Eulerian Momentum 3 Off Grid....: ', &
      EulerianMomentum3_OffGrid
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Eulerian Momentum 3 Change......: ', &
      EulerianMomentum3_Change
    WRITE(*,*)

  END SUBROUTINE DisplayTally


  SUBROUTINE FinalizeTally

  END SUBROUTINE FinalizeTally


END MODULE TwoMoment_TallyModule_FMC
