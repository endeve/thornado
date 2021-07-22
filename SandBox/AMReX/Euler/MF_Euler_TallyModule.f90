MODULE MF_Euler_TallyModule

  ! --- AMReX Modules ---

  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_geometry_module, ONLY: &
    amrex_geometry
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor, &
    amrex_parallel_reduce_sum

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    nNodesX, &
    nDimsX
  USE ReferenceElementModuleX, ONLY: &
    WeightsX_q
  USE UnitsModule, ONLY: &
    UnitsDisplay
  USE MeshModule, ONLY: &
    MeshType, &
    CreateMesh, &
    DestroyMesh
  USE GeometryFieldsModule, ONLY: &
    nGF, &
    iGF_SqrtGm, &
    CoordinateSystem
  USE FluidFieldsModule, ONLY: &
    nCF, &
    iCF_D, &
    iCF_E

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero, &
    FourPi
  USE InputParsingModule, ONLY: &
    nX, &
    nLevels, &
    ProgramName, &
    UseTiling, &
    xL, &
    xR
  USE TimersModule_AMReX_Euler, ONLY: &
    TimersStart_AMReX_Euler, &
    TimersStop_AMReX_Euler, &
    Timer_AMReX_Euler_Allocate
  USE MF_UtilitiesModule, ONLY: &
    amrex2thornado_X

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MF_InitializeTally_Euler
  PUBLIC :: MF_ComputeTally_Euler
  PUBLIC :: MF_IncrementOffGridTally_Euler
  PUBLIC :: MF_FinalizeTally_Euler

  LOGICAL :: SuppressTally

  CHARACTER(256) :: BaryonicMass_FileName
  REAL(DP), ALLOCATABLE :: BaryonicMass_Interior(:)
  REAL(DP), ALLOCATABLE :: BaryonicMass_Initial (:)
  REAL(DP), ALLOCATABLE :: BaryonicMass_OffGrid (:)
  REAL(DP), ALLOCATABLE :: BaryonicMass_Change  (:)

  CHARACTER(256) :: Energy_FileName
  REAL(DP), ALLOCATABLE :: Energy_Interior(:)
  REAL(DP), ALLOCATABLE :: Energy_Initial (:)
  REAL(DP), ALLOCATABLE :: Energy_OffGrid (:)
  REAL(DP), ALLOCATABLE :: Energy_Change  (:)


CONTAINS


  SUBROUTINE MF_InitializeTally_Euler &
    ( SuppressTally_Option, BaseFileName_Option )

    LOGICAL,  INTENT(in),         OPTIONAL :: &
      SuppressTally_Option
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: &
      BaseFileName_Option

    CHARACTER(256) :: BaseFileName
    INTEGER        :: FileUnit

    CHARACTER(256) :: TimeLabel
    CHARACTER(256) :: InteriorLabel, InitialLabel, OffGridLabel, ChangeLabel

    SuppressTally = .FALSE.
    IF( PRESENT( SuppressTally_Option ) ) &
      SuppressTally = SuppressTally_Option

    IF( SuppressTally ) RETURN

    ALLOCATE( BaryonicMass_Interior(0:nLevels-1) )
    ALLOCATE( BaryonicMass_Initial (0:nLevels-1) )
    ALLOCATE( BaryonicMass_OffGrid (0:nLevels-1) )
    ALLOCATE( BaryonicMass_Change  (0:nLevels-1) )

    ALLOCATE( Energy_Interior(0:nLevels-1) )
    ALLOCATE( Energy_Initial (0:nLevels-1) )
    ALLOCATE( Energy_OffGrid (0:nLevels-1) )
    ALLOCATE( Energy_Change  (0:nLevels-1) )

    IF( amrex_parallel_ioprocessor() )THEN

      BaseFileName = ''
      IF( PRESENT( BaseFileName_Option ) ) &
        BaseFileName = TRIM( BaseFileName_Option )

      BaseFileName = TRIM( BaseFileName ) // TRIM( ProgramName )

      ! --- Baryonic Mass ---

      BaryonicMass_FileName &
        = TRIM( BaseFileName ) // '_Tally_BaryonicMass.dat'

      TimeLabel     &
        = 'Time ['     // TRIM( UnitsDisplay % TimeLabel ) // ']'
      InteriorLabel &
        = 'Interior [' // TRIM( UnitsDisplay % MassLabel ) // ']'
      OffGridLabel  &
        = 'Off Grid [' // TRIM( UnitsDisplay % MassLabel ) // ']'
      InitialLabel  &
        = 'Initial ['  // TRIM( UnitsDisplay % MassLabel ) // ']'
      ChangeLabel   &
        = 'Change ['   // TRIM( UnitsDisplay % MassLabel ) // ']'

      OPEN( NEWUNIT = FileUnit, FILE = TRIM( BaryonicMass_FileName ) )

      WRITE(FileUnit,'(5(A25,x))') &
        TRIM( TimeLabel ), TRIM( InteriorLabel ), TRIM( OffGridLabel ), &
        TRIM( InitialLabel ), TRIM( ChangeLabel )

      CLOSE( FileUnit )

      ! --- Energy ---

      Energy_FileName &
        = TRIM( BaseFileName ) // '_Tally_Energy.dat'

      TimeLabel     &
        = 'Time ['     // TRIM( UnitsDisplay % TimeLabel ) // ']'
      InteriorLabel &
        = 'Interior [' // TRIM( UnitsDisplay % EnergyGlobalLabel ) // ']'
      OffGridLabel  &
        = 'Off Grid [' // TRIM( UnitsDisplay % EnergyGlobalLabel ) // ']'
      InitialLabel  &
        = 'Initial ['  // TRIM( UnitsDisplay % EnergyGlobalLabel ) // ']'
      ChangeLabel   &
        = 'Change ['   // TRIM( UnitsDisplay % EnergyGlobalLabel ) // ']'

      OPEN( NEWUNIT = FileUnit, FILE = TRIM( Energy_FileName ) )

      WRITE(FileUnit,'(5(A25,x))') &
        TRIM( TimeLabel ), TRIM( InteriorLabel ), TRIM( OffGridLabel ), &
        TRIM( InitialLabel ), TRIM( ChangeLabel )

      CLOSE( FileUnit )

    END IF

    BaryonicMass_Interior = Zero
    BaryonicMass_Initial  = Zero
    BaryonicMass_OffGrid  = Zero
    BaryonicMass_Change   = Zero

    Energy_Interior = Zero
    Energy_Initial  = Zero
    Energy_OffGrid  = Zero
    Energy_Change   = Zero

  END SUBROUTINE MF_InitializeTally_Euler


  SUBROUTINE MF_ComputeTally_Euler &
    ( GEOM, MF_uGF, MF_uCF, Time, SetInitialValues_Option, Verbose_Option )

    TYPE(amrex_geometry), INTENT(in) :: GEOM  (0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in) :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in) :: MF_uCF(0:nLevels-1)
    REAL(DP),             INTENT(in) :: Time
    LOGICAL,              INTENT(in), OPTIONAL :: SetInitialValues_Option
    LOGICAL,              INTENT(in), OPTIONAL :: Verbose_Option

    LOGICAL :: SetInitialValues
    LOGICAL :: Verbose

    INTEGER                       :: iX_B0(3), iX_E0(3)
    INTEGER                       :: iLevel, iLo_MF(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    REAL(DP), ALLOCATABLE         :: G(:,:,:,:,:)
    REAL(DP), ALLOCATABLE         :: U(:,:,:,:,:)

    IF( SuppressTally ) RETURN

    SetInitialValues = .FALSE.
    IF( PRESENT( SetInitialValues_Option ) ) &
      SetInitialValues = SetInitialValues_Option

    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      BaryonicMass_Interior(iLevel) = Zero
      Energy_Interior      (iLevel) = Zero

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCF => MF_uCF(iLevel) % DataPtr( MFI )

        iLo_MF = LBOUND( uGF )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi

        CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_Allocate )

        ALLOCATE( G(1:nDOFX,iX_B0(1):iX_E0(1), &
                            iX_B0(2):iX_E0(2), &
                            iX_B0(3):iX_E0(3),1:nGF) )

        ALLOCATE( U(1:nDOFX,iX_B0(1):iX_E0(1), &
                            iX_B0(2):iX_E0(2), &
                            iX_B0(3):iX_E0(3),1:nCF) )

        CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_Allocate )

        CALL amrex2thornado_X( nGF, iX_B0, iX_E0, iLo_MF, iX_B0, iX_E0, uGF, G )
        CALL amrex2thornado_X( nCF, iX_B0, iX_E0, iLo_MF, iX_B0, iX_E0, uCF, U )

        CALL ComputeTally_Euler( iX_B0, iX_E0, G, U, iLevel )

        CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_Allocate )

        DEALLOCATE( U )
        DEALLOCATE( G )

        CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_Allocate )

      END DO

      CALL amrex_mfiter_destroy( MFI )

    END DO

    IF( SetInitialValues )THEN

      DO iLevel = 0, nLevels-1

        BaryonicMass_Initial(iLevel) = BaryonicMass_Interior(iLevel)
        Energy_Initial      (iLevel) = Energy_Interior      (iLevel)

      END DO

      CALL amrex_parallel_reduce_sum( BaryonicMass_Initial, nLevels )
      CALL amrex_parallel_reduce_sum( Energy_Initial      , nLevels )

    END IF

    CALL amrex_parallel_reduce_sum( BaryonicMass_Interior, nLevels )
    CALL amrex_parallel_reduce_sum( Energy_Interior      , nLevels )

    DO iLevel = 0, nLevels-1

      BaryonicMass_Change(iLevel) &
        = BaryonicMass_Interior(iLevel) &
            - ( BaryonicMass_Initial(iLevel) + BaryonicMass_OffGrid(iLevel) )

      Energy_Change(iLevel) &
        = Energy_Interior(iLevel) &
            - ( Energy_Initial(iLevel) + Energy_OffGrid(iLevel) )

    END DO

    CALL WriteTally_Euler( Time )

    IF( Verbose )THEN

      CALL DisplayTally( Time )

    END IF

  END SUBROUTINE MF_ComputeTally_Euler


  SUBROUTINE MF_IncrementOffGridTally_Euler( dM )

    REAL(DP), INTENT(in) :: dM(0:nLevels-1,nCF)

    INTEGER :: iLevel

    IF( SuppressTally ) RETURN

    DO iLevel = 0, nLevels-1

      BaryonicMass_OffGrid(iLevel) &
        = BaryonicMass_OffGrid(iLevel) + dM(iLevel,iCF_D)

      Energy_OffGrid(iLevel) &
        = Energy_OffGrid(iLevel) + dM(iLevel,iCF_E)

    END DO

  END SUBROUTINE MF_IncrementOffGridTally_Euler


  SUBROUTINE MF_FinalizeTally_Euler

    IF( SuppressTally ) RETURN

    DEALLOCATE( Energy_Change )
    DEALLOCATE( Energy_OffGrid )
    DEALLOCATE( Energy_Initial )
    DEALLOCATE( Energy_Interior )

    DEALLOCATE( BaryonicMass_Change )
    DEALLOCATE( BaryonicMass_OffGrid )
    DEALLOCATE( BaryonicMass_Initial )
    DEALLOCATE( BaryonicMass_Interior )

  END SUBROUTINE MF_FinalizeTally_Euler


  ! --- PRIVATE Subroutines ---


  SUBROUTINE ComputeTally_Euler( iX_B0, iX_E0, G, U, iLevel )

    INTEGER,  INTENT(in) :: &
      iX_B0(3), iX_E0(3), iLevel
    REAL(DP), INTENT(in) :: &
      G(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)
    REAL(DP), INTENT(in) :: &
      U(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)

    TYPE(MeshType) :: MeshX(3)
    INTEGER        :: iNX, iX1, iX2, iX3, iDim
    REAL(DP)       :: d3X

    DO iDim = 1, 3

      CALL CreateMesh &
             ( MeshX(iDim), nX(iDim), nNodesX(iDim), 0, &
               xL(iDim), xR(iDim) )

    END DO

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1, nDOFX

      IF( TRIM( CoordinateSystem ) .EQ. 'SPHERICAL' .AND. nDimsX .EQ. 1 )THEN

        d3X = FourPi * MeshX(1) % Width(iX1)

      ELSE

        d3X = MeshX(1) % Width(iX1) &
                * MeshX(2) % Width(iX2) &
                * MeshX(3) % Width(iX3)

      END IF

      BaryonicMass_Interior(iLevel) &
        = BaryonicMass_Interior(iLevel) &
            + d3X &
                * WeightsX_q(iNX) &
                * G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                * U(iNX,iX1,iX2,iX3,iCF_D)

      Energy_Interior(iLevel) &
        = Energy_Interior(iLevel) &
            + d3X &
                * WeightsX_q(iNX) &
                * G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                * U(iNX,iX1,iX2,iX3,iCF_E)

    END DO
    END DO
    END DO
    END DO

    DO iDim = 1, 3

      CALL DestroyMesh( MeshX(iDim) )

    END DO

  END SUBROUTINE ComputeTally_Euler


  SUBROUTINE WriteTally_Euler( Time )

    REAL(DP), INTENT(in) :: Time

    INTEGER :: FileUnit

    IF( amrex_parallel_ioprocessor() )THEN

      ! --- Baryonic Mass ---

      OPEN( NEWUNIT = FileUnit, FILE = TRIM( BaryonicMass_FileName ), &
            POSITION = 'APPEND', ACTION = 'WRITE' )

      WRITE( FileUnit, '(5(ES25.16E3,1x))' ) &
        Time / UnitsDisplay % TimeUnit, &
        BaryonicMass_Interior(0) / UnitsDisplay % MassUnit, &
        BaryonicMass_OffGrid (0) / UnitsDisplay % MassUnit, &
        BaryonicMass_Initial (0) / UnitsDisplay % MassUnit, &
        BaryonicMass_Change  (0) / UnitsDisplay % MassUnit

      CLOSE( FileUnit )

      ! --- Energy ---

      OPEN( NEWUNIT = FileUnit, FILE = TRIM( Energy_FileName ), &
            POSITION = 'APPEND', ACTION = 'WRITE' )

      WRITE( FileUnit, '(5(ES25.16E3,1x))' ) &
        Time / UnitsDisplay % TimeUnit, &
        Energy_Interior(0) / UnitsDisplay % EnergyGlobalUnit, &
        Energy_OffGrid (0) / UnitsDisplay % EnergyGlobalUnit, &
        Energy_Initial (0) / UnitsDisplay % EnergyGlobalUnit, &
        Energy_Change  (0) / UnitsDisplay % EnergyGlobalUnit

      CLOSE( FileUnit )

    END IF

  END SUBROUTINE WriteTally_Euler


  SUBROUTINE DisplayTally( Time )

    REAL(DP), INTENT(in) :: Time

    IF( amrex_parallel_ioprocessor() )THEN

      WRITE(*,*)
      WRITE(*,'(A8,A,ES8.2E2,x,A)') &
        '', 'Euler Tally. t = ', &
        Time / UnitsDisplay % TimeUnit, &
        UnitsDisplay % TimeLabel
      WRITE(*,*)
      WRITE(*,'(A6,A40,ES14.7E2,x,A)') &
        '', 'Baryonic Mass Interior.: ', &
        BaryonicMass_Interior(0) / UnitsDisplay % MassUnit, &
        UnitsDisplay % MassLabel
      WRITE(*,'(A6,A40,ES14.7E2,x,A)') &
        '', 'Baryonic Mass Initial..: ', &
        BaryonicMass_Initial(0)  / UnitsDisplay % MassUnit, &
        UnitsDisplay % MassLabel
      WRITE(*,'(A6,A40,ES14.7E2,x,A)') &
        '', 'Baryonic Mass Off Grid.: ', &
        BaryonicMass_OffGrid(0)  / UnitsDisplay % MassUnit, &
        UnitsDisplay % MassLabel
      WRITE(*,'(A6,A40,ES14.7E2,x,A)') &
        '', 'Baryonic Mass Change...: ', &
        BaryonicMass_Change(0)   / UnitsDisplay % MassUnit, &
        UnitsDisplay % MassLabel

      WRITE(*,*)
      WRITE(*,'(A6,A40,ES14.7E2,x,A)') &
        '', 'Energy Interior.: ', &
        Energy_Interior(0) / UnitsDisplay % EnergyGlobalUnit, &
        UnitsDisplay % EnergyGlobalLabel
      WRITE(*,'(A6,A40,ES14.7E2,x,A)') &
        '', 'Energy Initial..: ', &
        Energy_Initial(0)  / UnitsDisplay % EnergyGlobalUnit, &
        UnitsDisplay % EnergyGlobalLabel
      WRITE(*,'(A6,A40,ES14.7E2,x,A)') &
        '', 'Energy Off Grid.: ', &
        Energy_OffGrid(0)  / UnitsDisplay % EnergyGlobalUnit, &
        UnitsDisplay % EnergyGlobalLabel
      WRITE(*,'(A6,A40,ES14.7E2,x,A)') &
        '', 'Energy Change...: ', &
        Energy_Change(0)   / UnitsDisplay % EnergyGlobalUnit, &
        UnitsDisplay % EnergyGlobalLabel

      WRITE(*,*)

    END IF

  END SUBROUTINE DisplayTally


END MODULE MF_Euler_TallyModule
