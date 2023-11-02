MODULE MF_Euler_TallyModule

  ! --- AMReX Modules ---

  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_geometry_module, ONLY: &
    amrex_geometry
  USE amrex_parmparse_module, ONLY: &
    amrex_parmparse, &
    amrex_parmparse_build, &
    amrex_parmparse_destroy
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_imultifab, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor, &
    amrex_parallel_reduce_sum
  USE amrex_amr_module, ONLY: &
    amrex_geom

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
    iGF_SqrtGm, &
    nGF, &
    CoordinateSystem
  USE FluidFieldsModule, ONLY: &
    iCF_D, &
    iCF_E, &
    iCF_Ne, &
    nCF

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero
  USE InputParsingModule, ONLY: &
    nX, &
    nLevels, &
    ProgramName, &
    UseTiling, &
    xL, &
    xR, &
    swX
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF
  USE MaskModule, ONLY: &
    CreateFineMask, &
    DestroyFineMask, &
    IsNotLeafElement
  USE MF_UtilitiesModule, ONLY: &
    amrex2thornado_X, &
    AllocateArray_X, &
    DeallocateArray_X

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

  USE ADM_Mass_Module, ONLY: &
    Calc_ADM_Mass

#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeTally_Euler_MF
  PUBLIC :: ComputeTally_Euler_MF
  PUBLIC :: IncrementOffGridTally_Euler_MF
  PUBLIC :: FinalizeTally_Euler_MF

  LOGICAL :: SuppressTally

  INTEGER, PARAMETER :: SL = 256

  CHARACTER(SL)    :: BaryonicMass_FileName
  REAL(DP), PUBLIC :: BaryonicMass_Initial
  REAL(DP), PUBLIC :: BaryonicMass_OffGrid
  REAL(DP)         :: BaryonicMass_Interior
  REAL(DP)         :: BaryonicMass_Interior_OMP
  REAL(DP)         :: BaryonicMass_Change

  CHARACTER(SL)    :: Energy_FileName
  REAL(DP), PUBLIC :: Energy_Initial
  REAL(DP), PUBLIC :: Energy_OffGrid
  REAL(DP)         :: Energy_Interior
  REAL(DP)         :: Energy_Interior_OMP
  REAL(DP)         :: Energy_Change

  CHARACTER(SL)    :: ElectronNumber_FileName
  REAL(DP), PUBLIC :: ElectronNumber_Initial
  REAL(DP), PUBLIC :: ElectronNumber_OffGrid
  REAL(DP)         :: ElectronNumber_Interior
  REAL(DP)         :: ElectronNumber_Interior_OMP
  REAL(DP)         :: ElectronNumber_Change

  CHARACTER(SL)    :: ADMMass_FileName
  REAL(DP), PUBLIC :: ADMMass_Initial
  REAL(DP), PUBLIC :: ADMMass_OffGrid
  REAL(DP)         :: ADMMass_Interior
  REAL(DP)         :: ADMMass_Change

CONTAINS


  SUBROUTINE InitializeTally_Euler_MF &
    ( SuppressTally_Option, InitializeFromCheckpoint_Option )

    LOGICAL, INTENT(in), OPTIONAL :: &
      SuppressTally_Option, &
      InitializeFromCheckpoint_Option

    CHARACTER(:), ALLOCATABLE :: TallyFileNameRoot_Euler

    CHARACTER(SL) :: FileNameRoot
    INTEGER       :: FileUnit

    LOGICAL :: InitializeFromCheckpoint

    TYPE(amrex_parmparse) :: PP

    CHARACTER(SL) :: TimeLabel
    CHARACTER(SL) :: InteriorLabel, InitialLabel, OffGridLabel, ChangeLabel

    SuppressTally = .FALSE.
    IF( PRESENT( SuppressTally_Option ) ) &
      SuppressTally = SuppressTally_Option

    InitializeFromCheckpoint = .FALSE.
    IF( PRESENT( InitializeFromCheckpoint_Option ) ) &
      InitializeFromCheckpoint = InitializeFromCheckpoint_Option

    TallyFileNameRoot_Euler = TRIM( ProgramName )
    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL pp % query( 'TallyFileNameRoot_Euler', &
                        TallyFileNameRoot_Euler )
    CALL amrex_parmparse_destroy( PP )

    IF( SuppressTally ) RETURN

    IF( amrex_parallel_ioprocessor() )THEN

      FileNameRoot = TRIM( TallyFileNameRoot_Euler )

      ! --- Baryonic Mass ---

      BaryonicMass_FileName &
        = TRIM( FileNameRoot ) // '_BaryonicMass.dat'

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

      CALL CheckFileExistenceAndAppend( BaryonicMass_FileName )

      OPEN( NEWUNIT = FileUnit, FILE = TRIM( BaryonicMass_FileName ) )

      WRITE(FileUnit,'(5(A25,x))') &
        TRIM( TimeLabel ), TRIM( InteriorLabel ), TRIM( OffGridLabel ), &
        TRIM( InitialLabel ), TRIM( ChangeLabel )

      CLOSE( FileUnit )

      ! --- Energy ---

      Energy_FileName &
        = TRIM( FileNameRoot ) // '_Energy.dat'

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

      CALL CheckFileExistenceAndAppend( Energy_FileName )

      OPEN( NEWUNIT = FileUnit, FILE = TRIM( Energy_FileName ) )

      WRITE(FileUnit,'(5(A25,x))') &
        TRIM( TimeLabel ), TRIM( InteriorLabel ), TRIM( OffGridLabel ), &
        TRIM( InitialLabel ), TRIM( ChangeLabel )

      CLOSE( FileUnit )

      ! --- Electron Number ---

      ElectronNumber_FileName &
        = TRIM( FileNameRoot ) // '_ElectronNumber.dat'

      TimeLabel     &
        = 'Time ['     // TRIM( UnitsDisplay % TimeLabel ) // ']'
      InteriorLabel &
        = 'Interior [' // '' // ']'
      OffGridLabel  &
        = 'Off Grid [' // '' // ']'
      InitialLabel  &
        = 'Initial ['  // '' // ']'
      ChangeLabel   &
        = 'Change ['   // '' // ']'

      CALL CheckFileExistenceAndAppend( ElectronNumber_FileName )

      OPEN( NEWUNIT = FileUnit, FILE = TRIM( ElectronNumber_FileName ) )

      WRITE(FileUnit,'(5(A25,x))') &
        TRIM( TimeLabel ), TRIM( InteriorLabel ), TRIM( OffGridLabel ), &
        TRIM( InitialLabel ), TRIM( ChangeLabel )

      CLOSE( FileUnit )

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

      ! --- ADM Mass ---

      ADMMass_FileName &
        = TRIM( FileNameRoot ) // '_ADMMass.dat'

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

      CALL CheckFileExistenceAndAppend( ADMMass_FileName )

      OPEN( NEWUNIT = FileUnit, FILE = TRIM( ADMMass_FileName ) )

      WRITE(FileUnit,'(5(A25,x))') &
        TRIM( TimeLabel ), TRIM( InteriorLabel ), TRIM( OffGridLabel ), &
        TRIM( InitialLabel ), TRIM( ChangeLabel )

      CLOSE( FileUnit )

#endif

    END IF

    IF( .NOT. InitializeFromCheckpoint )THEN

      BaryonicMass_Initial = Zero
      BaryonicMass_OffGrid = Zero

      Energy_Initial = Zero
      Energy_OffGrid = Zero

      ElectronNumber_Initial = Zero
      ElectronNumber_OffGrid = Zero

      ADMMass_Initial = Zero
      ADMMass_OffGrid = Zero

    END IF

    BaryonicMass_Interior = Zero
    BaryonicMass_Change   = Zero

    Energy_Interior = Zero
    Energy_Change   = Zero

    ElectronNumber_Interior = Zero
    ElectronNumber_Change   = Zero

    ADMMass_Interior = Zero
    ADMMass_Change   = Zero

  END SUBROUTINE InitializeTally_Euler_MF


  SUBROUTINE ComputeTally_Euler_MF &
    ( Time, MF_uGF, MF_uCF, SetInitialValues_Option, &
      WriteTally_Option, Verbose_Option )

    REAL(DP),             INTENT(in) :: Time  (0:)
    TYPE(amrex_multifab), INTENT(in) :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(in) :: MF_uCF(0:)
    LOGICAL,              INTENT(in), OPTIONAL :: SetInitialValues_Option
    LOGICAL,              INTENT(in), OPTIONAL :: WriteTally_Option
    LOGICAL,              INTENT(in), OPTIONAL :: Verbose_Option

    LOGICAL :: SetInitialValues
    LOGICAL :: Verbose
    LOGICAL :: WriteTally

    INTEGER                       :: iX_B0(3), iX_E0(3)
    INTEGER                       :: iX_B1(3), iX_E1(3)
    INTEGER                       :: iLevel, iLo_MF(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    INTEGER , CONTIGUOUS, POINTER :: FineMask(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uGF     (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF     (:,:,:,:)
    REAL(DP), ALLOCATABLE         :: G     (:,:,:,:,:)
    REAL(DP), ALLOCATABLE         :: U     (:,:,:,:,:)

    TYPE(amrex_imultifab) :: iMF_FineMask

    TYPE(MeshType) :: MeshX(3)
    INTEGER        :: iNX, iX1, iX2, iX3
    REAL(DP)       :: d3X

    IF( SuppressTally ) RETURN

    SetInitialValues = .FALSE.
    IF( PRESENT( SetInitialValues_Option ) ) &
      SetInitialValues = SetInitialValues_Option

    WriteTally = .TRUE.
    IF( PRESENT( WriteTally_Option ) ) &
      WriteTally = WriteTally_Option

    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    BaryonicMass_Interior   = Zero
    Energy_Interior         = Zero
    ElectronNumber_Interior = Zero
    ADMMass_Interior        = Zero

    DO iLevel = 0, nLevels-1

      CALL CreateFineMask( iLevel, iMF_FineMask, MF_uCF % BA, MF_uCF % DM )

      CALL CreateMesh_MF( iLevel, MeshX )

      BaryonicMass_Interior_OMP   = Zero
      Energy_Interior_OMP         = Zero
      ElectronNumber_Interior_OMP = Zero

#if defined( THORNADO_OMP )
      !$OMP PARALLEL &
      !$OMP PRIVATE( iX_B0, iX_E0, iX_B1, iX_E1, iLo_MF, &
      !$OMP          BX, MFI, FineMask, uGF, uCF, G, U, d3X ) &
      !$OMP REDUCTION( +:BaryonicMass_Interior_OMP, &
      !$OMP              Energy_Interior_OMP, &
      !$OMP              ElectronNumber_Interior_OMP )
#endif

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        FineMask => iMF_FineMask   % DataPtr( MFI )
        uGF      => MF_uGF(iLevel) % DataPtr( MFI )
        uCF      => MF_uCF(iLevel) % DataPtr( MFI )

        iLo_MF = LBOUND( uGF )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi

        iX_B1 = iX_B0 - swX
        iX_E1 = iX_E0 + swX

        CALL AllocateArray_X &
               ( [ 1    , iX_B0(1), iX_B0(2), iX_B0(3), 1   ], &
                 [ nDOFX, iX_E0(1), iX_E0(2), iX_E0(3), nGF ], &
                 G )

        CALL AllocateArray_X &
               ( [ 1    , iX_B0(1), iX_B0(2), iX_B0(3), 1   ], &
                 [ nDOFX, iX_E0(1), iX_E0(2), iX_E0(3), nCF ], &
                 U )

        CALL amrex2thornado_X( nGF, iX_B0, iX_E0, iLo_MF, iX_B0, iX_E0, uGF, G )
        CALL amrex2thornado_X( nCF, iX_B0, iX_E0, iLo_MF, iX_B0, iX_E0, uCF, U )

        d3X =   MeshX(1) % Width(iX_B0(1)) &
              * MeshX(2) % Width(iX_B0(2)) &
              * MeshX(3) % Width(iX_B0(3))

        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1       , nDOFX

          IF( IsNotLeafElement( FineMask(iX1,iX2,iX3,1) ) ) CYCLE

          BaryonicMass_Interior_OMP &
            = BaryonicMass_Interior_OMP &
                + d3X &
                    * WeightsX_q(iNX) &
                    * G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                    * U(iNX,iX1,iX2,iX3,iCF_D)

          Energy_Interior_OMP &
            = Energy_Interior_OMP &
                + d3X &
                    * WeightsX_q(iNX) &
                    * G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                    * U(iNX,iX1,iX2,iX3,iCF_E)

          ElectronNumber_Interior_OMP &
            = ElectronNumber_Interior_OMP &
                + d3X &
                    * WeightsX_q(iNX) &
                    * G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                    * U(iNX,iX1,iX2,iX3,iCF_Ne)

        END DO
        END DO
        END DO
        END DO

        CALL DeallocateArray_X &
               ( [ 1    , iX_B0(1), iX_B0(2), iX_B0(3), 1   ], &
                 [ nDOFX, iX_E0(1), iX_E0(2), iX_E0(3), nCF ], &
                 U )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B0(1), iX_B0(2), iX_B0(3), 1   ], &
                 [ nDOFX, iX_E0(1), iX_E0(2), iX_E0(3), nGF ], &
                 G )

      END DO ! WHILE( MFI % next() )

      CALL amrex_mfiter_destroy( MFI )

#if defined( THORNADO_OMP )
      !$OMP END PARALLEL
#endif

      BaryonicMass_Interior &
        = BaryonicMass_Interior   + BaryonicMass_Interior_OMP
      Energy_Interior &
        = Energy_Interior         + Energy_Interior_OMP
      ElectronNumber_Interior &
        = ElectronNumber_Interior + ElectronNumber_Interior_OMP

      CALL DestroyMesh_MF( MeshX )

      CALL DestroyFineMask( iMF_FineMask )

    END DO ! iLevel = 0, nLevels-1

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

    CALL Calc_ADM_Mass( ADMMass_Interior )

#endif

    CALL amrex_parallel_reduce_sum( BaryonicMass_Interior   )
    CALL amrex_parallel_reduce_sum( Energy_Interior         )
    CALL amrex_parallel_reduce_sum( ElectronNumber_Interior )

    IF( SetInitialValues )THEN

      BaryonicMass_Initial   = BaryonicMass_Interior
      Energy_Initial         = Energy_Interior
      ElectronNumber_Initial = ElectronNumber_Interior

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

      ADMMass_Initial = ADMMass_Interior

#endif

    END IF

    ! --- dM = Minterior - Minitial + ( OffGrid_Outer - OffGrid_Inner ) ---

    BaryonicMass_Change &
      = BaryonicMass_Interior &
          - BaryonicMass_Initial + BaryonicMass_OffGrid

    Energy_Change &
      = Energy_Interior &
          - Energy_Initial + Energy_OffGrid

    ElectronNumber_Change &
      = ElectronNumber_Interior &
          - ElectronNumber_Initial + ElectronNumber_OffGrid

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

    ADMMass_Change &
      = ADMMass_Interior &
          - ( ADMMass_Initial + ADMMass_OffGrid )

#endif

    IF( WriteTally ) &
      CALL WriteTally_Euler( Time(0) )

    IF( Verbose ) CALL DisplayTally( Time(0) )

  END SUBROUTINE ComputeTally_Euler_MF


  SUBROUTINE IncrementOffGridTally_Euler_MF( dM )

    REAL(DP), INTENT(in) :: dM(1:,0:)

    INTEGER :: iLevel

    IF( SuppressTally ) RETURN

    DO iLevel = 0, nLevels-1

      BaryonicMass_OffGrid &
        = BaryonicMass_OffGrid + dM(iCF_D,iLevel)

      Energy_OffGrid &
        = Energy_OffGrid + dM(iCF_E,iLevel)

      ElectronNumber_OffGrid &
        = ElectronNumber_OffGrid + dM(iCF_Ne,iLevel)

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

      ADMMass_OffGrid &
        = Zero

#endif

    END DO

  END SUBROUTINE IncrementOffGridTally_Euler_MF


  SUBROUTINE FinalizeTally_Euler_MF

    IF( SuppressTally ) RETURN

  END SUBROUTINE FinalizeTally_Euler_MF


  ! --- PRIVATE Subroutines ---


  SUBROUTINE WriteTally_Euler( Time )

    REAL(DP), INTENT(in) :: Time

    INTEGER :: FileUnit

    IF( amrex_parallel_ioprocessor() )THEN

      ! --- Baryonic Mass ---

      OPEN( NEWUNIT = FileUnit, FILE = TRIM( BaryonicMass_FileName ), &
            POSITION = 'APPEND', ACTION = 'WRITE' )

      WRITE( FileUnit, '(5(ES25.16E3,1x))' ) &
        Time / UnitsDisplay % TimeUnit, &
        BaryonicMass_Interior / UnitsDisplay % MassUnit, &
        BaryonicMass_OffGrid  / UnitsDisplay % MassUnit, &
        BaryonicMass_Initial  / UnitsDisplay % MassUnit, &
        BaryonicMass_Change   / UnitsDisplay % MassUnit

      CLOSE( FileUnit )

      ! --- Energy ---

      OPEN( NEWUNIT = FileUnit, FILE = TRIM( Energy_FileName ), &
            POSITION = 'APPEND', ACTION = 'WRITE' )

      WRITE( FileUnit, '(5(ES25.16E3,1x))' ) &
        Time / UnitsDisplay % TimeUnit, &
        Energy_Interior / UnitsDisplay % EnergyGlobalUnit, &
        Energy_OffGrid  / UnitsDisplay % EnergyGlobalUnit, &
        Energy_Initial  / UnitsDisplay % EnergyGlobalUnit, &
        Energy_Change   / UnitsDisplay % EnergyGlobalUnit

      CLOSE( FileUnit )

      CLOSE( FileUnit )

      ! --- Electron Number ---

      OPEN( NEWUNIT = FileUnit, FILE = TRIM( ElectronNumber_FileName ), &
            POSITION = 'APPEND', ACTION = 'WRITE' )

      WRITE( FileUnit, '(5(ES25.16E3,1x))' ) &
        Time / UnitsDisplay % TimeUnit, &
        ElectronNumber_Interior, &
        ElectronNumber_OffGrid , &
        ElectronNumber_Initial , &
        ElectronNumber_Change

      CLOSE( FileUnit )

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

      ! --- ADM Mass ---

      OPEN( NEWUNIT = FileUnit, FILE = TRIM( ADMMass_FileName ), &
            POSITION = 'APPEND', ACTION = 'WRITE' )

      WRITE( FileUnit, '(5(ES25.16E3,1x))' ) &
        Time / UnitsDisplay % TimeUnit, &
        ADMMass_Interior / UnitsDisplay % MassUnit, &
        ADMMass_OffGrid  / UnitsDisplay % MassUnit, &
        ADMMass_Initial  / UnitsDisplay % MassUnit, &
        ADMMass_Change   / UnitsDisplay % MassUnit

      CLOSE( FileUnit )

#endif

    END IF

  END SUBROUTINE WriteTally_Euler


  SUBROUTINE DisplayTally( Time )

    REAL(DP), INTENT(in) :: Time

    CHARACTER(32) :: FMT

    FMT = '(6x,A40,ES15.7E3,x,A)'

    IF( amrex_parallel_ioprocessor() )THEN

      WRITE(*,*)
      WRITE(*,'(A8,A,ES8.2E2,x,A)') &
        '', 'Euler Tally. t = ', &
        Time / UnitsDisplay % TimeUnit, &
        UnitsDisplay % TimeLabel
      WRITE(*,*)
      WRITE(*,TRIM(FMT)) &
        'Baryonic Mass Interior.: ', &
        BaryonicMass_Interior / UnitsDisplay % MassUnit, &
        UnitsDisplay % MassLabel
      WRITE(*,TRIM(FMT)) &
        'Baryonic Mass Initial..: ', &
        BaryonicMass_Initial  / UnitsDisplay % MassUnit, &
        UnitsDisplay % MassLabel
      WRITE(*,TRIM(FMT)) &
        'Baryonic Mass Off Grid.: ', &
        BaryonicMass_OffGrid  / UnitsDisplay % MassUnit, &
        UnitsDisplay % MassLabel
      WRITE(*,TRIM(FMT)) &
        'Baryonic Mass Change...: ', &
        BaryonicMass_Change   / UnitsDisplay % MassUnit, &
        UnitsDisplay % MassLabel

      WRITE(*,*)
      WRITE(*,TRIM(FMT)) &
        'Energy Interior.: ', &
        Energy_Interior / UnitsDisplay % EnergyGlobalUnit, &
        UnitsDisplay % EnergyGlobalLabel
      WRITE(*,TRIM(FMT)) &
        'Energy Initial..: ', &
        Energy_Initial  / UnitsDisplay % EnergyGlobalUnit, &
        UnitsDisplay % EnergyGlobalLabel
      WRITE(*,TRIM(FMT)) &
        'Energy Off Grid.: ', &
        Energy_OffGrid  / UnitsDisplay % EnergyGlobalUnit, &
        UnitsDisplay % EnergyGlobalLabel
      WRITE(*,TRIM(FMT)) &
        'Energy Change...: ', &
        Energy_Change   / UnitsDisplay % EnergyGlobalUnit, &
        UnitsDisplay % EnergyGlobalLabel

      WRITE(*,*)
      WRITE(*,TRIM(FMT)) &
        'Electron Number Interior.: ', &
        ElectronNumber_Interior, &
        ''
      WRITE(*,TRIM(FMT)) &
        'Electron Number Initial..: ', &
        ElectronNumber_Initial, &
        ''
      WRITE(*,TRIM(FMT)) &
        'Electron Number Off Grid.: ', &
        ElectronNumber_OffGrid, &
        ''
      WRITE(*,TRIM(FMT)) &
        'Electron Number Change...: ', &
        ElectronNumber_Change, &
        ''

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

      WRITE(*,*)
      WRITE(*,'(A6,A40,ES14.7E2,x,A)') &
        '', 'ADM Mass Interior.: ', &
        ADMMass_Interior / UnitsDisplay % MassUnit, &
        UnitsDisplay % MassLabel
      WRITE(*,'(A6,A40,ES14.7E2,x,A)') &
        '', 'ADM Mass Initial..: ', &
        ADMMass_Initial  / UnitsDisplay % MassUnit, &
        UnitsDisplay % MassLabel
      WRITE(*,'(A6,A40,ES14.7E2,x,A)') &
        '', 'ADM Mass Off Grid.: ', &
        ADMMass_OffGrid  / UnitsDisplay % MassUnit, &
        UnitsDisplay % MassLabel
      WRITE(*,'(A6,A40,ES14.7E2,x,A)') &
        '', 'ADM Mass Change...: ', &
        ADMMass_Change   / UnitsDisplay % MassUnit, &
        UnitsDisplay % MassLabel

#endif

      WRITE(*,*)

    END IF

  END SUBROUTINE DisplayTally


  RECURSIVE SUBROUTINE CheckFileExistenceAndAppend( FileName, IntSuffix_Option )

    CHARACTER(LEN=SL), INTENT(inout) :: FileName
    INTEGER          , INTENT(inout), OPTIONAL :: IntSuffix_Option

    LOGICAL :: IsFile
    INTEGER :: IntSuffix
    INTEGER :: SL_T

    IntSuffix = 1
    IF( PRESENT( IntSuffix_Option ) ) &
      IntSuffix = IntSuffix_Option

    SL_T = LEN( TRIM( FileName ) )

    INQUIRE( FILE = TRIM( FileName ), EXIST = IsFile )

    IF( IsFile )THEN

      IF( FileName(SL_T-3:SL_T) .EQ. '.dat' )THEN

        WRITE(FileName,'(A,A,I2.2)') TRIM( FileName ), '_', IntSuffix

      ELSE

        WRITE(FileName(SL_T-1:SL_T),'(I2.2)') IntSuffix

      END IF

      IntSuffix = IntSuffix + 1

      CALL CheckFileExistenceAndAppend &
             ( FileName, IntSuffix_Option = IntSuffix )

    END IF

  END SUBROUTINE CheckFileExistenceAndAppend


END MODULE MF_Euler_TallyModule
