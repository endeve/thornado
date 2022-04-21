MODULE MF_Euler_TallyModule

  ! --- AMReX Modules ---

  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_geometry_module, ONLY: &
    amrex_geometry
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
    Zero
  USE InputParsingModule, ONLY: &
    nX, &
    nLevels, &
    ProgramName, &
    UseTiling, &
    xL, &
    xR
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF
!  USE TimersModule_AMReX_Euler, ONLY: &
!    TimersStart_AMReX_Euler, &
!    TimersStop_AMReX_Euler, &
!    Timer_AMReX_Euler_Allocate
  USE MakeFineMaskModule, ONLY: &
    MakeFineMask, &
    DestroyFineMask, &
    iLeaf_MFM
  USE MF_UtilitiesModule, ONLY: &
    amrex2thornado_X

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeTally_Euler_MF
  PUBLIC :: ComputeTally_Euler_MF
  PUBLIC :: IncrementOffGridTally_Euler_MF
  PUBLIC :: FinalizeTally_Euler_MF

  LOGICAL :: SuppressTally

  CHARACTER(256) :: BaryonicMass_FileName
  REAL(DP) :: BaryonicMass_Interior
  REAL(DP) :: BaryonicMass_Initial
  REAL(DP) :: BaryonicMass_OffGrid
  REAL(DP) :: BaryonicMass_Change

  CHARACTER(256) :: Energy_FileName
  REAL(DP) :: Energy_Interior
  REAL(DP) :: Energy_Initial
  REAL(DP) :: Energy_OffGrid
  REAL(DP) :: Energy_Change

CONTAINS


  SUBROUTINE InitializeTally_Euler_MF &
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

    IF( nLevels .GT. 1 )THEN

      IF( amrex_parallel_ioprocessor() )THEN

        WRITE(*,*)
        WRITE(*,*) &
          'WARNING: Euler_TallyModule not accurate for multi-level mesh'
        WRITE(*,*)

      END IF

    END IF

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

  END SUBROUTINE InitializeTally_Euler_MF


  SUBROUTINE ComputeTally_Euler_MF &
    ( Time, MF_uGF, MF_uCF, SetInitialValues_Option, Verbose_Option )

    REAL(DP),             INTENT(in) :: Time  (0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in) :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in) :: MF_uCF(0:nLevels-1)
    LOGICAL,              INTENT(in), OPTIONAL :: SetInitialValues_Option
    LOGICAL,              INTENT(in), OPTIONAL :: Verbose_Option

    LOGICAL :: SetInitialValues
    LOGICAL :: Verbose

    INTEGER                       :: iX_B0(3), iX_E0(3)
    INTEGER                       :: iLevel, iLo_MF(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    INTEGER , CONTIGUOUS, POINTER :: Mask(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uGF (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF (:,:,:,:)
    REAL(DP), ALLOCATABLE         :: G (:,:,:,:,:)
    REAL(DP), ALLOCATABLE         :: U (:,:,:,:,:)

    TYPE(amrex_imultifab) :: iMF_Mask

    IF( SuppressTally ) RETURN

    SetInitialValues = .FALSE.
    IF( PRESENT( SetInitialValues_Option ) ) &
      SetInitialValues = SetInitialValues_Option

    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    BaryonicMass_Interior = Zero
    Energy_Interior       = Zero

    DO iLevel = 0, nLevels-1

      CALL MakeFineMask( iLevel, iMF_Mask, MF_uCF % BA, MF_uCF % DM )

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        Mask => iMF_Mask % DataPtr( MFI )
        uGF  => MF_uGF(iLevel) % DataPtr( MFI )
        uCF  => MF_uCF(iLevel) % DataPtr( MFI )

        iLo_MF = LBOUND( uGF )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi

!        CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_Allocate )

        ALLOCATE( G(1:nDOFX,iX_B0(1):iX_E0(1), &
                            iX_B0(2):iX_E0(2), &
                            iX_B0(3):iX_E0(3),1:nGF) )

        ALLOCATE( U(1:nDOFX,iX_B0(1):iX_E0(1), &
                            iX_B0(2):iX_E0(2), &
                            iX_B0(3):iX_E0(3),1:nCF) )

!        CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_Allocate )

        CALL amrex2thornado_X( nGF, iX_B0, iX_E0, iLo_MF, iX_B0, iX_E0, uGF, G )
        CALL amrex2thornado_X( nCF, iX_B0, iX_E0, iLo_MF, iX_B0, iX_E0, uCF, U )

        CALL ComputeTally_Euler( iX_B0, iX_E0, G, U, Mask, iLevel )

!        CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_Allocate )

        DEALLOCATE( U )
        DEALLOCATE( G )

!        CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_Allocate )

      END DO

      CALL amrex_mfiter_destroy( MFI )

      CALL DestroyFineMask( iLevel, iMF_Mask )

    END DO

    CALL amrex_parallel_reduce_sum( BaryonicMass_Interior )
    CALL amrex_parallel_reduce_sum( Energy_Interior       )

    IF( SetInitialValues )THEN

      BaryonicMass_Initial = BaryonicMass_Interior
      Energy_Initial       = Energy_Interior

    END IF

    BaryonicMass_Change &
      = BaryonicMass_Interior &
          - ( BaryonicMass_Initial + BaryonicMass_OffGrid )

    Energy_Change &
      = Energy_Interior &
          - ( Energy_Initial + Energy_OffGrid )

    CALL WriteTally_Euler( Time(0) )

    IF( Verbose )THEN

      CALL DisplayTally( Time(0) )

    END IF

  END SUBROUTINE ComputeTally_Euler_MF


  SUBROUTINE IncrementOffGridTally_Euler_MF( dM )

    REAL(DP), INTENT(in) :: dM(nCF,0:nLevels-1)

    INTEGER :: iLevel

    IF( SuppressTally ) RETURN

    DO iLevel = 0, nLevels-1

      BaryonicMass_OffGrid &
        = BaryonicMass_OffGrid + dM(iCF_D,iLevel)

      Energy_OffGrid &
        = Energy_OffGrid + dM(iCF_E,iLevel)

    END DO

    CALL amrex_parallel_reduce_sum( BaryonicMass_OffGrid )
    CALL amrex_parallel_reduce_sum( Energy_OffGrid       )

  END SUBROUTINE IncrementOffGridTally_Euler_MF


  SUBROUTINE FinalizeTally_Euler_MF

    IF( SuppressTally ) RETURN

  END SUBROUTINE FinalizeTally_Euler_MF


  ! --- PRIVATE Subroutines ---


  SUBROUTINE ComputeTally_Euler( iX_B0, iX_E0, G, U, Mask, iLevel )

    INTEGER,  INTENT(in) :: &
      iX_B0(3), iX_E0(3), iLevel
    REAL(DP), INTENT(in) :: &
      G(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)
    REAL(DP), INTENT(in) :: &
      U(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)
    INTEGER , INTENT(in) :: &
      Mask(iX_B0(1):,iX_B0(2):,iX_B0(3):,:)

    TYPE(MeshType) :: MeshX(3)
    INTEGER        :: iNX, iX1, iX2, iX3
    REAL(DP)       :: d3X

    CALL CreateMesh_MF( iLevel, MeshX )

    d3X =   MeshX(1) % Width(iX_B0(1)) &
          * MeshX(2) % Width(iX_B0(2)) &
          * MeshX(3) % Width(iX_B0(3))

    CALL DestroyMesh_MF( MeshX )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1, nDOFX

      IF( nLevels .GT. 1 .AND. iLevel .LT. nLevels-1 )THEN

        IF( Mask(iX1,iX2,iX3,1) .NE. iLeaf_MFM ) CYCLE

      END IF

      BaryonicMass_Interior &
        = BaryonicMass_Interior &
            + d3X &
                * WeightsX_q(iNX) &
                * G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                * U(iNX,iX1,iX2,iX3,iCF_D)

      Energy_Interior &
        = Energy_Interior &
            + d3X &
                * WeightsX_q(iNX) &
                * G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                * U(iNX,iX1,iX2,iX3,iCF_E)

    END DO
    END DO
    END DO
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
        BaryonicMass_Interior / UnitsDisplay % MassUnit, &
        UnitsDisplay % MassLabel
      WRITE(*,'(A6,A40,ES14.7E2,x,A)') &
        '', 'Baryonic Mass Initial..: ', &
        BaryonicMass_Initial  / UnitsDisplay % MassUnit, &
        UnitsDisplay % MassLabel
      WRITE(*,'(A6,A40,ES14.7E2,x,A)') &
        '', 'Baryonic Mass Off Grid.: ', &
        BaryonicMass_OffGrid  / UnitsDisplay % MassUnit, &
        UnitsDisplay % MassLabel
      WRITE(*,'(A6,A40,ES14.7E2,x,A)') &
        '', 'Baryonic Mass Change...: ', &
        BaryonicMass_Change   / UnitsDisplay % MassUnit, &
        UnitsDisplay % MassLabel

      WRITE(*,*)
      WRITE(*,'(A6,A40,ES14.7E2,x,A)') &
        '', 'Energy Interior.: ', &
        Energy_Interior / UnitsDisplay % EnergyGlobalUnit, &
        UnitsDisplay % EnergyGlobalLabel
      WRITE(*,'(A6,A40,ES14.7E2,x,A)') &
        '', 'Energy Initial..: ', &
        Energy_Initial  / UnitsDisplay % EnergyGlobalUnit, &
        UnitsDisplay % EnergyGlobalLabel
      WRITE(*,'(A6,A40,ES14.7E2,x,A)') &
        '', 'Energy Off Grid.: ', &
        Energy_OffGrid  / UnitsDisplay % EnergyGlobalUnit, &
        UnitsDisplay % EnergyGlobalLabel
      WRITE(*,'(A6,A40,ES14.7E2,x,A)') &
        '', 'Energy Change...: ', &
        Energy_Change   / UnitsDisplay % EnergyGlobalUnit, &
        UnitsDisplay % EnergyGlobalLabel

      WRITE(*,*)

    END IF

  END SUBROUTINE DisplayTally


END MODULE MF_Euler_TallyModule
