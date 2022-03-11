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
    iCoarse_MFM
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
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    INTEGER , CONTIGUOUS, POINTER :: Mask(:,:,:,:)
    REAL(DP), ALLOCATABLE         :: G(:,:,:,:,:)
    REAL(DP), ALLOCATABLE         :: U(:,:,:,:,:)

    TYPE(amrex_imultifab) :: iMF_Mask

    IF( SuppressTally ) RETURN

    SetInitialValues = .FALSE.
    IF( PRESENT( SetInitialValues_Option ) ) &
      SetInitialValues = SetInitialValues_Option

    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    DO iLevel = 0, nLevels-1

      IF( nLevels .GT. 1 .AND. iLevel .LT. nLevels-1 ) &
        CALL MakeFineMask &
               ( iMF_Mask, MF_uCF(iLevel) % BA, MF_uCF(iLevel) % DM, &
                 MF_uCF(iLevel+1) % BA )

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      BaryonicMass_Interior(iLevel) = Zero
      Energy_Interior      (iLevel) = Zero

      DO WHILE( MFI % next() )

        IF( nLevels-1 .GT. 0 .AND. iLevel .LT. nLevels-1 ) &
          Mask => iMF_Mask % DataPtr( MFI )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCF => MF_uCF(iLevel) % DataPtr( MFI )

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

    CALL amrex_parallel_reduce_sum( BaryonicMass_OffGrid, nLevels )
    CALL amrex_parallel_reduce_sum( Energy_OffGrid      , nLevels )

    DO iLevel = 0, nLevels-1

      BaryonicMass_Change(iLevel) &
        = BaryonicMass_Interior(iLevel) &
            - ( BaryonicMass_Initial(iLevel) + BaryonicMass_OffGrid(iLevel) )

      Energy_Change(iLevel) &
        = Energy_Interior(iLevel) &
            - ( Energy_Initial(iLevel) + Energy_OffGrid(iLevel) )

    END DO

    CALL WriteTally_Euler( Time(0) )

    IF( Verbose )THEN

      CALL DisplayTally( Time(0) )

    END IF

  END SUBROUTINE ComputeTally_Euler_MF


  SUBROUTINE IncrementOffGridTally_Euler_MF( dM )

    REAL(DP), INTENT(in) :: dM(0:nLevels-1,nCF)

    INTEGER :: iLevel

    IF( SuppressTally ) RETURN

    DO iLevel = 0, nLevels-1

      BaryonicMass_OffGrid(iLevel) &
        = BaryonicMass_OffGrid(iLevel) + dM(iLevel,iCF_D)

      Energy_OffGrid(iLevel) &
        = Energy_OffGrid(iLevel) + dM(iLevel,iCF_E)

    END DO

  END SUBROUTINE IncrementOffGridTally_Euler_MF


  SUBROUTINE FinalizeTally_Euler_MF

    IF( SuppressTally ) RETURN

    DEALLOCATE( Energy_Change )
    DEALLOCATE( Energy_OffGrid )
    DEALLOCATE( Energy_Initial )
    DEALLOCATE( Energy_Interior )

    DEALLOCATE( BaryonicMass_Change )
    DEALLOCATE( BaryonicMass_OffGrid )
    DEALLOCATE( BaryonicMass_Initial )
    DEALLOCATE( BaryonicMass_Interior )

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

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1, nDOFX

      IF( nLevels .GT. 1 .AND. iLevel .LT. nLevels-1 )THEN

        IF( Mask(iX1,iX2,iX3,1) .NE. iCoarse_MFM ) CYCLE

      END IF

      d3X = MeshX(1) % Width(iX1) &
              * MeshX(2) % Width(iX2) &
              * MeshX(3) % Width(iX3)

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

    CALL DestroyMesh_MF( MeshX )

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
        SUM( BaryonicMass_Interior ) / UnitsDisplay % MassUnit, &
        SUM( BaryonicMass_OffGrid  ) / UnitsDisplay % MassUnit, &
        SUM( BaryonicMass_Initial  ) / UnitsDisplay % MassUnit, &
        SUM( BaryonicMass_Change   ) / UnitsDisplay % MassUnit

      CLOSE( FileUnit )

      ! --- Energy ---

      OPEN( NEWUNIT = FileUnit, FILE = TRIM( Energy_FileName ), &
            POSITION = 'APPEND', ACTION = 'WRITE' )

      WRITE( FileUnit, '(5(ES25.16E3,1x))' ) &
        Time / UnitsDisplay % TimeUnit, &
        SUM( Energy_Interior ) / UnitsDisplay % EnergyGlobalUnit, &
        SUM( Energy_OffGrid  ) / UnitsDisplay % EnergyGlobalUnit, &
        SUM( Energy_Initial  ) / UnitsDisplay % EnergyGlobalUnit, &
        SUM( Energy_Change   ) / UnitsDisplay % EnergyGlobalUnit

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
        SUM( BaryonicMass_Interior ) / UnitsDisplay % MassUnit, &
        UnitsDisplay % MassLabel
      WRITE(*,'(A6,A40,ES14.7E2,x,A)') &
        '', 'Baryonic Mass Initial..: ', &
        SUM( BaryonicMass_Initial )  / UnitsDisplay % MassUnit, &
        UnitsDisplay % MassLabel
      WRITE(*,'(A6,A40,ES14.7E2,x,A)') &
        '', 'Baryonic Mass Off Grid.: ', &
        SUM( BaryonicMass_OffGrid )  / UnitsDisplay % MassUnit, &
        UnitsDisplay % MassLabel
      WRITE(*,'(A6,A40,ES14.7E2,x,A)') &
        '', 'Baryonic Mass Change...: ', &
        SUM( BaryonicMass_Change )   / UnitsDisplay % MassUnit, &
        UnitsDisplay % MassLabel

      WRITE(*,*)
      WRITE(*,'(A6,A40,ES14.7E2,x,A)') &
        '', 'Energy Interior.: ', &
        SUM( Energy_Interior ) / UnitsDisplay % EnergyGlobalUnit, &
        UnitsDisplay % EnergyGlobalLabel
      WRITE(*,'(A6,A40,ES14.7E2,x,A)') &
        '', 'Energy Initial..: ', &
        SUM( Energy_Initial )  / UnitsDisplay % EnergyGlobalUnit, &
        UnitsDisplay % EnergyGlobalLabel
      WRITE(*,'(A6,A40,ES14.7E2,x,A)') &
        '', 'Energy Off Grid.: ', &
        SUM( Energy_OffGrid )  / UnitsDisplay % EnergyGlobalUnit, &
        UnitsDisplay % EnergyGlobalLabel
      WRITE(*,'(A6,A40,ES14.7E2,x,A)') &
        '', 'Energy Change...: ', &
        SUM( Energy_Change )   / UnitsDisplay % EnergyGlobalUnit, &
        UnitsDisplay % EnergyGlobalLabel

      WRITE(*,*)

    END IF

  END SUBROUTINE DisplayTally


END MODULE MF_Euler_TallyModule
