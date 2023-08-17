MODULE InputOutputModuleAMReX

  USE ISO_C_BINDING

  ! --- AMReX Modules ---

  USE amrex_fort_module, ONLY: &
    amrex_spacedim
  USE amrex_plotfile_module, ONLY: &
    amrex_write_plotfile
  USE amrex_string_module, ONLY: &
    amrex_string, &
    amrex_string_build
  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_boxarray_module, ONLY: &
    amrex_boxarray, &
    amrex_boxarray_build, &
    amrex_boxarray_destroy
  USE amrex_distromap_module,  ONLY: &
    amrex_distromap,       &
    amrex_distromap_build, &
    amrex_distromap_destroy
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_multifab_build, &
    amrex_multifab_destroy, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy, &
    amrex_imultifab
  USE amrex_geometry_module, ONLY: &
    amrex_geometry, &
    amrex_geometry_build, &
    amrex_geometry_destroy
  USE amrex_amrcore_module, ONLY: &
    amrex_get_amrcore, &
    amrex_get_numlevels, &
    amrex_ref_ratio, &
    amrex_set_boxarray, &
    amrex_set_distromap, &
    amrex_set_geometry, &
    amrex_set_finest_level
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor, &
    amrex_parallel_myproc
  USE amrex_fluxregister_module, ONLY: &
    amrex_fluxregister_build
  USE amrex_amr_module, ONLY: &
    amrex_geom

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    nDOFZ, &
    nDOFE, &
    iZ_B0, &
    iZ_E0
  USE ReferenceElementModuleX, ONLY: &
    WeightsX_q, &
    nDOFX_X1
  USE ReferenceElementModuleZ, ONLY: &
    nDOFZ_Z2
  USE ReferenceElementModuleE, ONLY: &
    WeightsE
  USE ReferenceElementModule, ONLY: &
    Weights_q
  USE MeshModule, ONLY: &
    MeshType, &
    MeshE, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    ShortNamesGF, &
    unitsGF, &
    nGF, &
    iGF_SqrtGm
  USE FluidFieldsModule, ONLY: &
    ShortNamesCF, &
    unitsCF, &
    nCF, &
    ShortNamesPF, &
    unitsPF, &
    nPF, &
    ShortNamesAF, &
    unitsAF, &
    nAF, &
    ShortNamesDF, &
    unitsDF, &
    nDF
  USE RadiationFieldsModule, ONLY: &
    ShortNamesCR, &
    unitsCR, &
    nCR, &
    ShortNamesPR, &
    unitsPR, &
    nPR, &
    ShortNamesGR, &
    unitsGR, &
    nGR
  USE UnitsModule, ONLY: &
    UnitsDisplay

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero, &
    Two
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF
  USE MF_FieldsModule_Geometry, ONLY: &
    MF_uGF
  USE MF_FieldsModule_Euler, ONLY: &
    MF_uCF, &
    MF_uPF, &
    MF_uAF, &
    MF_uDF, &
    FluxRegister_Euler
  USE MF_FieldsModule_TwoMoment, ONLY: &
    MF_uCR, &
    MF_uPR, &
    MF_uGR, &
    FluxRegister_TwoMoment
  USE FillPatchModule, ONLY: &
    FillPatch
  USE InputParsingModule, ONLY: &
    nLevels, &
    nMaxLevels, &
    MaxGridSizeX, &
    dt, &
    StepNo, &
    swX, &
    t_new, &
    PlotFileNameRoot, &
    nX, &
    nE, &
    nSpecies, &
    iRestart, &
    UseTiling, &
    UseFluxCorrection_Euler, &
    UseFluxCorrection_TwoMoment, &
    iOS_CPP
  USE MF_Euler_TallyModule, ONLY: &
    BaryonicMass_Initial, &
    BaryonicMass_OffGrid, &
    Energy_Initial, &
    Energy_OffGrid, &
    ElectronNumber_Initial, &
    ElectronNumber_OffGrid, &
    ADMMass_Initial, &
    ADMMass_OffGrid

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ReadCheckpointFile
  PUBLIC :: WriteFieldsAMReX_Checkpoint
  PUBLIC :: WriteFieldsAMReX_PlotFile

  INTERFACE

    SUBROUTINE WriteFieldsAMReX_Checkpoint &
                 ( StepNo, nLevels, dt, time, &
                   BaryonicMassArr, &
                   EnergyArr, &
                   ElectronNumberArr, &
                   ADMMassArr, &
                   pBA, &
                   iWriteFields_uGF, iWriteFields_uCF, iWriteFields_uCR, &
                   pMF_uGF_Option, pMF_uCF_Option, pMF_uCR_Option ) BIND(c)
      IMPORT
      IMPLICIT NONE
      INTEGER(c_int),        INTENT(in) :: StepNo(*)
      INTEGER(c_int), VALUE, INTENT(in) :: nLevels
      REAL(DP)      ,        INTENT(in) :: dt(*), time(*)
      REAL(DP)      ,        INTENT(in) :: BaryonicMassArr(*)
      REAL(DP)      ,        INTENT(in) :: EnergyArr(*)
      REAL(DP)      ,        INTENT(in) :: ElectronNumberArr(*)
      REAL(DP)      ,        INTENT(in) :: ADMMassArr(*)
      TYPE(c_ptr)   ,        INTENT(in) :: pBA(*)
      INTEGER(c_int), VALUE, INTENT(in) :: iWriteFields_uGF
      INTEGER(c_int), VALUE, INTENT(in) :: iWriteFields_uCF
      INTEGER(c_int), VALUE, INTENT(in) :: iWriteFields_uCR
      TYPE(c_ptr)   ,        INTENT(in), OPTIONAL :: pMF_uGF_Option(*)
      TYPE(c_ptr)   ,        INTENT(in), OPTIONAL :: pMF_uCF_Option(*)
      TYPE(c_ptr)   ,        INTENT(in), OPTIONAL :: pMF_uCR_Option(*)
    END SUBROUTINE WriteFieldsAMReX_Checkpoint

    SUBROUTINE ReadHeaderAndBoxArrayData &
                 ( FinestLevelArr, StepNo, dt, Time, &
                   BaryonicMassArr, &
                   EnergyArr, &
                   ElectronNumberArr, &
                   ADMMassArr, &
                   pBA, pDM, iChkFile ) BIND(c)
      IMPORT
      IMPLICIT NONE
      INTEGER(c_int), INTENT(out) :: FinestLevelArr(*)
      INTEGER(c_int), INTENT(out) :: StepNo(*)
      REAL(DP)      , INTENT(out) :: dt(*), Time(*)
      REAL(DP)      , INTENT(out) :: BaryonicMassArr(*)
      REAL(DP)      , INTENT(out) :: EnergyArr(*)
      REAL(DP)      , INTENT(out) :: ElectronNumberArr(*)
      REAL(DP)      , INTENT(out) :: ADMMassArr(*)
      TYPE(c_ptr)   , INTENT(out) :: pBA(*), pDM(*)
      INTEGER(c_int), VALUE       :: iChkFile
    END SUBROUTINE ReadHeaderAndBoxArrayData

    SUBROUTINE ReadMultiFabData( FinestLevel, pMF, iMF, iChkFile ) BIND(c)
      IMPORT
      IMPLICIT NONE
      INTEGER(c_int), VALUE       :: FinestLevel
      TYPE(c_ptr),    INTENT(out) :: pMF(*)
      INTEGER(c_int), VALUE       :: iMF
      INTEGER(c_int), VALUE       :: iChkFile
    END SUBROUTINE ReadMultiFabData

  END INTERFACE

CONTAINS


  SUBROUTINE WriteFieldsAMReX_PlotFile &
    ( Time, StepNo, MF_uGF, &
      MF_uGF_Option, MF_uCF_Option, MF_uPF_Option, &
      MF_uAF_Option, MF_uDF_Option, &
      MF_uCR_Option, MF_uPR_Option, MF_uGR_Option, PlotFileNumber_Option )

    REAL(DP)            , INTENT(in) :: Time
    INTEGER             , INTENT(in) :: StepNo(0:)
    TYPE(amrex_multifab), INTENT(in) :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(in), OPTIONAL :: MF_uGF_Option(0:)
    TYPE(amrex_multifab), INTENT(in), OPTIONAL :: MF_uCF_Option(0:)
    TYPE(amrex_multifab), INTENT(in), OPTIONAL :: MF_uPF_Option(0:)
    TYPE(amrex_multifab), INTENT(in), OPTIONAL :: MF_uAF_Option(0:)
    TYPE(amrex_multifab), INTENT(in), OPTIONAL :: MF_uDF_Option(0:)
    TYPE(amrex_multifab), INTENT(in), OPTIONAL :: MF_uCR_Option(0:)
    TYPE(amrex_multifab), INTENT(in), OPTIONAL :: MF_uPR_Option(0:)
    TYPE(amrex_multifab), INTENT(in), OPTIONAL :: MF_uGR_Option(0:)
    INTEGER             , INTENT(in), OPTIONAL :: PlotFileNumber_Option

    CHARACTER(08)                   :: NumberString
    CHARACTER(64)                   :: PlotFileName
    CHARACTER(32)                   :: ShortNamesCR_Z
    CHARACTER(32)                   :: ShortNamesPR_Z
    CHARACTER(32)                   :: ShortNamesGR_Z
    CHARACTER(3)                    :: iSC, iZ1C
    LOGICAL                         :: WriteGF
    LOGICAL                         :: WriteFF_C
    LOGICAL                         :: WriteFF_P
    LOGICAL                         :: WriteFF_A
    LOGICAL                         :: WriteFF_D
    LOGICAL                         :: WriteRF_C
    LOGICAL                         :: WriteRF_P
    LOGICAL                         :: WriteRF_GR
    INTEGER                         :: iFd, iOS, iLevel, nF, iS, iZ1
    TYPE(amrex_multifab)            :: MF_plt(0:nLevels-1)
    TYPE(amrex_string), ALLOCATABLE :: VarNames(:)

    nF = 7 ! MPI proc, X1_C, X2_C, X3_C, dX1, dX2, dX3

    WriteGF = .FALSE.
    IF( PRESENT( MF_uGF_Option ) )THEN

      WriteGF = .TRUE.
      nF = nF + nGF

    END IF

    WriteFF_C = .FALSE.
    IF( PRESENT( MF_uCF_Option ) )THEN

      WriteFF_C = .TRUE.
      nF = nF + nCF

    END IF

    WriteFF_P = .FALSE.
    IF( PRESENT( MF_uPF_Option ) )THEN

      WriteFF_P = .TRUE.
      nF = nF + nPF

    END IF

    WriteFF_A = .FALSE.
    IF( PRESENT( MF_uAF_Option ) )THEN

      WriteFF_A = .TRUE.
      nF = nF + nAF

    END IF

    WriteFF_D = .FALSE.
    IF( PRESENT( MF_uDF_Option ) )THEN

      WriteFF_D = .TRUE.
      nF = nF + nDF

    END IF

    WriteRF_C = .FALSE.
    IF( PRESENT( MF_uCR_Option ) )THEN

      WriteRF_C = .TRUE.
      nF = nF + nCR * nE * nSpecies

    END IF

    WriteRF_P = .FALSE.
    IF( PRESENT( MF_uPR_Option ) )THEN

      WriteRF_P = .TRUE.
      nF = nF + nPR * nE * nSpecies

    END IF

    WriteRF_GR = .FALSE.
    IF( PRESENT( MF_uGR_Option ) )THEN

      WriteRF_GR = .TRUE.
      nF = nF + nGR *  nSpecies

    END IF

    IF( PRESENT( PlotFileNumber_Option ) )THEN

      WRITE(NumberString,'(I8.8)') PlotFileNumber_Option

    ELSE

      WRITE(NumberString,'(I8.8)') StepNo(0)

    END IF

    PlotFileName = TRIM( PlotFileNameRoot ) // NumberString

    IF( amrex_parallel_ioprocessor() )THEN

      WRITE(*,*)
      WRITE(*,'(6x,A,A)') 'Writing PlotFile ', PlotFileName

    END IF

    ALLOCATE( VarNames(nF) )

    CALL amrex_string_build( VarNames( 1 ), 'MPIProcess' )
    CALL amrex_string_build( VarNames( 2 ), 'X1_C' )
    CALL amrex_string_build( VarNames( 3 ), 'X2_C' )
    CALL amrex_string_build( VarNames( 4 ), 'X3_C' )
    CALL amrex_string_build( VarNames( 5 ), 'dX1' )
    CALL amrex_string_build( VarNames( 6 ), 'dX2' )
    CALL amrex_string_build( VarNames( 7 ), 'dX3' )

    iOS = 7

    IF( WriteGF )THEN

      DO iFd = 1, nGF

        CALL amrex_string_build &
               ( VarNames( iFd + iOS ), TRIM( ShortNamesGF(iFd) ) )

      END DO

      iOS = iOS + nGF

    END IF

    IF( WriteFF_C )THEN

      DO iFd = 1, nCF

        CALL amrex_string_build &
               ( VarNames( iFd + iOS ), TRIM( ShortNamesCF(iFd) ) )

      END DO

      iOS = iOS + nCF

    END IF

    IF( WriteFF_P )THEN

      DO iFd = 1, nPF

        CALL amrex_string_build &
               ( VarNames( iFd + iOS ), TRIM( ShortNamesPF(iFd) ) )

      END DO

      iOS = iOS + nPF

    END IF

    IF( WriteFF_A )THEN

      DO iFd = 1, nAF

        CALL amrex_string_build &
               ( VarNames( iFd + iOS ), TRIM( ShortNamesAF(iFd) ) )

      END DO

      iOS = iOS + nAF

    END IF

    IF( WriteFF_D )THEN

      DO iFd = 1, nDF

        CALL amrex_string_build &
               ( VarNames( iFd + iOS ), TRIM( ShortNamesDF(iFd) ) )

      END DO

      iOS = iOS + nDF

    END IF

    IF( WriteRF_C )THEN

      DO iS  = 1       , nSpecies
      DO iZ1 = iZ_B0(1), iZ_E0(1)
      DO iFd = 1       , nCR

        WRITE(iZ1C,'(I3.3)') iZ1
        WRITE(iSC ,'(I3.3)') iS

        ShortNamesCR_Z = TRIM( ShortNamesCR(iFd) ) // '_' // iZ1C // '_' // iSC

        CALL amrex_string_build &
               ( VarNames( iFd + ( iZ1 - 1 ) * nCR &
                   + ( iS - 1 ) * nE * nCR + iOS ), TRIM( ShortNamesCR_Z ) )

      END DO
      END DO
      END DO

      iOS = iOS + nCR * nE * nSpecies

    END IF

    IF( WriteRF_P )THEN

      DO iS  = 1       , nSpecies
      DO iZ1 = iZ_B0(1), iZ_E0(1)
      DO iFd = 1       , nPR

        WRITE(iZ1C,'(I3.3)') iZ1
        WRITE(iSC ,'(I3.3)') iS

        ShortNamesPR_Z = TRIM( ShortNamesPR(iFd) ) // '_' // iZ1C // '_' // iSC

        CALL amrex_string_build &
               ( VarNames( iFd + ( iZ1 - 1 ) * nPR &
                   + ( iS - 1 ) * nE * nPR + iOS ), TRIM( ShortNamesPR_Z ) )

      END DO
      END DO
      END DO

      iOS = iOS + nPR * nE * nSpecies

    END IF


    IF( WriteRF_GR )THEN

      DO iS  = 1       , nSpecies
      DO iFd = 1, nGR

        WRITE(iSC ,'(I3.3)') iS

        ShortNamesGR_Z = TRIM( ShortNamesGR(iFd) ) // '_' // iSC

        CALL amrex_string_build &
               ( VarNames( iFd &
                   + ( iS - 1 ) * nGR + iOS ), TRIM( ShortNamesGR_Z ) )

      END DO
      END DO

      iOS = iOS + nGR * nSpecies

    END IF



    DO iLevel = 0, nLevels-1

      CALL amrex_multifab_build &
             ( MF_plt(iLevel), MF_uGF(iLevel) % BA, &
                               MF_uGF(iLevel) % DM, &
               nF, 0 )
      CALL MF_plt(iLevel) % setVal( Zero )

      CALL WriteMPI( MF_uGF(iLevel), MF_plt(iLevel) )

      CALL WriteMesh( iLevel, MF_uGF(iLevel), MF_plt(iLevel) )

      iOS = 7

      IF( WriteGF )THEN

        CALL ComputeCellAverage_X_MF &
               ( nGF, MF_uGF(iLevel), MF_uGF_Option(iLevel), &
                 iOS, 'GF', MF_plt(iLevel) )

        iOS = iOS + nGF

      END IF

      IF( WriteFF_C )THEN

        CALL ComputeCellAverage_X_MF &
               ( nCF, MF_uGF(iLevel), MF_uCF_Option(iLevel), &
                 iOS, 'CF', MF_plt(iLevel) )

        iOS = iOS + nCF

      END IF

      IF( WriteFF_P )THEN

        CALL ComputeCellAverage_X_MF &
               ( nPF, MF_uGF(iLevel), MF_uPF_Option(iLevel), &
                 iOS, 'PF', MF_plt(iLevel) )

        iOS = iOS + nPF

      END IF

      IF( WriteFF_A )THEN

        CALL ComputeCellAverage_X_MF &
               ( nAF, MF_uGF(iLevel), MF_uAF_Option(iLevel), &
                 iOS, 'AF', MF_plt(iLevel) )

        iOS = iOS + nAF

      END IF

      IF( WriteFF_D )THEN

        CALL ComputeCellAverage_X_MF &
               ( nDF, MF_uGF(iLevel), MF_uDF_Option(iLevel), &
                 iOS, 'DF', MF_plt(iLevel) )

        iOS = iOS + nDF

      END IF

      IF( WriteRF_C )THEN

        CALL ComputeCellAverage_Z_MF &
               ( nCR, MF_uGF(iLevel), MF_uCR_Option(iLevel), &
                 iOS, 'CR', MF_plt(iLevel) )

        iOS = iOS + nCR * nE * nSpecies

      END IF

      IF( WriteRF_P )THEN

        CALL ComputeCellAverage_Z_MF &
               ( nPR, MF_uGF(iLevel), MF_uPR_Option(iLevel), &
                 iOS, 'PR', MF_plt(iLevel) )

        iOS = iOS + nPR * nE * nSpecies

      END IF

      IF( WriteRF_GR )THEN

        CALL ComputeCellAverage_Integral_MF &
               ( nGR, MF_uGF(iLevel), MF_uGR_Option(iLevel), &
                 iOS, 'GR', MF_plt(iLevel) )

        iOS = iOS + nGR * nSpecies

      END IF
    END DO ! iLevel = 0, nLevels-1

    CALL amrex_write_plotfile &
           ( PlotFileName, nLevels, MF_plt, VarNames, &
             amrex_geom, Time / UnitsDisplay % TimeUnit, StepNo, &
             amrex_ref_ratio )

    DO iLevel = 0, nLevels-1

      CALL amrex_multifab_destroy ( MF_plt(iLevel) )

    END DO

    DEALLOCATE( VarNames )

  END SUBROUTINE WriteFieldsAMReX_PlotFile


  SUBROUTINE ReadCheckpointFile &
    ( ReadFields_uCF_Option, ReadFields_uCR_Option )

    LOGICAL, INTENT(in), OPTIONAL :: ReadFields_uCF_Option
    LOGICAL, INTENT(in), OPTIONAL :: ReadFields_uCR_Option

    INTEGER     :: iLevel, FinestLevel
    TYPE(c_ptr) :: pBA(0:nMaxLevels-1)
    TYPE(c_ptr) :: pDM(0:nMaxLevels-1)
    TYPE(c_ptr) :: pGF(0:nMaxLevels-1)
    TYPE(c_ptr) :: pCF(0:nMaxLevels-1)
    TYPE(c_ptr) :: pCR(0:nMaxLevels-1)
    TYPE(c_ptr) :: amrcore

    TYPE(amrex_box)       :: BX
    TYPE(amrex_distromap) :: DM  (0:nMaxLevels-1)
    TYPE(amrex_boxarray)  :: BA  (0:nMaxLevels-1)
    TYPE(amrex_geometry)  :: GEOM(0:nMaxLevels-1)

    INTEGER :: nXX(3)
    INTEGER :: FinestLevelArr(0:0) ! Hack

    REAL(DP) :: BaryonicMassArr  (0:1)
    REAL(DP) :: EnergyArr        (0:1)
    REAL(DP) :: ElectronNumberArr(0:1)
    REAL(DP) :: ADMMassArr       (0:1)

    LOGICAL :: ReadFields_uCF
    LOGICAL :: ReadFields_uCR

    ReadFields_uCF = .FALSE.
    IF( PRESENT( ReadFields_uCF_Option ) ) &
      ReadFields_uCF = ReadFields_uCF_Option

    ReadFields_uCR = .FALSE.
    IF( PRESENT( ReadFields_uCR_Option ) ) &
      ReadFields_uCR = ReadFields_uCR_Option

    amrcore = amrex_get_amrcore()

    DO iLevel = 0, nMaxLevels-1

      nXX = nX

      nXX(1) = 2**( iLevel ) * nX(1)
      IF( amrex_spacedim .GT. 1 ) nXX(2) = 2**( iLevel ) * nX(2)
      IF( amrex_spacedim .GT. 2 ) nXX(3) = 2**( iLevel ) * nX(3)

      BX = amrex_box( 1 - iOS_CPP, nXX - iOS_CPP )

      CALL amrex_boxarray_build( BA(iLevel), BX )

      CALL BA(iLevel) % maxSize( MaxGridSizeX )

      CALL amrex_geometry_build( GEOM(iLevel), BX )

      CALL amrex_distromap_build( DM(iLevel), BA(iLevel) )

    END DO

    pBA(0:nMaxLevels-1) = BA(0:nMaxLevels-1) % P
    pDM(0:nMaxLevels-1) = DM(0:nMaxLevels-1) % P

    CALL ReadHeaderAndBoxArrayData &
           ( FinestLevelArr, StepNo, dt, t_new, &
             BaryonicMassArr, &
             EnergyArr, &
             ElectronNumberArr, &
             ADMMassArr, &
             pBA, pDM, iRestart )

    FinestLevel = FinestLevelArr(0) ! Hack
    CALL amrex_set_finest_level( FinestLevel )
    nLevels = amrex_get_numlevels()

    BaryonicMass_Initial   = BaryonicMassArr  (0)
    BaryonicMass_OffGrid   = BaryonicMassArr  (1)
    Energy_Initial         = EnergyArr        (0)
    Energy_OffGrid         = EnergyArr        (1)
    ElectronNumber_Initial = ElectronNumberArr(0)
    ElectronNumber_OffGrid = ElectronNumberArr(1)
    ADMMass_Initial        = ADMMassArr       (0)
    ADMMass_OffGrid        = ADMMassArr       (1)

    DO iLevel = 0, nLevels-1

      BA(iLevel) = pBA(iLevel)
      DM(iLevel) = pDM(iLevel)

      CALL amrex_set_boxarray ( iLevel, BA  (iLevel) )
      CALL amrex_set_distromap( iLevel, DM  (iLevel) )
      CALL amrex_set_geometry ( iLevel, GEOM(iLevel) )

      CALL amrex_multifab_build &
             ( MF_uGF(iLevel), BA(iLevel), DM(iLevel), nDOFX * nGF, swX )
      CALL MF_uGF(iLevel) % SetVal( Zero )

      IF( ReadFields_uCF )THEN

        CALL amrex_multifab_build &
               ( MF_uCF(iLevel), BA(iLevel), DM(iLevel), nDOFX * nCF, swX )
        CALL MF_uCF(iLevel) % SetVal( Zero )

        CALL amrex_multifab_build &
               ( MF_uPF(iLevel), BA(iLevel), DM(iLevel), nDOFX * nPF, swX )
        CALL MF_uPF(iLevel) % SetVal( Zero )

        CALL amrex_multifab_build &
               ( MF_uAF(iLevel), BA(iLevel), DM(iLevel), nDOFX * nAF, swX )
        CALL MF_uAF(iLevel) % SetVal( Zero )

        CALL amrex_multifab_build &
               ( MF_uDF(iLevel), BA(iLevel), DM(iLevel), nDOFX * nDF, swX )
        CALL MF_uDF(iLevel) % SetVal( Zero )

        ! Assume nDOFX_X2 = nDOFX_X3 = nDOFX_X1
        IF( iLevel .GT. 0 .AND. UseFluxCorrection_Euler )THEN

          CALL amrex_fluxregister_build &
                 ( FluxRegister_Euler(iLevel), BA(iLevel), DM(iLevel), &
                   amrex_ref_ratio(iLevel-1), iLevel, nDOFX_X1*nCF )

        END IF

      END IF

      IF( ReadFields_uCR )THEN

        CALL amrex_multifab_build &
               ( MF_uCR(iLevel), BA(iLevel), DM(iLevel), &
                 nDOFZ * nCR * nE * nSpecies, swX )
        CALL MF_uCR(iLevel) % SetVal( Zero )

        CALL amrex_multifab_build &
               ( MF_uPR(iLevel), BA(iLevel), DM(iLevel), &
                 nDOFZ * nPR * nE * nSpecies, swX )
        CALL MF_uPR(iLevel) % SetVal( Zero )


        CALL amrex_multifab_build &
               ( MF_uGR(iLevel), BA(iLevel), DM(iLevel), &
                 nDOFX * nGR * nSpecies, swX )
        CALL MF_uGR(iLevel) % SetVal( Zero )


        ! Assume nDOFZ_Z3 = nDOFZ_Z4 = nDOFZ_Z2
        IF( iLevel .GT. 0 .AND. UseFluxCorrection_TwoMoment )THEN

          CALL amrex_fluxregister_build &
                 ( FluxRegister_TwoMoment(iLevel), BA(iLevel), DM(iLevel), &
                   amrex_ref_ratio(iLevel-1), iLevel, &
                   nDOFZ_Z2 * nCR * nE * nSpecies )

        END IF

      END IF

    END DO

    pGF(0:nLevels-1) = MF_uGF(0:nLevels-1) % P
    CALL ReadMultiFabData( FinestLevel, pGF, 0, iRestart )

    IF( ReadFields_uCF )THEN

      pCF(0:nLevels-1) = MF_uCF(0:nLevels-1) % P
      CALL ReadMultiFabData( FinestLevel, pCF, 1, iRestart )

    END IF

    IF( ReadFields_uCR )THEN

      pCR(0:nLevels-1) = MF_uCR(0:nLevels-1) % P
      CALL ReadMultiFabData( FinestLevel, pCR, 2, iRestart )

    END IF

  END SUBROUTINE ReadCheckpointFile


  ! --- PRIVATE SUBROUTINES ---


  SUBROUTINE ComputeCellAverage_X_MF &
    ( nFd, MF_uGF, MF, iOS, Field, MF_plt )

    INTEGER              , INTENT(in)    :: nFd, iOS
    TYPE(amrex_multifab) , INTENT(in)    :: MF_uGF
    TYPE(amrex_multifab) , INTENT(in)    :: MF
    CHARACTER(2)         , INTENT(in)    :: Field
    TYPE(amrex_multifab) , INTENT(inout) :: MF_plt

    INTEGER                       :: iX1, iX2, iX3, iFd
    INTEGER                       :: iX_B0(3), iX_E0(3)
    INTEGER                       :: lo_G(4), hi_G(4)
    INTEGER                       :: lo_U(4), hi_U(4)
    REAL(DP)                      :: G_K(nDOFX,nGF)
    REAL(DP)                      :: U_K(nDOFX,nFd)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(DP), CONTIGUOUS, POINTER :: G    (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: U    (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: U_plt(:,:,:,:)

    CALL amrex_mfiter_build( MFI, MF_uGF, tiling = UseTiling )

    DO WHILE( MFI % next() )

      G     => MF_uGF   % DataPtr( MFI )
      U     => MF       % DataPtr( MFI )
      U_plt => MF_plt   % DataPtr( MFI )

      BX = MFI % TileBox()

      iX_B0 = BX % lo
      iX_E0 = BX % hi

      lo_G = LBOUND( G ); hi_G = UBOUND( G )
      lo_U = LBOUND( U ); hi_U = UBOUND( U )
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)

        G_K(1:nDOFX,1:nGF) &
          = RESHAPE( G(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )

        U_K(1:nDOFX,1:nFd) &
          = RESHAPE( U(iX1,iX2,iX3,lo_U(4):hi_U(4)), [ nDOFX, nFd ] )

        DO iFd = 1, nFd

          U_plt(iX1,iX2,iX3,iFd+iOS) &
            = SUM( WeightsX_q * U_K(:,iFd) * G_K(:,iGF_SqrtGm) ) &
                / SUM( WeightsX_q * G_K(:,iGF_SqrtGm) )

        END DO

      END DO
      END DO
      END DO

      CALL ConvertUnits( Field, nFd, iOS, U_plt )

    END DO

    CALL amrex_mfiter_destroy( MFI )

  END SUBROUTINE ComputeCellAverage_X_MF




  SUBROUTINE ComputeCellAverage_Z_MF &
    ( nFd, MF_uGF, MF, iOS, Field, MF_plt )

    INTEGER              , INTENT(in)    :: nFd, iOS
    TYPE(amrex_multifab) , INTENT(in)    :: MF_uGF
    TYPE(amrex_multifab) , INTENT(in)    :: MF
    CHARACTER(2)         , INTENT(in)    :: Field
    TYPE(amrex_multifab) , INTENT(inout) :: MF_plt

    INTEGER                       :: iX1, iX2, iX3, iFd, iS, iZ1, iNZ, iNE, iNX
    INTEGER                       :: iX_B0(3), iX_E0(3)
    INTEGER                       :: lo_G(4), hi_G(4)
    INTEGER                       :: lo_U(4), hi_U(4)
    REAL(DP)                      :: G_K(nDOFX,nGF)
    REAL(DP)                      :: U_K(nDOFZ,nE,nFd,nSpecies)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(DP), CONTIGUOUS, POINTER :: G    (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: U    (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: U_plt(:,:,:,:)
    REAL(DP)                      :: Eq(1:nDOFE), E(1:nDOFZ), SqrtGM(1:nDOFZ), V_K, SUM2

    CALL amrex_mfiter_build( MFI, MF_uGF, tiling = UseTiling )

    DO WHILE( MFI % next() )

      G     => MF_uGF   % DataPtr( MFI )
      U     => MF       % DataPtr( MFI )
      U_plt => MF_plt   % DataPtr( MFI )

      BX = MFI % tilebox()

      iX_B0 = BX % lo
      iX_E0 = BX % hi

      lo_G = LBOUND( G ); hi_G = UBOUND( G )
      lo_U = LBOUND( U ); hi_U = UBOUND( U )

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)

        G_K(1:nDOFX,1:nGF) &
          = RESHAPE( G(iX1,iX2,iX3,lo_G(4):hi_G(4)), &
                     [ nDOFX, nGF ] )

        U_K(1:nDOFZ,1:nE,1:nFd,1:nSpecies) &
          = RESHAPE( U(iX1,iX2,iX3,lo_U(4):hi_U(4)), &
                     [ nDOFZ, nE, nFd, nSpecies ] )

        DO iS  = 1, nSpecies
        DO iZ1 = 1, nE

          ! --- Compute cell-average ---

            V_K = 0.0_DP

          DO iNZ = 1, nDOFZ

            iNE = MOD( iNZ - 1, nDOFE ) + 1
            iNX = MOD( (iNZ-1) / nDOFE, nDOFX ) + 1

            E (iNZ) = NodeCoordinate( MeshE, iZ1, iNE )
            Eq(iNE) = NodeCoordinate( MeshE, iZ1, iNE )

            SqrtGM(iNZ) = G_K(iNX,iGF_SqrtGm)

            V_K = V_K + Weights_q(iNZ) * Eq(iNE)**2 * G_K(iNX,iGF_SqrtGm)

          END DO


          DO iFd = 1, nFd


            SUM2 = 0.0_DP

            DO iNZ = 1, nDOFZ

              iNE = MOD( iNZ - 1, nDOFE ) + 1
              iNX = MOD( (iNZ-1) / nDOFE, nDOFX ) + 1

              E (iNZ) = NodeCoordinate( MeshE, iZ1, iNE )
              Eq(iNE) = NodeCoordinate( MeshE, iZ1, iNE )

              SUM2 = SUM2 &
                   + Weights_q(iNZ) * Eq(iNE)**2 * G_K(iNX,iGF_SqrtGm) * U_K(iNZ,iZ1,iFd,iS)

            END DO

            U_plt(iX1,iX2,iX3, &
                  iFd + ( iZ1 - 1 ) * nFd + ( iS - 1 ) * nFd * nE + iOS) &
              = SUM2 / V_K

          END DO

        END DO
        END DO

      END DO
      END DO
      END DO

      CALL ConvertUnits( Field, nFd, iOS, U_plt )

    END DO

    CALL amrex_mfiter_destroy( MFI )

  END SUBROUTINE ComputeCellAverage_Z_MF


  SUBROUTINE ComputeCellAverage_Integral_MF &
    ( nFd, MF_uGF, MF, iOS, Field, MF_plt )

    INTEGER,              INTENT(in)    :: nFd, iOS
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF
    TYPE(amrex_multifab), INTENT(in)    :: MF
    CHARACTER(2),         INTENT(in)    :: Field
    TYPE(amrex_multifab), INTENT(inout) :: MF_plt

    INTEGER                       :: iX1, iX2, iX3, iFd, iS
    INTEGER                       :: lo_G(4), hi_G(4)
    INTEGER                       :: lo_U(4), hi_U(4)
    REAL(DP)                      :: G_K(nDOFX,nGF)
    REAL(DP)                      :: U_K(nDOFX,nFd,nSpecies)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(DP), CONTIGUOUS, POINTER :: G    (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: U    (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: U_plt(:,:,:,:)

    CALL amrex_mfiter_build( MFI, MF_uGF, tiling = UseTiling )

    DO WHILE( MFI % next() )

      G     => MF_uGF % DataPtr( MFI )
      U     => MF     % DataPtr( MFI )
      U_plt => MF_plt % DataPtr( MFI )

      BX = MFI % TileBox()

      lo_G = LBOUND( G ); hi_G = UBOUND( G )
      lo_U = LBOUND( U ); hi_U = UBOUND( U )

      DO iS = 1, nSpecies
      DO iX3 = BX % lo(3), BX % hi(3)
      DO iX2 = BX % lo(2), BX % hi(2)
      DO iX1 = BX % lo(1), BX % hi(1)

        G_K(1:nDOFX,1:nGF) &
          = RESHAPE( G(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )

        U_K(1:nDOFX,1:nFd,1:nSpecies) &
          = RESHAPE( U(iX1,iX2,iX3,lo_U(4):hi_U(4)), [ nDOFX, nFd, nSpecies ] )

        DO iFd = 1, nFd

          U_plt(iX1,iX2,iX3,iFd + ( iS - 1 ) * nFd + iOS) &
            = SUM( WeightsX_q * U_K(:,iFd,iS) * G_K(:,iGF_SqrtGm) ) &
                / SUM( WeightsX_q * G_K(:,iGF_SqrtGm) )

        END DO

      END DO
      END DO
      END DO
      END DO

      CALL ConvertUnits( Field, nFd, iOS, U_plt )

    END DO

    CALL amrex_mfiter_destroy( MFI )

  END SUBROUTINE ComputeCellAverage_Integral_MF

  SUBROUTINE ConvertUnits( Field, nFd, iOS, U_plt )

    CHARACTER(2), INTENT(in)    :: Field
    INTEGER     , INTENT(in)    :: nFd, iOS
    REAL(DP)    , INTENT(inout) :: U_plt(:,:,:,:)

    INTEGER :: iFd, iS

    SELECT CASE( Field )

      CASE( 'GF' )

        DO iFd = 1, nFd

          U_plt(:,:,:,iFd+iOS) = U_plt(:,:,:,iFd+iOS) / unitsGF(iFd)

        END DO

      CASE( 'CF' )

        DO iFd = 1, nFd

          U_plt(:,:,:,iFd+iOS) = U_plt(:,:,:,iFd+iOS) / unitsCF(iFd)

        END DO

      CASE( 'PF' )

        DO iFd = 1, nFd

          U_plt(:,:,:,iFd+iOS) = U_plt(:,:,:,iFd+iOS) / unitsPF(iFd)

        END DO

      CASE( 'AF' )

        DO iFd = 1, nFd

          U_plt(:,:,:,iFd+iOS) = U_plt(:,:,:,iFd+iOS) / unitsAF(iFd)

        END DO

      CASE( 'DF' )

        DO iFd = 1, nFd

          U_plt(:,:,:,iFd+iOS) = U_plt(:,:,:,iFd+iOS) / unitsDF(iFd)

        END DO

      CASE( 'CR' )

        DO iFd = 1, nFd

          U_plt(:,:,:,iFd+iOS) = U_plt(:,:,:,iFd+iOS) / unitsCR(iFd)

        END DO

      CASE( 'PR' )

        DO iFd = 1, nFd

          U_plt(:,:,:,iFd+iOS) = U_plt(:,:,:,iFd+iOS) / unitsPR(iFd)

        END DO

      CASE( 'GR' )

        DO iS = 1, nSpecies
        DO iFd = 1, nFd

          U_plt(:,:,:,iFd+(iS-1)*nFd+iOS) &
            = U_plt(:,:,:,iFd+(iS-1)*nFd+iOS) / unitsGR(iFd)

        END DO
        END DO

      CASE DEFAULT

        RETURN

    END SELECT

  END SUBROUTINE ConvertUnits


  SUBROUTINE WriteMPI( MF_uGF, MF_plt )

    TYPE(amrex_multifab) , INTENT(in)    :: MF_uGF
    TYPE(amrex_multifab) , INTENT(inout) :: MF_plt

    INTEGER            :: iX1, iX2, iX3
    INTEGER            :: iX_B0(3), iX_E0(3)
    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI

    REAL(DP), CONTIGUOUS, POINTER :: U_plt(:,:,:,:)

    CALL amrex_mfiter_build( MFI, MF_uGF, tiling = UseTiling )

    DO WHILE( MFI % next() )

      U_plt => MF_plt   % DataPtr( MFI )

      BX = MFI % TileBox()

      iX_B0 = BX % lo
      iX_E0 = BX % hi

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)

        U_plt(iX1,iX2,iX3,1) = amrex_parallel_myproc()

      END DO
      END DO
      END DO

    END DO

    CALL amrex_mfiter_destroy( MFI )

  END SUBROUTINE WriteMPI


  SUBROUTINE WriteMesh( iLevel, MF_uGF, MF_plt )

    INTEGER              , INTENT(in)    :: iLevel
    TYPE(amrex_multifab) , INTENT(in)    :: MF_uGF
    TYPE(amrex_multifab) , INTENT(inout) :: MF_plt

    INTEGER            :: iX1, iX2, iX3
    INTEGER            :: iX_B0(3), iX_E0(3)
    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI

    TYPE(MeshType) :: MeshX(3)

    REAL(DP), CONTIGUOUS, POINTER :: U_plt(:,:,:,:)

    CALL CreateMesh_MF( iLevel, MeshX )

    CALL amrex_mfiter_build( MFI, MF_uGF, tiling = UseTiling )

    ASSOCIATE( U => UnitsDisplay )

    DO WHILE( MFI % next() )

      U_plt => MF_plt   % DataPtr( MFI )

      BX = MFI % TileBox()

      iX_B0 = BX % lo
      iX_E0 = BX % hi

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)

        U_plt(iX1,iX2,iX3,2) = MeshX(1) % Center(iX1) / U % LengthX1Unit
        U_plt(iX1,iX2,iX3,3) = MeshX(2) % Center(iX2) / U % LengthX2Unit
        U_plt(iX1,iX2,iX3,4) = MeshX(3) % Center(iX3) / U % LengthX3Unit

        U_plt(iX1,iX2,iX3,5) = MeshX(1) % Width(iX1) / U % LengthX1Unit
        U_plt(iX1,iX2,iX3,6) = MeshX(2) % Width(iX2) / U % LengthX2Unit
        U_plt(iX1,iX2,iX3,7) = MeshX(3) % Width(iX3) / U % LengthX3Unit

      END DO
      END DO
      END DO

    END DO

    END ASSOCIATE ! U

    CALL amrex_mfiter_destroy( MFI )

    CALL DestroyMesh_MF( MeshX )

  END SUBROUTINE WriteMesh


END MODULE InputOutputModuleAMReX
