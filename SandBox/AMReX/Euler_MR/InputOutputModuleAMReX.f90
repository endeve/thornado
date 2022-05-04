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
    amrex_mfiter_destroy
  USE amrex_geometry_module, ONLY: &
    amrex_geometry, &
    amrex_geometry_build, &
    amrex_geometry_destroy
  USE amrex_amrcore_module, ONLY: &
    amrex_get_amrcore, &
    amrex_ref_ratio
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor, &
    amrex_parallel_myproc
  USE amrex_fluxregister_module, ONLY: &
    amrex_fluxregister_build
  USE amrex_amr_module, ONLY: &
    amrex_geom

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nDOFX
  USE ReferenceElementModuleX, ONLY: &
    WeightsX_q, &
    nDOFX_X1
  USE MeshModule, ONLY: &
    MeshX
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
  USE UnitsModule, ONLY: &
    UnitsDisplay

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF
  USE FillPatchModule, ONLY: &
    FillPatch
  USE InputParsingModule, ONLY: &
    nLevels, &
    MaxGridSizeX, &
    dt, &
    StepNo, &
    swX, &
    t_new, &
    PlotFileBaseName, &
    nX, &
    iRestart, &
    UseTiling, &
    do_reflux

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ReadCheckpointFile
  PUBLIC :: WriteFieldsAMReX_Checkpoint
  PUBLIC :: WriteFieldsAMReX_PlotFile

  INTERFACE

    SUBROUTINE WriteFieldsAMReX_Checkpoint &
                 ( StepNo, nLevels, dt, time, pBA, &
                   pMF_uGF, pMF_uCF ) BIND(c)
       IMPORT
       IMPLICIT NONE
       INTEGER(c_int), INTENT(in) :: StepNo(*)
       INTEGER(c_int), VALUE      :: nLevels
       REAL(DP),       INTENT(in) :: dt(*), time(*)
       TYPE(c_ptr),    INTENT(in) :: pBA(*)
       TYPE(c_ptr),    INTENT(in) :: pMF_uGF(*)
       TYPE(c_ptr),    INTENT(in) :: pMF_uCF(*)
    END SUBROUTINE WriteFieldsAMReX_Checkpoint

    SUBROUTINE ReadHeaderAndBoxArrayData &
                 ( FinestLevel, StepNo, dt, time, &
                   pBA, pDM, iChkFile ) BIND(c)
      IMPORT
      IMPLICIT NONE
      INTEGER(c_int), INTENT(out) :: FinestLevel
      INTEGER(c_int), INTENT(out) :: StepNo(*)
      REAL(DP),       INTENT(out) :: dt(*), time(*)
      TYPE(c_ptr),    INTENT(out) :: pBA(*), pDM(*)
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

    SUBROUTINE amrex_fi_set_boxarray( iLevel, pBA, amrcore ) BIND(c)
      IMPORT
      IMPLICIT NONE
      TYPE(c_ptr),    VALUE :: pBA
      INTEGER(c_int), VALUE :: iLevel
      TYPE(c_ptr),    VALUE :: amrcore
    END SUBROUTINE amrex_fi_set_boxarray

    SUBROUTINE amrex_fi_set_distromap  (lev, pdm, amrcore) BIND(c)
      IMPORT
      IMPLICIT NONE
      TYPE(c_ptr),    VALUE :: pdm
      INTEGER(c_int), VALUE :: lev
      TYPE(c_ptr),    VALUE :: amrcore
    END SUBROUTINE amrex_fi_set_distromap

    SUBROUTINE amrex_fi_clone_boxarray (bao, bai) BIND(c)
      IMPORT
      IMPLICIT NONE
      TYPE(c_ptr)        :: bao
      TYPE(c_ptr), VALUE :: bai
    END SUBROUTINE amrex_fi_clone_boxarray

    SUBROUTINE amrex_fi_set_finest_level (lev, amrcore) BIND(c)
      IMPORT
      IMPLICIT NONE
      INTEGER(c_int), VALUE :: lev
      TYPE(c_ptr),    VALUE :: amrcore
    END SUBROUTINE amrex_fi_set_finest_level

  END INTERFACE

CONTAINS


  SUBROUTINE WriteFieldsAMReX_PlotFile &
    ( Time, StepNo, MF_uGF, &
      MF_uGF_Option, MF_uCF_Option, MF_uPF_Option, &
      MF_uAF_Option, MF_uDF_Option )

    REAL(DP),             INTENT(in) :: Time
    INTEGER,              INTENT(in) :: StepNo(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in) :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in), OPTIONAL :: &
      MF_uGF_Option(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in), OPTIONAL :: &
      MF_uCF_Option(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in), OPTIONAL :: &
      MF_uPF_Option(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in), OPTIONAL :: &
      MF_uAF_Option(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in), OPTIONAL :: &
      MF_uDF_Option(0:nLevels-1)

    CHARACTER(08)                   :: NumberString
    CHARACTER(64)                   :: PlotFileName
    LOGICAL                         :: WriteGF
    LOGICAL                         :: WriteFF_C
    LOGICAL                         :: WriteFF_P
    LOGICAL                         :: WriteFF_A
    LOGICAL                         :: WriteFF_D
    INTEGER                         :: iFd, iOS, iLevel, nF
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

    WRITE(NumberString,'(I8.8)') StepNo(0)

    PlotFileName = TRIM( PlotFileBaseName ) // '_' // NumberString

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

        CALL ComputeCellAverage_MF &
               ( nGF, MF_uGF(iLevel), MF_uGF_Option(iLevel), &
                 iOS, 'GF', MF_plt(iLevel) )

        iOS = iOS + nGF

      END IF

      IF( WriteFF_C )THEN

        CALL ComputeCellAverage_MF &
               ( nCF, MF_uGF(iLevel), MF_uCF_Option(iLevel), &
                 iOS, 'CF', MF_plt(iLevel) )

        iOS = iOS + nCF

      END IF

      IF( WriteFF_P )THEN

        CALL ComputeCellAverage_MF &
               ( nPF, MF_uGF(iLevel), MF_uPF_Option(iLevel), &
                 iOS, 'PF', MF_plt(iLevel) )

        iOS = iOS + nPF

      END IF

      IF( WriteFF_A )THEN

        CALL ComputeCellAverage_MF &
               ( nAF, MF_uGF(iLevel), MF_uAF_Option(iLevel), &
                 iOS, 'AF', MF_plt(iLevel) )

        iOS = iOS + nAF

      END IF

      IF( WriteFF_D )THEN

        CALL ComputeCellAverage_MF &
               ( nDF, MF_uGF(iLevel), MF_uDF_Option(iLevel), &
                 iOS, 'DF', MF_plt(iLevel) )

        iOS = iOS + nDF

      END IF

    END DO ! End of loop over levels

    CALL amrex_write_plotfile &
           ( PlotFileName, nLevels, MF_plt, VarNames, &
             amrex_geom, Time / UnitsDisplay % TimeUnit, StepNo, &
             amrex_ref_ratio )

    DO iLevel = 0, nLevels-1

      CALL amrex_multifab_destroy ( MF_plt(iLevel) )

    END DO

    DEALLOCATE( VarNames )

  END SUBROUTINE WriteFieldsAMReX_PlotFile


  SUBROUTINE ReadCheckpointFile

    USE MF_FieldsModule, ONLY: &
      MF_uGF, &
      MF_uCF, &
      MF_uPF, &
      MF_uAF, &
      MF_uDF, &
      FluxRegister

    IMPLICIT NONE

    INTEGER     :: iLevel, FinestLevel
    TYPE(c_ptr) :: pBA(0:nLevels-1)
    TYPE(c_ptr) :: pDM(0:nLevels-1)
    TYPE(c_ptr) :: pGF(0:nLevels-1)
    TYPE(c_ptr) :: pCF(0:nLevels-1)
    TYPE(c_ptr) :: amrcore

    TYPE(amrex_box)       :: BX
    TYPE(amrex_distromap) :: DM  (0:nLevels-1)
    TYPE(amrex_boxarray)  :: BA  (0:nLevels-1)
    TYPE(amrex_geometry)  :: GEOM(0:nLevels-1)

    amrcore = amrex_get_amrcore()

    BX = amrex_box( [ 0, 0, 0 ], [ nX(1)-1, nX(2)-1, nX(3)-1 ] )

    DO iLevel = 0, nLevels-1

      CALL amrex_boxarray_build ( BA(iLevel), BX )

      CALL BA(iLevel) % maxSize( MaxGridSizeX )

      CALL amrex_geometry_build( GEOM(iLevel), BX )

      CALL amrex_distromap_build( DM(iLevel), BA(iLevel) )

    END DO

    pBA(0:nLevels-1) = BA(0:nLevels-1) % P
    pDM(0:nLevels-1) = DM(0:nLevels-1) % P

    FinestLevel = nLevels-1

    CALL ReadHeaderAndBoxArrayData &
           ( FinestLevel, StepNo, dt, t_new, pBA, pDM, iRestart )

    DO iLevel = 0, nLevels-1

      BA(iLevel) = pBA(iLevel)
      DM(iLevel) = pDM(iLevel)

    END DO

    DO iLevel = 0, nLevels-1

      CALL amrex_fi_set_boxarray ( iLevel, BA(iLevel) % P, amrcore )
      CALL amrex_fi_set_distromap( iLevel, DM(iLevel) % P, amrcore )

    END DO

    DO iLevel = 0, nLevels-1

      CALL amrex_multifab_build &
             ( MF_uGF(iLevel), BA(iLevel), DM(iLevel), nDOFX * nGF, swX )

      CALL amrex_multifab_build &
             ( MF_uCF(iLevel), BA(iLevel), DM(iLevel), nDOFX * nCF, swX )

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
      IF( iLevel .GT. 0 .AND. do_reflux ) &
        CALL amrex_fluxregister_build &
               ( FluxRegister(iLevel), BA(iLevel), DM(iLevel), &
                 amrex_ref_ratio(iLevel-1), iLevel, nDOFX_X1*nCF )

    END DO

    pGF(0:nLevels-1) = MF_uGF(0:nLevels-1) % P
    pCF(0:nLevels-1) = MF_uCF(0:nLevels-1) % P

    CALL ReadMultiFabData( nLevels-1, pGF, 0, iRestart )
    CALL ReadMultiFabData( nLevels-1, pCF, 1, iRestart )

    DO iLevel = 0, nLevels-1

      CALL FillPatch( iLevel, t_new(0), MF_uGF )
      CALL FillPatch( iLevel, t_new(0), MF_uCF )

    END DO

    CALL amrex_fi_set_finest_level( nLevels-1, amrcore )

  END SUBROUTINE ReadCheckpointFile


  ! --- PRIVATE SUBROUTINES ---


  SUBROUTINE ComputeCellAverage_MF &
    ( nFd, MF_uGF, MF, iOS, Field, MF_plt )

    INTEGER,              INTENT(in)    :: nFd, iOS
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF
    TYPE(amrex_multifab), INTENT(in)    :: MF
    CHARACTER(2),         INTENT(in)    :: Field
    TYPE(amrex_multifab), INTENT(inout) :: MF_plt

    INTEGER                       :: iX1, iX2, iX3, iFd
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

      G     => MF_uGF % DataPtr( MFI )
      U     => MF     % DataPtr( MFI )
      U_plt => MF_plt % DataPtr( MFI )

      BX = MFI % TileBox()

      lo_G = LBOUND( G ); hi_G = UBOUND( G )
      lo_U = LBOUND( U ); hi_U = UBOUND( U )

      DO iX3 = BX % lo(3), BX % hi(3)
      DO iX2 = BX % lo(2), BX % hi(2)
      DO iX1 = BX % lo(1), BX % hi(1)

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

  END SUBROUTINE ComputeCellAverage_MF


  SUBROUTINE ConvertUnits( Field, nFd, iOS, U_plt )

    CHARACTER(2), INTENT(in)    :: Field
    INTEGER,      INTENT(in)    :: nFd, iOS
    REAL(DP),     INTENT(inout) :: U_plt(:,:,:,:)

    INTEGER :: iFd

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

      CASE DEFAULT

        RETURN

    END SELECT

  END SUBROUTINE ConvertUnits


  SUBROUTINE WriteMPI( MF_uGF, MF_plt )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF
    TYPE(amrex_multifab), INTENT(inout) :: MF_plt

    INTEGER            :: iX1, iX2, iX3
    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI

    REAL(DP), CONTIGUOUS, POINTER :: U_plt(:,:,:,:)

    CALL amrex_mfiter_build( MFI, MF_uGF, tiling = UseTiling )

    DO WHILE( MFI % next() )

      U_plt => MF_plt % DataPtr( MFI )

      BX = MFI % TileBox()

      DO iX3 = BX % lo(3), BX % hi(3)
      DO iX2 = BX % lo(2), BX % hi(2)
      DO iX1 = BX % lo(1), BX % hi(1)

        U_plt(iX1,iX2,iX3,1) = amrex_parallel_myproc()

      END DO
      END DO
      END DO

    END DO

    CALL amrex_mfiter_destroy( MFI )

  END SUBROUTINE WriteMPI


  SUBROUTINE WriteMesh( iLevel, MF_uGF, MF_plt )

    INTEGER,              INTENT(in)    :: iLevel
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF
    TYPE(amrex_multifab), INTENT(inout) :: MF_plt

    INTEGER            :: iX1, iX2, iX3
    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI

    REAL(DP), CONTIGUOUS, POINTER :: U_plt(:,:,:,:)

    CALL CreateMesh_MF( iLevel, MeshX )

    CALL amrex_mfiter_build( MFI, MF_uGF, tiling = UseTiling )

    DO WHILE( MFI % next() )

      U_plt => MF_plt % DataPtr( MFI )

      BX = MFI % TileBox()

      DO iX3 = BX % lo(3), BX % hi(3)
      DO iX2 = BX % lo(2), BX % hi(2)
      DO iX1 = BX % lo(1), BX % hi(1)

        U_plt(iX1,iX2,iX3,2) = MeshX(1) % Center(iX1)
        U_plt(iX1,iX2,iX3,3) = MeshX(2) % Center(iX2)
        U_plt(iX1,iX2,iX3,4) = MeshX(3) % Center(iX3)

        U_plt(iX1,iX2,iX3,5) = MeshX(1) % Width(iX1)
        U_plt(iX1,iX2,iX3,6) = MeshX(2) % Width(iX2)
        U_plt(iX1,iX2,iX3,7) = MeshX(3) % Width(iX3)

      END DO
      END DO
      END DO

    END DO

    CALL amrex_mfiter_destroy( MFI )

    CALL DestroyMesh_MF( MeshX )

  END SUBROUTINE WriteMesh


END MODULE InputOutputModuleAMReX
