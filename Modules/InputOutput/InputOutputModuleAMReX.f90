MODULE InputOutputModuleAMReX

  ! --- AMReX Modules ---

  USE ISO_C_BINDING
  USE amrex_base_module
  USE amrex_amr_module

  ! --- For getting MPI process info (could also
  !     add this to amrex_base_module) ---
  USE amrex_paralleldescriptor_module

  ! --- thornado Modules ---

  USE KindModule,              ONLY: &
    DP
  USE ProgramHeaderModule,     ONLY: &
    nDOFX
  USE ReferenceElementModuleX, ONLY: &
    WeightsX_q
  USE GeometryFieldsModule,    ONLY: &
    nGF, ShortNamesGF
  USE FluidFieldsModule,       ONLY: &
    nCF, ShortNamesCF, &
    nPF, ShortNamesPF, &
    nAF, ShortNamesAF

  IMPLICIT NONE
  PRIVATE

  CHARACTER(8) :: BaseFileName = 'thornado'

  PUBLIC :: WriteFieldsAMReX_PlotFile
  PUBLIC :: WriteFieldsAMReX_Checkpoint
  PUBLIC :: ReadCheckpointFile

  INTERFACE

    SUBROUTINE WriteFieldsAMReX_Checkpoint &
                 ( StepNo, FinestLevel, dt, time, pBA, &
                   pMF_uGF, pMF_uCF, pMF_uPF, pMF_uAF ) BIND(c)
       IMPORT
       IMPLICIT NONE
       INTEGER(c_int),   INTENT(in) :: StepNo(*)
       INTEGER(c_int),   VALUE      :: FinestLevel
       REAL(amrex_real), INTENT(in) :: dt(*), time(*)
       TYPE(c_ptr),      INTENT(in) :: pBA(*)
       TYPE(c_ptr),      INTENT(in) :: pMF_uGF(*)
       TYPE(c_ptr),      INTENT(in) :: pMF_uCF(*)
       TYPE(c_ptr),      INTENT(in) :: pMF_uPF(*)
       TYPE(c_ptr),      INTENT(in) :: pMF_uAF(*)
    END SUBROUTINE WriteFieldsAMReX_Checkpoint

    SUBROUTINE ReadHeaderAndBoxArrayData &
                 ( FinestLevel, StepNo, dt, time, pBA, pDM, iChkFile ) BIND(c)
      IMPORT
      IMPLICIT NONE
      INTEGER(c_int),   INTENT(out) :: FinestLevel(*)
      INTEGER(c_int),   INTENT(out) :: StepNo(*)
      REAL(amrex_real), INTENT(out) :: dt(*),time(*)
      TYPE(c_ptr),      INTENT(out) :: pBA(*), pDM(*)
      INTEGER(c_int),   VALUE       :: iChkFile
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

  SUBROUTINE ReadCheckpointFile( iChkFile )

    USE amrex_amr_module
    USE MyAmrModule
    USE MyAmrDataModule

    IMPLICIT NONE

    INTEGER, INTENT(in) :: iChkFile

    INTEGER               :: iLevel, FinestLevel(1), nComps
    TYPE(c_ptr)           :: pBA(0:amrex_max_level)
    TYPE(c_ptr)           :: pDM(0:amrex_max_level)
    TYPE(c_ptr)           :: pGF(0:amrex_max_level)
    TYPE(c_ptr)           :: pCF(0:amrex_max_level)
    TYPE(c_ptr)           :: pPF(0:amrex_max_level)
    TYPE(c_ptr)           :: pAF(0:amrex_max_level)
    TYPE(amrex_boxarray)  :: BA(0:amrex_max_level)
    TYPE(amrex_distromap) :: DM(0:amrex_max_level)
    TYPE(c_ptr)           :: amrcore
    TYPE(amrex_box)       :: DOMAIN, DOM
    TYPE(amrex_geometry)  :: GEOM(0:amrex_max_level)

    ! --- Parse parameter file ---
    CALL MyAmrInit

    amrcore = amrex_get_amrcore()

    ! Dummy variables 
    DOMAIN = amrex_box( [0,0,0], [1,1,1] )

    DO iLevel = 0, nLevels
      CALL amrex_geometry_build( GEOM(iLevel), DOMAIN )
    END DO

    DOM = GEOM(0) % DOMAIN

    DO iLevel = 0, nLevels
      CALL amrex_boxarray_build ( BA(iLevel), DOM )
      CALL amrex_distromap_build( DM(iLevel), BA(iLevel) )
    END DO

    pBA(0:amrex_max_level) = BA(0:amrex_max_level) % P
    pDM(0:amrex_max_level) = DM(0:amrex_max_level) % P

    FinestLevel = nLevels

    CALL ReadHeaderAndBoxArrayData &
           ( FinestLevel, StepNo, dt, t, &
             pBA(0:amrex_max_level), pDM(0:amrex_max_level), iChkFile )

    DO iLevel = 0, nLevels
      BA(iLevel) = pBA(iLevel)
      DM(iLevel) = pDM(iLevel)
    END DO

    DO iLevel = 0, nLevels
      CALL amrex_fi_set_boxarray ( iLevel, BA(iLevel) % P, amrcore )
      CALL amrex_fi_set_distromap( iLevel, DM(iLevel) % P, amrcore )
    END DO

    nComps = 1

    DO iLevel = 0, nLevels
      CALL amrex_multifab_build &
             ( MF_uGF(iLevel), BA(iLevel), DM(iLevel), nDOFX * nGF, swX(1) )
      CALL amrex_multifab_build &
             ( MF_uCF(iLevel), BA(iLevel), DM(iLevel), nDOFX * nCF, swX(1) )
      CALL amrex_multifab_build &
             ( MF_uPF(iLevel), BA(iLevel), DM(iLevel), nDOFX * nPF, swX(1) )
      CALL amrex_multifab_build &
             ( MF_uAF(iLevel), BA(iLevel), DM(iLevel), nDOFX * nAF, swX(1) )
    END DO

    pGF(0:amrex_max_level) = MF_uGF(0:amrex_max_level) % P
    pCF(0:amrex_max_level) = MF_uCF(0:amrex_max_level) % P
    pPF(0:amrex_max_level) = MF_uPF(0:amrex_max_level) % P
    pAF(0:amrex_max_level) = MF_uAF(0:amrex_max_level) % P

    CALL ReadMultiFabData( nLevels, pGF(0:amrex_max_level), 0, iChkFile )
    CALL ReadMultiFabData( nLevels, pCF(0:amrex_max_level), 1, iChkFile )
    CALL ReadMultiFabData( nLevels, pPF(0:amrex_max_level), 2, iChkFile )
    CALL ReadMultiFabData( nLevels, pAF(0:amrex_max_level), 3, iChkFile )

    DO iLevel = 0, nLevels
      CALL MF_uGF(iLevel) % Fill_Boundary( GEOM(iLevel) )
      CALL MF_uCF(iLevel) % Fill_Boundary( GEOM(iLevel) )
      CALL MF_uPF(iLevel) % Fill_Boundary( GEOM(iLevel) )
      CALL MF_uAF(iLevel) % Fill_Boundary( GEOM(iLevel) )
    END DO

    CALL amrex_fi_set_finest_level( FinestLevel(1), amrcore )
	
  END SUBROUTINE ReadCheckpointFile


  SUBROUTINE WriteFieldsAMReX_PlotFile &
    ( Time, nLevels, GEOM, StepNo, &
      MF_uGF_Option, MF_uCF_Option, MF_uPF_Option, MF_uAF_Option )

    REAL(DP),             INTENT(in)           :: Time
    INTEGER,              INTENT(in)           :: nLevels
    TYPE(amrex_geometry), INTENT(in)           :: GEOM(0:nLevels)
    INTEGER,              INTENT(in)           :: StepNo(0:nLevels)
    TYPE(amrex_multifab), INTENT(in), OPTIONAL :: MF_uGF_Option(0:nLevels)
    TYPE(amrex_multifab), INTENT(in), OPTIONAL :: MF_uCF_Option(0:nLevels)
    TYPE(amrex_multifab), INTENT(in), OPTIONAL :: MF_uPF_Option(0:nLevels)
    TYPE(amrex_multifab), INTENT(in), OPTIONAL :: MF_uAF_Option(0:nLevels)

    CHARACTER(08)                   :: NumberString
    CHARACTER(32)                   :: PlotFileName
    LOGICAL                         :: WriteGF
    LOGICAL                         :: WriteFF_C, WriteFF_P, WriteFF_A
    INTEGER                         :: iComp, iOS, iLevel, nF
    INTEGER                         :: MyProc
    TYPE(amrex_multifab)            :: MF_PF(0:nLevels)
    TYPE(amrex_boxarray)            :: BA(0:nLevels)
    TYPE(amrex_distromap)           :: DM(0:nLevels)
    TYPE(amrex_string), ALLOCATABLE :: VarNames(:)

    ! --- Only needed to get and write BX % lo and BX % hi ---
    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    nF = 0
    WriteGF   = .FALSE.
    IF( PRESENT( MF_uGF_Option ) )THEN
      WriteGF   = .TRUE.
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

    MyProc = amrex_pd_myproc()

    IF( amrex_parallel_ioprocessor() )THEN

      WRITE(*,*)
      WRITE(*,'(A4,A18,I9.8)') '', 'Writing PlotFile: ', StepNo

    END IF

    IF( nLevels .EQ. 0 )THEN
      WRITE(NumberString,'(I8.8)') StepNo(:)
    ELSE
      WRITE(NumberString,'(I8.8)') StepNo(0)
    END IF

    PlotFileName = TRIM( BaseFileName ) // '_' // NumberString

    ALLOCATE( VarNames(nF) )

    iOS = 0
    IF( WriteGF )THEN
      DO iComp = 1, nGF
        CALL amrex_string_build &
               ( VarNames(iComp + iOS), TRIM( ShortNamesGF(iComp) ) )
      END DO
      iOS = iOS + nGF
      DO iLevel = 0, nLevels
        BA(iLevel) = MF_uGF_Option(iLevel) % BA
        DM(iLevel) = MF_uGF_Option(iLevel) % DM
      END DO
    END IF

    IF( WriteFF_C )THEN
      DO iComp = 1, nCF
        CALL amrex_string_build &
               ( VarNames(iComp + iOS), TRIM( ShortNamesCF(iComp) ) )
      END DO
      iOS = iOS + nCF
      DO iLevel = 0, nLevels
        BA(iLevel) = MF_uCF_Option(iLevel) % BA
        DM(iLevel) = MF_uCF_Option(iLevel) % DM
      END DO
    END IF

    IF( WriteFF_P )THEN
      DO iComp = 1, nPF
        CALL amrex_string_build &
               ( VarNames(iComp + iOS), TRIM( ShortNamesPF(iComp) ) )
      END DO
      iOS = iOS + nPF
      DO iLevel = 0, nLevels
        BA(iLevel) = MF_uPF_Option(iLevel) % BA
        DM(iLevel) = MF_uPF_Option(iLevel) % DM
      END DO
    END IF

    IF( WriteFF_A )THEN
      DO iComp = 1, nAF
        CALL amrex_string_build &
               ( VarNames(iComp + iOS), TRIM( ShortNamesAF(iComp) ) )
      END DO
      iOS = iOS + nAF
      DO iLevel = 0, nLevels
        BA(iLevel) = MF_uAF_Option(iLevel) % BA
        DM(iLevel) = MF_uAF_Option(iLevel) % DM
      END DO
    END IF

    DO iLevel = 0, nLevels

      CALL amrex_multifab_build &
             ( MF_PF(iLevel), BA(iLevel), DM(iLevel), nF, 0 )
      CALL MF_PF(iLevel) % setVal( 0.0d0 )

      CALL amrex_mfiter_build( MFI, MF_PF(iLevel), tiling = .TRUE. )

      BX = MFI % tilebox()
!!$      WRITE(*,'(A,I2.2,A,3I3.2,A,3I3.2)') &
!!$        'MyProc = ', MyProc, ': lo = ', BX % lo, ', hi = ', BX % hi
!!$      CALL amrex_print( BX )

      iOS = 0
      IF( WriteGF   )THEN
        CALL MF_ComputeCellAverage &
          ( nGF, MF_uGF_Option(iLevel), MF_PF(iLevel), iOS )
        iOS = iOS + nGF
      END IF
      IF( WriteFF_C )THEN
        CALL MF_ComputeCellAverage &
          ( nCF, MF_uCF_Option(iLevel), MF_PF(iLevel), iOS )
        iOS = iOS + nCF
      END IF
      IF( WriteFF_P )THEN
        CALL MF_ComputeCellAverage &
          ( nPF, MF_uPF_Option(iLevel), MF_PF(iLevel), iOS )
        iOS = iOS + nPF
      END IF
      IF( WriteFF_A )THEN
        CALL MF_ComputeCellAverage &
          ( nAF, MF_uAF_Option(iLevel), MF_PF(iLevel), iOS )
        iOS = iOS + nAF
      END IF

    END DO ! End of loop over levels

    CALL amrex_write_plotfile &
           ( PlotFileName, nLevels+1, MF_PF, VarNames, &
             GEOM, Time, StepNo, amrex_ref_ratio )

    DO iLevel = 0, nLevels
      CALL amrex_multifab_destroy( MF_PF(iLevel) )
    END DO

    DEALLOCATE( VarNames )

  END SUBROUTINE WriteFieldsAMReX_PlotFile


  SUBROUTINE MF_ComputeCellAverage( nComp, MF, MF_A, iOS )

    INTEGER,              INTENT(in   ) :: nComp, iOS
    TYPE(amrex_multifab), INTENT(in   ) :: MF
    TYPE(amrex_multifab), INTENT(inout) :: MF_A

    INTEGER            :: iX1, iX2, iX3, iComp
    INTEGER            :: lo(4), hi(4)
    REAL(amrex_real)   :: u_K(nDOFX,nComp)
    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI
    REAL(amrex_real), CONTIGUOUS, POINTER :: u  (:,:,:,:)
    REAL(amrex_real), CONTIGUOUS, POINTER :: u_A(:,:,:,:)

    CALL amrex_mfiter_build( MFI, MF, tiling = .TRUE. )

    DO WHILE( MFI % next() )

      u   => MF   % DataPtr( MFI )
      u_A => MF_A % DataPtr( MFI )

      BX = MFI % tilebox()

      lo = LBOUND( u ); hi = UBOUND( u )

      DO iX3 = BX % lo(3), BX % hi(3)
      DO iX2 = BX % lo(2), BX % hi(2)
      DO iX1 = BX % lo(1), BX % hi(1)

        u_K(1:nDOFX,1:nComp) &
          = RESHAPE( u(iX1,iX2,iX3,lo(4):hi(4)), [ nDOFX, nComp ] )

        ! --- Compute cell-average ---
        DO iComp = 1, nComp

          u_A(iX1,iX2,iX3,iComp + iOS) &
            = DOT_PRODUCT( WeightsX_q(:), u_K(:,iComp) )

        END DO

      END DO
      END DO
      END DO

    END DO

    CALL amrex_mfiter_destroy( MFI )

  END SUBROUTINE MF_ComputeCellAverage


END MODULE InputOutputModuleAMReX
