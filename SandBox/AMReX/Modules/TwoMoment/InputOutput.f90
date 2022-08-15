MODULE InputOutput

  ! --- AMReX Modules ---

  USE amrex_fort_module, ONLY: &
    AR => amrex_real

  USE ISO_C_BINDING
  USE amrex_base_module
  USE amrex_amr_module

  ! --- thornado Modules ---

  USE ProgramHeaderModule,     ONLY: &
    nDOFZ,                  &
    nDOFE,                  &
    iZ_B0,                  &
    iZ_E0
  USE ReferenceElementModule, ONLY: &
    Weights_q
  USE ReferenceElementModuleE, ONLY: &
    WeightsE
  USE RadiationFieldsModule,       ONLY: &
    nCR,    &
    uCR,    &
    nPR,    &
    uPR,    &
    nSpecies, &
    iPR_D,  &
    iPR_I1, &
    iPR_I2, &
    iPR_I3, &
    unitsPR, &
    iCR_N,  &
    iCR_G1, &
    iCR_G2, &
    iCR_G3, &
    unitsCR, &
    ShortNamesCR, &
    ShortNamesPR
  USE UnitsModule,             ONLY: &
    Joule,  &
    Kelvin, &
    UnitsDisplay


  USE MeshModule,              ONLY: &
    MeshE,          &
    NodeCoordinate

  USE InputParsingModule,                      ONLY: &
    t_end,                     &
    t,                         &
    dt,                        &
    t_wrt,                     &
    dt_wrt,                    &
    t_chk,                     &
    dt_chk,                    &
    CFL,                       &
    nNodes,                    &
    nStages,                   &
    nX,                        &
    nE,                        &
    swX,                       &
    swE,                       &
    bcX,                       &
    xL,                        &
    xR,                        &
    ProgramName,               &
    CoordSys,                  &
    StepNo,                    &
    nLevels,                   &
    iRestart,                  &
    MaxGridSizeX,              &
    BA,                        &
    DM,                        &
    GEOM,                      &
    StepNo,                    &
    UseTiling,                 &
    MyAmrInit
  USE MF_FieldsModule,    ONLY: &
    MF_uCR, &
    MF_uPR

  IMPLICIT NONE
  PRIVATE

  CHARACTER(8) :: BaseFileName = 'thornado'


  REAL(AR), PARAMETER :: One = 1.0_AR

  PUBLIC :: ReadCheckpointFile
  PUBLIC :: WriteFieldsAMReX_Checkpoint
  PUBLIC :: WriteFieldsAMReX_PlotFile

  INTERFACE

    SUBROUTINE WriteFieldsAMReX_Checkpoint &
                 ( StepNo, nLevels, dt, time, t_wrt, pBA, &
                   pMF_uCR, pMF_uPR ) BIND(c)
       IMPORT
       IMPLICIT NONE
       INTEGER(c_int), INTENT(in) :: StepNo(*)
       INTEGER(c_int), VALUE      :: nLevels
       REAL(AR),       INTENT(in) :: dt(*), time(*), t_wrt
       TYPE(c_ptr),    INTENT(in) :: pBA(*)
       TYPE(c_ptr),    INTENT(in) :: pMF_uCR(*)
       TYPE(c_ptr),    INTENT(in) :: pMF_uPR(*)

    END SUBROUTINE WriteFieldsAMReX_Checkpoint

    SUBROUTINE ReadHeaderAndBoxArrayData &
                 ( FinestLevel, StepNo, dt, time, t_wrt, &
                   pBA, pDM, iChkFile ) BIND(c)
      IMPORT
      IMPLICIT NONE
      INTEGER(c_int), INTENT(out) :: FinestLevel
      INTEGER(c_int), INTENT(out) :: StepNo(*)
      REAL(AR),       INTENT(out) :: dt(*), time(*), t_wrt
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

 SUBROUTINE ReadCheckpointFile( iChkFile )

    USE InputParsingModule

    IMPLICIT NONE

    INTEGER, INTENT(in) :: iChkFile

    INTEGER         :: iLevel, FinestLevel
    TYPE(c_ptr)     :: pBA(0:amrex_max_level)
    TYPE(c_ptr)     :: pDM(0:amrex_max_level)
    TYPE(c_ptr)     :: pCR(0:amrex_max_level)
    TYPE(c_ptr)     :: pPR(0:amrex_max_level)
    TYPE(c_ptr)     :: amrcore
    TYPE(amrex_box) :: BX

    amrcore = amrex_get_amrcore()

    BX = amrex_box( [ 1, 1, 1 ], [ nX(1), nX(2), nX(3) ] )

    ALLOCATE( BA(0:nLevels-1) )
    DO iLevel = 0, nLevels-1
      CALL amrex_boxarray_build ( BA(iLevel), BX )
    END DO

    DO iLevel = 0, nLevels-1
      CALL BA(iLevel) % maxSize( MaxGridSizeX )
    END DO

    ALLOCATE( DM  (0:nLevels-1) )
    ALLOCATE( GEOM(0:nLevels-1) )
    DO iLevel = 0, nLevels-1
      CALL amrex_geometry_build( GEOM(iLevel), BX )
      CALL amrex_distromap_build( DM (iLevel), BA(iLevel) )
    END DO

    pBA(0:amrex_max_level) = BA(0:amrex_max_level) % P
    pDM(0:amrex_max_level) = DM(0:amrex_max_level) % P

    FinestLevel = nLevels-1

    CALL ReadHeaderAndBoxArrayData &
           ( FinestLevel, StepNo, dt, t, t_wrt, pBA, pDM, iChkFile )
    print*, 'here'
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
             ( MF_uCR(iLevel), BA(iLevel), DM(iLevel), &
               nDOFZ * nCR * ( iZ_E0( 1 ) - iZ_B0( 1 ) + 1 ) * nSpecies, swX(1) )
      CALL amrex_multifab_build &
             ( MF_uPR(iLevel), BA(iLevel), DM(iLevel), &
               nDOFZ * nPR * ( iZ_E0( 1 ) - iZ_B0( 1 ) + 1 ) * nSpecies, swX(1) )
    END DO

    pCR(0:amrex_max_level) = MF_uCR(0:amrex_max_level) % P
    pPR(0:amrex_max_level) = MF_uPR(0:amrex_max_level) % P

    CALL ReadMultiFabData( nLevels-1, pCR, 1, iChkFile )
    CALL ReadMultiFabData( nLevels-1, pPR, 2, iChkFile )

    DO iLevel = 0, nLevels-1
      CALL MF_uCR(iLevel) % Fill_Boundary( GEOM(iLevel) )
      CALL MF_uPR(iLevel) % Fill_Boundary( GEOM(iLevel) )
    END DO

    CALL amrex_fi_set_finest_level( nLevels-1, amrcore )

  END SUBROUTINE ReadCheckpointFile

  SUBROUTINE WriteFieldsAMReX_PlotFile &
    ( Time, StepNo, MF_uCR_Option, MF_uPR_Option, num_Option )

    REAL(AR),             INTENT(in)           :: Time
    INTEGER,              INTENT(in)           :: StepNo(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in), OPTIONAL :: MF_uCR_Option(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in), OPTIONAL :: MF_uPR_Option(0:nLevels-1)
    INTEGER,              INTENT(in), OPTIONAL :: num_Option


    CHARACTER(08)                   :: NumberString
    CHARACTER(32)                   :: PlotFileName
    CHARACTER(32)                   :: Names
    CHARACTER(3)                    :: iSC, iZ1C
    CHARACTER(64)                   :: nm
    LOGICAL                         :: WriteFF_C, WriteFF_P
    INTEGER                         :: iComp, iOS, iLevel, nF, iOS_CPP(3), iS, iZ1
    TYPE(amrex_box)                 :: BX
    TYPE(amrex_multifab)            :: MF_PR(0:nLevels-1)
    TYPE(amrex_geometry)            :: GEOM (0:nLevels-1)
    TYPE(amrex_boxarray)            :: BA   (0:nLevels-1)
    TYPE(amrex_distromap)           :: DM   (0:nLevels-1)
    TYPE(amrex_string), ALLOCATABLE :: VarNames(:)

    ! --- Offset for C++ indexing ---
    iOS_CPP = 1
    IF     ( amrex_spacedim .EQ. 1 )THEN
      iOS_CPP(2) = 0
      iOS_CPP(3) = 0
    ELSE IF( amrex_spacedim .EQ. 2 )THEN
      iOS_CPP(3) = 0
    END IF

    BX = amrex_box( [ 0, 0, 0 ], [ nX(1)-1, nX(2)-1, nX(3)-1 ] )

    DO iLevel = 0, nLevels-1
      CALL amrex_boxarray_build( BA(iLevel), BX )
    END DO

    DO iLevel = 0, nLevels-1
      CALL BA(iLevel) % maxSize( MaxGridSizeX )
    END DO

    DO iLevel = 0, nLevels-1
      CALL amrex_geometry_build ( GEOM(iLevel), BX )
      CALL amrex_distromap_build( DM  (iLevel), BA(iLevel) )
    END DO

    nF = 0

    WriteFF_C = .FALSE.
    IF( PRESENT( MF_uCR_Option ) )THEN
      WriteFF_C = .TRUE.
      nF = nF + nCR
    END IF

    WriteFF_P = .FALSE.
    IF( PRESENT( MF_uPR_Option ) )THEN
      WriteFF_P = .TRUE.
      nF = nF + nPR
    END IF

    IF( amrex_parallel_ioprocessor() )THEN

      WRITE(*,*)
      WRITE(*,'(A6,A,I8.8)') '', 'Writing PlotFile: ', StepNo

    END IF

    IF( nLevels-1 .EQ. 0 )THEN
      WRITE(NumberString,'(I8.8)') StepNo(:)
    ELSE
      WRITE(NumberString,'(I8.8)') StepNo(0)
    END IF

    PlotFileName = TRIM( BaseFileName ) // '_' // NumberString

    IF( PRESENT( num_Option ) )THEN
      
      WRITE(nm,*) num_Option
      nm = ADJUSTL(nm)
      PlotFileName = TRIM( BaseFileName ) // '_' // TRIM(nm)

    END IF 

    nF = nF * nE * nSpecies

    ALLOCATE( VarNames(nF) )

    iOS = 0

    IF( WriteFF_C )THEN
      DO iS = 1, nSpecies
      DO iZ1 = iZ_B0(1), iZ_E0(1)
      DO iComp = 1, nCR
        write(iZ1C,'(i3.3)') iZ1
        write(iSC,'(i3.3)') iS
        Names = ShortNamesCR(iComp) // '_' // iZ1C // '_' // iSC
        CALL amrex_string_build &
               ( VarNames(iComp + ( iZ1 - 1 ) * nCR + ( iS - 1 ) * nE * nCR + iOS), TRIM( Names ) )

      END DO
      END DO
      END DO

        iOS = iOS + nCR * nE * nSpecies
      !iOS = iOS + nCR * nE * nSpecies

    END IF

    IF( WriteFF_P )THEN

      DO iS = 1, nSpecies
      DO iZ1 = iZ_B0(1), iZ_E0(1)
      DO iComp = 1, nPR

        write(iZ1C,'(i3.3)') iZ1
        write(iSC,'(i3.3)') iS

        Names = ShortNamesPR(iComp) // '_' // iZ1C // '_' // iSC
        CALL amrex_string_build &
               ( VarNames(iComp  + ( iZ1 - 1 ) * nPR + ( iS - 1 ) * nPR * nE + iOS ), TRIM( Names ) )

      END DO
      END DO
      END DO

      !iOS = iOS + nPR * nE * nSpecies

        iOS = iOS + nPR * nE * nSpecies
    END IF

    DO iLevel = 0, nLevels-1

      CALL amrex_multifab_build &
             ( MF_PR(iLevel), BA(iLevel), DM(iLevel), nF, 0 )
      CALL MF_PR(iLevel) % setVal( 0.0d0 )

      iOS = 0


      IF( WriteFF_C )THEN
        CALL MF_ComputeCellAverage &
          ( nCR, MF_uCR_Option(iLevel), MF_PR(iLevel), iOS, iOS_CPP, 'CR' )
        iOS = iOS + nCR * nE * nSpecies
      END IF

      IF( WriteFF_P )THEN
        CALL MF_ComputeCellAverage &
          ( nPR, MF_uPR_Option(iLevel), MF_PR(iLevel), iOS, iOS_CPP, 'PR' )
        iOS = iOS + nPR * nE * nSpecies
      END IF



    END DO ! End of loop over levels
    CALL amrex_write_plotfile &
           ( PlotFileName, nLevels, MF_PR, VarNames, &
             GEOM, Time / UnitsDisplay % TimeUnit, StepNo, amrex_ref_ratio )

    DO iLevel = 0, nLevels-1
      CALL amrex_multifab_destroy ( MF_PR(iLevel) )
      CALL amrex_distromap_destroy( DM   (iLevel) )
      CALL amrex_geometry_destroy ( GEOM (iLevel) )
      CALL amrex_boxarray_destroy ( BA   (iLevel) )
    END DO

    DEALLOCATE( VarNames )

  END SUBROUTINE WriteFieldsAMReX_PlotFile


  SUBROUTINE MF_ComputeCellAverage( nComp, MF, MF_A, iOS, iOS_CPP, Field )

    INTEGER,              INTENT(in)    :: nComp, iOS, iOS_CPP(3)
    TYPE(amrex_multifab), INTENT(in)    :: MF
    TYPE(amrex_multifab), INTENT(inout) :: MF_A
    CHARACTER(2),         INTENT(in)    :: Field

    INTEGER                       :: iX1, iX2, iX3, iComp, iS, iZ1, iNodeZ, iNodeE
    INTEGER                       :: lo(4), hi(4)
    REAL(AR)                      :: u_K( nDOFZ, nE, nComp, nSpecies )
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(AR), CONTIGUOUS, POINTER :: u  (:,:,:,:)
    REAL(AR), CONTIGUOUS, POINTER :: u_A(:,:,:,:)
    REAL(AR)                      :: Eq(1:nDOFE), E(1:nDOFZ), V_K


    CALL amrex_mfiter_build( MFI, MF, tiling = UseTiling )

    DO WHILE( MFI % next() )

      u   => MF   % DataPtr( MFI )
      u_A => MF_A % DataPtr( MFI )

      BX = MFI % tilebox()

      lo = LBOUND( u ); hi = UBOUND( u )

      DO iX3 = BX % lo(3), BX % hi(3)
      DO iX2 = BX % lo(2), BX % hi(2)
      DO iX1 = BX % lo(1), BX % hi(1)

        u_K(1:nDOFZ,1:nE,1:nComp,1:nSpecies) &
          = RESHAPE( u(iX1,iX2,iX3,lo(4):hi(4)), [ nDOFZ, nE, nComp, nSpecies ] )
      DO iS = 1, nSpecies
      DO iZ1 = 1, nE
        ! --- Compute cell-average ---

        DO iNodeZ = 1, nDOFZ

          iNodeE = MOD( (iNodeZ-1)        , nDOFE ) + 1

          E(iNodeZ) = NodeCoordinate( MeshE, iZ1, iNodeE )
          
          Eq(iNodeE)= NodeCoordinate( MeshE, iZ1, iNodeE )
 
        END DO
        
        V_K = SUM( WeightsE * Eq(:)**2 )
 
        DO iComp = 1, nComp

          u_A(iX1-iOS_CPP(1),iX2-iOS_CPP(2),iX3-iOS_CPP(3),iComp+( iZ1 - 1 ) * nComp + ( iS - 1 ) * nComp * nE + iOS) &
            = SUM( Weights_q(:) * u_K(:,iZ1,iComp,iS) * E(:)**2 ) / V_K

        END DO

      END DO
      END DO

      END DO
      END DO
      END DO

      CALL ConvertUnits( nComp, iOS, u_A, Field )

    END DO

    CALL amrex_mfiter_destroy( MFI )

  END SUBROUTINE MF_ComputeCellAverage


  SUBROUTINE ConvertUnits( nComp, iOS, U_plt, Field )

    INTEGER,      INTENT(in)    :: nComp, iOS
    REAL(AR),     INTENT(inout) :: U_plt(:,:,:,:)
    CHARACTER(2), INTENT(in)    :: Field

    INTEGER  :: iComp


    IF( Field .EQ. 'CR' )THEN

      DO iComp = 1, nComp

        U_plt(:,:,:,iComp+iOS) = U_plt(:,:,:,iComp+iOS) / unitsCR(iComp)

      END DO

    ELSE IF( Field .EQ. 'PR' )THEN

      DO iComp = 1, nComp

        U_plt(:,:,:,iComp+iOS) = U_plt(:,:,:,iComp+iOS) / unitsPR(iComp)

      END DO

    ELSE

      RETURN

    END IF

  END SUBROUTINE ConvertUnits


END MODULE InputOutput
