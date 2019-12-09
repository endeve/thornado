MODULE InputOutput

  ! --- AMReX Modules ---

  USE amrex_fort_module, ONLY: &
    AR => amrex_real

  USE ISO_C_BINDING
  USE amrex_base_module
  USE amrex_amr_module

  ! --- thornado Modules ---

  USE ProgramHeaderModule,     ONLY: &
    nDOFZ
  USE ReferenceElementModuleZ, ONLY: &
    WeightsZ_q
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
    iCR_N,  &
    iCR_G1, &
    iCR_G2, &
    iCR_G3, &
    ShortNamesCR, &
    ShortNamesPR
  USE UnitsModule,             ONLY: &
    Joule,  &
    Kelvin, &
    UnitsDisplay




  USE MyAmrModule,                      ONLY: &
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
    MyAmrInit


  IMPLICIT NONE
  PRIVATE

  CHARACTER(8) :: BaseFileName = 'thornado'


  PUBLIC :: WriteFieldsAMReX_Checkpoint
  REAL(AR), PARAMETER :: One = 1.0_AR

  PUBLIC :: WriteFieldsAMReX_PlotFile

  INTERFACE

    SUBROUTINE WriteFieldsAMReX_Checkpoint &
                 ( StepNo, nLevels, dt, time, t_wrt, pBA, &
                   pMF_uGF, pMF_uCF, pMF_uPF, pMF_uAF ) BIND(c)
       IMPORT
       IMPLICIT NONE
       INTEGER(c_int), INTENT(in) :: StepNo(*)
       INTEGER(c_int), VALUE      :: nLevels
       REAL(AR),       INTENT(in) :: dt(*), time(*), t_wrt
       TYPE(c_ptr),    INTENT(in) :: pBA(*)
       TYPE(c_ptr),    INTENT(in) :: pMF_uGF(*)
       TYPE(c_ptr),    INTENT(in) :: pMF_uCF(*)
       TYPE(c_ptr),    INTENT(in) :: pMF_uPF(*)
       TYPE(c_ptr),    INTENT(in) :: pMF_uAF(*)
    END SUBROUTINE WriteFieldsAMReX_Checkpoint

  END INTERFACE



  CONTAINS

  SUBROUTINE WriteFieldsAMReX_PlotFile &
    ( Time, StepNo, MF_uCR_Option, MF_uPR_Option )

    REAL(AR),             INTENT(in)           :: Time
    INTEGER,              INTENT(in)           :: StepNo(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in), OPTIONAL :: MF_uCR_Option(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in), OPTIONAL :: MF_uPR_Option(0:nLevels-1)

    CHARACTER(08)                   :: NumberString
    CHARACTER(32)                   :: PlotFileName
    LOGICAL                         :: WriteFF_C, WriteFF_P
    INTEGER                         :: iComp, iOS, iLevel, nF, iOS_CPP(3)
    TYPE(amrex_mfiter)              :: MFI
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

    ALLOCATE( VarNames(nF) )

    iOS = 0

    IF( WriteFF_C )THEN
      DO iComp = 1, nCR
        CALL amrex_string_build &
               ( VarNames(iComp + iOS), TRIM( ShortNamesCR(iComp) ) )
      END DO
      iOS = iOS + nCR
    END IF

    IF( WriteFF_P )THEN
      DO iComp = 1, nPR
        CALL amrex_string_build &
               ( VarNames(iComp + iOS), TRIM( ShortNamesPR(iComp) ) )
      END DO
      iOS = iOS + nPR
    END IF

    DO iLevel = 0, nLevels-1

      CALL amrex_multifab_build &
             ( MF_PR(iLevel), BA(iLevel), DM(iLevel), nF, 0 )
      CALL MF_PR(iLevel) % setVal( 0.0d0 )

      CALL amrex_mfiter_build( MFI, MF_PR(iLevel), tiling = .TRUE. )

      iOS = 0


      IF( WriteFF_C )THEN
        CALL MF_ComputeCellAverage &
          ( nCR, MF_uCR_Option(iLevel), MF_PR(iLevel), iOS, iOS_CPP, 'CR' )
        iOS = iOS + nCR
      END IF

      IF( WriteFF_P )THEN
        CALL MF_ComputeCellAverage &
          ( nPR, MF_uPR_Option(iLevel), MF_PR(iLevel), iOS, iOS_CPP, 'PR' )
        iOS = iOS + nPR
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

    INTEGER                       :: iX1, iX2, iX3, iComp, iS, iZ1
    INTEGER                       :: lo(4), hi(4)
    REAL(AR)                      :: u_K( nDOFZ, nE, nComp, nSpecies )
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(AR), CONTIGUOUS, POINTER :: u  (:,:,:,:)
    REAL(AR), CONTIGUOUS, POINTER :: u_A(:,:,:,:)

    CALL amrex_mfiter_build( MFI, MF, tiling = .TRUE. )

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
        DO iComp = 1, nComp

          u_A(iX1-iOS_CPP(1),iX2-iOS_CPP(2),iX3-iOS_CPP(3),iComp+iOS) &
            = DOT_PRODUCT( WeightsZ_q(:), u_K(:,iZ1,iComp,iS) )

        END DO

      END DO 
      END DO

      END DO
      END DO
      END DO


    END DO

    CALL amrex_mfiter_destroy( MFI )

  END SUBROUTINE MF_ComputeCellAverage




END MODULE InputOutput
