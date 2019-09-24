MODULE InputOutputModuleAMReX

  ! --- AMReX Modules ---

  USE ISO_C_BINDING
  USE amrex_base_module
  USE amrex_amr_module

  ! --- thornado Modules ---

  USE KindModule,              ONLY: &
    DP
  USE ProgramHeaderModule,     ONLY: &
    nDOFX
  USE ReferenceElementModuleX, ONLY: &
    WeightsX_q
  USE GeometryFieldsModule
  USE FluidFieldsModule

  ! --- thornado Modules ---
  USE ProgramHeaderModule,     ONLY: &
    InitializeProgramHeader
  USE ReferenceElementModuleX, ONLY: &
    InitializeReferenceElementX, &
    FinalizeReferenceElementX
  USE FluidFieldsModule,       ONLY: &
    nCF, nPF, nAF
  USE GeometryFieldsModule,    ONLY: &
    nGF
  USE UnitsModule,             ONLY: &
    Millisecond, &
    Gram, &
    Centimeter, &
    Kilometer, &
    Second, &
    Kelvin, &
    MeV, &
    Erg, &
    BoltzmannConstant

  ! --- Local Modules ---
  USE MyAmrModule
  USE MyAmrDataModule
  USE MF_UtilitiesModule, ONLY: &
    AMReX2thornado, &
    ShowVariableFromMultiFab

  IMPLICIT NONE
  PRIVATE

  CHARACTER(8) :: BaseFileName = 'thornado'

  PUBLIC :: ReadCheckpointFile
  PUBLIC :: WriteFieldsAMReX_Checkpoint
  PUBLIC :: WriteFieldsAMReX_PlotFile

  INTERFACE

    SUBROUTINE WriteFieldsAMReX_Checkpoint &
                 ( StepNo, FinestLevel, dt, time, t_wrt, pBA, &
                   pMF_uGF, pMF_uCF, pMF_uPF, pMF_uAF ) BIND(c)
       IMPORT
       IMPLICIT NONE
       INTEGER(c_int),   INTENT(in) :: StepNo(*)
       INTEGER(c_int),   VALUE      :: FinestLevel
       REAL(amrex_real), INTENT(in) :: dt(*), time(*), t_wrt
       TYPE(c_ptr),      INTENT(in) :: pBA(*)
       TYPE(c_ptr),      INTENT(in) :: pMF_uGF(*)
       TYPE(c_ptr),      INTENT(in) :: pMF_uCF(*)
       TYPE(c_ptr),      INTENT(in) :: pMF_uPF(*)
       TYPE(c_ptr),      INTENT(in) :: pMF_uAF(*)
    END SUBROUTINE WriteFieldsAMReX_Checkpoint

    SUBROUTINE ReadHeaderAndBoxArrayData &
                 ( FinestLevel, StepNo, dt, time, t_wrt, &
                   pBA, pDM, iChkFile ) BIND(c)
      IMPORT
      IMPLICIT NONE
      INTEGER(c_int),   INTENT(out) :: FinestLevel(*)
      INTEGER(c_int),   INTENT(out) :: StepNo(*)
      REAL(amrex_real), INTENT(out) :: dt(*), time(*), t_wrt
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

    INTEGER         :: iLevel, FinestLevel(1)
    TYPE(c_ptr)     :: pBA(0:amrex_max_level)
    TYPE(c_ptr)     :: pDM(0:amrex_max_level)
    TYPE(c_ptr)     :: pGF(0:amrex_max_level)
    TYPE(c_ptr)     :: pCF(0:amrex_max_level)
    TYPE(c_ptr)     :: pPF(0:amrex_max_level)
    TYPE(c_ptr)     :: pAF(0:amrex_max_level)
    TYPE(c_ptr)     :: amrcore
    TYPE(amrex_box) :: BX

    amrcore = amrex_get_amrcore()

    BX = amrex_box( [ 1, 1, 1 ], [ nX(1), nX(2), nX(3) ] )

    ALLOCATE( BA(0:nLevels) )
    DO iLevel = 0, nLevels
      CALL amrex_boxarray_build ( BA(iLevel), BX )
    END DO

    DO iLevel = 0, nLevels
      CALL BA(iLevel) % maxSize( MaxGridSize )
    END DO

    ALLOCATE( DM  (0:nLevels) )
    ALLOCATE( GEOM(0:nLevels) )
    DO iLevel = 0, nLevels
      CALL amrex_geometry_build( GEOM(iLevel), BX )
      CALL amrex_distromap_build( DM (iLevel), BA(iLevel) )
    END DO

    pBA(0:amrex_max_level) = BA(0:amrex_max_level) % P
    pDM(0:amrex_max_level) = DM(0:amrex_max_level) % P

    FinestLevel = nLevels

    CALL ReadHeaderAndBoxArrayData &
           ( FinestLevel, StepNo, dt, t, t_wrt, &
             pBA(0:nLevels), pDM(0:nLevels), iChkFile )

    DO iLevel = 0, nLevels
      BA(iLevel) = pBA(iLevel)
      DM(iLevel) = pDM(iLevel)
    END DO

    DO iLevel = 0, nLevels
      CALL amrex_fi_set_boxarray ( iLevel, BA(iLevel) % P, amrcore )
      CALL amrex_fi_set_distromap( iLevel, DM(iLevel) % P, amrcore )
    END DO

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

    CALL amrex_fi_set_finest_level( nLevels, amrcore )

  END SUBROUTINE ReadCheckpointFile


  SUBROUTINE WriteFieldsAMReX_PlotFile &
    ( Time, StepNo, MF_uGF_Option, MF_uCF_Option, MF_uPF_Option, MF_uAF_Option )

    REAL(amrex_real),     INTENT(in)           :: Time
    INTEGER,              INTENT(in)           :: StepNo(0:nLevels)
    TYPE(amrex_multifab), INTENT(in), OPTIONAL :: MF_uGF_Option(0:nLevels)
    TYPE(amrex_multifab), INTENT(in), OPTIONAL :: MF_uCF_Option(0:nLevels)
    TYPE(amrex_multifab), INTENT(in), OPTIONAL :: MF_uPF_Option(0:nLevels)
    TYPE(amrex_multifab), INTENT(in), OPTIONAL :: MF_uAF_Option(0:nLevels)

    CHARACTER(08)                   :: NumberString
    CHARACTER(32)                   :: PlotFileName
    LOGICAL                         :: WriteGF
    LOGICAL                         :: WriteFF_C, WriteFF_P, WriteFF_A
    INTEGER                         :: iComp, iOS, iLevel, nF, iOS_CPP(3)
    TYPE(amrex_mfiter)              :: MFI
    TYPE(amrex_box)                 :: BX
    TYPE(amrex_multifab)            :: MF_PF(0:nLevels)
    TYPE(amrex_geometry)            :: GEOM (0:nLevels)
    TYPE(amrex_boxarray)            :: BA   (0:nLevels)
    TYPE(amrex_distromap)           :: DM   (0:nLevels)
    TYPE(amrex_string), ALLOCATABLE :: VarNames(:)
    REAL(amrex_real)                :: TimeUnit

    ! --- Offset for C++ indexing ---
    iOS_CPP = 1
    IF     ( amrex_spacedim .EQ. 1 )THEN
      iOS_CPP(2) = 0
      iOS_CPP(3) = 0
    ELSE IF( amrex_spacedim .EQ. 2 )THEN
      iOS_CPP(3) = 0
    END IF

    BX = amrex_box( [ 0, 0, 0 ], [ nX(1)-1, nX(2)-1, nX(3)-1 ] )

    DO iLevel = 0, nLevels
      CALL amrex_boxarray_build( BA(iLevel), BX )
    END DO

    DO iLevel = 0, nLevels
      CALL BA(iLevel) % maxSize( MaxGridSize )
    END DO

    DO iLevel = 0, nLevels
      CALL amrex_geometry_build ( GEOM(iLevel), BX )
      CALL amrex_distromap_build( DM  (iLevel), BA(iLevel) )
    END DO

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

    IF( amrex_parallel_ioprocessor() )THEN

      WRITE(*,*)
      WRITE(*,'(A6,A,I8.8)') '', 'Writing PlotFile: ', StepNo

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
    END IF

    IF( WriteFF_C )THEN
      DO iComp = 1, nCF
        CALL amrex_string_build &
               ( VarNames(iComp + iOS), TRIM( ShortNamesCF(iComp) ) )
      END DO
      iOS = iOS + nCF
    END IF

    IF( WriteFF_P )THEN
      DO iComp = 1, nPF
        CALL amrex_string_build &
               ( VarNames(iComp + iOS), TRIM( ShortNamesPF(iComp) ) )
      END DO
      iOS = iOS + nPF
    END IF

    IF( WriteFF_A )THEN
      DO iComp = 1, nAF
        CALL amrex_string_build &
               ( VarNames(iComp + iOS), TRIM( ShortNamesAF(iComp) ) )
      END DO
      iOS = iOS + nAF
    END IF

    DO iLevel = 0, nLevels

      CALL amrex_multifab_build &
             ( MF_PF(iLevel), BA(iLevel), DM(iLevel), nF, 0 )
      CALL MF_PF(iLevel) % setVal( 0.0d0 )

      CALL amrex_mfiter_build( MFI, MF_PF(iLevel), tiling = .TRUE. )

      iOS = 0
      IF( WriteGF   )THEN
        CALL MF_ComputeCellAverage &
          ( nGF, MF_uGF_Option(iLevel), MF_PF(iLevel), iOS, iOS_CPP, 'GF' )
        iOS = iOS + nGF
      END IF
      IF( WriteFF_C )THEN
        CALL MF_ComputeCellAverage &
          ( nCF, MF_uCF_Option(iLevel), MF_PF(iLevel), iOS, iOS_CPP, 'CF' )
        iOS = iOS + nCF
      END IF
      IF( WriteFF_P )THEN
        CALL MF_ComputeCellAverage &
          ( nPF, MF_uPF_Option(iLevel), MF_PF(iLevel), iOS, iOS_CPP, 'PF' )
        iOS = iOS + nPF
      END IF
      IF( WriteFF_A )THEN
        CALL MF_ComputeCellAverage &
          ( nAF, MF_uAF_Option(iLevel), MF_PF(iLevel), iOS, iOS_CPP, 'AF' )
        iOS = iOS + nAF
      END IF

    END DO ! End of loop over levels

    IF( UsePhysicalUnits )THEN
      TimeUnit = Millisecond
    ELSE
      TimeUnit = 1.0_amrex_real
    END IF

    CALL amrex_write_plotfile &
           ( PlotFileName, nLevels+1, MF_PF, VarNames, &
             GEOM, Time / TimeUnit, StepNo, amrex_ref_ratio )

    DO iLevel = 0, nLevels
      CALL amrex_multifab_destroy ( MF_PF(iLevel) )
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

          u_A(iX1-iOS_CPP(1),iX2-iOS_CPP(2),iX3-iOS_CPP(3),iComp+iOS) &
            = DOT_PRODUCT( WeightsX_q(:), u_K(:,iComp) )

        END DO

      END DO
      END DO
      END DO

      IF( UsePhysicalUnits ) &
        CALL ConvertToPhysicalUnits( nComp, iOS, u_A, Field )

    END DO

    CALL amrex_mfiter_destroy( MFI )

  END SUBROUTINE MF_ComputeCellAverage


  SUBROUTINE ConvertToPhysicalUnits( nComp, iOS, u_A, Field )

    INTEGER,          INTENT(in)    :: nComp, iOS
    REAL(amrex_real), INTENT(inout) :: u_A(:,:,:,:)
    CHARACTER(2),     INTENT(in)    :: Field

    INTEGER          :: iComp
    REAL(amrex_real) :: unitsGF(nGF), unitsCF(nCF), unitsPF(nPF), unitsAF(nAF)

    CALL SetUnitsFields( unitsGF, unitsCF, unitsPF, unitsAF )

    IF     ( Field .EQ. 'GF' )THEN

      DO iComp = 1, nComp
        u_A(:,:,:,iComp+iOS) = u_A(:,:,:,iComp+iOS) / unitsGF(iComp)
      END DO

    ELSE IF( Field .EQ. 'CF' )THEN

      DO iComp = 1, nComp
        u_A(:,:,:,iComp+iOS) = u_A(:,:,:,iComp+iOS) / unitsCF(iComp)
      END DO

    ELSE IF( Field .EQ. 'PF' )THEN

      DO iComp = 1, nComp
        u_A(:,:,:,iComp+iOS) = u_A(:,:,:,iComp+iOS) / unitsPF(iComp)
      END DO

    ELSE IF( Field .EQ. 'AF' )THEN

      DO iComp = 1, nComp
        u_A(:,:,:,iComp+iOS) = u_A(:,:,:,iComp+iOS) / unitsAF(iComp)
      END DO

    ELSE

      RETURN

    END IF

  END SUBROUTINE ConvertToPhysicalUnits


  SUBROUTINE SetUnitsFields( unitsGF, unitsCF, unitsPF, unitsAF )

    REAL(amrex_real), INTENT(out) :: unitsGF(nGF)
    REAL(amrex_real), INTENT(out) :: unitsCF(nCF)
    REAL(amrex_real), INTENT(out) :: unitsPF(nPF)
    REAL(amrex_real), INTENT(out) :: unitsAF(nAF)

    ! --- Geometry ---

    unitsGF(iGF_Phi_N)    = Erg / Gram
    unitsGF(iGF_h_1)      = 1.0_amrex_real
    unitsGF(iGF_h_2)      = 1.0_amrex_real
    unitsGF(iGF_h_3)      = 1.0_amrex_real
    unitsGF(iGF_Gm_dd_11) = 1.0_amrex_real
    unitsGF(iGF_Gm_dd_22) = 1.0_amrex_real
    unitsGF(iGF_Gm_dd_33) = 1.0_amrex_real
    unitsGF(iGF_SqrtGm)   = 1.0_amrex_real
    unitsGF(iGF_Alpha)    = 1.0_amrex_real
    unitsGF(iGF_Beta_1)   = 1.0_amrex_real
    unitsGF(iGF_Beta_2)   = 1.0_amrex_real
    unitsGF(iGF_Beta_3)   = 1.0_amrex_real
    unitsGF(iGF_Psi)      = 1.0_amrex_real

    ! --- Conserved ---

    unitsCF(iCF_D)  = Gram / Centimeter**3
    unitsCF(iCF_S1) = Gram / Centimeter**2 / Second
    unitsCF(iCF_S2) = Gram / Centimeter**2 / Second
    unitsCF(iCF_S3) = Gram / Centimeter**2 / Second
    unitsCF(iCF_E)  = Erg / Centimeter**3
    unitsCF(iCF_Ne) = 1.0_amrex_real / Centimeter**3

    ! --- Primitive ---

    unitsPF(iPF_D)  = Gram / Centimeter**3
    unitsPF(iPF_V1) = Kilometer / Second
    unitsPF(iPF_V2) = Kilometer / Second
    unitsPF(iPF_V3) = Kilometer / Second
    unitsPF(iPF_E)  = Erg / Centimeter**3
    unitsPF(iPF_Ne) = 1.0_amrex_real / Centimeter**3

    ! --- Auxiliary ---

    unitsAF(iAF_P)  = Erg / Centimeter**3
    unitsAF(iAF_T)  = Kelvin
    unitsAF(iAF_Ye) = 1.0_amrex_real
    unitsAF(iAF_S)  = BoltzmannConstant
    unitsAF(iAF_E)  = Erg / Gram
    unitsAF(iAF_Me) = MeV
    unitsAF(iAF_Mp) = MeV
    unitsAF(iAF_Mn) = MeV
    unitsAF(iAF_Xp) = 1.0_amrex_real
    unitsAF(iAF_Xn) = 1.0_amrex_real
    unitsAF(iAF_Xa) = 1.0_amrex_real
    unitsAF(iAF_Xh) = 1.0_amrex_real
    unitsAF(iAF_Gm) = 1.0_amrex_real
    unitsAF(iAF_Cs) = Kilometer / Second

  END SUBROUTINE SetUnitsFields


END MODULE InputOutputModuleAMReX
