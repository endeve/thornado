MODULE InitializationModule

  USE ISO_C_BINDING

  ! --- AMReX Modules ---

  USE amrex_init_module, ONLY: &
    amrex_init
  USE amrex_parmparse_module, ONLY: &
    amrex_parmparse, &
    amrex_parmparse_build, &
    amrex_parmparse_destroy
  USE amrex_amrcore_module, ONLY: &
    amrex_amrcore_init, &
    amrex_init_virtual_functions, &
    amrex_init_from_scratch, &
    amrex_ref_ratio, &
    amrex_get_numlevels
  USE amrex_boxarray_module, ONLY: &
    amrex_boxarray
  USE amrex_distromap_module, ONLY: &
    amrex_distromap
  USE amrex_multifab_module, ONLY: &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy, &
    amrex_imultifab, &
    amrex_multifab, &
    amrex_multifab_build, &
    amrex_multifab_destroy
  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor, &
    amrex_parallel_myproc
  USE amrex_fluxregister_module, ONLY: &
    amrex_fluxregister_build, &
    amrex_fluxregister_destroy
  USE amrex_tagbox_module, ONLY: &
    amrex_tagboxarray

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    swX, &
    DescribeProgramHeaderX
  USE PolynomialBasisModule_Lagrange, ONLY: &
    InitializePolynomialBasis_Lagrange
  USE PolynomialBasisModule_Legendre, ONLY: &
    InitializePolynomialBasis_Legendre
  USE PolynomialBasisModuleX_Lagrange, ONLY: &
    InitializePolynomialBasisX_Lagrange
  USE PolynomialBasisModuleX_Legendre, ONLY: &
    InitializePolynomialBasisX_Legendre
  USE PolynomialBasisMappingModule, ONLY: &
    InitializePolynomialBasisMapping
  USE ReferenceElementModuleX, ONLY: &
    InitializeReferenceElementX, &
    nDOFX_X1, &
    WeightsX_q
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    InitializeReferenceElementX_Lagrange
  USE UnitsModule, ONLY: &
    DescribeUnitsDisplay
  USE MeshModule, ONLY: &
    MeshX, &
    MeshType
  USE GeometryFieldsModule, ONLY: &
    iGF_SqrtGm, &
    nGF
  USE FluidFieldsModule, ONLY: &
    unitsCF, &
    iCF_D, &
    nCF, &
    nPF, &
    nAF, &
    nDF

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero
  USE MaskModule, ONLY: &
    IsNotLeafElement, &
    CreateFineMask, &
    DestroyFineMask
  USE MF_FieldsModule_Geometry, ONLY: &
    CreateFields_Geometry_MF, &
    MF_uGF
  USE MF_FieldsModule_Euler, ONLY: &
    CreateFields_Euler_MF, &
    MF_uCF, &
    MF_uPF, &
    MF_uAF, &
    MF_uDF, &
    FluxRegister_Euler
  USE MF_EquationOfStateModule, ONLY: &
    InitializeEquationOfState_MF
  USE MF_Euler_SlopeLimiterModule, ONLY: &
    InitializeSlopeLimiter_Euler_MF, &
    ApplySlopeLimiter_Euler_MF
  USE MF_Euler_PositivityLimiterModule, ONLY: &
    InitializePositivityLimiter_Euler_MF, &
    ApplyPositivityLimiter_Euler_MF
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF
  USE FillPatchModule, ONLY: &
    FillPatch, &
    FillCoarsePatch
  USE RegridModule, ONLY: &
    ReGrid
  USE InputParsingModule, ONLY: &
    InitializeParameters, &
    nLevels, &
    nMaxLevels, &
    StepNo, &
    iRestart, &
    t_old, &
    t_new, &
    UseTiling, &
    UseFluxCorrection_Euler, &
    TagCriteria, &
    DescribeProgramHeader_AMReX
  USE InputOutputModuleAMReX, ONLY: &
    ReadCheckpointFile
  USE AverageDownModule, ONLY: &
    AverageDown
  USE Euler_MeshRefinementModule, ONLY: &
    InitializeMeshRefinement_Euler

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeProgram

CONTAINS


  SUBROUTINE InitializeProgram

    INTEGER :: iLevel
    TYPE(amrex_multifab) :: MF_dU

    CALL amrex_init()

    CALL amrex_amrcore_init()

    CALL InitializeParameters

    IF( amrex_parallel_ioprocessor() )THEN

      CALL DescribeUnitsDisplay
      CALL DescribeProgramHeaderX

    END IF

    CALL CreateFields_Geometry_MF
    CALL CreateFields_Euler_MF

    CALL InitializePolynomialBasisX_Lagrange
    CALL InitializePolynomialBasisX_Legendre

    CALL InitializePolynomialBasis_Lagrange
    CALL InitializePolynomialBasis_Legendre

    CALL CreateMesh_MF( 0, MeshX )

    CALL InitializePolynomialBasisMapping &
           ( [Zero], MeshX(1) % Nodes, MeshX(2) % Nodes, MeshX(3) % Nodes )

    CALL DestroyMesh_MF( MeshX )

    ! --- Ordering of calls is important here ---
    CALL InitializeReferenceElementX
    CALL InitializeReferenceElementX_Lagrange

    CALL InitializeMeshRefinement_Euler

    CALL InitializeEquationOfState_MF

    CALL InitializePositivityLimiter_Euler_MF

    CALL InitializeSlopeLimiter_Euler_MF

    CALL amrex_init_virtual_functions &
           ( MakeNewLevelFromScratch, &
             MakeNewLevelFromCoarse, &
             RemakeLevel, &
             ClearLevel, &
             ErrorEstimate )

    ALLOCATE( StepNo(0:nMaxLevels-1) )
    ALLOCATE( t_old (0:nMaxLevels-1) )
    ALLOCATE( t_new (0:nMaxLevels-1) )

    StepNo = 0
    t_new  = 0.0_DP

    IF( iRestart .LT. 0 )THEN

      CALL amrex_init_from_scratch( 0.0_DP )
      nLevels = amrex_get_numlevels()

      CALL AverageDown( MF_uGF )
      CALL AverageDown( MF_uGF, MF_uCF )

      CALL ApplySlopeLimiter_Euler_MF &
             ( MF_uGF, MF_uCF, MF_uDF )

      CALL ApplyPositivityLimiter_Euler_MF &
             ( MF_uGF, MF_uCF, MF_uDF )

    ELSE

      CALL ReadCheckpointFile( ReadFields_uCF_Option = .TRUE. )

    END IF

    CALL AverageDown( MF_uGF )
    CALL AverageDown( MF_uGF, MF_uCF )
    CALL ApplyPositivityLimiter_Euler_MF &
           ( MF_uGF, MF_uCF, MF_uDF )

    DO iLevel = 0, nLevels-1

      CALL amrex_multifab_build &
             ( MF_dU, MF_uCF(iLevel) % BA, MF_uCF(iLevel) % DM, &
               nCF * nDOFX, swX )

      CALL MF_dU % SetVal( 1.0_DP )

      CALL MF_uCF(iLevel) % Add( MF_dU, 1, 1, nCF * nDOFX, swX )

      CALL amrex_multifab_destroy( MF_dU )

    END DO

    CALL ReGrid

    t_old = t_new

    CALL DescribeProgramHeader_AMReX

    CALL WriteCellAverageData

  END SUBROUTINE InitializeProgram


  SUBROUTINE WriteCellAverageData

    INTEGER  :: iX1, iX2, iX3, iCF, iX_B0(3), iX_E0(3), iLevel
    REAL(DP) :: Volume, SqrtGm(nDOFX), UN(nDOFX)

    TYPE(amrex_box)      :: BX
    TYPE(amrex_mfiter)   :: MFI
    TYPE(amrex_multifab) :: MF_uCF_K(0:nMaxLevels-1)

    REAL(DP), CONTIGUOUS, POINTER :: G  (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: U  (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: U_K(:,:,:,:)

    DO iLevel = 0, nLevels-1

      CALL amrex_multifab_build &
             ( MF_uCF_K(iLevel), MF_uCF(iLevel) % BA, MF_uCF(iLevel) % DM, &
               nCF, swX )

      CALL amrex_mfiter_build( MFI, MF_uCF_K(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        G   => MF_uGF  (iLevel) % DataPtr( MFI )
        U   => MF_uCF  (iLevel) % DataPtr( MFI )
        U_K => MF_uCF_K(iLevel) % DataPtr( MFI )

        BX = MFI % TileBox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi

        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)

          SqrtGm &
            = G(iX1,iX2,iX3,nDOFX*(iGF_SqrtGm-1)+1:nDOFX*(iGF_SqrtGm-1)+nDOFX)

          Volume = SUM( WeightsX_q * SqrtGm )

          DO iCF = 1, nCF

            UN = U(iX1,iX2,iX3,nDOFX*(iCF-1)+1:nDOFX*(iCF-1)+nDOFX)

            U_K(iX1,iX2,iX3,iCF) &
              = SUM( WeightsX_q * UN * SqrtGm ) / Volume

            U_K(iX1,iX2,iX3,iCF) = U_K(iX1,iX2,iX3,iCF) / unitsCF(iCF)

          END DO

        END DO
        END DO
        END DO

      END DO

      CALL amrex_mfiter_destroy( MFI )

    END DO

    CALL ShowVariableFromMultiFab_Vector &
           ( MF_uCF_K, iCF_D, WriteToFile_Option = .TRUE., &
             FileNameBase_Option = 'CF_D' )

    DO iLevel = 0, nLevels-1

      CALL amrex_multifab_destroy( MF_uCF_K(iLevel) )

    END DO

  END SUBROUTINE WriteCellAverageData


  SUBROUTINE ShowVariableFromMultiFab_Single &
    ( iLevel, MF, iField, iMF_FineMask_Option, &
      swXX_Option, WriteToFile_Option, FileNameBase_Option )

    INTEGER              , INTENT(in) :: iLevel, iField
    TYPE(amrex_multifab) , INTENT(in) :: MF
    TYPE(amrex_imultifab), INTENT(in), OPTIONAL :: iMF_FineMask_Option
    INTEGER              , INTENT(in), OPTIONAL :: swXX_Option(3)
    LOGICAL              , INTENT(in), OPTIONAL :: WriteToFile_Option
    CHARACTER(*)         , INTENT(in), OPTIONAL :: FileNameBase_Option

    INTEGER                       :: iX1, iX2, iX3
    INTEGER                       :: lo(4), hi(4), iX_B1(3), iX_E1(3)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(DP), CONTIGUOUS, POINTER :: F       (:,:,:,:)
    INTEGER , CONTIGUOUS, POINTER :: FineMask(:,:,:,:)
    INTEGER                       :: swXX(3)
    INTEGER                       :: iFileNo
    LOGICAL                       :: WriteToFile
    CHARACTER(128)                :: FMT
    CHARACTER(128)                :: FileNameBase, FileName

    TYPE(MeshType) :: MeshXX(3)

    swXX = 0
    IF( PRESENT( swXX_Option ) ) swXX = swXX_Option

    WriteToFile = .FALSE.
    IF( PRESENT( WriteToFile_Option ) ) WriteToFile = WriteToFile_Option

    WRITE(FMT,'(A,I2.2,A,I2.2,A,I2.2,A,I3.3,A)') &
      '(I2.2,3I8.6,SP3ES25.16E3,SP', &
      1, 'ES25.16E3,SP', &
      1, 'ES25.16E3,SP', &
      1, 'ES25.16E3,SP', &
      1, 'ES25.16E3)'

    CALL amrex_mfiter_build( MFI, MF, tiling = UseTiling )

    CALL CreateMesh_MF( iLevel, MeshXX )

    DO WHILE( MFI % next() )

      IF( WriteToFile )THEN

        iFileNo = 100 + amrex_parallel_myproc()

        IF( PRESENT( FileNameBase_Option ) ) &
          FileNameBase = TRIM( FileNameBase_Option )

        WRITE(FileName,'(A,A,I3.3,A,I8.8,A)') &
          TRIM( FileNameBase ), '_proc', &
          amrex_parallel_myproc(), '_', StepNo(0), '.dat'

        OPEN( iFileNo, FILE = TRIM( FileName ), POSITION = 'APPEND' )

      END IF

      IF( PRESENT( iMF_FineMask_Option ) ) &
        FineMask => iMF_FineMask_Option % DataPtr( MFI )

      F => MF % DataPtr( MFI )

      BX = MFI % tilebox()

      lo = LBOUND( F ); hi = UBOUND( F )

      iX_B1 = BX % lo - swXX
      iX_E1 = BX % hi + swXX

      DO iX3 = iX_B1(3), iX_E1(3)
      DO iX2 = iX_B1(2), iX_E1(2)
      DO iX1 = iX_B1(1), iX_E1(1)

!        IF( PRESENT( iMF_FineMask_Option ) )THEN
!
!          IF( IsNotLeafElement( FineMask(iX1,iX2,iX3,1) ) ) CYCLE
!
!        END IF

        IF( WriteToFile )THEN

          WRITE(iFileNo,TRIM(FMT)) &
            iLevel, iX1, iX2, iX3, &
            MeshXX(1) % Width (iX1), &
            MeshXX(2) % Width (iX2), &
            MeshXX(3) % Width (iX3), &
            MeshXX(1) % Center(iX1), &
            MeshXX(2) % Center(iX2), &
            MeshXX(3) % Center(iX3), &
            F(iX1,iX2,iX3,iField)

        ELSE

          WRITE(*,TRIM(FMT)) &
            iLevel, iX1, iX2, iX3, &
            F(iX1,iX2,iX3,iField)

        END IF

      END DO
      END DO
      END DO

      IF( WriteToFile )THEN

        WRITE(iFileNo,*)

        CLOSE( iFileNo )

      ELSE IF( ANY( FineMask(:,:,:,1) .EQ. 0 ) )THEN

        WRITE(*,*)

      END IF

    END DO ! WHILE( MFI % next() )

    CALL amrex_mfiter_destroy( MFI )

    CALL DestroyMesh_MF( MeshXX )

    IF( .NOT. WriteToFile ) WRITE(*,*)

  END SUBROUTINE ShowVariableFromMultiFab_Single


  SUBROUTINE ShowVariableFromMultiFab_Vector &
    ( MF, iField, swXX_Option, WriteToFile_Option, FileNameBase_Option )

    INTEGER             , INTENT(in) :: iField
    TYPE(amrex_multifab), INTENT(in) :: MF(0:)
    INTEGER             , INTENT(in), OPTIONAL :: swXX_Option(3)
    LOGICAL             , INTENT(in), OPTIONAL :: WriteToFile_Option
    CHARACTER(*)        , INTENT(in), OPTIONAL :: FileNameBase_Option

    INTEGER :: iLevel

    INTEGER        :: swXX(3)
    LOGICAL        :: WriteToFile
    CHARACTER(128) :: FileNameBase

    TYPE(amrex_imultifab) :: iMF_FineMask

    swXX = 0
    IF( PRESENT( swXX_Option ) ) swXX = swXX_Option

    WriteToFile = .FALSE.
    IF( PRESENT( WriteToFile_Option ) ) WriteToFile = WriteToFile_Option

    FileNameBase = ''
    IF( PRESENT( FileNameBase_Option ) ) &
      FileNameBase = TRIM( FileNameBase_Option )

    DO iLevel = 0, nLevels-1

      CALL CreateFineMask( iLevel, iMF_FineMask, MF % BA, MF % DM )

      CALL ShowVariableFromMultiFab_Single &
             ( iLevel, MF(iLevel), iField, &
               iMF_FineMask_Option = iMF_FineMask, &
               swXX_Option = swXX, &
               WriteToFile_Option = WriteToFile, &
               FileNameBase_Option = TRIM( FileNameBase ) )

      CALL DestroyFineMask( iMF_FineMask )

    END DO

  END SUBROUTINE ShowVariableFromMultiFab_Vector


  SUBROUTINE MakeNewLevelFromScratch( iLevel, Time, pBA, pDM ) BIND(c)

    USE MF_GeometryModule, ONLY: &
      ComputeGeometryX_MF

    USE MF_InitializationModule, ONLY: &
      InitializeFields_MF

    INTEGER,     INTENT(in), VALUE :: iLevel
    REAL(DP),    INTENT(in), VALUE :: Time
    TYPE(c_ptr), INTENT(in), VALUE :: pBA, pDM

    TYPE(amrex_boxarray)  :: BA
    TYPE(amrex_distromap) :: DM

    BA = pBA
    DM = pDM

    t_new(iLevel) = Time
    t_old(iLevel) = Time - 1.0e200_DP

    CALL ClearLevel( iLevel )

    CALL amrex_multifab_build( MF_uGF(iLevel), BA, DM, nDOFX * nGF, swX )
    CALL MF_uGF(iLevel) % SetVal( Zero )

    CALL amrex_multifab_build( MF_uCF(iLevel), BA, DM, nDOFX * nCF, swX )
    CALL MF_uCF(iLevel) % SetVal( Zero )

    CALL amrex_multifab_build( MF_uPF(iLevel), BA, DM, nDOFX * nPF, swX )
    CALL MF_uPF(iLevel) % SetVal( Zero )

    CALL amrex_multifab_build( MF_uAF(iLevel), BA, DM, nDOFX * nAF, swX )
    CALL MF_uAF(iLevel) % SetVal( Zero )

    CALL amrex_multifab_build( MF_uDF(iLevel), BA, DM, nDOFX * nDF, swX )
    CALL MF_uDF(iLevel) % SetVal( Zero )

    ! Assume nDOFX_X2 = nDOFX_X3 = nDOFX_X1
    IF( iLevel .GT. 0 .AND. UseFluxCorrection_Euler ) &
      CALL amrex_fluxregister_build &
             ( FluxRegister_Euler(iLevel), BA, DM, &
               amrex_ref_ratio(iLevel-1), iLevel, nDOFX_X1*nCF )

    CALL CreateMesh_MF( iLevel, MeshX )

    CALL ComputeGeometryX_MF( MF_uGF(iLevel) )

    CALL InitializeFields_MF( iLevel, MF_uGF(iLevel), MF_uCF(iLevel) )

    CALL DestroyMesh_MF( MeshX )

  END SUBROUTINE MakeNewLevelFromScratch


  SUBROUTINE MakeNewLevelFromCoarse( iLevel, Time, pBA, pDM ) BIND(c)

    INTEGER,     INTENT(in), VALUE :: iLevel
    REAL(DP),    INTENT(in), VALUE :: Time
    TYPE(c_ptr), INTENT(in), VALUE :: pBA, pDM

    TYPE(amrex_boxarray)  :: BA
    TYPE(amrex_distromap) :: DM

    BA = pBA
    DM = pDM

    CALL ClearLevel( iLevel )

    t_new( iLevel ) = Time
    t_old( iLevel ) = Time - 1.0e200_DP

    CALL amrex_multifab_build( MF_uGF(iLevel), BA, DM, nDOFX * nGF, swX )
    CALL amrex_multifab_build( MF_uCF(iLevel), BA, DM, nDOFX * nCF, swX )
    CALL amrex_multifab_build( MF_uPF(iLevel), BA, DM, nDOFX * nPF, swX )
    CALL amrex_multifab_build( MF_uAF(iLevel), BA, DM, nDOFX * nAF, swX )
    CALL amrex_multifab_build( MF_uDF(iLevel), BA, DM, nDOFX * nDF, swX )

    IF( iLevel .GT. 0 .AND. UseFluxCorrection_Euler ) &
      CALL amrex_fluxregister_build &
             ( FluxRegister_Euler(iLevel), BA, DM, amrex_ref_ratio(iLevel-1), &
               iLevel, nDOFX_X1 * nCF )

    CALL FillCoarsePatch( iLevel, MF_uGF, &
                          ApplyBoundaryConditions_Geometry_Option = .TRUE. )

    CALL FillCoarsePatch( iLevel, MF_uDF )

    CALL FillCoarsePatch( iLevel, MF_uGF, MF_uCF, &
                          ApplyBoundaryConditions_Euler_Option = .TRUE. )

    CALL ApplyPositivityLimiter_Euler_MF &
           ( iLevel, MF_uGF(iLevel), MF_uCF(iLevel), MF_uDF(iLevel) )

  END SUBROUTINE MakeNewLevelFromCoarse


  SUBROUTINE ClearLevel( iLevel ) BIND(c)

    INTEGER, INTENT(in), VALUE :: iLevel

    CALL amrex_multifab_destroy( MF_uDF(iLevel) )
    CALL amrex_multifab_destroy( MF_uAF(iLevel) )
    CALL amrex_multifab_destroy( MF_uPF(iLevel) )
    CALL amrex_multifab_destroy( MF_uCF(iLevel) )
    CALL amrex_multifab_destroy( MF_uGF(iLevel) )

    IF( iLevel .GT. 0 .AND. UseFluxCorrection_Euler ) &
      CALL amrex_fluxregister_destroy( FluxRegister_Euler(iLevel) )

  END SUBROUTINE ClearLevel


  SUBROUTINE RemakeLevel( iLevel, Time, pBA, pDM ) BIND(c)

    INTEGER,     INTENT(in), VALUE :: iLevel
    REAL(DP),    INTENT(in), VALUE :: Time
    TYPE(c_ptr), INTENT(in), VALUE :: pBA, pDM

    TYPE(amrex_boxarray)  :: BA
    TYPE(amrex_distromap) :: DM
    TYPE(amrex_multifab)  :: MF_uGF_tmp, MF_uCF_tmp, MF_uPF_tmp, &
                             MF_uAF_tmp, MF_uDF_tmp

    BA = pBA
    DM = pDM

    CALL amrex_multifab_build( MF_uGF_tmp, BA, DM, nDOFX * nGF, swX )
    CALL amrex_multifab_build( MF_uCF_tmp, BA, DM, nDOFX * nCF, swX )
    CALL amrex_multifab_build( MF_uPF_tmp, BA, DM, nDOFX * nPF, swX )
    CALL amrex_multifab_build( MF_uAF_tmp, BA, DM, nDOFX * nAF, swX )
    CALL amrex_multifab_build( MF_uDF_tmp, BA, DM, nDOFX * nDF, swX )

    CALL FillPatch( iLevel, MF_uGF, MF_uGF_tmp, &
                    ApplyBoundaryConditions_Geometry_Option = .TRUE. )

    CALL FillPatch( iLevel, MF_uDF, MF_uDF_tmp )

    CALL FillPatch( iLevel, MF_uGF, MF_uGF_tmp, MF_uCF, MF_uCF_tmp, &
                    ApplyBoundaryConditions_Euler_Option = .TRUE. )

    CALL ApplyPositivityLimiter_Euler_MF &
           ( iLevel, MF_uGF(iLevel), MF_uCF(iLevel), MF_uDF(iLevel) )

    CALL ClearLevel( iLevel )

    CALL amrex_multifab_build( MF_uGF(iLevel), BA, DM, nDOFX * nGF, swX )
    CALL amrex_multifab_build( MF_uCF(iLevel), BA, DM, nDOFX * nCF, swX )
    CALL amrex_multifab_build( MF_uPF(iLevel), BA, DM, nDOFX * nPF, swX )
    CALL amrex_multifab_build( MF_uAF(iLevel), BA, DM, nDOFX * nAF, swX )
    CALL amrex_multifab_build( MF_uDF(iLevel), BA, DM, nDOFX * nDF, swX )

    IF( iLevel .GT. 0 .AND. UseFluxCorrection_Euler ) &
      CALL amrex_fluxregister_build &
             ( FluxRegister_Euler(iLevel), BA, DM, amrex_ref_ratio(iLevel-1), &
               iLevel, nDOFX_X1 * nCF )

    CALL MF_uGF(iLevel) % COPY( MF_uGF_tmp, 1, 1, nDOFX * nGF, swX )
    CALL MF_uCF(iLevel) % COPY( MF_uCF_tmp, 1, 1, nDOFX * nCF, swX )
    CALL MF_uDF(iLevel) % COPY( MF_uDF_tmp, 1, 1, nDOFX * nDF, swX )

    CALL amrex_multifab_destroy( MF_uDF_tmp )
    CALL amrex_multifab_destroy( MF_uAF_tmp )
    CALL amrex_multifab_destroy( MF_uPF_tmp )
    CALL amrex_multifab_destroy( MF_uCF_tmp )
    CALL amrex_multifab_destroy( MF_uGF_tmp )

  END SUBROUTINE RemakeLevel


  SUBROUTINE ErrorEstimate( iLevel, cp, Time, SetTag, ClearTag ) BIND(c)

    USE TaggingModule, ONLY: &
      TagElements

    INTEGER,                INTENT(in), VALUE :: iLevel
    TYPE(c_ptr),            INTENT(in), VALUE :: cp
    REAL(DP),               INTENT(in), VALUE :: Time
    CHARACTER(KIND=c_char), INTENT(in), VALUE :: SetTag, ClearTag

    TYPE(amrex_parmparse)   :: PP
    TYPE(amrex_tagboxarray) :: Tag
    TYPE(amrex_mfiter)      :: MFI
    TYPE(amrex_box)         :: BX
    REAL(DP),               CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    CHARACTER(KIND=c_char), CONTIGUOUS, POINTER :: TagArr(:,:,:,:)

    IF( .NOT. ALLOCATED( TagCriteria ) )THEN

       CALL amrex_parmparse_build( PP, "amr" )

         CALL PP % getarr( "TagCriteria", TagCriteria )

       CALL amrex_parmparse_destroy( PP )

    END IF

    Tag = cp

    CALL CreateMesh_MF( iLevel, MeshX )

    !$OMP PARALLEL PRIVATE( MFI, BX, uCF, TagArr )
    CALL amrex_mfiter_build( MFI, MF_uCF( iLevel ), Tiling = UseTiling )

    DO WHILE( MFI % next() )

      BX = MFI % TileBox()

      uCF    => MF_uCF( iLevel ) % DataPtr( MFI )
      TagArr => Tag              % DataPtr( MFI )

      ! TagCriteria(iLevel+1) because iLevel starts at 0 but
      ! TagCriteria starts with 1

      CALL TagElements &
             ( iLevel, BX % lo, BX % hi, LBOUND( uCF ), UBOUND( uCF ), &
               uCF, TagCriteria(iLevel+1), SetTag, ClearTag, &
               LBOUND( TagArr ), UBOUND( TagArr ), TagArr )

    END DO

    CALL amrex_mfiter_destroy( MFI )
    !$OMP END PARALLEL

    CALL DestroyMesh_MF( MeshX )

  END SUBROUTINE ErrorEstimate


END MODULE InitializationModule
