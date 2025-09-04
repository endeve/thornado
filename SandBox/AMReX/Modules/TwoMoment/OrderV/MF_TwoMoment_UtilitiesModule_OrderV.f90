MODULE MF_TwoMoment_UtilitiesModule

  ! --- AMReX Modules ---

  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_reduce_min, &
    amrex_parallel_ioprocessor

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    swX, &
    swE, &
    nDOFX, &
    nDOFZ, &
    nDOFE, &
    iE_B0, &
    iE_E0, &
    iE_B1, &
    iE_E1
  USE RadiationFieldsModule, ONLY: &
    nCR, &
    iCR_N, &
    iCR_G1, &
    iCR_G2, &
    iCR_G3, &
    nPR, &
    nAR, &
    iPR_D, &
    iPR_I1,  &
    iPR_I2,  &
    iPR_I3,  &
    iGR_N,   &      
    iGR_D,   & 
    iGR_I1,  & 
    iGR_I2,  &
    iGR_I3,  & 
    iGR_J,   & 
    iGR_H1,  & 
    iGR_H2,  & 
    iGR_H3,  & 
    iGR_RMS, &
    iGR_F,   & 
    iGR_K,   & 
    iGR_Q,   & 
    nGR,     & 
    LeptonNumber
  USE FluidFieldsModule, ONLY: &
    nCF, &
    iCF_D, &
    iCF_S1, &
    iCF_S2, &
    iCF_S3, &
    iCF_E, &
    iCF_Ne, &
    nPF, &
    iPF_D, &
    iPF_V1, &
    iPF_V2, &
    iPF_V3, &
    iPF_E, &
    iPF_Ne, &
    nAF, &
    nDF
  USE GeometryFieldsModuleE, ONLY: &
    uGE
  USE GeometryFieldsModule, ONLY: &
    nGF, &
    iGF_h_1, &
    iGF_h_2, &
    iGF_h_3, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33
  USE Euler_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_Euler
  USE TwoMoment_UtilitiesModule, ONLY: &
    ComputeFromConserved_TwoMoment, &
    ComputeTimeStep_TwoMoment, &
    ComputeTimeStep_TwoMoment_Realizability
  USE Euler_UtilitiesModule_NonRelativistic, ONLY: &
    ComputePrimitive_Euler_NonRelativistic, &
    ComputeFromConserved_Euler_NonRelativistic
  USE Euler_ErrorModule, ONLY: &
    DescribeError_Euler
  USE MeshModule, ONLY: &
    MeshX, &
    MeshE, &
    NodeCoordinate
  USE UnitsModule, ONLY: &
    UnitsActive, &
    AtomicMassUnit, &
    SpeedOfLight, &
    PlanckConstant
  USE TwoMoment_ClosureModule, ONLY: &
    FluxFactor, &
    EddingtonFactor, &
    HeatFluxFactor

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero, &
    One, &
    Two, &
    FourPi
  USE InputParsingModule, ONLY: &
    nLevels, &
    nSpecies, &
    nE, &
    UseTiling
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF
  USE MF_UtilitiesModule, ONLY: &
    amrex2thornado_X, &
    thornado2amrex_X, &
    amrex2thornado_Z, &
    thornado2amrex_Z, &
    amrex2thornado_Integrated, &
    thornado2amrex_Integrated, &
    AllocateArray_X, &
    DeallocateArray_X, &
    AllocateArray_Z, &
    DeallocateArray_Z, &
    AllocateArray_Integrated, &
    DeallocateArray_Integrated
  IMPLICIT NONE
  PRIVATE


  PUBLIC :: ComputeFromConserved_TwoMoment_MF
  PUBLIC :: ComputeTimeStep_TwoMoment_MF
  PUBLIC :: ComputeTimeStep_TwoMoment_Realizability_MF


CONTAINS

  SUBROUTINE ComputeTimeStep_TwoMoment_MF( MF_uGF, CFL, TimeStep )

    TYPE(amrex_multifab),  INTENT(in)  :: MF_uGF(0:nLevels-1)
    REAL(DP)            ,  INTENT(in)  :: CFL
    REAL(DP)            ,  INTENT(out) :: TimeStep(0:nLevels-1)


    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)

    REAL(DP), ALLOCATABLE :: G(:,:,:,:,:)

    REAL(DP) :: myTimeStep(0:nLevels-1)
    INTEGER  :: iLevel, iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    INTEGER  :: iLo_MF(4)

    myTimeStep = HUGE( One )
    TimeStep = HUGE( One )

    DO iLevel = 0, nLevels-1

      CALL CreateMesh_MF( iLevel, MeshX )

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )

        iLo_MF = LBOUND( uGF )

        BX = MFI % tilebox()

        iX_B0(1:3) = BX % lo(1:3)
        iX_E0(1:3) = BX % hi(1:3)
        iX_B1(1:3) = BX % lo(1:3) - swX(1:3)
        iX_E1(1:3) = BX % hi(1:3) + swX(1:3)

        CALL AllocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
                 G )

        CALL amrex2thornado_X &
               ( nGF, iX_B1, iX_E1, iLo_MF, iX_B0, iX_E0, uGF, G )

        !PRINT *, G

        CALL ComputeTimeStep_TwoMoment &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G, CFL, myTimeStep( iLevel ) )
        
        !PRINT *, 'TIME 1 STEP HERE:'
        !PRINT *, myTimeStep(iLevel)

        TimeStep( iLevel ) = MIN( TimeStep( iLevel ), myTimeStep( iLevel ) )
        !PRINT *, 'TIME STEP HERE:'
        !PRINT *, TimeStep(iLevel)
        !PRINT *, myTimeStep(iLevel)
        CALL DeallocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
                 G )

      END DO ! --- Loop over grids (boxes) ---

      CALL amrex_mfiter_destroy( MFI )

      CALL DestroyMesh_MF( MeshX )

    END DO ! --- Loop over levels ---

    CALL amrex_parallel_reduce_min( TimeStep, SIZE(TimeStep) )

  END SUBROUTINE ComputeTimeStep_TwoMoment_MF


  SUBROUTINE ComputeTimeStep_TwoMoment_Realizability_MF( MF_uGF, MF_uCF, CFL, TimeStep, Verbose_Option)

    TYPE(amrex_multifab),  INTENT(in)  :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab),  INTENT(in)  :: MF_uCF(0:nLevels-1)
    REAL(DP)            ,  INTENT(in)  :: CFL
    REAL(DP)            ,  INTENT(out) :: TimeStep(0:nLevels-1)


    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    REAL(DP), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: CF(:,:,:,:,:)

    REAL(DP) :: myTimeStep(0:nLevels-1)
    INTEGER  :: iLevel, iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    INTEGER  :: iLo_MF(4), iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    LOGICAL, INTENT(in), OPTIONAL :: Verbose_Option

    myTimeStep = HUGE( One )
    TimeStep = HUGE( One )

!PRINT *, 'ENTERING LOOP'

    DO iLevel = 0, nLevels-1

!PRINT *, 'Create MEsh'
      CALL CreateMesh_MF( iLevel, MeshX )

!PRINT *, 'CALL ITERATOR'
      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCF => MF_uCF(iLevel) % DataPtr( MFI )

        iLo_MF = LBOUND( uGF )

        BX = MFI % tilebox()

        iX_B0(1:3) = BX % lo(1:3)
        iX_E0(1:3) = BX % hi(1:3)
        iX_B1(1:3) = BX % lo(1:3) - swX(1:3)
        iX_E1(1:3) = BX % hi(1:3) + swX(1:3)
      
        iZ_B0(1) = iE_B0
        iZ_E0(1) = iE_E0
        iZ_B1(1) = iE_B1
        iZ_E1(1) = iE_E1

        iZ_B0(2:4) = iX_B0
        iZ_E0(2:4) = iX_E0
        iZ_B1(2:4) = iX_B1
        iZ_E1(2:4) = iX_E1

        CALL AllocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
                 G )

        CALL amrex2thornado_X &
               ( nGF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uGF, G )

        CALL AllocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCF ], &
                 CF )

        CALL amrex2thornado_X &
               ( nCF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uCF, CF )

        CALL ComputeTimeStep_TwoMoment_Realizability &
                ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, G, CF, CFL, myTimeStep( iLevel ), Verbose_Option=.FALSE. )

        TimeStep( iLevel ) = MIN( TimeStep( iLevel ), myTimeStep( iLevel ) )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
                 G )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCF ], &
                 CF )

      END DO ! --- Loop over grids (boxes) ---

      CALL amrex_mfiter_destroy( MFI )

      CALL DestroyMesh_MF( MeshX )

    END DO ! --- Loop over levels ---

    CALL amrex_parallel_reduce_min( TimeStep, SIZE(TimeStep) )

  END SUBROUTINE ComputeTimeStep_TwoMoment_Realizability_MF

  SUBROUTINE ComputeFromConserved_TwoMoment_MF &
    ( MF_uGF, MF_uCF, MF_uCR, MF_uPR, MF_uAR, MF_uGR )


    TYPE(amrex_multifab), INTENT(in)    :: &
      MF_uGF(0:nLevels-1), MF_uCF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in)    :: &
      MF_uCR(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout)    :: &
      MF_uAR(0:nLevels-1), MF_uPR(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: &
      MF_uGR(0:nLevels-1)

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCR(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uPR(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uGR(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uAR(:,:,:,:)

    REAL(DP), ALLOCATABLE :: G (:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: CF(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: CR(:,:,:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: PR(:,:,:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: AR(:,:,:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: GR(:,:,:,:,:,:)


    INTEGER :: iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    INTEGER :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    INTEGER :: iLevel, iLo_MF(4)

    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )


        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCF => MF_uCF(iLevel) % DataPtr( MFI )
        uCR => MF_uCR(iLevel) % DataPtr( MFI )
        uPR => MF_uPR(iLevel) % DataPtr( MFI )
        uGR => MF_uGR(iLevel) % DataPtr( MFI )
        uAR => MF_uAR(iLevel) % DataPtr( MFI )

        iLo_MF = LBOUND( uGF )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = BX % lo - swX
        iX_E1 = BX % hi + swX

        iZ_B0(1) = iE_B0
        iZ_E0(1) = iE_E0
        iZ_B1(1) = iE_B1
        iZ_E1(1) = iE_E1

        iZ_B0(2:4) = iX_B0
        iZ_E0(2:4) = iX_E0
        iZ_B1(2:4) = iX_B1
        iZ_E1(2:4) = iX_E1

        CALL AllocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
                 G )
        CALL AllocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCF ], &
                 CF )
        CALL AllocateArray_Z &
               ( [ 1       , &
                   iZ_B1(1), &
                   iZ_B1(2), &
                   iZ_B1(3), &
                   iZ_B1(4), &
                   1       , &
                   1        ], &
                 [ nDOFZ   , &
                   iZ_E1(1), &
                   iZ_E1(2), &
                   iZ_E1(3), &
                   iZ_E1(4), &
                   nCR     , &
                   nSpecies ], &
                 CR )

        CALL AllocateArray_Z &
               ( [ 1       , &
                   iZ_B1(1), &
                   iZ_B1(2), &
                   iZ_B1(3), &
                   iZ_B1(4), &
                   1       , &
                   1        ], &
                 [ nDOFZ   , &
                   iZ_E1(1), &
                   iZ_E1(2), &
                   iZ_E1(3), &
                   iZ_E1(4), &
                   nPR     , &
                   nSpecies ], &
                 PR )

        CALL AllocateArray_Z &
               ( [ 1       , &
                   iZ_B1(1), &
                   iZ_B1(2), &
                   iZ_B1(3), &
                   iZ_B1(4), &
                   1       , &
                   1        ], &
                 [ nDOFZ   , &
                   iZ_E1(1), &
                   iZ_E1(2), &
                   iZ_E1(3), &
                   iZ_E1(4), &
                   nAR     , &
                   nSpecies ], &
                 AR )

        CALL AllocateArray_Integrated &
               ( [ 1       , &
                   iX_B1(1), &
                   iX_B1(2), &
                   iX_B1(3), &
                   1       , &
                   1        ], &
                 [ nDOFX   , &
                   iX_E1(1), &
                   iX_E1(2), &
                   iX_E1(3), &
                   nGR     , &
                   nSpecies ], &
                 GR )

        CALL amrex2thornado_X &
               ( nGF, iX_B1, iX_E1, iLo_MF, iX_B0, iX_E0, uGF, G )

        CALL amrex2thornado_X &
               ( nCF, iX_B1, iX_E1, iLo_MF, iX_B0, iX_E0, uCF, CF )

        CALL amrex2thornado_Z &
               ( nCR, nSpecies, nE, iE_B0, iE_E0, &
                 iZ_B1, iZ_E1, iLo_MF, iZ_B0, iZ_E0, uCR, CR )

        CALL amrex2thornado_Z &
               ( nPR, nSpecies, nE, iE_B0, iE_E0, &
                 iZ_B1, iZ_E1, iLo_MF, iZ_B0, iZ_E0, uPR, PR )

        CALL amrex2thornado_Z &
               ( nAR, nSpecies, nE, iE_B0, iE_E0, &
                 iZ_B1, iZ_E1, iLo_MF, iZ_B0, iZ_E0, uAR, AR )

        CALL amrex2thornado_Integrated &
               ( nGR, nSpecies, iX_B1, iX_E1, &
                 iLo_MF, iX_B0, iX_E0, uGR, GR )

        CALL ComputeFromConserved_TwoMoment &
              ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, G, CF, CR, PR, AR, GR )

        CALL thornado2amrex_Integrated &
               ( nGR, nSpecies, iX_B1, iX_E1, &
                 iLo_MF, iX_B0, iX_E0, uGR, GR )

        CALL thornado2amrex_Z &
               ( nCR, nSpecies, nE, iE_B0, iE_E0, &
                 iZ_B1, iZ_E1, iLo_MF, iZ_B0, iZ_E0, uCR, CR )

        CALL thornado2amrex_Z &
               ( nPR, nSpecies, nE, iE_B0, iE_E0, &
                 iZ_B1, iZ_E1, iLo_MF, iZ_B0, iZ_E0, uPR, PR )

        CALL thornado2amrex_Z &
               ( nAR, nSpecies, nE, iE_B0, iE_E0, &
                 iZ_B1, iZ_E1, iLo_MF, iZ_B0, iZ_E0, uAR, AR )


        CALL DeallocateArray_Integrated &
               ( [ 1       , &
                   iX_B1(1), &
                   iX_B1(2), &
                   iX_B1(3), &
                   1       , &
                   1        ], &
                 [ nDOFX   , &
                   iX_E1(1), &
                   iX_E1(2), &
                   iX_E1(3), &
                   nGR     , &
                   nSpecies ], &
                 GR )

        CALL DeallocateArray_Z &
               ( [ 1       , &
                   iZ_B1(1), &
                   iZ_B1(2), &
                   iZ_B1(3), &
                   iZ_B1(4), &
                   1       , &
                   1        ], &
                 [ nDOFZ   , &
                   iZ_E1(1), &
                   iZ_E1(2), &
                   iZ_E1(3), &
                   iZ_E1(4), &
                   3       , &
                   nSpecies ], &
                 AR )

        CALL DeallocateArray_Z &
               ( [ 1       , &
                   iZ_B1(1), &
                   iZ_B1(2), &
                   iZ_B1(3), &
                   iZ_B1(4), &
                   1       , &
                   1        ], &
                 [ nDOFZ   , &
                   iZ_E1(1), &
                   iZ_E1(2), &
                   iZ_E1(3), &
                   iZ_E1(4), &
                   nPR     , &
                   nSpecies ], &
                 PR )

        CALL DeallocateArray_Z &
               ( [ 1       , &
                   iZ_B1(1), &
                   iZ_B1(2), &
                   iZ_B1(3), &
                   iZ_B1(4), &
                   1       , &
                   1        ], &
                 [ nDOFZ   , &
                   iZ_E1(1), &
                   iZ_E1(2), &
                   iZ_E1(3), &
                   iZ_E1(4), &
                   nCR     , &
                   nSpecies ], &
                 CR )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCF ], &
                 CF )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
                 G )



      END DO

      CALL amrex_mfiter_destroy( MFI )

    END DO


  END SUBROUTINE ComputeFromConserved_TwoMoment_MF
 
END MODULE MF_TwoMoment_UtilitiesModule
