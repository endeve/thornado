MODULE MF_TwoMoment_OpacityModule

  ! --- AMReX Modules ---

  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_geometry_module, ONLY: &
    amrex_geometry
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor
  USE amrex_parmparse_module, ONLY: &
    amrex_parmparse, &
    amrex_parmparse_build, &
    amrex_parmparse_destroy

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    swX, &
    nDOFX, &
    nDOFZ, &
    iE_B0, &
    iE_E0, &
    iE_B1, &
    iE_E1
  USE MeshModule, ONLY: &
    MeshX
  USE GeometryFieldsModule, ONLY: &
    nGF
  USE GeometryFieldsModuleE, ONLY: &
    uGE
  USE RadiationFieldsModule, ONLY: &
    nCR, &
    nSpecies
  USE FluidFieldsModule, ONLY: &
    nCF
  USE TwoMoment_SlopeLimiterModule, ONLY: &
    InitializeSlopeLimiter_TwoMoment, &
    FinalizeSlopeLimiter_TwoMoment, &
    ApplySlopeLimiter_TwoMoment
  USE TwoMoment_TroubledCellIndicatorModule, ONLY: &
    InitializeTroubledCellIndicator_TwoMoment, &
    FinalizeTroubledCellIndicator_TwoMoment

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP
  USE MF_UtilitiesModule, ONLY: &
    amrex2thornado_X, &
    amrex2thornado_Z, &
    thornado2amrex_Z, &
    AllocateArray_X, &
    DeallocateArray_X, &
    AllocateArray_Z, &
    DeallocateArray_Z
  USE InputParsingModule, ONLY: &
    nLevels, &
    UseTiling, &
    nE
  USE MF_TwoMoment_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_TwoMoment_MF
  USE MF_EdgeMapModule, ONLY: &
    EdgeMap,          &
    ConstructEdgeMap
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: SetOpacities_TwoMoment_MF
  !PUBLIC :: CreateOpacities_TwoMoment_MF
  !PUBLIC :: DestroyOpacities_TwoMoment_MF


CONTAINS


  SUBROUTINE SetOpacities_TwoMoment_MF &
    ( GEOM, MF_uGF, Verbose_Option )

    TYPE(amrex_geometry), INTENT(in)    :: GEOM  (0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:nLevels-1)

    LOGICAL             , INTENT(in), OPTIONAL :: Verbose_Option

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)

    REAL(DP), ALLOCATABLE :: G(:,:,:,:,:)

    INTEGER :: iLevel
    INTEGER :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    INTEGER :: iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4), iLo_MF(4)

    LOGICAL :: Verbose

    Verbose = .FALSE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    DO iLevel = 0, nLevels-1

      ! --- Apply boundary conditions to interior domains ---

      CALL MF_uGF(iLevel) % Fill_Boundary( GEOM(iLevel) )

      CALL CreateMesh_MF( iLevel, MeshX )

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )

        iLo_MF = LBOUND( uGF )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = BX % lo - swX
        iX_E1 = BX % hi + swX

        iZ_B0(1) = iE_B0
        iZ_B1(1) = iE_B1
        iZ_E0(1) = iE_E0
        iZ_E1(1) = iE_E1

        iZ_B0(2:4) = iX_B0
        iZ_B1(2:4) = iX_B1
        iZ_E0(2:4) = iX_E0
        iZ_E1(2:4) = iX_E1

        CALL AllocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
                 G )

        CALL amrex2thornado_X &
               ( nGF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uGF, G )

        CALL SetOpacities &
               ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, D0, Chi, Sigma, Verbose_Option )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
                 G )

      END DO ! WHILE( MFI % next() )

      CALL amrex_mfiter_destroy( MFI )

      CALL DestroyMesh_MF( MeshX )

    END DO ! iLevel = 0, nLevels-1

  END SUBROUTINE SetOpacities_TwoMoment_MF

