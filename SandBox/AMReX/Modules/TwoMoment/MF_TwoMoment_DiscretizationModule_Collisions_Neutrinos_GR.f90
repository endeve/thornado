MODULE MF_TwoMoment_DiscretizationModule_Collisions_Neutrinos_GR

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

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    swX, &
    swE, &
    nDOFX, &
    nDOFE, &
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
    nGE, &
    uGE
  USE RadiationFieldsModule, ONLY: &
    nCR
  USE FluidFieldsModule, ONLY: &
    nCF
  USE TwoMoment_DiscretizationModule_Collisions_Neutrinos, ONLY: &
    ComputeIncrement_TwoMoment_Implicit


  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero
  USE MF_UtilitiesModule, ONLY: &
    amrex2thornado_X, &
    amrex2thornado_Z, &
    thornado2amrex_Z, &
    thornado2amrex_X, &
    AllocateArray_X, &
    DeallocateArray_X, &
    AllocateArray_Z, &
    DeallocateArray_Z
  USE InputParsingModule,ONLY: &
    nLevels, &
    nSpecies, &
    UseTiling, &
    nE
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF
  USE MF_TwoMoment_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_TwoMoment_MF
  USE MF_Euler_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_Euler_MF
  USE MF_EdgeMapModule, ONLY: &
    EdgeMap,          &
    ConstructEdgeMap

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeIncrement_TwoMoment_Implicit_Neutrinos_MF


CONTAINS

  SUBROUTINE ComputeIncrement_TwoMoment_Implicit_Neutrinos_MF &
    ( GEOM, MF_uGF, MF_uCF, MF_duCF, MF_uCR, MF_duCR, dt, Verbose_Option )

    TYPE(amrex_geometry), INTENT(in)    :: GEOM   (0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF (0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in)    :: MF_uCF (0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_duCF(0:nLevels-1)
    REAL(DP)            , INTENT(in)    :: dt
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR (0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_duCR(0:nLevels-1)
    LOGICAL             , INTENT(in), OPTIONAL :: Verbose_Option

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    REAL(DP), CONTIGUOUS, POINTER :: uGF (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: duCF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCR (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: duCR(:,:,:,:)

    REAL(DP), ALLOCATABLE :: G (:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: C (:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: dC(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: U (:,:,:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: dU(:,:,:,:,:,:,:)

    INTEGER :: iLevel
    INTEGER :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), iLo_MF(4)
    INTEGER :: iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    LOGICAL :: Verbose

    TYPE(EdgeMap) :: Edge_Map

    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    DO iLevel = 0, nLevels-1

      ! --- Apply boundary conditions to interior domains ---

      CALL MF_uCF(iLevel) % Fill_Boundary( GEOM(iLevel) )

      CALL MF_uCR(iLevel) % Fill_Boundary( GEOM(iLevel) )

      CALL MF_duCF(iLevel) % setval( 0.0_DP )

      CALL MF_duCR(iLevel) % setval( 0.0_DP )

      CALL CreateMesh_MF( iLevel, MeshX )

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF  => MF_uGF (iLevel) % DataPtr( MFI )
        uCF  => MF_uCF (iLevel) % DataPtr( MFI )
        duCF => MF_duCF(iLevel) % DataPtr( MFI )
        uCR  => MF_uCR (iLevel) % DataPtr( MFI )
        duCR => MF_duCR(iLevel) % DataPtr( MFI )

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

        CALL AllocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCF ], &
                 C )

        CALL AllocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCF ], &
                 dC )

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
                 U )

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
                 dU )

        CALL amrex2thornado_X &
               ( nGF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uGF, G )

        CALL amrex2thornado_X &
               ( nCF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uCF, C )

        CALL amrex2thornado_X &
               ( nCF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, duCF, dC )

        CALL amrex2thornado_Z &
               ( nCR, nSpecies, nE, iE_B0, iE_E0, &
                 iZ_B1, iZ_E1, iLo_MF, iZ_B1, iZ_E1, uCR, U )

        CALL ConstructEdgeMap( iLevel, BX, Edge_Map )

        CALL ApplyBoundaryConditions_Euler_MF &
               ( iX_B0, iX_E0, iX_B1, iX_E1, C, Edge_Map )

        CALL ApplyBoundaryConditions_TwoMoment_MF &
               ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U, Edge_Map )

        CALL ComputeIncrement_TwoMoment_Implicit &
               ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, dt, uGE, G, C, dC, U, dU )
        CALL thornado2amrex_X &
               ( nCF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, duCF, dC )

        CALL thornado2amrex_Z &
               ( nCR, nSpecies, nE, iE_B0, iE_E0, &
                 iZ_B1, iZ_E1, iLo_MF, iZ_B0, iZ_E0, duCR, dU )

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
                 dU )

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
                 U )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCF ], &
                 dC )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCF ], &
                 C )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
                 G )

      END DO ! DO WHILE( MFI % next() )

      CALL amrex_mfiter_destroy( MFI )

      CALL DestroyMesh_MF( MeshX )

    END DO ! iLevel = 0, nLevels-1

  END SUBROUTINE ComputeIncrement_TwoMoment_Implicit_Neutrinos_MF


END MODULE  MF_TwoMoment_DiscretizationModule_Collisions_Neutrinos_GR
