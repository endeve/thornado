MODULE MF_MHD_BoundaryConditionsModule

  ! --- AMReX Modules ---

  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    swX, &
    nDOFX
  USE MeshModule, ONLY: &
    MeshX
  USE MagnetofluidFieldsModule, ONLY: &
    nCM, &
    nDM
  USE MHD_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_MHD

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP
  USE MF_EdgeMapModule_MHD, ONLY: &
    EdgeMap, &
    ConstructEdgeMap
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF
  USE MF_UtilitiesModule_MHD, ONLY: &
    AllocateArray_X, &
    DeallocateArray_X, &
    amrex2thornado_X, &
    thornado2amrex_X
  USE InputParsingModule, ONLY: &
    UseTiling, &
    nLevels

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ApplyBoundaryConditions_MHD_MF

  INTERFACE ApplyBoundaryConditions_MHD_MF
    MODULE PROCEDURE ApplyBoundaryConditions_MHD_MF_MultiLevel
    MODULE PROCEDURE ApplyBoundaryConditions_MHD_MF_SingleLevel
    MODULE PROCEDURE ApplyBoundaryConditions_MHD_MF_SingleLevel_Box
  END INTERFACE ApplyBoundaryConditions_MHD_MF

CONTAINS


  SUBROUTINE ApplyBoundaryConditions_MHD_MF_MultiLevel( t, MF_uCM, MF_uDM )

    REAL(DP),             INTENT(in   ) :: t(0:)
    TYPE(amrex_multifab), INTENT(in   ) :: MF_uDM(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCM(0:)

    INTEGER :: iLevel

    DO iLevel = 0, nLevels-1

      CALL ApplyBoundaryConditions_MHD_MF( t(iLevel), iLevel, MF_uCM(iLevel), MF_uDM(iLevel) )

    END DO

  END SUBROUTINE ApplyBoundaryConditions_MHD_MF_MultiLevel


  SUBROUTINE ApplyBoundaryConditions_MHD_MF_SingleLevel( t, iLevel, MF_uCM, MF_uDM )

    REAL(DP),             INTENT(in   ) :: t
    INTEGER             , INTENT(in)    :: iLevel
    TYPE(amrex_multifab), INTENT(in)    :: MF_uDM
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCM

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX
    TYPE(EdgeMap)      :: Edge_Map

    REAL(DP), CONTIGUOUS, POINTER :: uCM(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uDM(:,:,:,:)
    REAL(DP), ALLOCATABLE         :: U  (:,:,:,:,:)
    REAL(DP), ALLOCATABLE         :: D  (:,:,:,:,:)


    INTEGER :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), iLo_MF(4)

    CALL CreateMesh_MF( iLevel, MeshX )

    CALL amrex_mfiter_build( MFI, MF_uCM, tiling = UseTiling )

    CALL amrex_mfiter_build( MFI, MF_uDM, tiling = UseTiling )

    DO WHILE( MFI % next() )

      uCM => MF_uCM % DataPtr( MFI )
      uDM => MF_uDM % DataPtr( MFI )

      iLo_MF = LBOUND( uCM )

      BX = MFI % tilebox()

      iX_B0 = BX % lo
      iX_E0 = BX % hi
      iX_B1 = BX % lo - swX
      iX_E1 = BX % hi + swX

      CALL AllocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCM ], &
               U )

      CALL AllocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nDM ], &
               D )

      CALL amrex2thornado_X( nCM, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uCM, U )

      CALL amrex2thornado_X( nDM, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uDM, D )

      ! --- Apply boundary conditions to physical boundaries ---

      CALL ConstructEdgeMap( iLevel, BX, Edge_Map )

      CALL ApplyBoundaryConditions_MHD_MF &
             ( t, iX_B0, iX_E0, iX_B1, iX_E1, U, D, Edge_Map )

      CALL thornado2amrex_X( nCM, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uCM, U )

      CALL DeallocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCM ], &
               U )

      CALL DeallocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nDM ], &
               D )

    END DO

    CALL amrex_mfiter_destroy( MFI )

  END SUBROUTINE ApplyBoundaryConditions_MHD_MF_SingleLevel


  SUBROUTINE ApplyBoundaryConditions_MHD_MF_SingleLevel_Box &
    ( t, iX_B0, iX_E0, iX_B1, iX_E1, U, D, Edge_Map )

    REAL(DP),             INTENT(in   ) :: t
    INTEGER,       INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP),      INTENT(in)    :: &
      D(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP),      INTENT(inout) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    TYPE(EdgeMap), INTENT(in)    :: &
      Edge_Map

    INTEGER :: iApplyBC(3)

    CALL Edge_Map % GetBC( iApplyBC )

    CALL ApplyBoundaryConditions_MHD &
           ( t, iX_B0, iX_E0, iX_B1, iX_E1, U, D, iApplyBC )

  END SUBROUTINE ApplyBoundaryConditions_MHD_MF_SingleLevel_Box


END MODULE MF_MHD_BoundaryConditionsModule
