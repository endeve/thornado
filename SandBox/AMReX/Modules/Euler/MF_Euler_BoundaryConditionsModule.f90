MODULE MF_Euler_BoundaryConditionsModule

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
  USE FluidFieldsModule, ONLY: &
    nCF
  USE Euler_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_Euler

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP
  USE MF_EdgeMapModule_Euler, ONLY: &
    EdgeMap, &
    ConstructEdgeMap
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF
  USE MF_UtilitiesModule, ONLY: &
    AllocateArray_X, &
    DeallocateArray_X, &
    amrex2thornado_X, &
    thornado2amrex_X
  USE InputParsingModule, ONLY: &
    UseTiling, &
    nLevels

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ApplyBoundaryConditions_Euler_MF

  INTERFACE ApplyBoundaryConditions_Euler_MF
    MODULE PROCEDURE ApplyBoundaryConditions_Euler_MF_MultiLevel
    MODULE PROCEDURE ApplyBoundaryConditions_Euler_MF_SingleLevel
    MODULE PROCEDURE ApplyBoundaryConditions_Euler_MF_SingleLevel_Box
  END INTERFACE ApplyBoundaryConditions_Euler_MF

CONTAINS


  SUBROUTINE ApplyBoundaryConditions_Euler_MF_MultiLevel( MF_uCF )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:)

    INTEGER :: iLevel

    DO iLevel = 0, nLevels-1

      CALL ApplyBoundaryConditions_Euler_MF( iLevel, MF_uCF(iLevel) )

    END DO

  END SUBROUTINE ApplyBoundaryConditions_Euler_MF_MultiLevel


  SUBROUTINE ApplyBoundaryConditions_Euler_MF_SingleLevel( iLevel, MF_uCF )

    INTEGER             , INTENT(in)    :: iLevel
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX
    TYPE(EdgeMap)      :: Edge_Map

    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    REAL(DP), ALLOCATABLE         :: U  (:,:,:,:,:)

    INTEGER :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), iLo_MF(4)

    CALL CreateMesh_MF( iLevel, MeshX )

    CALL amrex_mfiter_build( MFI, MF_uCF, tiling = UseTiling )

    DO WHILE( MFI % next() )

      uCF => MF_uCF % DataPtr( MFI )

      iLo_MF = LBOUND( uCF )

      BX = MFI % tilebox()

      iX_B0 = BX % lo
      iX_E0 = BX % hi
      iX_B1 = BX % lo - swX
      iX_E1 = BX % hi + swX

      CALL AllocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCF ], &
               U )

      CALL amrex2thornado_X( nCF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uCF, U )

      ! --- Apply boundary conditions to physical boundaries ---

      CALL ConstructEdgeMap( iLevel, BX, Edge_Map )

      CALL ApplyBoundaryConditions_Euler_MF &
             ( iX_B0, iX_E0, iX_B1, iX_E1, U, Edge_Map )

      CALL thornado2amrex_X( nCF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uCF, U )

      CALL DeallocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCF ], &
               U )

    END DO

    CALL amrex_mfiter_destroy( MFI )

  END SUBROUTINE ApplyBoundaryConditions_Euler_MF_SingleLevel


  SUBROUTINE ApplyBoundaryConditions_Euler_MF_SingleLevel_Box &
    ( iX_B0, iX_E0, iX_B1, iX_E1, U, Edge_Map )

    INTEGER,       INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP),      INTENT(inout) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    TYPE(EdgeMap), INTENT(in)    :: &
      Edge_Map

    INTEGER :: iApplyBC(3)

    CALL Edge_Map % GetBC( iApplyBC )

    CALL ApplyBoundaryConditions_Euler &
           ( iX_B0, iX_E0, iX_B1, iX_E1, U, iApplyBC )

  END SUBROUTINE ApplyBoundaryConditions_Euler_MF_SingleLevel_Box


END MODULE MF_Euler_BoundaryConditionsModule
