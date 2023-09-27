MODULE MF_Euler_BoundaryConditionsModule

  ! --- AMReX Modules ---

  USE amrex_fort_module, ONLY: &
    amrex_spacedim
  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_geometry_module, ONLY: &
    amrex_geometry
  USE amrex_amrcore_module, ONLY: &
    amrex_geom

  ! --- thornado Modules ---

  USE Euler_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_Euler

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP
  USE MF_EdgeMapModule,               ONLY: &
    EdgeMap

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ApplyBoundaryConditions_Euler_MF

CONTAINS


  SUBROUTINE ApplyBoundaryConditions_Euler_MF &
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

  END SUBROUTINE ApplyBoundaryConditions_Euler_MF



END MODULE MF_Euler_BoundaryConditionsModule
