 MODULE MF_TwoMoment_BoundaryConditionsModule

  ! --- AMReX Modules ---

  USE amrex_fort_module,     ONLY: &
    AR => amrex_real, &
    amrex_spacedim
  USE amrex_box_module,      ONLY: &
    amrex_box
  USE amrex_geometry_module, ONLY: &
    amrex_geometry

  ! --- thornado Modules ---

  USE TwoMoment_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_TwoMoment
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    nCR
  USE ProgramHeaderModule, ONLY: &
    nDOF, nDOFE, nDOFX
  USE GeometryFieldsModuleE, ONLY: &
    nGE
  USE GeometryFieldsModule, ONLY: &
    nGF
  ! --- Local Modules ---

  USE InputParsingModule,              ONLY: &
    DEBUG
  USE MF_EdgeMapModule,                ONLY: &
    EdgeMap

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ApplyBoundaryConditions_TwoMoment_MF


CONTAINS


  SUBROUTINE ApplyBoundaryConditions_TwoMoment_MF &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U, Edge_Map )

    INTEGER,       INTENT(in   ) :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(AR),      INTENT(inout) :: &
      U(1:nDOF, &
        iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
        iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4), &
        1:nCR,1:nSpecies)
    TYPE(EdgeMap), INTENT(in   ) :: &
      Edge_Map

    INTEGER :: iApplyBC(3)

    CALL Edge_Map % GetBC( iApplyBC )

    IF( DEBUG ) WRITE(*,'(A)') '      CALL ApplyBoundaryConditions_TwoMoment'

    CALL ApplyBoundaryConditions_TwoMoment &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U, &
             iApplyBC_Option = iApplyBC )

  END SUBROUTINE ApplyBoundaryConditions_TwoMoment_MF




END MODULE MF_TwoMoment_BoundaryConditionsModule
