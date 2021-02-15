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

  USE MyAmrModule,              ONLY: &
    DEBUG

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MF_ApplyBoundaryConditions_TwoMoment
  PUBLIC :: ConstructEdgeMap

  TYPE, PUBLIC :: EdgeMap
    LOGICAL :: IsLowerBoundary(3)
    LOGICAL :: IsUpperBoundary(3)
  CONTAINS
    PROCEDURE :: TwoMoment_GetBC => EdgeMap_TwoMoment_GetBC
  END TYPE EdgeMap

  ! --- Hack to get iApplyBC_Euler_XXX. DO NOT CHANGE THESE VALUES ---
  INTEGER, PARAMETER :: iApplyBC_Euler_Both  = 0
  INTEGER, PARAMETER :: iApplyBC_Euler_Inner = 1
  INTEGER, PARAMETER :: iApplyBC_Euler_Outer = 2
  INTEGER, PARAMETER :: iApplyBC_Euler_None  = 3


CONTAINS


  SUBROUTINE MF_ApplyBoundaryConditions_TwoMoment &
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

    CALL Edge_Map % TwoMoment_GetBC( iApplyBC )

    IF( DEBUG ) WRITE(*,'(A)') '      CALL ApplyBoundaryConditions_TwoMoment'

    CALL ApplyBoundaryConditions_TwoMoment &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U )

  END SUBROUTINE MF_ApplyBoundaryConditions_TwoMoment


  SUBROUTINE ConstructEdgeMap( GEOM, BX, Edge_Map )

    TYPE(amrex_geometry), INTENT(in   ) :: GEOM
    TYPE(amrex_box),      INTENT(in   ) :: BX
    TYPE(EdgeMap),        INTENT(inout) :: Edge_Map

    INTEGER :: iDim


    Edge_Map % IsLowerBoundary = .FALSE.
    Edge_Map % IsUpperBoundary = .FALSE.

    DO iDim = 1, 3

      IF( iDim .LE. amrex_spacedim )THEN

        IF( BX % lo( iDim ) .LE. GEOM % DOMAIN % lo( iDim ) ) &
          Edge_Map % IsLowerBoundary( iDim ) = .TRUE.

        IF( BX % hi( iDim ) .GE. GEOM % DOMAIN % hi( iDim ) ) &
          Edge_Map % IsUpperBoundary( iDim ) = .TRUE.

      END IF

    END DO


  END SUBROUTINE ConstructEdgeMap




  SUBROUTINE EdgeMap_TwoMoment_GetBC( this, iApplyBC )

    CLASS(EdgeMap), INTENT(in)  :: this
    INTEGER,        INTENT(out) :: iApplyBC(3)

    INTEGER :: iDim


    DO iDim = 1, 3

      IF     ( this % IsLowerBoundary( iDim ) .AND. &
               this % IsUpperBoundary( iDim ) )THEN

        iApplyBC(iDim) = iApplyBC_Euler_Both

      ELSE IF( this % IsLowerBoundary( iDim ) )THEN

        iApplyBC(iDim) = iApplyBC_Euler_Inner

      ELSE IF( this % IsUpperBoundary( iDim ) )THEN

        iApplyBC(iDim) = iApplyBC_Euler_Outer

      ELSE

        iApplyBC(iDim) = iApplyBC_Euler_None

      END IF

    END DO


  END SUBROUTINE EdgeMap_TwoMoment_GetBC




END MODULE MF_TwoMoment_BoundaryConditionsModule
