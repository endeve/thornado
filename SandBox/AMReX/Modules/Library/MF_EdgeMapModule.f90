MODULE MF_EdgeMapModule

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
    iApplyBC_Euler_Both,  &
    iApplyBC_Euler_Inner, &
    iApplyBC_Euler_Outer, &
    iApplyBC_Euler_None

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ConstructEdgeMap

  TYPE, PUBLIC :: EdgeMap
    LOGICAL :: IsLowerBoundary(3)
    LOGICAL :: IsUpperBoundary(3)
  CONTAINS
    PROCEDURE :: GetBC => EdgeMap_GetBC
  END TYPE EdgeMap

CONTAINS


  SUBROUTINE ConstructEdgeMap( iLevel, BX, Edge_Map )

    INTEGER        , INTENT(in)    :: iLevel
    TYPE(amrex_box), INTENT(in)    :: BX
    TYPE(EdgeMap)  , INTENT(inout) :: Edge_Map

    INTEGER :: iDimX

    Edge_Map % IsLowerBoundary = .FALSE.
    Edge_Map % IsUpperBoundary = .FALSE.

    DO iDimX = 1, 3

      IF( iDimX .LE. amrex_spacedim )THEN

        IF( BX % lo( iDimX ) .LE. amrex_geom(iLevel) % DOMAIN % lo( iDimX ) ) &
          Edge_Map % IsLowerBoundary( iDimX ) = .TRUE.

        IF( BX % hi( iDimX ) .GE. amrex_geom(iLevel) % DOMAIN % hi( iDimX ) ) &
          Edge_Map % IsUpperBoundary( iDimX ) = .TRUE.

      END IF

    END DO

  END SUBROUTINE ConstructEdgeMap


  SUBROUTINE EdgeMap_GetBC( this, iApplyBC )

    CLASS(EdgeMap), INTENT(in)  :: this
    INTEGER       , INTENT(out) :: iApplyBC(3)

    INTEGER :: iDimX

    DO iDimX = 1, 3

      IF     ( this % IsLowerBoundary( iDimX ) .AND. &
               this % IsUpperBoundary( iDimX ) )THEN

        iApplyBC(iDimX) = iApplyBC_Euler_Both

      ELSE IF( this % IsLowerBoundary( iDimX ) )THEN

        iApplyBC(iDimX) = iApplyBC_Euler_Inner

      ELSE IF( this % IsUpperBoundary( iDimX ) )THEN

        iApplyBC(iDimX) = iApplyBC_Euler_Outer

      ELSE

        iApplyBC(iDimX) = iApplyBC_Euler_None

      END IF

    END DO

  END SUBROUTINE EdgeMap_GetBC


END MODULE MF_EdgeMapModule
