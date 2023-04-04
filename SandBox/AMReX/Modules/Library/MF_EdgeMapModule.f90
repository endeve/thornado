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

    INTEGER,         INTENT(in)    :: iLevel
    TYPE(amrex_box), INTENT(in)    :: BX
    TYPE(EdgeMap),   INTENT(inout) :: Edge_Map

    INTEGER :: iDim

    Edge_Map % IsLowerBoundary = .FALSE.
    Edge_Map % IsUpperBoundary = .FALSE.

    DO iDim = 1, 3

      IF( iDim .LE. amrex_spacedim )THEN

        IF( BX % lo( iDim ) .LE. amrex_geom(iLevel) % DOMAIN % lo( iDim ) ) &
          Edge_Map % IsLowerBoundary( iDim ) = .TRUE.

        IF( BX % hi( iDim ) .GE. amrex_geom(iLevel) % DOMAIN % hi( iDim ) ) &
          Edge_Map % IsUpperBoundary( iDim ) = .TRUE.

      END IF

    END DO

  END SUBROUTINE ConstructEdgeMap


  SUBROUTINE EdgeMap_GetBC( this, iApplyBC )

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

  END SUBROUTINE EdgeMap_GetBC


END MODULE MF_EdgeMapModule
