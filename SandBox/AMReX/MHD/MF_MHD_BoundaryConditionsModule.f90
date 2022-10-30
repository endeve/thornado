MODULE MF_MHD_BoundaryConditionsModule

  ! --- AMReX Modules ---

  USE amrex_fort_module,              ONLY: &
    amrex_spacedim
  USE amrex_box_module,               ONLY: &
    amrex_box
  USE amrex_geometry_module,          ONLY: &
    amrex_geometry

  ! --- thornado Modules ---

  USE MHD_BoundaryConditionsModule, ONLY: &
    iApplyBC_MHD_Both,  &
    iApplyBC_MHD_Inner, &
    iApplyBC_MHD_Outer, &
    iApplyBC_MHD_None,  &
    ApplyBoundaryConditions_MHD

  ! --- Local Modules ---

  USE MF_KindModule,                  ONLY: &
    DP
  USE InputParsingModule,             ONLY: &
    DEBUG

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MF_ApplyBoundaryConditions_MHD
  PUBLIC :: ConstructEdgeMap

  TYPE, PUBLIC :: EdgeMap
    LOGICAL :: IsLowerBoundary(3)
    LOGICAL :: IsUpperBoundary(3)
  CONTAINS
    PROCEDURE :: MHD_GetBC => EdgeMap_MHD_GetBC
  END TYPE EdgeMap


CONTAINS


  SUBROUTINE MF_ApplyBoundaryConditions_MHD &
    ( iX_B0, iX_E0, iX_B1, iX_E1, U, Edge_Map )

    INTEGER,       INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP),      INTENT(inout) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    TYPE(EdgeMap), INTENT(in)    :: &
      Edge_Map

    INTEGER :: iApplyBC(3)

    CALL Edge_Map % MHD_GetBC( iApplyBC )

    IF( DEBUG ) WRITE(*,'(A)') '      CALL ApplyBoundaryConditions_MHD'

    CALL ApplyBoundaryConditions_MHD &
           ( iX_B0, iX_E0, iX_B1, iX_E1, U, iApplyBC )

  END SUBROUTINE MF_ApplyBoundaryConditions_MHD


  SUBROUTINE ConstructEdgeMap( GEOM, BX, Edge_Map )

    TYPE(amrex_geometry), INTENT(in)    :: GEOM
    TYPE(amrex_box),      INTENT(in)    :: BX
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


  SUBROUTINE EdgeMap_MHD_GetBC( this, iApplyBC )

    CLASS(EdgeMap), INTENT(in)  :: this
    INTEGER,        INTENT(out) :: iApplyBC(3)

    INTEGER :: iDim

    DO iDim = 1, 3

      IF     ( this % IsLowerBoundary( iDim ) .AND. &
               this % IsUpperBoundary( iDim ) )THEN

        iApplyBC(iDim) = iApplyBC_MHD_Both

      ELSE IF( this % IsLowerBoundary( iDim ) )THEN

        iApplyBC(iDim) = iApplyBC_MHD_Inner

      ELSE IF( this % IsUpperBoundary( iDim ) )THEN

        iApplyBC(iDim) = iApplyBC_MHD_Outer

      ELSE

        iApplyBC(iDim) = iApplyBC_MHD_None

      END IF

    END DO

  END SUBROUTINE EdgeMap_MHD_GetBC


END MODULE MF_MHD_BoundaryConditionsModule
