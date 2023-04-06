MODULE MF_Euler_BoundaryConditionsModule

  ! --- AMReX Modules ---

  USE amrex_fort_module,              ONLY: &
    amrex_spacedim
  USE amrex_box_module,               ONLY: &
    amrex_box
  USE amrex_geometry_module,          ONLY: &
    amrex_geometry

  ! --- thornado Modules ---

  USE Euler_BoundaryConditionsModule, ONLY: &
    iApplyBC_Euler_Both,  &
    iApplyBC_Euler_Inner, &
    iApplyBC_Euler_Outer, &
    iApplyBC_Euler_None,  &
    ApplyBoundaryConditions_Euler

  ! --- Local Modules ---

  USE MF_KindModule,                  ONLY: &
    DP
  USE InputParsingModule,             ONLY: &
    DEBUG
  USE TimersModule_AMReX_Euler,       ONLY: &
    TimersStart_AMReX_Euler,            &
    TimersStop_AMReX_Euler,             &
    Timer_AMReX_Euler_ConstructEdgeMap, &
    Timer_AMReX_Euler_GetBC
  USE MF_EdgeMapModule,               ONLY: &
    EdgeMap
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MF_ApplyBoundaryConditions_Euler
  PUBLIC :: ConstructEdgeMap

  TYPE, PUBLIC :: EdgeMap
    LOGICAL :: IsLowerBoundary(3)
    LOGICAL :: IsUpperBoundary(3)
  CONTAINS
    PROCEDURE :: Euler_GetBC => EdgeMap_Euler_GetBC
  END TYPE EdgeMap


CONTAINS


  SUBROUTINE MF_ApplyBoundaryConditions_Euler &
    ( iX_B0, iX_E0, iX_B1, iX_E1, U, Edge_Map )

    INTEGER,       INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP),      INTENT(inout) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    TYPE(EdgeMap), INTENT(in)    :: &
      Edge_Map

    INTEGER :: iApplyBC(3)

    CALL Edge_Map % Euler_GetBC( iApplyBC )

    IF( DEBUG ) WRITE(*,'(A)') '      CALL ApplyBoundaryConditions_Euler'

    CALL ApplyBoundaryConditions_Euler &
           ( iX_B0, iX_E0, iX_B1, iX_E1, U, iApplyBC )

  END SUBROUTINE MF_ApplyBoundaryConditions_Euler


  SUBROUTINE ConstructEdgeMap( GEOM, BX, Edge_Map )

    TYPE(amrex_geometry), INTENT(in)    :: GEOM
    TYPE(amrex_box),      INTENT(in)    :: BX
    TYPE(EdgeMap),        INTENT(inout) :: Edge_Map

    INTEGER :: iDim

    CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_ConstructEdgeMap )

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

    CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_ConstructEdgeMap )

  END SUBROUTINE ConstructEdgeMap


  SUBROUTINE EdgeMap_Euler_GetBC( this, iApplyBC )

    CLASS(EdgeMap), INTENT(in)  :: this
    INTEGER,        INTENT(out) :: iApplyBC(3)

    INTEGER :: iDim

    CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_GetBC )

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

    CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_GetBC )

  END SUBROUTINE EdgeMap_Euler_GetBC


END MODULE MF_Euler_BoundaryConditionsModule
