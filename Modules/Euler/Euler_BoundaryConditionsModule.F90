MODULE Euler_BoundaryConditionsModule

  USE KindModule, ONLY: &
    DP

#ifdef HYDRO_NONRELATIVISTIC

  USE Euler_BoundaryConditionsModule_NonRelativistic, ONLY: &
    Euler_ApplyBoundaryConditions_NonRelativistic, &
    iEuler_ApplyBC_Both, &
    iEuler_ApplyBC_Inner, &
    iEuler_ApplyBC_Outer, &
    iEuler_ApplyBC_None

#elif HYDRO_RELATIVISTIC

  USE Euler_BoundaryConditionsModule_Relativistic, ONLY: &
    Euler_ApplyBoundaryConditions_Relativistic, &
    iEuler_ApplyBC_Both, &
    iEuler_ApplyBC_Inner, &
    iEuler_ApplyBC_Outer, &
    iEuler_ApplyBC_None

#endif

  IMPLICIT NONE
  PRIVATE


  PUBLIC :: Euler_ApplyBoundaryConditions


CONTAINS


  SUBROUTINE Euler_ApplyBoundaryConditions &
    ( iX_B0, iX_E0, iX_B1, iX_E1, U, iApplyBC_Option )

    INTEGER,  INTENT(in)            :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout)        :: &
      U (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    INTEGER,  INTENT(in), OPTIONAL :: &
      iApplyBC_Option(3)

    INTEGER :: iApplyBC(3)

    iApplyBC = iEuler_ApplyBC_Both
    IF( PRESENT( iApplyBC_Option ) ) &
      iApplyBC = iApplyBC_Option

#ifdef HYDRO_NONRELATIVISTIC

    CALL Euler_ApplyBoundaryConditions_NonRelativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, U, iApplyBC )

#elif HYDRO_RELATIVISTIC

    CALL Euler_ApplyBoundaryConditions_Relativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, U, iApplyBC )

#endif

  END SUBROUTINE Euler_ApplyBoundaryConditions


END MODULE Euler_BoundaryConditionsModule
