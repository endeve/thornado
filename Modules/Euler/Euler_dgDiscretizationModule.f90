MODULE Euler_dgDiscretizationModule

  USE KindModule, ONLY: &
    DP

#ifdef HYDRO_NONRELATIVISTIC

  USE Euler_dgDiscretizationModule_NonRelativistic, ONLY: &
    Euler_ComputeIncrement_DG_Explicit_NonRelativistic

#elif HYDRO_RELATIVISTIC

  USE Euler_dgDiscretizationModule_Relativistic, ONLY: &
    Euler_ComputeIncrement_DG_Explicit_Relativistic

#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Euler_ComputeIncrement_DG_Explicit


CONTAINS


  SUBROUTINE Euler_ComputeIncrement_DG_Explicit &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU, SuppressBC_Option )

    INTEGER,  INTENT(in)           :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)           :: &
      G (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    REAL(DP), INTENT(inout)        :: &
      U (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    REAL(DP), INTENT(out)          :: &
      dU(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)
    LOGICAL,  INTENT(in), OPTIONAL :: &
      SuppressBC_Option

    LOGICAL :: SuppressBC

    SuppressBC = .FALSE.
    IF( PRESENT( SuppressBC_Option ) ) &
      SuppressBC = SuppressBC_Option

#ifdef HYDRO_NONRELATIVISTIC

    CALL Euler_ComputeIncrement_DG_Explicit_NonRelativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU, SuppressBC )

#elif HYDRO_RELATIVISTIC

    CALL Euler_ComputeIncrement_DG_Explicit_Relativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU, SuppressBC )

#endif

  END SUBROUTINE Euler_ComputeIncrement_DG_Explicit


END MODULE Euler_dgDiscretizationModule
