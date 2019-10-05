MODULE UtilitiesModule_NuclearEOS

  USE KindModule, ONLY: &
    DP

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeFVEC
  PUBLIC :: ComputeFJAC

CONTAINS


  SUBROUTINE ComputeFVEC &
               ( D_R, V_R, P_R, E_R, Y_R, D_L, V_L, P_L, E_L, Y_L, V_Sh, FVEC )

    REAL(DP),               INTENT(in)  :: D_R, V_R, P_R, E_R, Y_R
    REAL(DP),               INTENT(in)  :: D_L, V_L, P_L, E_L, Y_L, V_Sh
    REAL(DP), DIMENSION(4), INTENT(out) :: FVEC

    FVEC(1) = D_L * ( V_Sh - V_L ) - V_Sh * D_R
    FVEC(2) = P_L - V_L * V_Sh * D_R - P_R
    FVEC(3) = E_L + 0.5_DP * V_L**2 - P_L * V_L / ( V_Sh * D_R ) - E_R
    FVEC(4) = Y_L - Y_R

  END SUBROUTINE ComputeFVEC


  SUBROUTINE ComputeFJAC &
               ( D_R, V_R, P_R, E_R, Y_R, D_L, V_L, P_L, E_L, Y_L, &
                 dPdD, dPdT, dPdY, dEdD, dEdT, dEdY, V_Sh, FJAC )

    REAL(DP),                 INTENT(in)  :: D_R, V_R, P_R, E_R, Y_R
    REAL(DP),                 INTENT(in)  :: D_L, V_L, P_L, E_L, Y_L, V_Sh
    REAL(DP),                 INTENT(in)  :: dPdD, dPdT, dPdY
    REAL(DP),                 INTENT(in)  :: dEdD, dEdT, dEdY
    REAL(DP), DIMENSION(4,4), INTENT(out) :: FJAC

    FJAC(1,1) = ( V_Sh - V_L )
    FJAC(1,2) = - D_L
    FJAC(1,3) = 0.0_DP
    FJAC(1,4) = 0.0_DP

    FJAC(2,1) = dPdD
    FJAC(2,2) = - V_Sh * D_R
    FJAC(2,3) = dPdT
    FJAC(2,4) = dPdY

    FJAC(3,1) = dEdD - dPdD * V_L / ( V_Sh * D_R )
    FJAC(3,2) = V_L - P_L / ( V_Sh * D_R )
    FJAC(3,3) = dEdT - dPdT * V_L / ( V_Sh * D_R )
    FJAC(3,4) = dEdY - dPdY * V_L / ( V_Sh * D_R )

    FJAC(4,1) = 0.0_DP
    FJAC(4,2) = 0.0_DP
    FJAC(4,3) = 0.0_DP
    FJAC(4,4) = 1.0_DP

  END SUBROUTINE ComputeFJAC


END MODULE UtilitiesModule_NuclearEOS
