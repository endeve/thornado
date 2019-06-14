MODULE Euler_PositivityLimiterModule

  USE KindModule, ONLY: &
    DP







  USE Euler_PositivityLimiterModule_Relativistic




  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Euler_InitializePositivityLimiter
  PUBLIC :: Euler_FinalizePositivityLimiter
  PUBLIC :: Euler_ApplyPositivityLimiter


CONTAINS


  SUBROUTINE Euler_InitializePositivityLimiter &
    ( Min_1_Option, Min_2_Option, UsePositivityLimiter_Option, Verbose_Option )

    REAL(DP), INTENT(in), OPTIONAL :: Min_1_Option
    REAL(DP), INTENT(in), OPTIONAL :: Min_2_Option
    LOGICAL,  INTENT(in), OPTIONAL :: UsePositivityLimiter_Option
    LOGICAL,  INTENT(in), OPTIONAL :: Verbose_Option









    CALL Euler_InitializePositivityLimiter_Relativistic &
           ( Min_1_Option, Min_2_Option, &
             UsePositivityLimiter_Option, Verbose_Option )



  END SUBROUTINE Euler_InitializePositivityLimiter


  SUBROUTINE Euler_FinalizePositivityLimiter







    CALL Euler_FinalizePositivityLimiter_Relativistic



  END SUBROUTINE Euler_FinalizePositivityLimiter


  SUBROUTINE Euler_ApplyPositivityLimiter &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)








    CALL Euler_ApplyPositivityLimiter_Relativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U )



  END SUBROUTINE Euler_ApplyPositivityLimiter


END MODULE Euler_PositivityLimiterModule

