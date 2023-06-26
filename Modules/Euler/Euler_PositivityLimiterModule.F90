MODULE Euler_PositivityLimiterModule

  USE KindModule, ONLY: &
    DP, &
    Zero, &
    One

#if   defined( MICROPHYSICS_WEAKLIB ) && defined( HYDRO_RELATIVISTIC    )

  USE Euler_PositivityLimiterModule_Relativistic_TABLE

#elif defined( MICROPHYSICS_WEAKLIB ) && defined( HYDRO_NONRELATIVISTIC )

  USE Euler_PositivityLimiterModule_NonRelativistic_TABLE

#elif defined( HYDRO_RELATIVISTIC )

  USE Euler_PositivityLimiterModule_Relativistic_IDEAL

#else

  USE Euler_PositivityLimiterModule_NonRelativistic_IDEAL

#endif


  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializePositivityLimiter_Euler
  PUBLIC :: FinalizePositivityLimiter_Euler
  PUBLIC :: ApplyPositivityLimiter_Euler


CONTAINS


  SUBROUTINE InitializePositivityLimiter_Euler &
    ( UsePositivityLimiter_Option, Verbose_Option, &
      Min_1_Option, Min_2_Option, Min_3_Option, &
      Max_1_Option, Max_2_Option, Max_3_Option, &
      D_Min_Euler_PL_Option, IntE_Min_Euler_PL_Option )

    LOGICAL , INTENT(in), OPTIONAL :: UsePositivityLimiter_Option
    LOGICAL , INTENT(in), OPTIONAL :: Verbose_Option
    REAL(DP), INTENT(in), OPTIONAL :: Min_1_Option, Max_1_Option
    REAL(DP), INTENT(in), OPTIONAL :: Min_2_Option, Max_2_Option
    REAL(DP), INTENT(in), OPTIONAL :: Min_3_Option, Max_3_Option
    REAL(DP), INTENT(in), OPTIONAL :: D_Min_Euler_PL_Option
    REAL(DP), INTENT(in), OPTIONAL :: IntE_Min_Euler_PL_Option

    LOGICAL  :: UsePositivityLimiter
    LOGICAL  :: Verbose
    REAL(DP) :: Min_1, Max_1
    REAL(DP) :: Min_2, Max_2
    REAL(DP) :: Min_3, Max_3
    REAL(DP) :: D_Min_Euler_PL
    REAL(DP) :: IntE_Min_Euler_PL

    UsePositivityLimiter = .TRUE.
    IF( PRESENT( UsePositivityLimiter_Option ) ) &
      UsePositivityLimiter = UsePositivityLimiter_Option

    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

#if   defined( MICROPHYSICS_WEAKLIB ) && defined( HYDRO_RELATIVISTIC    )

    Min_1 = - HUGE( One )
    IF( PRESENT( Min_1_Option ) ) &
      Min_1 = Min_1_Option

    Min_2 = - HUGE( One )
    IF( PRESENT( Min_2_Option ) ) &
      Min_2 = Min_2_Option

    Min_3 = - HUGE( One )
    IF( PRESENT( Min_3_Option ) ) &
      Min_3 = Min_3_Option

    Max_1 = - HUGE( One )
    IF( PRESENT( Max_1_Option ) ) &
      Max_1 = Max_1_Option

    Max_2 = - HUGE( One )
    IF( PRESENT( Max_2_Option ) ) &
      Max_2 = Max_2_Option

    Max_3 = - HUGE( One )
    IF( PRESENT( Max_3_Option ) ) &
      Max_3 = Max_3_Option

    CALL InitializePositivityLimiter_Euler_Relativistic_TABLE &
           ( UsePositivityLimiter_Option = UsePositivityLimiter, &
             Verbose_Option              = Verbose, &
             Min_1_Option                = Min_1, &
             Min_2_Option                = Min_2, &
             Min_3_Option                = Min_3, &
             Max_1_Option                = Max_1, &
             Max_2_Option                = Max_2, &
             Max_3_Option                = Max_3 )

#elif defined( MICROPHYSICS_WEAKLIB ) && defined( HYDRO_NONRELATIVISTIC )

    Min_1 = - HUGE( One )
    IF( PRESENT( Min_1_Option ) ) &
      Min_1 = Min_1_Option

    Min_2 = - HUGE( One )
    IF( PRESENT( Min_2_Option ) ) &
      Min_2 = Min_2_Option

    Min_3 = - HUGE( One )
    IF( PRESENT( Min_3_Option ) ) &
      Min_3 = Min_3_Option

    Max_1 = - HUGE( One )
    IF( PRESENT( Max_1_Option ) ) &
      Max_1 = Max_1_Option

    Max_2 = - HUGE( One )
    IF( PRESENT( Max_2_Option ) ) &
      Max_2 = Max_2_Option

    Max_3 = - HUGE( One )
    IF( PRESENT( Max_3_Option ) ) &
      Max_3 = Max_3_Option

    CALL InitializePositivityLimiter_Euler_NonRelativistic_TABLE &
           ( UsePositivityLimiter_Option = UsePositivityLimiter, &
             Verbose_Option              = Verbose, &
             Min_1_Option                = Min_1, &
             Min_2_Option                = Min_2, &
             Min_3_Option                = Min_3, &
             Max_1_Option                = Max_1, &
             Max_2_Option                = Max_2, &
             Max_3_Option                = Max_3 )

#elif defined( HYDRO_RELATIVISTIC )

    Min_1 = - HUGE( One )
    IF( PRESENT( Min_1_Option ) ) &
      Min_1 = Min_1_Option

    Min_2 = - HUGE( One )
    IF( PRESENT( Min_2_Option ) ) &
      Min_2 = Min_2_Option

    CALL InitializePositivityLimiter_Euler_Relativistic_IDEAL &
           ( UsePositivityLimiter_Option = UsePositivityLimiter, &
             Verbose_Option              = Verbose, &
             Min_1_Option                = Min_1, &
             Min_2_Option                = Min_2 )

#else

    Min_1 = - HUGE( One )
    IF( PRESENT( Min_1_Option ) ) &
      Min_1 = Min_1_Option

    Min_2 = - HUGE( One )
    IF( PRESENT( Min_2_Option ) ) &
      Min_2 = Min_2_Option

    D_Min_Euler_PL = Zero
    IF( PRESENT( D_Min_Euler_PL_Option ) ) &
      D_Min_Euler_PL = D_Min_Euler_PL_Option

    IntE_Min_Euler_PL = Zero
    IF( PRESENT( IntE_Min_Euler_PL_Option ) ) &
      IntE_Min_Euler_PL = IntE_Min_Euler_PL_Option

    CALL InitializePositivityLimiter_Euler_NonRelativistic_IDEAL &
           ( UsePositivityLimiter_Option = UsePositivityLimiter, &
             Verbose_Option              = Verbose, &
             Min_1_Option                = Min_1, &
             Min_2_Option                = Min_2, &
             D_Min_Euler_PL_Option       = D_Min_Euler_PL , &
             IntE_Min_Euler_PL_Option    = IntE_Min_Euler_PL )

#endif

  END SUBROUTINE InitializePositivityLimiter_Euler


  SUBROUTINE FinalizePositivityLimiter_Euler

#if   defined( MICROPHYSICS_WEAKLIB ) && defined( HYDRO_RELATIVISTIC    )

    CALL FinalizePositivityLimiter_Euler_Relativistic_TABLE

#elif defined( MICROPHYSICS_WEAKLIB ) && defined( HYDRO_NONRELATIVISTIC )

    CALL FinalizePositivityLimiter_Euler_NonRelativistic_TABLE

#elif defined( HYDRO_RELATIVISTIC )

    CALL FinalizePositivityLimiter_Euler_Relativistic_IDEAL

#else

    CALL FinalizePositivityLimiter_Euler_NonRelativistic_IDEAL

#endif

  END SUBROUTINE FinalizePositivityLimiter_Euler


  SUBROUTINE ApplyPositivityLimiter_Euler &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D )

    INTEGER,  INTENT(in)              :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)              :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout)           :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout), OPTIONAL :: &
      D(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

#if   defined( MICROPHYSICS_WEAKLIB ) && defined( HYDRO_RELATIVISTIC    )

    CALL ApplyPositivityLimiter_Euler_Relativistic_TABLE &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U )

#elif defined( MICROPHYSICS_WEAKLIB ) && defined( HYDRO_NONRELATIVISTIC )

    CALL ApplyPositivityLimiter_Euler_NonRelativistic_TABLE &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D )

#elif defined( HYDRO_RELATIVISTIC )

    CALL ApplyPositivityLimiter_Euler_Relativistic_IDEAL &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U )

#else

    CALL ApplyPositivityLimiter_Euler_NonRelativistic_IDEAL &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U )

#endif

  END SUBROUTINE ApplyPositivityLimiter_Euler


END MODULE Euler_PositivityLimiterModule
