MODULE MHD_TallyModule

  USE KindModule, ONLY: &
    DP
  USE MagnetofluidFieldsModule, ONLY: &
    nCM
  USE MHD_TallyModule_Relativistic

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeTally_MHD
  PUBLIC :: FinalizeTally_MHD
  PUBLIC :: IncrementOffGridTally_MHD
  PUBLIC :: ComputeTally_MHD


CONTAINS


  SUBROUTINE InitializeTally_MHD &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, &
      SuppressTally_Option, BaseFileName_Option )

    INTEGER,  INTENT(in)           :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)           :: &
      G(:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    REAL(DP), INTENT(in)           :: &
      U(:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    LOGICAL,  INTENT(in), OPTIONAL :: &
      SuppressTally_Option
    CHARACTER(LEN=*),  INTENT(in), OPTIONAL :: &
      BaseFileName_Option

    LOGICAL        :: SuppressTally
    CHARACTER(256) :: BaseFileName

    SuppressTally = .FALSE.
    IF( PRESENT( SuppressTally_Option ) ) &
      SuppressTally = SuppressTally_Option

    BaseFileName = '../Output/'
    IF( PRESENT( BaseFileName_Option ) ) &
      BaseFileName = TRIM( BaseFileName_Option )

    CALL InitializeTally_MHD_Relativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, &
             SuppressTally_Option = SuppressTally, &
             BaseFileName_Option = TRIM( BaseFileName ) )

  END SUBROUTINE InitializeTally_MHD


  SUBROUTINE FinalizeTally_MHD

    CALL FinalizeTally_MHD_Relativistic

  END SUBROUTINE FinalizeTally_MHD


  SUBROUTINE ComputeTally_MHD &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, Time, &
      SetInitialValues_Option, Verbose_Option )

    INTEGER,  INTENT(in) :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in) :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(in) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(in) :: &
      Time
    LOGICAL,  INTENT(in), OPTIONAL :: &
      SetInitialValues_Option, Verbose_Option

    LOGICAL :: SetInitialValues, Verbose

    SetInitialValues = .FALSE.
    IF( PRESENT( SetInitialValues_Option ) ) &
      SetInitialValues = SetInitialValues_Option

    Verbose = .FALSE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    CALL ComputeTally_MHD_Relativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, Time, &
             SetInitialValues_Option = SetInitialValues, &
             Verbose_Option = Verbose )

  END SUBROUTINE ComputeTally_MHD


  SUBROUTINE IncrementOffGridTally_MHD( dM )

    REAL(DP), INTENT(in) :: dM(nCM)

    CALL IncrementOffGridTally_MHD_Relativistic( dM )

  END SUBROUTINE IncrementOffGridTally_MHD


END MODULE MHD_TallyModule
