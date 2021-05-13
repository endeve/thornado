MODULE MF_Euler_TallyModule

  ! --- AMReX Modules ---

  USE amrex_fort_module, ONLY: &
    AR => amrex_real
  USE amrex_geometry_module, ONLY: &
    amrex_geometry
  USE amrex_multifab_module, ONLY: &
    amrex_multifab
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nDOFX
  USE GeometryFieldsModule, ONLY: &
    nGF
  USE FluidFieldsModule, ONLY: &
    nCF
  USE Euler_TallyModule, ONLY: &
    InitializeTally_Euler, &
    ComputeTally_Euler, &
    IncrementOffGridTally_Euler, &
    FinalizeTally_Euler

  ! --- Local Modules ---

  USE InputParsingModule, ONLY: &
    nX, &
    swX, &
    nLevels
  USE MF_UtilitiesModule, ONLY: &
    amrex2thornado_X_Global

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MF_InitializeTally_Euler
  PUBLIC :: MF_ComputeTally_Euler
  PUBLIC :: MF_IncrementOffGridTally_Euler
  PUBLIC :: MF_FinalizeTally_Euler


CONTAINS


  SUBROUTINE MF_InitializeTally_Euler( BaseFileName_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: BaseFileName_Option

    CHARACTER(256) :: BaseFileName
    INTEGER        :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)

    ! --- Dummy arrays ---

    REAL(AR) :: G(nDOFX,1-swX(1):nX(1)+swX(1), &
                        1-swX(2):nX(2)+swX(2), &
                        1-swX(3):nX(3)+swX(3), &
                  nGF)
    REAL(AR) :: U(nDOFX,1-swX(1):nX(1)+swX(1), &
                        1-swX(2):nX(2)+swX(2), &
                        1-swX(3):nX(3)+swX(3), &
                  nCF)

    BaseFileName = './'
    IF( PRESENT( BaseFileName_Option ) ) &
      BaseFileName = TRIM( BaseFileName_Option )

    iX_B0 = [ 1, 1, 1 ]
    iX_E0 = nX
    iX_B1 = 1 - swX
    iX_E1 = nX + swX

    IF( amrex_parallel_ioprocessor() ) &
      CALL InitializeTally_Euler &
             ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, &
               BaseFileName_Option = TRIM( BaseFileName ) )

  END SUBROUTINE MF_InitializeTally_Euler


  SUBROUTINE MF_ComputeTally_Euler &
    ( GEOM, MF_uGF, MF_uCF, Time, SetInitialValues_Option, Verbose_Option )

    TYPE(amrex_geometry), INTENT(in) :: GEOM  (0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in) :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in) :: MF_uCF(0:nLevels-1)
    REAL(AR),             INTENT(in) :: Time
    LOGICAL,              INTENT(in), OPTIONAL :: SetInitialValues_Option
    LOGICAL,              INTENT(in), OPTIONAL :: Verbose_Option

    LOGICAL :: SetInitialValues, Verbose
    INTEGER :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)

    REAL(AR) :: G(nDOFX,1-swX(1):nX(1)+swX(1), &
                        1-swX(2):nX(2)+swX(2), &
                        1-swX(3):nX(3)+swX(3), &
                  nGF)
    REAL(AR) :: U(nDOFX,1-swX(1):nX(1)+swX(1), &
                        1-swX(2):nX(2)+swX(2), &
                        1-swX(3):nX(3)+swX(3), &
                  nCF)

    SetInitialValues = .FALSE.
    IF( PRESENT( SetInitialValues_Option ) ) &
      SetInitialValues = SetInitialValues_Option

    Verbose = .FALSE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    iX_B0 = [ 1, 1, 1 ]
    iX_E0 = nX
    iX_B1 = 1 - swX
    iX_E1 = nX + swX

    CALL amrex2thornado_X_Global( GEOM, MF_uGF, nGF, G )
    CALL amrex2thornado_X_Global( GEOM, MF_uCF, nCF, U )

    IF( amrex_parallel_ioprocessor() ) &
      CALL ComputeTally_Euler &
             ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, Time, &
               SetInitialValues_Option = SetInitialValues, &
               Verbose_Option = Verbose )

  END SUBROUTINE MF_ComputeTally_Euler


  SUBROUTINE MF_IncrementOffGridTally_Euler( dM )

    REAL(AR), INTENT(in) :: dM(nCF)

    IF( amrex_parallel_ioprocessor() ) &
      CALL IncrementOffGridTally_Euler( dM )

  END SUBROUTINE MF_IncrementOffGridTally_Euler


  SUBROUTINE MF_FinalizeTally_Euler
  END SUBROUTINE MF_FinalizeTally_Euler


END MODULE MF_Euler_TallyModule
