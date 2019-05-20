MODULE TwoMoment_SlopeLimiterModule

  USE KindModule, ONLY: &
    DP, Zero, One
  USE ProgramHeaderModule, ONLY: &
    nNodesZ, nDOF, nDOFE, nDOFX
  USE GeometryFieldsModule, ONLY: &
    nGF
  USE GeometryFieldsModuleE, ONLY: &
    nGE
  USE RadiationFieldsModule, ONLY: &
    nSpecies, nCR, &
    iCR_N, iCR_G1, iCR_G2, iCR_G3

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: TwoMoment_InitializeSlopeLimiter
  PUBLIC :: TwoMoment_FinalizeSlopeLimiter
  PUBLIC :: TwoMoment_ApplySlopeLimiter

  LOGICAL  :: UseSlopeLimiter
  LOGICAL  :: UseCharacteristicLimiting
  REAL(DP) :: SlopeTolerance
  REAL(DP) :: I_4x4(1:4,1:4)

CONTAINS


  SUBROUTINE TwoMoment_InitializeSlopeLimiter &
    ( SlopeTolerance_Option, &
      UseSlopeLimiter_Option, &
      UseCharacteristicLimiting_Option, &
      Verbose_Option ) 

    REAL(DP), INTENT(in), OPTIONAL :: &
      SlopeTolerance_Option
    LOGICAL, INTENT(in), OPTIONAL :: &
      UseSlopeLimiter_Option, &
      UseCharacteristicLimiting_Option, &
      Verbose_Option

    LOGICAL  :: Verbose
    INTEGER  :: i

    IF( PRESENT( SlopeTolerance_Option ) )THEN
      SlopeTolerance = SlopeTolerance_Option
    ELSE
      SlopeTolerance = 1.0d-3
    END IF

    IF( PRESENT( UseSlopeLimiter_Option ) )THEN
      UseSlopeLimiter = UseSlopeLimiter_Option
    ELSE
      UseSlopeLimiter = .TRUE.
    END IF

    IF( PRESENT( UseCharacteristicLimiting_Option ) )THEN
      UseCharacteristicLimiting = UseCharacteristicLimiting_Option
    ELSE
      UseCharacteristicLimiting = .TRUE.
    END IF

    IF( PRESENT( Verbose_Option ) )THEN
      Verbose = Verbose_Option
    ELSE
      Verbose = .TRUE.
    END IF

    I_4x4 = Zero
    DO i = 1, 4
      I_4x4(i,i) = One
    END DO

    IF( Verbose )THEN
      WRITE(*,*)
      WRITE(*,'(A)') '  INFO: TwoMoment_InitializeSlopeLimiter:'
      WRITE(*,'(A)') '  -----------------------------------'
      WRITE(*,*)
      WRITE(*,'(A4,A27,L1)'      ) '', 'UseSlopeLimiter: ' , &
        UseSlopeLimiter
      WRITE(*,*)
      WRITE(*,'(A4,A27,ES9.3E2)' ) '', 'SlopeTolerance: ' , &
        SlopeTolerance
      WRITE(*,'(A4,A27,L1)'      ) '', 'UseCharacteristicLimiting: ' , &
        UseCharacteristicLimiting
      WRITE(*,*)
    END IF
  
  END SUBROUTINE TwoMoment_InitializeSlopeLimiter


  SUBROUTINE TwoMoment_FinalizeSlopeLimiter

  END SUBROUTINE TwoMoment_FinalizeSlopeLimiter


  SUBROUTINE TwoMoment_ApplySlopeLimiter &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U, SuppressBC_Option )

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)    :: &
      GE(1:nDOFE,iZ_B1(1):iZ_E1(1),1:nGE)
    REAL(DP), INTENT(in)    :: &
      GX(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(inout) :: &
      U (1:nDOF ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    LOGICAL,  INTENT(in), OPTIONAL :: &
      SuppressBC_Option

    LOGICAL  :: SuppressBC

    IF( .NOT. UseSlopeLimiter ) RETURN

    IF( nDOFE == 1 ) RETURN

  END SUBROUTINE TwoMoment_ApplySlopeLimiter


END MODULE TwoMoment_SlopeLimiterModule
