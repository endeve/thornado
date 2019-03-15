MODULE TimeSteppingModule_Castro

  USE KindModule, ONLY: &
    DP, Zero, Half, One
  USE ProgramHeaderModule, ONLY: &
    nDOFX, nDOFE, nDOF, &
    iX_B0, iX_B1, iX_E0, iX_E1, &
    iZ_B0, iZ_B1, iZ_E0, iZ_E1
  USE TimersModule, ONLY: &
    TimersStart, &
    TimersStop, &
    Timer_AddFieldsF, &
    Timer_AddFieldsR
  USE FluidFieldsModule, ONLY: &
    nCF, iCF_Ne
  USE RadiationFieldsModule, ONLY: &
    nCR, nSpecies
  USE TwoMoment_DiscretizationModule_Streaming, ONLY: &
    ComputeIncrement_TwoMoment_Explicit
  USE TwoMoment_DiscretizationModule_Collisions_Neutrinos, ONLY: &
    ComputeIncrement_TwoMoment_Implicit, &
    ComputeIncrement_TwoMoment_Implicit_New, &
    ComputeIncrement_TwoMoment_Implicit_DGFV
  USE TwoMoment_PositivityLimiterModule, ONLY: &
    ApplyPositivityLimiter_TwoMoment

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeTimeStep_TwoMoment
  PUBLIC :: Update_IMEX_PDARS

CONTAINS


  SUBROUTINE ComputeTimeStep_TwoMoment( dX_CGS, dt_CGS )

    use PhysicalConstantsModule, only : SpeedOfLightCGS
    use ProgramHeaderModule,     only : nDimsX, nNodes

    REAL(DP), INTENT(in)  :: dX_CGS(3)
    REAL(DP), INTENT(out) :: dt_CGS

    INTEGER  :: i
    REAL(DP) :: CFL, InverseLenghtScale

    CFL = One / ( SpeedOfLightCGS * DBLE( 2 * nNodes - 1 ) )

    InverseLenghtScale = 0.0_DP
    DO i = 1, nDimsX
      InverseLenghtScale &
        = InverseLenghtScale + 1.0_DP / dX_CGS(i)
    END DO

    dt_CGS = CFL / InverseLenghtScale

  END SUBROUTINE ComputeTimeStep_TwoMoment


  SUBROUTINE Update_IMEX_PDARS &
    ( dt, U_F, U_R, Explicit_Option, Implicit_Option, &
      SingleStage_Option, CallFromThornado_Option )

    use GeometryFieldsModuleE, only : uGE
    use GeometryFieldsModule,  only : uGF
    use ProgramHeaderModule  , only : nDimsX

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    REAL(DP), INTENT(in)    :: &
      dt
    REAL(DP), INTENT(inout) :: &
      U_F(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCF)
    REAL(DP), INTENT(inout) :: &
      U_R(1:nDOF ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    LOGICAL,  INTENT(in), OPTIONAL :: &
      Explicit_Option, &
      Implicit_Option, &
      SingleStage_Option, &
      CallFromThornado_Option

    LOGICAL  :: &
      Explicit, &
      Implicit, &
      SingleStage, &
      CallFromThornado
    INTEGER  :: &
      iX_SW(3), iZ_SW(4)
    REAL(DP) :: &
      U0_F &
        (1:nDOFX, &
         iZ_B1(2):iZ_E1(2), &
         iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4), &
         1:nCF), &
      Q1_F &
        (1:nDOFX, &
         iZ_B1(2):iZ_E1(2), &
         iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4), &
         1:nCF)
    REAL(DP) :: &
      U0_R &
        (1:nDOF, &
         iZ_B1(1):iZ_E1(1), &
         iZ_B1(2):iZ_E1(2), &
         iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4), &
         1:nCR,1:nSpecies), &
      T0_R &
        (1:nDOF, &
         iZ_B1(1):iZ_E1(1), &
         iZ_B1(2):iZ_E1(2), &
         iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4), &
         1:nCR,1:nSpecies), &
      T1_R &
        (1:nDOF, &
         iZ_B1(1):iZ_E1(1), &
         iZ_B1(2):iZ_E1(2), &
         iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4), &
         1:nCR,1:nSpecies), &
      Q1_R &
        (1:nDOF, &
         iZ_B1(1):iZ_E1(1), &
         iZ_B1(2):iZ_E1(2), &
         iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4), &
         1:nCR,1:nSpecies)

    IF( PRESENT( Explicit_Option ) )THEN
      Explicit = Explicit_Option
    ELSE
      Explicit = .TRUE.
    END IF

    IF( PRESENT( Implicit_Option ) )THEN
      Implicit = Implicit_Option
    ELSE
      Implicit = .TRUE.
    END IF

    IF( PRESENT( SingleStage_Option ) )THEN
      SingleStage = SingleStage_Option
    ELSE
      SingleStage = .FALSE.
    END IF

    IF( PRESENT( CallFromThornado_Option ) )THEN
      CallFromThornado = CallFromThornado_Option
    ELSE
      CallFromThornado = .FALSE.
    END IF

    ! ----------------------------------------------------------------
    ! --- Positive, Diffusion Accurate IMEX Scheme from Chu et al. ---
    ! --- arXiv:1809.06949 -------------------------------------------
    ! ----------------------------------------------------------------

    CALL AddFields_Fluid &
           ( iX_B1, iX_E1, One, Zero, U_F, U_F, U0_F )

    CALL AddFields_Radiation &
           ( iZ_B1, iZ_E1, One, Zero, U_R, U_R, U0_R )

    ! ---------------
    ! --- Stage 1 ---
    ! ---------------

    ! --- Include One Layer of Spatial Ghost Cells in Update

    if (nDimsX .eq. 3) then
       iX_SW = [    1, 1, 1 ]
       iZ_SW = [ 0, 1, 1, 1 ]
    else
       iX_SW = [    1, 1, 0 ]
       iZ_SW = [ 0, 1, 1, 0 ]
    end if

    IF( CallFromThornado )THEN
      iX_SW = [ 0, 0, 0 ]     ! --- For Debugging within thornado
      iZ_SW = [ 0, 0, 0, 0 ]  ! --- For Debugging within thornado
    END IF

    ! --- Explicit Step (Radiation Only) ---

    IF( Explicit )THEN

      ! --- Apply Positivity Limiter ---

      CALL ApplyPositivityLimiter_TwoMoment &
             ( iZ_B0-iZ_SW, iZ_E0+iZ_SW, iZ_B1, iZ_E1, uGE, uGF, U_R )

      CALL ComputeIncrement_TwoMoment_Explicit &
             ( iZ_B0-iZ_SW, iZ_E0+iZ_SW, iZ_B1, iZ_E1, &
               uGE, uGF, &
               U_R, T0_R &
                      (1:nDOF, &
                       iZ_B0(1)-iZ_SW(1):iZ_E0(1)+iZ_SW(1), &
                       iZ_B0(2)-iZ_SW(2):iZ_E0(2)+iZ_SW(2), &
                       iZ_B0(3)-iZ_SW(3):iZ_E0(3)+iZ_SW(3), &
                       iZ_B0(4)-iZ_SW(4):iZ_E0(4)+iZ_SW(4), &
                       1:nCR,1:nSpecies) )

    ELSE

      T0_R = Zero

    END IF

    ! --- Apply Increment ---

    CALL AddFields_Radiation &
           ( iZ_B0-iZ_SW, iZ_E0+iZ_SW, One, dt, U0_R, T0_R, U_R )

    ! --- Apply Positivity Limiter ---

    CALL ApplyPositivityLimiter_TwoMoment &
           ( iZ_B0-iZ_SW, iZ_E0+iZ_SW, iZ_B1, iZ_E1, uGE, uGF, U_R )

    ! --- Implicit Step ---

    IF( Implicit )THEN

      CALL ComputeIncrement_TwoMoment_Implicit_New &
             ( iZ_B0-iZ_SW, iZ_E0+iZ_SW, iZ_B1, iZ_E1, dt, &
               uGE, uGF, &
               U_F, Q1_F &
                      (1:nDOFX, &
                       iX_B0(1)-iX_SW(1):iX_E0(1)+iX_SW(1), &
                       iX_B0(2)-iX_SW(2):iX_E0(2)+iX_SW(2), &
                       iX_B0(3)-iX_SW(3):iX_E0(3)+iX_SW(3), &
                       1:nCF), &
               U_R, Q1_R &
                      (1:nDOF, &
                       iZ_B0(1)-iZ_SW(1):iZ_E0(1)+iZ_SW(1), &
                       iZ_B0(2)-iZ_SW(2):iZ_E0(2)+iZ_SW(2), &
                       iZ_B0(3)-iZ_SW(3):iZ_E0(3)+iZ_SW(3), &
                       iZ_B0(4)-iZ_SW(4):iZ_E0(4)+iZ_SW(4), &
                       1:nCR,1:nSpecies) )

    ELSE

      Q1_F = Zero
      Q1_R = Zero

    END IF

    CALL AddFields_Fluid &
           ( iX_B0-iX_SW, iX_E0+iX_SW, One, dt, U_F, Q1_F, U_F )

    CALL AddFields_Radiation &
           ( iZ_B0-iZ_SW, iZ_E0+iZ_SW, One, dt, U_R, Q1_R, U_R )

    ! --- Apply Positivity Limiter ---

    CALL ApplyPositivityLimiter_TwoMoment &
           ( iZ_B0-iZ_SW, iZ_E0+iZ_SW, iZ_B1, iZ_E1, uGE, uGF, U_R )


    IF( SingleStage ) RETURN

    ! ---------------
    ! --- Stage 2 ---
    ! ---------------

    iX_SW = [ 0, 0, 0 ]
    iZ_SW = [ 0, 0, 0, 0 ]

    ! --- Explicit Step (Radiation Only) ---

    IF( Explicit )THEN

      CALL ComputeIncrement_TwoMoment_Explicit &
             ( iZ_B0-iZ_SW, iZ_E0+iZ_SW, iZ_B1, iZ_E1, &
               uGE, uGF, &
               U_R, T1_R &
                      (1:nDOF, &
                       iZ_B0(1)-iZ_SW(1):iZ_E0(1)+iZ_SW(1), &
                       iZ_B0(2)-iZ_SW(2):iZ_E0(2)+iZ_SW(2), &
                       iZ_B0(3)-iZ_SW(3):iZ_E0(3)+iZ_SW(3), &
                       iZ_B0(4)-iZ_SW(4):iZ_E0(4)+iZ_SW(4), &
                       1:nCR, 1:nSpecies) )

    ELSE

      T1_R = 0.0_DP

    END IF

    ! --- Apply Increment ---

    CALL AddFields_Fluid &
           ( iX_B0-iX_SW, iX_E0+iX_SW, One, Half * dt, U0_F, Q1_F, U_F )

    CALL AddFields_Radiation &
           ( iZ_B0-iZ_SW, iZ_E0+iZ_SW, One, Half * dt, U0_R, T0_R, U_R )

    CALL AddFields_Radiation &
           ( iZ_B0-iZ_SW, iZ_E0+iZ_SW, One, Half * dt, U_R,  T1_R, U_R )

    CALL AddFields_Radiation &
           ( iZ_B0-iZ_SW, iZ_E0+iZ_SW, One, Half * dt, U_R,  Q1_R, U_R )

    ! --- Apply Positivity Limiter ---

    CALL ApplyPositivityLimiter_TwoMoment &
           ( iZ_B0-iZ_SW, iZ_E0+iZ_SW, iZ_B1, iZ_E1, uGE, uGF, U_R )

    ! --- Implicit Step ---

    IF( Implicit )THEN

      CALL ComputeIncrement_TwoMoment_Implicit_New &
             ( iZ_B0-iZ_SW, iZ_E0+iZ_SW, iZ_B1, iZ_E1, Half * dt, &
               uGE, uGF, &
               U_F, Q1_F &
                      (1:nDOFX, &
                       iX_B0(1)-iX_SW(1):iX_E0(1)+iX_SW(1), &
                       iX_B0(2)-iX_SW(2):iX_E0(2)+iX_SW(2), &
                       iX_B0(3)-iX_SW(3):iX_E0(3)+iX_SW(3), &
                       1:nCF), &
               U_R, Q1_R &
                      (1:nDOF, &
                       iZ_B0(1)-iZ_SW(1):iZ_E0(1)+iZ_SW(1), &
                       iZ_B0(2)-iZ_SW(2):iZ_E0(2)+iZ_SW(2), &
                       iZ_B0(3)-iZ_SW(3):iZ_E0(3)+iZ_SW(3), &
                       iZ_B0(4)-iZ_SW(4):iZ_E0(4)+iZ_SW(4), &
                       1:nCR, 1:nSpecies) )

    ELSE

      Q1_F = Zero
      Q1_R = Zero

    END IF

    CALL AddFields_Fluid &
           ( iX_B0-iX_SW, iX_E0+iX_SW, One, Half * dt, U_F, Q1_F, U_F )

    CALL AddFields_Radiation &
           ( iZ_B0-iZ_SW, iZ_E0+iZ_SW, One, Half * dt, U_R, Q1_R, U_R )

    ! --- Apply Positivity Limiter ---

    CALL ApplyPositivityLimiter_TwoMoment &
           ( iZ_B0-iZ_SW, iZ_E0+iZ_SW, iZ_B1, iZ_E1, uGE, uGF, U_R )

  END SUBROUTINE Update_IMEX_PDARS


  SUBROUTINE AddFields_Fluid( iX_B, iX_E, alpha, beta, A, B, C )

    ! --- C = alpha * A + beta * B

    INTEGER,  INTENT(in)    :: iX_B(3), iX_E(3)
    REAL(DP), INTENT(in)    :: alpha, beta
    REAL(DP), INTENT(inout) :: &
      A(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: &
      B(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(out)   :: &
      C(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER :: iX1, iX2, iX3, iFF

    CALL TimersStart( Timer_AddFieldsF )

    DO iFF = 1, nCF
    DO iX3 = iX_B(3), iX_E(3)
    DO iX2 = iX_B(2), iX_E(2)
    DO iX1 = iX_B(1), iX_E(1)

      C(:,iX1,iX2,iX3,iFF) &
        = alpha * A(:,iX1,iX2,iX3,iFF) &
            + beta * B(:,iX1,iX2,iX3,iFF)

    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_AddFieldsF )

  END SUBROUTINE AddFields_Fluid


  SUBROUTINE AddFields_Radiation( iZ_B, iZ_E, alpha, beta, A, B, C )

    ! --- C = alpha * A + beta * B

    INTEGER,  INTENT(in)    :: iZ_B(4), iZ_E(4)
    REAL(DP), INTENT(in)    :: alpha, beta
    REAL(DP), INTENT(inout) :: &
      A(1:,iZ_B1(1):,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:,1:)
    REAL(DP), INTENT(inout) :: &
      B(1:,iZ_B1(1):,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:,1:)
    REAL(DP), INTENT(out)   :: &
      C(1:,iZ_B1(1):,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:,1:)

    INTEGER :: iZ1, iZ2, iZ3, iZ4, iRF, iS

    CALL TimersStart( Timer_AddFieldsR )

    DO iS = 1, nSpecies
    DO iRF = 1, nCR
    DO iZ4 = iZ_B(4), iZ_E(4)
    DO iZ3 = iZ_B(3), iZ_E(3)
    DO iZ2 = iZ_B(2), iZ_E(2)
    DO iZ1 = iZ_B(1), iZ_E(1)

      C(:,iZ1,iZ2,iZ3,iZ4,iRF,iS) &
        = alpha * A(:,iZ1,iZ2,iZ3,iZ4,iRF,iS) &
            + beta * B(:,iZ1,iZ2,iZ3,iZ4,iRF,iS)

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_AddFieldsR )

  END SUBROUTINE AddFields_Radiation


END MODULE TimeSteppingModule_Castro
