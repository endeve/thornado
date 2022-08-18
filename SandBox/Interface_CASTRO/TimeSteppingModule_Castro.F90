#ifdef THORNADO_DEBUG
#define THORNADO_DEBUG_IMEX
#endif
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
    ComputeIncrement_TwoMoment_Implicit_New
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
      iS, iCR, iZ4, iZ3, iZ2, iZ1, iNode, iCF, iNodeX
    INTEGER  :: &
      iX_SW(3), iZ_SW(4), iZ_SW_P(4)
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

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: U_F, U_R, uGE, uGF ) &
    !$OMP MAP( alloc: U0_F, Q1_F, U0_R, T0_R, T1_R, Q1_R )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( U_F, U_R, uGE, uGF ) &
    !$ACC CREATE( U0_F, Q1_F, U0_R, T0_R, T1_R, Q1_R )
#endif

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
       iX_SW   = [    1, 1, 1 ]
       iZ_SW   = [ 0, 1, 1, 1 ] 
       iZ_SW_P = [ 0, 2, 2, 2 ]
    else
       iX_SW   = [    1, 1, 0 ]
       iZ_SW   = [ 0, 1, 1, 0 ]
       iZ_SW_P = [ 0, 2, 2, 0 ]
    end if

    IF( CallFromThornado )THEN
      iX_SW   = [ 0, 0, 0 ]     ! --- For Debugging within thornado
      iZ_SW   = [ 0, 0, 0, 0 ]  ! --- For Debugging within thornado
      iZ_SW_P = [ 0, 0, 0, 0 ]  ! --- For Debugging within thornado
    END IF

    ! --- Explicit Step (Radiation Only) ---

    IF( Explicit )THEN

      ! --- Apply Positivity Limiter ---

      CALL ApplyPositivityLimiter_TwoMoment &
             ( iZ_B0-iZ_SW_P, iZ_E0+iZ_SW_P, iZ_B1, iZ_E1, uGE, uGF, U_R )

      CALL ComputeIncrement_TwoMoment_Explicit &
             ( iZ_B0-iZ_SW, iZ_E0+iZ_SW, iZ_B1, iZ_E1, &
               uGE, uGF, &
               U_R, T0_R )

    ELSE

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7)
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) &
      !$ACC PRESENT( T0_R, iZ_B1, iZ_E1 )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO COLLAPSE(7)
#endif
      DO iS = 1, nSpecies
        DO iCR = 1, nCR
          DO iZ4 = iZ_B1(4), iZ_E1(4)
            DO iZ3 = iZ_B1(3), iZ_E1(3)
              DO iZ2 = iZ_B1(2), iZ_E1(2)
                DO iZ1 = iZ_B1(1), iZ_E1(1)
                  DO iNode = 1, nDOF
                    T0_R(iNode,iZ1,iZ2,iZ3,iZ4,iCR,iS) = Zero
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO

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
               U_F, Q1_F, &
               U_R, Q1_R )

    ELSE

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7)
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) &
      !$ACC PRESENT( Q1_R, iZ_B1, iZ_E1 )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO COLLAPSE(7)
#endif
      DO iS = 1, nSpecies
        DO iCR = 1, nCR
          DO iZ4 = iZ_B1(4), iZ_E1(4)
            DO iZ3 = iZ_B1(3), iZ_E1(3)
              DO iZ2 = iZ_B1(2), iZ_E1(2)
                DO iZ1 = iZ_B1(1), iZ_E1(1)
                  DO iNode = 1, nDOF
                    Q1_R(iNode,iZ1,iZ2,iZ3,iZ4,iCR,iS) = Zero
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
      !$ACC PRESENT( Q1_F, iZ_B1, iZ_E1 )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO COLLAPSE(7)
#endif
      DO iCF = 1, nCF
        DO iZ4 = iZ_B1(4), iZ_E1(4)
          DO iZ3 = iZ_B1(3), iZ_E1(3)
            DO iZ2 = iZ_B1(2), iZ_E1(2)
              DO iNodeX = 1, nDOFX
                Q1_F(iNodeX,iZ2,iZ3,iZ4,iCF) = Zero
              END DO
            END DO
          END DO
        END DO
      END DO

    END IF

    CALL AddFields_Fluid &
           ( iX_B0-iX_SW, iX_E0+iX_SW, One, dt, U_F, Q1_F, U_F )

    CALL AddFields_Radiation &
           ( iZ_B0-iZ_SW, iZ_E0+iZ_SW, One, dt, U_R, Q1_R, U_R )

    ! --- Apply Positivity Limiter ---

    CALL ApplyPositivityLimiter_TwoMoment &
           ( iZ_B0-iZ_SW, iZ_E0+iZ_SW, iZ_B1, iZ_E1, uGE, uGF, U_R )

    IF( .NOT. SingleStage ) THEN

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
                 U_R, T1_R )

      ELSE

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) &
        !$ACC PRESENT( T1_R, iZ_B1, iZ_E1 )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(7)
#endif
        DO iS = 1, nSpecies
          DO iCR = 1, nCR
            DO iZ4 = iZ_B1(4), iZ_E1(4)
              DO iZ3 = iZ_B1(3), iZ_E1(3)
                DO iZ2 = iZ_B1(2), iZ_E1(2)
                  DO iZ1 = iZ_B1(1), iZ_E1(1)
                    DO iNode = 1, nDOF
                      T1_R(iNode,iZ1,iZ2,iZ3,iZ4,iCR,iS) = Zero
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO

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
                 U_F, Q1_F, &
                 U_R, Q1_R )

      ELSE

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) &
        !$ACC PRESENT( Q1_R, iZ_B1, iZ_E1 )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(7)
#endif
        DO iS = 1, nSpecies
          DO iCR = 1, nCR
            DO iZ4 = iZ_B1(4), iZ_E1(4)
              DO iZ3 = iZ_B1(3), iZ_E1(3)
                DO iZ2 = iZ_B1(2), iZ_E1(2)
                  DO iZ1 = iZ_B1(1), iZ_E1(1)
                    DO iNode = 1, nDOF
                      Q1_R(iNode,iZ1,iZ2,iZ3,iZ4,iCR,iS) = Zero
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
        !$ACC PRESENT( Q1_F, iZ_B1, iZ_E1 )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(7)
#endif
        DO iCF = 1, nCF
          DO iZ4 = iZ_B1(4), iZ_E1(4)
            DO iZ3 = iZ_B1(3), iZ_E1(3)
              DO iZ2 = iZ_B1(2), iZ_E1(2)
                DO iNodeX = 1, nDOFX
                  Q1_F(iNodeX,iZ2,iZ3,iZ4,iCF) = Zero
                END DO
              END DO
            END DO
          END DO
        END DO

      END IF

      CALL AddFields_Fluid &
             ( iX_B0-iX_SW, iX_E0+iX_SW, One, Half * dt, U_F, Q1_F, U_F )

      CALL AddFields_Radiation &
             ( iZ_B0-iZ_SW, iZ_E0+iZ_SW, One, Half * dt, U_R, Q1_R, U_R )

      ! --- Apply Positivity Limiter ---

      CALL ApplyPositivityLimiter_TwoMoment &
             ( iZ_B0-iZ_SW, iZ_E0+iZ_SW, iZ_B1, iZ_E1, uGE, uGF, U_R )

    END IF

#ifdef THORNADO_DEBUG_IMEX
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE FROM( Q1_F, T0_R, T1_R, Q1_R )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE HOST( Q1_F, T0_R, T1_R, Q1_R )
#endif
    WRITE(*,'(a,8x,5i4,es23.15)') 'MINLOC(Q1_F), MINVAL(Q1_F)', MINLOC(Q1_F), MINVAL(Q1_F)
    WRITE(*,'(a,8x,5i4,es23.15)') 'MAXLOC(Q1_F), MAXVAL(Q1_F)', MAXLOC(Q1_F), MAXVAL(Q1_F)
    WRITE(*,'(a,7i4,es23.15)')    'MINLOC(T0_R), MINVAL(T0_R)', MINLOC(T0_R), MINVAL(T0_R)
    WRITE(*,'(a,7i4,es23.15)')    'MAXLOC(T0_R), MAXVAL(T0_R)', MAXLOC(T0_R), MAXVAL(T0_R)
    WRITE(*,'(a,7i4,es23.15)')    'MINLOC(T1_R), MINVAL(T1_R)', MINLOC(T1_R), MINVAL(T1_R)
    WRITE(*,'(a,7i4,es23.15)')    'MAXLOC(T1_R), MAXVAL(T1_R)', MAXLOC(T1_R), MAXVAL(T1_R)
    WRITE(*,'(a,7i4,es23.15)')    'MINLOC(Q1_R), MINVAL(Q1_R)', MINLOC(Q1_R), MINVAL(Q1_R)
    WRITE(*,'(a,7i4,es23.15)')    'MAXLOC(Q1_R), MAXVAL(Q1_R)', MAXLOC(Q1_R), MAXVAL(Q1_R)
#endif

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from: U_F, U_R ) &
    !$OMP MAP( release: U0_F, Q1_F, U0_R, T0_R, T1_R, Q1_R, uGE, uGF )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC COPYOUT( U_F, U_R ) &
    !$ACC DELETE( U0_F, Q1_F, U0_R, T0_R, T1_R, Q1_R, uGE, uGF )
#endif

  END SUBROUTINE Update_IMEX_PDARS


  SUBROUTINE AddFields_Fluid( iX_B, iX_E, alpha, beta, A, B, C )

    ! --- C = alpha * A + beta * B

    INTEGER,  INTENT(in)    :: iX_B(3), iX_E(3)
    REAL(DP), INTENT(in)    :: alpha, beta
    REAL(DP), INTENT(inout) :: &
      A(1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nCF)
    REAL(DP), INTENT(inout) :: &
      B(1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nCF)
    REAL(DP), INTENT(out)   :: &
      C(1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nCF)

    INTEGER :: iNodeX, iX1, iX2, iX3, iFF

    CALL TimersStart( Timer_AddFieldsF )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: A, B, C, iX_B, iX_E )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( A, B, C, iX_B, iX_E )
#endif

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( A, B, C, iX_B, iX_E )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iFF = 1, nCF
    DO iX3 = iX_B(3), iX_E(3)
    DO iX2 = iX_B(2), iX_E(2)
    DO iX1 = iX_B(1), iX_E(1)
    DO iNodeX = 1, nDOFX

      C(iNodeX,iX1,iX2,iX3,iFF) &
        = alpha * A(iNodeX,iX1,iX2,iX3,iFF) &
            + beta * B(iNodeX,iX1,iX2,iX3,iFF)

    END DO
    END DO
    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from: A, B, C ) &
    !$OMP MAP( release: iX_B, iX_E )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC COPYOUT( A, B, C ) &
    !$ACC DELETE( iX_B, iX_E )
#endif

    CALL TimersStop( Timer_AddFieldsF )

  END SUBROUTINE AddFields_Fluid


  SUBROUTINE AddFields_Radiation( iZ_B, iZ_E, alpha, beta, A, B, C )

    ! --- C = alpha * A + beta * B

    INTEGER,  INTENT(in)    :: iZ_B(4), iZ_E(4)
    REAL(DP), INTENT(in)    :: alpha, beta
    REAL(DP), INTENT(inout) :: &
      A(1:nDOF,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(inout) :: &
      B(1:nDOF,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(out)   :: &
      C(1:nDOF,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)

    INTEGER :: iNode, iZ1, iZ2, iZ3, iZ4, iRF, iS

    CALL TimersStart( Timer_AddFieldsR )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: A, B, C, iZ_B, iZ_E )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( A, B, C, iZ_B, iZ_E )
#endif

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) &
    !$ACC PRESENT( A, B, C, iZ_B, iZ_E )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(7)
#endif
    DO iS = 1, nSpecies
    DO iRF = 1, nCR
    DO iZ4 = iZ_B(4), iZ_E(4)
    DO iZ3 = iZ_B(3), iZ_E(3)
    DO iZ2 = iZ_B(2), iZ_E(2)
    DO iZ1 = iZ_B(1), iZ_E(1)
    DO iNode = 1, nDOF

      C(iNode,iZ1,iZ2,iZ3,iZ4,iRF,iS) &
        = alpha * A(iNode,iZ1,iZ2,iZ3,iZ4,iRF,iS) &
            + beta * B(iNode,iZ1,iZ2,iZ3,iZ4,iRF,iS)

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from: A, B, C ) &
    !$OMP MAP( release: iZ_B, iZ_E )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC COPYOUT( A, B, C ) &
    !$ACC DELETE( iZ_B, iZ_E )
#endif

    CALL TimersStop( Timer_AddFieldsR )

  END SUBROUTINE AddFields_Radiation


END MODULE TimeSteppingModule_Castro
