#ifdef THORNADO_DEBUG
#define THORNADO_DEBUG_IMEX
#endif
MODULE TimeSteppingModule_Flash

  USE KindModule, ONLY: &
    DP, Zero, Half, One
  USE ProgramHeaderModule, ONLY: &
    nDOFX, nDOFE, nDOF, nDOFZ, nNodesZ, &
    iX_B0, iX_B1, iX_E0, iX_E1, nNodesX, &
    iZ_B0, iZ_B1, iZ_E0, iZ_E1
  USE ReferenceElementModule, ONLY: &
    NodeNumberTable4D
  USE UtilitiesModule, ONLY: &
    NodeNumberX
  USE TimersModule, ONLY: &
    TimersStart, &
    TimersStop, &
    Timer_AddFieldsF, &
    Timer_AddFieldsR
  USE FluidFieldsModule, ONLY: &
    nCF, uDF, iCF_D, iCF_S1, iCF_S2, iCF_S3, &
    iCF_E, iCF_Ne
  USE RadiationFieldsModule, ONLY: &
    nCR, nSpecies, iCR_N, iCR_G1, iCR_G2, iCR_G3
  USE TwoMoment_DiscretizationModule_Streaming, ONLY: &
    ComputeIncrement_TwoMoment_Explicit
  USE TwoMoment_DiscretizationModule_Collisions_Neutrinos, ONLY: &
    ComputeIncrement_TwoMoment_Implicit
  USE TwoMoment_PositivityLimiterModule, ONLY: &
    ApplyPositivityLimiter_TwoMoment
  USE TwoMoment_SlopeLimiterModule, ONLY : &
    ApplySlopeLimiter_TwoMoment
  USE Euler_PositivityLimiterModule_NonRelativistic_TABLE, ONLY: &
    ApplyPositivityLimiter_Euler_NonRelativistic_TABLE
#ifdef TWOMOMENT_ORDER_V
  USE TwoMoment_DiscretizationModule_Streaming, ONLY: &
    OffGridFlux_TwoMoment
#endif

  USE, INTRINSIC :: ieee_arithmetic, ONLY: &
    IEEE_IS_NAN

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeTimeStep_TwoMoment
  PUBLIC :: Update_IMEX_PDARS
  PUBLIC :: ApplyBoundaryConditions_Radiation
  PUBLIC :: ApplyBoundaryConditions_Fluid

  LOGICAL :: DEBUG = .FALSE.

  INTEGER, PARAMETER, PUBLIC :: iApplyBC_Both  = 0
  INTEGER, PARAMETER, PUBLIC :: iApplyBC_Inner = 1
  INTEGER, PARAMETER, PUBLIC :: iApplyBC_Outer = 2
  INTEGER, PARAMETER, PUBLIC :: iApplyBC_None  = 3

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
      SingleStage_Option, CallFromThornado_Option, &
      bcX_Option, iApplyBC_Option, OffGridFluxR_Option )

    use GeometryFieldsModuleE, only : uGE
    use GeometryFieldsModule,  only : uGF
    use ProgramHeaderModule  , only : nDimsX

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    REAL(DP), INTENT(in)    :: &
      dt
    REAL(DP), INTENT(inout) :: &
      U_F(1:nDOFX, &
          iZ_B1(2):iZ_E1(2), &
          iZ_B1(3):iZ_E1(3), &
          iZ_B1(4):iZ_E1(4), &
          1:nCF)
    REAL(DP), INTENT(inout) :: &
      U_R(1:nDOFZ, &
          iZ_B1(1):iZ_E1(1), &
          iZ_B1(2):iZ_E1(2), &
          iZ_B1(3):iZ_E1(3), &
          iZ_B1(4):iZ_E1(4), &
          1:nCR,1:nSpecies)
    LOGICAL,  INTENT(in), OPTIONAL :: &
      Explicit_Option, &
      Implicit_Option, &
      SingleStage_Option, &
      CallFromThornado_Option
    INTEGER, INTENT(in), OPTIONAL :: &
      bcX_Option(3), iApplyBC_Option(3)
    REAL(DP), INTENT(out), OPTIONAL :: &
      OffGridFluxR_Option(2*nCR)

    LOGICAL  :: &
      Explicit, &
      Implicit, &
      SingleStage, &
      CallFromThornado
    INTEGER  :: &
      iS, iCR, iZ4, iZ3, iZ2, iZ1, iNode, iCF, iNodeX, bcX(3), iApplyBC(3)
    INTEGER  :: &
      iX_SW(3), iZ_SW(4), iZ_SW_P(4)
    INTEGER  :: &
      iX_B0_SW(3), iX_E0_SW(3), iZ_B0_SW(4), iZ_E0_SW(4), iZ_B0_SW_P(4), iZ_E0_SW_P(4)
    REAL(DP) :: &
      OffGridFluxR   (2*nCR), &
      OffGridFluxR_T0(2*nCR), &
      OffGridFluxR_T1(2*nCR)
    REAL(DP), ALLOCATABLE, DIMENSION(:,:,:,:,:)     :: U0_F, Q1_F
    REAL(DP), ALLOCATABLE, DIMENSION(:,:,:,:,:,:,:) :: U0_R, T0_R, T1_R, Q1_R

    ALLOCATE( U0_F(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCF) )
    ALLOCATE( Q1_F(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCF) )

    ALLOCATE( U0_R(1:nDOF,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies) )
    ALLOCATE( T0_R(1:nDOF,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies) )
    ALLOCATE( T1_R(1:nDOF,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies) )
    ALLOCATE( Q1_R(1:nDOF,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies) )

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

    IF( PRESENT( bcX_Option ) )THEN
      bcX = bcX_Option
    ELSE
      bcX = 0 ! No Boundary Condition
    END IF

    IF( PRESENT( iApplyBC_Option ) )THEN
      iApplyBC = iApplyBC_Option
    ELSE
      iApplyBC = iApplyBC_None
    END IF

#ifdef THORNADO_DEBUG
    IF( ANY(IEEE_IS_NAN(U_R)) ) PRINT*, 'NaN when enter TimeStep'
#endif

    U0_F = Zero; Q1_F = Zero

    U0_R = Zero; T0_R = Zero; T1_R = Zero; Q1_R = Zero

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: U_F, U_R, uGE, uGF, &
    !$OMP          U0_F, Q1_F, U0_R, T0_R, T1_R, Q1_R )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( U_F, U_R, uGE, uGF, &
    !$ACC         U0_F, Q1_F, U0_R, T0_R, T1_R, Q1_R )
#endif

    OffGridFluxR = Zero

#ifdef TWOMOMENT_ORDER_V
#if defined MICROPHYSICS_WEAKLIB
    CALL ApplyPositivityLimiter_Euler_NonRelativistic_TABLE &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, U_F, uDF )
#endif
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
    else if (nDimsX .eq. 2) then
       iX_SW   = [    1, 1, 0 ]
       iZ_SW   = [ 0, 1, 1, 0 ]
       iZ_SW_P = [ 0, 2, 2, 0 ]
    else if (nDimsX .eq. 1) then
       iX_SW   = [    1, 0, 0 ]
       iZ_SW   = [ 0, 1, 0, 0 ]
       iZ_SW_P = [ 0, 2, 0, 0 ]
    end if

    IF( CallFromThornado )THEN
      iX_SW   = [ 0, 0, 0 ]     ! --- For Debugging within thornado
      iZ_SW   = [ 0, 0, 0, 0 ]  ! --- For Debugging within thornado
      iZ_SW_P = [ 0, 0, 0, 0 ]  ! --- For Debugging within thornado
    END IF

    iX_B0_SW = iX_B0 - iX_SW
    iX_E0_SW = iX_E0 + iX_SW
    iZ_B0_SW = iZ_B0 - iZ_SW
    iZ_E0_SW = iZ_E0 + iZ_SW
    iZ_B0_SW_P = iZ_B0 - iZ_SW_P
    iZ_E0_SW_P = iZ_E0 + iZ_SW_P

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: iX_B0_SW, iX_E0_SW, iZ_B0_SW, iZ_E0_SW, iZ_B0_SW_P, iZ_E0_SW_P, iZ_SW_P )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( iX_B0_SW, iX_E0_SW, iZ_B0_SW, iZ_E0_SW, iZ_B0_SW_P, iZ_E0_SW_P, iZ_SW_P )
#endif

    ! --- Explicit Step (Radiation Only) ---

    IF( Explicit )THEN

      ! --- Apply Limiter ---

#ifdef TWOMOMENT_ORDER_1
      CALL ApplyPositivityLimiter_TwoMoment &
             ( iZ_B0_SW_P, iZ_E0_SW_P, iZ_B1, iZ_E1, uGE, uGF, U_R )

#elif TWOMOMENT_ORDER_V
      CALL ApplySlopeLimiter_TwoMoment &
             ( iZ_B0_SW, iZ_E0_SW, iZ_B1, iZ_E1, uGE, uGF, U_F, U_R )

      CALL ApplyPositivityLimiter_TwoMoment &
             ( iZ_B0_SW_P, iZ_E0_SW_P, iZ_B1, iZ_E1, uGE, uGF, U_F, U_R )
#endif

      ! --- Apply Boundary Condition ---

      CALL ApplyBoundaryConditions_Radiation &
             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U_R, iZ_SW_P, bcX, iApplyBC )

#ifdef TWOMOMENT_ORDER_V
      CALL ApplyBoundaryConditions_Fluid &
             ( iX_B0, iX_E0, iX_B1, iX_E1, U_F, iZ_SW_P, bcX, iApplyBC )
#endif

      ! --- Explicit Solver ---

#ifdef TWOMOMENT_ORDER_1
      CALL ComputeIncrement_TwoMoment_Explicit &
             ( iZ_B0_SW, iZ_E0_SW, iZ_B1, iZ_E1, &
               uGE, uGF, U_R, T0_R )

#elif TWOMOMENT_ORDER_V
      CALL ComputeIncrement_TwoMoment_Explicit &
             ( iZ_B0_SW, iZ_E0_SW, iZ_B1, iZ_E1, &
               uGE, uGF, U_F, U_R, T0_R )
      OffGridFluxR_T0 = OffGridFlux_TwoMoment
#endif
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
      OffGridFluxR_T0 = Zero

    END IF

    ! --- Apply Increment ---

    CALL AddFields_Radiation &
           ( iZ_B0_SW, iZ_E0_SW, One, dt, U0_R, T0_R, U_R )

#ifdef TWOMOMENT_ORDER_V
    OffGridFluxR = dt * OffGridFluxR_T0
#endif

    ! --- Apply Limiter ---

#ifdef TWOMOMENT_ORDER_1
    CALL ApplyPositivityLimiter_TwoMoment &
           ( iZ_B0_SW, iZ_E0_SW, iZ_B1, iZ_E1, uGE, uGF, U_R )

#elif TWOMOMENT_ORDER_V
    CALL ApplySlopeLimiter_TwoMoment &
           ( iZ_B0_SW, iZ_E0_SW, iZ_B1, iZ_E1, uGE, uGF, U_F, U_R )

    CALL ApplyPositivityLimiter_TwoMoment &
           ( iZ_B0_SW, iZ_E0_SW, iZ_B1, iZ_E1, uGE, uGF, U_F, U_R )
#endif

    ! --- Implicit Step ---

    IF( Implicit )THEN
      CALL ComputeIncrement_TwoMoment_Implicit &
             ( iZ_B0_SW, iZ_E0_SW, iZ_B1, iZ_E1, dt, &
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
      !$OMP PARALLEL DO COLLAPSE(5)
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
           ( iX_B0_SW, iX_E0_SW, One, dt, U_F, Q1_F, U_F )

    CALL AddFields_Radiation &
           ( iZ_B0_SW, iZ_E0_SW, One, dt, U_R, Q1_R, U_R )

    ! --- Apply Positivity Limiter ---

#ifdef TWOMOMENT_ORDER_1
    CALL ApplyPositivityLimiter_TwoMoment &
           ( iZ_B0_SW, iZ_E0_SW, iZ_B1, iZ_E1, uGE, uGF, U_R )

#elif TWOMOMENT_ORDER_V
#if defined MICROPHYSICS_WEAKLIB
    CALL ApplyPositivityLimiter_Euler_NonRelativistic_TABLE &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, U_F, uDF )
#endif

    CALL ApplySlopeLimiter_TwoMoment &
           ( iZ_B0_SW, iZ_E0_SW, iZ_B1, iZ_E1, uGE, uGF, U_F, U_R )

    CALL ApplyPositivityLimiter_TwoMoment &
           ( iZ_B0_SW, iZ_E0_SW, iZ_B1, iZ_E1, uGE, uGF, U_F, U_R )

#endif

    IF( .NOT. SingleStage ) THEN

      ! ---------------
      ! --- Stage 2 ---
      ! ---------------

      iX_SW = [ 0, 0, 0 ]
      iZ_SW = [ 0, 0, 0, 0 ]
      if (nDimsX .eq. 3) then
         iZ_SW_P = [ 0, 1, 1, 1 ]
      else if (nDimsX .eq. 2) then
         iZ_SW_P = [ 0, 1, 1, 0 ]
      else if (nDimsX .eq. 1) then
         iZ_SW_P = [ 0, 1, 0, 0 ]
      end if

      iX_B0_SW = iX_B0 - iX_SW
      iX_E0_SW = iX_E0 + iX_SW
      iZ_B0_SW = iZ_B0 - iZ_SW
      iZ_E0_SW = iZ_E0 + iZ_SW

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET UPDATE &
      !$OMP TO( iX_B0_SW, iX_E0_SW, iZ_B0_SW, iZ_E0_SW, iZ_SW_P )
#elif defined(THORNADO_OACC)
      !$ACC UPDATE &
      !$ACC DEVICE( iX_B0_SW, iX_E0_SW, iZ_B0_SW, iZ_E0_SW, iZ_SW_P )
#endif

      ! --- Explicit Step (Radiation Only) ---

      IF( Explicit )THEN

        CALL ApplyBoundaryConditions_Radiation &
               ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U_R, iZ_SW_P, bcX, iApplyBC )

#ifdef TWOMOMENT_ORDER_V
        CALL ApplyBoundaryConditions_Fluid &
               ( iX_B0, iX_E0, iX_B1, iX_E1, U_F, iZ_SW_P, bcX, iApplyBC )
#endif

#ifdef TWOMOMENT_ORDER_1
        CALL ComputeIncrement_TwoMoment_Explicit &
               ( iZ_B0_SW, iZ_E0_SW, iZ_B1, iZ_E1, &
                 uGE, uGF, U_R, T1_R )

#elif TWOMOMENT_ORDER_V
        CALL ComputeIncrement_TwoMoment_Explicit &
               ( iZ_B0_SW, iZ_E0_SW, iZ_B1, iZ_E1, &
                 uGE, uGF, U_F, U_R, T1_R )
        OffGridFluxR_T1 = OffGridFlux_TwoMoment
#endif

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
        OffGridFluxR_T1 = Zero

      END IF

      ! --- Apply Increment ---

      CALL AddFields_Fluid &
             ( iX_B0_SW, iX_E0_SW, One, Half * dt, U0_F, Q1_F, U_F )

      CALL AddFields_Radiation &
             ( iZ_B0_SW, iZ_E0_SW, One, Half * dt, U0_R, T0_R, U_R )

      CALL AddFields_Radiation &
             ( iZ_B0_SW, iZ_E0_SW, One, Half * dt, U_R,  T1_R, U_R )

      CALL AddFields_Radiation &
             ( iZ_B0_SW, iZ_E0_SW, One, Half * dt, U_R,  Q1_R, U_R )

#ifdef TWOMOMENT_ORDER_V
      OffGridFluxR = Half * dt * OffGridFluxR_T0 + Half * dt * OffGridFluxR_T1
#endif

      ! --- Apply Limiter ---

#ifdef TWOMOMENT_ORDER_1
      CALL ApplyPositivityLimiter_TwoMoment &
             ( iZ_B0_SW, iZ_E0_SW, iZ_B1, iZ_E1, uGE, uGF, U_R )

#elif TWOMOMENT_ORDER_V
      CALL ApplySlopeLimiter_TwoMoment &
             ( iZ_B0_SW, iZ_E0_SW, iZ_B1, iZ_E1, uGE, uGF, U_F, U_R )

      CALL ApplyPositivityLimiter_TwoMoment &
             ( iZ_B0_SW, iZ_E0_SW, iZ_B1, iZ_E1, uGE, uGF, U_F, U_R )

#endif

      ! --- Implicit Step ---

      IF( Implicit )THEN

        CALL ComputeIncrement_TwoMoment_Implicit &
               (iZ_B0_SW, iZ_E0_SW, iZ_B1, iZ_E1, Half * dt, &
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
        !$OMP PARALLEL DO COLLAPSE(5)
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
             ( iX_B0_SW, iX_E0_SW, One, Half * dt, U_F, Q1_F, U_F )

      CALL AddFields_Radiation &
             ( iZ_B0_SW, iZ_E0_SW, One, Half * dt, U_R, Q1_R, U_R )

      ! --- Apply Limiter ---
#ifdef TWOMOMENT_ORDER_1
      CALL ApplyPositivityLimiter_TwoMoment &
             ( iZ_B0_SW, iZ_E0_SW, iZ_B1, iZ_E1, uGE, uGF, U_R )

#elif TWOMOMENT_ORDER_V
#if defined MICROPHYSICS_WEAKLIB
      CALL ApplyPositivityLimiter_Euler_NonRelativistic_TABLE &
             ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, U_F, uDF )
#endif

      CALL ApplySlopeLimiter_TwoMoment &
             ( iZ_B0_SW, iZ_E0_SW, iZ_B1, iZ_E1, uGE, uGF, U_F, U_R )

      CALL ApplyPositivityLimiter_TwoMoment &
             ( iZ_B0_SW, iZ_E0_SW, iZ_B1, iZ_E1, uGE, uGF, U_F, U_R )
#endif

    END IF

#ifdef THORNADO_DEBUG
    IF( ANY(IEEE_IS_NAN(U_R)) ) PRINT*, 'NaN @ end of [Update_IMEX_PDARS]'
#endif

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
    !$OMP MAP( release: U0_F, Q1_F, U0_R, T0_R, T1_R, Q1_R, uGE, uGF, &
    !$OMP               iX_B0_SW, iX_E0_SW, iZ_B0_SW, iZ_E0_SW, iZ_B0_SW_P, iZ_E0_SW_P, iZ_SW_P )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC COPYOUT( U_F, U_R ) &
    !$ACC DELETE( U0_F, Q1_F, U0_R, T0_R, T1_R, Q1_R, uGE, uGF, &
    !$ACC         iX_B0_SW, iX_E0_SW, iZ_B0_SW, iZ_E0_SW, iZ_B0_SW_P, iZ_E0_SW_P, iZ_SW_P )
#endif

    IF( PRESENT( OffGridFluxR_Option ) )THEN
      OffGridFluxR_Option = OffGridFluxR
    END IF

  END SUBROUTINE Update_IMEX_PDARS


  SUBROUTINE AddFields_Fluid( iX_B, iX_E, alpha, beta, A, B, C )

    ! --- C = alpha * A + beta * B

    INTEGER,  INTENT(in)    :: iX_B(3), iX_E(3)
    REAL(DP), INTENT(in)    :: alpha, beta
    REAL(DP), INTENT(inout) :: &
      A(1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nCF)
    REAL(DP), INTENT(inout) :: &
      B(1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nCF)
    REAL(DP), INTENT(inout) :: &
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
    REAL(DP), INTENT(inout) :: &
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


  LOGICAL FUNCTION ApplyInnerBC( iApplyBC )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    INTEGER, INTENT(in) :: iApplyBC

    ApplyInnerBC = .FALSE.
    IF( iApplyBC .EQ. iApplyBC_Inner .OR. &
        iApplyBC .EQ. iApplyBC_Both ) &
    ApplyInnerBC = .TRUE.

  END FUNCTION ApplyInnerBC


  LOGICAL FUNCTION ApplyOuterBC( iApplyBC )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    INTEGER, INTENT(in) :: iApplyBC

    ApplyOuterBC = .FALSE.
    IF( iApplyBC .EQ. iApplyBC_Outer .OR. &
        iApplyBC .EQ. iApplyBC_Both ) &
    ApplyOuterBC = .TRUE.

  END FUNCTION ApplyOuterBC


  SUBROUTINE ApplyBoundaryConditions_Radiation &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U, swZ, bcX, iApplyBC )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4), swZ(4), bcX(3), iApplyBC(3)
    REAL(DP), INTENT(inout) :: &
      U(1:nDOF, &
        iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
        iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4), &
        1:nCR,1:nSpecies)

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: U, iZ_B0, iZ_E0, iZ_B1, iZ_E1, swZ, bcX, iApplyBC )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA ASYNC &
    !$ACC COPYIN( U, iZ_B0, iZ_E0, iZ_B1, iZ_E1, swZ, bcX, iApplyBC )
#endif

    CALL ApplyBC_Radiation_X1 &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U, swZ(2), bcX(1), iApplyBC(1) )

    CALL ApplyBC_Radiation_X2 &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U, swZ(3), bcX(2), iApplyBC(2) )

    CALL ApplyBC_Radiation_X3 &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U, swZ(4), bcX(3), iApplyBC(3) )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from: U ) &
    !$OMP MAP( release: iZ_B0, iZ_E0, iZ_B1, iZ_E1, swZ, bcX, iApplyBC )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA ASYNC &
    !$ACC COPYOUT( U ) &
    !$ACC DELETE( iZ_B0, iZ_E0, iZ_B1, iZ_E1, swZ, bcX, iApplyBC )
#endif

  END SUBROUTINE ApplyBoundaryConditions_Radiation


  SUBROUTINE ApplyBoundaryConditions_Fluid &
    ( iX_B0, iX_E0, iX_B1, iX_E1, U, swZ, bcX, iApplyBC )

    INTEGER,  INTENT(in)           :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), swZ(4), bcX(3), iApplyBC(3)
    REAL(DP), INTENT(inout)        :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: U, iX_B0, iX_E0, iX_B1, iX_E1, swZ, bcX, iApplyBC )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN(  U, iX_B0, iX_E0, iX_B1, iX_E1, swZ, bcX, iApplyBC )
#endif

    CALL ApplyBC_Fluid_X1( iX_B0, iX_E0, iX_B1, iX_E1, U, swZ(2), bcX(1), iApplyBC(1) )

    CALL ApplyBC_Fluid_X2( iX_B0, iX_E0, iX_B1, iX_E1, U, swZ(3), bcX(2), iApplyBC(2) )

    CALL ApplyBC_Fluid_X3( iX_B0, iX_E0, iX_B1, iX_E1, U, swZ(4), bcX(3), iApplyBC(3) )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from:    U ) &
    !$OMP MAP( release: iX_B0, iX_E0, iX_B1, iX_E1, swZ, bcX, iApplyBC )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC COPYOUT(      U ) &
    !$ACC DELETE(       iX_B0, iX_E0, iX_B1, iX_E1, swZ, bcX, iApplyBC )
#endif

  END SUBROUTINE ApplyBoundaryConditions_Fluid


  SUBROUTINE ApplyBC_Fluid_X1( iX_B0, iX_E0, iX_B1, iX_E1, U, swX1, bcX1, iApplyBC )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), &
      swX1, bcX1, iApplyBC
    REAL(DP), INTENT(inout) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER  :: iCF, iX1, iX2, iX3
    INTEGER  :: iNX
    INTEGER  :: iNX1, iNX2, iNX3, jNX, jNX1

    SELECT CASE ( bcX1 )

    CASE ( 0 ) ! No Boundary Condition

    CASE ( 2 ) ! Homogeneous

      ! --- Inner Boundary ---

      IF( ApplyInnerBC( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
        !$ACC PRESENT( U, iX_B0, iX_E0 )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(5)
#endif
        DO iCF = 1, nCF
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = 1, swX1
        DO iNX = 1, nDOFX

          IF( iCF == iCF_S1 )THEN

            U(iNX,iX_B0(1)-iX1,iX2,iX3,iCF) = MIN( U(iNX,iX_B0(1),iX2,iX3,iCF), Zero )

          ELSE

            U(iNX,iX_B0(1)-iX1,iX2,iX3,iCF) = U(iNX,iX_B0(1),iX2,iX3,iCF)

          END IF

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

      ! --- Outer Boundary ---

      IF( ApplyOuterBC( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
        !$ACC PRESENT( U, iX_B0, iX_E0 )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(5)
#endif
        DO iCF = 1, nCF
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = 1, swX1
        DO iNX = 1, nDOFX

          IF( iCF == iCF_S1 )THEN

            U(iNX,iX_E0(1)+iX1,iX2,iX3,iCF) = MAX( U(iNX,iX_E0(1),iX2,iX3,iCF), Zero )

          ELSE

            U(iNX,iX_E0(1)+iX1,iX2,iX3,iCF) = U(iNX,iX_E0(1),iX2,iX3,iCF)

          END IF

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

    CASE ( 3 ) ! Reflecting

      ! --- Inner Boundary ---

      IF( ApplyInnerBC( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
        !$OMP PRIVATE( iNX, jNX, jNX1 )
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) &
        !$ACC PRIVATE( iNX, jNX, jNX1 ) &
        !$ACC PRESENT( U, iX_B0, iX_E0, nNodesX )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(6) &
        !$OMP PRIVATE( iNX, jNX, jNX1 )
#endif
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = 1, swX1

          DO iNX3 = 1, nNodesX(3)
          DO iNX2 = 1, nNodesX(2)
          DO iNX1 = 1, nNodesX(1)

            jNX1 = ( nNodesX(1) - iNX1 ) + 1

            iNX = NodeNumberX( iNX1, iNX2, iNX3 )
            jNX = NodeNumberX( jNX1, iNX2, iNX3 )

            U(iNX,iX_B0(1)-iX1,iX2,iX3,iCF_D ) = + U(jNX,iX_B0(1)+iX1-1,iX2,iX3,iCF_D )
            U(iNX,iX_B0(1)-iX1,iX2,iX3,iCF_S1) = - U(jNX,iX_B0(1)+iX1-1,iX2,iX3,iCF_S1)
            U(iNX,iX_B0(1)-iX1,iX2,iX3,iCF_S2) = + U(jNX,iX_B0(1)+iX1-1,iX2,iX3,iCF_S2)
            U(iNX,iX_B0(1)-iX1,iX2,iX3,iCF_S3) = + U(jNX,iX_B0(1)+iX1-1,iX2,iX3,iCF_S3)
            U(iNX,iX_B0(1)-iX1,iX2,iX3,iCF_E ) = + U(jNX,iX_B0(1)+iX1-1,iX2,iX3,iCF_E )
            U(iNX,iX_B0(1)-iX1,iX2,iX3,iCF_Ne) = + U(jNX,iX_B0(1)+iX1-1,iX2,iX3,iCF_Ne)

          END DO
          END DO
          END DO

        END DO
        END DO
        END DO

      END IF

      ! --- Outer Boundary ---

      IF( ApplyOuterBC( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
        !$OMP PRIVATE( iNX, jNX, jNX1 )
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) &
        !$ACC PRIVATE( iNX, jNX, jNX1 ) &
        !$ACC PRESENT( U, iX_B0, iX_E0, nNodesX )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(6) &
        !$OMP PRIVATE( iNX, jNX, jNX1 )
#endif
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = 1, swX1

          DO iNX3 = 1, nNodesX(3)
          DO iNX2 = 1, nNodesX(2)
          DO iNX1 = 1, nNodesX(1)

            jNX1 = ( nNodesX(1) - iNX1 ) + 1

            iNX = NodeNumberX( iNX1, iNX2, iNX3 )
            jNX = NodeNumberX( jNX1, iNX2, iNX3 )

            U(iNX,iX_E0(1)+iX1,iX2,iX3,iCF_D ) = + U(jNX,iX_E0(1)-iX1+1,iX2,iX3,iCF_D )
            U(iNX,iX_E0(1)+iX1,iX2,iX3,iCF_S1) = - U(jNX,iX_E0(1)-iX1+1,iX2,iX3,iCF_S1)
            U(iNX,iX_E0(1)+iX1,iX2,iX3,iCF_S2) = + U(jNX,iX_E0(1)-iX1+1,iX2,iX3,iCF_S2)
            U(iNX,iX_E0(1)+iX1,iX2,iX3,iCF_S3) = + U(jNX,iX_E0(1)-iX1+1,iX2,iX3,iCF_S3)
            U(iNX,iX_E0(1)+iX1,iX2,iX3,iCF_E ) = + U(jNX,iX_E0(1)-iX1+1,iX2,iX3,iCF_E )
            U(iNX,iX_E0(1)+iX1,iX2,iX3,iCF_Ne) = + U(jNX,iX_E0(1)-iX1+1,iX2,iX3,iCF_Ne)

          END DO
          END DO
          END DO

        END DO
        END DO
        END DO

      END IF

    CASE ( 32 ) ! Reflecting (Inner), Homogeneous (Outer)

      ! --- Inner Boundary ---

      IF( ApplyInnerBC( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
        !$OMP PRIVATE( iNX, jNX, jNX1 )
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) &
        !$ACC PRIVATE( iNX, jNX, jNX1 ) &
        !$ACC PRESENT( U, iX_B0, iX_E0, nNodesX )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(6) &
        !$OMP PRIVATE( iNX, jNX, jNX1 )
#endif
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = 1, swX1

          DO iNX3 = 1, nNodesX(3)
          DO iNX2 = 1, nNodesX(2)
          DO iNX1 = 1, nNodesX(1)

            jNX1 = ( nNodesX(1) - iNX1 ) + 1

            iNX = NodeNumberX( iNX1, iNX2, iNX3 )
            jNX = NodeNumberX( jNX1, iNX2, iNX3 )

            U(iNX,iX_B0(1)-iX1,iX2,iX3,iCF_D ) = + U(jNX,iX_B0(1)+iX1-1,iX2,iX3,iCF_D )
            U(iNX,iX_B0(1)-iX1,iX2,iX3,iCF_S1) = - U(jNX,iX_B0(1)+iX1-1,iX2,iX3,iCF_S1)
            U(iNX,iX_B0(1)-iX1,iX2,iX3,iCF_S2) = + U(jNX,iX_B0(1)+iX1-1,iX2,iX3,iCF_S2)
            U(iNX,iX_B0(1)-iX1,iX2,iX3,iCF_S3) = + U(jNX,iX_B0(1)+iX1-1,iX2,iX3,iCF_S3)
            U(iNX,iX_B0(1)-iX1,iX2,iX3,iCF_E ) = + U(jNX,iX_B0(1)+iX1-1,iX2,iX3,iCF_E )
            U(iNX,iX_B0(1)-iX1,iX2,iX3,iCF_Ne) = + U(jNX,iX_B0(1)+iX1-1,iX2,iX3,iCF_Ne)

          END DO
          END DO
          END DO

        END DO
        END DO
        END DO

      END IF

      ! --- Outer Boundary ---

      IF( ApplyOuterBC( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
        !$ACC PRESENT( U, iX_B0, iX_E0 )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(5)
#endif
        DO iCF = 1, nCF
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = 1, swX1
        DO iNX = 1, nDOFX

          IF( iCF == iCF_S1 )THEN

            U(iNX,iX_E0(1)+iX1,iX2,iX3,iCF) = MAX( U(iNX,iX_E0(1),iX2,iX3,iCF), Zero )

          ELSE

            U(iNX,iX_E0(1)+iX1,iX2,iX3,iCF) = U(iNX,iX_E0(1),iX2,iX3,iCF)

          END IF

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

    CASE DEFAULT

    END SELECT

  END SUBROUTINE ApplyBC_Fluid_X1


  SUBROUTINE ApplyBC_Fluid_X2( iX_B0, iX_E0, iX_B1, iX_E1, U, swX2, bcX2, iApplyBC )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), &
      swX2, bcX2, iApplyBC
    REAL(DP), INTENT(inout) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER  :: iCF, iX1, iX2, iX3
    INTEGER  :: iNX
    INTEGER  :: iNX1, iNX2, iNX3, jNX, jNX2

    SELECT CASE ( bcX2 )

    CASE ( 0 ) ! No Boundary Condition

    CASE ( 2 ) ! Homogeneous

      ! --- Inner Boundary ---

      IF( ApplyInnerBC( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
        !$ACC PRESENT( U, iX_B0, iX_E0 )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(5)
#endif
        DO iCF = 1, nCF
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = 1, swX2
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1, nDOFX

          IF( iCF == iCF_S2 )THEN

            U(iNX,iX1,iX_B0(2)-iX2,iX3,iCF) = MIN( U(iNX,iX1,iX_B0(2),iX3,iCF), Zero )

          ELSE

            U(iNX,iX1,iX_B0(2)-iX2,iX3,iCF) = U(iNX,iX1,iX_B0(2),iX3,iCF)

          END IF

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

      ! --- Outer Boundary ---

      IF( ApplyOuterBC( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
        !$ACC PRESENT( U, iX_B0, iX_E0 )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(5)
#endif
        DO iCF = 1, nCF
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = 1, swX2
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1, nDOFX

          IF( iCF == iCF_S2 )THEN

            U(iNX,iX1,iX_E0(2)+iX2,iX3,iCF) = MAX( U(iNX,iX1,iX_E0(2),iX3,iCF), Zero )

          ELSE

            U(iNX,iX1,iX_E0(2)+iX2,iX3,iCF) = U(iNX,iX1,iX_E0(2),iX3,iCF)

          END IF

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

    CASE ( 3 ) ! Reflecting

      ! --- Inner Boundary ---

      IF( ApplyInnerBC( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
        !$OMP PRIVATE( iNX, jNX, jNX2 )
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) &
        !$ACC PRIVATE( iNX, jNX, jNX2 ) &
        !$ACC PRESENT( U, iX_B0, iX_E0, nNodesX )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(6) &
        !$OMP PRIVATE( iNX, jNX, jNX2 )
#endif
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = 1, swX2
        DO iX1 = iX_B0(1), iX_E0(1)

          DO iNX3 = 1, nNodesX(3)
          DO iNX2 = 1, nNodesX(2)
          DO iNX1 = 1, nNodesX(1)

            jNX2 = ( nNodesX(2) - iNX2 ) + 1

            iNX = NodeNumberX( iNX1, iNX2, iNX3 )
            jNX = NodeNumberX( iNX1, jNX2, iNX3 )

            U(iNX,iX1,iX_B0(2)-iX2,iX3,iCF_D ) = + U(jNX,iX1,iX_B0(2)+iX2-1,iX3,iCF_D )
            U(iNX,iX1,iX_B0(2)-iX2,iX3,iCF_S1) = + U(jNX,iX1,iX_B0(2)+iX2-1,iX3,iCF_S1)
            U(iNX,iX1,iX_B0(2)-iX2,iX3,iCF_S2) = - U(jNX,iX1,iX_B0(2)+iX2-1,iX3,iCF_S2)
            U(iNX,iX1,iX_B0(2)-iX2,iX3,iCF_S3) = + U(jNX,iX1,iX_B0(2)+iX2-1,iX3,iCF_S3)
            U(iNX,iX1,iX_B0(2)-iX2,iX3,iCF_E ) = + U(jNX,iX1,iX_B0(2)+iX2-1,iX3,iCF_E )
            U(iNX,iX1,iX_B0(2)-iX2,iX3,iCF_Ne) = + U(jNX,iX1,iX_B0(2)+iX2-1,iX3,iCF_Ne)

          END DO
          END DO
          END DO

        END DO
        END DO
        END DO

      END IF

      ! --- Outer Boundary ---

      IF( ApplyOuterBC( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
        !$OMP PRIVATE( iNX, jNX, jNX2 )
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) &
        !$ACC PRIVATE( iNX, jNX, jNX2 ) &
        !$ACC PRESENT( U, iX_B0, iX_E0, nNodesX )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(6) &
        !$OMP PRIVATE( iNX, jNX, jNX2 )
#endif
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = 1, swX2
        DO iX1 = iX_B0(1), iX_E0(1)

          DO iNX3 = 1, nNodesX(3)
          DO iNX2 = 1, nNodesX(2)
          DO iNX1 = 1, nNodesX(1)

            jNX2 = ( nNodesX(2) - iNX2 ) + 1

            iNX = NodeNumberX( iNX1, iNX2, iNX3 )
            jNX = NodeNumberX( iNX1, jNX2, iNX3 )

            U(iNX,iX1,iX_E0(2)+iX2,iX3,iCF_D ) = + U(jNX,iX1,iX_E0(2)-iX2+1,iX3,iCF_D )
            U(iNX,iX1,iX_E0(2)+iX2,iX3,iCF_S1) = + U(jNX,iX1,iX_E0(2)-iX2+1,iX3,iCF_S1)
            U(iNX,iX1,iX_E0(2)+iX2,iX3,iCF_S2) = - U(jNX,iX1,iX_E0(2)-iX2+1,iX3,iCF_S2)
            U(iNX,iX1,iX_E0(2)+iX2,iX3,iCF_S3) = + U(jNX,iX1,iX_E0(2)-iX2+1,iX3,iCF_S3)
            U(iNX,iX1,iX_E0(2)+iX2,iX3,iCF_E ) = + U(jNX,iX1,iX_E0(2)-iX2+1,iX3,iCF_E )
            U(iNX,iX1,iX_E0(2)+iX2,iX3,iCF_Ne) = + U(jNX,iX1,iX_E0(2)-iX2+1,iX3,iCF_Ne)

          END DO
          END DO
          END DO

        END DO
        END DO
        END DO

      END IF

    CASE ( 32 ) ! Reflecting (Inner), Homogeneous (Outer)

      ! --- Inner Boundary ---

      IF( ApplyInnerBC( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
        !$OMP PRIVATE( iNX, jNX, jNX2 )
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) &
        !$ACC PRIVATE( iNX, jNX, jNX2 ) &
        !$ACC PRESENT( U, iX_B0, iX_E0, nNodesX )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(6) &
        !$OMP PRIVATE( iNX, jNX, jNX2 )
#endif
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = 1, swX2
        DO iX1 = iX_B0(1), iX_E0(1)

          DO iNX3 = 1, nNodesX(3)
          DO iNX2 = 1, nNodesX(2)
          DO iNX1 = 1, nNodesX(1)

            jNX2 = ( nNodesX(2) - iNX2 ) + 1

            iNX = NodeNumberX( iNX1, iNX2, iNX3 )
            jNX = NodeNumberX( iNX1, jNX2, iNX3 )

            U(iNX,iX1,iX_B0(2)-iX2,iX3,iCF_D ) = + U(jNX,iX1,iX_B0(2)+iX2-1,iX3,iCF_D )
            U(iNX,iX1,iX_B0(2)-iX2,iX3,iCF_S1) = + U(jNX,iX1,iX_B0(2)+iX2-1,iX3,iCF_S1)
            U(iNX,iX1,iX_B0(2)-iX2,iX3,iCF_S2) = - U(jNX,iX1,iX_B0(2)+iX2-1,iX3,iCF_S2)
            U(iNX,iX1,iX_B0(2)-iX2,iX3,iCF_S3) = + U(jNX,iX1,iX_B0(2)+iX2-1,iX3,iCF_S3)
            U(iNX,iX1,iX_B0(2)-iX2,iX3,iCF_E ) = + U(jNX,iX1,iX_B0(2)+iX2-1,iX3,iCF_E )
            U(iNX,iX1,iX_B0(2)-iX2,iX3,iCF_Ne) = + U(jNX,iX1,iX_B0(2)+iX2-1,iX3,iCF_Ne)

          END DO
          END DO
          END DO

        END DO
        END DO
        END DO

      END IF

      ! --- Outer Boundary ---

      IF( ApplyOuterBC( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
        !$ACC PRESENT( U, iX_B0, iX_E0 )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(5)
#endif
        DO iCF = 1, nCF
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = 1, swX2
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1, nDOFX

          IF( iCF == iCF_S2 )THEN

            U(iNX,iX1,iX_E0(2)+iX2,iX3,iCF) = MAX( U(iNX,iX1,iX_E0(2),iX3,iCF), Zero )

          ELSE

            U(iNX,iX1,iX_E0(2)+iX2,iX3,iCF) = U(iNX,iX1,iX_E0(2),iX3,iCF)

          END IF

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

    CASE DEFAULT

    END SELECT

  END SUBROUTINE ApplyBC_Fluid_X2


  SUBROUTINE ApplyBC_Fluid_X3( iX_B0, iX_E0, iX_B1, iX_E1, U, swX3, bcX3, iApplyBC )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), &
      swX3, bcX3, iApplyBC
    REAL(DP), INTENT(inout) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER  :: iCF, iX1, iX2, iX3
    INTEGER  :: iNX
    INTEGER  :: iNX1, iNX2, iNX3, jNX, jNX3

    SELECT CASE ( bcX3 )

    CASE ( 0 ) ! No Boundary Condition

    CASE ( 2 ) ! Homogeneous

      ! --- Inner Boundary ---

      IF( ApplyInnerBC( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
        !$ACC PRESENT( U, iX_B0, iX_E0 )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(5)
#endif
        DO iCF = 1, nCF
        DO iX3 = 1, swX3
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1, nDOFX

          IF( iCF == iCF_S3 )THEN

            U(iNX,iX1,iX2,iX_B0(3)-iX3,iCF) = MIN( U(iNX,iX1,iX2,iX_B0(3),iCF), Zero )

          ELSE

            U(iNX,iX1,iX2,iX_B0(3)-iX3,iCF) = U(iNX,iX1,iX2,iX_B0(3),iCF)

          END IF

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

      ! --- Outer Boundary ---

      IF( ApplyOuterBC( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
        !$ACC PRESENT( U, iX_B0, iX_E0 )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(5)
#endif
        DO iCF = 1, nCF
        DO iX3 = 1, swX3
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1, nDOFX

          IF( iCF == iCF_S3 )THEN

            U(iNX,iX1,iX2,iX_E0(3)+iX3,iCF) = MAX( U(iNX,iX1,iX2,iX_E0(3),iCF), Zero )

          ELSE

            U(iNX,iX1,iX2,iX_E0(3)+iX3,iCF) = U(iNX,iX1,iX2,iX_E0(3),iCF)

          END IF

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

    CASE ( 3 ) ! Reflecting

      ! --- Inner Boundary ---

      IF( ApplyInnerBC( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
        !$OMP PRIVATE( iNX, jNX, jNX3 )
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) &
        !$ACC PRIVATE( iNX, jNX, jNX3 ) &
        !$ACC PRESENT( U, iX_B0, iX_E0, nNodesX )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(6) &
        !$OMP PRIVATE( iNX, jNX, jNX3 )
#endif
        DO iX3 = 1, swX3
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)

          DO iNX3 = 1, nNodesX(3)
          DO iNX2 = 1, nNodesX(2)
          DO iNX1 = 1, nNodesX(1)

            jNX3 = ( nNodesX(3) - iNX3 ) + 1

            iNX = NodeNumberX( iNX1, iNX2, iNX3 )
            jNX = NodeNumberX( iNX1, iNX2, jNX3 )

            U(iNX,iX1,iX2,iX_B0(3)-iX3,iCF_D ) = + U(jNX,iX1,iX2,iX_B0(3)+iX3-1,iCF_D )
            U(iNX,iX1,iX2,iX_B0(3)-iX3,iCF_S1) = + U(jNX,iX1,iX2,iX_B0(3)+iX3-1,iCF_S1)
            U(iNX,iX1,iX2,iX_B0(3)-iX3,iCF_S2) = + U(jNX,iX1,iX2,iX_B0(3)+iX3-1,iCF_S2)
            U(iNX,iX1,iX2,iX_B0(3)-iX3,iCF_S3) = - U(jNX,iX1,iX2,iX_B0(3)+iX3-1,iCF_S3)
            U(iNX,iX1,iX2,iX_B0(3)-iX3,iCF_E ) = + U(jNX,iX1,iX2,iX_B0(3)+iX3-1,iCF_E )
            U(iNX,iX1,iX2,iX_B0(3)-iX3,iCF_Ne) = + U(jNX,iX1,iX2,iX_B0(3)+iX3-1,iCF_Ne)

          END DO
          END DO
          END DO

        END DO
        END DO
        END DO

      END IF

      ! --- Outer Boundary ---

      IF( ApplyOuterBC( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
        !$OMP PRIVATE( iNX, jNX, jNX3 )
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) &
        !$ACC PRIVATE( iNX, jNX, jNX3 ) &
        !$ACC PRESENT( U, iX_B0, iX_E0, nNodesX )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(6) &
        !$OMP PRIVATE( iNX, jNX, jNX3 )
#endif
        DO iX3 = 1, swX3
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)

          DO iNX3 = 1, nNodesX(3)
          DO iNX2 = 1, nNodesX(2)
          DO iNX1 = 1, nNodesX(1)

            jNX3 = ( nNodesX(3) - iNX3 ) + 1

            iNX = NodeNumberX( iNX1, iNX2, iNX3 )
            jNX = NodeNumberX( iNX1, iNX2, jNX3 )

            U(iNX,iX1,iX2,iX_E0(3)+iX3,iCF_D ) = + U(jNX,iX1,iX2,iX_E0(3)-iX3+1,iCF_D )
            U(iNX,iX1,iX2,iX_E0(3)+iX3,iCF_S1) = + U(jNX,iX1,iX2,iX_E0(3)-iX3+1,iCF_S1)
            U(iNX,iX1,iX2,iX_E0(3)+iX3,iCF_S2) = + U(jNX,iX1,iX2,iX_E0(3)-iX3+1,iCF_S2)
            U(iNX,iX1,iX2,iX_E0(3)+iX3,iCF_S3) = - U(jNX,iX1,iX2,iX_E0(3)-iX3+1,iCF_S3)
            U(iNX,iX1,iX2,iX_E0(3)+iX3,iCF_E ) = + U(jNX,iX1,iX2,iX_E0(3)-iX3+1,iCF_E )
            U(iNX,iX1,iX2,iX_E0(3)+iX3,iCF_Ne) = + U(jNX,iX1,iX2,iX_E0(3)-iX3+1,iCF_Ne)

          END DO
          END DO
          END DO

        END DO
        END DO
        END DO

      END IF

    CASE ( 32 ) ! Reflecting (Inner), Homogeneous (Outer)

      ! --- Inner Boundary ---

      IF( ApplyInnerBC( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
        !$OMP PRIVATE( iNX, jNX, jNX3 )
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) &
        !$ACC PRIVATE( iNX, jNX, jNX3 ) &
        !$ACC PRESENT( U, iX_B0, iX_E0, nNodesX )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(6) &
        !$OMP PRIVATE( iNX, jNX, jNX3 )
#endif
        DO iX3 = 1, swX3
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)

          DO iNX3 = 1, nNodesX(3)
          DO iNX2 = 1, nNodesX(2)
          DO iNX1 = 1, nNodesX(1)

            jNX3 = ( nNodesX(3) - iNX3 ) + 1

            iNX = NodeNumberX( iNX1, iNX2, iNX3 )
            jNX = NodeNumberX( iNX1, iNX2, jNX3 )

            U(iNX,iX1,iX2,iX_B0(3)-iX3,iCF_D ) = + U(jNX,iX1,iX2,iX_B0(3)+iX3-1,iCF_D )
            U(iNX,iX1,iX2,iX_B0(3)-iX3,iCF_S1) = + U(jNX,iX1,iX2,iX_B0(3)+iX3-1,iCF_S1)
            U(iNX,iX1,iX2,iX_B0(3)-iX3,iCF_S2) = + U(jNX,iX1,iX2,iX_B0(3)+iX3-1,iCF_S2)
            U(iNX,iX1,iX2,iX_B0(3)-iX3,iCF_S3) = - U(jNX,iX1,iX2,iX_B0(3)+iX3-1,iCF_S3)
            U(iNX,iX1,iX2,iX_B0(3)-iX3,iCF_E ) = + U(jNX,iX1,iX2,iX_B0(3)+iX3-1,iCF_E )
            U(iNX,iX1,iX2,iX_B0(3)-iX3,iCF_Ne) = + U(jNX,iX1,iX2,iX_B0(3)+iX3-1,iCF_Ne)

          END DO
          END DO
          END DO

        END DO
        END DO
        END DO

      END IF

      ! --- Outer Boundary ---

      IF( ApplyOuterBC( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
        !$ACC PRESENT( U, iX_B0, iX_E0 )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(5)
#endif
        DO iCF = 1, nCF
        DO iX3 = 1, swX3
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1, nDOFX

          IF( iCF == iCF_S3 )THEN

            U(iNX,iX1,iX2,iX_E0(3)+iX3,iCF) = MAX( U(iNX,iX1,iX2,iX_E0(3),iCF), Zero )

          ELSE

            U(iNX,iX1,iX2,iX_E0(3)+iX3,iCF) = U(iNX,iX1,iX2,iX_E0(3),iCF)

          END IF

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

    CASE DEFAULT

    END SELECT

  END SUBROUTINE ApplyBC_Fluid_X3


  SUBROUTINE ApplyBC_Radiation_X1( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U, swX1, bcX1, iApplyBC )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4), &
      swX1, bcX1, iApplyBC
    REAL(DP), INTENT(inout) :: &
      U(1:nDOF, &
        iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
        iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4), &
        1:nCR,1:nSpecies)

    INTEGER :: iNode, iS, iCR, iZ1, iZ2, iZ3, iZ4
    INTEGER :: iNodeZ1, iNodeZ2, iNodeZ3, iNodeZ4
    INTEGER :: jNodeZ2, iNodeZ, jNodeZ

    SELECT CASE ( bcX1 )

    CASE ( 0 ) ! No Boundary Condition

    CASE ( 2 ) ! Homogeneous

      IF( ApplyInnerBC( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) ASYNC &
        !$ACC PRESENT( U, iZ_B0, iZ_E0 )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(7)
#endif
        DO iS = 1, nSpecies
          DO iCR = 1, nCR
            DO iZ4 = iZ_B0(4), iZ_E0(4)
              DO iZ3 = iZ_B0(3), iZ_E0(3)
                DO iZ2 = 1, swX1
                  DO iZ1 = iZ_B0(1), iZ_E0(1)
                    DO iNode = 1, nDOF

                      ! --- Inner Boundary ---

                      IF( iCR == iCR_G1 )THEN

                        U(iNode,iZ1,iZ_B0(2)-iZ2,iZ3,iZ4,iCR,iS) &
                          = MIN( U(iNode,iZ1,iZ_B0(2),iZ3,iZ4,iCR,iS), Zero )

                      ELSE

                        U(iNode,iZ1,iZ_B0(2)-iZ2,iZ3,iZ4,iCR,iS) &
                          = U(iNode,iZ1,iZ_B0(2),iZ3,iZ4,iCR,iS)

                      END IF

                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO

      END IF

      IF( ApplyOuterBC( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) ASYNC &
        !$ACC PRESENT( U, iZ_B0, iZ_E0 )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(7)
#endif
        DO iS = 1, nSpecies
          DO iCR = 1, nCR
            DO iZ4 = iZ_B0(4), iZ_E0(4)
              DO iZ3 = iZ_B0(3), iZ_E0(3)
                DO iZ2 = 1, swX1
                  DO iZ1 = iZ_B0(1), iZ_E0(1)
                    DO iNode = 1, nDOF

                      ! --- Outer Boundary ---

                      IF( iCR == iCR_G1 )THEN

                        U(iNode,iZ1,iZ_E0(2)+iZ2,iZ3,iZ4,iCR,iS) &
                          = MAX( U(iNode,iZ1,iZ_E0(2),iZ3,iZ4,iCR,iS), Zero )

                      ELSE

                        U(iNode,iZ1,iZ_E0(2)+iZ2,iZ3,iZ4,iCR,iS) &
                          = U(iNode,iZ1,iZ_E0(2),iZ3,iZ4,iCR,iS)

                      END IF

                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO

      END IF

    CASE ( 3 ) ! Reflecting

      IF( ApplyInnerBC( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(10) &
        !$OMP PRIVATE( jNodeZ2, iNodeZ, jNodeZ )
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(10) ASYNC &
        !$ACC PRIVATE( jNodeZ2, iNodeZ, jNodeZ ) &
        !$ACC PRESENT( U, iZ_B0, iZ_E0, nNodesZ, NodeNumberTable4D )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(10) &
        !$OMP PRIVATE( jNodeZ2, iNodeZ, jNodeZ )
#endif
        DO iS = 1, nSpecies
          DO iCR = 1, nCR
            DO iZ4 = iZ_B0(4), iZ_E0(4)
              DO iZ3 = iZ_B0(3), iZ_E0(3)
                DO iZ2 = 1, swX1
                  DO iZ1 = iZ_B0(1), iZ_E0(1)

                    DO iNodeZ4 = 1, nNodesZ(4)
                    DO iNodeZ3 = 1, nNodesZ(3)
                    DO iNodeZ2 = 1, nNodesZ(2)
                    DO iNodeZ1 = 1, nNodesZ(1)

                      jNodeZ2 = (nNodesZ(2)-iNodeZ2) + 1

                      iNodeZ = NodeNumberTable4D &
                                 ( iNodeZ1, iNodeZ2, iNodeZ3, iNodeZ4 )
                      jNodeZ = NodeNumberTable4D &
                                 ( iNodeZ1, jNodeZ2, iNodeZ3, iNodeZ4 )

                      ! --- Inner Boundary ---

                      IF( iCR == iCR_G1 )THEN

                        U(iNodeZ,iZ1,iZ_B0(2)-iZ2,iZ3,iZ4,iCR,iS) &
                          = - U(jNodeZ,iZ1,iZ_B0(2)+iZ2-1,iZ3,iZ4,iCR,iS)

                      ELSE

                        U(iNodeZ,iZ1,iZ_B0(2)-iZ2,iZ3,iZ4,iCR,iS) &
                          = U(jNodeZ,iZ1,iZ_B0(2)+iZ2-1,iZ3,iZ4,iCR,iS)

                      END IF

                    END DO
                    END DO
                    END DO
                    END DO

                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO

      END IF

      IF( ApplyOuterBC( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(10) &
        !$OMP PRIVATE( jNodeZ2, iNodeZ, jNodeZ )
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(10) ASYNC &
        !$ACC PRIVATE( jNodeZ2, iNodeZ, jNodeZ ) &
        !$ACC PRESENT( U, iZ_B0, iZ_E0, nNodesZ, NodeNumberTable4D )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(10) &
        !$OMP PRIVATE( jNodeZ2, iNodeZ, jNodeZ )
#endif
        DO iS = 1, nSpecies
          DO iCR = 1, nCR
            DO iZ4 = iZ_B0(4), iZ_E0(4)
              DO iZ3 = iZ_B0(3), iZ_E0(3)
                DO iZ2 = 1, swX1
                  DO iZ1 = iZ_B0(1), iZ_E0(1)

                    DO iNodeZ4 = 1, nNodesZ(4)
                    DO iNodeZ3 = 1, nNodesZ(3)
                    DO iNodeZ2 = 1, nNodesZ(2)
                    DO iNodeZ1 = 1, nNodesZ(1)

                      jNodeZ2 = (nNodesZ(2)-iNodeZ2) + 1

                      iNodeZ = NodeNumberTable4D &
                                 ( iNodeZ1, iNodeZ2, iNodeZ3, iNodeZ4 )
                      jNodeZ = NodeNumberTable4D &
                                 ( iNodeZ1, jNodeZ2, iNodeZ3, iNodeZ4 )

                      ! --- Outer Boundary ---

                      IF( iCR == iCR_G1 )THEN

                        U(iNodeZ,iZ1,iZ_E0(2)+iZ2,iZ3,iZ4,iCR,iS) &
                          = - U(jNodeZ,iZ1,iZ_E0(2)-iZ2+1,iZ3,iZ4,iCR,iS)

                      ELSE

                        U(iNodeZ,iZ1,iZ_E0(2)+iZ2,iZ3,iZ4,iCR,iS) &
                          = U(jNodeZ,iZ1,iZ_E0(2)-iZ2+1,iZ3,iZ4,iCR,iS)

                      END IF

                    END DO
                    END DO
                    END DO
                    END DO

                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO

      END IF

    CASE ( 32 ) ! Reflecting (Inner), Homogeneous (Outer)

      IF( ApplyInnerBC( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(10) &
        !$OMP PRIVATE( jNodeZ2, iNodeZ, jNodeZ )
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(10) ASYNC &
        !$ACC PRIVATE( jNodeZ2, iNodeZ, jNodeZ ) &
        !$ACC PRESENT( U, iZ_B0, iZ_E0, nNodesZ, NodeNumberTable4D )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(10) &
        !$OMP PRIVATE( jNodeZ2, iNodeZ, jNodeZ )
#endif
        DO iS = 1, nSpecies
          DO iCR = 1, nCR
            DO iZ4 = iZ_B0(4), iZ_E0(4)
              DO iZ3 = iZ_B0(3), iZ_E0(3)
                DO iZ2 = 1, swX1
                  DO iZ1 = iZ_B0(1), iZ_E0(1)

                    DO iNodeZ4 = 1, nNodesZ(4)
                    DO iNodeZ3 = 1, nNodesZ(3)
                    DO iNodeZ2 = 1, nNodesZ(2)
                    DO iNodeZ1 = 1, nNodesZ(1)

                      jNodeZ2 = (nNodesZ(2)-iNodeZ2) + 1

                      iNodeZ = NodeNumberTable4D &
                                 ( iNodeZ1, iNodeZ2, iNodeZ3, iNodeZ4 )
                      jNodeZ = NodeNumberTable4D &
                                 ( iNodeZ1, jNodeZ2, iNodeZ3, iNodeZ4 )

                      ! --- Inner Boundary ---

                      IF( iCR == iCR_G1 )THEN

                        U(iNodeZ,iZ1,iZ_B0(2)-iZ2,iZ3,iZ4,iCR,iS) &
                          = - U(jNodeZ,iZ1,iZ_B0(2)+iZ2-1,iZ3,iZ4,iCR,iS)

                      ELSE

                        U(iNodeZ,iZ1,iZ_B0(2)-iZ2,iZ3,iZ4,iCR,iS) &
                          = U(jNodeZ,iZ1,iZ_B0(2)+iZ2-1,iZ3,iZ4,iCR,iS)

                      END IF

                    END DO
                    END DO
                    END DO
                    END DO

                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO

      END IF

      IF( ApplyOuterBC( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) ASYNC &
        !$ACC PRESENT( U, iZ_B0, iZ_E0 )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(7)
#endif
        DO iS = 1, nSpecies
          DO iCR = 1, nCR
            DO iZ4 = iZ_B0(4), iZ_E0(4)
              DO iZ3 = iZ_B0(3), iZ_E0(3)
                DO iZ2 = 1, swX1
                  DO iZ1 = iZ_B0(1), iZ_E0(1)
                    DO iNode = 1, nDOF

                      ! --- Outer Boundary ---

                      IF( iCR == iCR_G1 )THEN

                        U(iNode,iZ1,iZ_E0(2)+iZ2,iZ3,iZ4,iCR,iS) &
                          = MAX( U(iNode,iZ1,iZ_E0(2),iZ3,iZ4,iCR,iS), Zero )

                      ELSE

                        U(iNode,iZ1,iZ_E0(2)+iZ2,iZ3,iZ4,iCR,iS) &
                          = U(iNode,iZ1,iZ_E0(2),iZ3,iZ4,iCR,iS)

                      END IF

                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO

      END IF

    CASE DEFAULT

    END SELECT

  END SUBROUTINE ApplyBC_Radiation_X1


  SUBROUTINE ApplyBC_Radiation_X2( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U, swX2, bcX2, iApplyBC )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4), &
      swX2, bcX2, iApplyBC
    REAL(DP), INTENT(inout) :: &
      U(1:nDOF, &
        iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
        iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4), &
        1:nCR,1:nSpecies)

    INTEGER :: iNode, iS, iCR, iZ1, iZ2, iZ3, iZ4
    INTEGER :: iNodeZ1, iNodeZ2, iNodeZ3, iNodeZ4
    INTEGER :: jNodeZ3, iNodeZ, jNodeZ, iNodeE

    SELECT CASE ( bcX2 )

    CASE ( 0 ) ! No Boundary Condition

    CASE ( 2 ) ! Homogeneous

      IF( ApplyInnerBC( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) ASYNC &
        !$ACC PRESENT( U, iZ_B0, iZ_E0 )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(7)
#endif
        DO iS = 1, nSpecies
          DO iCR = 1, nCR
            DO iZ4 = iZ_B0(4), iZ_E0(4)
              DO iZ3 = 1, swX2
                DO iZ2 = iZ_B0(2), iZ_E0(2)
                  DO iZ1 = iZ_B0(1), iZ_E0(1)
                    DO iNode = 1, nDOF

                      ! --- Inner Boundary ---

                      IF( iCR == iCR_G2 )THEN

                        U(iNode,iZ1,iZ2,iZ_B0(3)-iZ3,iZ4,iCR,iS) &
                          = MIN( U(iNode,iZ1,iZ2,iZ_B0(3),iZ4,iCR,iS), Zero )

                      ELSE

                        U(iNode,iZ1,iZ2,iZ_B0(3)-iZ3,iZ4,iCR,iS) &
                          = U(iNode,iZ1,iZ2,iZ_B0(3),iZ4,iCR,iS)

                      END IF

                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO

      END IF

      IF( ApplyOuterBC( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) ASYNC &
        !$ACC PRESENT( U, iZ_B0, iZ_E0 )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(7)
#endif
        DO iS = 1, nSpecies
          DO iCR = 1, nCR
            DO iZ4 = iZ_B0(4), iZ_E0(4)
              DO iZ3 = 1, swX2
                DO iZ2 = iZ_B0(2), iZ_E0(2)
                  DO iZ1 = iZ_B0(1), iZ_E0(1)
                    DO iNode = 1, nDOF

                      ! --- Outer Boundary ---

                      IF( iCR == iCR_G2 )THEN

                        U(iNode,iZ1,iZ2,iZ_E0(3)+iZ3,iZ4,iCR,iS) &
                          = MAX( U(iNode,iZ1,iZ2,iZ_E0(3),iZ4,iCR,iS), Zero )

                      ELSE

                        U(iNode,iZ1,iZ2,iZ_E0(3)+iZ3,iZ4,iCR,iS) &
                          = U(iNode,iZ1,iZ2,iZ_E0(3),iZ4,iCR,iS)

                      END IF

                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO

      END IF

    CASE ( 3 ) ! Reflecting

      IF( ApplyInnerBC( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(10) &
        !$OMP PRIVATE( jNodeZ3, iNodeZ, jNodeZ )
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(10) ASYNC &
        !$ACC PRIVATE( jNodeZ3, iNodeZ, jNodeZ ) &
        !$ACC PRESENT( U, iZ_B0, iZ_E0, nNodesZ, NodeNumberTable4D )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(10) &
        !$OMP PRIVATE( jNodeZ3, iNodeZ, jNodeZ )
#endif
        DO iS = 1, nSpecies
          DO iCR = 1, nCR
            DO iZ4 = iZ_B0(4), iZ_E0(4)
              DO iZ3 = 1, swX2
                DO iZ2 = iZ_B0(2), iZ_E0(2)
                  DO iZ1 = iZ_B0(1), iZ_E0(1)

                    DO iNodeZ4 = 1, nNodesZ(4)
                    DO iNodeZ3 = 1, nNodesZ(3)
                    DO iNodeZ2 = 1, nNodesZ(2)
                    DO iNodeZ1 = 1, nNodesZ(1)

                      jNodeZ3 = (nNodesZ(3)-iNodeZ3) + 1

                      iNodeZ = NodeNumberTable4D &
                                 ( iNodeZ1, iNodeZ2, iNodeZ3, iNodeZ4 )
                      jNodeZ = NodeNumberTable4D &
                                 ( iNodeZ1, iNodeZ2, jNodeZ3, iNodeZ4 )

                      ! --- Inner Boundary ---

                      IF( iCR == iCR_G2 )THEN

                        U(iNodeZ,iZ1,iZ2,iZ_B0(3)-iZ3,iZ4,iCR,iS) &
                          = - U(jNodeZ,iZ1,iZ2,iZ_B0(3)+iZ3-1,iZ4,iCR,iS)

                      ELSE

                        U(iNodeZ,iZ1,iZ2,iZ_B0(3)-iZ3,iZ4,iCR,iS) &
                          = U(jNodeZ,iZ1,iZ2,iZ_B0(3)+iZ3-1,iZ4,iCR,iS)

                      END IF

                    END DO
                    END DO
                    END DO
                    END DO

                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO

      END IF

      IF( ApplyOuterBC( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(10) &
        !$OMP PRIVATE( jNodeZ3, iNodeZ, jNodeZ )
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(10) ASYNC &
        !$ACC PRIVATE( jNodeZ3, iNodeZ, jNodeZ ) &
        !$ACC PRESENT( U, iZ_B0, iZ_E0, nNodesZ, NodeNumberTable4D )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(10) &
        !$OMP PRIVATE( jNodeZ3, iNodeZ, jNodeZ )
#endif
        DO iS = 1, nSpecies
          DO iCR = 1, nCR
            DO iZ4 = iZ_B0(4), iZ_E0(4)
              DO iZ3 = 1, swX2
                DO iZ2 = iZ_B0(2), iZ_E0(2)
                  DO iZ1 = iZ_B0(1), iZ_E0(1)

                    DO iNodeZ4 = 1, nNodesZ(4)
                    DO iNodeZ3 = 1, nNodesZ(3)
                    DO iNodeZ2 = 1, nNodesZ(2)
                    DO iNodeZ1 = 1, nNodesZ(1)

                      jNodeZ3 = (nNodesZ(3)-iNodeZ3) + 1

                      iNodeZ = NodeNumberTable4D &
                                 ( iNodeZ1, iNodeZ2, iNodeZ3, iNodeZ4 )
                      jNodeZ = NodeNumberTable4D &
                                 ( iNodeZ1, iNodeZ2, jNodeZ3, iNodeZ4 )


                      ! --- Outer Boundary ---

                      IF( iCR == iCR_G2 )THEN

                        U(iNodeZ,iZ1,iZ2,iZ_E0(3)+iZ3,iZ4,iCR,iS) &
                          = - U(jNodeZ,iZ1,iZ2,iZ_E0(3)-iZ3+1,iZ4,iCR,iS)

                      ELSE

                        U(iNodeZ,iZ1,iZ2,iZ_E0(3)+iZ3,iZ4,iCR,iS) &
                          = U(jNodeZ,iZ1,iZ2,iZ_E0(3)-iZ3+1,iZ4,iCR,iS)

                      END IF

                    END DO
                    END DO
                    END DO
                    END DO

                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO

      END IF

    CASE ( 32 ) ! Reflecting (Inner), Homogeneous (Outer)

      IF( ApplyInnerBC( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(10) &
        !$OMP PRIVATE( jNodeZ3, iNodeZ, jNodeZ )
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(10) ASYNC &
        !$ACC PRIVATE( jNodeZ3, iNodeZ, jNodeZ ) &
        !$ACC PRESENT( U, iZ_B0, iZ_E0, nNodesZ, NodeNumberTable4D )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(10) &
        !$OMP PRIVATE( jNodeZ3, iNodeZ, jNodeZ )
#endif
        DO iS = 1, nSpecies
          DO iCR = 1, nCR
            DO iZ4 = iZ_B0(4), iZ_E0(4)
              DO iZ3 = 1, swX2
                DO iZ2 = iZ_B0(2), iZ_E0(2)
                  DO iZ1 = iZ_B0(1), iZ_E0(1)

                    DO iNodeZ4 = 1, nNodesZ(4)
                    DO iNodeZ3 = 1, nNodesZ(3)
                    DO iNodeZ2 = 1, nNodesZ(2)
                    DO iNodeZ1 = 1, nNodesZ(1)

                      jNodeZ3 = (nNodesZ(3)-iNodeZ3) + 1

                      iNodeZ = NodeNumberTable4D &
                                 ( iNodeZ1, iNodeZ2, iNodeZ3, iNodeZ4 )
                      jNodeZ = NodeNumberTable4D &
                                 ( iNodeZ1, iNodeZ2, jNodeZ3, iNodeZ4 )

                      ! --- Inner Boundary ---

                      IF( iCR == iCR_G2 )THEN

                        U(iNodeZ,iZ1,iZ2,iZ_B0(3)-iZ3,iZ4,iCR,iS) &
                          = - U(jNodeZ,iZ1,iZ2,iZ_B0(3)+iZ3-1,iZ4,iCR,iS)

                      ELSE

                        U(iNodeZ,iZ1,iZ2,iZ_B0(3)-iZ3,iZ4,iCR,iS) &
                          = U(jNodeZ,iZ1,iZ2,iZ_B0(3)+iZ3-1,iZ4,iCR,iS)

                      END IF

                    END DO
                    END DO
                    END DO
                    END DO

                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO

      END IF

      IF( ApplyOuterBC( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) ASYNC &
        !$ACC PRESENT( U, iZ_B0, iZ_E0 )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(7)
#endif
        DO iS = 1, nSpecies
          DO iCR = 1, nCR
            DO iZ4 = iZ_B0(4), iZ_E0(4)
              DO iZ3 = 1, swX2
                DO iZ2 = iZ_B0(2), iZ_E0(2)
                  DO iZ1 = iZ_B0(1), iZ_E0(1)
                    DO iNode = 1, nDOF

                      ! --- Outer Boundary ---

                      IF( iCR == iCR_G2 )THEN

                        U(iNode,iZ1,iZ2,iZ_E0(3)+iZ3,iZ4,iCR,iS) &
                          = MAX( U(iNode,iZ1,iZ2,iZ_E0(3),iZ4,iCR,iS), Zero )

                      ELSE

                        U(iNode,iZ1,iZ2,iZ_E0(3)+iZ3,iZ4,iCR,iS) &
                          = U(iNode,iZ1,iZ2,iZ_E0(3),iZ4,iCR,iS)

                      END IF

                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO

      END IF

    CASE DEFAULT

    END SELECT

  END SUBROUTINE ApplyBC_Radiation_X2


  SUBROUTINE ApplyBC_Radiation_X3( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U, swX3, bcX3, iApplyBC )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4), &
      swX3, bcX3, iApplyBC
    REAL(DP), INTENT(inout) :: &
      U(1:nDOF, &
        iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
        iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4), &
        1:nCR,1:nSpecies)

    INTEGER :: iNode, iS, iCR, iZ1, iZ2, iZ3, iZ4
    INTEGER :: iNodeZ1, iNodeZ2, iNodeZ3, iNodeZ4
    INTEGER :: jNodeZ4, iNodeZ, jNodeZ, iNodeE

    SELECT CASE ( bcX3 )

    CASE ( 0 ) ! No Boundary Condition

    CASE ( 2 ) ! Homogeneous

      IF( ApplyInnerBC( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) ASYNC &
        !$ACC PRESENT( U, iZ_B0, iZ_E0 )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(7)
#endif
        DO iS = 1, nSpecies
          DO iCR = 1, nCR
            DO iZ4 = 1, swX3
              DO iZ3 = iZ_B0(3), iZ_E0(3)
                DO iZ2 = iZ_B0(2), iZ_E0(2)
                  DO iZ1 = iZ_B0(1), iZ_E0(1)
                    DO iNode = 1, nDOF

                      ! --- Inner Boundary ---

                      IF( iCR == iCR_G3 )THEN

                        U(iNode,iZ1,iZ2,iZ3,iZ_B0(4)-iZ4,iCR,iS) &
                          = MIN( U(iNode,iZ1,iZ2,iZ3,iZ_B0(4),iCR,iS), Zero )

                      ELSE

                        U(iNode,iZ1,iZ2,iZ3,iZ_B0(4)-iZ4,iCR,iS) &
                          = U(iNode,iZ1,iZ2,iZ3,iZ_B0(4),iCR,iS)

                      END IF

                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO

      END IF

      IF( ApplyOuterBC( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) ASYNC &
        !$ACC PRESENT( U, iZ_B0, iZ_E0 )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(7)
#endif
        DO iS = 1, nSpecies
          DO iCR = 1, nCR
            DO iZ4 = 1, swX3
              DO iZ3 = iZ_B0(3), iZ_E0(3)
                DO iZ2 = iZ_B0(2), iZ_E0(2)
                  DO iZ1 = iZ_B0(1), iZ_E0(1)
                    DO iNode = 1, nDOF

                      ! --- Outer Boundary ---

                      IF( iCR == iCR_G3 )THEN

                        U(iNode,iZ1,iZ2,iZ3,iZ_E0(4)+iZ4,iCR,iS) &
                          = MAX( U(iNode,iZ1,iZ2,iZ3,iZ_E0(4),iCR,iS), Zero )

                      ELSE

                        U(iNode,iZ1,iZ2,iZ3,iZ_E0(4)+iZ4,iCR,iS) &
                          = U(iNode,iZ1,iZ2,iZ3,iZ_E0(4),iCR,iS)

                      END IF

                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO

      END IF

    CASE ( 3 ) ! Reflecting

      IF( ApplyInnerBC( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(10) &
        !$OMP PRIVATE( jNodeZ4, iNodeZ, jNodeZ )
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(10) ASYNC &
        !$ACC PRIVATE( jNodeZ4, iNodeZ, jNodeZ ) &
        !$ACC PRESENT( U, iZ_B0, iZ_E0, nNodesZ, NodeNumberTable4D )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(10) &
        !$OMP PRIVATE( jNodeZ4, iNodeZ, jNodeZ )
#endif
        DO iS = 1, nSpecies
          DO iCR = 1, nCR
            DO iZ4 = 1, swX3
              DO iZ3 = iZ_B0(3), iZ_E0(3)
                DO iZ2 = iZ_B0(2), iZ_E0(2)
                  DO iZ1 = iZ_B0(1), iZ_E0(1)

                    DO iNodeZ4 = 1, nNodesZ(4)
                    DO iNodeZ3 = 1, nNodesZ(3)
                    DO iNodeZ2 = 1, nNodesZ(2)
                    DO iNodeZ1 = 1, nNodesZ(1)

                      jNodeZ4 = (nNodesZ(4)-iNodeZ4) + 1

                      iNodeZ = NodeNumberTable4D &
                                 ( iNodeZ1, iNodeZ2, iNodeZ3, iNodeZ4 )
                      jNodeZ = NodeNumberTable4D &
                                 ( iNodeZ1, iNodeZ2, iNodeZ3, jNodeZ4 )

                      ! --- Inner Boundary ---

                      IF( iCR == iCR_G3 )THEN

                        U(iNodeZ,iZ1,iZ2,iZ3,iZ_B0(4)-iZ4,iCR,iS) &
                          = - U(jNodeZ,iZ1,iZ2,iZ3,iZ_B0(4)+iZ4-1,iCR,iS)

                      ELSE

                        U(iNodeZ,iZ1,iZ2,iZ3,iZ_B0(4)-iZ4,iCR,iS) &
                          = U(jNodeZ,iZ1,iZ2,iZ3,iZ_B0(4)+iZ4-1,iCR,iS)

                      END IF

                    END DO
                    END DO
                    END DO
                    END DO

                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO

      END IF

      IF( ApplyOuterBC( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(10) &
        !$OMP PRIVATE( jNodeZ4, iNodeZ, jNodeZ )
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(10) ASYNC &
        !$ACC PRIVATE( jNodeZ4, iNodeZ, jNodeZ ) &
        !$ACC PRESENT( U, iZ_B0, iZ_E0, nNodesZ, NodeNumberTable4D )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(10) &
        !$OMP PRIVATE( jNodeZ4, iNodeZ, jNodeZ )
#endif
        DO iS = 1, nSpecies
          DO iCR = 1, nCR
            DO iZ4 = 1, swX3
              DO iZ3 = iZ_B0(3), iZ_E0(3)
                DO iZ2 = iZ_B0(2), iZ_E0(2)
                  DO iZ1 = iZ_B0(1), iZ_E0(1)

                    DO iNodeZ4 = 1, nNodesZ(4)
                    DO iNodeZ3 = 1, nNodesZ(3)
                    DO iNodeZ2 = 1, nNodesZ(2)
                    DO iNodeZ1 = 1, nNodesZ(1)

                      jNodeZ4 = (nNodesZ(4)-iNodeZ4) + 1

                      iNodeZ = NodeNumberTable4D &
                                 ( iNodeZ1, iNodeZ2, iNodeZ3, iNodeZ4 )
                      jNodeZ = NodeNumberTable4D &
                                 ( iNodeZ1, iNodeZ2, iNodeZ3, jNodeZ4 )


                      ! --- Outer Boundary ---

                      IF( iCR == iCR_G3 )THEN

                        U(iNodeZ,iZ1,iZ2,iZ3,iZ_E0(4)+iZ4,iCR,iS) &
                          = - U(jNodeZ,iZ1,iZ2,iZ3,iZ_E0(4)-iZ4+1,iCR,iS)

                      ELSE

                        U(iNodeZ,iZ1,iZ2,iZ3,iZ_E0(4)+iZ4,iCR,iS) &
                          = U(jNodeZ,iZ1,iZ2,iZ3,iZ_E0(4)-iZ4+1,iCR,iS)

                      END IF

                    END DO
                    END DO
                    END DO
                    END DO

                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO

      END IF

    CASE ( 32 ) ! Reflecting (Inner), Homogeneous (Outer)

      IF( ApplyInnerBC( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(10) &
        !$OMP PRIVATE( jNodeZ4, iNodeZ, jNodeZ )
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(10) ASYNC &
        !$ACC PRIVATE( jNodeZ4, iNodeZ, jNodeZ ) &
        !$ACC PRESENT( U, iZ_B0, iZ_E0, nNodesZ, NodeNumberTable4D )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(10) &
        !$OMP PRIVATE( jNodeZ4, iNodeZ, jNodeZ )
#endif
        DO iS = 1, nSpecies
          DO iCR = 1, nCR
            DO iZ4 = 1, swX3
              DO iZ3 = iZ_B0(3), iZ_E0(3)
                DO iZ2 = iZ_B0(2), iZ_E0(2)
                  DO iZ1 = iZ_B0(1), iZ_E0(1)

                    DO iNodeZ4 = 1, nNodesZ(4)
                    DO iNodeZ3 = 1, nNodesZ(3)
                    DO iNodeZ2 = 1, nNodesZ(2)
                    DO iNodeZ1 = 1, nNodesZ(1)

                      jNodeZ4 = (nNodesZ(4)-iNodeZ4) + 1

                      iNodeZ = NodeNumberTable4D &
                                 ( iNodeZ1, iNodeZ2, iNodeZ3, iNodeZ4 )
                      jNodeZ = NodeNumberTable4D &
                                 ( iNodeZ1, iNodeZ2, iNodeZ3, jNodeZ4 )

                      ! --- Inner Boundary ---

                      IF( iCR == iCR_G3 )THEN

                        U(iNodeZ,iZ1,iZ2,iZ3,iZ_B0(4)-iZ4,iCR,iS) &
                          = - U(jNodeZ,iZ1,iZ2,iZ3,iZ_B0(4)+iZ4-1,iCR,iS)

                      ELSE

                        U(iNodeZ,iZ1,iZ2,iZ3,iZ_B0(4)-iZ4,iCR,iS) &
                          = U(jNodeZ,iZ1,iZ2,iZ3,iZ_B0(4)+iZ4-1,iCR,iS)

                      END IF

                    END DO
                    END DO
                    END DO
                    END DO

                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO

      END IF

      IF( ApplyOuterBC( iApplyBC ) )THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) ASYNC &
        !$ACC PRESENT( U, iZ_B0, iZ_E0 )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO COLLAPSE(7)
#endif
        DO iS = 1, nSpecies
          DO iCR = 1, nCR
            DO iZ4 = 1, swX3
              DO iZ3 = iZ_B0(3), iZ_E0(3)
                DO iZ2 = iZ_B0(2), iZ_E0(2)
                  DO iZ1 = iZ_B0(1), iZ_E0(1)
                    DO iNode = 1, nDOF

                      ! --- Outer Boundary ---

                      IF( iCR == iCR_G3 )THEN

                        U(iNode,iZ1,iZ2,iZ3,iZ_E0(4)+iZ4,iCR,iS) &
                          = MAX( U(iNode,iZ1,iZ2,iZ3,iZ_E0(4),iCR,iS), Zero )

                      ELSE

                        U(iNode,iZ1,iZ2,iZ3,iZ_E0(4)+iZ4,iCR,iS) &
                          = U(iNode,iZ1,iZ2,iZ3,iZ_E0(4),iCR,iS)

                      END IF

                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO

      END IF

    CASE DEFAULT

    END SELECT

  END SUBROUTINE ApplyBC_Radiation_X3

END MODULE TimeSteppingModule_Flash
