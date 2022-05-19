MODULE TwoMoment_TimeSteppingModule_OrderV

  USE KindModule, ONLY: &
    DP, Zero, One
  USE ProgramHeaderModule, ONLY: &
    nDOFZ, nDOFX, nDOFE, &
    iZ_B0, iZ_B1, iZ_E0, iZ_E1, &
    iX_B0, iX_B1, iX_E0, iX_E1, &
    iE_B0, iE_B1, iE_E0, iE_E1
  USE GeometryFieldsModuleE, ONLY: &
    nGE
  USE GeometryFieldsModule, ONLY: &
    nGF
  USE FluidFieldsModule, ONLY: &
    nCF, iCF_D, uDF
  USE Euler_SlopeLimiterModule_NonRelativistic_TABLE, ONLY: &
    ApplySlopeLimiter_Euler_NonRelativistic_TABLE
  USE Euler_PositivityLimiterModule_NonRelativistic_TABLE, ONLY: &
    ApplyPositivityLimiter_Euler_NonRelativistic_TABLE
  USE Euler_dgDiscretizationModule, ONLY: &
    ComputeIncrement_Euler_DG_Explicit, &
    OffGridFlux_Euler_X1_Inner, &
    OffGridFlux_Euler_X1_Outer, &
    OffGridFlux_Euler_X2_Inner, &
    OffGridFlux_Euler_X2_Outer, &
    OffGridFlux_Euler_X3_Inner, &
    OffGridFlux_Euler_X3_Outer
  USE RadiationFieldsModule, ONLY: &
    nCR, nSpecies
  USE TwoMoment_TimersModule_OrderV, ONLY: &
    TimersStart, &
    TimersStop, &
    Timer_IMEX, &
    Timer_TimeStepper
  USE TwoMoment_SlopeLimiterModule_OrderV, ONLY: &
    ApplySlopeLimiter_TwoMoment
  USE TwoMoment_PositivityLimiterModule_OrderV, ONLY: &
    ApplyPositivityLimiter_TwoMoment, &
    dEnergyMomentum_PL_TwoMoment
  USE TwoMoment_DiscretizationModule_Streaming_OrderV, ONLY: &
    ComputeIncrement_TwoMoment_Explicit, &
    OffGridFlux_TwoMoment
  USE TwoMoment_TallyModule_OrderV, ONLY: &
    IncrementOffGridTally_Euler, &
    IncrementOffGridTally_TwoMoment, &
    IncrementPositivityLimiterTally_TwoMoment

  IMPLICIT NONE
  PRIVATE

  TYPE :: StageDataType
    REAL(DP)              :: OffGridFlux_U(nCF)
    REAL(DP)              :: OffGridFlux_M(2*nCR)
    REAL(DP), ALLOCATABLE :: dU_IM(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: dU_EX(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: dM_IM(:,:,:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: dM_EX(:,:,:,:,:,:,:)
  END TYPE StageDataType

  LOGICAL                          :: EvolveEuler
  LOGICAL                          :: EvolveTwoMoment
  INTEGER                          :: nStages
  REAL(DP),            ALLOCATABLE :: c_IM(:), w_IM(:), a_IM(:,:)
  REAL(DP),            ALLOCATABLE :: c_EX(:), w_EX(:), a_EX(:,:)
  REAL(DP),            ALLOCATABLE :: U0(:,:,:,:,:)
  REAL(DP),            ALLOCATABLE :: Ui(:,:,:,:,:)
  REAL(DP),            ALLOCATABLE :: M0(:,:,:,:,:,:,:)
  REAL(DP),            ALLOCATABLE :: Mi(:,:,:,:,:,:,:)
  TYPE(StageDataType), ALLOCATABLE :: StageData(:)

  PUBLIC :: Initialize_IMEX_RK
  PUBLIC :: Finalize_IMEX_RK
  PUBLIC :: Update_IMEX_RK

#if   defined( THORNADO_OMP_OL )
#elif defined( THORNADO_OACC   )
#endif

  INTERFACE
    SUBROUTINE ImplicitIncrement &
      ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, dt, GE, GX, U_F, dU_F, U_R, dU_R )
      USE KindModule           , ONLY: DP
      USE ProgramHeaderModule  , ONLY: nDOFX, nDOFE, nDOFZ
      USE GeometryFieldsModuleE, ONLY: nGE
      USE GeometryFieldsModule , ONLY: nGF
      USE FluidFieldsModule    , ONLY: nCF
      USE RadiationFieldsModule, ONLY: nCR, nSpecies
      INTEGER,  INTENT(in)    :: &
        iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
      REAL(DP), INTENT(in)    :: &
        dt
      REAL(DP), INTENT(in)    :: &
        GE  (1:nDOFE, &
             iZ_B1(1):iZ_E1(1), &
             1:nGE)
      REAL(DP), INTENT(in)    :: &
        GX  (1:nDOFX, &
             iZ_B1(2):iZ_E1(2), &
             iZ_B1(3):iZ_E1(3), &
             iZ_B1(4):iZ_E1(4), &
             1:nGF)
      REAL(DP), INTENT(in) :: &
        U_F (1:nDOFX, &
             iZ_B1(2):iZ_E1(2), &
             iZ_B1(3):iZ_E1(3), &
             iZ_B1(4):iZ_E1(4), &
             1:nCF)
      REAL(DP), INTENT(out) :: &
        dU_F(1:nDOFX, &
             iZ_B1(2):iZ_E1(2), &
             iZ_B1(3):iZ_E1(3), &
             iZ_B1(4):iZ_E1(4), &
             1:nCF)
      REAL(DP), INTENT(in) :: &
        U_R (1:nDOFZ, &
             iZ_B1(1):iZ_E1(1), &
             iZ_B1(2):iZ_E1(2), &
             iZ_B1(3):iZ_E1(3), &
             iZ_B1(4):iZ_E1(4), &
             1:nCR, &
             1:nSpecies)
      REAL(DP), INTENT(out) :: &
        dU_R(1:nDOFZ, &
             iZ_B1(1):iZ_E1(1), &
             iZ_B1(2):iZ_E1(2), &
             iZ_B1(3):iZ_E1(3), &
             iZ_B1(4):iZ_E1(4), &
             1:nCR, &
             1:nSpecies)
    END SUBROUTINE ImplicitIncrement
  END INTERFACE

  INTERFACE
    SUBROUTINE GravitySolver( iX_B0, iX_E0, iX_B1, iX_E1, GX, D )
      USE KindModule          , ONLY: DP
      USE ProgramHeaderModule , ONLY: nDOFX
      USE GeometryFieldsModule, ONLY: nGF
      INTEGER,  INTENT(in)    :: &
        iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
      REAL(DP), INTENT(inout) :: &
        GX(1:, &
           iX_B1(1):, &
           iX_B1(2):, &
           iX_B1(3):, &
           1:)
      REAL(DP), INTENT(in)    :: &
        D (1:, &
           iX_B1(1):, &
           iX_B1(2):, &
           iX_B1(3):)
    END SUBROUTINE GravitySolver
  END INTERFACE

  INTERFACE AllocateArray
    MODULE PROCEDURE AllocateArray5D
    MODULE PROCEDURE AllocateArray7D
  END INTERFACE AllocateArray

  INTERFACE DeallocateArray
    MODULE PROCEDURE DeallocateArray5D
    MODULE PROCEDURE DeallocateArray7D
  END INTERFACE DeallocateArray

  INTERFACE CopyArray
    MODULE PROCEDURE CopyArray5D
    MODULE PROCEDURE CopyArray7D
  END INTERFACE CopyArray

  INTERFACE AddToArray
    MODULE PROCEDURE AddToArray5D
    MODULE PROCEDURE AddToArray7D
  END INTERFACE AddToArray

CONTAINS


  SUBROUTINE Update_IMEX_RK &
    ( dt, GE, GX, U, M, ComputeIncrement_TwoMoment_Implicit, SolveGravity )

    REAL(DP), INTENT(in) :: &
      dt
    REAL(DP), INTENT(in) :: &
      GE(1:nDOFE, &
         iE_B1:iE_E1, &
         1:nGE)
    REAL(DP), INTENT(inout) :: &
      GX(1:nDOFX, &
         iX_B1(1):iX_E1(1), &
         iX_B1(2):iX_E1(2), &
         iX_B1(3):iX_E1(3), &
         1:nGF)
    REAL(DP), INTENT(inout) :: &
      U(1:nDOFX, &
        iX_B1(1):iX_E1(1), &
        iX_B1(2):iX_E1(2), &
        iX_B1(3):iX_E1(3), &
        1:nCF)
    REAL(DP), INTENT(inout) :: &
      M(1:nDOFZ, &
        iZ_B1(1):iZ_E1(1), &
        iZ_B1(2):iZ_E1(2), &
        iZ_B1(3):iZ_E1(3), &
        iZ_B1(4):iZ_E1(4), &
        1:nCR,1:nSpecies)
    PROCEDURE(ImplicitIncrement) :: &
      ComputeIncrement_TwoMoment_Implicit
    PROCEDURE(GravitySolver), OPTIONAL :: &
      SolveGravity

    INTEGER  :: iS, jS
    REAL(DP) :: dU_OffGrid(nCF)
    REAL(DP) :: dM_OffGrid(2*nCR)
    REAL(DP) :: dM_PL(nCR)

    CALL TimersStart( Timer_IMEX )

    dM_PL = Zero

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: GE, GX, U, M )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA &
    !$ACC COPYIN( GE, GX, U, M )
#endif

    CALL CopyArray( U0, One, U )
    CALL CopyArray( M0, One, M )

    DO iS = 1, nStages

      CALL CopyArray( Ui, One, U0 )
      CALL CopyArray( Mi, One, M0 )

      DO jS = 1, iS - 1

        IF( a_IM(iS,jS) .NE. Zero )THEN

          CALL AddToArray( One, Ui, dt * a_IM(iS,jS), StageData(jS) % dU_IM )
          CALL AddToArray( One, Mi, dt * a_IM(iS,jS), StageData(jS) % dM_IM )

        END IF

        IF( a_EX(iS,jS) .NE. Zero )THEN

          CALL AddToArray( One, Ui, dt * a_EX(iS,jS), StageData(jS) % dU_EX )
          CALL AddToArray( One, Mi, dt * a_EX(iS,jS), StageData(jS) % dM_EX )

        END IF

        IF( jS == iS - 1 )THEN

          ! --- Apply Limiters ---

          IF( EvolveEuler )THEN

            CALL ApplySlopeLimiter_Euler_NonRelativistic_TABLE &
                   ( iX_B0, iX_E0, iX_B1, iX_E1, GX, Ui, uDF )

            CALL ApplyPositivityLimiter_Euler_NonRelativistic_TABLE &
                   ( iX_B0, iX_E0, iX_B1, iX_E1, GX, Ui, uDF )

            IF( PRESENT( SolveGravity ) )THEN

              CALL SolveGravity &
                     ( iX_B0, iX_E0, iX_B1, iX_E1, GX, Ui(:,:,:,:,iCF_D) )

            END IF

          END IF

          IF( EvolveTwoMoment )THEN

            CALL ApplySlopeLimiter_TwoMoment &
                   ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, Ui, Mi )

            CALL ApplyPositivityLimiter_TwoMoment &
                   ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, Ui, Mi )

            dM_PL = dM_PL + dEnergyMomentum_PL_TwoMoment

          END IF

        END IF

      END DO ! jS = 1, iS - 1

      IF( ANY( a_IM(:,iS) .NE. Zero ) .OR. ( w_IM(iS) .NE. Zero ) )THEN

        ! --- Implicit Solve ---

        CALL ComputeIncrement_TwoMoment_Implicit &
               ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, dt * a_IM(iS,iS), GE, GX, &
                 Ui, StageData(iS) % dU_IM, Mi, StageData(iS) % dM_IM )

        CALL AddToArray( One, Ui, dt * a_IM(iS,iS), StageData(iS) % dU_IM )
        CALL AddToArray( One, Mi, dt * a_IM(iS,iS), StageData(iS) % dM_IM )

        IF( EvolveEuler )THEN

          CALL ApplyPositivityLimiter_Euler_NonRelativistic_TABLE &
                 ( iX_B0, iX_E0, iX_B1, iX_E1, GX, Ui, uDF )

        END IF

        IF( EvolveTwoMoment )THEN

          CALL ApplyPositivityLimiter_TwoMoment &
                 ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, Ui, Mi )

          dM_PL = dM_PL + dEnergyMomentum_PL_TwoMoment

        END IF

      END IF

      IF( ANY( a_EX(:,iS) .NE. Zero ) .OR. ( w_EX(iS) .NE. Zero ) )THEN

        ! --- Explicit Solve ---

        IF( EvolveEuler )THEN

          CALL ComputeIncrement_Euler_DG_Explicit &
                 ( iX_B0, iX_E0, iX_B1, iX_E1, GX, &
                   Ui, uDF, StageData(iS) % dU_EX )

          StageData(iS) % OffGridFlux_U &
                  =   OffGridFlux_Euler_X1_Outer &
                    - OffGridFlux_Euler_X1_Inner &
                    + OffGridFlux_Euler_X2_Outer &
                    - OffGridFlux_Euler_X2_Inner &
                    + OffGridFlux_Euler_X3_Outer &
                    - OffGridFlux_Euler_X3_Inner

        END IF

        IF( EvolveTwoMoment )THEN

          CALL ComputeIncrement_TwoMoment_Explicit &
                 ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, &
                   Ui, Mi, StageData(iS) % dM_EX )

          StageData(iS) % OffGridFlux_M = OffGridFlux_TwoMoment

        END IF

      END IF

    END DO ! iS = 1, nStages

    ! --- Assembly Step ---

    IF( ANY( a_IM(nStages,:) .NE. w_IM(:) ) .OR. &
        ANY( a_EX(nStages,:) .NE. w_EX(:) ) )THEN

      CALL CopyArray( Ui, One, U0 )
      CALL CopyArray( Mi, One, M0 )

      DO iS = 1, nStages

        IF( w_IM(iS) .NE. Zero )THEN

          CALL AddToArray( One, Ui, dt * w_IM(iS), StageData(iS) % dU_IM )
          CALL AddToArray( One, Mi, dt * w_IM(iS), StageData(iS) % dM_IM )

        END IF

        IF( w_EX(iS) .NE. Zero )THEN

          CALL AddToArray( One, Ui, dt * w_EX(iS), StageData(iS) % dU_EX )
          CALL AddToArray( One, Mi, dt * w_EX(iS), StageData(iS) % dM_EX )

        END IF

      END DO

      IF( EvolveEuler )THEN

        CALL ApplySlopeLimiter_Euler_NonRelativistic_TABLE &
               ( iX_B0, iX_E0, iX_B1, iX_E1, GX, Ui, uDF )

        CALL ApplyPositivityLimiter_Euler_NonRelativistic_TABLE &
               ( iX_B0, iX_E0, iX_B1, iX_E1, GX, Ui, uDF )

      END IF

      IF( EvolveTwoMoment )THEN

        CALL ApplySlopeLimiter_TwoMoment &
               ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, Ui, Mi )

        CALL ApplyPositivityLimiter_TwoMoment &
               ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, Ui, Mi )

        dM_PL = dM_PL + dEnergyMomentum_PL_TwoMoment

      END IF

    END IF

    IF( EvolveEuler .AND. PRESENT( SolveGravity ) )THEN

      CALL SolveGravity &
             ( iX_B0, iX_E0, iX_B1, iX_E1, GX, Ui(:,:,:,:,iCF_D) )

    END IF

    CALL CopyArray( U, One, Ui )
    CALL CopyArray( M, One, Mi )

    IF( EvolveEuler )THEN

      dU_OffGrid = Zero
      DO iS = 1, nStages

        dU_OffGrid &
          = dU_OffGrid + dt * w_EX(iS) * StageData(iS) % OffGridFlux_U

      END DO

      CALL IncrementOffGridTally_Euler( dU_OffGrid )

    END IF

    IF( EvolveTwoMoment )THEN

      dM_OffGrid = Zero
      DO iS = 1, nStages

        dM_OffGrid &
          = dM_OffGrid + dt * w_EX(iS) * StageData(iS) % OffGridFlux_M

      END DO

      CALL IncrementOffGridTally_TwoMoment( dM_OffGrid )

      CALL IncrementPositivityLimiterTally_TwoMoment( dM_PL )

    END IF

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from: U, M ) &
    !$OMP MAP( release: GE, GX )
#elif defined( THORNADO_OACC   )
    !$ACC EXIT DATA &
    !$ACC COPYOUT( U, M ) &
    !$ACC DELETE( GE, GX )
#endif

    CALL TimersStop( Timer_IMEX )

  END SUBROUTINE Update_IMEX_RK


  SUBROUTINE Initialize_IMEX_RK &
    ( Scheme, EvolveEuler_Option, EvolveTwoMoment_Option )

    CHARACTER(LEN=*), INTENT(in)           :: Scheme
    LOGICAL         , INTENT(in), OPTIONAL :: EvolveEuler_Option
    LOGICAL         , INTENT(in), OPTIONAL :: EvolveTwoMoment_Option

    INTEGER :: i

    EvolveEuler = .FALSE.
    IF( PRESENT( EvolveEuler_Option ) )THEN
      EvolveEuler = EvolveEuler_Option
    END IF

    EvolveTwoMoment = .TRUE.
    IF( PRESENT( EvolveTwoMoment_Option ) )THEN
      EvolveTwoMoment = EvolveTwoMoment_Option
    END IF

    WRITE(*,*)
    WRITE(*,'(A6,A19,A)' ) '', 'IMEX-RK Scheme: '   , TRIM( Scheme )
    WRITE(*,'(A6,A19,L1)') '', 'Evolve Euler: '     , EvolveEuler
    WRITE(*,'(A6,A19,L1)') '', 'Evolve Two-Moment: ', EvolveTwoMoment

    SELECT CASE ( TRIM( Scheme ) )

      CASE ( 'SSPRK1' )

        nStages = 1

        CALL AllocateButcherTables

        a_EX(1,1) = 0.0_DP
        w_EX(1)   = 1.0_DP

      CASE ( 'SSPRK2' )

        nStages = 2

        CALL AllocateButcherTables

        a_EX(1,1:2) = [ 0.0_DP, 0.0_DP ]
        a_EX(2,1:2) = [ 1.0_DP, 0.0_DP ]
        w_EX(1:2)   = [ 0.5_DP, 0.5_DP ]

      CASE ( 'SSPRK3' )

        nStages = 3

        CALL AllocateButcherTables

        a_EX(1,1:3) = [ 0.00_DP, 0.00_DP, 0.00_DP ]
        a_EX(2,1:3) = [ 1.00_DP, 0.00_DP, 0.00_DP ]
        a_EX(3,1:3) = [ 0.25_DP, 0.25_DP, 0.00_DP ]
        w_EX(1:3)   = [ 1.00_DP, 1.00_DP, 4.00_DP ] / 6.0_DP

      CASE ( 'BackwardEuler' )

        nStages = 1

        CALL AllocateButcherTables

        a_IM(1,1) = 1.0_DP
        w_IM(1)   = 1.0_DP

      CASE ( 'IMEX_ARS_111' )

        nStages = 2

        CALL AllocateButcherTables

        ! --- Coefficients from Ascher et al. (1997) ---

        a_EX(2,1) = 1.0_DP
        w_EX(1)   = a_EX(2,1)

        a_IM(2,2) = 1.0_DP
        w_IM(2)   = a_IM(2,2)

      CASE ( 'IMEX_PDARS' )

        nStages = 3

        CALL AllocateButcherTables

        a_EX(2,1) = 1.0_DP
        a_EX(3,1) = 0.5_DP
        a_EX(3,2) = 0.5_DP

        w_EX(1)   = a_EX(3,1)
        w_EX(2)   = a_EX(3,2)

        a_IM(2,2) = 1.0_DP
        a_IM(3,2) = 0.5_DP
        a_IM(3,3) = 0.5_DP

        w_IM(2)   = a_IM(3,2)
        w_IM(3)   = a_IM(3,3)

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,'(A6,A,A)') &
          '', 'Unknown Time Stepping Scheme: ', TRIM( Scheme )
        WRITE(*,*)
        WRITE(*,'(A6,A)') &
          '', 'Available Options:'
        WRITE(*,*)
        WRITE(*,'(A6,A)') '', 'SSPRK1'
        WRITE(*,'(A6,A)') '', 'SSPRK2'
        WRITE(*,'(A6,A)') '', 'SSPRK3'
        WRITE(*,'(A6,A)') '', 'BackwardEuler'
        WRITE(*,'(A6,A)') '', 'IMEX_ARS_111'
        WRITE(*,'(A6,A)') '', 'IMEX_PDARS'
        WRITE(*,*)
        STOP

    END SELECT

    DO i = 1, nStages
      c_IM(i) = SUM( a_IM(i,1:i) )
      c_EX(i) = SUM( a_EX(i,1:i-1) )
    END DO

    WRITE(*,*)
    WRITE(*,'(A6,A)') '', 'Implicit Butcher Table:'
    WRITE(*,'(A6,A)') '', '-----------------------'
    DO i = 1, nStages
      WRITE(*,'(A6,5ES14.4E3)') '', c_IM(i), a_IM(i,1:nStages)
    END DO
    WRITE(*,'(A6,A14,4ES14.4E3)') '', '', w_IM(1:nStages)

    WRITE(*,*)
    WRITE(*,'(A6,A)') '', 'Explicit Butcher Table:'
    WRITE(*,'(A6,A)') '', '-----------------------'
    DO i = 1, nStages
      WRITE(*,'(A6,5ES14.4E3)') '', c_EX(i), a_EX(i,1:nStages)
    END DO
    WRITE(*,'(A6,A14,4ES14.4E3)') '', '', w_EX(1:nStages)

    CALL AllocateArray( U0 )
    CALL AllocateArray( Ui )
    CALL AllocateArray( M0 )
    CALL AllocateArray( Mi )

    ALLOCATE( StageData(nStages) )

    DO i = 1, nStages

      CALL AllocateArray( StageData(i) % dU_IM )
      CALL AllocateArray( StageData(i) % dU_EX )
      CALL AllocateArray( StageData(i) % dM_IM )
      CALL AllocateArray( StageData(i) % dM_EX )

      StageData(i) % OffGridFlux_U = Zero
      StageData(i) % OffGridFlux_M = Zero

    END DO

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: c_IM, w_IM, a_IM, c_EX, w_EX, a_EX )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA &
    !$ACC COPYIN( c_IM, w_IM, a_IM, c_EX, w_EX, a_EX )
#endif

  END SUBROUTINE Initialize_IMEX_RK


  SUBROUTINE Finalize_IMEX_RK

    INTEGER :: i

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: c_IM, w_IM, a_IM, c_EX, w_EX, a_EX )
#elif defined( THORNADO_OACC   )
    !$ACC EXIT DATA &
    !$ACC DELETE( c_IM, w_IM, a_IM, c_EX, w_EX, a_EX )
#endif

    DEALLOCATE( c_IM, w_IM, a_IM )
    DEALLOCATE( c_EX, w_EX, a_EX )

    CALL DeallocateArray( U0 )
    CALL DeallocateArray( Ui )
    CALL DeallocateArray( M0 )
    CALL DeallocateArray( Mi )

    DO i = 1, nStages

      CALL DeallocateArray( StageData(i) % dU_IM )
      CALL DeallocateArray( StageData(i) % dU_EX )
      CALL DeallocateArray( StageData(i) % dM_IM )
      CALL DeallocateArray( StageData(i) % dM_EX )

    END DO

    DEALLOCATE( StageData )

  END SUBROUTINE Finalize_IMEX_RK


  SUBROUTINE AllocateButcherTables

    ! --- Implicit Coefficients ---

    ALLOCATE( c_IM(nStages) )
    ALLOCATE( w_IM(nStages) )
    ALLOCATE( a_IM(nStages,nStages) )

    c_IM = Zero
    w_IM = Zero
    a_IM = Zero

    ! --- Explicit Coefficients ---

    ALLOCATE( c_EX(nStages) )
    ALLOCATE( w_EX(nStages) )
    ALLOCATE( a_EX(nStages,nStages) )

    c_EX = Zero
    w_EX = Zero
    a_EX = Zero

  END SUBROUTINE AllocateButcherTables


  SUBROUTINE AllocateArray5D( Array5D )

    REAL(DP), ALLOCATABLE, INTENT(inout) :: Array5D(:,:,:,:,:)

    CALL TimersStart( Timer_TimeStepper )

    ALLOCATE &
      ( Array5D(nDOFX, &
                iX_B1(1):iX_E1(1), &
                iX_B1(2):iX_E1(2), &
                iX_B1(3):iX_E1(3), &
                nCF) )

    Array5D = Zero

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: Array5D )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA &
    !$ACC COPYIN( Array5D )
#endif

    CALL TimersStop( Timer_TimeStepper )

  END SUBROUTINE AllocateArray5D


  SUBROUTINE DeallocateArray5D( Array5D )

    REAL(DP), ALLOCATABLE, INTENT(inout) :: Array5D(:,:,:,:,:)

    CALL TimersStart( Timer_TimeStepper )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: Array5D )
#elif defined( THORNADO_OACC   )
    !$ACC EXIT DATA &
    !$ACC DELETE( Array5D )
#endif

    DEALLOCATE( Array5D )

    CALL TimersStop( Timer_TimeStepper )

  END SUBROUTINE DeallocateArray5D


  SUBROUTINE AllocateArray7D( Array7D )

    REAL(DP), ALLOCATABLE, INTENT(inout) :: Array7D(:,:,:,:,:,:,:)

    CALL TimersStart( Timer_TimeStepper )

    ALLOCATE &
      ( Array7D(nDOFZ, &
                iZ_B1(1):iZ_E1(1), &
                iZ_B1(2):iZ_E1(2), &
                iZ_B1(3):iZ_E1(3), &
                iZ_B1(4):iZ_E1(4), &
                nCR,nSpecies) )

    Array7D = Zero

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: Array7D )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA &
    !$ACC COPYIN( Array7D )
#endif

    CALL TimersStop( Timer_TimeStepper )

  END SUBROUTINE AllocateArray7D


  SUBROUTINE DeallocateArray7D( Array7D )

    REAL(DP), ALLOCATABLE, INTENT(inout) :: Array7D(:,:,:,:,:,:,:)

    CALL TimersStart( Timer_TimeStepper )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: Array7D )
#elif defined( THORNADO_OACC   )
    !$ACC EXIT DATA &
    !$ACC DELETE( Array7D )
#endif

    DEALLOCATE( Array7D )

    CALL TimersStop( Timer_TimeStepper )

  END SUBROUTINE DeallocateArray7D


  SUBROUTINE CopyArray5D( U, alpha, V )

    REAL(DP), INTENT(inout) :: &
      U(1:nDOFX, &
        iX_B1(1):iX_E1(1), &
        iX_B1(2):iX_E1(2), &
        iX_B1(3):iX_E1(3), &
        1:nCF)
    REAL(DP), INTENT(in)    :: &
      alpha
    REAL(DP), INTENT(in)    :: &
      V(1:nDOFX, &
        iX_B1(1):iX_E1(1), &
        iX_B1(2):iX_E1(2), &
        iX_B1(3):iX_E1(3), &
        1:nCF)

    INTEGER :: iNodeX, iX1, iX2, iX3, iCF

    CALL TimersStart( Timer_TimeStepper )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( U, V, iX_B1, iX_E1 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iCF = 1, nCF
    DO iX3 = iX_B1(3), iX_E1(3)
    DO iX2 = iX_B1(2), iX_E1(2)
    DO iX1 = iX_B1(1), iX_E1(1)

      DO iNodeX = 1, nDOFX

        U(iNodeX,iX1,iX2,iX3,iCF) &
          = alpha * V(iNodeX,iX1,iX2,iX3,iCF)

      END DO

    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_TimeStepper )

  END SUBROUTINE CopyArray5D


  SUBROUTINE CopyArray7D( U, alpha, V )

    REAL(DP), INTENT(inout) :: &
      U(1:nDOFZ, &
        iZ_B1(1):iZ_E1(1), &
        iZ_B1(2):iZ_E1(2), &
        iZ_B1(3):iZ_E1(3), &
        iZ_B1(4):iZ_E1(4), &
        1:nCR,1:nSpecies)
    REAL(DP), INTENT(in)    :: &
      alpha
    REAL(DP), INTENT(in)    :: &
      V(1:nDOFZ, &
        iZ_B1(1):iZ_E1(1), &
        iZ_B1(2):iZ_E1(2), &
        iZ_B1(3):iZ_E1(3), &
        iZ_B1(4):iZ_E1(4), &
        1:nCR,1:nSpecies)

    INTEGER :: iNodeZ, iZ1, iZ2, iZ3, iZ4, iCR, iS

    CALL TimersStart( Timer_TimeStepper )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) &
    !$ACC PRESENT( U, V, iZ_B1, iZ_E1 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(7)
#endif
    DO iS  = 1, nSpecies
    DO iCR = 1, nCR
    DO iZ4 = iZ_B1(4), iZ_E1(4)
    DO iZ3 = iZ_B1(3), iZ_E1(3)
    DO iZ2 = iZ_B1(2), iZ_E1(2)
    DO iZ1 = iZ_B1(1), iZ_E1(1)

      DO iNodeZ = 1, nDOFZ

        U(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS) &
          = alpha * V(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS)

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_TimeStepper )

  END SUBROUTINE CopyArray7D


  SUBROUTINE AddToArray5D( alpha, U, beta, V )

    REAL(DP), INTENT(in)    :: &
      alpha
    REAL(DP), INTENT(inout) :: &
      U(1:nDOFX, &
        iX_B1(1):iX_E1(1), &
        iX_B1(2):iX_E1(2), &
        iX_B1(3):iX_E1(3), &
        1:nCF)
    REAL(DP), INTENT(in)    :: &
      beta
    REAL(DP), INTENT(in)    :: &
      V(1:nDOFX, &
        iX_B1(1):iX_E1(1), &
        iX_B1(2):iX_E1(2), &
        iX_B1(3):iX_E1(3), &
        1:nCF)

    INTEGER :: iNodeX, iX1, iX2, iX3, iCF

    CALL TimersStart( Timer_TimeStepper )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( U, V, iX_B1, iX_E1 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iCF = 1, nCF
    DO iX3 = iX_B1(3), iX_E1(3)
    DO iX2 = iX_B1(2), iX_E1(2)
    DO iX1 = iX_B1(1), iX_E1(1)

      DO iNodeX = 1, nDOFX

        U(iNodeX,iX1,iX2,iX3,iCF) &
          = alpha * U(iNodeX,iX1,iX2,iX3,iCF) &
            + beta * V(iNodeX,iX1,iX2,iX3,iCF)

      END DO

    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_TimeStepper )

  END SUBROUTINE AddToArray5D


  SUBROUTINE AddToArray7D( alpha, U, beta, V )

    REAL(DP), INTENT(in)    :: &
      alpha
    REAL(DP), INTENT(inout) :: &
      U(1:nDOFZ, &
        iZ_B1(1):iZ_E1(1), &
        iZ_B1(2):iZ_E1(2), &
        iZ_B1(3):iZ_E1(3), &
        iZ_B1(4):iZ_E1(4), &
        1:nCR,1:nSpecies)
    REAL(DP), INTENT(in)    :: &
      beta
    REAL(DP), INTENT(in)    :: &
      V(1:nDOFZ, &
        iZ_B1(1):iZ_E1(1), &
        iZ_B1(2):iZ_E1(2), &
        iZ_B1(3):iZ_E1(3), &
        iZ_B1(4):iZ_E1(4), &
        1:nCR,1:nSpecies)

    INTEGER :: iNodeZ, iZ1, iZ2, iZ3, iZ4, iCR, iS

    CALL TimersStart( Timer_TimeStepper )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) &
    !$ACC PRESENT( U, V, iZ_B1, iZ_E1 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(7)
#endif
    DO iS  = 1, nSpecies
    DO iCR = 1, nCR
    DO iZ4 = iZ_B1(4), iZ_E1(4)
    DO iZ3 = iZ_B1(3), iZ_E1(3)
    DO iZ2 = iZ_B1(2), iZ_E1(2)
    DO iZ1 = iZ_B1(1), iZ_E1(1)

      DO iNodeZ = 1, nDOFZ

        U(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS) &
          = alpha * U(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS) &
            + beta * V(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS)

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_TimeStepper )

  END SUBROUTINE AddToArray7D


END MODULE TwoMoment_TimeSteppingModule_OrderV
