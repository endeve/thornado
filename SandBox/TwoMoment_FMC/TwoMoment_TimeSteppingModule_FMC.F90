MODULE TwoMoment_TimeSteppingModule_FMC

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
    nCF, iCF_D, uDF, &
    nPF
  USE TwoMoment_FieldsModule_FMC, ONLY: &
    nCM, nSpecies
  USE TwoMoment_TimersModule, ONLY: &
    TimersStart, &
    TimersStop, &
    Timer_IMEX, &
    Timer_TimeStepper, &
    Timer_Euler, &
    Timer_Poisson
  USE TwoMoment_PositivityLimiterModule_FMC, ONLY: &
    ApplyPositivityLimiter_TwoMoment
  USE TwoMoment_DiscretizationModule_Streaming_FMC, ONLY: &
    ComputeIncrement_TwoMoment_Explicit

  IMPLICIT NONE
  PRIVATE

  TYPE :: StageDataType
    REAL(DP)              :: OffGridFlux_M(2*nCM)
    REAL(DP), ALLOCATABLE :: dM_IM(:,:,:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: dM_EX(:,:,:,:,:,:,:)
  END TYPE StageDataType

  LOGICAL                          :: EvolveEuler
  LOGICAL                          :: EvolveTwoMoment
  INTEGER                          :: nStages
  REAL(DP),            ALLOCATABLE :: c_IM(:), w_IM(:), a_IM(:,:)
  REAL(DP),            ALLOCATABLE :: c_EX(:), w_EX(:), a_EX(:,:)
  REAL(DP),            ALLOCATABLE :: M0(:,:,:,:,:,:,:)
  REAL(DP),            ALLOCATABLE :: Mi(:,:,:,:,:,:,:)
  TYPE(StageDataType), ALLOCATABLE :: StageData(:)

  PUBLIC :: Initialize_IMEX_RK
  PUBLIC :: Finalize_IMEX_RK
  PUBLIC :: Update_IMEX_RK

  INTERFACE
    SUBROUTINE ImplicitIncrement &
      ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, dt, GE, GX, U_M, dU_M )
      USE KindModule           , ONLY: DP
      USE ProgramHeaderModule  , ONLY: nDOFX, nDOFE, nDOFZ
      USE GeometryFieldsModuleE, ONLY: nGE
      USE GeometryFieldsModule , ONLY: nGF
      USE FluidFieldsModule    , ONLY: nCF
      USE TwoMoment_FieldsModule_FMC, ONLY: nCM, nSpecies
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
        U_M (1:nDOFZ, &
             iZ_B1(1):iZ_E1(1), &
             iZ_B1(2):iZ_E1(2), &
             iZ_B1(3):iZ_E1(3), &
             iZ_B1(4):iZ_E1(4), &
             1:nCM, &
             1:nSpecies)
      REAL(DP), INTENT(out) :: &
        dU_M(1:nDOFZ, &
             iZ_B1(1):iZ_E1(1), &
             iZ_B1(2):iZ_E1(2), &
             iZ_B1(3):iZ_E1(3), &
             iZ_B1(4):iZ_E1(4), &
             1:nCM, &
             1:nSpecies)
    END SUBROUTINE ImplicitIncrement
  END INTERFACE

  INTERFACE AllocateArray
    MODULE PROCEDURE AllocateArray7D
  END INTERFACE AllocateArray

  INTERFACE DeallocateArray
    MODULE PROCEDURE DeallocateArray7D
  END INTERFACE DeallocateArray

  INTERFACE CopyArray
    MODULE PROCEDURE CopyArray7D
  END INTERFACE CopyArray

  INTERFACE AddToArray
    MODULE PROCEDURE AddToArray7D
  END INTERFACE AddToArray

CONTAINS

  SUBROUTINE Update_IMEX_RK &
    ( dt, GE, GX, U, M, ComputeIncrement_TwoMoment_Implicit)

    !--- Input/Output variable ---
    REAL(DP), INTENT(in) :: dt
    REAL(DP), INTENT(in) :: GE(1:nDOFE, &
                               iE_B1:iE_E1, &
                               1:nGE)
    REAL(DP), INTENT(inout) :: GX(1:nDOFX, &
                                  iX_B1(1):iX_E1(1), &
                                  iX_B1(2):iX_E1(2), &
                                  iX_B1(3):iX_E1(3), &
                                  1:nGF)
    REAL(DP), INTENT(inout) :: U(1:nDOFX, &
                                 iX_B1(1):iX_E1(1), &
                                 iX_B1(2):iX_E1(2), &
                                 iX_B1(3):iX_E1(3), &
                                 1:nPF)
    REAL(DP), INTENT(inout) :: M(1:nDOFZ, &
                                 iZ_B1(1):iZ_E1(1), &
                                 iZ_B1(2):iZ_E1(2), &
                                 iZ_B1(3):iZ_E1(3), &
                                 iZ_B1(4):iZ_E1(4), &
                                 1:nCM, 1:nSpecies)
    PROCEDURE(ImplicitIncrement) :: ComputeIncrement_TwoMoment_Implicit

    ! --- Local variables ---
    INTEGER :: iS, jS

    Write(*,*)
    print *,'Update_IMEX_RK'
    
    CALL CopyArray( M0, One, M)

    DO iS = 1, nStages

      CALL CopyArray(Mi, One, M0)

      DO jS = 1, iS - 1

        IF( a_IM(iS,jS) .NE. Zero )THEN

          CALL AddToArray( One, Mi, dt * a_IM(iS,jS), StageData(jS) % dM_IM )

        END IF

        IF (a_EX(iS,jS) .NE. Zero)THEN

          CALL AddToArray( One, Mi, dt * a_EX(iS,jS), StageData(jS) % dM_EX )

        END IF

        IF( jS==iS - 1)THEN

          ! --- Apply Limiters ---
          CALL ApplyPositivityLimiter_TwoMoment &
                 ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U, Mi )

        END IF

      END DO !jS = 1, iS -1

      ! --- Implicit Solve ---
      IF( ANY( a_IM(:,iS) .NE. Zero ) .OR. ( w_IM(iS) .NE. Zero ) )THEN

        CALL ComputeIncrement_TwoMoment_Implicit &
          ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, dt * a_IM(iS,iS), GE, GX, &
            Mi, StageData(iS) % dM_IM)

        CALL AddToArray ( One, Mi, dt * a_IM(iS,iS), StageData(iS) % dM_IM )

        ! --- Apply Limiters ---
        CALL ApplyPositivityLimiter_TwoMoment &
               ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U, Mi )

      END IF

      ! --- Explicit Solve ---

      CALL ComputeIncrement_TwoMoment_Explicit &
        ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, &
          U, Mi, StageData(iS) % dM_EX )

      ! StageData(iS) % OffGridFlux_M = OffGridFlux_TwoMoment

    END DO !iS = 1, nStages

    ! --- Assembly Step ---

    IF( ANY( a_IM(nStages,:) .NE. w_IM(:) ) .OR. &
        ANY( a_EX(nStages,:) .NE. w_EX(:) ) )THEN

      CALL CopyArray( Mi, One, M0 )

      DO iS = 1, nStages

        IF( w_IM(iS) .NE. Zero )THEN

          CALL AddToArray( One, Mi, dt * w_IM(iS), StageData(iS) % dM_IM )

        END IF

        IF( w_EX(iS) .NE. Zero )THEN

          CALL AddToArray( One, Mi, dt * w_EX(iS), StageData(iS) % dM_EX )

        END IF

      END DO

      ! --- Apply Limiters ---
      CALL ApplyPositivityLimiter_TwoMoment &
             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U, Mi )

      CALL CopyArray( M, One, Mi )

      ! --- Offgrid? ---

    END IF

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

    CALL AllocateArray( M0 )
    CALL AllocateArray( Mi )

    ALLOCATE( StageData(nStages) )

    DO i = 1, nStages

      CALL AllocateArray( StageData(i) % dM_IM )
      CALL AllocateArray( StageData(i) % dM_EX )

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

    CALL DeallocateArray( M0 )
    CALL DeallocateArray( Mi )

    DO i = 1, nStages

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

  SUBROUTINE AllocateArray7D( Array7D )

    REAL(DP), ALLOCATABLE, INTENT(inout) :: Array7D(:,:,:,:,:,:,:)

    CALL TimersStart( Timer_TimeStepper )

    ALLOCATE &
      ( Array7D(nDOFZ, &
                iZ_B1(1):iZ_E1(1), &
                iZ_B1(2):iZ_E1(2), &
                iZ_B1(3):iZ_E1(3), &
                iZ_B1(4):iZ_E1(4), &
                nCM,nSpecies) )

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

  SUBROUTINE CopyArray7D( U, alpha, V )

    REAL(DP), INTENT(inout) :: &
      U(1:nDOFZ, &
        iZ_B1(1):iZ_E1(1), &
        iZ_B1(2):iZ_E1(2), &
        iZ_B1(3):iZ_E1(3), &
        iZ_B1(4):iZ_E1(4), &
        1:nCM,1:nSpecies)
    REAL(DP), INTENT(in)    :: &
      alpha
    REAL(DP), INTENT(in)    :: &
      V(1:nDOFZ, &
        iZ_B1(1):iZ_E1(1), &
        iZ_B1(2):iZ_E1(2), &
        iZ_B1(3):iZ_E1(3), &
        iZ_B1(4):iZ_E1(4), &
        1:nCM,1:nSpecies)

    INTEGER :: iNodeZ, iZ1, iZ2, iZ3, iZ4, iCM, iS

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
    DO iCM = 1, nCM
    DO iZ4 = iZ_B1(4), iZ_E1(4)
    DO iZ3 = iZ_B1(3), iZ_E1(3)
    DO iZ2 = iZ_B1(2), iZ_E1(2)
    DO iZ1 = iZ_B1(1), iZ_E1(1)

      DO iNodeZ = 1, nDOFZ

        U(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM,iS) &
          = alpha * V(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM,iS)

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_TimeStepper )

  END SUBROUTINE CopyArray7D

  SUBROUTINE AddToArray7D( alpha, U, beta, V )

    REAL(DP), INTENT(in)    :: &
      alpha
    REAL(DP), INTENT(inout) :: &
      U(1:nDOFZ, &
        iZ_B1(1):iZ_E1(1), &
        iZ_B1(2):iZ_E1(2), &
        iZ_B1(3):iZ_E1(3), &
        iZ_B1(4):iZ_E1(4), &
        1:nCM,1:nSpecies)
    REAL(DP), INTENT(in)    :: &
      beta
    REAL(DP), INTENT(in)    :: &
      V(1:nDOFZ, &
        iZ_B1(1):iZ_E1(1), &
        iZ_B1(2):iZ_E1(2), &
        iZ_B1(3):iZ_E1(3), &
        iZ_B1(4):iZ_E1(4), &
        1:nCM,1:nSpecies)

    INTEGER :: iNodeZ, iZ1, iZ2, iZ3, iZ4, iCM, iS

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
    DO iCM = 1, nCM
    DO iZ4 = iZ_B1(4), iZ_E1(4)
    DO iZ3 = iZ_B1(3), iZ_E1(3)
    DO iZ2 = iZ_B1(2), iZ_E1(2)
    DO iZ1 = iZ_B1(1), iZ_E1(1)

      DO iNodeZ = 1, nDOFZ

        U(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM,iS) &
          = alpha * U(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM,iS) &
            + beta * V(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM,iS)

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_TimeStepper )

  END SUBROUTINE AddToArray7D

END MODULE TwoMoment_TimeSteppingModule_FMC