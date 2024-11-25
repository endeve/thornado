MODULE TimeSteppingModule_CCSN

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
  USE RadiationFieldsModule, ONLY: &
    nCR, nSpecies
  USE Euler_SlopeLimiterModule_NonRelativistic_TABLE, ONLY: &
    ApplySlopeLimiter_Euler_NonRelativistic_TABLE
  USE Euler_PositivityLimiterModule_NonRelativistic_TABLE, ONLY: &
    ApplyPositivityLimiter_Euler_NonRelativistic_TABLE
  USE TwoMoment_PositivityLimiterModule_Old, ONLY: &
    ApplyPositivityLimiter_TwoMoment
  USE GravitySolutionModule_Newtonian_Poseidon, ONLY: &
    SolveGravity_Newtonian_Poseidon
  USE Euler_dgDiscretizationModule, ONLY: &
    ComputeIncrement_Euler_DG_Explicit
  USE TwoMoment_DiscretizationModule_Streaming, ONLY: &
    ComputeIncrement_TwoMoment_Explicit
  USE TwoMoment_DiscretizationModule_Collisions_Neutrinos, ONLY: &
    ComputeIncrement_TwoMoment_Implicit_New

  IMPLICIT NONE
  PRIVATE

  TYPE :: StageDataType
    REAL(DP), ALLOCATABLE :: dU_IM(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: dU_EX(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: dM_IM(:,:,:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: dM_EX(:,:,:,:,:,:,:)
  END type StageDataType

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

  PUBLIC :: InitializeTimeStepping
  PUBLIC :: FinalizeTimeStepping
  PUBLIC :: UpdateFields

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


  SUBROUTINE UpdateFields( dt, GE, GX, U, M )

    REAL(DP), INTENT(in)    :: &
      dt
    REAL(DP), INTENT(in)    :: &
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

    INTEGER :: iS, jS

    CALL CopyArray( U0, One, U )
    CALL CopyArray( M0, One, M )

    ! --- IMEX Stages ---

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

        IF( jS .EQ. iS - 1 )THEN

          ! --- Apply Limiters ---

          IF( EvolveEuler )THEN

            CALL ApplySlopeLimiter_Euler_NonRelativistic_TABLE &
                   ( iX_B0, iX_E0, iX_B1, iX_E1, GX, Ui, uDF )

            CALL ApplyPositivityLimiter_Euler_NonRelativistic_TABLE &
                   ( iX_B0, iX_E0, iX_B1, iX_E1, GX, Ui, uDF )

          END IF

          IF( EvolveTwoMoment )THEN

            CALL ApplyPositivityLimiter_TwoMoment &
                   ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, Mi )

          END IF

          ! --- Solve Gravity ---

          IF( EvolveEuler )THEN

            CALL SolveGravity_Newtonian_Poseidon &
                   ( iX_B0, iX_E0, iX_B1, iX_E1, GX, Ui(:,:,:,:,iCF_D) )

          END IF

        END IF

      END DO ! jS = 1, iS - 1

      IF( ANY( a_IM(:,iS) .NE. Zero ) .OR. ( w_IM(iS) .NE. Zero ) )THEN

        ! --- Implicit Solve ---

        CALL ComputeIncrement_TwoMoment_Implicit_New &
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
                 ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, Mi )

        END IF

      END IF

      IF( ANY( a_EX(:,iS) .NE. Zero ) .OR. ( w_EX(iS) .NE. Zero ) )THEN

        ! --- Explicit Solve ---

        IF( EvolveEuler )THEN

          CALL ComputeIncrement_Euler_DG_Explicit &
                 ( iX_B0, iX_E0, iX_B1, iX_E1, GX, Ui, uDF, &
                   StageData(iS) % dU_EX )

        END IF

        IF( EvolveTwoMoment )THEN

          CALL ComputeIncrement_TwoMoment_Explicit &
                 ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, Mi, &
                   StageData(iS) % dM_EX )

        END IF

      END IF

    END DO ! iS = 1, nStages

    ! --- IMEX Assembly ---

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

      ! --- Apply Limiters ---

      IF( EvolveEuler )THEN

        CALL ApplySlopeLimiter_Euler_NonRelativistic_TABLE &
               ( iX_B0, iX_E0, iX_B1, iX_E1, GX, Ui, uDF )

        CALL ApplyPositivityLimiter_Euler_NonRelativistic_TABLE &
               ( iX_B0, iX_E0, iX_B1, iX_E1, GX, Ui, uDF )

      END IF

      IF( EvolveTwoMoment )THEN

        CALL ApplyPositivityLimiter_TwoMoment &
               ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, Mi )

      END IF

    END IF

    ! --- Solve Gravity ---

    IF( EvolveEuler )THEN

      CALL SolveGravity_Newtonian_Poseidon &
             ( iX_B0, iX_E0, iX_B1, iX_E1, GX, Ui(:,:,:,:,iCF_D) )

    END IF

    CALL CopyArray( U, One, Ui )
    CALL CopyArray( M, One, Mi )

  END SUBROUTINE UpdateFields


  SUBROUTINE InitializeTimeStepping &
    ( Scheme_Option, EvolveEuler_Option, EvolveTwoMoment_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: Scheme_Option
    LOGICAL,          INTENT(in), OPTIONAL :: EvolveEuler_Option
    LOGICAL,          INTENT(in), OPTIONAL :: EvolveTwoMoment_Option

    CHARACTER(16) :: Scheme
    INTEGER       :: i

    Scheme = 'IMEX_PDARS'
    IF( PRESENT( Scheme_Option ) ) &
      Scheme = TRIM( Scheme_Option )

    EvolveEuler = .TRUE.
    IF( PRESENT( EvolveEuler_Option ) ) &
      EvolveEuler = EvolveEuler_Option

    EvolveTwoMoment = .TRUE.
    IF( PRESENT( EvolveTwoMoment_Option ) ) &
      EvolveTwoMoment = EvolveTwoMoment_Option

    WRITE(*,*)
    WRITE(*,'(A6,A22,A)' ) '', 'Time Stepping Scheme: ', TRIM( Scheme )
    WRITE(*,*)
    WRITE(*,'(A8,A18,L1)') '', 'Evolve Fluid: ', EvolveEuler
    WRITE(*,'(A8,A18,L1)') '', 'Evolve Neutrinos: ', EvolveTwoMoment

    SELECT CASE( TRIM( Scheme ) )

      CASE ( 'SSPRK2' )

        nStages = 2

        CALL AllocateButcherTables

        a_EX(1,1:2) = [ 0.0_DP, 0.0_DP ]
        a_EX(2,1:2) = [ 1.0_DP, 0.0_DP ]
        w_EX(1:2)   = [ 0.5_DP, 0.5_DP ]

       
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
        WRITE(*,'(A6,A)') '', 'SSPRK2, IMEX_PDARS'
        WRITE(*,*)
        STOP

    END SELECT

    DO i = 1, nStages
      c_IM(i) = SUM( a_IM(i,1:i) )
      c_EX(i) = SUM( a_EX(i,1:i-1) )
    END DO

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

    END DO

  END SUBROUTINE InitializeTimeStepping


  SUBROUTINE FinalizeTimeStepping

    INTEGER :: i

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

  END SUBROUTINE FinalizeTimeStepping


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

    ALLOCATE &
      ( Array5D(nDOFX, &
                iX_B1(1):iX_E1(1), &
                iX_B1(2):iX_E1(2), &
                iX_B1(3):iX_E1(3), &
                nCF) )

    Array5D = Zero

  END SUBROUTINE AllocateArray5D


  SUBROUTINE DeallocateArray5D( Array5D )

    REAL(DP), ALLOCATABLE, INTENT(inout) :: Array5D(:,:,:,:,:)

    DEALLOCATE( Array5D )

  END SUBROUTINE DeallocateArray5D


  SUBROUTINE AllocateArray7D( Array7D )

    REAL(DP), ALLOCATABLE, INTENT(inout) :: Array7D(:,:,:,:,:,:,:)

    ALLOCATE &
      ( Array7D(nDOFZ, &
                iZ_B1(1):iZ_E1(1), &
                iZ_B1(2):iZ_E1(2), &
                iZ_B1(3):iZ_E1(3), &
                iZ_B1(4):iZ_E1(4), &
                nCR, &
                nSpecies) )

    Array7D = Zero

  END SUBROUTINE AllocateArray7D


  SUBROUTINE DeallocateArray7D( Array7D )

    REAL(DP), ALLOCATABLE, INTENT(inout) :: Array7D(:,:,:,:,:,:,:)

    DEALLOCATE( Array7D )

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

    IF( .NOT. EvolveEuler ) RETURN

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

    IF( .NOT. EvolveTwoMoment ) RETURN

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

  END SUBROUTINE AddToArray7D


END MODULE TimeSteppingModule_CCSN
