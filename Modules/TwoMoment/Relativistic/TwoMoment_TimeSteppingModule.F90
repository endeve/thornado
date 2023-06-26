MODULE TwoMoment_TimeSteppingModule

  USE KindModule, ONLY: &
    DP, Zero
  USE ProgramHeaderModule, ONLY: &
    nDOFZ, nDOFX, nDOFE, &
    iZ_B0, iZ_B1, iZ_E0, iZ_E1
  USE GeometryFieldsModuleE, ONLY: &
    nGE
  USE GeometryFieldsModule, ONLY: &
    nGF
  USE FluidFieldsModule, ONLY: &
    nCF
  USE RadiationFieldsModule, ONLY: &
    nCR, nSpecies
!!$  USE TwoMoment_TroubledCellIndicatorModule, ONLY: &
!!$    DetectTroubledCells_TwoMoment
!!$  USE TwoMoment_SlopeLimiterModule_OrderV, ONLY: &
!!$    ApplySlopeLimiter_TwoMoment
  USE TwoMoment_PositivityLimiterModule, ONLY: &
    ApplyPositivityLimiter_TwoMoment
  USE TwoMoment_DiscretizationModule_Streaming, ONLY: &
    ComputeIncrement_TwoMoment_Explicit
!!$  USE TwoMoment_DiscretizationModule_Collisions_OrderV, ONLY: &
!!$    ComputeIncrement_TwoMoment_Implicit

  IMPLICIT NONE
  PRIVATE

  TYPE :: StageDataType
    REAL(DP), ALLOCATABLE :: dU_IM(:,:,:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: dU_EX(:,:,:,:,:,:,:)
  END TYPE StageDataType

  INTEGER                          :: nStages
  REAL(DP),            ALLOCATABLE :: c_IM(:), w_IM(:), a_IM(:,:)
  REAL(DP),            ALLOCATABLE :: c_EX(:), w_EX(:), a_EX(:,:)
  REAL(DP),            ALLOCATABLE :: U0(:,:,:,:,:,:,:)
  REAL(DP),            ALLOCATABLE :: Ui(:,:,:,:,:,:,:)
  TYPE(StageDataType), ALLOCATABLE :: StageData(:)

  PUBLIC :: Initialize_IMEX_RK
  PUBLIC :: Finalize_IMEX_RK
  PUBLIC :: Update_IMEX_RK

CONTAINS


  SUBROUTINE Update_IMEX_RK( dt, GE, GX, UF, UR )

    REAL(DP), INTENT(in) :: &
      dt
    REAL(DP), INTENT(in)    :: &
      GE(1:nDOFE,iZ_B1(1):iZ_E1(1),1:nGE)
    REAL(DP), INTENT(in)    :: &
      GX(1:nDOFX, &
         iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4), &
         1:nGF)
    REAL(DP), INTENT(inout) :: &
      UF(1:nDOFX, &
         iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4), &
         1:nCF)
    REAL(DP), INTENT(inout) :: &
      UR(1:nDOFZ, &
         iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
         iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4), &
         1:nCR,1:nSpecies)

    INTEGER :: iS, jS

    U0 = UR

    DO iS = 1, nStages

      Ui = U0

      DO jS = 1, iS - 1

        IF( a_IM(iS,jS) .NE. Zero )THEN

          Ui = Ui + dt * a_IM(iS,jS) * StageData(jS) % dU_IM

        END IF

        IF( a_EX(iS,jS) .NE. Zero )THEN

          Ui = Ui + dt * a_EX(iS,jS) * StageData(jS) % dU_EX

        END IF

        IF( jS == iS - 1 )THEN

!!$          CALL DetectTroubledCells_TwoMoment &
!!$                 ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, Ui )

          ! --- Apply Limiters ---

!!$          CALL ApplySLopeLimiter_TwoMoment &
!!$                 ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, UF, Ui )

          CALL ApplyPositivityLimiter_TwoMoment &
                 ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, UF, Ui )

        END IF

      END DO ! jS = 1, iS - 1

      IF( ANY( a_IM(:,iS) .NE. Zero ) .OR. ( w_IM(iS) .NE. Zero ) )THEN

        ! --- Implicit Solve ---

        PRINT*, "    IMPLICIT: ", iS

!!$        CALL ComputeIncrement_TwoMoment_Implicit &
!!$               ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, dt * a_IM(iS,iS), &
!!$                 GE, GX, UF, Ui, StageData(iS) % dU_IM )

        Ui = Ui + dt * a_IM(iS,iS) * StageData(iS) % dU_IM

      END IF

      IF( ANY( a_EX(:,iS) .NE. Zero ) .OR. ( w_EX(iS) .NE. Zero ) )THEN

        ! --- Explicit Solve ---

        PRINT*, "    EXPLICIT: ", iS

        CALL ComputeIncrement_TwoMoment_Explicit &
               ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, &
                 GE, GX, UF, Ui, StageData(iS) % dU_EX )

      END IF

    END DO ! iS = 1, nStages

    ! --- Assembly Step ---

    IF( ANY( a_IM(nStages,:) .NE. w_IM(:) ) .OR. &
        ANY( a_EX(nStages,:) .NE. w_EX(:) ) )THEN

      PRINT*, "    ASSEMBLY:"

      Ui = U0

      DO iS = 1, nStages

        IF( w_IM(iS) .NE. Zero )THEN

          Ui = Ui + dt * w_IM(iS) * StageData(iS) % dU_IM

        END IF

        IF( w_EX(iS) .NE. Zero )THEN

          Ui = Ui + dt * w_EX(iS) * StageData(iS) % dU_EX

        END IF

      END DO

!!$      CALL DetectTroubledCells_TwoMoment &
!!$             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, Ui )

!!$      CALL ApplySlopeLimiter_TwoMoment &
!!$             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, UF, Ui )

!!$      CALL ApplyPositivityLimiter_TwoMoment &
!!$             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, UF, Ui )

    END IF

    UR = Ui

  END SUBROUTINE Update_IMEX_RK


  SUBROUTINE Initialize_IMEX_RK( Scheme )

    CHARACTER(LEN=*), INTENT(in) :: Scheme

    INTEGER :: i

    WRITE(*,*)
    WRITE(*,'(A6,A,A)') '', 'IMEX-RK Scheme: ', TRIM( Scheme )

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

    CALL AllocateArray7D( U0 )
    CALL AllocateArray7D( Ui )

    ALLOCATE( StageData(nStages) )

    DO i = 1, nStages

      CALL AllocateArray7D( StageData(i) % dU_IM )
      CALL AllocateArray7D( StageData(i) % dU_EX )

    END DO

  END SUBROUTINE Initialize_IMEX_RK


  SUBROUTINE Finalize_IMEX_RK

    INTEGER :: i

    DEALLOCATE( c_IM, w_IM, a_IM )
    DEALLOCATE( c_EX, w_EX, a_EX )

    CALL DeallocateArray7D( U0 )
    CALL DeallocateArray7D( Ui )

    DO i = 1, nStages

      CALL DeallocateArray7D( StageData(i) % dU_IM )
      CALL DeallocateArray7D( StageData(i) % dU_EX )

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

    ALLOCATE &
      ( Array7D(nDOFZ, &
                iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
                iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4), &
                nCR,nSpecies) )

    Array7D = Zero

  END SUBROUTINE AllocateArray7D


  SUBROUTINE DeallocateArray7D( Array7D )

    REAL(DP), ALLOCATABLE, INTENT(inout) :: Array7D(:,:,:,:,:,:,:)

    DEALLOCATE( Array7D )

  END SUBROUTINE DeallocateArray7D


END MODULE TwoMoment_TimeSteppingModule
