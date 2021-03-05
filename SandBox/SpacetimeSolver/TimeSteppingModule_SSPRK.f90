MODULE TimeSteppingModule_SSPRK

  USE KindModule, ONLY: &
    DP, Zero, One
  USE ProgramHeaderModule, ONLY: &
    iX_B0, iX_B1, iX_E0, iX_E1, &
    nDOFX
  USE ScalarFieldsModule, ONLY: &
    nSF, iSF_u, iSF_v

  IMPLICIT NONE
  PRIVATE

  INTEGER :: nStages_SSPRK
  REAL(DP), DIMENSION(:),   ALLOCATABLE :: c_SSPRK
  REAL(DP), DIMENSION(:),   ALLOCATABLE :: w_SSPRK
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: a_SSPRK

  REAL(DP), DIMENSION(:,:,:,:,:),   ALLOCATABLE :: U_SSPRK
  REAL(DP), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: D_SSPRK

  LOGICAL :: Verbose

  PUBLIC :: InitializeScalarWave_SSPRK
  PUBLIC :: UpdateScalarWave_SSPRK
  PUBLIC :: FinalizeScalarWave_SSPRK

  INTERFACE
    SUBROUTINE ScalarWaveIncrement &
      ( iX_B0, iX_E0, iX_B1, iX_E1, U, dU )
      USE KindModule, ONLY: DP
      INTEGER,  INTENT(in)           :: &
        iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
      REAL(DP), INTENT(inout)        :: &
        U (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
      REAL(DP), INTENT(out)          :: &
        dU(:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    END SUBROUTINE ScalarWaveIncrement
  END INTERFACE

CONTAINS


  SUBROUTINE InitializeScalarWave_SSPRK( nStages, Verbose_Option )

    INTEGER, INTENT(in)           :: nStages
    LOGICAL, INTENT(in), OPTIONAL :: Verbose_Option

    INTEGER :: i

    IF( PRESENT( Verbose_Option ) )THEN
      Verbose = Verbose_Option
    ELSE
       Verbose = .TRUE.
    END IF

    nStages_SSPRK = nStages

    CALL InitializeSSPRK( nStages )

    IF( Verbose )THEN
      WRITE(*,*)
      WRITE(*,'(A5,A,I1)') '', 'SSP RK Scheme: ', nStages

      WRITE(*,*)
      WRITE(*,'(A5,A)') '', 'Butcher Table:'
      WRITE(*,'(A5,A)') '', '--------------'
      DO i = 1, nStages
        WRITE(*,'(A5,4ES14.4E3)') '', c_SSPRK(i), a_SSPRK(i,1:nStages)
      END DO
      WRITE(*,'(A5,A14,3ES14.4E3)') '', '', w_SSPRK(1:nStages)
      WRITE(*,*)
    END IF

    ALLOCATE( U_SSPRK &
                (1:nDOFX, &
                 iX_B1(1):iX_E1(1), &
                 iX_B1(2):iX_E1(2), &
                 iX_B1(3):iX_E1(3), &
                 1:nSF) )

    ALLOCATE( D_SSPRK &
                (1:nDOFX, &
                 iX_B1(1):iX_E1(1), &
                 iX_B1(2):iX_E1(2), &
                 iX_B1(3):iX_E1(3), &
                 1:nSF,1:nStages) )

  END SUBROUTINE InitializeScalarWave_SSPRK


  SUBROUTINE FinalizeScalarWave_SSPRK

    DEALLOCATE( a_SSPRK, c_SSPRK, w_SSPRK )

    DEALLOCATE( U_SSPRK, D_SSPRK )

  END SUBROUTINE FinalizeScalarWave_SSPRK


  SUBROUTINE InitializeSSPRK( nStages )

    INTEGER, INTENT(in) :: nStages

    INTEGER :: iS

    CALL AllocateButcherTables_SSPRK( nStages )

    SELECT CASE ( nStages )
      CASE ( 1 )

        a_SSPRK(1,1) = 0.0_DP
        w_SSPRK(1)   = 1.0_DP

      CASE ( 2 )

        a_SSPRK(1,1:2) = [ 0.0_DP, 0.0_DP ]
        a_SSPRK(2,1:2) = [ 1.0_DP, 0.0_DP ]
        w_SSPRK(1:2)   = [ 0.5_DP, 0.5_DP ]

      CASE ( 3 )

        a_SSPRK(1,1:3) = [ 0.00_DP, 0.00_DP, 0.00_DP ]
        a_SSPRK(2,1:3) = [ 1.00_DP, 0.00_DP, 0.00_DP ]
        a_SSPRK(3,1:3) = [ 0.25_DP, 0.25_DP, 0.00_DP ]
        w_SSPRK(1:3)   = [ 1.0_DP / 6.0_DP, &
                           1.0_DP / 6.0_DP, &
                           2.0_DP / 3.0_DP ]

    END SELECT

    DO iS = 1, nStages
      c_SSPRK(iS) = SUM( a_SSPRK(iS,1:iS-1) )
    END DO

  END SUBROUTINE InitializeSSPRK


  SUBROUTINE AllocateButcherTables_SSPRK( nStages )

    INTEGER, INTENT(in) :: nStages

    ALLOCATE( a_SSPRK(nStages,nStages) )
    ALLOCATE( c_SSPRK(nStages) )
    ALLOCATE( w_SSPRK(nStages) )

    a_SSPRK = Zero
    c_SSPRK = Zero
    w_SSPRK = Zero

  END SUBROUTINE AllocateButcherTables_SSPRK


  SUBROUTINE UpdateScalarWave_SSPRK &
    ( t, dt, U, ComputeIncrement_ScalarWave )

    REAL(DP), INTENT(in) :: &
      t, dt
    REAL(DP), INTENT(inout) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    PROCEDURE(ScalarWaveIncrement) :: &
      ComputeIncrement_ScalarWave

    INTEGER :: iS, jS

    U_SSPRK = Zero ! --- State
    D_SSPRK = Zero ! --- Increment

    DO iS = 1, nStages_SSPRK

      U_SSPRK = U

      DO jS = 1, iS - 1

        IF( a_SSPRK(iS,jS) .NE. Zero )THEN

          CALL AddIncrement_ScalarWave &
                 ( One, U_SSPRK, dt * a_SSPRK(iS,jS), D_SSPRK(:,:,:,:,:,jS) )
        END IF

      END DO

      IF( ANY( a_SSPRK(:,iS) .NE. Zero ) &
          .OR. ( w_SSPRK(iS) .NE. Zero ) )THEN

        CALL ComputeIncrement_ScalarWave &
               ( iX_B0, iX_E0, iX_B1, iX_E1, U_SSPRK, D_SSPRK(:,:,:,:,:,iS) )

      END IF

    END DO

    DO iS = 1, nStages_SSPRK

      IF( w_SSPRK(iS) .NE. Zero )THEN

        CALL AddIncrement_ScalarWave &
               ( One, U, dt * w_SSPRK(iS), D_SSPRK(:,:,:,:,:,iS) )

      END IF

    END DO

  END SUBROUTINE UpdateScalarWave_SSPRK


  SUBROUTINE AddIncrement_ScalarWave( alpha, U, beta, D )

    REAL(DP), INTENT(in)    :: &
      alpha, beta
    REAL(DP), INTENT(inout) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(in)    :: &
      D(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER :: iSF, iX1, iX2, iX3

    DO iSF = 1, nSF
      DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
          DO iX1 = iX_B0(1), iX_E0(1)

            U(:,iX1,iX2,iX3,iSF) &
              = alpha * U(:,iX1,iX2,iX3,iSF) &
                  + beta * D(:,iX1,iX2,iX3,iSF)

          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE AddIncrement_ScalarWave


END MODULE TimeSteppingModule_SSPRK
