MODULE TimeSteppingModule_SSPRK

  USE KindModule, ONLY: &
    DP, &
    Zero, &
    One
  USE ProgramHeaderModule, ONLY: &
    iX_B0, &
    iX_B1, &
    iX_E0, &
    iX_E1, &
    nDOFX
  USE MagnetofluidFieldsModule, ONLY: &
    nCM
  USE MHD_DiscretizationModule_Relativistic, ONLY: &
    OffGridFlux_MHD
  !USE Euler_TallyModule_Relativistic, ONLY: &
    !IncrementOffGridTally_MHD_Relativistic

  IMPLICIT NONE
  PRIVATE

  INTEGER :: nStages_SSPRK
  REAL(DP), DIMENSION(:),   ALLOCATABLE :: c_SSPRK
  REAL(DP), DIMENSION(:),   ALLOCATABLE :: w_SSPRK
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: a_SSPRK

  REAL(DP), DIMENSION(:,:,:,:,:),   ALLOCATABLE :: U_SSPRK
  REAL(DP), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: D_SSPRK

  PUBLIC :: InitializeMagnetofluid_SSPRK
  PUBLIC :: UpdateMagnetofluid_SSPRK
  PUBLIC :: FinalizeMagnetofluid_SSPRK

  INTERFACE
    SUBROUTINE MagnetofluidIncrement &
      ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dU, &
        SuppressBC_Option, UseXCFC_Option )
      USE KindModule, ONLY: DP
      INTEGER, INTENT(in)     :: &
        iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
      REAL(DP), INTENT(in)    :: &
        G (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
      REAL(DP), INTENT(inout) :: &
        U (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
      REAL(DP), INTENT(inout) :: &
        D (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
      REAL(DP), INTENT(out)   :: &
        dU(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
      LOGICAL, INTENT(in), OPTIONAL :: &
        SuppressBC_Option
      LOGICAL, INTENT(in), OPTIONAL :: &
        UseXCFC_Option
    END SUBROUTINE MagnetofluidIncrement
  END INTERFACE

CONTAINS


  SUBROUTINE InitializeMagnetofluid_SSPRK( nStages )

    INTEGER, INTENT(in) :: nStages

    INTEGER :: i

    nStages_SSPRK = nStages

    CALL InitializeSSPRK( nStages )

    WRITE(*,*)
    WRITE(*,*)
    WRITE(*,'(A)') '    INFO: TimeSteppingModule_SSPRK'
    WRITE(*,'(A)') '    ------------------------------'

    WRITE(*,*)
    WRITE(*,'(A5,A,I1)') '', 'SSP RK Scheme: ', nStages

    WRITE(*,*)
    WRITE(*,'(A5,A)') '', 'Butcher Table:'
    WRITE(*,'(A5,A)') '', '--------------'
    DO i = 1, nStages
      WRITE(*,'(A5,4ES14.4E3)') '', c_SSPRK(i), a_SSPRK(i,1:nStages)
    END DO
    WRITE(*,'(A5,A14,3ES14.4E3)') '', '', w_SSPRK(1:nStages)

    ALLOCATE( U_SSPRK &
                (1:nDOFX, &
                 iX_B1(1):iX_E1(1), &
                 iX_B1(2):iX_E1(2), &
                 iX_B1(3):iX_E1(3), &
                 1:nCM) )

    ALLOCATE( D_SSPRK &
                (1:nDOFX, &
                 iX_B1(1):iX_E1(1), &
                 iX_B1(2):iX_E1(2), &
                 iX_B1(3):iX_E1(3), &
                 1:nCM,1:nStages) )

  END SUBROUTINE InitializeMagnetofluid_SSPRK


  SUBROUTINE FinalizeMagnetofluid_SSPRK

    DEALLOCATE( a_SSPRK, c_SSPRK, w_SSPRK )

    DEALLOCATE( U_SSPRK, D_SSPRK )

  END SUBROUTINE FinalizeMagnetofluid_SSPRK


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


  SUBROUTINE UpdateMagnetofluid_SSPRK( t, dt, G, U, D, ComputeIncrement_Magnetofluid )

    REAL(DP), INTENT(in) :: &
      t, dt
    REAL(DP), INTENT(in) :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: &
      D(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    PROCEDURE (MagnetofluidIncrement) :: &
      ComputeIncrement_Magnetofluid
    LOGICAL :: DEBUG = .FALSE.

    INTEGER :: iNX, iX1, iX2, iX3, iCM
    INTEGER :: iS, jS

    REAL(DP) :: dM_OffGrid_MHD(nCM)

    dM_OffGrid_MHD = Zero

    DO iS = 1, nStages_SSPRK

      DO iCM = 1, nCM
      DO iX3 = iX_B1(3), iX_E1(3)
      DO iX2 = iX_B1(2), iX_E1(2)
      DO iX1 = iX_B1(1), iX_E1(1)
      DO iNX = 1, nDOFX

        U_SSPRK(iNX,iX1,iX2,iX3,iCM) = U(iNX,iX1,iX2,iX3,iCM)

      END DO
      END DO
      END DO
      END DO
      END DO

      DO jS = 1, iS - 1

        IF( a_SSPRK(iS,jS) .NE. Zero )THEN

          CALL AddIncrement_Magnetofluid &
                 ( jS, One, U_SSPRK, dt * a_SSPRK(iS,jS), D_SSPRK )

        END IF

      END DO

      IF( ANY( a_SSPRK(:,iS) .NE. Zero ) &
          .OR. ( w_SSPRK(iS) .NE. Zero ) )THEN

        CALL ComputeIncrement_Magnetofluid &
               ( iX_B0, iX_E0, iX_B1, iX_E1, &
                 G, U_SSPRK, D, D_SSPRK(:,:,:,:,:,iS) )

        dM_OffGrid_MHD &
          = dM_OffGrid_MHD &
              + dt * w_SSPRK(iS) * OffGridFlux_MHD

      END IF

    END DO

    DO iS = 1, nStages_SSPRK

      IF( w_SSPRK(iS) .NE. Zero )THEN

        CALL AddIncrement_Magnetofluid &
               ( iS, One, U, dt * w_SSPRK(iS), D_SSPRK )

      END IF

    END DO

    !CALL IncrementOffGridTally_MHD_Relativistic( dM_OffGrid_MHD )

  END SUBROUTINE UpdateMagnetofluid_SSPRK


  SUBROUTINE AddIncrement_Magnetofluid( iS, alpha, U, beta, D )

    INTEGER,  INTENT(in)    :: &
      iS
    REAL(DP), INTENT(in)    :: &
      alpha, beta
    REAL(DP), INTENT(inout) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(in)    :: &
      D(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:,1:)

    INTEGER :: iNX, iX1, iX2, iX3, iCM

    DO iCM = 1, nCM
    DO iX3 = iX_B1(3), iX_E1(3)
    DO iX2 = iX_B1(2), iX_E1(2)
    DO iX1 = iX_B1(1), iX_E1(1)
    DO iNX = 1, nDOFX

      U(iNX,iX1,iX2,iX3,iCM) &
        = alpha * U(iNX,iX1,iX2,iX3,iCM) &
            + beta * D(iNX,iX1,iX2,iX3,iCM,iS)

    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE AddIncrement_Magnetofluid


END MODULE TimeSteppingModule_SSPRK
