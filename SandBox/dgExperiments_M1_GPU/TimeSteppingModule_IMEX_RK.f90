MODULE TimeSteppingModule_IMEX_RK

  USE KindModule, ONLY: &
    DP, Zero
  USE ProgramHeaderModule, ONLY: &
    iZ_B0, iZ_B1, iZ_E0, iZ_E1, nDOF
  USE RadiationFieldsModule, ONLY: &
    nCR, nSpecies

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  INTEGER               :: nStages
  REAL(DP)              :: alpha
  REAL(DP), ALLOCATABLE :: c_IM(:), w_IM(:), a_IM(:,:)
  REAL(DP), ALLOCATABLE :: c_EX(:), w_EX(:), a_EX(:,:)
  REAL(DP), ALLOCATABLE :: U_IMEX(:), dU_IM(:,:), dU_EX(:,:)
  REAL(DP), ALLOCATABLE :: U0(:,:,:,:,:,:,:)
  REAL(DP), ALLOCATABLE :: dU(:,:,:,:,:,:,:)

  PUBLIC :: Initialize_IMEX_RK
  PUBLIC :: Finalize_IMEX_RK
  PUBLIC :: Update_IMEX_RK

  INTERFACE
    SUBROUTINE IncrementExplicit &
      ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U, dU )
      USE KindModule, ONLY: DP
      INTEGER, INTENT(in)     :: &
        iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
      REAL(DP), INTENT(in)    :: &
        GE(1:,iZ_B1(1):,1:)
      REAL(DP), INTENT(in)    :: &
        GX(1:,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:)
      REAL(DP), INTENT(inout) :: &
        U (1:,iZ_B1(1):,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:,1:)
      REAL(DP), INTENT(inout) :: &
        dU(1:,iZ_B0(1):,iZ_B0(2):,iZ_B0(3):,iZ_B0(4):,1:,1:)
    END SUBROUTINE IncrementExplicit
  END INTERFACE

  INTERFACE
    SUBROUTINE IncrementImplicit &
      ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, dt, GE, GX, U, dU )
      USE KindModule, ONLY: DP
      INTEGER, INTENT(in)     :: &
        iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
      REAL(DP), INTENT(in)    :: &
        dt
      REAL(DP), INTENT(in)    :: &
        GE(1:,iZ_B1(1):,1:)
      REAL(DP), INTENT(in)    :: &
        GX(1:,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:)
      REAL(DP), INTENT(inout) :: &
        U (1:,iZ_B1(1):,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:,1:)
      REAL(DP), INTENT(inout) :: &
        dU(1:,iZ_B0(1):,iZ_B0(2):,iZ_B0(3):,iZ_B0(4):,1:,1:)
    END SUBROUTINE IncrementImplicit
  END INTERFACE

CONTAINS


  SUBROUTINE Initialize_IMEX_RK( Scheme )

    CHARACTER(LEN=*), INTENT(in) :: Scheme

    INTEGER :: i, nDOF_IMEX

    WRITE(*,*)
    WRITE(*,'(A6,A,A)') '', 'IMEX-RK Scheme: ', TRIM( Scheme )

    SELECT CASE ( TRIM( Scheme ) )

      CASE ( 'SSPRK1' )

        nStages = 1
        CALL AllocateButcherTables( nStages )

        a_EX(1,1) = 0.0_DP
        w_EX(1)   = 1.0_DP

      CASE ( 'SSPRK2' )

        nStages = 2
        CALL AllocateButcherTables( nStages )

        a_EX(1,1:2) = [ 0.0_DP, 0.0_DP ]
        a_EX(2,1:2) = [ 1.0_DP, 0.0_DP ]
        w_EX(1:2)   = [ 0.5_DP, 0.5_DP ]

      CASE ( 'SSPRK3' )

        nStages = 3
        CALL AllocateButcherTables( nStages )

        a_EX(1,1:3) = [ 0.00_DP, 0.00_DP, 0.00_DP ]
        a_EX(2,1:3) = [ 1.00_DP, 0.00_DP, 0.00_DP ]
        a_EX(3,1:3) = [ 0.25_DP, 0.25_DP, 0.00_DP ]
        w_EX(1:3)   = [ 1.0_DP / 6.0_DP, &
                        1.0_DP / 6.0_DP, &
                        2.0_DP / 3.0_DP ]

      CASE ( 'IMEX_P_A2' )

        nStages = 3
        CALL AllocateButcherTables( nStages )

        ! --- Coefficients from Hu et al. (2017) ---
        ! --- arXiv:1708.06279v1, Section 2.6.1 ----

        a_EX(2,1) = 0.7369502715_DP
        a_EX(3,1) = 0.3215281691_DP
        a_EX(3,2) = 0.6784718309_DP

        w_EX(1)   = a_EX(3,1)
        w_EX(2)   = a_EX(3,2)

        a_IM(1,1) = 0.6286351712_DP
        a_IM(2,1) = 0.2431004655_DP
        a_IM(2,2) = 0.1959392570_DP
        a_IM(3,1) = 0.4803651051_DP
        a_IM(3,2) = 0.0746432814_DP
        a_IM(3,3) = 0.4449916135_DP

        w_IM(1)   = a_IM(3,1)
        w_IM(2)   = a_IM(3,2)
        w_IM(3)   = a_IM(3,3)

        alpha = 0.2797373792_DP

      CASE ( 'IMEX_P_ARS2' )

        nStages = 4
        CALL AllocateButcherTables( nStages )

        ! --- Coefficients from Hu et al. (2017) ---
        ! --- arXiv:1708.06279v1, Section 2.6.2 ----

        a_EX(3,1) = 1.0_DP
        a_EX(4,1) = 0.5_DP
        a_EX(4,3) = 0.5_DP

        w_EX(1)   = a_EX(4,1)
        w_EX(3)   = a_EX(4,3)

        a_IM(2,2) = 1.6_DP
        a_IM(3,2) = 0.3_DP
        a_IM(3,3) = 0.7_DP
        a_IM(4,2) = 0.5_DP
        a_IM(4,3) = 0.3_DP
        a_IM(4,4) = 0.2_DP

        w_IM(2)   = a_IM(4,2)
        w_IM(3)   = a_IM(4,3)
        w_IM(4)   = a_IM(4,4)

        alpha = 0.8_DP

      CASE ( 'IMEX_SSP2332' )

        ! --- Coefficients from Pareschi & Russo (2005) ---
        ! --- J. Sci. Comput. 25, 129-154 -----------------

        nStages = 3
        CALL AllocateButcherTables( nStages )

        a_EX(2,1) = 0.5_DP
        a_EX(3,1) = 0.5_DP
        a_EX(3,2) = 0.5_DP

        w_EX(1)   = 1.0_DP / 3.0_DP
        w_EX(2)   = 1.0_DP / 3.0_DP
        w_EX(3)   = 1.0_DP / 3.0_DP

        a_IM(1,1) = 0.25_DP
        a_IM(2,2) = 0.25_DP
        a_IM(3,1) = 1.0_DP / 3.0_DP
        a_IM(3,2) = 1.0_DP / 3.0_DP
        a_IM(3,3) = 1.0_DP / 3.0_DP

        w_IM(1)   = a_IM(3,1)
        w_IM(2)   = a_IM(3,2)
        w_IM(3)   = a_IM(3,3)

        alpha = 0.0_DP

      CASE ( 'IMEX_RKCB2' )

        ! --- Coefficients from Cavaglieri & Bewley (2015) ---
        ! --- JCP, 286, 172-193 ------------------------------

        nStages = 3
        CALL AllocateButcherTables( nStages )

        a_EX(2,1) = 0.4_DP
        a_EX(3,2) = 1.0_DP

        w_EX(2)   = 5.0_DP / 6.0_DP
        w_EX(3)   = 1.0_DP / 6.0_DP

        a_IM(2,2) = 0.4_DP
        a_IM(3,2) = 5.0_DP / 6.0_DP
        a_IM(3,3) = 1.0_DP / 6.0_DP

        w_IM(2)   = 5.0_DP / 6.0_DP
        w_IM(3)   = 1.0_DP / 6.0_DP

        alpha = 0.0_DP

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,'(A6,A,A)'), &
          '', 'Unknown Time Stepping Scheme: ', TRIM( Scheme )
        WRITE(*,*)
        WRITE(*,'(A6,A)'), &
          '', 'Available Options:'
        WRITE(*,*)
        WRITE(*,'(A6,A)') '', 'SSPRK1'
        WRITE(*,'(A6,A)') '', 'SSPRK2'
        WRITE(*,'(A6,A)') '', 'SSPRK3'
        WRITE(*,'(A6,A)') '', 'IMEX_P_A2'
        WRITE(*,'(A6,A)') '', 'IMEX_P_ARS2'
        WRITE(*,'(A6,A)') '', 'IMEX_SSP2332'
        WRITE(*,'(A6,A)') '', 'IMEX_RKCB2'

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
    WRITE(*,'(A6,A8,ES11.4E3)') '', 'alpha = ', alpha

    WRITE(*,*)
    WRITE(*,'(A6,A)') '', 'Explicit Butcher Table:'
    WRITE(*,'(A6,A)') '', '-----------------------'
    DO i = 1, nStages
      WRITE(*,'(A6,5ES14.4E3)') '', c_EX(i), a_EX(i,1:nStages)
    END DO
    WRITE(*,'(A6,A14,4ES14.4E3)') '', '', w_EX(1:nStages)

    ALLOCATE &
      ( U0(1:nDOF, &
           iZ_B0(1):iZ_E0(1), &
           iZ_B0(2):iZ_E0(2), &
           iZ_B0(3):iZ_E0(3), &
           iZ_B0(4):iZ_E0(4), &
           1:nCR, 1:nSpecies) )

    ALLOCATE &
      ( dU(1:nDOF, &
           iZ_B0(1):iZ_E0(1), &
           iZ_B0(2):iZ_E0(2), &
           iZ_B0(3):iZ_E0(3), &
           iZ_B0(4):iZ_E0(4), &
           1:nCR, 1:nSpecies) )

    nDOF_IMEX = nDOF * PRODUCT( iZ_E0 ) * nCR * nSpecies

    ALLOCATE( U_IMEX(1:nDOF_IMEX) )
    ALLOCATE( dU_IM (1:nDOF_IMEX,1:nStages) )
    ALLOCATE( dU_EX (1:nDOF_IMEX,1:nStages) )

  END SUBROUTINE Initialize_IMEX_RK


  SUBROUTINE Finalize_IMEX_RK

    DEALLOCATE( c_IM, w_IM, a_IM )
    DEALLOCATE( c_EX, w_EX, a_EX )
    DEALLOCATE( U0, dU )
    DEALLOCATE( U_IMEX, dU_IM, dU_EX )

  END SUBROUTINE Finalize_IMEX_RK


  SUBROUTINE AllocateButcherTables( nStages )

    INTEGER, INTENT(in) :: nStages

    alpha = Zero

    ALLOCATE( c_IM(nStages) )
    ALLOCATE( w_IM(nStages) )
    ALLOCATE( a_IM(nStages,nStages) )

    c_IM = Zero
    w_IM = Zero
    a_IM = Zero

    ALLOCATE( c_EX(nStages) )
    ALLOCATE( w_EX(nStages) )
    ALLOCATE( a_EX(nStages,nStages) )

    c_EX = Zero
    w_EX = Zero
    a_EX = Zero

  END SUBROUTINE AllocateButcherTables


  SUBROUTINE Update_IMEX_RK &
    ( dt, GE, GX, U, ComputeIncrement_Explicit, &
      ComputeIncrement_Implicit, ComputeCorrection_Implicit )

    REAL(DP), INTENT(in)         :: &
      dt
    REAL(DP), INTENT(in)         :: &
      GE(1:,iZ_B1(1):,1:)
    REAL(DP), INTENT(in)         :: &
      GX(1:,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:)
    REAL(DP), INTENT(inout)      :: &
      U (1:,iZ_B1(1):,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:,1:)
    PROCEDURE(IncrementExplicit) :: &
      ComputeIncrement_Explicit
    PROCEDURE(IncrementImplicit) :: &
      ComputeIncrement_Implicit, &
      ComputeCorrection_Implicit

    INTEGER :: iS, jS

    CALL InitializeStep_IMEX_RK( U )

    DO iS = 1, nStages

      CALL MapToStage( iZ_B0, U0, U_IMEX )

      DO jS = 1, iS - 1

        IF( a_IM(iS,jS) .NE. Zero )THEN

          U_IMEX(:) = U_IMEX(:) + dt * a_IM(iS,jS) * dU_IM(:,jS)

        END IF

        IF( a_EX(iS,jS) .NE. Zero )THEN

          U_IMEX(:) = U_IMEX(:) + dt * a_EX(iS,jS) * dU_EX(:,jS)

        END IF

      END DO

      IF( ANY( a_IM(:,iS) .NE. Zero ) .OR. ( w_IM(iS) .NE. Zero ) )THEN

        ! --- Implicit Solve ---

        CALL MapFromStage( iZ_B1, U, U_IMEX )

        CALL ComputeIncrement_Implicit &
               ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, &
                 dt * a_IM(iS,iS), GE, GX, U, dU )

        CALL MapToStage( iZ_B0, dU, dU_IM(:,iS) )

        U_IMEX(:) = U_IMEX(:) + dt * a_IM(iS,iS) * dU_IM(:,iS)

      END IF

      IF( ANY( a_EX(:,iS) .NE. Zero ) .OR. ( w_EX(iS) .NE. Zero ) )THEN

        ! --- Explicit Solve ---

        CALL MapFromStage( iZ_B1, U, U_IMEX )

        CALL ComputeIncrement_Explicit &
               ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, &
                 GE, GX, U, dU )

        CALL MapToStage( iZ_B0, dU, dU_EX(:,iS) )

      END IF

    END DO

    IF( ANY( a_IM(nStages,:) .NE. w_IM(:) ) &
          .OR. ANY( a_EX(nStages,:) .NE. w_EX(:) ) )THEN

      CALL MapToStage( iZ_B0, U0, U_IMEX )

      DO iS = 1, nStages

        IF( w_IM(iS) .NE. Zero )THEN

          U_IMEX(:) = U_IMEX(:) + dt * w_IM(iS) * dU_IM(:,iS)

        END IF

        IF( w_EX(iS) .NE. Zero )THEN

          U_IMEX(:) = U_IMEX(:) + dt * w_EX(iS) * dU_EX(:,iS)

        END IF

      END DO

    END IF

    CALL MapFromStage( iZ_B1, U, U_IMEX )

    IF( alpha > Zero )THEN

      CALL ComputeCorrection_Implicit &
             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, &
               alpha * dt**2, GE, GX, U, dU )

      U(:,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
          iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),:,:) &
        = U(:,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
              iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),:,:) &
          - alpha * dt**2 &
            * dU(:,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
                   iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),:,:)

    END IF

  END SUBROUTINE Update_IMEX_RK


  SUBROUTINE InitializeStep_IMEX_RK( U )

    REAL(DP), INTENT(in) :: &
      U(1:,iZ_B1(1):,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:,1:)

    U0(1:nDOF, &
       iZ_B0(1):iZ_E0(1), &
       iZ_B0(2):iZ_E0(2), &
       iZ_B0(3):iZ_E0(3), &
       iZ_B0(4):iZ_E0(4), &
       1:nCR,1:nSpecies) &
    = U(1:nDOF, &
        iZ_B0(1):iZ_E0(1), &
        iZ_B0(2):iZ_E0(2), &
        iZ_B0(3):iZ_E0(3), &
        iZ_B0(4):iZ_E0(4), &
        1:nCR,1:nSpecies)

    U_IMEX = Zero
    dU_IM  = Zero
    dU_EX  = Zero

  END SUBROUTINE InitializeStep_IMEX_RK


  SUBROUTINE MapToStage( iZ_B, U_7D, U_1D )

    INTEGER,  INTENT(in)  :: &
      iZ_B(4)
    REAL(DP), INTENT(in)  :: &
      U_7D(1:,iZ_B(1):,iZ_B(2):,iZ_B(3):,iZ_B(4):,1:,1:)
    REAL(DP), INTENT(out) :: &
      U_1D(1:)

    INTEGER :: i, iNode, iZ1, iZ2, iZ3, iZ4, iCR, iS

    i = 1
    DO iS = 1, nSpecies
      DO iCR = 1, nCR
        DO iZ4 = iZ_B0(4), iZ_E0(4)
          DO iZ3 = iZ_B0(3), iZ_E0(3)
            DO iZ2 = iZ_B0(2), iZ_E0(2)
              DO iZ1 = iZ_B0(1), iZ_E0(1)
                DO iNode = 1, nDOF

                  U_1D(i) = U_7D(iNode,iZ1,iZ2,iZ3,iZ4,iCR,iS)

                  i = i + 1

                END DO
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE MapToStage


  SUBROUTINE MapFromStage( iZ_B, U_7D, U_1D )

    INTEGER,  INTENT(in)  :: &
      iZ_B(4)
    REAL(DP), INTENT(out) :: &
      U_7D(1:,iZ_B(1):,iZ_B(2):,iZ_B(3):,iZ_B(4):,1:,1:)
    REAL(DP), INTENT(in)  :: &
      U_1D(1:)

    INTEGER :: i, iNode, iZ1, iZ2, iZ3, iZ4, iCR, iS

    i = 1
    DO iS = 1, nSpecies
      DO iCR = 1, nCR
        DO iZ4 = iZ_B0(4), iZ_E0(4)
          DO iZ3 = iZ_B0(3), iZ_E0(3)
            DO iZ2 = iZ_B0(2), iZ_E0(2)
              DO iZ1 = iZ_B0(1), iZ_E0(1)
                DO iNode = 1, nDOF

                  U_7D(iNode,iZ1,iZ2,iZ3,iZ4,iCR,iS) = U_1D(i)

                  i = i + 1

                END DO
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE MapFromStage


END MODULE TimeSteppingModule_IMEX_RK
