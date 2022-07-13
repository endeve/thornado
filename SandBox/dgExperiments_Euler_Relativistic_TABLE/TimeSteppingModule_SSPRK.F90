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
  USE FluidFieldsModule, ONLY: &
    nCF
  USE Euler_SlopeLimiterModule_Relativistic_TABLE, ONLY: &
    ApplySlopeLimiter_Euler_Relativistic_TABLE
  USE Euler_PositivityLimiterModule_Relativistic_TABLE, ONLY: &
    ApplyPositivityLimiter_Euler_Relativistic_TABLE
  USE TimersModule_Euler, ONLY: &
    TimersStart_Euler, &
    TimersStop_Euler, &
    Timer_Euler_UpdateFluid
  USE Euler_dgDiscretizationModule, ONLY: &
    OffGridFlux_Euler_X1_Inner, &
    OffGridFlux_Euler_X1_Outer, &
    OffGridFlux_Euler_X2_Inner, &
    OffGridFlux_Euler_X2_Outer, &
    OffGridFlux_Euler_X3_Inner, &
    OffGridFlux_Euler_X3_Outer
  USE Euler_TallyModule_Relativistic, ONLY: &
    IncrementOffGridTally_Euler_Relativistic

  IMPLICIT NONE
  PRIVATE

  INTEGER :: nStages_SSPRK
  REAL(DP), DIMENSION(:),   ALLOCATABLE :: c_SSPRK
  REAL(DP), DIMENSION(:),   ALLOCATABLE :: w_SSPRK
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: a_SSPRK

  REAL(DP), DIMENSION(:,:,:,:,:),   ALLOCATABLE :: U_SSPRK
  REAL(DP), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: D_SSPRK

  PUBLIC :: InitializeFluid_SSPRK
  PUBLIC :: UpdateFluid_SSPRK
  PUBLIC :: FinalizeFluid_SSPRK

  INTERFACE
    SUBROUTINE FluidIncrement &
      ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dU, &
        SuppressBC_Option, UseXCFC_Option, &
        SurfaceFlux_X1_Option, &
        SurfaceFlux_X2_Option, &
        SurfaceFlux_X3_Option )
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
      REAL(DP), INTENT(out), OPTIONAL :: &
        SurfaceFlux_X1_Option(:,:,:,:,:), &
        SurfaceFlux_X2_Option(:,:,:,:,:), &
        SurfaceFlux_X3_Option(:,:,:,:,:)
    END SUBROUTINE FluidIncrement
  END INTERFACE


CONTAINS


  SUBROUTINE InitializeFluid_SSPRK( nStages )

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
                 1:nCF) )

    ALLOCATE( D_SSPRK &
                (1:nDOFX, &
                 iX_B1(1):iX_E1(1), &
                 iX_B1(2):iX_E1(2), &
                 iX_B1(3):iX_E1(3), &
                 1:nCF,1:nStages) )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to:    a_SSPRK, c_SSPRK, w_SSPRK ) &
    !$OMP MAP( alloc: U_SSPRK, D_SSPRK )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN(     a_SSPRK, c_SSPRK, w_SSPRK ) &
    !$ACC CREATE(     U_SSPRK, D_SSPRK )
#endif

  END SUBROUTINE InitializeFluid_SSPRK


  SUBROUTINE FinalizeFluid_SSPRK

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: U_SSPRK, D_SSPRK, a_SSPRK, c_SSPRK, w_SSPRK )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE(       U_SSPRK, D_SSPRK, a_SSPRK, c_SSPRK, w_SSPRK )
#endif

    DEALLOCATE( a_SSPRK, c_SSPRK, w_SSPRK )

    DEALLOCATE( U_SSPRK, D_SSPRK )

  END SUBROUTINE FinalizeFluid_SSPRK


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


  SUBROUTINE UpdateFluid_SSPRK &
    ( t, dt, G, U, D, ComputeIncrement_Fluid )

    REAL(DP), INTENT(in)    :: &
      t, dt
    REAL(DP), INTENT(inout) :: &
      G(:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    REAL(DP), INTENT(inout) :: &
      U(:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    REAL(DP), INTENT(out)   :: &
      D(:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    PROCEDURE (FluidIncrement) :: &
      ComputeIncrement_Fluid

    LOGICAL :: DEBUG = .FALSE.

    INTEGER :: iNX, iX1, iX2, iX3, iCF
    INTEGER :: iS, jS

    REAL(DP) :: dM_OffGrid_Euler(nCF)

    dM_OffGrid_Euler = Zero

    CALL TimersStart_Euler( Timer_Euler_UpdateFluid )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: iX_B0, iX_E0, iX_B1, iX_E1, G, U, D )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN(  iX_B0, iX_E0, iX_B1, iX_E1, G, U, D )
#endif

    DO iS = 1, nStages_SSPRK

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
      !$ACC PRESENT( iX_B1, iX_E1, U_SSPRK, U )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO COLLAPSE(5)
#endif
      DO iCF = 1, nCF
      DO iX3 = iX_B1(3), iX_E1(3)
      DO iX2 = iX_B1(2), iX_E1(2)
      DO iX1 = iX_B1(1), iX_E1(1)
      DO iNX = 1, nDOFX

        U_SSPRK(iNX,iX1,iX2,iX3,iCF) = U(iNX,iX1,iX2,iX3,iCF)

      END DO
      END DO
      END DO
      END DO
      END DO

      DO jS = 1, iS - 1

        IF( a_SSPRK(iS,jS) .NE. Zero )THEN

          IF( DEBUG ) WRITE(*,'(A)') 'CALL AddIncrement_Fluid (1)'
          CALL AddIncrement_Fluid &
                 ( jS, One, U_SSPRK, dt * a_SSPRK(iS,jS), D_SSPRK )

        END IF

      END DO

      IF( ANY( a_SSPRK(:,iS) .NE. Zero ) &
          .OR. ( w_SSPRK(iS) .NE. Zero ) )THEN

        CALL ApplySlopeLimiter_Euler_Relativistic_TABLE &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G, U_SSPRK, D )

        CALL ApplyPositivityLimiter_Euler_Relativistic_TABLE &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G, U_SSPRK )

        CALL ComputeIncrement_Fluid &
               ( iX_B0, iX_E0, iX_B1, iX_E1, &
                 G, U_SSPRK, D, D_SSPRK(:,:,:,:,:,iS) )

        dM_OffGrid_Euler &
          = dM_OffGrid_Euler &
              + dt * w_SSPRK(iS) &
                  * (   OffGridFlux_Euler_X1_Outer &
                      - OffGridFlux_Euler_X1_Inner &
                      + OffGridFlux_Euler_X2_Outer &
                      - OffGridFlux_Euler_X2_Inner &
                      + OffGridFlux_Euler_X3_Outer &
                      - OffGridFlux_Euler_X3_Inner )

      END IF

    END DO

    DO iS = 1, nStages_SSPRK

      IF( w_SSPRK(iS) .NE. Zero )THEN

        CALL AddIncrement_Fluid &
               ( iS, One, U, dt * w_SSPRK(iS), D_SSPRK )

      END IF

    END DO

    CALL ApplySlopeLimiter_Euler_Relativistic_TABLE &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D )

    CALL ApplyPositivityLimiter_Euler_Relativistic_TABLE &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from:    U, D ) &
    !$OMP MAP( release: iX_B0, iX_E0, iX_B1, iX_E1, G )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC COPYOUT(      U, D ) &
    !$ACC DELETE(       iX_B0, iX_E0, iX_B1, iX_E1, G )
#endif

    CALL IncrementOffGridTally_Euler_Relativistic( dM_OffGrid_Euler )

    CALL TimersStop_Euler( Timer_Euler_UpdateFluid )

  END SUBROUTINE UpdateFluid_SSPRK


  SUBROUTINE AddIncrement_Fluid( iS, alpha, U, beta, D )

    INTEGER,  INTENT(in)    :: &
      iS
    REAL(DP), INTENT(in)    :: &
      alpha, beta
    REAL(DP), INTENT(inout) :: &
      U(:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    REAL(DP), INTENT(in)    :: &
      D(:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:,:)

    INTEGER :: iNX, iX1, iX2, iX3, iCF

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( iX_B1, iX_E1, U, D )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iCF = 1, nCF
    DO iX3 = iX_B1(3), iX_E1(3)
    DO iX2 = iX_B1(2), iX_E1(2)
    DO iX1 = iX_B1(1), iX_E1(1)
    DO iNX = 1, nDOFX

      U(iNX,iX1,iX2,iX3,iCF) &
        = alpha * U(iNX,iX1,iX2,iX3,iCF) &
            + beta * D(iNX,iX1,iX2,iX3,iCF,iS)

    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE AddIncrement_Fluid


END MODULE TimeSteppingModule_SSPRK
