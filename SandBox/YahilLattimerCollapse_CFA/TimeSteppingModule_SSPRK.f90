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
  USE GeometryFieldsModule, ONLY: &
    iGF_Psi
  USE GravitySolutionModule_CFA_Poseidon, ONLY: &
    ComputeConformalFactor_Poseidon, &
    ComputeGeometry_Poseidon
  USE FluidFieldsModule, ONLY: &
    nCF
  USE Euler_SlopeLimiterModule_Relativistic_TABLE, ONLY: &
    ApplySlopeLimiter_Euler_Relativistic_TABLE
  USE Euler_PositivityLimiterModule_Relativistic_TABLE, ONLY: &
    ApplyPositivityLimiter_Euler_Relativistic_TABLE
  USE Poseidon_UtilitiesModule, ONLY: &
    ComputeMatterSources_Poseidon, &
    ComputePressureTensorTrace_Poseidon, &
    MultiplyByPsi6, &
    DivideByPsi6
  USE TimersModule_Euler, ONLY: &
    TimersStart_Euler, &
    TimersStop_Euler,  &
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

  REAL(DP), DIMENSION(:,:,:,:,:),   ALLOCATABLE :: Ustar
  REAL(DP), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: Dstar

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

    ALLOCATE( Ustar &
                (1:nDOFX, &
                 iX_B1(1):iX_E1(1), &
                 iX_B1(2):iX_E1(2), &
                 iX_B1(3):iX_E1(3), &
                 1:nCF) )

    ALLOCATE( Dstar &
                (1:nDOFX, &
                 iX_B1(1):iX_E1(1), &
                 iX_B1(2):iX_E1(2), &
                 iX_B1(3):iX_E1(3), &
                 1:nCF,1:nStages) )

  END SUBROUTINE InitializeFluid_SSPRK


  SUBROUTINE FinalizeFluid_SSPRK

    DEALLOCATE( a_SSPRK, c_SSPRK, w_SSPRK )

    DEALLOCATE( Ustar, Dstar )

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
      U(:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:), &
      D(:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    PROCEDURE (FluidIncrement) :: &
      ComputeIncrement_Fluid

    INTEGER :: iS, jS, iX1, iX2, iX3

    REAL(DP) :: E (nDOFX,iX_B0(1):iX_E0(1), &
                         iX_B0(2):iX_E0(2), &
                         iX_B0(3):iX_E0(3))
    REAL(DP) :: Si(nDOFX,iX_B0(1):iX_E0(1), &
                         iX_B0(2):iX_E0(2), &
                         iX_B0(3):iX_E0(3),3)
    REAL(DP) :: S (nDOFX,iX_B0(1):iX_E0(1), &
                         iX_B0(2):iX_E0(2), &
                         iX_B0(3):iX_E0(3))
    REAL(DP) :: Mg(nDOFX,iX_B0(1):iX_E0(1), &
                         iX_B0(2):iX_E0(2), &
                         iX_B0(3):iX_E0(3))

    REAL(DP) :: dM_OffGrid_Euler(nCF)

    dM_OffGrid_Euler = Zero

    CALL TimersStart_Euler( Timer_Euler_UpdateFluid )

    CALL MultiplyByPsi6( iX_B1, iX_E1, G, U ) ! Ustar = psi^6 * U

    Dstar = Zero ! --- Increment

    DO iS = 1, nStages_SSPRK

      Ustar = U

      DO jS = 1, iS - 1

        IF( a_SSPRK(iS,jS) .NE. Zero )THEN

          CALL AddIncrement_Fluid &
                 ( One, Ustar, dt * a_SSPRK(iS,jS), Dstar(:,:,:,:,:,jS) )

        END IF

      END DO

      IF( ANY( a_SSPRK(:,iS) .NE. Zero ) &
          .OR. ( w_SSPRK(iS) .NE. Zero ) )THEN

        IF( iS .NE. 1 )THEN ! At first stage, U and psi are known

          CALL ComputeMatterSources_Poseidon &
                 ( iX_B0, iX_E0, iX_B1, iX_E1, G, Ustar, E, Si, Mg )

          CALL ComputeConformalFactor_Poseidon &
                 ( iX_B0, iX_E0, iX_B1, iX_E1, E, Si, Mg, G )

          CALL DivideByPsi6( iX_B1, iX_E1, G, Ustar )

          CALL ApplySlopeLimiter_Euler_Relativistic_TABLE &
                 ( iX_B0, iX_E0, iX_B1, iX_E1, G, Ustar, D )

          CALL ApplyPositivityLimiter_Euler_Relativistic_TABLE &
                 ( iX_B0, iX_E0, iX_B1, iX_E1, G, Ustar )

          CALL MultiplyByPsi6( iX_B1, iX_E1, G, Ustar )

          CALL ComputeMatterSources_Poseidon &
                 ( iX_B0, iX_E0, iX_B1, iX_E1, G, Ustar, E, Si, Mg )

          CALL ComputeConformalFactor_Poseidon &
                 ( iX_B0, iX_E0, iX_B1, iX_E1, E, Si, Mg, G )

          CALL ComputePressureTensorTrace_Poseidon &
                 ( iX_B0, iX_E0, iX_B1, iX_E1, G, Ustar, S )

          CALL ComputeGeometry_Poseidon &
                 ( iX_B0, iX_E0, iX_B1, iX_E1, E, S, Si, G )

        END IF

        CALL DivideByPsi6( iX_B1, iX_E1, G, Ustar )

        CALL ComputeIncrement_Fluid &
               ( iX_B0, iX_E0, iX_B1, iX_E1, &
                 G, Ustar, D, Dstar(:,:,:,:,:,iS), &
                 UseXCFC_Option = .TRUE. )

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
               ( One, U, dt * w_SSPRK(iS), Dstar(:,:,:,:,:,iS) )

      END IF

    END DO

    CALL ComputeMatterSources_Poseidon &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, E, Si, Mg )

    CALL ComputeConformalFactor_Poseidon &
           ( iX_B0, iX_E0, iX_B1, iX_E1, E, Si, Mg, G )

    CALL DivideByPsi6( iX_B1, iX_E1, G, U )

    CALL ApplySlopeLimiter_Euler_Relativistic_TABLE &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D )

    CALL ApplyPositivityLimiter_Euler_Relativistic_TABLE &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U )

    CALL MultiplyByPsi6( iX_B1, iX_E1, G, U )

    CALL ComputeMatterSources_Poseidon &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, E, Si, Mg )

    CALL ComputeConformalFactor_Poseidon &
           ( iX_B0, iX_E0, iX_B1, iX_E1, E, Si, Mg, G )

    CALL ComputePressureTensorTrace_Poseidon &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, S )

    CALL ComputeGeometry_Poseidon &
           ( iX_B0, iX_E0, iX_B1, iX_E1, E, S, Si, G )

    CALL DivideByPsi6( iX_B1, iX_E1, G, U )

    CALL IncrementOffGridTally_Euler_Relativistic( dM_OffGrid_Euler )

    CALL TimersStop_Euler( Timer_Euler_UpdateFluid )

  END SUBROUTINE UpdateFluid_SSPRK


  ! --- PRIVATE Subroutines ---


  SUBROUTINE AddIncrement_Fluid( alpha, U, beta, D )

    REAL(DP), INTENT(in)    :: &
      alpha, beta
    REAL(DP), INTENT(inout) :: &
      U(:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    REAL(DP), INTENT(in)    :: &
      D(:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)

    INTEGER :: iCF, iNX, iX1, iX2, iX3

    DO iCF = 1       , nCF
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1       , nDOFX

      U(iNX,iX1,iX2,iX3,iCF) &
        = alpha * U(iNX,iX1,iX2,iX3,iCF) &
            + beta * D(iNX,iX1,iX2,iX3,iCF)

    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE AddIncrement_Fluid


END MODULE TimeSteppingModule_SSPRK
