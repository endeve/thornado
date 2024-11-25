MODULE TimeSteppingModule

  USE KindModule, ONLY: &
    DP, Zero, Fourth, Third
  USE UnitsModule, ONLY: &
    UnitsDisplay
  USE ProgramHeaderModule, ONLY: &
    nX, nDOFX, nNodes, &
    nE, nDOF, &
    iX_B0, iX_E0, &
    iX_B1, iX_E1
  USE MeshModule, ONLY: &
    MeshX
  USE FluidFieldsModule, ONLY: &
    uCF, rhsCF, nCF, &
    uPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne, nPF, &
    uAF, iAF_P, iAF_T,  iAF_Ye, iAF_S,  iAF_E, iAF_Gm, iAF_Cs, nAF
  USE RadiationFieldsModule, ONLY: &
    uCR, rhsCR, nCR, nSpecies
  USE EquationOfStateModule, ONLY: &
    ComputeAuxiliary_Fluid
  USE InputOutputModule, ONLY: &
    WriteFields1D, &
    WriteFieldsRestart1D
  USE TallyModule, ONLY: &
    InitializeGlobalTally, &
    ComputeGlobalTally
  USE GravitySolutionModule, ONLY: &
    SolveGravity
  USE EulerEquationsUtilitiesModule, ONLY: &
    ComputePrimitive, &
    Eigenvalues
  USE FluidEvolutionModule, ONLY: &
    ComputeRHS_Fluid, &
    ApplySlopeLimiter_Fluid, &
    ApplyPositivityLimiter_Fluid
  USE RadiationEvolutionModule, ONLY: &
    ComputeRHS_Radiation, &
    ComputeExplicitIncrement_Radiation, &
    ApplySlopeLimiter_Radiation, &
    ApplyPositivityLimiter_Radiation
  USE FluidRadiationCouplingModule, ONLY: &
    CoupleFluidRadiation, &
    ComputeImplicitIncrement_FluidRadiation
  USE BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_Fluid, &
    ApplyBoundaryConditions_Radiation

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  LOGICAL :: EvolveGravity   = .FALSE.
  LOGICAL :: EvolveFluid     = .FALSE.
  LOGICAL :: EvolveRadiation = .FALSE.
  INTEGER :: nStages_SSP_RK  = 1
  INTEGER :: nStages_SI_RK   = 1
  INTEGER :: nDOF_SSP        = 0
  INTEGER :: nStages_IMEX    = 1
  INTEGER :: nDOF_IMEX       = 0
  REAL(DP), DIMENSION(:),             ALLOCATABLE :: c_SSP
  REAL(DP), DIMENSION(:),             ALLOCATABLE :: w_SSP
  REAL(DP), DIMENSION(:),             ALLOCATABLE :: U_SSP
  REAL(DP), DIMENSION(:,:),           ALLOCATABLE :: a_SSP
  REAL(DP), DIMENSION(:,:),           ALLOCATABLE :: F_SSP
  REAL(DP), DIMENSION(:),             ALLOCATABLE :: c_IM, c_EX
  REAL(DP), DIMENSION(:),             ALLOCATABLE :: w_IM, w_EX
  REAL(DP), DIMENSION(:),             ALLOCATABLE :: U_IMEX
  REAL(DP), DIMENSION(:,:),           ALLOCATABLE :: a_IM, a_EX
  REAL(DP), DIMENSION(:,:),           ALLOCATABLE :: F_IM, F_EX
  REAL(DP), DIMENSION(:,:,:,:,:),     ALLOCATABLE :: uCF_0
  REAL(DP), DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE :: uCR_0

  PROCEDURE (TimeStepper), POINTER, PUBLIC :: &
    SSP_RK => NULL(), &
    SI_RK  => NULL(), &
    SSPRKs => NULL(), &
    IMEX   => NULL()

  INTERFACE
    SUBROUTINE TimeStepper( t, dt )
      USE KindModule, ONLY: DP
      REAL(DP), INTENT(in) :: t, dt
    END SUBROUTINE TimeStepper
  END INTERFACE

  PUBLIC :: InitializeTimeStepping
  PUBLIC :: FinalizeTimeStepping
  PUBLIC :: EvolveFields
  PUBLIC :: BackwardEuler

CONTAINS


  SUBROUTINE InitializeTimeStepping &
               ( EvolveGravity_Option, EvolveFluid_Option, &
                 EvolveRadiation_Option, nStages_SSP_RK_Option, &
                 nStages_SI_RK_Option, IMEX_Scheme_Option )

    LOGICAL,      INTENT(in), OPTIONAL :: EvolveGravity_Option
    LOGICAL,      INTENT(in), OPTIONAL :: EvolveFluid_Option
    LOGICAL,      INTENT(in), OPTIONAL :: EvolveRadiation_Option
    INTEGER,      INTENT(in), OPTIONAL :: nStages_SSP_RK_Option
    INTEGER,      INTENT(in), OPTIONAL :: nStages_SI_RK_Option
    CHARACTER(*), INTENT(in), OPTIONAL :: IMEX_Scheme_Option

    CHARACTER(32) :: IMEX_Scheme

    IF( PRESENT( EvolveGravity_Option ) )THEN
      EvolveGravity = EvolveGravity_Option
    END IF
    WRITE(*,*)
    WRITE(*,'(A5,A19,L1)') '', &
      'Evolve Gravity = ', EvolveGravity

    IF( PRESENT( EvolveFluid_Option ) )THEN
      EvolveFluid = EvolveFluid_Option
    END IF
    WRITE(*,'(A5,A19,L1)') '', &
      'Evolve Fluid = ', EvolveFluid

    IF( PRESENT( EvolveRadiation_Option ) )THEN
      EvolveRadiation = EvolveRadiation_Option
    END IF
    WRITE(*,'(A5,A19,L1)') '', &
      'Evolve Radiation = ', EvolveRadiation
    WRITE(*,*)

    nStages_SSP_RK = 0

    IF( PRESENT( nStages_SSP_RK_Option ) )THEN

      nStages_SSP_RK = nStages_SSP_RK_Option

    END IF

    SELECT CASE ( nStages_SSP_RK )
      CASE ( 1 )

        SSP_RK => SSP_RK1

      CASE ( 2 )

        SSP_RK => SSP_RK2

      CASE ( 3 )

        SSP_RK => SSP_RK3

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,'(A4,A43,I2.2)') &
          '', 'SSP_RK not implemented for nStatgesSSPRK = ', nStages_SSP_RK
        STOP

    END SELECT

    IF( PRESENT( nStages_SI_RK_Option ) )THEN

      nStages_SI_RK = nStages_SI_RK_Option

    END IF

    SELECT CASE ( nStages_SI_RK )
      CASE ( 1 )

        SI_RK => SI_RK1

      CASE ( 2 )

        SI_RK => SI_RK2

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,'(A4,A43,I2.2)') &
          '', 'SI_RK not implemented for nStatges_SI_RK = ', nStages_SI_RK
        STOP

    END SELECT

    CALL InitializeTimeStepper_SSPRK( nStages_SSP_RK )

    IMEX_Scheme = 'DUMMY'
    IF( PRESENT( IMEX_Scheme_Option ) )THEN

      IMEX_Scheme = TRIM( IMEX_Scheme_Option )

    END IF

    CALL InitializeTimeStepper_IMEX( IMEX_Scheme )

  END SUBROUTINE InitializeTimeStepping


  SUBROUTINE FinalizeTimeStepping

    NULLIFY( SSP_RK, SI_RK, IMEX )

  END SUBROUTINE FinalizeTimeStepping


  SUBROUTINE TimeStepper_DUMMY( t, dt )

    REAL(DP), INTENT(in) :: t, dt

    WRITE(*,*)
    WRITE(*,'(A5,A)') '', 'TimeStepper_DUMMY'
    WRITE(*,*)

  END SUBROUTINE TimeStepper_DUMMY


  SUBROUTINE EvolveFields &
               ( t_begin, t_end, dt_write, UpdateFields, dt_fixed_Option )

    REAL(DP), INTENT(in) :: t_begin, t_end, dt_write
    PROCEDURE (TimeStepper) :: UpdateFields
    REAL(DP), INTENT(in), OPTIONAL :: dt_fixed_Option

    LOGICAL  :: WriteOutput = .FALSE.
    LOGICAL  :: FixedTimeStep
    INTEGER  :: iCycle
    REAL(DP) :: t, t_write, dt, dt_fixed
    REAL(DP), DIMENSION(0:1) :: WallTime

    FixedTimeStep = .FALSE.
    IF( PRESENT( dt_fixed_Option ) )THEN
      FixedTimeStep = .TRUE.
      dt_fixed = dt_fixed_Option
    END IF

    ASSOCIATE( UD => UnitsDisplay )

    WRITE(*,*)
    WRITE(*,'(A4,A21)') '', 'INFO: Evolving Fields'
    WRITE(*,*)
    WRITE(*,'(A6,A10,ES10.4E2,A1,2A2,A8,ES10.4E2,&
              &A1,2A2,A11,ES10.4E2,A1,A2)') &
      '', 't_begin = ',  &
      t_begin  / UD % TimeUnit, '', TRIM( UD % TimeLabel ), &
      '', 't_end = ',    &
      t_end    / UD % TimeUnit, '', TRIM( UD % TimeLabel ), &
      '', 'dt_write = ', &
      dt_write / UD % TimeUnit, '', TRIM( UD % TimeLabel )
    WRITE(*,*)

    iCycle  = 0
    t       = t_begin
    t_write = t_begin + dt_write

    IF( EvolveFluid )THEN

      CALL ApplyPositivityLimiter_Fluid

    END IF

    IF( EvolveRadiation )THEN

      CALL ApplyPositivityLimiter_Radiation &
             ( iX_B0 = iX_B0, iX_E0 = iX_E0, &
               iX_B1 = iX_B1, iX_E1 = iX_E1, &
               U = uCR )

    END IF

    IF( EvolveGravity )THEN

      CALL SolveGravity

    END IF

    CALL InitializeGlobalTally &
           ( Time = t, &
             TallyGravity_Option = EvolveGravity, &
             TallyFluid_Option = EvolveFluid,  &
             TallyRadiation_Option = EvolveRadiation )

    CALL WriteFields1D &
           ( Time = t, &
             WriteGeometryFields_Option = EvolveGravity, &
             WriteFluidFields_Option = EvolveFluid, &
             WriteRadiationFields_Option = EvolveRadiation )

    CALL WriteFieldsRestart1D &
           ( Time = t, &
             WriteGeometryFields_Option = EvolveGravity, &
             WriteFluidFields_Option = EvolveFluid, &
             WriteRadiationFields_Option = EvolveRadiation )

    WallTime(0) = MPI_WTIME( )

    DO WHILE( t < t_end )

      iCycle = iCycle + 1

      IF( FixedTimeStep )THEN
        dt = dt_fixed
      ELSE
        CALL ComputeTimeStep( dt )
      END IF

      IF( t + dt > t_end )THEN

        dt = t_end - t

      END IF

      IF( t + dt > t_write )THEN

        dt          = t_write - t
        t_write     = t_write + dt_write
        WriteOutput = .TRUE.

      END IF

      IF( MOD( iCycle, 100 ) == 0 )THEN

        WRITE(*,'(A8,A8,I8.8,A2,A4,ES12.6E2,A1,A2,A2,A5,ES12.6E2,A1,A2)') &
          '', 'Cycle = ', iCycle, &
          '', 't = ',  t  / UD % TimeUnit, '', TRIM( UD % TimeLabel ), &
          '', 'dt = ', dt / UD % TimeUnit, '', TRIM( UD % TimeLabel )

        CALL ComputeGlobalTally( Time = t )

      END IF

      CALL UpdateFields( t, dt )

      t = t + dt

      IF( WriteOutput )THEN

        CALL WriteFields1D &
               ( Time = t, &
                 WriteGeometryFields_Option = EvolveGravity, &
                 WriteFluidFields_Option = EvolveFluid, &
                 WriteRadiationFields_Option = EvolveRadiation )

        CALL WriteFieldsRestart1D &
               ( Time = t, &
                 WriteGeometryFields_Option = EvolveGravity, &
                 WriteFluidFields_Option = EvolveFluid, &
                 WriteRadiationFields_Option = EvolveRadiation )

        WriteOutput = .FALSE.

      END IF

    END DO

    WallTime(1) = MPI_WTIME( )

    WRITE(*,*)
    WRITE(*,'(A6,A15,ES10.4E2,A1,A2,A6,I8.8,A7,A4,ES10.4E2,A2)') &
      '', 'Evolved to t = ', t / UD % TimeUnit, '', TRIM( UD % TimeLabel ), &
      ' with ', iCycle, ' cycles', &
      ' in ', WallTime(1)-WallTime(0), ' s'
    WRITE(*,*)

    CALL WriteFields1D &
           ( Time = t, &
             WriteGeometryFields_Option = EvolveGravity, &
             WriteFluidFields_Option = EvolveFluid, &
             WriteRadiationFields_Option = EvolveRadiation )

    CALL WriteFieldsRestart1D &
           ( Time = t, &
             WriteGeometryFields_Option = EvolveGravity, &
             WriteFluidFields_Option = EvolveFluid, &
             WriteRadiationFields_Option = EvolveRadiation )

    END ASSOCIATE ! UD

  END SUBROUTINE EvolveFields


  SUBROUTINE ComputeTimeStep( dt )

    REAL(DP), INTENT(out) :: dt

    REAL(DP) :: dt_Fluid, dt_Radiation

    dt_Fluid     = HUGE( 1.0_DP )
    dt_Radiation = HUGE( 1.0_DP )

    IF( EvolveFluid )THEN

      CALL ComputeTimeStep_Fluid( dt_Fluid )

    END IF

    IF( EvolveRadiation )THEN

      CALL ComputeTimeStep_Radiation( dt_Radiation )

    END IF

    dt = MIN( dt_Fluid, dt_Radiation )

  END SUBROUTINE ComputeTimeStep


  SUBROUTINE ComputeTimeStep_Fluid( dt )

    REAL(DP), INTENT(out) :: dt

    INTEGER                    :: iX1, iX2, iX3, iNodeX
    REAL(DP)                   :: CFL, dt_X1, dt_X2, dt_X3
    REAL(DP), DIMENSION(1:nPF) :: uPF_N

    dt = HUGE( 1.0_DP )

    ASSOCIATE &
      ( dX1 => MeshX(1) % Width (1:nX(1)), &
        dX2 => MeshX(2) % Width (1:nX(2)), &
        dX3 => MeshX(3) % Width (1:nX(3)), &
        X1  => MeshX(1) % Center(1:nX(1)), &
        X2  => MeshX(2) % Center(1:nX(2)), &
        X3  => MeshX(3) % Center(1:nX(3))  )

    CFL = 0.2_DP / ( 2.0_DP * DBLE( nNodes - 1 ) + 1.0_DP ) ! For Debugging

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          CALL ComputePrimitive &
                 ( uCF(:,iX1,iX2,iX3,1:nCF), uPF(:,iX1,iX2,iX3,1:nPF) )

          CALL ComputeAuxiliary_Fluid &
                 ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_E ), &
                   uPF(:,iX1,iX2,iX3,iPF_Ne), uAF(:,iX1,iX2,iX3,iAF_P ), &
                   uAF(:,iX1,iX2,iX3,iAF_T ), uAF(:,iX1,iX2,iX3,iAF_Ye), &
                   uAF(:,iX1,iX2,iX3,iAF_S ), uAF(:,iX1,iX2,iX3,iAF_E ), &
                   uAF(:,iX1,iX2,iX3,iAF_Gm), uAF(:,iX1,iX2,iX3,iAF_Cs) )

          DO iNodeX = 1, nDOFX

            dt_X1 &
              = CFL * dX1(iX1) &
                   / MAXVAL(ABS(Eigenvalues &
                     ( uPF(iNodeX,iX1,iX2,iX3,iPF_V1), &
                       uAF(iNodeX,iX1,iX2,iX3,iAF_Cs) )))

            dt_X2 &
              = CFL * dX2(iX2) &
                    / MAXVAL(ABS(Eigenvalues &
                      ( uPF(iNodeX,iX1,iX2,iX3,iPF_V2), &
                        uAF(iNodeX,iX1,iX2,iX3,iAF_Cs) )))

            dt_X3 &
              = CFL * dX3(iX3) &
                   / MAXVAL(ABS(Eigenvalues &
                     ( uPF(iNodeX,iX1,iX2,iX3,iPF_V3), &
                       uAF(iNodeX,iX1,iX2,iX3,iAF_Cs) )))

            dt = MIN( dt, dt_x1, dt_X2, dt_X3 )

          END DO
        END DO
      END DO
    END DO

    END ASSOCIATE ! dX1, etc.

  END SUBROUTINE ComputeTimeStep_Fluid


  SUBROUTINE ComputeTimeStep_Radiation( dt )

    REAL(DP), INTENT(out) :: dt

    INTEGER  :: iX1, iX2, iX3, iNodeX
    REAL(DP) :: CFL, dt_X1, dt_X2, dt_X3

    dt = HUGE( 1.0_DP )

    ASSOCIATE( dX1 => MeshX(1) % Width (1:nX(1)), &
               dX2 => MeshX(2) % Width (1:nX(2)), &
               dX3 => MeshX(3) % Width (1:nX(3)), &
               X1  => MeshX(1) % Center(1:nX(1)), &
               X2  => MeshX(2) % Center(1:nX(2)), &
               X3  => MeshX(3) % Center(1:nX(3)) )

    CFL = 0.2_DP / ( 2.0_DP * DBLE( nNodes - 1 ) + 1.0_DP ) ! For Debugging

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          dt_X1 = CFL * dX1(iX1)
          dt_X2 = CFL * dX2(iX2)
          dt_X3 = CFL * dX3(iX3)

          dt = MIN( dt, dt_x1, dt_X2, dt_X3 )

        END DO
      END DO
    END DO

    END ASSOCIATE ! dX1, etc.

  END SUBROUTINE ComputeTimeStep_Radiation


  SUBROUTINE SSP_RK1( t, dt )

    REAL(DP), INTENT(in) :: t, dt

    CALL Initialize_SSP_RK

    IF( EvolveFluid )THEN

      CALL ApplyBoundaryConditions_Fluid

      CALL ComputeRHS_Fluid &
             ( iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ] )

    END IF

    IF( EvolveRadiation )THEN

      CALL ApplyBoundaryConditions_Radiation( t )

      CALL ComputeRHS_Radiation &
             ( iX_B0 = iX_B0, iX_E0 = iX_E0, &
               iX_B1 = iX_B1, iX_E1 = iX_E1, &
               U = uCR, dU = rhsCR )

    END IF

    IF( EvolveFluid )THEN

      CALL ApplyRHS_Fluid &
             ( iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ], &
               dt = dt, alpha = 0.0_DP, beta = 1.0_DP )

      CALL ApplyBoundaryConditions_Fluid

      CALL ApplySlopeLimiter_Fluid

      CALL ApplyPositivityLimiter_Fluid

    END IF

    IF( EvolveRadiation )THEN

      CALL ApplyRHS_Radiation &
             ( iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ], &
               dt = dt, alpha = 0.0_DP, beta = 1.0_DP )

      CALL ApplyBoundaryConditions_Radiation &
             ( t + dt, LimiterBC_Option = .TRUE. )

      CALL ApplySlopeLimiter_Radiation &
             ( iX_B0 = iX_B0, iX_E0 = iX_E0, &
               iX_B1 = iX_B1, iX_E1 = iX_E1, &
               U = uCR )

      CALL ApplyPositivityLimiter_Radiation &
             ( iX_B0 = iX_B0, iX_E0 = iX_E0, &
               iX_B1 = iX_B1, iX_E1 = iX_E1, &
               U = uCR )

    END IF

    CALL Finalize_SSP_RK

  END SUBROUTINE SSP_RK1


  SUBROUTINE SSP_RK2( t, dt )

    REAL(DP), INTENT(in) :: t, dt

    CALL Initialize_SSP_RK

    ! -- RK Stage 1 --

    IF( EvolveFluid )THEN

      CALL ApplyBoundaryConditions_Fluid

      CALL ComputeRHS_Fluid &
             ( iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ] )

    END IF

    IF( EvolveRadiation )THEN

      CALL ApplyBoundaryConditions_Radiation( t )

      CALL ComputeRHS_Radiation &
             ( iX_B0 = iX_B0, iX_E0 = iX_E0, &
               iX_B1 = iX_B1, iX_E1 = iX_E1, &
               U = uCR, dU = rhsCR )

    END IF

    IF( EvolveFluid )THEN

      CALL ApplyRHS_Fluid &
             ( iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ], &
               dt = dt, alpha = 0.0_DP, beta = 1.0_DP )

      CALL ApplyBoundaryConditions_Fluid

      CALL ApplySlopeLimiter_Fluid

      CALL ApplyPositivityLimiter_Fluid

    END IF

    IF( EvolveGravity )THEN

      CALL SolveGravity

    END IF

    IF( EvolveRadiation )THEN

      CALL ApplyRHS_Radiation &
             ( iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ], &
               dt = dt, alpha = 0.0_DP, beta = 1.0_DP )

      CALL ApplyBoundaryConditions_Radiation &
             ( t + dt, LimiterBC_Option = .TRUE. )

      CALL ApplySlopeLimiter_Radiation &
             ( iX_B0 = iX_B0, iX_E0 = iX_E0, &
               iX_B1 = iX_B1, iX_E1 = iX_E1, &
               U = uCR )

      CALL ApplyPositivityLimiter_Radiation &
             ( iX_B0 = iX_B0, iX_E0 = iX_E0, &
               iX_B1 = iX_B1, iX_E1 = iX_E1, &
               U = uCR )

    END IF

    ! -- RK Stage 2 --

    IF( EvolveFluid )THEN

      CALL ApplyBoundaryConditions_Fluid

      CALL ComputeRHS_Fluid &
             ( iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ] )

    END IF

    IF( EvolveRadiation )THEN

      CALL ApplyBoundaryConditions_Radiation( t + dt )

      CALL ComputeRHS_Radiation &
             ( iX_B0 = iX_B0, iX_E0 = iX_E0, &
               iX_B1 = iX_B1, iX_E1 = iX_E1, &
               U = uCR, dU = rhsCR )

    END IF

    IF( EvolveFluid )THEN

      CALL ApplyRHS_Fluid &
             ( iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ], &
               dt = dt, alpha = 0.5_DP, beta = 0.5_DP )

      CALL ApplyBoundaryConditions_Fluid

      CALL ApplySlopeLimiter_Fluid

      CALL ApplyPositivityLimiter_Fluid

    END IF

    IF( EvolveGravity )THEN

      CALL SolveGravity

    END IF

    IF( EvolveRadiation )THEN

      CALL ApplyRHS_Radiation &
             ( iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ], &
               dt = dt, alpha = 0.5_DP, beta = 0.5_DP )

      CALL ApplyBoundaryConditions_Radiation &
             ( t + dt, LimiterBC_Option = .TRUE. )

      CALL ApplySlopeLimiter_Radiation &
             ( iX_B0 = iX_B0, iX_E0 = iX_E0, &
               iX_B1 = iX_B1, iX_E1 = iX_E1, &
               U = uCR )

      CALL ApplyPositivityLimiter_Radiation &
             ( iX_B0 = iX_B0, iX_E0 = iX_E0, &
               iX_B1 = iX_B1, iX_E1 = iX_E1, &
               U = uCR )

    END IF

    CALL Finalize_SSP_RK

  END SUBROUTINE SSP_RK2


  SUBROUTINE SSP_RK3( t, dt )

    REAL(DP), INTENT(in) :: t, dt

    CALL Initialize_SSP_RK

    ! -- RK Stage 1 -- 

    IF( EvolveFluid )THEN

      CALL ApplyBoundaryConditions_Fluid

      CALL ComputeRHS_Fluid &
             ( iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ] )

    END IF

    IF( EvolveRadiation )THEN

      CALL ApplyBoundaryConditions_Radiation( t )

      CALL ComputeRHS_Radiation &
             ( iX_B0 = iX_B0, iX_E0 = iX_E0, &
               iX_B1 = iX_B1, iX_E1 = iX_E1, &
               U = uCR, dU = rhsCR )

    END IF

    IF( EvolveFluid )THEN

      CALL ApplyRHS_Fluid &
             ( iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ], &
               dt = dt, alpha = 0.0_DP, beta = 1.0_DP )

      CALL ApplyBoundaryConditions_Fluid

      CALL ApplySlopeLimiter_Fluid

      CALL ApplyPositivityLimiter_Fluid

    END IF

    IF( EvolveGravity )THEN

      CALL SolveGravity

    END IF

    IF( EvolveRadiation )THEN

      CALL ApplyRHS_Radiation &
             ( iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ], &
               dt = dt, alpha = 0.0_DP, beta = 1.0_DP )

      CALL ApplyBoundaryConditions_Radiation &
             ( t + dt, LimiterBC_Option = .TRUE. )

      CALL ApplySlopeLimiter_Radiation &
             ( iX_B0 = iX_B0, iX_E0 = iX_E0, &
               iX_B1 = iX_B1, iX_E1 = iX_E1, &
               U = uCR )

      CALL ApplyPositivityLimiter_Radiation &
             ( iX_B0 = iX_B0, iX_E0 = iX_E0, &
               iX_B1 = iX_B1, iX_E1 = iX_E1, &
               U = uCR )

    END IF

    ! -- RK Stage 2 --

    IF( EvolveFluid )THEN

      CALL ApplyBoundaryConditions_Fluid

      CALL ComputeRHS_Fluid &
             ( iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ] )

    END IF

    IF( EvolveRadiation )THEN

      CALL ApplyBoundaryConditions_Radiation( t + dt )

      CALL ComputeRHS_Radiation &
             ( iX_B0 = iX_B0, iX_E0 = iX_E0, &
               iX_B1 = iX_B1, iX_E1 = iX_E1, &
               U = uCR, dU = rhsCR )

    END IF

    IF( EvolveFluid )THEN

      CALL ApplyRHS_Fluid &
             ( iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ], &
               dt = dt, alpha = 0.75_DP, beta = 0.25_DP )

      CALL ApplyBoundaryConditions_Fluid

      CALL ApplySlopeLimiter_Fluid

      CALL ApplyPositivityLimiter_Fluid

    END IF

    IF( EvolveGravity )THEN

      CALL SolveGravity

    END IF

    IF( EvolveRadiation )THEN

      CALL ApplyRHS_Radiation &
             ( iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ], &
               dt = dt, alpha = 0.75_DP, beta = 0.25_DP )

      CALL ApplyBoundaryConditions_Radiation &
             ( t + 0.5_DP * dt, LimiterBC_Option = .TRUE. )

      CALL ApplySlopeLimiter_Radiation &
             ( iX_B0 = iX_B0, iX_E0 = iX_E0, &
               iX_B1 = iX_B1, iX_E1 = iX_E1, &
               U = uCR )

      CALL ApplyPositivityLimiter_Radiation &
             ( iX_B0 = iX_B0, iX_E0 = iX_E0, &
               iX_B1 = iX_B1, iX_E1 = iX_E1, &
               U = uCR )

    END IF

    ! -- RK Stage 3 --

    IF( EvolveFluid )THEN

      CALL ApplyBoundaryConditions_Fluid

      CALL ComputeRHS_Fluid &
             ( iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ] )

    END IF

    IF( EvolveRadiation )THEN

      CALL ApplyBoundaryConditions_Radiation( t + 0.5_DP * dt )

      CALL ComputeRHS_Radiation &
             ( iX_B0 = iX_B0, iX_E0 = iX_E0, &
               iX_B1 = iX_B1, iX_E1 = iX_E1, &
               U = uCR, dU = rhsCR )

    END IF

    IF( EvolveFluid )THEN

      CALL ApplyRHS_Fluid &
             ( iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ], &
               dt = dt, alpha = 1.0_DP / 3.0_DP, beta = 2.0_DP / 3.0_DP )

      CALL ApplyBoundaryConditions_Fluid

      CALL ApplySlopeLimiter_Fluid

      CALL ApplyPositivityLimiter_Fluid

    END IF

    IF( EvolveGravity )THEN

      CALL SolveGravity

    END IF

    IF( EvolveRadiation )THEN

      CALL ApplyRHS_Radiation &
             ( iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ], &
               dt = dt, alpha = 1.0_DP / 3.0_DP , beta = 2.0_DP / 3.0_DP )

      CALL ApplyBoundaryConditions_Radiation &
             ( t + dt, LimiterBC_Option = .TRUE. )

      CALL ApplySlopeLimiter_Radiation &
             ( iX_B0 = iX_B0, iX_E0 = iX_E0, &
               iX_B1 = iX_B1, iX_E1 = iX_E1, &
               U = uCR )

      CALL ApplyPositivityLimiter_Radiation &
             ( iX_B0 = iX_B0, iX_E0 = iX_E0, &
               iX_B1 = iX_B1, iX_E1 = iX_E1, &
               U = uCR )

    END IF

    CALL Finalize_SSP_RK

  END SUBROUTINE SSP_RK3


  SUBROUTINE Initialize_SSP_RK

    IF( EvolveFluid )THEN

      ALLOCATE( uCF_0(1:nDOFX,1:nX(1),1:nX(2),1:nX(3),1:nCF) )

      uCF_0(1:nDOFX,1:nX(1),1:nX(2),1:nX(3),1:nCF) &
        = uCF(1:nDOFX,1:nX(1),1:nX(2),1:nX(3),1:nCF)

    END IF

    IF( EvolveRadiation )THEN

      ALLOCATE( uCR_0(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),1:nCR,1:nSpecies) )

      uCR_0(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),1:nCR,1:nSpecies) &
        = uCR(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),1:nCR,1:nSpecies)

    END IF

  END SUBROUTINE Initialize_SSP_RK


  SUBROUTINE Finalize_SSP_RK

    IF( EvolveFluid )THEN

      DEALLOCATE( uCF_0 )

    END IF

    IF( EvolveRadiation )THEN

      DEALLOCATE( uCR_0 )

    END IF

  END SUBROUTINE Finalize_SSP_RK


  SUBROUTINE SI_RK1( t, dt )

    REAL(DP), INTENT(in) :: t, dt

    CALL Initialize_SI_RK

    IF( EvolveRadiation )THEN

      CALL ApplyBoundaryConditions_Radiation( t )

      CALL ComputeRHS_Radiation &
             ( iX_B0 = iX_B0, iX_E0 = iX_E0, &
               iX_B1 = iX_B1, iX_E1 = iX_E1, &
               U = uCR, dU = rhsCR )

    END IF

    IF( EvolveRadiation )THEN

      CALL ApplyRHS_Radiation &
             ( iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ], &
               dt = dt, alpha = 0.0_DP, beta = 1.0_DP )

      CALL ApplyBoundaryConditions_Radiation &
             ( t + dt, LimiterBC_Option = .TRUE. )

      CALL ApplySlopeLimiter_Radiation &
             ( iX_B0 = iX_B0, iX_E0 = iX_E0, &
               iX_B1 = iX_B1, iX_E1 = iX_E1, &
               U = uCR )

      CALL ApplyPositivityLimiter_Radiation &
             ( iX_B0 = iX_B0, iX_E0 = iX_E0, &
               iX_B1 = iX_B1, iX_E1 = iX_E1, &
               U = uCR )

    END IF

    CALL CoupleFluidRadiation &
           ( dt = dt, &
             iX_B0 = iX_B0, iX_E0 = iX_E0, &
             iX_B1 = iX_B1, iX_E1 = iX_E1, &
             U_F = uCF, dU_F = rhsCF, &
             U_R = uCR, dU_R = rhsCR, &
             EvolveFluid_Option = EvolveFluid )

    CALL Finalize_SI_RK

  END SUBROUTINE SI_RK1


  SUBROUTINE SI_RK2( t, dt )

    REAL(DP), INTENT(in) :: t, dt

    INTEGER :: iE, iX1, iX2, iX3, iCF, iCR, iS

    CALL Initialize_SI_RK

    ! -- SI-RK Stage 1 --

    IF( EvolveRadiation )THEN

      CALL ApplyBoundaryConditions_Radiation( t )

      CALL ComputeRHS_Radiation &
             ( iX_B0 = iX_B0, iX_E0 = iX_E0, &
               iX_B1 = iX_B1, iX_E1 = iX_E1, &
               U = uCR, dU = rhsCR )

    END IF

    IF( EvolveRadiation )THEN

      CALL ApplyRHS_Radiation &
             ( iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ], &
               dt = dt, alpha = 0.0_DP, beta = 1.0_DP )

      CALL ApplyBoundaryConditions_Radiation &
             ( t + dt, LimiterBC_Option = .TRUE. )

      CALL ApplySlopeLimiter_Radiation &
             ( iX_B0 = iX_B0, iX_E0 = iX_E0, &
               iX_B1 = iX_B1, iX_E1 = iX_E1, &
               U = uCR )

      CALL ApplyPositivityLimiter_Radiation &
             ( iX_B0 = iX_B0, iX_E0 = iX_E0, &
               iX_B1 = iX_B1, iX_E1 = iX_E1, &
               U = uCR )

    END IF

    CALL CoupleFluidRadiation &
           ( dt = dt, &
             iX_B0 = iX_B0, iX_E0 = iX_E0, &
             iX_B1 = iX_B1, iX_E1 = iX_E1, &
             U_F = uCF, dU_F = rhsCF, &
             U_R = uCR, dU_R = rhsCR, &
             EvolveFluid_Option = EvolveFluid )

    IF( EvolveRadiation )THEN

      CALL ApplyPositivityLimiter_Radiation &
             ( iX_B0 = iX_B0, iX_E0 = iX_E0, &
               iX_B1 = iX_B1, iX_E1 = iX_E1, &
               U = uCR )

    END IF

    ! -- SI-RK Stage 2 --

    IF( EvolveRadiation )THEN

      CALL ApplyBoundaryConditions_Radiation( t + dt )

      CALL ComputeRHS_Radiation &
             ( iX_B0 = iX_B0, iX_E0 = iX_E0, &
               iX_B1 = iX_B1, iX_E1 = iX_E1, &
               U = uCR, dU = rhsCR )

    END IF

    IF( EvolveRadiation )THEN

      CALL ApplyRHS_Radiation &
             ( iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ], &
               dt = dt, alpha = 0.0_DP, beta = 1.0_DP )

      CALL ApplyBoundaryConditions_Radiation &
             ( t + dt, LimiterBC_Option = .TRUE. )

      CALL ApplySlopeLimiter_Radiation &
             ( iX_B0 = iX_B0, iX_E0 = iX_E0, &
               iX_B1 = iX_B1, iX_E1 = iX_E1, &
               U = uCR )

      CALL ApplyPositivityLimiter_Radiation &
             ( iX_B0 = iX_B0, iX_E0 = iX_E0, &
               iX_B1 = iX_B1, iX_E1 = iX_E1, &
               U = uCR )

    END IF

    CALL CoupleFluidRadiation &
           ( dt = dt, &
             iX_B0 = iX_B0, iX_E0 = iX_E0, &
             iX_B1 = iX_B1, iX_E1 = iX_E1, &
             U_F = uCF, dU_F = rhsCF, &
             U_R = uCR, dU_R = rhsCR, &
             EvolveFluid_Option = EvolveFluid )

    IF( EvolveRadiation )THEN

      CALL ApplyPositivityLimiter_Radiation &
             ( iX_B0 = iX_B0, iX_E0 = iX_E0, &
               iX_B1 = iX_B1, iX_E1 = iX_E1, &
               U = uCR )

    END IF

    ! -- Combine Steps --

    IF( EvolveFluid )THEN

      DO iCF = 1, nCF
        DO iX3 = 1, nX(3)
          DO iX2 = 1, nX(2)
            DO iX1 = 1, nX(1)

              uCF(:,iX1,iX2,iX3,iCF) &
                = 0.5_DP * ( uCF_0(:,iX1,iX2,iX3,iCF) &
                             + uCF(:,iX1,iX2,iX3,iCF) )

            END DO
          END DO
        END DO
      END DO

    END IF

    IF( EvolveRadiation )THEN

      DO iS = 1, nSpecies
        DO iCR = 1, nCR
          DO iX3 = 1, nX(3)
            DO iX2 = 1, nX(2)
              DO iX1 = 1, nX(1)
                DO iE = 1, nE

                  uCR(:,iE,iX1,iX2,iX3,iCR,iS) &
                    = 0.5_DP * ( uCR_0(:,iE,iX1,iX2,iX3,iCR,iS) &
                                 + uCR(:,iE,iX1,iX2,iX3,iCR,iS) )

                END DO
              END DO
            END DO
          END DO
        END DO
      END DO

    END IF

    CALL Finalize_SI_RK

  END SUBROUTINE SI_RK2


  SUBROUTINE Initialize_SI_RK

    IF( EvolveFluid )THEN

      ALLOCATE( uCF_0(1:nDOFX,1:nX(1),1:nX(2),1:nX(3),1:nCF) )

      uCF_0(1:nDOFX,1:nX(1),1:nX(2),1:nX(3),1:nCF) &
        = uCF(1:nDOFX,1:nX(1),1:nX(2),1:nX(3),1:nCF)

    END IF

    IF( EvolveRadiation )THEN

      ALLOCATE( uCR_0(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),1:nCR,1:nSpecies) )

      uCR_0(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),1:nCR,1:nSpecies) &
        = uCR(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),1:nCR,1:nSpecies)

    END IF

  END SUBROUTINE Initialize_SI_RK


  SUBROUTINE Finalize_SI_RK

    IF( EvolveFluid )THEN

      DEALLOCATE( uCF_0 )

    END IF

    IF( EvolveRadiation )THEN

      DEALLOCATE( uCR_0 )

    END IF

  END SUBROUTINE Finalize_SI_RK


  SUBROUTINE ApplyRHS_Fluid( iX_Begin, iX_End, dt, alpha, beta )

    INTEGER, DIMENSION(3), INTENT(in) :: iX_Begin, iX_End
    REAL(DP),              INTENT(in) :: dt, alpha, beta

    INTEGER :: iX1, iX2, iX3, iCF

    DO iCF = 1, nCF
      DO iX3 = iX_Begin(3), iX_End(3)
        DO iX2 = iX_Begin(2), iX_End(2)
          DO iX1 = iX_Begin(1), iX_End(1)

            uCF(:,iX1,iX2,iX3,iCF) &
              = alpha * uCF_0(:,iX1,iX2,iX3,iCF) &
                  + beta * ( uCF(:,iX1,iX2,iX3,iCF) &
                             + dt * rhsCF(:,iX1,iX2,iX3,iCF) )

          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE ApplyRHS_Fluid


  SUBROUTINE ApplyRHS_Radiation( iX_Begin, iX_End, dt, alpha, beta )

    INTEGER, DIMENSION(3), INTENT(in) :: iX_Begin, iX_End
    REAL(DP),              INTENT(in) :: dt, alpha, beta

    INTEGER :: iE, iX1, iX2, iX3, iCR, iS

    DO iS = 1, nSpecies
      DO iCR = 1, nCR
        DO iX3 = iX_Begin(3), iX_End(3)
          DO iX2 = iX_Begin(2), iX_End(2)
            DO iX1 = iX_Begin(1), iX_End(1)
              DO iE = 1, nE

                uCR(:,iE,iX1,iX2,iX3,iCR,iS) &
                  = alpha * uCR_0(:,iE,iX1,iX2,iX3,iCR,iS) &
                      + beta * ( uCR(:,iE,iX1,iX2,iX3,iCR,iS) &
                                 + dt * rhsCR(:,iE,iX1,iX2,iX3,iCR,iS) )

              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE ApplyRHS_Radiation


  SUBROUTINE BackwardEuler( t, dt )

    REAL(DP), INTENT(in) :: t, dt

    CALL CoupleFluidRadiation &
           ( dt = dt, &
             iX_B0 = iX_B0, iX_E0 = iX_E0, &
             iX_B1 = iX_B1, iX_E1 = iX_E1, &
             U_F = uCF, dU_F = rhsCF, &
             U_R = uCR, dU_R = rhsCR, &
             EvolveFluid_Option = EvolveFluid )

  END SUBROUTINE BackwardEuler


  SUBROUTINE InitializeTimeStepper_SSPRK( nStages )

    INTEGER, INTENT(in) :: nStages

    INTEGER :: i

    IF( nStages == 0 )THEN

      SSPRKs => TimeStepper_DUMMY
      RETURN

    ELSE

      CALL InitializeSSPRK_Scheme( nStages )

      WRITE(*,*)
      WRITE(*,'(A5,A,I1)') '', 'SSP RK Scheme: ', nStages

      WRITE(*,*)
      WRITE(*,'(A5,A)') '', 'Butcher Table:'
      WRITE(*,'(A5,A)') '', '--------------'
      DO i = 1, nStages
        WRITE(*,'(A5,4ES14.4E3)') '', c_SSP(i), a_SSP(i,1:nStages)
      END DO
      WRITE(*,'(A5,A14,3ES14.4E3)') '', '', w_SSP(1:nStages)

      nDOF_SSP = 0

      IF( EvolveRadiation )THEN

        nDOF_SSP &
          = nDOF_SSP &
              + nDOF * nE * PRODUCT( nX ) * nCR * nSpecies

      END IF

      IF( EvolveFluid )THEN

        nDOF_SSP &
          = nDOF_SSP &
              + nDOFX * PRODUCT( nX ) * nCF

      END IF

      ALLOCATE( U_SSP(nDOF_SSP) )
      ALLOCATE( F_SSP(nDOF_SSP,nStages) )

      SSPRKs => SSPRKs_Radiation

    END IF

  END SUBROUTINE InitializeTimeStepper_SSPRK


  SUBROUTINE InitializeSSPRK_Scheme( nStages )

    INTEGER, INTENT(in) :: nStages

    INTEGER :: iS

    CALL AllocateButcherTables_SSPRK( nStages )

    SELECT CASE ( nStages )
      CASE ( 1 )

        a_SSP(1,1) = 0.0_DP
        w_SSP(1)   = 1.0_DP

      CASE ( 2 )

        a_SSP(1,1:2) = [ 0.0_DP, 0.0_DP ]
        a_SSP(2,1:2) = [ 1.0_DP, 0.0_DP ]
        w_SSP(1:2)   = [ 0.5_DP, 0.5_DP ]

      CASE ( 3 )

        a_SSP(1,1:3) = [ 0.00_DP, 0.00_DP, 0.00_DP ]
        a_SSP(2,1:3) = [ 1.00_DP, 0.00_DP, 0.00_DP ]
        a_SSP(3,1:3) = [ 0.25_DP, 0.25_DP, 0.00_DP ]
        w_SSP(1:3)   = [ 1.0_DP / 6.0_DP, &
                         1.0_DP / 6.0_DP, &
                         2.0_DP / 3.0_DP ]

    END SELECT

    DO iS = 1, nStages
      c_SSP(iS) = SUM( a_SSP(iS,1:iS-1) )
    END DO

  END SUBROUTINE InitializeSSPRK_Scheme


  SUBROUTINE AllocateButcherTables_SSPRK( nStages )

    INTEGER, INTENT(in) :: nStages

    ALLOCATE( a_SSP(nStages,nStages) )
    ALLOCATE( c_SSP(nStages) )
    ALLOCATE( w_SSP(nStages) )

    a_SSP = Zero
    c_SSP = Zero
    w_SSP = Zero

  END SUBROUTINE AllocateButcherTables_SSPRK


  SUBROUTINE InitializeSSPRK_Step

    ALLOCATE( uCR_0(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),1:nCR,1:nSpecies) )

    uCR_0(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),1:nCR,1:nSpecies) &
      = uCR(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),1:nCR,1:nSpecies)

    U_SSP = Zero
    F_SSP = Zero

  END SUBROUTINE InitializeSSPRK_Step


  SUBROUTINE FinalizeSSPRK_Step

    DEALLOCATE( uCR_0 )

  END SUBROUTINE FinalizeSSPRK_Step


  SUBROUTINE SSPRKs_Radiation( t, dt )

    REAL(DP), INTENT(in) :: t, dt

    INTEGER :: iS, jS

    CALL InitializeSSPRK_Step

    DO iS = 1, nStages_SSP_RK

      CALL MapToStage_Radiation( iX_B0, uCR_0, U_SSP )

      DO jS = 1, iS - 1

        IF( a_SSP(iS,jS) .NE. Zero )THEN

          U_SSP(:) = U_SSP(:) + dt * a_SSP(iS,jS) * F_SSP(:,jS)

        END IF

      END DO

      IF( iS > 1 )THEN

        ! --- Apply Limiters ---

        CALL MapFromStage_Radiation( iX_B1, uCR, U_SSP )

        CALL ApplyBoundaryConditions_Radiation &
               ( t + c_SSP(iS) * dt )

        CALL ApplySlopeLimiter_Radiation &
             ( iX_B0 = iX_B0, iX_E0 = iX_E0, &
               iX_B1 = iX_B1, iX_E1 = iX_E1, &
               U = uCR )

        CALL ApplyPositivityLimiter_Radiation &
             ( iX_B0 = iX_B0, iX_E0 = iX_E0, &
               iX_B1 = iX_B1, iX_E1 = iX_E1, &
               U = uCR )

        CALL MapToStage_Radiation( iX_B1, uCR, U_SSP )

      END IF

      IF( ANY( a_SSP(:,iS) .NE. Zero ) .OR. ( w_SSP(iS) .NE. Zero ) )THEN

        ! --- Explicit Solve ---

        CALL MapFromStage_Radiation( iX_B1, uCR, U_SSP )

        CALL ApplyBoundaryConditions_Radiation &
               ( t + c_SSP(iS) * dt )

        CALL ComputeExplicitIncrement_Radiation &
               ( iX_B0 = iX_B0, iX_E0 = iX_E0, &
                 iX_B1 = iX_B1, iX_E1 = iX_E1, &
                 U = uCR, dU = rhsCR )

        CALL MapToStage_Radiation( iX_B0, rhsCR, F_SSP(:,iS) )

      END IF

    END DO

    CALL MapToStage_Radiation( iX_B0, uCR_0, U_SSP )

    DO iS = 1, nStages_SSP_RK

      IF( w_SSP(iS) .NE. Zero )THEN

        U_SSP(:) = U_SSP(:) + dt * w_SSP(iS) * F_SSP(:,iS)

      END IF

    END DO

    CALL MapFromStage_Radiation( iX_B1, uCR, U_SSP )

    ! --- Apply Limiters ---

    CALL ApplyBoundaryConditions_Radiation &
           ( t + dt )

    CALL ApplySlopeLimiter_Radiation &
           ( iX_B0 = iX_B0, iX_E0 = iX_E0, &
             iX_B1 = iX_B1, iX_E1 = iX_E1, &
             U = uCR )

    CALL ApplyPositivityLimiter_Radiation &
           ( iX_B0 = iX_B0, iX_E0 = iX_E0, &
             iX_B1 = iX_B1, iX_E1 = iX_E1, &
             U = uCR )

    CALL FinalizeSSPRK_Step

  END SUBROUTINE SSPRKs_Radiation


  SUBROUTINE InitializeTimeStepper_IMEX( Scheme )

    CHARACTER(LEN=*), INTENT(in) :: Scheme

    INTEGER :: i

    IF( TRIM( Scheme ) == 'DUMMY' )THEN

      IMEX => TimeStepper_DUMMY
      RETURN

    ELSE

      CALL InitializeIMEX_Scheme( TRIM( Scheme ) )

      WRITE(*,*)
      WRITE(*,'(A5,A,A)') '', 'IMEX Scheme: ', TRIM( Scheme )

      WRITE(*,*)
      WRITE(*,'(A5,A)') '', 'Implicit Butcher Table:'
      WRITE(*,'(A5,A)') '', '-----------------------'
      DO i = 1, nStages_IMEX
        WRITE(*,'(A5,4ES14.4E3)') '', c_IM(i), a_IM(i,1:nStages_IMEX)
      END DO
      WRITE(*,'(A5,A14,3ES14.4E3)') '', '', w_IM(1:nStages_IMEX)

      WRITE(*,*)
      WRITE(*,'(A5,A)') '', 'Explicit Butcher Table:'
      WRITE(*,'(A5,A)') '', '-----------------------'
      DO i = 1, nStages_IMEX
        WRITE(*,'(A5,4ES14.4E3)') '', c_EX(i), a_EX(i,1:nStages_IMEX)
      END DO
      WRITE(*,'(A5,A14,3ES14.4E3)') '', '', w_EX(1:nStages_IMEX)

      nDOF_IMEX = 0

      IF( EvolveRadiation )THEN

        nDOF_IMEX &
          = nDOF_IMEX &
              + nDOF * nE * PRODUCT( nX ) * nCR * nSpecies

      END IF

      IF( EvolveFluid )THEN

        nDOF_IMEX &
          = nDOF_IMEX &
              + nDOFX * PRODUCT( nX ) * nCF

      END IF

      ALLOCATE( U_IMEX(nDOF_IMEX) )
      ALLOCATE( F_IM  (nDOF_IMEX,nStages_IMEX) )
      ALLOCATE( F_EX  (nDOF_IMEX,nStages_IMEX) )

      IMEX => IMEX_Radiation

    END IF

  END SUBROUTINE InitializeTimeStepper_IMEX


  SUBROUTINE InitializeIMEX_Scheme( Scheme )

    CHARACTER(LEN=*), INTENT(in) :: Scheme

    INTEGER :: i

    SELECT CASE( TRIM( Scheme ) )
      CASE ( 'IMEX_SIRK2' )

        nStages_IMEX = 3

        CALL AllocateButcherTables_IMEX( nStages_IMEX )

        a_IM(1,1:3) = [ 0.0_DP, 0.0_DP, 0.0_DP ]
        a_IM(2,1:3) = [ 0.0_DP, 1.0_DP, 0.0_DP ]
        a_IM(3,1:3) = [ 0.0_DP, 1.0_DP, 1.0_DP ]
        w_IM(1:3)   = [ 0.0_DP, 0.5_DP, 0.5_DP ]
        DO i = 1, nStages_IMEX
          c_IM(i) = SUM( a_IM(i,1:i) )
        END DO

        a_EX(1,1:3) = [ 0.0_DP, 0.0_DP, 0.0_DP ]
        a_EX(2,1:3) = [ 1.0_DP, 0.0_DP, 0.0_DP ]
        a_EX(3,1:3) = [ 1.0_DP, 1.0_DP, 0.0_DP ]
        w_EX(1:3)   = [ 0.5_DP, 0.5_DP, 0.0_DP ]
        DO i = 1, nStages_IMEX
          c_EX(i) = SUM( a_EX(i,1:i-1) )
        END DO

      CASE( 'IMEX_2332' )

        nStages_IMEX = 3

        CALL AllocateButcherTables_IMEX( nStages_IMEX )

        ! --- 2332: Pareschi & Russo (2005) ---

        a_IM(1,1:3) = [ Fourth, 0.0_DP, 0.0_DP ]
        a_IM(2,1:3) = [ 0.0_DP, Fourth, 0.0_DP ]
        a_IM(3,1:3) = [ Third,  Third,  Third  ]
        w_IM(1:3)   = [ Third,  Third,  Third  ]
        DO i = 1, nStages_IMEX
          c_IM(i) = SUM( a_IM(i,1:i) )
        END DO

        a_EX(1,1:3) = [ 0.0_DP, 0.0_DP, 0.0_DP ]
        a_EX(2,1:3) = [ 0.5_DP, 0.0_DP, 0.0_DP ]
        a_EX(3,1:3) = [ 0.5_DP, 0.5_DP, 0.0_DP ]
        w_EX(1:3)   = [ Third,  Third,  Third  ]
        DO i = 1, nStages_IMEX
          c_EX(i) = SUM( a_EX(i,1:i-1) )
        END DO

      CASE( 'IMEX_RKCB2' )

        nStages_IMEX = 3

        CALL AllocateButcherTables_IMEX( nStages_IMEX )

        ! --- RKCB2: Cavaglieri & Bewley (2015) ---

        a_IM(1,1:3) = [ 0.0_DP, 0.0_DP,          0.0_DP ]
        a_IM(2,1:3) = [ 0.0_DP, 0.4_DP,          0.0_DP ]
        a_IM(3,1:3) = [ 0.0_DP, 5.0_DP / 6.0_DP, 1.0_DP / 6.0_DP ]
        w_IM(1:3)   = [ 0.0_DP, 5.0_DP / 6.0_DP, 1.0_DP / 6.0_DP ]
        DO i = 1, nStages_IMEX
          c_IM(i) = SUM( a_IM(i,1:i) )
        END DO

        a_EX(1,1:3) = [ 0.0_DP, 0.0_DP,          0.0_DP ]
        a_EX(2,1:3) = [ 0.4_DP, 0.0_DP,          0.0_DP ]
        a_EX(3,1:3) = [ 0.0_DP, 1.0_DP,          0.0_DP ]
        w_EX(1:3)   = [ 0.0_DP, 5.0_DP / 6.0_DP, 1.0_DP / 6.0_DP ]
        DO i = 1, nStages_IMEX
          c_EX(i) = SUM( a_EX(i,1:i-1) )
        END DO

    END SELECT

  END SUBROUTINE InitializeIMEX_Scheme


  SUBROUTINE AllocateButcherTables_IMEX( nStages )

    INTEGER, INTENT(in) :: nStages

    ALLOCATE( a_IM(nStages,nStages) )
    ALLOCATE( a_EX(nStages,nStages) )

    ALLOCATE( c_IM(nStages) )
    ALLOCATE( c_EX(nStages) )

    ALLOCATE( w_IM(nStages) )
    ALLOCATE( w_EX(nStages) )

    a_IM = Zero; a_EX = Zero
    c_IM = Zero; c_EX = Zero
    w_IM = Zero; w_EX = Zero

  END SUBROUTINE AllocateButcherTables_IMEX


  SUBROUTINE InitializeIMEX_Step

    ALLOCATE( uCR_0(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),1:nCR,1:nSpecies) )

    uCR_0(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),1:nCR,1:nSpecies) &
      = uCR(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),1:nCR,1:nSpecies)

    U_IMEX = Zero
    F_IM   = Zero
    F_EX   = Zero

  END SUBROUTINE InitializeIMEX_Step


  SUBROUTINE FinalizeIMEX_Step

    DEALLOCATE( uCR_0 )

  END SUBROUTINE FinalizeIMEX_Step


  SUBROUTINE IMEX_Radiation( t, dt )

    REAL(DP), INTENT(in) :: t, dt

    INTEGER :: iS, jS

    CALL InitializeIMEX_Step

    DO iS = 1, nStages_IMEX

      CALL MapToStage_Radiation( iX_B0, uCR_0, U_IMEX )

      DO jS = 1, iS - 1

        IF( a_IM(iS,jS) .NE. Zero )THEN

          U_IMEX(:) = U_IMEX(:) + dt * a_IM(iS,jS) * F_IM(:,jS)

        END IF

        IF( a_EX(iS,jS) .NE. Zero )THEN

          U_IMEX(:) = U_IMEX(:) + dt * a_EX(iS,jS) * F_EX(:,jS)

        END IF

      END DO

      IF( ANY( a_IM(:,iS) .NE. Zero ) .OR. ( w_IM(iS) .NE. Zero ) )THEN

        ! --- Implicit Solve ---

        CALL MapFromStage_Radiation( iX_B1, uCR, U_IMEX )

        CALL ComputeImplicitIncrement_FluidRadiation &
               ( dt = a_IM(iS,iS) * dt, &
                 iX_B0 = iX_B0, iX_E0 = iX_E0, &
                 iX_B1 = iX_B1, iX_E1 = iX_E1, &
                 U_F = uCF, dU_F = rhsCF, &
                 U_R = uCR, dU_R = rhsCR )

        CALL MapToStage_Radiation( iX_B0, rhsCR, F_IM(:,iS) )

        U_IMEX(:) = U_IMEX(:) + dt * a_IM(iS,iS) * F_IM(:,iS)

      END IF

      IF( ANY( a_EX(:,iS) .NE. Zero ) .OR. ( w_EX(iS) .NE. Zero ) )THEN

        ! --- Explicit Solve ---

        CALL MapFromStage_Radiation( iX_B1, uCR, U_IMEX )

        CALL ApplyBoundaryConditions_Radiation &
               ( t + c_EX(iS) * dt )

        CALL ComputeExplicitIncrement_Radiation &
               ( iX_B0 = iX_B0, iX_E0 = iX_E0, &
                 iX_B1 = iX_B1, iX_E1 = iX_E1, &
                 U = uCR, dU = rhsCR )

        CALL MapToStage_Radiation( iX_B0, rhsCR, F_EX(:,iS) )

      END IF

    END DO

    CALL MapToStage_Radiation( iX_B0, uCR_0, U_IMEX )

    DO iS = 1, nStages_IMEX

      IF( w_IM(iS) .NE. Zero )THEN

        U_IMEX(:) = U_IMEX(:) + dt * w_IM(iS) * F_IM(:,iS)

      END IF

      IF( w_EX(iS) .NE. Zero )THEN

        U_IMEX(:) = U_IMEX(:) + dt * w_EX(iS) * F_EX(:,iS)

      END IF

    END DO

    CALL MapFromStage_Radiation( iX_B1, uCR, U_IMEX )

    CALL FinalizeIMEX_Step

  END SUBROUTINE IMEX_Radiation


  SUBROUTINE MapToStage_Radiation( iX_B, U_R, U )

    INTEGER,  INTENT(in)  :: iX_B(3)
    REAL(DP), INTENT(in)  :: U_R(:,:,iX_B(1):,iX_B(2):,iX_B(3):,:,:)
    REAL(DP), INTENT(out) :: U(:)

    INTEGER :: i, iNode, iE, iX1, iX2, iX3, iCR, iS

    i = 1

    DO iS = 1, nSpecies
      DO iCR = 1, nCR
        DO iX3 = 1, nX(3)
          DO iX2 = 1, nX(2)
            DO iX1 = 1, nX(1)
              DO iE = 1, nE
                DO iNode = 1, nDOF

                  U(i) = U_R(iNode,iE,iX1,iX2,iX3,iCR,iS)

                  i = i + 1

                END DO
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE MapToStage_Radiation


  SUBROUTINE MapFromStage_Radiation( iX_B, U_R, U )

    INTEGER,  INTENT(in)  :: iX_B(3)
    REAL(DP), INTENT(out) :: U_R(:,:,iX_B(1):,iX_B(2):,iX_B(3):,:,:)
    REAL(DP), INTENT(in)  :: U(:)

    INTEGER :: i, iNode, iE, iX1, iX2, iX3, iCR, iS

    i = 1

    DO iS = 1, nSpecies
      DO iCR = 1, nCR
        DO iX3 = 1, nX(3)
          DO iX2 = 1, nX(2)
            DO iX1 = 1, nX(1)
              DO iE = 1, nE
                DO iNode = 1, nDOF

                  U_R(iNode,iE,iX1,iX2,iX3,iCR,iS) = U(i)

                  i = i + 1

                END DO
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE MapFromStage_Radiation


END MODULE TimeSteppingModule
