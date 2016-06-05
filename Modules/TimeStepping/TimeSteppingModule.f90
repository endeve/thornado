MODULE TimeSteppingModule

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    nX, nDOFX, nNodes
  USE MeshModule, ONLY: &
    MeshX
  USE FluidFieldsModule, ONLY: &
    uCF, rhsCF, nCF, nPF, nAF, iPF_V1, iPF_V2, iPF_V3, iAF_Cs
  USE EquationOfStateModule, ONLY: &
    Auxiliary_Fluid
  USE InputOutputModule, ONLY: &
    WriteFields1D
  USE EulerEquationsUtilitiesModule, ONLY: &
    Primitive, &
    Eigenvalues
  USE FluidEvolutionModule, ONLY: &
    ComputeRHS_Fluid, &
    ApplySlopeLimiter_Fluid, &
    ApplyPositivityLimiter_Fluid
  USE BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_Fluid

  IMPLICIT NONE
  PRIVATE

  REAL(DP) :: wtR, wtS, wtP

  LOGICAL :: EvolveFluid     = .FALSE.
  LOGICAL :: EvolveRadiation = .FALSE.
  INTEGER :: nStagesSSPRK    = 1
  REAL(DP), DIMENSION(:,:,:,:,:), ALLOCATABLE :: uCF_0

  PROCEDURE (TimeStepExplicit), POINTER, PUBLIC :: &
    SSP_RK => NULL()

  INTERFACE
    SUBROUTINE TimeStepExplicit( dt )
      USE KindModule, ONLY: DP
      REAL(DP), INTENT(in) :: dt
    END SUBROUTINE TimeStepExplicit
  END INTERFACE

  PUBLIC :: InitializeTimeStepping
  PUBLIC :: FinalizeTimeStepping
  PUBLIC :: EvolveFields

CONTAINS


  SUBROUTINE InitializeTimeStepping &
               ( EvolveFluid_Option, EvolveRadiation_Option, &
                 nStagesSSPRK_Option )

    LOGICAL, INTENT(in), OPTIONAL :: EvolveFluid_Option
    LOGICAL, INTENT(in), OPTIONAL :: EvolveRadiation_Option
    INTEGER, INTENT(in), OPTIONAL :: nStagesSSPRK_Option

    IF( PRESENT( EvolveFluid_Option ) )THEN
      EvolveFluid = EvolveFluid_Option
    END IF

    IF( PRESENT( EvolveRadiation_Option ) )THEN
      EvolveRadiation = EvolveRadiation_Option
    END IF

    IF( PRESENT( nStagesSSPRK_Option ) )THEN
      nStagesSSPRK = nStagesSSPRK_Option
    END IF

    SELECT CASE ( nStagesSSPRK )
      CASE ( 1 )

        SSP_RK => SSP_RK1

      CASE ( 2 )

        SSP_RK => SSP_RK2

      CASE ( 3 )

        SSP_RK => SSP_RK3

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,'(A4,A43,I2.2)') &
          '', 'SSP_RK not implemented for nStatgesSSPRK = ', nStagesSSPRK
        STOP

    END SELECT

  END SUBROUTINE InitializeTimeStepping


  SUBROUTINE FinalizeTimeStepping

    NULLIFY( SSP_RK )

  END SUBROUTINE FinalizeTimeStepping


  SUBROUTINE EvolveFields( t_begin, t_end, dt_write, UpdateFields )

    REAL(DP), INTENT(in) :: t_begin, t_end, dt_write
    INTERFACE
      SUBROUTINE UpdateFields( dt )
        USE KindModule, ONLY: DP
        REAL(DP), INTENT(in) :: dt
      END SUBROUTINE UpdateFields
    END INTERFACE

    LOGICAL  :: WriteOutput = .FALSE.
    INTEGER  :: iCycle
    REAL(DP) :: t, t_write, dt
    REAL(DP), DIMENSION(0:1) :: WallTime

    WRITE(*,*)
    WRITE(*,'(A4,A21)') '', 'INFO: Evolving Fields'
    WRITE(*,*)
    WRITE(*,'(A6,A10,ES10.4E2,A2,A8,ES10.4E2,A2,A11,ES10.4E2)') &
      '  ', 't_begin = ', t_begin, ', ', 't_end = ', t_end, &
      ', ', 'dt_write = ', dt_write
    WRITE(*,*)

    iCycle  = 0
    t       = t_begin
    t_write = dt_write

    CALL WriteFields1D

    CALL CPU_TIME( WallTime(0) )

    wtR = 0.0_DP; wtS = 0.0_DP; wtP = 0.0_DP

    DO WHILE( t < t_end )

      iCycle = iCycle + 1

      CALL ComputeTimeStep( dt )

      IF( t + dt > t_end )THEN

        dt = t_end - t

      END IF

      IF( t + dt > t_write )THEN

        dt          = t_write - t
        t_write     = t_write + dt_write
        WriteOutput = .TRUE.

      END IF

      IF( MOD( iCycle, 100 ) == 0 )THEN

        WRITE(*,'(A8,A8,I8.8,A2,A4,ES10.4E2,A2,A5,ES10.4E2)') &
          '', 'Cycle = ', iCycle, ' ', 't = ', t, ' ', 'dt = ', dt

      END IF

      CALL UpdateFields( dt )

      t = t + dt

      IF( WriteOutput )THEN

        CALL WriteFields1D

        WriteOutput = .FALSE.

      END IF

    END DO

    CALL CPU_TIME( WallTime(1) )

    WRITE(*,*)
    WRITE(*,'(A6,A15,ES10.4E2,A6,I8.8,A7,A4,ES10.4E2,A2)') &
      '', 'Evolved to t = ', t, ' with ', iCycle, ' cycles', &
      ' in ', WallTime(1)-WallTime(0), ' s'
    WRITE(*,*)
    WRITE(*,'(A6,A12,ES10.4E2,A2,A12,ES10.4E2,A2,A12,ES10.4E2)') &
      '', '  wt(RHS) = ', wtR / (WallTime(1)-WallTime(0)), &
      '', '  wt(SLM) = ', wtS / (WallTime(1)-WallTime(0)), &
      '', '  wt(PLM) = ', wtP / (WallTime(1)-WallTime(0))
    WRITE(*,*)

    CALL WriteFields1D

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
    REAL(DP), DIMENSION(1:nAF) :: uAF_N

    dt = HUGE( 1.0_DP )

    ASSOCIATE( dX1 => MeshX(1) % Width(1:nX(1)), &
               dX2 => MeshX(2) % Width(1:nX(2)), &
               dX3 => MeshX(3) % Width(1:nX(3)) )

    CFL = 0.2_DP / ( 2.0_DP * DBLE( nNodes - 1 ) + 1.0_DP ) ! For Debugging

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)
          DO iNodeX = 1, nDOFX

            uPF_N = Primitive( uCF(iNodeX,iX1,ix2,ix3,1:nCF) )

            uAF_N = Auxiliary_Fluid( uPF_N )

            dt_X1 &
              = CFL * dX1(iX1) &
                  / MAXVAL(ABS(Eigenvalues( uPF_N(iPF_V1), uAF_N(iAF_Cs) )))

            dt_X2 &
              = CFL * dX2(iX2) &
                  / MAXVAL(ABS(Eigenvalues( uPF_N(iPF_V2), uAF_N(iAF_Cs) )))

            dt_X3 &
              = CFL * dX3(iX3) &
                  / MAXVAL(ABS(Eigenvalues( uPF_N(iPF_V3), uAF_N(iAF_Cs) )))

            dt = MIN( dt, dt_x1, dt_X2, dt_X3 )

          END DO
        END DO
      END DO
    END DO

    END ASSOCIATE ! dX1

  END SUBROUTINE ComputeTimeStep_Fluid


  SUBROUTINE ComputeTimeStep_Radiation( dt )

    REAL(DP), INTENT(out) :: dt

  END SUBROUTINE ComputeTimeStep_Radiation


  SUBROUTINE SSP_RK1( dt )

    REAL(DP), INTENT(in) :: dt

    IF( EvolveFluid )THEN

      CALL Initialize_SSP_RK

      CALL ApplyBoundaryConditions_Fluid

      CALL ComputeRHS_Fluid &
             ( iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ] )

      CALL ApplyRHS_Fluid &
             ( iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ], &
               dt = dt, a = 0.0_DP, b = 1.0_DP )

      CALL ApplySlopeLimiter_Fluid

      CALL ApplyPositivityLimiter_Fluid

      CALL Finalize_SSP_RK

    END IF

  END SUBROUTINE SSP_RK1


  SUBROUTINE SSP_RK2( dt )

    REAL(DP), INTENT(in) :: dt

    REAL(DP), DIMENSION(0:1) :: WT

    IF( EvolveFluid )THEN

      CALL Initialize_SSP_RK

      CALL ApplyBoundaryConditions_Fluid

      CALL CPU_TIME( WT(0) )

      CALL ComputeRHS_Fluid &
             ( iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ] )

      CALL CPU_TIME( WT(1) )
      wtR = wtR + ( WT(1) - WT(0) )

      CALL ApplyRHS_Fluid &
             ( iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ], &
               dt = dt, a = 0.0_DP, b = 1.0_DP )

      CALL CPU_TIME( WT(0) )

      CALL ApplySlopeLimiter_Fluid

      CALL CPU_TIME( WT(1) )
      wtS = wtS + ( WT(1) - WT(0) )

      CALL CPU_TIME( WT(0) )

      CALL ApplyPositivityLimiter_Fluid

      CALL CPU_TIME( WT(1) )
      wtP = wtP + ( WT(1) - WT(0) )

      CALL ApplyBoundaryConditions_Fluid

      CALL CPU_TIME( WT(0) )

      CALL ComputeRHS_Fluid &
             ( iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ] )

      CALL CPU_TIME( WT(1) )
      wtR = wtR + ( WT(1) - WT(0) )

      CALL ApplyRHS_Fluid &
             ( iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ], &
               dt = dt, a = 0.5_DP, b = 0.5_DP )

      CALL CPU_TIME( WT(0) )

      CALL ApplySlopeLimiter_Fluid

      CALL CPU_TIME( WT(1) )
      wtS = wtS + ( WT(1) - WT(0) )

      CALL CPU_TIME( WT(0) )

      CALL ApplyPositivityLimiter_Fluid

      CALL CPU_TIME( WT(1) )
      wtP = wtP + ( WT(1) - WT(0) )

      CALL Finalize_SSP_RK

    END IF

  END SUBROUTINE SSP_RK2


  SUBROUTINE SSP_RK3( dt )

    REAL(DP), INTENT(in) :: dt

    IF( EvolveFluid )THEN

      CALL Initialize_SSP_RK

      CALL ApplyBoundaryConditions_Fluid

      CALL ComputeRHS_Fluid &
             ( iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ] )

      CALL ApplyRHS_Fluid &
             ( iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ], &
               dt = dt, a = 0.0_DP, b = 1.0_DP )

      CALL ApplySlopeLimiter_Fluid

      CALL ApplyPositivityLimiter_Fluid

      CALL ApplyBoundaryConditions_Fluid

      CALL ComputeRHS_Fluid &
             ( iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ] )

      CALL ApplyRHS_Fluid &
             ( iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ], &
               dt = dt, a = 0.75_DP, b = 0.25_DP )

      CALL ApplySlopeLimiter_Fluid

      CALL ApplyPositivityLimiter_Fluid

      CALL ApplyBoundaryConditions_Fluid

      CALL ComputeRHS_Fluid &
             ( iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ] )

      CALL ApplyRHS_Fluid &
             ( iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ], &
               dt = dt, a = 1.0_DP / 3.0_DP, b = 2.0_DP / 3.0_DP )

      CALL ApplySlopeLimiter_Fluid

      CALL ApplyPositivityLimiter_Fluid

      CALL Finalize_SSP_RK

    END IF    

  END SUBROUTINE SSP_RK3


  SUBROUTINE Initialize_SSP_RK

    ALLOCATE( uCF_0(1:nDOFX,1:nX(1),1:nX(2),1:nX(3),1:nCF) )

    uCF_0(1:nDOFX,1:nX(1),1:nX(2),1:nX(3),1:nCF) &
      = uCF(1:nDOFX,1:nX(1),1:nX(2),1:nX(3),1:nCF)

  END SUBROUTINE Initialize_SSP_RK


  SUBROUTINE ApplyRHS_Fluid( iX_Begin, iX_End, dt, a, b )

    INTEGER, DIMENSION(3), INTENT(in) :: iX_Begin, iX_End
    REAL(DP),              INTENT(in) :: dt, a, b

    INTEGER :: iX1, iX2, iX3, iCF

    DO iCF = 1, nCF
      DO iX3 = iX_Begin(3), iX_End(3)
        DO iX2 = iX_Begin(2), iX_End(2)
          DO iX1 = iX_Begin(1), iX_End(1)

            uCF(:,iX1,iX2,iX3,iCF) &
              = a * uCF_0(:,iX1,iX2,iX3,iCF) &
                  + b * ( uCF(:,iX1,iX2,iX3,iCF) &
                          + dt * rhsCF(:,iX1,iX2,iX3,iCF) )

          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE ApplyRHS_Fluid


  SUBROUTINE Finalize_SSP_RK

    DEALLOCATE( uCF_0 )

  END SUBROUTINE Finalize_SSP_RK


END MODULE TimeSteppingModule
