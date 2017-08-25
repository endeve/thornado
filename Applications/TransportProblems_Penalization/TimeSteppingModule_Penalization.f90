MODULE TimeSteppingModule_Penalization

  USE KindModule, ONLY: &
    DP
  USE UnitsModule, ONLY: &
    UnitsDisplay, &
    Millisecond
  USE ProgramHeaderModule, ONLY: &
    nX, nNodesX, nDOFX, &
    nE, nNodesE, nDOFE, &
    nDOF
  USE RadiationFieldsModule, ONLY: &
    iPR_D
  USE InputOutputModule, ONLY: &
    WriteFields1D
  USE FluidRadiationCouplingUtilitiesModule, ONLY: &
    FinalizeFluidFields, &
    FinalizeRadiationFields
  USE FluidRadiationCouplingSolutionModule_Penalization, ONLY: &
    absLambda, rhsCR_C, &
    InitializeFluidRadiationCoupling, &
    FinalizeFluidRadiationCoupling, &
    ComputeRHS_C

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: EvolveFields

CONTAINS


  SUBROUTINE EvolveFields &
               ( t_begin, t_end, dt_write, dt_fixed_Option, &
                 iDisplayCycle_Option )

    REAL(DP), INTENT(in) :: t_begin, t_end, dt_write
    REAL(DP), INTENT(in), OPTIONAL :: dt_fixed_Option
    INTEGER,  INTENT(in), OPTIONAL :: iDisplayCycle_Option

    LOGICAL  :: WriteOutput = .FALSE.
    LOGICAL  :: FixedTimeStep
    INTEGER  :: iCycle, iDisplayCycle
    REAL(DP) :: t, t_write, dt, dt_fixed

    FixedTimeStep = .FALSE.
    IF( PRESENT( dt_fixed_Option ) )THEN
      FixedTimeStep = .TRUE.
      dt_fixed = dt_fixed_Option
    END IF

    iDisplayCycle = 10
    IF( PRESENT( iDisplayCycle_Option ) ) &
      iDisplayCycle = iDisplayCycle_Option

    ASSOCIATE( U => UnitsDisplay )

    iCycle  = 0
    t       = t_begin
    t_write = t_begin + dt_write

    CALL WriteFields1D &
           ( Time = t, WriteFluidFields_Option = .TRUE., &
             WriteRadiationFields_Option = .TRUE. )

    CALL InitializeFluidRadiationCoupling
  
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

      IF( MOD( iCycle, iDisplayCycle ) == 0 )THEN

        WRITE(*,'(A8,A8,I8.8,A2,A4,ES12.6E2,A1,A2,A2,A5,ES12.6E2,A1,A2)') &
          '', 'Cycle = ', iCycle, &
          '', 't = ', t / U % TimeUnit, '', TRIM( U % TimeLabel ), &
          '', 'dt = ', dt / U % TimeUnit, '', TRIM( U % TimeLabel )

      END IF
 
      CALL UpdateFields( dt )

      t = t + dt

      IF( WriteOutput )THEN

        CALL WriteFields1D &
               ( Time = t, WriteFluidFields_Option = .TRUE., &
                 WriteRadiationFields_Option = .TRUE. )

        WriteOutput = .FALSE.

      END IF

    END DO ! WHILE

    CALL FinalizeFluidRadiationCoupling

    WRITE(*,*)
    WRITE(*,'(A6,A15,ES10.4E2,A1,A2,A6,I8.8,A7,A4,ES10.4E2,A2)') &
      '', 'Evolved to t = ', t / U % TimeUnit, '', TRIM( U % TimeLabel ), &
      ' with ', iCycle, ' cycles'
    WRITE(*,*)

    END ASSOCIATE ! U
 
     CALL WriteFields1D &
           ( Time = t, WriteFluidFields_Option = .TRUE., &
             WriteRadiationFields_Option = .TRUE. )

  END SUBROUTINE EvolveFields


  SUBROUTINE ComputeTimestep( dt )

    REAL(DP), INTENT(out) :: dt
   
    REAL(DP)  :: LAMP
    REAL(DP)  :: dt_Fluid, dt_Radiation

    dt = 1.0d-2 * Millisecond
    dt_Fluid = 1.0d-2 * Millisecond    

    CALL ComputeTimestepPenalization &
         ( dt_Radiation, dt_max = 1.0d-2 * Millisecond, Diff = 1.0d-5 )

    dt = MIN( dt_Fluid, dt_Radiation )

  END SUBROUTINE ComputeTimestep


  SUBROUTINE ComputeTimestepPenalization( dt, dt_max, Diff )

    REAL(DP),                 INTENT(inout) :: dt
    REAL(DP),                 INTENT(in)  :: dt_max, Diff
   
    REAL(DP)                  :: LAMP, invLAMP, lim, dt_0, LLN, dt_1
    REAL(DP), DIMENSION(nNodesE) :: NN, Collision, &
                                lim1sign, lim1, lim1p, &
                                lim2sign, lim2, lim2p

    INTEGER :: iX

    dt = 1.0d-2 * Millisecond

!    DO iX = 1, nNodesX_G
!
!      NN = uPR_N(:,iPR_D,iX)
!
!      LAMP = RHSLAMP &
!             ( nNodesE_G, R0_In(:,:,iX), R0_Out(:,:,iX), Neq(:,iX) )
!
!      invLAMP = 1.0 / LAMP
!
!      ! --- Boundary Limit ---
!    
!      lim = MINVAL( lim1p + lim2p )    
!
!      IF( lim < invLAMP )THEN
!        dt_0 = lim / ( 1.0 - lim * LAMP )
!      ELSE
!        dt_0 = dt_max
!      END IF
!
!   
!      dt = MIN( dt_0, dt*1.05, dt_max )
!   
!      ! --- Accuracy Limit ---
!    
!      LLN = MAXVAL( ABS( Collision ) ) ! NEEDS CHANGE
!  
!      dt_1 = Diff / ( 2.0 * LLN ) + &
!             SQRT( Diff / ( LAMP * LLN ) + &
!                   Diff**2 / ( 4.0 * LLN**2 ) )
!                                         
!      dt = MIN( dt, dt_1 )
!
!      PRINT*, 'new dt: ', dt 
!
!   END DO

  END SUBROUTINE ComputeTimestepPenalization


  SUBROUTINE UpdateFields( dt )

    REAL(DP), INTENT(in) :: dt

    CALL ComputeRHS_C

  END SUBROUTINE UpdateFields


END MODULE TimeSteppingModule_Penalization
