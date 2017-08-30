MODULE TimeSteppingModule_Penalization

  USE KindModule, ONLY: &
    DP
  USE UnitsModule, ONLY: &
    UnitsDisplay, &
    Millisecond
  USE UtilitiesModule, ONLY: &
    NodeNumberX, &
    NodeNumber
  USE ProgramHeaderModule, ONLY: &
    nX, nNodesX, nDOFX, &
    nE, nNodesE, nDOFE, &
    nDOF
  USE RadiationFieldsModule, ONLY: &
    iCR_N, &
    iPR_D, uPR
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
   
    REAL(DP)                  :: dt_old, dt_i, dt_1, dt_2, &
                                 LAMB, LNMax, lim, lim1max, lim2max
    REAL(DP), DIMENSION(nNodesE,nE) :: NN, Coll, OnesMat, &
                                       lim1, lim2

    INTEGER :: iE, iX1, iX2, iX3
    INTEGER :: iNodeE, iNode
    INTEGER :: iNodeX1, iNodeX2, iNodeX3, iNodeX

    dt = 1.0d-2 * Millisecond
    
    dt_old = dt

    dt_i = dt_old

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)
          DO iNodeX3 = 1, nNodesX(3)
            DO iNodeX2 = 1, nNodesX(2)
              DO iNodeX1 = 1, nNodesX(1)

                iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

                LAMB   = MAX( absLambda(iNodeX,iX1,iX2,iX3), 1.0d-100 )

                DO iE = 1, nE
                  DO iNodeE = 1, nNodesE

                    iNode = NodeNumber( iNodeE, iNodeX1, iNodeX2, iNodeX3 )

                    Coll(iNodeE,iE) = rhsCR_C(iNode,iE,iX1,iX2,iX3,iCR_N,1) 
                     NN (iNodeE,iE) =  uPR   (iNode,iE,iX1,iX2,iX3,iPR_D,1)

                  END DO
                END DO

!               ! --- Boundary Limit ---
   
                OnesMat = 1.0

                lim1 = - Coll / NN
                lim2 =   Coll / ( OnesMat - NN )

                lim1max = MAXVAL( lim1 ) 
                lim2max = MAXVAL( lim2 ) 
                    lim = MAX( lim1max, lim2max ) 

                IF( lim > LAMB )THEN
                  dt_1 = 1.0 / ( lim - LAMB )
                ELSE
                  dt_1 = dt_max
                END IF
             
                ! --- Accuracy Limit ---
              
                LNMax = MAXVAL( ABS( Coll / NN ) )
            
                dt_2 = Diff / ( 2.0 * LNMax ) + &
                       SQRT( Diff / ( LAMB * LNMax ) + &
                             Diff**2 / ( 4.0 * LNMax**2 ) )
                                                   
                dt_i = MIN( dt_i, dt_max, dt_1, dt_2 )

              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

    dt = MIN( dt_i, 1.05*dt_old )

  END SUBROUTINE ComputeTimestepPenalization


  SUBROUTINE UpdateFields( dt )

    REAL(DP), INTENT(in) :: dt

    CALL ComputeRHS_C

  END SUBROUTINE UpdateFields


END MODULE TimeSteppingModule_Penalization
