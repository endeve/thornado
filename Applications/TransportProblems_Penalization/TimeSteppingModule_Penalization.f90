MODULE TimeSteppingModule_Penalization

  USE KindModule, ONLY: &
    DP, FourPi
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
    iCR_N, uCR, nCR, &
    iPR_D, uPR, &
    nSpecies
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

  REAL(DP), PARAMETER :: Diff_FE = 1.0d-5

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

      CALL ComputeRHS_C

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
    REAL(DP)  :: dt_Radiation, dt_max

    dt_max = 1.0d-2 * Millisecond

    CALL ComputeTimestepPenalization &
           ( dt_Radiation, dt_max = dt_max )

    dt = MIN( dt_Radiation, dt_max )

  END SUBROUTINE ComputeTimestep


  SUBROUTINE ComputeTimestepPenalization( dt, dt_max )

    REAL(DP), INTENT(out) :: dt
    REAL(DP), INTENT(in)  :: dt_max

    INTEGER  :: iE, iX1, iX2, iX3
    INTEGER  :: iNodeE, iNode
    INTEGER  :: iNodeX1, iNodeX2, iNodeX3, iNodeX
    REAL(DP) :: dt_1, dt_2, LAMB, LNMax, lim
    REAL(DP) :: NN, Coll, lim1, lim2

    dt = HUGE( 1.0_DP )

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)
          DO iE = 1, nE

            DO iNodeX3 = 1, nNodesX(3)
              DO iNodeX2 = 1, nNodesX(2)
                DO iNodeX1 = 1, nNodesX(1)

                  iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )
                  LAMB   = MAX( absLambda(iNodeX,iX1,iX2,iX3), 1.0d-100 )

                  DO iNodeE = 1, nNodesE

                    iNode = NodeNumber( iNodeE, iNodeX1, iNodeX2, iNodeX3 )

                    Coll = rhsCR_C(iNode,iE,iX1,iX2,iX3,iCR_N,1) 
                    NN   = uCR    (iNode,iE,iX1,iX2,iX3,iCR_N,1)

                    ! --- Boundary Limit ---

                    lim1 = - Coll / NN
                    lim2 =   Coll / ( FourPi - NN )

                    lim  = MAX( lim1, lim2 ) 

                    IF( lim > LAMB )THEN
                      dt_1 = 1.0 / ( lim - LAMB )
                    ELSE
                      dt_1 = dt_max
                    END IF
             
                    ! --- Accuracy Limit ---
              
                    LNMax = ABS( Coll / NN )

                    dt_2 = Diff_FE / ( 2.0 * LNMax ) + &
                           SQRT( Diff_FE / ( LAMB * LNMax ) + &
                                 Diff_FE**2 / ( 4.0 * LNMax**2 ) )

                    dt = MIN( dt, dt_1, dt_2 )

                  END DO

                END DO
              END DO
            END DO

          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE ComputeTimestepPenalization


  SUBROUTINE UpdateFields( dt )

    REAL(DP), INTENT(in) :: dt

    INTEGER  :: iE, iX1, iX2, iX3, iS, iCR
    INTEGER  :: iNodeE, iNodeX1, iNodeX2, iNodeX3
    INTEGER  :: iNode, iNodeX

    DO iS = 1, nSpecies

      DO iCR = 1, nCR

        DO iX3 = 1, nX(3)
          DO iX2 = 1, nX(2)
            DO iX1 = 1, nX(1)
              DO iE = 1, nE

                DO iNodeX3 = 1, nNodesX(3)
                  DO iNodeX2 = 1, nNodesX(2)
                    DO iNodeX1 = 1, nNodesX(1)

                      iNodeX &
                        = NodeNumberX &
                            ( iNodeX1, iNodeX2, iNodeX3 )

                      DO iNodeE = 1, nNodesE

                        iNode &
                          = NodeNumber &
                              ( iNodeE, iNodeX1, iNodeX2, iNodeX3 )

                        uCR(iNode,iE,iX1,iX2,iX3,iCR,iS) &
                          = uCR(iNode,iE,iX1,iX2,iX3,iCR,iS) &
                              + dt * rhsCR_C(iNode,iE,iX1,iX2,iX3,iCR,iS)

                      END DO
                    END DO
                  END DO
                END DO

              END DO
            END DO
          END DO
        END DO

      END DO

    END DO

  END SUBROUTINE UpdateFields


END MODULE TimeSteppingModule_Penalization
