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

  REAL(DP), PARAMETER :: Diff_FE = 1.0d-3

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
    INTEGER, PARAMETER :: out_unit = 20
    INTEGER, DIMENSION(4) :: SmallestPosition 
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
  
    OPEN(unit = out_unit, file = "tvsdt", action = "write", status = "replace" )

    dt = 1.d-8 * Millisecond 

    DO WHILE( t < t_end )

      iCycle = iCycle + 1

      CALL ComputeRHS_C_J

      IF( FixedTimeStep )THEN

        dt = dt_fixed

      ELSE

        CALL ComputeTimeStep( t, dt, SmallestPosition )

      END IF

      IF( t + dt > t_end )THEN

        dt = t_end - t

      END IF

      IF( t + dt > t_write )THEN

        t_write     = t_write + dt_write
        WriteOutput = .TRUE.

      END IF

      IF( MOD( iCycle, iDisplayCycle ) == 0 )THEN

        WRITE(*,'(A8,A8,I8.8,A2,A4,ES12.6E2,A1,A2,A2,A5,ES12.6E2,A1,A2,A1,4I4)') &
          '', 'Cycle = ', iCycle, &
          '', 't = ', t / U % TimeUnit, '', TRIM( U % TimeLabel ), &
          '', 'dt = ', dt / U % TimeUnit, '', TRIM( U % TimeLabel ), &
          '', SmallestPosition

      END IF

      CALL ComputeRHS_C_H( dt )
 
      CALL UpdateFields( dt, .true. )

      t = t + dt

      IF( WriteOutput )THEN

        WRITE( out_unit, '(2E15.6,4I4)' ) t/MilliSecond, dt/MilliSecond, SmallestPosition

        CALL WriteFields1D &
               ( Time = t, WriteFluidFields_Option = .TRUE., &
                 WriteRadiationFields_Option = .TRUE. )

        WriteOutput = .FALSE.

      END IF

    END DO ! WHILE

    CLOSE( out_unit )

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


  SUBROUTINE ComputeTimestep( t, dt, SmallestPosition )

    REAL(DP), INTENT(in) :: t
    INTEGER, DIMENSION(4), INTENT(out) :: SmallestPosition 
    REAL(DP), INTENT(inout) :: dt
   
    REAL(DP)  :: dt_Radiation, dt_max

    dt_max = 1.0d-2 * Millisecond

    CALL ComputeTimestepPenalization &
           ( dt_Radiation, dt_max, t, SmallestPosition )

    dt = MIN( dt_Radiation, dt_max, dt*1.1d0 )

  END SUBROUTINE ComputeTimestep


  SUBROUTINE ComputeTimestepPenalization( dt, dt_max, t, SmallestPosition )

    REAL(DP), INTENT(in)  :: dt_max, t
    REAL(DP), INTENT(out) :: dt
    INTEGER, DIMENSION(4), INTENT(out) :: SmallestPosition 

    INTEGER  :: out_unit_debug
    INTEGER  :: iE, iX1, iX2, iX3
    INTEGER  :: iNodeE, iNode
    INTEGER  :: iNodeX1, iNodeX2, iNodeX3, iNodeX
    CHARACTER(12)  :: FileNumberString
    CHARACTER(19)   :: FileName
    REAL(DP) :: dt_1, dt_2, LAMB, LNMax, lim
    REAL(DP) :: NN, Coll, lim1, lim2

    dt = HUGE( 1.0_DP )
    SmallestPosition = (/0,0,0,0/)

  WRITE( FileNumberString, FMT='(i12.12)') INT(t/MilliSecond*1.d7)
  FileName = 'o.dtdis'// FileNumberString
  OPEN( NEWUNIT = out_unit_debug, file = FileName, action = "write", status = "replace" )
    
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

                    IF( dt >= MIN( dt_1, dt_2 ) )THEN 
                      SmallestPosition = (/iE, iX1, iX2, iX3/)
                    END IF
 
                    dt = MIN( dt, dt_1, dt_2 )
                    WRITE(out_unit_debug,'(5I4,E15.6)'), iNode, iE, iX1, iX2, iX3, dt/ MilliSecond                   

                  END DO

                END DO
              END DO
            END DO

          END DO
        END DO
      END DO
    END DO

    CLOSE( out_unit_debug )

    IF( dt < 1.0d-15 * Millisecond ) THEN
      PRINT*,"dt too small: ", dt / MilliSecond, "ms"
      STOP
    END IF

  END SUBROUTINE ComputeTimestepPenalization


  SUBROUTINE ComputeTimestepForwardEuler( dt )

    REAL(DP), INTENT(out) :: dt

    INTEGER  :: iE, iX1, iX2, iX3
    INTEGER  :: iNodeE, iNode
    INTEGER  :: iNodeX1, iNodeX2, iNodeX3, iNodeX
    REAL(DP) :: Coll, NN

    dt = HUGE( 1.0_DP )

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)
          DO iE = 1, nE

            DO iNodeX3 = 1, nNodesX(3)
              DO iNodeX2 = 1, nNodesX(2)
                DO iNodeX1 = 1, nNodesX(1)

                  iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

                  DO iNodeE = 1, nNodesE

                    iNode = NodeNumber( iNodeE, iNodeX1, iNodeX2, iNodeX3 )

                    Coll = rhsCR_C(iNode,iE,iX1,iX2,iX3,iCR_N,1)

                    NN   = uCR    (iNode,iE,iX1,iX2,iX3,iCR_N,1)

                    IF( Coll > 0.d0 ) THEN
                      dt = MIN( dt, ( 1.d0 - NN ) / Coll )
                    ELSE IF( Coll < 0.d0 ) THEN
                      dt = MIN( dt, - NN / Coll )
                    END IF

                  END DO

                END DO
              END DO
            END DO

          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE ComputeTimestepForwardEuler


  SUBROUTINE UpdateFields( dt, PenalizationMethod_flag )

    REAL(DP), INTENT(in) :: dt
    LOGICAL,  INTENT(in) :: PenalizationMethod_flag

    INTEGER  :: iE, iX1, iX2, iX3, iS, iCR
    INTEGER  :: iNodeE, iNodeX1, iNodeX2, iNodeX3
    INTEGER  :: iNode, iNodeX
    REAL(DP) :: LAMB

    IF( PenalizationMethod_flag )THEN
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
                      LAMB   &
                        = MAX( absLambda(iNodeX,iX1,iX2,iX3), 1.0d-100 )

                      DO iNodeE = 1, nNodesE

                        iNode &
                          = NodeNumber &
                              ( iNodeE, iNodeX1, iNodeX2, iNodeX3 )

                        uCR(iNode,iE,iX1,iX2,iX3,iCR,iS) &
                          = uCR(iNode,iE,iX1,iX2,iX3,iCR,iS) &
                              + (dt/(1.d0 +dt*LAMB)) * rhsCR_C(iNode,iE,iX1,iX2,iX3,iCR,iS)

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

    ELSE

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
                      LAMB   &
                        = MAX( absLambda(iNodeX,iX1,iX2,iX3), 1.0d-100 )

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
    END IF

  END SUBROUTINE UpdateFields


END MODULE TimeSteppingModule_Penalization
