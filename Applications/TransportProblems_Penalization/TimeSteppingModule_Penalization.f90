MODULE TimeSteppingModule_Penalization

  USE KindModule, ONLY: &
    DP, FourPi
  USE UnitsModule, ONLY: &
    UnitsDisplay, &
    Kilometer, &
    Second, &
    Millisecond
  USE UtilitiesModule, ONLY: &
    NodeNumberX, &
    NodeNumber
  USE ProgramHeaderModule, ONLY: &
    nX, nDOFX, nNodes, &
    nE, nDOF
  USE MeshModule, ONLY: &
    MeshX
  USE GeometryFieldsModule, ONLY: &
    a, b, c
  USE ProgramHeaderModule, ONLY: &
    nX, nNodesX, nDOFX, &
    nE, nNodesE, nDOFE, &
    iX_B0, iX_B1, &
    iX_E0, iX_E1, &
    nDOF
  USE RadiationFieldsModule, ONLY: &
    iCR_N, iCR_G1, iCR_G2, iCR_G3, uCR, nCR, &
    iPR_D, uPR, &
    nSpecies, rhsCR
  USE MomentEquationsUtilitiesModule, ONLY: &
    ComputePrimitiveMoments
  USE RadiationEvolutionModule, ONLY: &
    ApplyPositivityLimiter_Radiation
  USE InputOutputModule, ONLY: &
    WriteFields1D
  USE MomentEquationsSolutionModule_M1_DG, ONLY: &
    ComputeRHS_M1_DG
  USE FluidRadiationCouplingUtilitiesModule, ONLY: &
    FinalizeFluidFields, &
    FinalizeRadiationFields
  USE FluidRadiationCouplingSolutionModule_Penalization, ONLY: &
    absLambda, C_J, Kappa, &
    InitializeFluidRadiationCoupling, &
    FinalizeFluidRadiationCoupling, &
    ComputeRHS_C_J, &
    ComputeRHS_C_H
  USE BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_Radiation

  IMPLICIT NONE
  PRIVATE

  REAL(DP), PARAMETER :: Diff_FE = 1.0d-2
  REAL(DP), PARAMETER :: Min_dt  = 1.0d-10 * MilliSecond
  REAL(DP), PARAMETER :: Min_N   = 0.0_DP 
  REAL(DP), PARAMETER :: dt_incf = 1.1_DP
  REAL(DP), PARAMETER :: dt_shrf = 1.0_DP
  REAL(DP), PARAMETER :: dt_init = 1.0d-5  * MilliSecond

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
    REAL(DP) :: t, t_write, dt, dt_fixed, &
                dt_stream, dt_accur, dt_boundary, dt_lower, dt_upper

    FixedTimeStep = .FALSE.
    IF( PRESENT( dt_fixed_Option ) )THEN
      FixedTimeStep = .TRUE.
      dt_fixed = dt_fixed_Option
    END IF

    iDisplayCycle = 10
    IF( PRESENT( iDisplayCycle_Option ) ) &
      iDisplayCycle = iDisplayCycle_Option

    ASSOCIATE( U => UnitsDisplay )

    CALL ApplyPositivityLimiter_Radiation &
	   ( iX_B0, iX_E0, iX_B1, iX_E1, uCR )

    iCycle  = 0
    t       = t_begin
    t_write = t_begin + dt_write

    CALL WriteFields1D &
           ( Time = t, WriteFluidFields_Option = .TRUE., &
             WriteRadiationFields_Option = .TRUE. )

    CALL InitializeFluidRadiationCoupling

    WRITE(*,'(A10,ES10.3)') 'Diff_FE =', Diff_FE 
    WRITE(*,'(A10,ES10.3)') 'Min_dt  =', Min_dt
    WRITE(*,'(A10,ES10.3)') 'Min_N   =', Min_N   
    WRITE(*,'(A10,ES10.3)') 'dt_incf =', dt_incf 
    WRITE(*,'(A10,ES10.3)') 'dt_shrf =', dt_shrf 
    WRITE(*,*) '' 
 
    OPEN( unit = out_unit, file = "tvsdt", action = "write", &
          status = "replace" )
    WRITE( out_unit, '(A1,8A15)') '#', 'time', 'dt', 'dt_stream', &
      'dt_accur','dt_boundary', 'dt_lower','dt_upper','Position'

    dt =  dt_init

    DO WHILE( t < t_end )

      iCycle = iCycle + 1

      CALL ComputeRHS_C_J  ! compute absLambda and C_J

      CALL ApplyBoundaryConditions_Radiation( t )

      CALL ComputeRHS_M1_DG  &
             ( iX_B0, iX_E0, iX_B1, iX_E1, uCR, rhsCR )

      CALL ComputeRHS_C_H ! compute Kappa

      IF( FixedTimeStep )THEN

        CALL ComputeTimeStep &
               ( t, dt, SmallestPosition, dt_stream, dt_accur, &
                 dt_boundary, dt_lower, dt_upper )

        dt = dt_fixed

      ELSE

        CALL ComputeTimeStep &
               ( t, dt, SmallestPosition, dt_stream, dt_accur, &
                 dt_boundary, dt_lower, dt_upper )

      END IF

      IF( t + dt > t_end )THEN

        dt = t_end - t

      END IF

      IF( t + dt > t_write )THEN

        t_write     = t_write + dt_write
        WriteOutput = .TRUE.

      END IF

      IF( MOD( iCycle, iDisplayCycle ) == 0 )THEN

        WRITE(*,'(A8,A8,I8.8,A2,A4,ES12.6E2,A1,A2,A2,A5,ES13.6E3,A1,A2,A1,4I4)') &
          '', 'Cycle = ', iCycle, &
          '', 't = ', t / U % TimeUnit, '', TRIM( U % TimeLabel ), &
          '', 'dt = ', dt / U % TimeUnit, '', TRIM( U % TimeLabel ), &
          '', SmallestPosition

        WRITE( out_unit, '(7E15.6,4I4)' ) &
          t / MilliSecond, dt / MilliSecond, &
          dt_stream / MilliSecond, dt_accur / (Second*Kilometer), &
          dt_boundary / MilliSecond, dt_lower / MilliSecond, &
          dt_upper / MilliSecond, SmallestPosition

      END IF

      CALL UpdateFields( dt, iCycle )

      CALL ApplyPositivityLimiter_Radiation &
	     ( iX_B0, iX_E0, iX_B1, iX_E1, uCR )

      CALL ComputePrimitiveMoments &
           ( [ 1, 1, 1 ], [ nX(1), nX(2), nX(3) ] )

      t = t + dt

      IF( dt < Min_dt ) THEN
         WriteOutput = .TRUE.
         WRITE(*,*)
         WRITE(*,'(A1,A42,ES10.4E2,A1,A2,A6,ES13.4E3,A1,A2,A6,I8,A10)' ) &
           '','ERROR: dt too small! Impose to end at t = ', t / U % TimeUnit, &
           '', TRIM( U % TimeLabel ),' dt = ',dt / U % TimeUnit,&
           '', TRIM( U % TimeLabel ),' with ', iCycle, ' cycles.'
         WRITE(*,*) 'Write out and end the loop. '
         WRITE(*,*)

         WRITE( out_unit, '(6E15.6,4I4)' ) &
           t / MilliSecond, dt / MilliSecond, &
           dt_stream / MilliSecond, dt_accur / (Second*Kilometer), &
           dt_boundary / MilliSecond, dt_lower / MilliSecond, &
           dt_upper / MilliSecond, SmallestPosition
         CLOSE( out_unit )

        CALL WriteFields1D &
               ( Time = t, WriteFluidFields_Option = .TRUE., &
                 WriteRadiationFields_Option = .TRUE. )

         RETURN
      END IF

      IF( WriteOutput )THEN

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


  SUBROUTINE ComputeTimestep &
               ( t, dt, SmallestPosition, dt_streams, dt_accur, &
                 dt_boundary, dt_lower, dt_upper )

    REAL(DP), INTENT(in) :: t
    INTEGER, DIMENSION(4), INTENT(out) :: SmallestPosition 
    REAL(DP), INTENT(inout) :: dt
    REAL(DP), INTENT(out) ::  dt_streams, dt_accur, dt_boundary,&
                              dt_lower, dt_upper
   
    REAL(DP)  :: dt_Radiation, dt_max, dt_Stream, dt_buffer

    dt_max = 1.0d-2 * Millisecond

    CALL ComputeTimestepPenalization &
           ( dt_Radiation, dt_max, t, SmallestPosition, dt_accur, &
             dt_boundary, dt_lower, dt_upper )

    CALL ComputeTimeStep_Streaming( dt_Stream )
    dt_streams = dt_Stream

    dt_buffer =  dt * dt_incf
    dt_Radiation = dt_Radiation * dt_shrf
    dt = MIN( dt_Radiation, dt_max, dt_Stream, dt_buffer )

    IF( dt < Min_dt ) THEN
      PRINT*,''
      PRINT*,'In ComputeTimestep, dt too small: ', dt / MilliSecond, 'ms'
      PRINT*,'dt_Stream ', dt_Stream / Millisecond, 'ms'
      PRINT*,'dt_Radiation ', dt_Radiation / Millisecond, 'ms'
      PRINT*,'    dt_accur ', dt_accur     / Millisecond, 'ms'
      PRINT*,'    dt_bound ', dt_boundary  / Millisecond, 'ms'
      PRINT*,'    dt_lower ', dt_lower     / Millisecond, 'ms'
      PRINT*,'    dt_upper ', dt_upper     / Millisecond, 'ms'
      PRINT*,'Smallest pos ', SmallestPosition
      RETURN
    END IF

  END SUBROUTINE ComputeTimestep


  SUBROUTINE ComputeTimestepPenalization &
  ( dt, dt_max, t, SmallestPosition, Cri_kappadx, dt_boundary, dt_lower, dt_upper )
!  ( dt, dt_max, t, SmallestPosition, dt_accur, dt_boundary, dt_lower, dt_upper )

    REAL(DP), INTENT(in)  :: dt_max, t
    REAL(DP), INTENT(out) :: dt, Cri_kappadx, dt_lower, dt_upper, dt_boundary
  !  REAL(DP), INTENT(out) :: dt, dt_accur, dt_lower, dt_upper, dt_boundary
    INTEGER, DIMENSION(4), INTENT(out) :: SmallestPosition 

    INTEGER  :: out_unit_debug
    INTEGER  :: iE, iX1, iX2, iX3
    INTEGER  :: iNodeE, iNode
    INTEGER  :: iNodeX1, iNodeX2, iNodeX3, iNodeX
    CHARACTER(12)  :: FileNumberString
    CHARACTER(19)   :: FileName
    REAL(DP) :: LAMB, LNMax, lim, dt_dx, lim_W, dt_buffer
    REAL(DP) :: NN, Coll, CpS
    REAL(DP) :: Invmindx, maxdx, dt_accur
    !REAL(DP) :: Invmindx

    dt = HUGE( 1.0_DP )
    dt_accur = dt
    dt_boundary = dt
    dt_lower = dt
    dt_upper = dt
   
    Cri_kappadx = 0.0

    SmallestPosition = (/0,0,0,0/)
    ASSOCIATE( dX1 => MeshX(1) % Width (1:nX(1)) )

    Invmindx = 1.d0 / MINVAL( dX1 )
    maxdx    = MAXVAL( dX1 )
 
    END ASSOCIATE

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)
          DO iE = 1, nE

            DO iNodeX3 = 1, nNodesX(3)
              DO iNodeX2 = 1, nNodesX(2)
                DO iNodeX1 = 1, nNodesX(1)

                  iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )
                  LAMB   = absLambda(iNodeX,iX1,iX2,iX3)

                  Cri_kappadx = MAX( LAMB*maxdx, Cri_kappadx )

                  DO iNodeE = 1, nNodesE

                    iNode = NodeNumber( iNodeE, iNodeX1, iNodeX2, iNodeX3 )

                    Coll = C_J(iNode,iE,iX1,iX2,iX3,1)
                    CpS  = Coll + rhsCR(iNode,iE,iX1,iX2,iX3,iCR_N,1)
                    NN   = uCR(iNode,iE,iX1,iX2,iX3,iCR_N,1)

                    ! --- Boundary Limit (lower) ---

                    lim_W = Invmindx - LAMB - Coll / NN

                    dt_buffer = HUGE( 1.0_DP )
                    IF( lim_W > 0.0d0 ) dt_buffer = 1.0d0 / lim_W

                    dt_lower = MIN( dt_lower, dt_buffer )

                    ! --- Boundary Limit (upper) ---
                    lim_W = Invmindx - LAMB + Coll / ( FourPi - NN )

                    dt_buffer = HUGE( 1.0_DP )
                    IF( lim_W > 0.0d0 ) dt_buffer = 1.0d0 / lim_W

                    dt_upper = MIN( dt_upper, dt_buffer )

                    ! --- Boundary Limit ---

                    dt_boundary = MIN( dt_boundary, dt_upper, dt_lower )

                    ! --- Accuracy Limit ---
              
!                    LNMax = ABS( CpS ) / MAX( NN, 1.0d-16 )
!
!                    dt_buffer = Diff_FE / ( 2.0d0 * LNMax ) + &
!                           SQRT( Diff_FE / ( LAMB * LNMax ) + &
!                                 Diff_FE * Diff_FE &
!                                 / ( 4.0d0 * LNMax * LNMax ) )
!
!                    dt_accur = MIN( dt_accur, dt_buffer )

                    dt_buffer = MIN( dt_boundary, dt_accur )

                    IF( dt > dt_buffer )THEN
                      dt = dt_buffer 
                      SmallestPosition = (/iE, iX1, iX2, iX3/)
                    END IF
 
                  END DO

                END DO
              END DO
            END DO

          END DO
        END DO
      END DO
    END DO

   
  END SUBROUTINE ComputeTimestepPenalization


  SUBROUTINE ComputeTimeStep_Streaming( dt )

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

    CFL = 0.4_DP / ( 2.0_DP * DBLE( nNodes - 1 ) + 1.0_DP ) ! For Debugging

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          dt_X1 = CFL * dX1(iX1)
          dt_X2 = CFL * a( [ X1(iX1), X2(iX2), X3(iX3) ] ) &
                      * dX2(iX2)
          dt_X3 = CFL * b( [ X1(iX1), X2(iX2), X3(iX3) ] ) &
                      * c( [ X1(iX1), X2(iX2), X3(iX3) ] ) &
                      * dX3(iX3)

          dt = MIN( dt, dt_x1, dt_X2, dt_X3 )

        END DO
      END DO
    END DO

    END ASSOCIATE ! dX1, etc.

  END SUBROUTINE ComputeTimeStep_Streaming


  SUBROUTINE UpdateFields( dt, iCycle )

    INTEGER,  INTENT(in)    :: iCycle
    REAL(DP), INTENT(inout) :: dt

    INTEGER  :: iE, iX1, iX2, iX3, iS, iCR
    INTEGER  :: iNodeE, iNodeX1, iNodeX2, iNodeX3
    INTEGER  :: iNode, iNodeX
    REAL(DP) :: LAMB, temp

    DO iS = 1, nSpecies
  
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

                    LAMB = MAX( absLambda(iNodeX,iX1,iX2,iX3), 1.d-100 )

                    DO iNodeE = 1, nNodesE

                      iNode &
                        = NodeNumber( iNodeE, iNodeX1, iNodeX2, iNodeX3 )

                      temp &
                        = uCR(iNode,iE,iX1,iX2,iX3,iCR_N,iS) &
                            + ( dt / (1.d0+dt*LAMB) ) &
                              * ( rhsCR(iNode,iE,iX1,iX2,iX3,iCR_N,iS) &
                                  + C_J(iNode,iE,iX1,iX2,iX3,iS) )

                      uCR(iNode,iE,iX1,iX2,iX3,iCR_N,iS) &
                        = MAX( temp, Min_N )

                      temp = 1.d0 / ( 1.d0 + dt * LAMB )

                      uCR(iNode,iE,iX1,iX2,iX3,iCR_G1,iS) &
                        = temp &
                          * ( uCR(iNode,iE,iX1,iX2,iX3,iCR_G1,iS) &
                               * ( 1.d0 + dt * LAMB - &
                                   dt * Kappa(iNode,iE,iX1,iX2,iX3) ) &
                              + dt * rhsCR(iNode,iE,iX1,iX2,iX3,iCR_G1,iS) )

                      uCR(iNode,iE,iX1,iX2,iX3,iCR_G2,iS) &
                        = temp &
                          * ( uCR(iNode,iE,iX1,iX2,iX3,iCR_G2,iS) &
                               * ( 1.d0 + dt * LAMB - &
                                   dt * Kappa(iNode,iE,iX1,iX2,iX3) ) &
                              + dt * rhsCR(iNode,iE,iX1,iX2,iX3,iCR_G2,iS) )

                      uCR(iNode,iE,iX1,iX2,iX3,iCR_G3,iS) &
                        = temp &
                          * ( uCR(iNode,iE,iX1,iX2,iX3,iCR_G3,iS) &
                               * ( 1.d0 + dt * LAMB - &
                                   dt * Kappa(iNode,iE,iX1,iX2,iX3) ) &
                              + dt * rhsCR(iNode,iE,iX1,iX2,iX3,iCR_G3,iS) )

                    END DO

                  END DO
                END DO
              END DO

            END DO
          END DO
        END DO
      END DO

    END DO

    IF( MINVAL( uCR(:,:,:,:,:,iCR_N,:) ) < Min_N )THEN
      PRINT*,''
      PRINT*,'ERROR in UpdateFields'
      PRINT*,'NEGATIVE or ZERO uCR(:,iCR_N):',MINVAL( uCR(:,:,:,:,:,iCR_N,:) )
      PRINT*,'dt =', dt
      PRINT*,'iCycle = ', iCycle
      PRINT*,''
      STOP
    END IF

    IF( MAXVAL( uCR(:,:,:,:,:,iCR_N,:) ) > FourPi )THEN
      PRINT*,''
      PRINT*,'ERROR in UpdateFields'
      PRINT*,'uCR(:,iCR_N) beyond 4Pi',MAXVAL( uCR(:,:,:,:,:,iCR_N,:) )
      PRINT*,'dt =', dt / Millisecond, 'ms'
      PRINT*,'iCycle = ', iCycle
      PRINT*,''
    !  STOP
    END IF


  END SUBROUTINE UpdateFields

END MODULE TimeSteppingModule_Penalization
