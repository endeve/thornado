PROGRAM ScalarWaveTest

  USE KindModule, ONLY: &
    DP, One, Two, Zero, TwoPi
  USE ProgramHeaderModule, ONLY: &
    iX_B0, iX_B1, iX_E0, iX_E1, &
    nDimsX
  USE ReferenceElementModuleX, ONLY: &
    InitializeReferenceElementX, &
    FinalizeReferenceElementx
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    InitializeReferenceElementX_Lagrange, &
    FinalizeReferenceElementX_Lagrange
  USE ScalarWave_ProgramInitializationModule, ONLY: &
    InitializeProgram, &
    FinalizeProgram
  USE InitializationModule, ONLY: &
    InitializeFields
  USE ScalarWave_InputOutputModuleHDF, ONLY: &
    WriteFieldsHDF, &
    ReadFieldsHDF  
  USE TimeSteppingModule_SSPRK, ONLY: &
    InitializeScalarWave_SSPRK, &
    FinalizeScalarWave_SSPRK, &
    UpdateScalarWave_SSPRK
  USE ScalarFieldsModule, ONLY : &
    uSF, &
    ComputeTimeStep_ScalarWave
  USE ScalarWave_dgDiscretizationModule, Only : &
    ComputeIncrement_ScalarWave_DG_Explicit

  IMPLICIT NONE

  LOGICAL  :: wrt
  INTEGER  :: iCycle, nNodes
  REAL(DP) :: t, dt, t_end, dt_wrt, t_wrt

  t = 0.0_DP
  t_end = 1.0_DP
  dt_wrt = 1.0d-2 * t_end

  nNodes = 2

  CALL InitializeProgram &
         ( ProgramName_Option = 'SpacetimeTest', &
           nX_Option = [ 16, 1, 1 ], swX_Option = [ 1, 0, 0 ], &
           bcX_Option = [ 1, 0, 0 ], xL_Option = [ Zero, Zero, Zero ], &
           xR_Option = [ TwoPi, One, One ], nNodes_Option = 1 )

  CALL InitializeReferenceElementX

  CALL InitializeReferenceElementX_Lagrange

  Call InitializeFields

  CALL WriteFieldsHDF( 0.0_DP, WriteSF_Option=.true.)

  CALL InitializeScalarWave_SSPRK( nStages = 2 )

  t_wrt = t + dt_wrt
  wrt   = .FALSE.

  iCycle = 0
  DO WHILE ( t < t_end )
    
    iCycle = iCycle + 1
print*,iCycle
    CALL ComputeTimeStep_ScalarWave &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uSF, &
             CFL = 0.5_DP / (nDimsX * ( Two * DBLE ( nNodes ) - One ) ), &
             Timestep = dt )

    IF( t + dt > t_end )THEN

      dt = t_end - t

    END IF

    IF( t + dt > t_wrt )THEN

      dt    = t_wrt - t
      t_wrt = t_wrt + dt_wrt
      wrt   = .TRUE.

    END IF

    CALL UpdateScalarWave_SSPRK &
           ( t, dt, uSF, ComputeIncrement_ScalarWave_DG_Explicit )

    t = t + dt

    IF( wrt )THEN
!print*,"Start write"
!      CALL WriteFieldsHDF &
!             ( t, WriteSF_Option = .TRUE. )
!print*,"End write"
      wrt = .FALSE.

    END IF

  END DO

  CALL WriteFieldsHDF &
         ( t, WriteSF_Option = .TRUE. )


  WRITE(*,*)
  WRITE(*,'(A6,A,I6.6,A)') &
    '', 'Finished ', iCycle, ' Cycles.'
  WRITE(*,*)

print*,"Start Finalize SW_SSPRK"
  CALL FinalizeScalarWave_SSPRK
print*,"Finish Finalize SW_SSPRK"
print*,"Start Finalize REX_L"
  CALL FinalizeReferenceElementX_Lagrange
print*,"Finish Finalize REX_L"  
print*,"Start Finalize REX"
  CALL FinalizeReferenceElementX
print*,"Finish Finalize REX"
print*,"Start Finalize Program"
  CALL FinalizeProgram
print*,"Finish Finalize Program"

END PROGRAM ScalarWaveTest
