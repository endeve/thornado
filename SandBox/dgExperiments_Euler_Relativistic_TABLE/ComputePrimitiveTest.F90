PROGRAM ComputePrimitiveTest

  USE KindModule, ONLY: &
    DP, &
    Zero, &
    Half, &
    One, &
    Four, &
    Pi, &
    TwoPi, &
    SqrtTiny
  USE ProgramInitializationModule, ONLY: &
    InitializeProgram, &
    FinalizeProgram
  USE ProgramHeaderModule, ONLY: &
    iX_B0, &
    iX_B1, &
    iX_E0, &
    iX_E1, &
    nDOFX
  USE UnitsModule, ONLY: &
    Kilometer, &
    Gram, &
    Second, &
    Erg, &
    Centimeter, &
    SpeedOfLight
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    uGF, &
    nGF, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX
  USE FluidFieldsModule, ONLY: &
    uCF, &
    nCF, &
    iCF_D, &
    iCF_S1, &
    iCF_S2, &
    iCF_S3, &
    iCF_E, &
    iCF_Ne, &
    uPF, &
    nPF, &
    iPF_D, &
    iPF_V1, &
    iPF_V2, &
    iPF_V3, &
    iPF_E, &
    iPF_Ne, &
    uAF, &
    nAF, &
    iAF_P, &
    iAF_T, &
    iAF_E, &
    iAF_Ye
  USE EquationOfStateModule, ONLY: &
    InitializeEquationOfState, &
    FinalizeEquationOfState, &
    ComputeTemperatureFromPressure, &
    ComputeThermodynamicStates_Primitive
  USE Euler_UtilitiesModule_Relativistic, ONLY: &
    ComputeConserved_Euler_Relativistic, &
    ComputePrimitive_Euler_Relativistic, &
    ComputePrimitive_Euler_Relativistic_GPU
  USE Euler_ErrorModule, ONLY: &
    DescribeError_Euler

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  REAL(DP) :: xL(3), xU(3), X1, L
  INTEGER  :: nNodesX(3), nX(3)
  INTEGER  :: swX(3)

  INTEGER :: iNX, iX1, iX2, iX3, iNX1, iNX2, iNX3, iPF

  REAL(DP), ALLOCATABLE :: uPFo(:,:,:,:,:)
  REAL(DP), ALLOCATABLE :: uPFn(:,:,:,:,:)
  INTEGER,  ALLOCATABLE :: iErr(:,:,:,:)

  REAL(DP), PARAMETER :: UnitsD    = Gram / Centimeter**3
  REAL(DP), PARAMETER :: UnitsP    = Erg  / Centimeter**3
  REAL(DP), PARAMETER :: UnitsE    = Erg  / Centimeter**3
  REAL(DP), PARAMETER :: UnitsNe   = One  / Centimeter**3
  REAL(DP), PARAMETER :: UnitsV1   = Kilometer / Second
  REAL(DP), PARAMETER :: UnitsV2   = Kilometer / Second
  REAL(DP), PARAMETER :: UnitsV3   = Kilometer / Second
  REAL(DP), PARAMETER :: UnitsX1   = Kilometer
  REAL(DP), PARAMETER :: UnitsX2   = Kilometer
  REAL(DP), PARAMETER :: UnitsX3   = Kilometer
  REAL(DP), PARAMETER :: UnitsGm11 = One
  REAL(DP), PARAMETER :: UnitsGm22 = One
  REAL(DP), PARAMETER :: UnitsGm33 = One

  REAL(DP) :: D0, Amp, V1, V2, V3, P, Ye

  REAL(DP) :: Timer_old, Timer_new

  CHARACTER(64) :: EosTableName

  L = 1.0e2_DP * Kilometer

  nNodesX = [ 3, 1, 1 ]
  nX      = [ 128, 1, 1 ]
  swX     = [ 1,   0, 0 ]
  xL      = [ Zero, Zero, Zero ]
  xU      = [ One, One, One ] * L

  CALL InitializeProgram &
         ( ProgramName_Option &
             = 'ComputePrimitiveTest', &
           nX_Option &
             = nX, &
           swX_Option &
             = swX, &
           xL_Option &
             = xL, &
           xR_Option &
             = xU, &
           nNodes_Option &
             = nNodesX(1), &
           ActivateUnits_Option &
             = .TRUE., &
           BasicInitialization_Option &
             = .TRUE. )

  CALL ComputeGeometryX &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF )

  EosTableName = 'wl-EOS-SFHo-25-50-100.h5'

  CALL InitializeEquationOfState &
         ( EquationOfState_Option = 'TABLE', &
           EquationOfStateTableName_Option &
             = TRIM( EosTableName ), &
           Verbose_Option = .TRUE. )

  ! --- Define primitives ---

  D0  = 1.0e12_DP * UnitsD
  Amp = 0.1_DP * D0
  V1  = 3.0e4_DP * UnitsV1
  V2  = 0.0_DP * UnitsV2
  V3  = 0.0_DP * UnitsV3
  P   = 0.01_DP * D0 * SpeedOfLight**2
  Ye  = 0.3_DP

  ASSOCIATE( X1C => MeshX(1) % Center, dX1 => MeshX(1) % Width )

  ALLOCATE( iErr(nDOFX,iX_B0(1):iX_E0(1), &
                       iX_B0(2):iX_E0(2), &
                       iX_B0(3):iX_E0(3)) )
  ALLOCATE( uPFo(nDOFX,iX_B1(1):iX_E1(1), &
                       iX_B1(2):iX_E1(2), &
                       iX_B1(3):iX_E1(3),nPF) )

#if defined(THORNADO_OMP_OL)
  !$OMP TARGET ENTER DATA &
  !$OMP MAP( to:    X1C, dX1 ) &
  !$OMP MAP( alloc: iErr, uPFo )
#elif defined(THORNADO_OACC)
  !$ACC ENTER DATA &
  !$ACC COPYIN(     X1C, dX1 ) &
  !$ACC CREATE(     iErr, uPFo )
#endif

#if defined(THORNADO_OMP_OL)
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
  !$OMP PRIVATE( iNX, X1 )
#elif defined(THORNADO_OACC)
  !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
  !$ACC PRIVATE( iNX, X1 ) &
  !$ACC PRESENT( iX_B0, iX_E0, uGF, uCF, uPF, uAF, X1C, dX1 )
#elif defined(THORNADO_OMP)
  !$OMP PARALLEL DO COLLAPSE(3) &
  !$OMP PRIVATE( iNX, X1 )
#endif
  DO iX3  = iX_B0(3), iX_E0(3)
  DO iX2  = iX_B0(2), iX_E0(2)
  DO iX1  = iX_B0(1), iX_E0(1)

    iNX = 0

    DO iNX3 = 1, nNodesX(3)
    DO iNX2 = 1, nNodesX(2)
    DO iNX1 = 1, nNodesX(1)

      X1 = NodeCoordinate( X1C(iX1), dX1(iX1), DBLE( iNX1 ) )

      iNX = iNX + 1

      uPF(iNX,iX1,iX2,iX3,iPF_D ) = ( D0 + Amp * SIN( TwoPi * X1 / L ) )
      uPF(iNX,iX1,iX2,iX3,iPF_V1) = V1
      uPF(iNX,iX1,iX2,iX3,iPF_V2) = V2
      uPF(iNX,iX1,iX2,iX3,iPF_V3) = V3

      uAF(iNX,iX1,iX2,iX3,iAF_P ) = P
      uAF(iNX,iX1,iX2,iX3,iAF_Ye) = Ye

      CALL ComputeTemperatureFromPressure &
             ( uPF(iNX,iX1,iX2,iX3,iPF_D ), &
               uAF(iNX,iX1,iX2,iX3,iAF_P ), &
               uAF(iNX,iX1,iX2,iX3,iAF_Ye), &
               uAF(iNX,iX1,iX2,iX3,iAF_T ) )

      CALL ComputeThermodynamicStates_Primitive &
             ( uPF(iNX,iX1,iX2,iX3,iPF_D ), &
               uAF(iNX,iX1,iX2,iX3,iAF_T ), &
               uAF(iNX,iX1,iX2,iX3,iAF_Ye), &
               uPF(iNX,iX1,iX2,iX3,iPF_E ), &
               uAF(iNX,iX1,iX2,iX3,iAF_E ), &
               uPF(iNX,iX1,iX2,iX3,iPF_Ne) )

      CALL ComputeConserved_Euler_Relativistic &
             ( uPF(iNX,iX1,iX2,iX3,iPF_D ), &
               uPF(iNX,iX1,iX2,iX3,iPF_V1), &
               uPF(iNX,iX1,iX2,iX3,iPF_V2), &
               uPF(iNX,iX1,iX2,iX3,iPF_V3), &
               uPF(iNX,iX1,iX2,iX3,iPF_E ), &
               uPF(iNX,iX1,iX2,iX3,iPF_Ne), &
               uCF(iNX,iX1,iX2,iX3,iCF_D ), &
               uCF(iNX,iX1,iX2,iX3,iCF_S1), &
               uCF(iNX,iX1,iX2,iX3,iCF_S2), &
               uCF(iNX,iX1,iX2,iX3,iCF_S3), &
               uCF(iNX,iX1,iX2,iX3,iCF_E ), &
               uCF(iNX,iX1,iX2,iX3,iCF_Ne), &
               uGF(iNX,iX1,iX2,iX3,iGF_Gm_dd_11), &
               uGF(iNX,iX1,iX2,iX3,iGF_Gm_dd_22), &
               uGF(iNX,iX1,iX2,iX3,iGF_Gm_dd_33), &
               uAF(iNX,iX1,iX2,iX3,iAF_P) )

    END DO
    END DO
    END DO

  END DO
  END DO
  END DO

  ! --- Compute primitives with old version of subroutine ---

  Timer_old = MPI_WTIME()

#if defined(THORNADO_OMP_OL)
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
  !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
  !$ACC PRESENT( iX_B0, iX_E0, uGF, uCF, uPFo, iErr )
#elif defined(THORNADO_OMP)
  !$OMP PARALLEL DO COLLAPSE(4)
#endif
  DO iX3 = iX_B0(3), iX_E0(3)
  DO iX2 = iX_B0(2), iX_E0(2)
  DO iX1 = iX_B0(1), iX_E0(1)
  DO iNX = 1, nDOFX

    iErr(iNX,iX1,iX2,iX3) = 0

    CALL ComputePrimitive_Euler_Relativistic &
           ( uCF (iNX,iX1,iX2,iX3,iCF_D ), &
             uCF (iNX,iX1,iX2,iX3,iCF_S1), &
             uCF (iNX,iX1,iX2,iX3,iCF_S2), &
             uCF (iNX,iX1,iX2,iX3,iCF_S3), &
             uCF (iNX,iX1,iX2,iX3,iCF_E ), &
             uCF (iNX,iX1,iX2,iX3,iCF_Ne), &
             uPFo(iNX,iX1,iX2,iX3,iPF_D ), &
             uPFo(iNX,iX1,iX2,iX3,iPF_V1), &
             uPFo(iNX,iX1,iX2,iX3,iPF_V2), &
             uPFo(iNX,iX1,iX2,iX3,iPF_V3), &
             uPFo(iNX,iX1,iX2,iX3,iPF_E ), &
             uPFo(iNX,iX1,iX2,iX3,iPF_Ne), &
             uGF (iNX,iX1,iX2,iX3,iGF_Gm_dd_11), &
             uGF (iNX,iX1,iX2,iX3,iGF_Gm_dd_22), &
             uGF (iNX,iX1,iX2,iX3,iGF_Gm_dd_33), &
             iErr(iNX,iX1,iX2,iX3) )

  END DO
  END DO
  END DO
  END DO

#if defined(THORNADO_OMP_OL)
  !$OMP TARGET EXIT DATA &
  !$OMP MAP( from:    iErr, uPFo ) &
  !$OMP MAP( release: X1C, dX1 )
#elif defined(THORNADO_OACC)
  !$ACC EXIT DATA &
  !$ACC COPYOUT(      iErr, uPFo ) &
  !$ACC DELETE(       X1C, dX1 )
#endif

  END ASSOCIATE ! X1C, dX1

  DO iX3 = iX_B0(3), iX_E0(3)
  DO iX2 = iX_B0(2), iX_E0(2)
  DO iX1 = iX_B0(1), iX_E0(1)
  DO iNX = 1, nDOFX

    CALL DescribeError_Euler( iErr(iNX,iX1,iX2,iX3) )

  END DO
  END DO
  END DO
  END DO

  Timer_old = MPI_WTIME() - Timer_old

  ! --- Compute primitives with new version of subroutine ---

  ALLOCATE( uPFn(nDOFX,iX_B1(1):iX_E1(1), &
                       iX_B1(2):iX_E1(2), &
                       iX_B1(3):iX_E1(3),nPF) )

  Timer_new = MPI_WTIME()

  CALL ComputePrimitive_Euler_Relativistic_GPU &
         ( nDOFX, iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uPFn )

  Timer_new = MPI_WTIME() - Timer_new

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE FROM( uPF )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE HOST       ( uPF )
#endif

  WRITE(*,'(A,ES11.3E3,A)') 'Timer_old: ', Timer_old, ' s'
  WRITE(*,'(A,ES11.3E3,A)') 'Timer_new: ', Timer_new, ' s'

  WRITE(*,*) 'Maximum relative difference'
  DO iPF = 1, nPF

    WRITE(*,'(A,I2.2,ES11.3E3,1x,ES11.3E3)') &
      'iPF = ', iPF, &
      MAXVAL( ABS(   uPF (:,iX_B0(1):iX_E0(1), &
                            iX_B0(2):iX_E0(2), &
                            iX_B0(3):iX_E0(3),iPF) &
                   - uPFo(:,iX_B0(1):iX_E0(1), &
                            iX_B0(2):iX_E0(2), &
                            iX_B0(3):iX_E0(3),iPF) ) &
               / ABS( uPF(:,iX_B0(1):iX_E0(1), &
                            iX_B0(2):iX_E0(2), &
                            iX_B0(3):iX_E0(3),iPF)+SqrtTiny ) ), &
      MAXVAL( ABS(   uPF (:,iX_B0(1):iX_E0(1), &
                            iX_B0(2):iX_E0(2), &
                            iX_B0(3):iX_E0(3),iPF) &
                   - uPFn(:,iX_B0(1):iX_E0(1), &
                            iX_B0(2):iX_E0(2), &
                            iX_B0(3):iX_E0(3),iPF) ) &
               / ABS( uPF(:,iX_B0(1):iX_E0(1), &
                            iX_B0(2):iX_E0(2), &
                            iX_B0(3):iX_E0(3),iPF)+SqrtTiny ) )

  END DO

  DEALLOCATE( uPFn )
  DEALLOCATE( uPFo )
  DEALLOCATE( iErr )

  CALL FinalizeEquationOfState

  CALL FinalizeProgram

END PROGRAM ComputePrimitiveTest
