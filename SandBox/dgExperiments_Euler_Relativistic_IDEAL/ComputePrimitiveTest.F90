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
    ComputePrimitive_Euler_Relativistic
  USE Euler_ErrorModule, ONLY: &
    DescribeError_Euler
  USE TimersModule_Euler, ONLY: &
    TimeIt_Euler, &
    InitializeTimers_Euler, &
    FinalizeTimers_Euler

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  REAL(DP) :: xL(3), xU(3), X1, L
  INTEGER  :: nNodesX(3), nX(3)
  INTEGER  :: swX(3)

  INTEGER :: iNX, iX1, iX2, iX3, iNX1, iNX2, iNX3

  REAL(DP), POINTER, CONTIGUOUS :: uD(:), uS1(:), uS2(:), uS3(:), uE(:), uNe(:)
  REAL(DP), POINTER, CONTIGUOUS :: Gm11(:), Gm22(:), Gm33(:)
  REAL(DP), POINTER, CONTIGUOUS :: p1D(:), p1V1(:), p1V2(:), &
                                   p1V3(:), p1E(:), p1Ne(:)
  REAL(DP), ALLOCATABLE :: p2D(:), p2V1(:), p2V2(:), p2V3(:), p2E(:), p2Ne(:)
  REAL(DP), ALLOCATABLE :: p3D(:), p3V1(:), p3V2(:), p3V3(:), p3E(:), p3Ne(:)

  REAL(DP), ALLOCATABLE, TARGET :: uGF_K(:,:,:,:,:)
  REAL(DP), ALLOCATABLE, TARGET :: uCF_K(:,:,:,:,:)
  REAL(DP), ALLOCATABLE, TARGET :: uPF_K(:,:,:,:,:)

  REAL(DP) :: Gamma_IDEAL

  INTEGER :: nX_K, nNodesX_K

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

  CHARACTER(64) :: EosTableName
  CHARACTER(5)  :: EOS

  REAL(DP) :: Timer_old, Timer_new, Timer_new_new

  TimeIt_Euler = .TRUE.
  CALL InitializeTimers_Euler

  EOS = 'IDEAL'

  L = 1.0e2_DP * Kilometer

  nNodesX = [ 3, 3, 3 ]
  nX      = [ 32, 32, 32 ]
  swX     = [ 1, 1, 1 ]
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

  IF( TRIM( EOS ) .EQ. 'IDEAL' )THEN

    Gamma_IDEAL = 4.0_DP / 3.0_DP

    CALL InitializeEquationOfState &
           ( EquationOfState_Option = TRIM( EOS ), &
             Gamma_IDEAL_Option = Gamma_IDEAL, &
             Verbose_Option = .TRUE. )

  ELSE

    EosTableName = 'wl-EOS-SFHo-25-50-100.h5'

    CALL InitializeEquationOfState &
           ( EquationOfState_Option = TRIM( EOS ), &
             EquationOfStateTableName_Option &
               = TRIM( EosTableName ), &
             Verbose_Option = .TRUE. )

  END IF

  ! --- Define primitives ---

  D0  = 1.0e12_DP * UnitsD
  Amp = 0.1_DP * D0
  V1  = 3.0e4_DP * UnitsV1
  V2  = 0.0_DP * UnitsV2
  V3  = 0.0_DP * UnitsV3
  P   = 0.01_DP * D0 * SpeedOfLight**2
  Ye  = 0.0_DP

  ASSOCIATE( X1C => MeshX(1) % Center, dX1 => MeshX(1) % Width )

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

      IF( TRIM( EOS ) .EQ. 'TABLE' )THEN

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

      ELSE

        uPF(iNX,iX1,iX2,iX3,iPF_E) = P / ( Gamma_IDEAL - One )

      END IF

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

  END ASSOCIATE ! X1_C, dX1

  nX_K      = PRODUCT( nX )
  nNodesX_K = nDOFX * nX_K

  ALLOCATE( uGF_K(nDOFX,nX(1),nX(2),nX(3),nGF) )
  ALLOCATE( uCF_K(nDOFX,nX(1),nX(2),nX(3),nCF) )
  ALLOCATE( uPF_K(nDOFX,nX(1),nX(2),nX(3),nPF) )

  uGF_K = uGF(:,1:nX(1),1:nX(2),1:nX(3),:)
  uCF_K = uCF(:,1:nX(1),1:nX(2),1:nX(3),:)
  uPF_K = uPF(:,1:nX(1),1:nX(2),1:nX(3),:)

  Gm11(1:nNodesX_K) => uGF_K(:,:,:,:,iGF_Gm_dd_11)
  Gm22(1:nNodesX_K) => uGF_K(:,:,:,:,iGF_Gm_dd_22)
  Gm33(1:nNodesX_K) => uGF_K(:,:,:,:,iGF_Gm_dd_33)

  uD (1:nNodesX_K) => uCF_K(:,:,:,:,iCF_D )
  uS1(1:nNodesX_K) => uCF_K(:,:,:,:,iCF_S1)
  uS2(1:nNodesX_K) => uCF_K(:,:,:,:,iCF_S2)
  uS3(1:nNodesX_K) => uCF_K(:,:,:,:,iCF_S3)
  uE (1:nNodesX_K) => uCF_K(:,:,:,:,iCF_E )
  uNe(1:nNodesX_K) => uCF_K(:,:,:,:,iCF_Ne)

  p1D (1:nNodesX_K) => uPF_K(:,:,:,:,iPF_D )
  p1V1(1:nNodesX_K) => uPF_K(:,:,:,:,iPF_V1)
  p1V2(1:nNodesX_K) => uPF_K(:,:,:,:,iPF_V2)
  p1V3(1:nNodesX_K) => uPF_K(:,:,:,:,iPF_V3)
  p1E (1:nNodesX_K) => uPF_K(:,:,:,:,iPF_E )
  p1Ne(1:nNodesX_K) => uPF_K(:,:,:,:,iPF_Ne)

  ALLOCATE( p2D (nNodesX_K) )
  ALLOCATE( p2V1(nNodesX_K) )
  ALLOCATE( p2V2(nNodesX_K) )
  ALLOCATE( p2V3(nNodesX_K) )
  ALLOCATE( p2E (nNodesX_K) )
  ALLOCATE( p2Ne(nNodesX_K) )

  ALLOCATE( p3D (nNodesX_K) )
  ALLOCATE( p3V1(nNodesX_K) )
  ALLOCATE( p3V2(nNodesX_K) )
  ALLOCATE( p3V3(nNodesX_K) )
  ALLOCATE( p3E (nNodesX_K) )
  ALLOCATE( p3Ne(nNodesX_K) )

#if defined(THORNADO_OMP_OL)
  !$OMP TARGET ENTER DATA &
  !$OMP MAP( to:    uD, uS1, uS2, uS3, uE, uNe, &
  !$OMP             Gm11, Gm22, Gm33 ) &
  !$OMP MAP( alloc: p2D, p2V1, p2V2, p2V3, p2E, p2Ne, &
  !$OMP             p3D, p3V1, p3V2, p3V3, p3E, p3Ne )
#elif defined(THORNADO_OACC)
  !$ACC ENTER DATA &
  !$ACC COPYIN(     uD, uS1, uS2, uS3, uE, uNe, &
  !$ACC             Gm11, Gm22, Gm33 ) &
  !$ACC CREATE(     p2D, p2V1, p2V2, p2V3, p2E, p2Ne, &
  !$ACC             p3D, p3V1, p3V2, p3V3, p3E, p3Ne )
#endif

  Timer_old = MPI_WTIME()

  CALL ComputePrimitive_Euler_Relativistic &
    ( uD , uS1 , uS2 , uS3 , uE , uNe , &
      p2D, p2V1, p2V2, p2V3, p2E, p2Ne, &
      Gm11, Gm22, Gm33 )

  Timer_old = MPI_WTIME() - Timer_old

  Timer_new = MPI_WTIME()

  CALL ComputePrimitive_Euler_Relativistic &
    ( uD , uS1 , uS2 , uS3 , uE , uNe , &
      p3D, p3V1, p3V2, p3V3, p3E, p3Ne, &
      Gm11, Gm22, Gm33 )

  Timer_new = MPI_WTIME() - Timer_new

#if defined(THORNADO_OMP_OL)
  !$OMP TARGET EXIT DATA &
  !$OMP MAP( from:    p2D, p2V1, p2V2, p2V3, p2E, p2Ne, &
  !$OMP               p3D, p3V1, p3V2, p3V3, p3E, p3Ne ) &
  !$OMP MAP( release: uD, uS1, uS2, uS3, uE, uNe, &
  !$OMP               Gm11, Gm22, Gm33 )
#elif defined(THORNADO_OACC)
  !$ACC EXIT DATA &
  !$ACC COPYOUT(      p2D, p2V1, p2V2, p2V3, p2E, p2Ne, &
  !$ACC               p3D, p3V1, p3V2, p3V3, p3E, p3Ne ) &
  !$ACC DELETE(       uD, uS1, uS2, uS3, uE, uNe, &
  !$ACC               Gm11, Gm22, Gm33 )
#endif

  PRINT*
  PRINT '(A,ES14.6E3,A)', 'Timer_old     = ', Timer_old, ' s'
  PRINT '(A,ES14.6E3,A)', 'Timer_new     = ', Timer_new, ' s'
  PRINT*
  PRINT*, 'Maximum relative error:'
  PRINT*, '-----------------------'
  PRINT*
  PRINT*, 'PF_D  (old): ', MAXVAL( ABS( ( p1D - p2D ) / p1D ) )
  PRINT*, 'PF_D  (new): ', MAXVAL( ABS( ( p1D - p3D ) / p1D ) )
  PRINT*, 'PF_V1 (old): ', MAXVAL( ABS( ( p1V1 - p2V1 ) / p1V1 ) )
  PRINT*, 'PF_V1 (new): ', MAXVAL( ABS( ( p1V1 - p3V1 ) / p1V1 ) )
  IF( ABS( V2 ) .GT. Zero )THEN
    PRINT*, 'PF_V2 (old): ', MAXVAL( ABS( ( p1V2 - p2V2 ) / p1V2 ) )
    PRINT*, 'PF_V2 (new): ', MAXVAL( ABS( ( p1V2 - p3V2 ) / p1V2 ) )
  END IF
  IF( ABS( V3 ) .GT. Zero )THEN
    PRINT*, 'PF_V3 (old): ', MAXVAL( ABS( ( p1V3 - p2V3 ) / p1V3 ) )
    PRINT*, 'PF_V3 (new): ', MAXVAL( ABS( ( p1V3 - p3V3 ) / p1V3 ) )
  END IF
  PRINT*, 'PF_E  (old): ', MAXVAL( ABS( ( p1E - p2E ) / p1E ) )
  PRINT*, 'PF_E  (new): ', MAXVAL( ABS( ( p1E - p3E ) / p1E ) )
  IF( TRIM( EOS ) .EQ. 'TABLE' )THEN
    PRINT*, 'PF_Ne (old): ', MAXVAL( ABS( ( p1Ne - p2Ne ) / p1Ne ) )
    PRINT*, 'PF_Ne (new): ', MAXVAL( ABS( ( p1Ne - p3Ne ) / p1Ne ) )
  END IF
  PRINT*

  DEALLOCATE( p3Ne )
  DEALLOCATE( p3E  )
  DEALLOCATE( p3V3 )
  DEALLOCATE( p3V2 )
  DEALLOCATE( p3V1 )
  DEALLOCATE( p3D  )

  DEALLOCATE( p2Ne )
  DEALLOCATE( p2E  )
  DEALLOCATE( p2V3 )
  DEALLOCATE( p2V2 )
  DEALLOCATE( p2V1 )
  DEALLOCATE( p2D  )

  NULLIFY( p1Ne  )
  NULLIFY( p1E   )
  NULLIFY( p1V3  )
  NULLIFY( p1V2  )
  NULLIFY( p1V1  )
  NULLIFY( p1D   )

  NULLIFY( uNe  )
  NULLIFY( uE   )
  NULLIFY( uS3  )
  NULLIFY( uS2  )
  NULLIFY( uS1  )
  NULLIFY( uD   )

  NULLIFY( Gm33 )
  NULLIFY( Gm22 )
  NULLIFY( Gm11 )

  DEALLOCATE( uPF_K )
  DEALLOCATE( uCF_K )
  DEALLOCATE( uGF_K )

  CALL FinalizeEquationOfState

  CALL FinalizeProgram

  CALL FinalizeTimers_Euler

END PROGRAM ComputePrimitiveTest
