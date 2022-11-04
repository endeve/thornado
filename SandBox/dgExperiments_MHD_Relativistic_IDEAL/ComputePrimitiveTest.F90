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
    InitializeProgram_Basic, &
    FinalizeProgram_Basic
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
    SpeedOfLight, &
    Gauss
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    uGF, &
    nGF, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_Alpha,    &
    iGF_Beta_1,   &
    iGF_Beta_2,   &
    iGF_Beta_3
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX
  USE MagnetofluidFieldsModule, ONLY: &
    uCM, &
    nCM, &
    iCM_D, &
    iCM_S1, &
    iCM_S2, &
    iCM_S3, &
    iCM_E, &
    iCM_Ne, &
    iCM_B1, &
    iCM_B2, &
    iCM_B3, &
    iCM_Chi, &
    uPM, &
    nPM, &
    iPM_D, &
    iPM_V1, &
    iPM_V2, &
    iPM_V3, &
    iPM_E, &
    iPM_Ne, &
    iPM_B1, &
    iPM_B2, &
    iPM_B3, &
    iPM_Chi, &
    uAM, &
    nAM, &
    iAM_P, &
    iAM_T, &
    iAM_E, &
    iAM_Ye
  USE EquationOfStateModule, ONLY: &
    InitializeEquationOfState, &
    FinalizeEquationOfState, &
    ComputeTemperatureFromPressure, &
    ComputeThermodynamicStates_Primitive
  USE MHD_UtilitiesModule_Relativistic, ONLY: &
    ComputeConserved_MHD_Relativistic, &
    ComputePrimitive_MHD_Relativistic

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  REAL(DP) :: xL(3), xU(3), X1, L
  INTEGER  :: nNodesX(3), nX(3)
  INTEGER  :: swX(3)

  INTEGER :: iNX, iX1, iX2, iX3, iNX1, iNX2, iNX3

  REAL(DP), POINTER, CONTIGUOUS :: uD(:), uS1(:), uS2(:), uS3(:), uE(:), uNe(:), uB1(:), uB2(:), uB3(:), uChi(:)
  REAL(DP), POINTER, CONTIGUOUS :: Gm11(:), Gm22(:), Gm33(:), Alpha(:), Beta1(:), Beta2(:), Beta3(:)
  REAL(DP), POINTER, CONTIGUOUS :: p1D(:), p1V1(:), p1V2(:), &
                                   p1V3(:), p1E(:), p1Ne(:), &
                                   p1B1(:), p1B2(:), p1B3(:), p1Chi(:)
  REAL(DP), ALLOCATABLE :: p2D(:), p2V1(:), p2V2(:), p2V3(:), p2E(:), p2Ne(:), p2B1(:), p2B2(:), p2B3(:), p2Chi(:)
  REAL(DP), ALLOCATABLE :: p3D(:), p3V1(:), p3V2(:), p3V3(:), p3E(:), p3Ne(:), p3B1(:), p3B2(:), p3B3(:), p3Chi(:)

  REAL(DP), ALLOCATABLE, TARGET :: uGF_K(:,:,:,:,:)
  REAL(DP), ALLOCATABLE, TARGET :: uCM_K(:,:,:,:,:)
  REAL(DP), ALLOCATABLE, TARGET :: uPM_K(:,:,:,:,:)

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
  REAL(DP) :: B1, B2, B3, Chi

  CHARACTER(64) :: EosTableName
  CHARACTER(5)  :: EOS

  EOS = 'IDEAL'

  L = 1.0e2_DP * Kilometer

  nNodesX = [ 1, 1, 1 ]
  nX      = [ 8, 1, 1 ]
  swX     = [ 1, 0, 0 ]
  xL      = [ Zero, Zero, Zero ]
  xU      = [ One,  One,  One  ] * L

  CALL InitializeProgram_Basic &
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
             = .FALSE., &
           UseMHD_Option &
             = .TRUE. )

  CALL ComputeGeometryX &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF )

  IF( TRIM( EOS ) .EQ. 'IDEAL' )THEN

    Gamma_IDEAL = 4.0_DP / 3.0_DP

    CALL InitializeEquationOfState &
           ( EquationOfState_Option = TRIM( EOS ), &
             Gamma_IDEAL_Option = Gamma_IDEAL, &
             Verbose_Option = .TRUE. )

  END IF

  ! --- Define primitives ---

  D0  = 1.0_DP
  V1  = 0.5_DP
  V2  = 0.0_DP
  V3  = 0.0_DP
  P   = 0.1_DP
  Ye  = 0.0_DP
  B1  = 0.0_DP
  B2  = 0.0_DP
  B3  = 0.0_DP
  Chi = 0.0_DP

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

      PRINT*, 'iNX', iNX

      uPM(iNX,iX1,iX2,iX3,iPM_D ) = D0
      uPM(iNX,iX1,iX2,iX3,iPM_V1) = V1
      uPM(iNX,iX1,iX2,iX3,iPM_V2) = V2
      uPM(iNX,iX1,iX2,iX3,iPM_V3) = V3

      uAM(iNX,iX1,iX2,iX3,iAM_P ) = P
      uAM(iNX,iX1,iX2,iX3,iAM_Ye) = Ye

      uPM(iNX,iX1,iX2,iX3,iPM_E) = P / ( Gamma_IDEAL - One )

      uPM(iNX,iX1,iX2,iX3,iPM_B1) = B1
      uPM(iNX,iX1,iX2,iX3,iPM_B2) = B2
      uPM(iNX,iX1,iX2,iX3,iPM_B3) = B3

      CALL ComputeConserved_MHD_Relativistic &
             ( uPM(iNX,iX1,iX2,iX3,iPM_D  ), &
               uPM(iNX,iX1,iX2,iX3,iPM_V1 ), &
               uPM(iNX,iX1,iX2,iX3,iPM_V2 ), &
               uPM(iNX,iX1,iX2,iX3,iPM_V3 ), &
               uPM(iNX,iX1,iX2,iX3,iPM_E  ), &
               uPM(iNX,iX1,iX2,iX3,iPM_Ne ), &
               uPM(iNX,iX1,iX2,iX3,iPM_B1 ), &
               uPM(iNX,iX1,iX2,iX3,iPM_B2 ), &
               uPM(iNX,iX1,iX2,iX3,iPM_B3 ), &
               uPM(iNX,iX1,iX2,iX3,iPM_Chi), &
               uCM(iNX,iX1,iX2,iX3,iCM_D  ), &
               uCM(iNX,iX1,iX2,iX3,iCM_S1 ), &
               uCM(iNX,iX1,iX2,iX3,iCM_S2 ), &
               uCM(iNX,iX1,iX2,iX3,iCM_S3 ), &
               uCM(iNX,iX1,iX2,iX3,iCM_E  ), &
               uCM(iNX,iX1,iX2,iX3,iCM_Ne ), &
               uCM(iNX,iX1,iX2,iX3,iCM_B1 ), &
               uCM(iNX,iX1,iX2,iX3,iCM_B2 ), &
               uCM(iNX,iX1,iX2,iX3,iCM_B3 ), &
               uCM(iNX,iX1,iX2,iX3,iCM_Chi), &
               uGF(iNX,iX1,iX2,iX3,iGF_Gm_dd_11), &
               uGF(iNX,iX1,iX2,iX3,iGF_Gm_dd_22), &
               uGF(iNX,iX1,iX2,iX3,iGF_Gm_dd_33), &
               uGF(iNX,iX1,iX2,iX3,iGF_Alpha), &
               uGF(iNX,iX1,iX2,iX3,iGF_Beta_1), &
               uGF(iNX,iX1,iX2,iX3,iGF_Beta_2), &
               uGF(iNX,iX1,iX2,iX3,iGF_Beta_3), &
               uAM(iNX,iX1,iX2,iX3,iAM_P), .FALSE. )

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
  ALLOCATE( uCM_K(nDOFX,nX(1),nX(2),nX(3),nCM) )
  ALLOCATE( uPM_K(nDOFX,nX(1),nX(2),nX(3),nPM) )

  uGF_K = uGF(:,1:nX(1),1:nX(2),1:nX(3),:)
  uCM_K = uCM(:,1:nX(1),1:nX(2),1:nX(3),:)
  uPM_K = uPM(:,1:nX(1),1:nX(2),1:nX(3),:)

  Gm11(1:nNodesX_K)  => uGF_K(:,:,:,:,iGF_Gm_dd_11)
  Gm22(1:nNodesX_K)  => uGF_K(:,:,:,:,iGF_Gm_dd_22)
  Gm33(1:nNodesX_K)  => uGF_K(:,:,:,:,iGF_Gm_dd_33)
  Alpha(1:nNodesX_K) => uGF_K(:,:,:,:,iGF_Alpha   )
  Beta1(1:nNodesX_K) => uGF_K(:,:,:,:,iGF_Beta_1  )
  Beta2(1:nNodesX_K) => uGF_K(:,:,:,:,iGF_Beta_2  )
  Beta3(1:nNodesX_K) => uGF_K(:,:,:,:,iGF_Beta_3  )

  uD  (1:nNodesX_K) => uCM_K(:,:,:,:,iCM_D )
  uS1 (1:nNodesX_K) => uCM_K(:,:,:,:,iCM_S1)
  uS2 (1:nNodesX_K) => uCM_K(:,:,:,:,iCM_S2)
  uS3 (1:nNodesX_K) => uCM_K(:,:,:,:,iCM_S3)
  uE  (1:nNodesX_K) => uCM_K(:,:,:,:,iCM_E )
  uNe (1:nNodesX_K) => uCM_K(:,:,:,:,iCM_Ne)
  uB1 (1:nNodesX_K) => uCM_K(:,:,:,:,iCM_B1)
  uB2 (1:nNodesX_K) => uCM_K(:,:,:,:,iCM_B2)
  uB3 (1:nNodesX_K) => uCM_K(:,:,:,:,iCM_B3)
  uChi(1:nNodesX_K) => uCM_K(:,:,:,:,iCM_Chi)

  p1D  (1:nNodesX_K) => uPM_K(:,:,:,:,iPM_D )
  p1V1 (1:nNodesX_K) => uPM_K(:,:,:,:,iPM_V1)
  p1V2 (1:nNodesX_K) => uPM_K(:,:,:,:,iPM_V2)
  p1V3 (1:nNodesX_K) => uPM_K(:,:,:,:,iPM_V3)
  p1E  (1:nNodesX_K) => uPM_K(:,:,:,:,iPM_E )
  p1Ne (1:nNodesX_K) => uPM_K(:,:,:,:,iPM_Ne)
  p1B1 (1:nNodesX_K) => uPM_K(:,:,:,:,iPM_B1)
  p1B2 (1:nNodesX_K) => uPM_K(:,:,:,:,iPM_B2)
  p1B3 (1:nNodesX_K) => uPM_K(:,:,:,:,iPM_B3)
  p1Chi(1:nNodesX_K) => uPM_K(:,:,:,:,iPM_Chi)

  ALLOCATE( p2D  (nNodesX_K) )
  ALLOCATE( p2V1 (nNodesX_K) )
  ALLOCATE( p2V2 (nNodesX_K) )
  ALLOCATE( p2V3 (nNodesX_K) )
  ALLOCATE( p2E  (nNodesX_K) )
  ALLOCATE( p2Ne (nNodesX_K) )
  ALLOCATE( p2B1 (nNodesX_K) )
  ALLOCATE( p2B2 (nNodesX_K) )
  ALLOCATE( p2B3 (nNodesX_K) )
  ALLOCATE( p2Chi(nNodesX_K) )

  ALLOCATE( p3D  (nNodesX_K) )
  ALLOCATE( p3V1 (nNodesX_K) )
  ALLOCATE( p3V2 (nNodesX_K) )
  ALLOCATE( p3V3 (nNodesX_K) )
  ALLOCATE( p3E  (nNodesX_K) )
  ALLOCATE( p3Ne (nNodesX_K) )
  ALLOCATE( p3B1 (nNodesX_K) )
  ALLOCATE( p3B2 (nNodesX_K) )
  ALLOCATE( p3B3 (nNodesX_K) )
  ALLOCATE( p3Chi(nNodesX_K) )

  CALL ComputePrimitive_MHD_Relativistic &
    ( uD , uS1 , uS2 , uS3 , uE , uNe , &
      uB1, uB2, uB3, uChi, &
      p2D, p2V1, p2V2, p2V3, p2E, p2Ne, &
      p2B1, p2B2, p2B3, p2Chi, &
      Gm11, Gm22, Gm33, Alpha, Beta1, Beta2, Beta3, &
      .FALSE. )

  CALL ComputePrimitive_MHD_Relativistic &
    ( uD , uS1 , uS2 , uS3 , uE , uNe , &
      uB1, uB2, uB3, uChi, &
      p3D, p3V1, p3V2, p3V3, p3E, p3Ne, &
      p3B1, p3B2, p3B3, p3Chi, &
      Gm11, Gm22, Gm33, Alpha, Beta1, Beta2, Beta3, &
      .FALSE. )

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

  DEALLOCATE( p3Chi )
  DEALLOCATE( p3B3  )
  DEALLOCATE( p3B2  )
  DEALLOCATE( p3B1  )
  DEALLOCATE( p3Ne  )
  DEALLOCATE( p3E   )
  DEALLOCATE( p3V3  )
  DEALLOCATE( p3V2  )
  DEALLOCATE( p3V1  )
  DEALLOCATE( p3D   )

  DEALLOCATE( p2Chi )
  DEALLOCATE( p2B3  )
  DEALLOCATE( p2B2  )
  DEALLOCATE( p2B1  )
  DEALLOCATE( p2Ne )
  DEALLOCATE( p2E  )
  DEALLOCATE( p2V3 )
  DEALLOCATE( p2V2 )
  DEALLOCATE( p2V1 )
  DEALLOCATE( p2D  )

  NULLIFY( p1Chi )
  NULLIFY( p1B3  )
  NULLIFY( p1B2  )
  NULLIFY( p1B1  )
  NULLIFY( p1Ne  )
  NULLIFY( p1E   )
  NULLIFY( p1V3  )
  NULLIFY( p1V2  )
  NULLIFY( p1V1  )
  NULLIFY( p1D   )

  NULLIFY( uChi )
  NULLIFY( uB3  )
  NULLIFY( uB2  )
  NULLIFY( uB1  )
  NULLIFY( uNe  )
  NULLIFY( uE   )
  NULLIFY( uS3  )
  NULLIFY( uS2  )
  NULLIFY( uS1  )
  NULLIFY( uD   )

  NULLIFY( Gm33 )
  NULLIFY( Gm22 )
  NULLIFY( Gm11 )

  DEALLOCATE( uPM_K )
  DEALLOCATE( uCM_K )
  DEALLOCATE( uGF_K )

  CALL FinalizeEquationOfState

  CALL FinalizeProgram_Basic( UseMHD_Option = .TRUE. )

END PROGRAM ComputePrimitiveTest
