PROGRAM PrimitiveConserved

  USE KindModule, ONLY: &
    DP, Zero, One
  USE ProgramHeaderModule, ONLY: &
    iX_B0, iX_E0, iX_B1, iX_E1, &
    iE_B0, iE_E0, iE_B1, iE_E1, &
    nDOF
  USE ProgramInitializationModule, ONLY: &
    InitializeProgram, &
    FinalizeProgram
  USE ReferenceElementModuleX, ONLY: &
    InitializeReferenceElementX, &
    FinalizeReferenceElementX
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    InitializeReferenceElementX_Lagrange, &
    FinalizeReferenceElementX_Lagrange
  USE ReferenceElementModuleE, ONLY: &
    InitializeReferenceElementE, &
    FinalizeReferenceElementE
  USE ReferenceElementModuleE_Lagrange, ONLY: &
    InitializeReferenceElementE_Lagrange, &
    FinalizeReferenceElementE_Lagrange
  USE ReferenceElementModule, ONLY: &
    InitializeReferenceElement, &
    FinalizeReferenceElement
  USE ReferenceElementModule_Lagrange, ONLY: &
    InitializeReferenceElement_Lagrange, &
    FinalizeReferenceElement_Lagrange
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX
  USE GeometryComputationModuleE, ONLY: &
    ComputeGeometryE
  USE GeometryFieldsModule, ONLY: &
    uGF
  USE GeometryFieldsModuleE, ONLY: &
    uGE
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    uPR, iPR_D, iPR_I1, iPR_I2, iPR_I3, &
    uCR, iCR_N, iCR_G1, iCR_G2, iCR_G3
  USE InputOutputModuleHDF, ONLY: &
    WriteFieldsHDF
  USE UtilitiesModule, ONLY: &
    WriteVector
  USE TwoMoment_ClosureModule, ONLY: &
    InitializeClosure_TwoMoment
  USE TwoMoment_UtilitiesModule_Relativistic, ONLY: &
    ComputePrimitive_TwoMoment, &
    ComputeConserved_TwoMoment, &
    Flux_X1,                    &
    Flux_X2,                    &
    Flux_X3                    

  IMPLICIT NONE

  INTEGER  :: nNodes, nPoints, iPoint
  INTEGER  :: nE, nX(3)
  INTEGER  :: iNode, iE, iX1, iX2, iX3, iS
  REAL(DP) :: eL, eR, xL(3), xR(3)
  REAL(DP) :: D, I1, I2, I3, V1, V2, V3
  REAL(DP) :: absI, nVec(3), Vmax, absV, Vvec(3), Vsq, W
  REAL(DP) :: Bmax, absB, Bvec(3), B1, B2, B3, alp, F1(4), F2(4), F3(4), pep
  INTEGER, ALLOCATABLE :: nIterations(:)

  nNodes = 2

  nX = [ 24, 24, 24 ]
  xL = [ Zero, Zero, Zero ]
  xR = [ One,  One,  One  ]


  nE = 8
  eL = Zero
  eR = One

  Bmax=0.0_DP
  CALL RANDOM_NUMBER(absB)
  CALL RANDOM_NUMBER(Bvec)
  !absV=absV*Vmax
  absB=Bmax
  Bvec = 2.0_DP * ( Bvec - 0.5_DP )
  Bvec = Bvec / SQRT( DOT_PRODUCT( Bvec, Bvec ) )
  B1=absB*Bvec(1)
  B2=absB*Bvec(2)
  B3=absB*Bvec(3)
  alp = 1.0_DP

  Vmax=0.9_DP
  CALL RANDOM_NUMBER(absV)
  CALL RANDOM_NUMBER(Vvec(1))
  !CALL RANDOM_NUMBER(Vvec)
  !absV=absV*Vmax
  !absV=Vmax
  !Vvec = 2.0_DP * ( Vvec - 0.5_DP )
  !Vvec = Vvec / SQRT( DOT_PRODUCT( Vvec, Vvec ) )

  Vvec(1) = 2.0_DP * ( Vvec(1) - 0.5_DP )
  Vvec(1) = Vvec(1) / SQRT( Vvec(1) * Vvec(1) )
  V1 = absV * Vmax * Vvec(1)
  V2 = 0.0_DP
  V3 = 0.0_DP
  print*, V1
  !V1=absV*Vvec(1)
  !V2=absV*Vvec(2)
  !V3=absV*Vvec(3)
  Vsq = V1**2 + V2**2 + V3**2
  W = 1.0_DP / SQRT( 1.0_DP - Vsq )
  CALL InitializeProgram &
         ( ProgramName_Option &
             = 'PrimitiveConserved', &
           nX_Option &
             = nX, &
           swX_Option &
             = [ 0, 0, 0 ], &
           bcX_Option &
             = [ 0, 0, 0 ], &
           xL_Option &
             = xL, &
           xR_Option &
             = xR, &
           nE_Option &
             = nE, &
           eL_Option &
             = eL, &
           eR_Option &
             = eR, &
           nNodes_Option &
             = nNodes, &
           CoordinateSystem_Option &
             = 'CARTESIAN', &
           nSpecies_Option &
             = 1, &
           BasicInitialization_Option &
             = .TRUE. )

  ! --- Position Space Reference Element and Geometry ---

  CALL InitializeReferenceElementX

  CALL InitializeReferenceElementX_Lagrange

  CALL ComputeGeometryX &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF )

  ! --- Energy Space Reference Element and Geometry ---

  CALL InitializeReferenceElementE

  CALL InitializeReferenceElementE_Lagrange

  CALL ComputeGeometryE &
         ( iE_B0, iE_E0, iE_B1, iE_E1, uGE )

  ! --- Phase Space Reference Element ---

  CALL InitializeReferenceElement

  CALL InitializeReferenceElement_Lagrange

  ! --- Initialize Moment Closure ---

  CALL InitializeClosure_TwoMoment

  ! --- Initialize Primitive Radiation Fields ---


  DO iS  = 1, nSpecies
  DO iX3 = iX_B0(3), iX_E0(3)
  DO iX2 = iX_B0(2), iX_E0(2)
  DO iX1 = iX_B0(1), iX_E0(1)

    DO iE = iE_B0, iE_E0

      DO iNode = 1, nDOF

        ! --- Number Density: Random Number in [0,1) ---

        CALL RANDOM_NUMBER( D )

        uPR(iNode,iE,iX1,iX2,iX3,iPR_D,iS) = D

       ! uPR(iNode,iE,iX1,iX2,iX3,iPR_D,iS) = 1.2_DP * D
        ! --- Number Flux: Realizable with Random Direction ---

        CALL RANDOM_NUMBER( absI )

        !absI = absI * ( One - D ) * D
        absI = absI * D

        !CALL RANDOM_NUMBER( nVec )

        !nVec = 2.0_DP * ( nVec - 0.5_DP )
!        nVec = nVec / SQRT( DOT_PRODUCT( nVec, nVec ) )

        CALL RANDOM_NUMBER( nVec(1) )

        nVec(1) = 2.0_DP * ( nVec(1) - 0.5_DP )
        nVec(1) = nVec(1) / SQRT( nVec(1) * nVec(1)  )

        I1 = absI * nVec(1)
        I2 = 0.0_DP
        I3 = 0.0_DP
!        I2 = absI * nVec(2)
 !       I3 = absI * nVec(3)

        uPR(iNode,iE,iX1,iX2,iX3,iPR_I1,iS) = ( 1.0_DP + ( ( W - 1.0_DP ) / Vsq ) * V1 * V1 ) * I1 &
                                            + ( ( W - 1.0_DP ) / Vsq ) * V1 * V2 * I2 &
                                            + ( ( W - 1.0_DP ) / Vsq ) * V1 * V3 * I3  
        uPR(iNode,iE,iX1,iX2,iX3,iPR_I2,iS) = ( 1.0_DP + ( ( W - 1.0_DP ) / Vsq ) * V2 * V2 ) * I2 &
                                            + ( ( W - 1.0_DP ) / Vsq ) * V1 * V2 * I1 &
                                            + ( ( W - 1.0_DP ) / Vsq ) * V2 * V3 * I3  
        uPR(iNode,iE,iX1,iX2,iX3,iPR_I3,iS) = ( 1.0_DP + ( ( W - 1.0_DP ) / Vsq ) * V3 * V3 ) * I3 &
                                            + ( ( W - 1.0_DP ) / Vsq ) * V1 * V3 * I1 &
                                            + ( ( W - 1.0_DP ) / Vsq ) * V2 * V3 * I2  

 
        CALL ComputeConserved_TwoMoment &
               ( uPR(iNode,iE,iX1,iX2,iX3,iPR_D ,iS), &
                 uPR(iNode,iE,iX1,iX2,iX3,iPR_I1,iS), &
                 uPR(iNode,iE,iX1,iX2,iX3,iPR_I2,iS), &
                 uPR(iNode,iE,iX1,iX2,iX3,iPR_I3,iS), &
                 uCR(iNode,iE,iX1,iX2,iX3,iCR_N ,iS), &
                 uCR(iNode,iE,iX1,iX2,iX3,iCR_G1,iS), &
                 uCR(iNode,iE,iX1,iX2,iX3,iCR_G2,iS), &
                 uCR(iNode,iE,iX1,iX2,iX3,iCR_G3,iS), &
                 V1, V2, V3,                          &
                 1.0_DP, 1.0_DP, 1.0_DP ,             &
                 0.0_DP, 0.0_DP, 0.0_DP,              &
                 alp, B1, B2, B3 )

      END DO

    END DO

  END DO
  END DO
  END DO
  END DO

  CALL WriteFieldsHDF &
         ( Time = 0.0_DP, &
           WriteGF_Option = .TRUE., &
           WriteFF_Option = .TRUE., &
           WriteRF_Option = .TRUE. )


  ! --- Compute Primitive From Conserved ---

  nPoints = nSpecies * PRODUCT( iX_E0 - iX_B0 + 1 ) &
              * ( iE_E0 - iE_B0 + 1 ) * nDOF

  ALLOCATE( nIterations(nPoints) )

  iPoint = 0
  
  DO iS  = 1, nSpecies
  DO iX3 = iX_B0(3), iX_E0(3)
  DO iX2 = iX_B0(2), iX_E0(2)
  DO iX1 = iX_B0(1), iX_E0(1)

    DO iE = iE_B0, iE_E0

      DO iNode = 1, nDOF

        iPoint = iPoint + 1

        CALL ComputePrimitive_TwoMoment &
               ( uCR(iNode,iE,iX1,iX2,iX3,iCR_N ,iS), &
                 uCR(iNode,iE,iX1,iX2,iX3,iCR_G1,iS), &
                 uCR(iNode,iE,iX1,iX2,iX3,iCR_G2,iS), &
                 uCR(iNode,iE,iX1,iX2,iX3,iCR_G3,iS), &
                 uPR(iNode,iE,iX1,iX2,iX3,iPR_D ,iS), &
                 uPR(iNode,iE,iX1,iX2,iX3,iPR_I1,iS), &
                 uPR(iNode,iE,iX1,iX2,iX3,iPR_I2,iS), &
                 uPR(iNode,iE,iX1,iX2,iX3,iPR_I3,iS), &
                 V1, V2, V3, 1.0_DP, 1.0_DP, 1.0_DP, &
                 0.0_DP, 0.0_DP, 0.0_DP,              &
                 alp, B1, B2, B3, nIterations(iPoint) )
        


      END DO

    END DO

  END DO
  END DO
  END DO
  END DO




   

  CALL WriteFieldsHDF &
         ( Time = 0.0_DP, &
           WriteGF_Option = .TRUE., &
           WriteFF_Option = .TRUE., &
           WriteRF_Option = .TRUE. )

  CALL WriteVector( nPoints, DBLE( nIterations ), 'nIterations.dat' )

  DEALLOCATE( nIterations )

  CALL FinalizeReferenceElementX

  CALL FinalizeReferenceElementX_Lagrange

  CALL FinalizeReferenceElementE

  CALL FinalizeReferenceElementE_Lagrange

  CALL FinalizeReferenceElement

  CALL FinalizeReferenceElement_Lagrange


  CALL FinalizeProgram

END PROGRAM PrimitiveConserved
