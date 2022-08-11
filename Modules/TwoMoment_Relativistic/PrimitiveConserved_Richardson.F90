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
  USE UtilitiesModule, ONLY: &
    WriteVector
  USE TwoMoment_ClosureModule, ONLY: &
    InitializeClosure_TwoMoment
  USE TwoMoment_UtilitiesModule_Relativistic, ONLY: &
    ComputePrimitive_TwoMoment_Vector_Richardson, &
    ComputePrimitive_TwoMoment, &
    ComputeConserved_TwoMoment, &
    Flux_X1,                    &
    Flux_X2,                    &
    Flux_X3                    
  USE TwoMoment_TimersModule_Relativistic, ONLY: &
    InitializeTimers, &
    FinalizeTimers, &
    TimersStart, &
    TimersStop, &
    Timer_ComputePrimative


  IMPLICIT NONE

  INTEGER  :: nNodes, nPoints, iPoint
  INTEGER  :: nE, nX(3), nZ_G, nX_G, nE_G, nSpecies, nNodesZ(4)
  INTEGER  :: iNode, iE, iX1, iX2, iX3, iS, iN_E, iN_X, iZ, i, j
  REAL(DP) :: eL, eR, xL(3), xR(3)
  REAL(DP) , DIMENSION(:), ALLOCATABLE :: V1, V2, V3, Ones, Zeros, B1, B2, B3, alp
  REAL(DP) :: absI, nVec(3), Vmax, absV, Vvec(3), Vsq, W
  REAL(DP) :: Bmax, absB, Bvec(3), F1(4), F2(4), F3(4), pep
  REAL(DP) :: I1T, I2T, I3T, absIi, Ii(1:3), Vi(1:3), gg(1:3,1:3)
  REAL(DP), DIMENSION(:), ALLOCATABLE :: D, I1, I2, I3, N, G1, G2, G3 

  INTEGER, ALLOCATABLE :: nIterations(:)


  INTEGER,  DIMENSION(:), ALLOCATABLE :: PositionIndexZ


    CALL InitializeTimers

  nNodes = 2
 gg= 0.0_DP
gg(1,1) = 1.0_DP
gg(2,2) = 1.0_DP
gg(3,3) = 1.0_DP
  nSpecies = 6
  nNodesZ = [2,2,2,2]
  nX = [ 16, 16, 16 ]
  xL = [ Zero, Zero, Zero ]
  xR = [ One,  One,  One  ]


  nE = 16
  eL = Zero
  eR = One

  nX_G = PRODUCT( nX(1:3) * nNodesZ(2:4))
  nE_G =  nE * nNodesZ(1)
  nZ_G =  nX_G * nE_G * nSpecies

    ALLOCATE( V1(nX_G) )
    ALLOCATE( V2(nX_G) )
    ALLOCATE( V3(nX_G) )
    ALLOCATE( Ones(nX_G)) 
    ALLOCATE(Zeros(nX_G))
    ALLOCATE( alp(nX_G) )
    ALLOCATE( B1(nX_G) )
    ALLOCATE( B2(nX_G) )
    ALLOCATE( B3(nX_G) )

  absV=0.2_DP
  CALL RANDOM_NUMBER(Vvec)
  Vvec = 2.0_DP * ( Vvec - 0.5_DP )
  Vvec = Vvec / SQRT( DOT_PRODUCT( Vvec, Vvec ) )
  V1=absV*Vvec(1)
  V2=absV*Vvec(2)
  V3=absV*Vvec(3)

  Bmax=0.0_DP
  CALL RANDOM_NUMBER(absB)
  CALL RANDOM_NUMBER(Bvec)
  alp = 1.0_DP
  absB=Bmax
  Bvec = 2.0_DP * ( Bvec - 0.5_DP )
  Bvec = Bvec / SQRT( DOT_PRODUCT( Bvec, Bvec ) )
  B1=absB*Bvec(1)
  B2=absB*Bvec(2)
  B3=absB*Bvec(3)
  Ones=1.0_DP
  Zeros = 0.0_DP




    ALLOCATE( PositionIndexZ(nZ_G) )

    iZ = 0
    DO iS   = 1, nSpecies
    DO iN_X = 1, nX_G
    DO iN_E = 1, nE_G

      iZ = iZ + 1

      PositionIndexZ(iZ) = iN_X

    END DO
    END DO
    END DO

    ALLOCATE( D(nZ_G) )
    ALLOCATE( I1(nZ_G) )
    ALLOCATE( I2(nZ_G) )
    ALLOCATE( I3(nZ_G) )
    ALLOCATE( N(nZ_G) )
    ALLOCATE( G1(nZ_G) )
    ALLOCATE( G2(nZ_G) )
    ALLOCATE( G3(nZ_G) )
!absIi = 0.0_DP

    DO iZ = 1, nZ_G
      Vsq = V1(PositionIndexZ(iZ))**2 &
      + V2(PositionIndexZ(iZ))**2 &
      + V3(PositionIndexZ(iZ))**2
      W = 1.0_DP / SQRT( 1.0_DP - Vsq )
      CALL RANDOM_NUMBER( D(iZ) )

      CALL RANDOM_NUMBER( absI )

      absI = absI * D(iZ)

      CALL RANDOM_NUMBER( nVec )

      nVec = 2.0_DP * ( nVec - 0.5_DP )
      nVec = nVec/ SQRT( DOT_PRODUCT( nVec, nVec )  )

      I1T = absI * nVeC(1)
      I2T = absI * nVeC(2)
      I3T = absI * nVeC(3)


      I1(iZ) = ( 1.0_DP + ( ( W - 1.0_DP ) / Vsq ) * V1(PositionIndexZ(iZ)) * V1(PositionIndexZ(iZ)) ) * I1T &
             + ( ( W - 1.0_DP ) / Vsq ) * V1(PositionIndexZ(iZ)) * V2(PositionIndexZ(iZ)) * I2T &
             + ( ( W - 1.0_DP ) / Vsq ) * V1(PositionIndexZ(iZ)) * V3(PositionIndexZ(iZ)) * I3T  
      I2(iZ) = ( 1.0_DP + ( ( W - 1.0_DP ) / Vsq ) * V2(PositionIndexZ(iZ)) * V2(PositionIndexZ(iZ)) ) * I2T &
             + ( ( W - 1.0_DP ) / Vsq ) * V1(PositionIndexZ(iZ)) * V2(PositionIndexZ(iZ)) * I1T &
             + ( ( W - 1.0_DP ) / Vsq ) * V2(PositionIndexZ(iZ)) * V3(PositionIndexZ(iZ)) * I3T  
      I3(iZ) = ( 1.0_DP + ( ( W - 1.0_DP ) / Vsq ) * V3(PositionIndexZ(iZ)) * V3(PositionIndexZ(iZ)) ) * I3T &
             + ( ( W - 1.0_DP ) / Vsq ) * V1(PositionIndexZ(iZ)) * V3(PositionIndexZ(iZ)) * I1T &
             + ( ( W - 1.0_DP ) / Vsq ) * V2(PositionIndexZ(iZ)) * V3(PositionIndexZ(iZ)) * I2T  
  
      CALL ComputeConserved_TwoMoment &
               ( D(iZ), I1(iZ), I2(iZ), I3(iZ), &
                 N(iZ), G1(iZ), G2(iZ), G3(iZ), &
                 V1(PositionIndexZ(iZ)), V2(PositionIndexZ(iZ)), V3(PositionIndexZ(iZ)),&
                 1.0_DP, 1.0_DP, 1.0_DP, &
                 0.0_DP, 0.0_DP,0.0_DP, &
                 alp(PositionIndexZ(iZ)), B1(PositionIndexZ(iZ)), &
                 B2(PositionIndexZ(iZ)),B3(PositionIndexZ(iZ)) )


    END DO

  CALL WriteVector(nZ_G,DBLE(D), 'Dinitial.dat')
  CALL WriteVector(nZ_G,DBLE(I1), 'I1initial.dat')
  CALL WriteVector(nZ_G,DBLE(I2), 'I2initial.dat')
  CALL WriteVector(nZ_G,DBLE(i3), 'I3initial.dat')
  ! --- Compute Primitive From Conserved ---


  ALLOCATE( nIterations(nZ_G) )
 
#if defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( N, G1, G2, G3, &
    !$ACC         D, I1, I2, I3, &
    !$ACC         V1, V2, V3, Ones, Zeros, &
    !$ACC         alp, B1, B2, B3, PositionIndexZ, nIterations  )
#endif





    CALL TimersStart( Timer_ComputePrimative )



 
       CALL ComputePrimitive_TwoMoment_Vector_Richardson &
           ( N, G1, G2, G3, &
             D, I1, I2, I3, &
             V1, &
             V2, &
             V3, &
             Ones, &
             Ones, &
             Ones, &
             Zeros, &
             Zeros, &
             Zeros, &
             alp, B1, B2, B3, &
             PositionIndexZ, &
             nIterations )     

    CALL TimersStop( Timer_ComputePrimative )


#if defined(THORNADO_OACC)
    !$ACC UPDATE HOST( N, G1, G2, G3, &
    !$ACC                D, I1, I2, I3, &
    !$ACC                V1, V2, V3, Ones, Zeros, &
    !$ACC                alp, B1, B2, B3, PositionIndexZ, nIterations  )
#endif
   


  CALL WriteVector( nZ_G, DBLE( nIterations ), 'nIterations.dat' )

  CALL WriteVector(nZ_G,DBLE(D), 'Dfinal.dat')
  CALL WriteVector(nZ_G,DBLE(I1), 'I1final.dat')
  CALL WriteVector(nZ_G,DBLE(I2), 'I2final.dat')
  CALL WriteVector(nZ_G,DBLE(I3), 'I3final.dat')
  
OPEN(UNIT=14,FILE="vel.txt",FORM="FORMATTED",STATUS="UNKNOWN",ACTION="WRITE")

  WRITE(UNIT=14, FMT=*) absV, V1(1), V2(1), V3(1)

print*,absV, V1(1), V2(1), V3(1)
print*,N(1), G1(1), G2(1), G3(1)

OPEN(UNIT=14,FILE="N1.txt",FORM="FORMATTED",STATUS="UNKNOWN",ACTION="WRITE")

  WRITE(UNIT=14, FMT=*) N(1), G1(1), G2(1), G3(1)
  DEALLOCATE( nIterations )

    CALL FinalizeTimers
END PROGRAM PrimitiveConserved
