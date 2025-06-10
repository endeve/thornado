PROGRAM MeshRefinementTest_TwoMoment

  USE KindModule, ONLY: &
    DP, Zero, One, Two, Pi, TwoPi, SqrtTiny
  USE UnitsModule, ONLY: &
    Kilometer
  USE ProgramHeaderModule, ONLY: &
    nDOFX, nDOFE, nNodesX, nDimsX
  USE ProgramInitializationModule, ONLY: &
    InitializeProgram, &
    FinalizeProgram
  USE ReferenceElementModuleX, ONLY: &
    InitializeReferenceElementX, &
    FinalizeReferenceElementX, &
    NodeNumberTableX
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
    FinalizeReferenceElement, &
    NodeNumberTable
  USE ReferenceElementModule_Lagrange, ONLY: &
    InitializeReferenceElement_Lagrange, &
    FinalizeReferenceElement_Lagrange
  USE MeshModule, ONLY : &
    MeshType, &
    CreateMesh, &
    DestroyMesh, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
     nGF, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33, iGF_SqrtGm
  USE GeometryComputationModule, ONLY: &
     ComputeGeometryX
  USE RadiationFieldsModule, ONLY: &
    iPR_D, iPR_I1, iPR_I2, iPR_I3, nPR, &
    iCR_N, iCR_G1, iCR_G2, iCR_G3, nCR
  USE UtilitiesModule, ONLY: &
    WriteVector
  USE TimersModule, ONLY: &
    TimersStart, &
    TimersStop
  USE TwoMoment_UtilitiesModule, ONLY: &
    ComputeConserved_TwoMoment
  USE TwoMoment_MeshRefinementModule, ONLY: &
    InitializeMeshRefinement_TwoMoment, &
    FinalizeMeshRefinement_TwoMoment, &
    RefineX_TwoMoment_SIMPLE, &
    RefineX_TwoMoment_CURVILINEAR, &
    CoarsenX_TwoMoment_SIMPLE, &
    CoarsenX_TwoMoment_CURVILINEAR

  IMPLICIT NONE

  INTEGER, PARAMETER :: refine_factor = 2

  REAL(DP) :: &
    Timer_Refine, &
    Timer_Coarsen, &
    Timer_Total

  CHARACTER(32) :: ProgramName
  CHARACTER(32) :: CoordinateSystem
  CHARACTER(32) :: CoordinateSystem_lc
  LOGICAL       :: UseSimpleMeshRefinement

  INTEGER       :: nNodes
  INTEGER       :: nX(3), bcX(3), swX(3)
  REAL(DP)      :: xL(3), xR(3), ZoomX(3)
  INTEGER       :: nE, bcE, swE
  REAL(DP)      :: eL, eR, ZoomE
  INTEGER       :: nSpecies

  CHARACTER(2)  :: MeshString
  CHARACTER(32) :: VectorName
  INTEGER       :: i, k
  INTEGER       :: iNodeX, iX1, iX2, iX3
  INTEGER       :: iNodeE, iE, iCR, iS
  INTEGER       :: iP_X1, iP_X2, iP_X3
  INTEGER       :: iNodeX1, iNodeX2, iNodeX3
  INTEGER       :: iFineX(3), iFine, nFine
  INTEGER       :: iX1_Fine, iX2_Fine, iX3_Fine
  INTEGER       :: iFineX1, iFineX2, iFineX3, nFineX(3)
  REAL(DP)      :: X1, X2, X3, R0
  REAL(DP), ALLOCATABLE :: X1_0(:), X2_0(:), X3_0(:)
  REAL(DP)      :: V_0(3)
  REAL(DP)      :: uPR(nPR), uCR(nCR)

  REAL(DP), ALLOCATABLE :: U_0    (:,:,:,:,:)
  REAL(DP), ALLOCATABLE :: U_Crse (:,:,:,:,:)
  REAL(DP), ALLOCATABLE :: U_1    (:,:,:,:,:,:)
  REAL(DP), ALLOCATABLE :: U_Fine (:,:,:,:,:,:)

  REAL(DP), ALLOCATABLE :: G_Crse(:,:,:,:,:)
  REAL(DP), ALLOCATABLE :: G_Fine(:,:,:,:,:)

  INTEGER :: nX_Crse(3)
  INTEGER :: nX_Fine(3)

  INTEGER :: iX_Crse_B1(3), iX_Crse_E1(3)
  INTEGER :: iX_Fine_B1(3), iX_Fine_E1(3)

  TYPE(MeshType) :: MeshX_Crse(3)
  TYPE(MeshType) :: MeshX_Fine(3)

  REAL(DP), ALLOCATABLE :: xL_Sub(:,:), xR_Sub(:,:)
  INTEGER :: iSub, iSubX1, iSubX2, iSubX3, nSub
  INTEGER :: iSubX(3), nSubX(3), nX_Sub(3)
  TYPE(MeshType), ALLOCATABLE :: MeshX_Sub(:,:)

  INTEGER :: iVar, nVar
  REAL(DP) :: U_T, U_A
  REAL(DP) :: AbsErr, RelErr, MaxError(nPR)

  CoordinateSystem = 'SPHERICAL'
  CoordinateSystem_lc = CoordinateSystem
  CALL string_lc( CoordinateSystem_lc )
  UseSimpleMeshRefinement = ( TRIM( CoordinateSystem ) == 'CARTESIAN' )

  ProgramName = 'MeshRefinementTest_TwoMoment'

  R0 = 4.0_DP * Kilometer

  nX = [ 16, 1, 1 ]
  swX = 0
  bcX = 0
  xL = [ 0.0_DP, 0.0_DP, 0.0_DP ]
  ZoomX = One

  SELECT CASE( TRIM( CoordinateSystem ) )
  CASE( 'SPHERICAL' )
    X1 = R0
    X2 = Pi
    X3 = TwoPi
  CASE( 'CYLINDRICAL' )
    X1 = R0 / SQRT( 2.0_DP )
    X2 = X1
    X3 = TwoPi
  CASE( 'CARTESIAN' )
    X1 = R0 / SQRT( 3.0_DP )
    X2 = X1
    X3 = X1
  CASE DEFAULT
    WRITE(*,*)
    WRITE(*,'(A5,A27,A)') &
      '', 'Invalid Coordinate System: ', TRIM( CoordinateSystem )
    STOP
  END SELECT
  xR = MERGE( [ X1,     X2,     X3 ], &
              [ X1, 1.0_DP, 1.0_DP ], &
              MASK = ( nX > 1 ) )

  nE = 2
  bcE = 0
  swE = 0
  eL = 0.0e0_DP
  eR = 1.0e2_DP
  ZoomE = One

  nNodes = 2
  nSpecies = 1
  nVar = nSpecies * nCR * nE * nNodes

  V_0 = [ 0.0_DP, 0.0_DP, 0.0_DP ]

  CALL InitializeDriver

  ! -- Initialize Meshes

  nFineX = 1
  nFineX(1:nDimsX) = refine_factor

  nFine = PRODUCT( nFineX )

  nX_Crse = 1
  nX_Crse(1:nDimsX) = nX(1:nDimsX)

  iX_Crse_B1 = 1
  iX_Crse_E1 = nX_Crse

  nX_Fine = 1
  nX_Fine(1:nDimsX) = nX_Crse(1:nDimsX) * nFineX(1:nDimsX)

  iX_Fine_B1 = 1
  iX_Fine_E1 = nX_Fine

  nSubX  = nFineX
  nSub   = PRODUCT( nSubX )
  nX_Sub = nX_Fine / nSubX

  ALLOCATE( xL_Sub(3,nSub) )
  ALLOCATE( xR_Sub(3,nSub) )
  ALLOCATE( MeshX_Sub(3,nSub) )

  ALLOCATE( U_0   (nDOFX,nVar,      nX_Crse(1),nX_Crse(2),nX_Crse(3) ) )
  ALLOCATE( U_Crse(nDOFX,nVar,      nX_Crse(1),nX_Crse(2),nX_Crse(3) ) )
  ALLOCATE( U_1   (nDOFX,nFine,nVar,nX_Crse(1),nX_Crse(2),nX_Crse(3) ) )
  ALLOCATE( U_Fine(nDOFX,nFine,nVar,nX_Crse(1),nX_Crse(2),nX_Crse(3) ) )

  PRINT*, "  SHAPE( U_Crse ) = ", SHAPE( U_Crse )
  PRINT*, "  SHAPE( U_Fine ) = ", SHAPE( U_Fine )
  PRINT*, "  nFine  = ", nFine
  PRINT*, "  nFineX = ", nFineX

  ALLOCATE( G_Crse(nDOFX,nX_Crse(1),nX_Crse(2),nX_Crse(3),nGF) )
  ALLOCATE( G_Fine(nDOFX,nX_Fine(1),nX_Fine(2),nX_Fine(3),nGF) )

  G_Crse = 0.0_DP
  G_Fine = 0.0_DP

  DO k = 1, 3
    CALL CreateMesh( MeshX_Crse(k), nX_Crse(k), nNodes, 0, xL(k), xR(k) )
    CALL CreateMesh( MeshX_Fine(k), nX_Fine(k), nNodes, 0, xL(k), xR(k) )
  END DO

  DO iSub = 1, nSub
    iSubX1 = MOD( (iSub-1)                          , nSubX(1) ) + 1
    iSubX2 = MOD( (iSub-1) / ( nSubX(1)            ), nSubX(2) ) + 1
    iSubX3 = MOD( (iSub-1) / ( nSubX(1) * nSubX(2) ), nSubX(3) ) + 1
    iSubX = [ iSubX1, iSubX2, iSubX3 ]
    DO k = 1, 3
      xL_Sub(k,iSub) = xL(k) + ( iSubX(k) - 1 ) * ( xR(k) - xL(k) ) / nSubX(k)
      xR_Sub(k,iSub) = xL(k) + ( iSubX(k)     ) * ( xR(k) - xL(k) ) / nSubX(k)
      CALL CreateMesh &
             ( MeshX_Sub(k,iSub), nX_Sub(k), &
               nNodes, 0, xL_Sub(k,iSub), xR_Sub(k,iSub) )
    END DO
  END DO

  ALLOCATE( X1_0(nX_Crse(1)*nNodesX(1)) )
  ALLOCATE( X2_0(nX_Crse(2)*nNodesX(2)) )
  ALLOCATE( X3_0(nX_Crse(3)*nNodesX(3)) )

#if   defined( THORNADO_OMP_OL )
  !$OMP TARGET ENTER DATA &
  !$OMP MAP( to:    nX_Crse, iX_Crse_B1, iX_Crse_E1, G_Crse, &
  !$OMP             nX_Fine, iX_Fine_B1, iX_Fine_E1, G_Fine, &
  !$OMP             nFineX, V_0 ) &
  !$OMP MAP( alloc: U_0, U_Crse, &
  !$OMP             U_1, U_Fine )
#elif defined( THORNADO_OACC   )
  !$ACC ENTER DATA &
  !$ACC COPYIN(     nX_Crse, iX_Crse_B1, iX_Crse_E1, G_Crse, &
  !$ACC             nX_Fine, iX_Fine_B1, iX_Fine_E1, G_Fine, &
  !$ACC             nFineX, V_0 ) &
  !$ACC CREATE(     U_0, U_Crse, &
  !$ACC             U_1, U_Fine )
#endif

  ! Calculate sqrt(Gamma) for geometry corrections
  call ComputeGeometryX( iX_Crse_B1, iX_Crse_E1, iX_Crse_B1, iX_Crse_E1, G_Crse, &
     MeshX_Option = MeshX_Crse, &
     CoordinateSystem_Option = TRIM( CoordinateSystem_lc ) )

  call ComputeGeometryX( iX_Fine_B1, iX_Fine_E1, iX_Fine_B1, iX_Fine_E1, G_Fine, &
     MeshX_Option = MeshX_Fine, &
     CoordinateSystem_Option = TRIM( CoordinateSystem_lc ) )

#if   defined( THORNADO_OMP_OL )
  !$OMP TARGET UPDATE &
  !$OMP FROM( G_Crse, G_Fine )
#elif defined( THORNADO_OACC   )
  !$ACC UPDATE &
  !$ACC HOST( G_Crse, G_Fine )
#endif

  ! --- Initialize Data on Coarse Level ---

  DO iS = 1, nSpecies
  DO iE = 1, nE
  DO iNodeE = 1, nDOFE

  DO iX3 = 1, nX_Crse(3)
  DO iX2 = 1, nX_Crse(2)
  DO iX1 = 1, nX_Crse(1)
  DO iNodeX = 1, nDOFX

    iNodeX1 = MOD( (iNodeX-1)                              , nNodesX(1) ) + 1
    iNodeX2 = MOD( (iNodeX-1) / ( nNodesX(1)              ), nNodesX(2) ) + 1
    iNodeX3 = MOD( (iNodeX-1) / ( nNodesX(1) * nNodesX(2) ), nNodesX(3) ) + 1

    X1 = NodeCoordinate( MeshX_Crse(1), iX1, iNodeX1 )
    X2 = NodeCoordinate( MeshX_Crse(2), iX2, iNodeX2 )
    X3 = NodeCoordinate( MeshX_Crse(3), iX3, iNodeX3 )

    iP_X1 = ( iX1 - 1 ) * nNodesX(1) + iNodeX1
    iP_X2 = ( iX2 - 1 ) * nNodesX(2) + iNodeX2
    iP_X3 = ( iX3 - 1 ) * nNodesX(3) + iNodeX3

    X1_0(iP_X1) = X1
    X2_0(iP_X2) = X2
    X3_0(iP_X3) = X3

    uPR = f_0( X1, X2, X3 )
    CALL ComputeConserved_TwoMoment &
           ( uPR(iPR_D ), uPR(iPR_I1), uPR(iPR_I2), uPR(iPR_I3), &
             uCR(iCR_N ), uCR(iCR_G1), uCR(iCR_G2), uCR(iCR_G3), &
             V_0(1), V_0(2), V_0(3), &
             G_Crse(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11),  &
             G_Crse(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22),  &
             G_Crse(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33) )

    DO iCR = 1, nCR

      iVar =   ( iS  - 1 ) * nDOFE * nE * nCR &
             + ( iCR - 1 ) * nDOFE * nE &
             + ( iE  - 1 ) * nDOFE &
             +   iNodeE

      U_0   (iNodeX,iVar,iX1,iX2,iX3) = uCR(iCR)
      U_Crse(iNodeX,iVar,iX1,iX2,iX3) = uCR(iCR)

    END DO

  END DO
  END DO
  END DO
  END DO

  END DO
  END DO
  END DO


  DO iS = 1, nSpecies
  DO iE = 1, nE
  DO iNodeE = 1, nDOFE

  DO iX3 = 1, nX_Crse(3)
  DO iX2 = 1, nX_Crse(2)
  DO iX1 = 1, nX_Crse(1)
  DO iFine = 1, nFine
  DO iNodeX = 1, nDOFX

    iNodeX1 = MOD( (iNodeX-1)                              , nNodesX(1) ) + 1
    iNodeX2 = MOD( (iNodeX-1) / ( nNodesX(1)              ), nNodesX(2) ) + 1
    iNodeX3 = MOD( (iNodeX-1) / ( nNodesX(1) * nNodesX(2) ), nNodesX(3) ) + 1

    iFineX1 = MOD( (iFine-1)                            , nFineX(1) ) + 1
    iFineX2 = MOD( (iFine-1) / ( nFineX(1)             ), nFineX(2) ) + 1
    iFineX3 = MOD( (iFine-1) / ( nFineX(1) * nFineX(2) ), nFineX(3) ) + 1

    iX1_Fine = ( iX1 - 1 ) * nFineX(1) + iFineX1
    iX2_Fine = ( iX2 - 1 ) * nFineX(2) + iFineX2
    iX3_Fine = ( iX3 - 1 ) * nFineX(3) + iFineX3

    X1 = NodeCoordinate( MeshX_Fine(1), iX1_Fine, iNodeX1 )
    X2 = NodeCoordinate( MeshX_Fine(2), iX2_Fine, iNodeX2 )
    X3 = NodeCoordinate( MeshX_Fine(3), iX3_Fine, iNodeX3 )

    uPR = f_0( X1, X2, X3 )
    CALL ComputeConserved_TwoMoment &
           ( uPR(iPR_D ), uPR(iPR_I1), uPR(iPR_I2), uPR(iPR_I3), &
             uCR(iCR_N ), uCR(iCR_G1), uCR(iCR_G2), uCR(iCR_G3), &
             V_0(1), V_0(2), V_0(3), &
             G_Fine(iNodeX,iX1_Fine,iX2_Fine,iX3_Fine,iGF_Gm_dd_11),  &
             G_Fine(iNodeX,iX1_Fine,iX2_Fine,iX3_Fine,iGF_Gm_dd_22),  &
             G_Fine(iNodeX,iX1_Fine,iX2_Fine,iX3_Fine,iGF_Gm_dd_33) )

    DO iCR = 1, nCR

      iVar =   ( iS  - 1 ) * nDOFE * nE * nCR &
             + ( iCR - 1 ) * nDOFE * nE &
             + ( iE  - 1 ) * nDOFE &
             +   iNodeE

      U_1(iNodeX,iFine,iVar,iX1,iX2,iX3) = uCR(iCR)

    END DO

  END DO
  END DO
  END DO
  END DO
  END DO

  END DO
  END DO
  END DO

#if   defined( THORNADO_OMP_OL )
  !$OMP TARGET UPDATE &
  !$OMP TO( U_0, U_1, U_Crse )
#elif defined( THORNADO_OACC   )
  !$ACC UPDATE &
  !$ACC DEVICE( U_0, U_1, U_Crse )
#endif

  CALL WriteRF_Crse( U_0 )
  CALL WriteVector( nDOFX*nVar*PRODUCT(nX_Crse),  &
                    RESHAPE( U_0(:,:,:,:,:), [nDOFX*nVar*PRODUCT(nX_Crse)] ), &
                    'U_Coarse_0.dat' )
  CALL WriteRF_Fine( U_1 )

  PRINT*, ""
  PRINT*, "Before Refinement: "
  PRINT*, "  MIN/MAX/SUM U_Crse = ", MINVAL( U_Crse ), MAXVAL( U_Crse ), SUM( U_Crse )

  ! --- Refine ---

  Timer_Refine = 0.0_DP
  CALL TimersStart( Timer_Refine )
  IF ( UseSimpleMeshRefinement ) THEN
    CALL RefineX_TwoMoment_SIMPLE( nX_Crse, nVar, U_Crse, U_Fine )
  ELSE
    CALL RefineX_TwoMoment_CURVILINEAR( nX_Crse, nX_Fine, nVar, G_Crse, U_Crse, G_Fine, U_Fine )
  END IF
  CALL TimersStop( Timer_Refine )

#if   defined( THORNADO_OMP_OL )
  !$OMP TARGET UPDATE &
  !$OMP FROM( U_Fine )
#elif defined( THORNADO_OACC   )
  !$ACC UPDATE &
  !$ACC HOST( U_Fine )
#endif

  PRINT*, ""
  PRINT*, "After Refinement: "
  PRINT*, "  MIN/MAX/SUM U_Fine = ", MINVAL( U_Fine ), MAXVAL( U_Fine ), SUM( U_Fine )

  MaxError = 0.0_DP
  DO iCR = 1, nCR

    DO iX3 = 1, nX_Crse(3)
    DO iX2 = 1, nX_Crse(2)
    DO iX1 = 1, nX_Crse(1)
    DO iFine = 1, nFine

    DO iS = 1, nSpecies
    DO iE = 1, nE
    DO iNodeE = 1, nDOFE

    DO iNodeX = 1, nDOFX

      iVar =   ( iS  - 1 ) * nDOFE * nE * nCR &
             + ( iCR - 1 ) * nDOFE * nE &
             + ( iE  - 1 ) * nDOFE &
             +   iNodeE

      U_T = U_Fine(iNodeX,iFine,iVar,iX1,iX2,iX3)
      U_A = U_1   (iNodeX,iFine,iVar,iX1,iX2,iX3)
      AbsErr = ABS( U_T - U_A )
      IF ( ABS( U_A ) > 0.0_DP ) THEN
        RelErr = AbsErr / ABS( U_A )
      ELSE IF ( ABS( U_T ) > 0.0_DP ) THEN
        RelErr = AbsErr / SqrtTiny
      ELSE
        RelErr = 0.0_DP
      END IF

      MaxError(iCR) = MAX( MaxError(iCR), RelErr )

    END DO

    END DO
    END DO
    END DO

    END DO
    END DO
    END DO
    END DO

  END DO

  WRITE(*,*)
  WRITE(*,'(A2,A)') '', 'INFO: Refine Relative Error'
  WRITE(*,*)
  WRITE(*,'(A4,4A12)') '', 'N', 'G1', 'G2', 'G3'
  WRITE(*,'(A4,4ES12.4E2)') '', MaxError

  CALL WriteVector( nX_Crse(1)*nNodesX(1), X1_0, 'X1_Coarse.dat' )
  CALL WriteVector( nX_Crse(2)*nNodesX(2), X2_0, 'X2_Coarse.dat' )
  CALL WriteVector( nX_Crse(3)*nNodesX(3), X3_0, 'X3_Coarse.dat' )

  CALL WriteVector( nX_Crse(1)*nNodesX(1), (xL(1)+X1_0)/Two, 'X1_Fine_01.dat' )
  CALL WriteVector( nX_Crse(1)*nNodesX(1), (xR(1)+X1_0)/Two, 'X1_Fine_02.dat' )
  CALL WriteVector( nX_Crse(2)*nNodesX(2), (xL(2)+X2_0)/Two, 'X2_Fine_01.dat' )
  CALL WriteVector( nX_Crse(2)*nNodesX(2), (XR(2)+X2_0)/Two, 'X2_Fine_02.dat' )
  CALL WriteVector( nX_Crse(3)*nNodesX(3), (xL(3)+X3_0)/Two, 'X3_Fine_01.dat' )
  CALL WriteVector( nX_Crse(3)*nNodesX(3), (XR(3)+X3_0)/Two, 'X3_Fine_02.dat' )

  CALL WriteRF_Fine( U_Fine )

  DO iFine = 1, nFine

    WRITE( MeshString, FMT='(i2.2)') iFine

    VectorName = 'U_Fine_' // MeshString // '.dat'

    CALL WriteVector( nDOFX*nVar*PRODUCT(nX_Crse),  &
                      RESHAPE( U_Fine(:,iFine,:,:,:,:), [nDOFX*nVar*PRODUCT(nX_Crse)] ), &
                      TRIM( VectorName ) )

  END DO

  ! --- Coarsen ---

#if   defined( THORNADO_OMP_OL )
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   )
  !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5)
#elif defined( THORNADO_OMP    )
  !$OMP PARALLEL DO COLLAPSE(5)
#endif
  DO iX3 = 1, nX_Crse(3)
  DO iX2 = 1, nX_Crse(2)
  DO iX1 = 1, nX_Crse(1)
  DO iVar = 1, nVar

    ! --- Reset U_Crse
    DO iNodeX = 1, nDOFX
      U_Crse(iNodeX,iVar,iX1,iX2,iX3) = Zero
    END DO

  END DO
  END DO
  END DO
  END DO

  PRINT*, ""
  PRINT*, "Before Coarsening: "
  PRINT*, "  MIN/MAX/SUM U_Fine = ", MINVAL( U_Fine ), MAXVAL( U_Fine ), SUM( U_Fine )

  Timer_Coarsen = 0.0_DP
  CALL TimersStart( Timer_Coarsen )
  IF ( UseSimpleMeshRefinement ) THEN
    CALL CoarsenX_TwoMoment_SIMPLE( nX_Crse, nVar, U_Fine, U_Crse )
  ELSE
    CALL CoarsenX_TwoMoment_CURVILINEAR( nX_Fine, nX_Crse, nVar, G_Fine, U_Fine, G_Crse, U_Crse )
  END IF
  CALL TimersStop( Timer_Coarsen )

  Timer_Total = Timer_Refine + Timer_Coarsen

#if   defined( THORNADO_OMP_OL )
  !$OMP TARGET UPDATE &
  !$OMP FROM( U_Crse )
#elif defined( THORNADO_OACC   )
  !$ACC UPDATE &
  !$ACC HOST( U_Crse )
#endif

  PRINT*, ""
  PRINT*, "After Coarsening: "
  PRINT*, "  MIN/MAX/SUM U_Crse = ", MINVAL( U_Crse ), MAXVAL( U_Crse ), SUM( U_Crse )

  MaxError = 0.0_DP
  DO iCR = 1, nCR

    DO iX3 = 1, nX_Crse(3)
    DO iX2 = 1, nX_Crse(2)
    DO iX1 = 1, nX_Crse(1)

    DO iS = 1, nSpecies
    DO iE = 1, nE
    DO iNodeE = 1, nDOFE

    DO iNodeX = 1, nDOFX

      iVar =   ( iS  - 1 ) * nDOFE * nE * nCR &
             + ( iCR - 1 ) * nDOFE * nE &
             + ( iE  - 1 ) * nDOFE &
             +   iNodeE

      U_T = U_Crse(iNodeX,iVar,iX1,iX2,iX3)
      U_A = U_0   (iNodeX,iVar,iX1,iX2,iX3)
      AbsErr = ABS( U_T - U_A )
      IF ( ABS( U_A ) > 0.0_DP ) THEN
        RelErr = AbsErr / ABS( U_A )
      ELSE IF ( ABS( U_T ) > 0.0_DP ) THEN
        RelErr = AbsErr / SqrtTiny
      ELSE
        RelErr = 0.0_DP
      END IF

      MaxError(iCR) = MAX( MaxError(iCR), RelErr )

    END DO

    END DO
    END DO
    END DO

    END DO
    END DO
    END DO

  END DO

  WRITE(*,*)
  WRITE(*,'(A2,A)') '', 'INFO: Coarsen Relative Error'
  WRITE(*,*)
  WRITE(*,'(A4,4A12)') '', 'N', 'G1', 'G2', 'G3'
  WRITE(*,'(A4,4ES12.4E2)') '', MaxError

  CALL WriteRF_Crse( U_Crse )
  CALL WriteVector( nDOFX*nVar*PRODUCT(nX_Crse),  &
                    RESHAPE( U_Crse(:,:,:,:,:), [nDOFX*nVar*PRODUCT(nX_Crse)] ), &
                    'U_Coarse_1.dat' )

  DO k = 1, 3
    CALL DestroyMesh( MeshX_Crse(k) )
    CALL DestroyMesh( MeshX_Fine(k) )
  END DO

  DO iFine = 1, nFine
  DO k = 1, 3
    CALL DestroyMesh( MeshX_Sub(k,iFine) )
  END DO
  END DO

#if   defined( THORNADO_OMP_OL )
  !$OMP TARGET EXIT DATA &
  !$OMP MAP( release: nX_Crse, iX_Crse_B1, iX_Crse_E1, G_Crse, &
  !$OMP               nX_Fine, iX_Fine_B1, iX_Fine_E1, G_Fine, &
  !$OMP               nFineX, V_0, &
  !$OMP               U_0, U_Crse, &
  !$OMP               U_1, U_Fine )
#elif defined( THORNADO_OACC   )
  !$ACC EXIT DATA &
  !$ACC DELETE(       nX_Crse, iX_Crse_B1, iX_Crse_E1, G_Crse, &
  !$ACC               nX_Fine, iX_Fine_B1, iX_Fine_E1, G_Fine, &
  !$ACC               nFineX, V_0, &
  !$ACC               U_0, U_Crse, &
  !$ACC               U_1, U_Fine )
#endif

  DEALLOCATE( X1_0, X2_0, X3_0 )
  DEALLOCATE( G_Crse, G_Fine, U_0, U_Crse, U_1, U_Fine )
  DEALLOCATE( xL_Sub, xR_Sub )
  DEALLOCATE( MeshX_Sub )

  WRITE(*,*)
  WRITE(*,*)
  WRITE(*,'(A4,A22,1ES10.2E2)') '', 'Refine = ',       &
    Timer_Refine
  WRITE(*,'(A4,A22,1ES10.2E2)') '', 'Coarsen = ',       &
    Timer_Coarsen
  WRITE(*,'(A4,A22,1ES10.2E2)') '', 'Total = ',       &
    Timer_Total

  CALL FinalizeDriver

CONTAINS


  FUNCTION f_0( X1, X2, X3 ) RESULT( uPR )

    REAL(DP), INTENT(in)  :: X1, X2, X3
    REAL(DP)              :: uPR(nPR)

    REAL(DP), PARAMETER :: Sigma = 1.0e2_DP
    REAL(DP), PARAMETER :: t_0   = 1.0e0_DP

    REAL(DP) :: D_0, I_0, R, Theta, Phi
    REAL(DP) :: C1, C2, C3

    SELECT CASE( TRIM( CoordinateSystem ) )
    CASE( 'SPHERICAL' )
      R = X1
      Theta = X2
      Phi = X3
      C1 = 1.0_DP
      C2 = 0.0_DP
      C3 = 0.0_DP
    CASE( 'CYLINDRICAL' )
      R = SQRT( X1*X1 + X2*X2 )
      IF ( X2 /= 0.0_DP ) THEN
        Theta = ATAN( X1 / X2 )
      ELSE
        Theta = 0.5 * Pi
      END IF
      Phi = X3
      C1 = SIN( Theta )
      C2 = COS( Theta )
      C3 = 0.0
    CASE( 'CARTESIAN' )
      R = SQRT( X1*X1 + X2*X2 + X3*X3 )
      Theta = ACOS( X3 / R )
      Phi = MOD( ATAN2( X2, X1 ), TwoPi )
      C1 = SIN( Theta ) * COS( Phi )
      C2 = SIN( Theta ) * SIN( Phi )
      C3 = COS( Theta )
    CASE DEFAULT
      WRITE(*,*)
      WRITE(*,'(A5,A27,A)') &
        '', 'Invalid Coordinate System: ', TRIM( CoordinateSystem )
    END SELECT

    !D_0 = ( Sigma / t_0 )**( 1.5_DP ) &
    !        * EXP( - 3.0_DP * Sigma * R**2 / ( 4.0_DP *t_0 ) )
    !I_0 = D_0 * R / ( 2.0_DP * t_0 )
    D_0 = SIN( TwoPi * ( R / R0 ) )
    I_0 = D_0

    uPR(iPR_D ) = D_0
    uPR(iPR_I1) = I_0 * C1
    uPR(iPR_I2) = I_0 * C2
    uPR(iPR_I3) = I_0 * C3

  END FUNCTION f_0


  SUBROUTINE InitializeDriver

    CALL InitializeProgram &
           ( ProgramName_Option &
               = TRIM( ProgramName ), &
             nX_Option &
               = nX, &
             swX_Option &
               = swX, &
             bcX_Option &
               = bcX, &
             xL_Option &
               = xL, &
             xR_Option &
               = xR, &
             ZoomX_Option &
               = ZoomX, &
             nE_Option &
               = nE, &
             swE_Option &
               = swE, &
             bcE_Option &
               = bcE, &
             eL_Option &
               = eL, &
             eR_Option &
               = eR, &
             ZoomE_Option &
               = ZoomE, &
             nNodes_Option &
               = nNodes, &
             CoordinateSystem_Option &
               = TRIM( CoordinateSystem ), &
             nSpecies_Option &
               = nSpecies, &
             BasicInitialization_Option &
               = .TRUE. )

    ! --- Position Space Reference Element ---

    CALL InitializeReferenceElementX

    CALL InitializeReferenceElementX_Lagrange

    ! --- Energy Space Reference Element ---

    CALL InitializeReferenceElementE

    CALL InitializeReferenceElementE_Lagrange

    ! --- Phase Space Reference Element ---

    CALL InitializeReferenceElement

    CALL InitializeReferenceElement_Lagrange

    ! --- Refine / Coarsen ---

    CALL InitializeMeshRefinement_TwoMoment &
           ( UseSimpleMeshRefinement_Option &
               = UseSimpleMeshRefinement, &
             Verbose_Option &
               = .TRUE. )

  END SUBROUTINE InitializeDriver


  SUBROUTINE FinalizeDriver

    CALL FinalizeMeshRefinement_TwoMoment

    ! --- Finalize ---

    CALL FinalizeReferenceElementX

    CALL FinalizeReferenceElementX_Lagrange

    CALL FinalizeReferenceElementE

    CALL FinalizeReferenceElementE_Lagrange

    CALL FinalizeReferenceElement

    CALL FinalizeReferenceElement_Lagrange

    CALL FinalizeProgram

  END SUBROUTINE FinalizeDriver


  SUBROUTINE string_lc(string)

    CHARACTER(*), INTENT(inout) :: string

    INTEGER, PARAMETER :: lc_a_ascii=IACHAR('a')
    INTEGER, PARAMETER :: uc_a_ascii=IACHAR('A')
    INTEGER, PARAMETER :: uc_z_ascii=IACHAR('Z')
    INTEGER :: i, x

    DO i = 1, LEN_TRIM(string)
      x = IACHAR(string(i:i))
      IF ( x >= uc_a_ascii .and. x <= uc_z_ascii ) THEN
        x = x + (lc_a_ascii - uc_a_ascii)
        string(i:i) = ACHAR(x)
      END IF
    END DO

  END SUBROUTINE string_lc


  SUBROUTINE WriteRF_Crse( U_Crse )

    Use MeshModule, Only : &
      MeshX
    USE RadiationFieldsModule, ONLY: &
      uPR, uCR
    USE InputOutputModuleHDF, ONLY: &
      WriteFieldsHDF

    REAL(DP), INTENT(in) :: U_Crse(:,:,:,:,:)

    TYPE(MeshType) :: MeshX_save(3)
    INTEGER :: iS, iCR, iX3, iX2, iX1, iE, iNodeX, iNodeE, iNode, iVar

    IF (      SIZE(U_Crse,3) /= SIZE(uCR,3) &
         .OR. SIZE(U_Crse,4) /= SIZE(uCR,4) &
         .OR. SIZE(U_Crse,5) /= SIZE(uCR,5) ) THEN
      WRITE(*,*)
      WRITE(*,'(A5,A27)') &
        '', 'WriteRF Shape Error'
      STOP
    END IF

    MeshX_save = MeshX
    MeshX = MeshX_crse

    DO iS = 1, nSpecies
    DO iCR = 1, nCR

    DO iX3 = 1, nX_Crse(3)
    DO iX2 = 1, nX_Crse(2)
    DO iX1 = 1, nX_Crse(1)

      DO iE = 1, nE
      DO iNodeX = 1, nDOFX
      DO iNodeE = 1, nDOFE

          iNode = ( iNodeX - 1 ) * nDOFE + iNodeE

          iVar =   ( iS  - 1 ) * nDOFE * nE * nCR &
                 + ( iCR - 1 ) * nDOFE * nE &
                 + ( iE  - 1 ) * nDOFE &
                 +   iNodeE

          uPR(iNode,iE,iX1,iX2,iX3,iCR,iS) = U_Crse(iNodeX,iVar,iX1,iX2,iX3)
          uCR(iNode,iE,iX1,iX2,iX3,iCR,iS) = U_Crse(iNodeX,iVar,iX1,iX2,iX3)

        END DO
        END DO
        END DO

    END DO
    END DO
    END DO

    END DO
    END DO

    CALL WriteFieldsHDF &
           ( Time = 0.0_DP, &
             WriteRF_Option = .TRUE. )

    MeshX = MeshX_save

  END SUBROUTINE WriteRF_Crse


  SUBROUTINE WriteRF_Fine( U_Fine )

    Use MeshModule, Only : &
      MeshX
    USE RadiationFieldsModule, ONLY: &
      uPR, uCR
    USE InputOutputModuleHDF, ONLY: &
      WriteFieldsHDF

    REAL(DP), INTENT(in) :: U_Fine(:,:,:,:,:,:)

    TYPE(MeshType) :: MeshX_save(3)
    INTEGER :: iSub, iSubX1, iSubX2, iSubX3
    INTEGER :: iX1_Sub, iX2_Sub, iX3_Sub
    INTEGER :: iS, iCR, iX3, iX2, iX1, iE, iNodeX, iNodeE, iNode, iVar

    IF (      SIZE(U_Fine,4) /= SIZE(uCR,3) &
         .OR. SIZE(U_Fine,5) /= SIZE(uCR,4) &
         .OR. SIZE(U_Fine,6) /= SIZE(uCR,5) ) THEN
      WRITE(*,*)
      WRITE(*,'(A5,A27)') &
        '', 'WriteRF Shape Error'
      STOP
    END IF

    MeshX_save = MeshX

    DO iSubX3 = 1, nSubX(3)
    DO iSubX2 = 1, nSubX(2)
    DO iSubX1 = 1, nSubX(1)

      iSub = ( iSubX3 -1 ) * nSubX(1) * nSubX(2) &
           + ( iSubX2 -1 ) * nSubX(1) &
           +   iSubX1

      MeshX = MeshX_Sub(:,iSub)

      DO iS = 1, nSpecies
      DO iCR = 1, nCR

      DO iX3_Sub = 1, nX_Sub(3)
      DO iX2_Sub = 1, nX_Sub(2)
      DO iX1_Sub = 1, nX_Sub(1)

        iX1_Fine = ( iSubX1 - 1 ) * nX_Sub(1) + iX1_Sub
        iX2_Fine = ( iSubX2 - 1 ) * nX_Sub(2) + iX2_Sub
        iX3_Fine = ( iSubX3 - 1 ) * nX_Sub(3) + iX3_Sub

        iFineX1 = MOD( ( iX1_Fine - 1 ), nFineX(1) ) + 1
        iFineX2 = MOD( ( iX2_Fine - 1 ), nFineX(2) ) + 1
        iFineX3 = MOD( ( iX3_Fine - 1 ), nFineX(3) ) + 1

        iX1 = MOD( ( iX1_Fine - 1 ) / nFineX(1), nX_Sub(1) ) + 1
        iX2 = MOD( ( iX2_Fine - 1 ) / nFineX(2), nX_Sub(2) ) + 1
        iX3 = MOD( ( iX3_Fine - 1 ) / nFineX(3), nX_Sub(3) ) + 1

        iFine = ( iFineX3 -1 ) * nFineX(1) * nFineX(2) &
              + ( iFineX2 -1 ) * nFineX(1) &
              +   iFineX1

        DO iE = 1, nE
        DO iNodeX = 1, nDOFX
        DO iNodeE = 1, nDOFE

          iVar =   ( iS  - 1 ) * nDOFE * nE * nCR &
                 + ( iCR - 1 ) * nDOFE * nE &
                 + ( iE  - 1 ) * nDOFE &
                 +   iNodeE

          iNode = ( iNodeX - 1 ) * nDOFE + iNodeE

          uPR(iNode,iE,iX1_Sub,iX2_Sub,iX3_Sub,iCR,iS) = U_Fine(iNodeX,iFine,iVar,iX1,iX2,iX3)
          uCR(iNode,iE,iX1_Sub,iX2_Sub,iX3_Sub,iCR,iS) = U_Fine(iNodeX,iFine,iVar,iX1,iX2,iX3)

        END DO
        END DO
        END DO

      END DO
      END DO
      END DO

      END DO
      END DO

      CALL WriteFieldsHDF &
             ( Time = 0.0_DP, &
               WriteRF_Option = .TRUE. )

    END DO
    END DO
    END DO

  END SUBROUTINE WriteRF_Fine

END PROGRAM MeshRefinementTest_TwoMoment
