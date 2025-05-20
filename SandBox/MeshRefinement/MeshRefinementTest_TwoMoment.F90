PROGRAM MeshRefinementTest_TwoMoment

  USE KindModule, ONLY: &
    DP, Zero, One, Two, TwoPi
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
     nGF, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33
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
    RefineX_TwoMoment, &
    CoarsenX_TwoMoment

  IMPLICIT NONE

  INTEGER, PARAMETER :: refine_factor = 2

  REAL(DP) :: &
    Timer_Refine, &
    Timer_Coarsen, &
    Timer_Total

  CHARACTER(32) :: ProgramName
  CHARACTER(32) :: CoordinateSystem
  CHARACTER(32) :: CoordinateSystem_lc

  INTEGER       :: nNodes
  INTEGER       :: nX(3), bcX(3), swX(3)
  REAL(DP)      :: xL(3), xR(3), ZoomX(3)
  INTEGER       :: nE, bcE, swE
  REAL(DP)      :: eL, eR, ZoomE
  INTEGER       :: nSpecies

  REAL(DP), ALLOCATABLE :: xL_Fine(:,:), xR_Fine(:,:)

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
  REAL(DP)      :: R, X1, X2, X3
  REAL(DP), ALLOCATABLE :: X1_0(:), X2_0(:), X3_0(:)
  REAL(DP)      :: V_0(3)
  REAL(DP)      :: uPR(nPR), uCR(nCR)

  REAL(DP), ALLOCATABLE :: U_0    (:,:,:,:,:)
  REAL(DP), ALLOCATABLE :: U_Crse (:,:,:,:,:)
  REAL(DP), ALLOCATABLE :: U_1    (:,:,:,:,:,:)
  REAL(DP), ALLOCATABLE :: U_Fine (:,:,:,:,:,:)
  REAL(DP), ALLOCATABLE :: U_Fine2(:,:,:,:,:,:)

  REAL(DP), ALLOCATABLE :: G_Crse(:,:,:,:,:)
  REAL(DP), ALLOCATABLE :: G_Fine(:,:,:,:,:)

  INTEGER :: nX_Crse(3)
  INTEGER :: nX_Fine(3)

  INTEGER :: iX_Crse_B1(3), iX_Crse_E1(3)
  INTEGER :: iX_Fine_B1(3), iX_Fine_E1(3)

  TYPE(MeshType) :: MeshX_Crse(3)
  TYPE(MeshType) :: MeshX_Fine(3)
  TYPE(MeshType), ALLOCATABLE :: SubMeshX_Fine(:,:)

  INTEGER :: iVar, nVar
  REAL(DP) :: MaxError(nPR)

  CoordinateSystem = 'SPHERICAL'
  CoordinateSystem_lc = CoordinateSystem
  CALL string_lc( CoordinateSystem_lc )

  ProgramName = 'MeshRefinementTest_TwoMoment'

  nX = [ 16, 1, 1 ]
  swX = 0
  bcX = 0
  xL = 0.0_DP
  xR = 1.0_DP
  ZoomX = One

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

  ALLOCATE( xL_Fine(3,nFine) )
  ALLOCATE( xR_Fine(3,nFine) )
  ALLOCATE( SubMeshX_Fine(3,nFine) )

  ALLOCATE( U_0    (nDOFX,      nX_Crse(1),nX_Crse(2),nX_Crse(3),nVar ) )
  ALLOCATE( U_Crse (nDOFX,      nX_Crse(1),nX_Crse(2),nX_Crse(3),nVar ) )
  ALLOCATE( U_1    (nDOFX,nFine,nX_Crse(1),nX_Crse(2),nX_Crse(3),nVar ) )
  ALLOCATE( U_Fine (nDOFX,nFine,nX_Crse(1),nX_Crse(2),nX_Crse(3),nVar ) )
  ALLOCATE( U_Fine2(nDOFX,nFine,nX_Crse(1),nX_Crse(2),nX_Crse(3),nVar) )

  PRINT*, "  SHAPE( U_Crse ) = ", SHAPE( U_Crse )
  PRINT*, "  SHAPE( U_Fine ) = ", SHAPE( U_Fine )
  PRINT*, "  nFine  = ", nFine
  PRINT*, "  nFineX = ", nFineX

  ALLOCATE( G_Crse(nDOFX,nX_Crse(1),nX_Crse(2),nX_Crse(3),nGF) )
  ALLOCATE( G_Fine(nDOFX,nX_Fine(1),nX_Fine(2),nX_Fine(3),nGF) )

  G_Crse = 0.0
  G_Fine = 0.0

  DO k = 1, 3
    CALL CreateMesh( MeshX_Crse(k), nX_Crse(k), nNodes, 0, xL(k), xR(k) )
    CALL CreateMesh( MeshX_Fine(k), nX_Fine(k), nNodes, 0, xL(k), xR(k) )
  END DO

  DO iFineX3 = 1, nFineX(3)
  DO iFineX2 = 1, nFineX(2)
  DO iFineX1 = 1, nFineX(1)

    iFine = (iFineX3-1)*nFineX(1)*nFineX(2) &
          + (iFineX2-1)*nFineX(1) &
          +  iFineX1

    xL_Fine(1,iFine) = xL(1) + ( iFineX1 - 1 ) * ( xR(1) - xL(1) ) / nFineX(1)
    xR_Fine(1,iFine) = xL(1) + ( iFineX1     ) * ( xR(1) - xL(1) ) / nFineX(1)
    xL_Fine(2,iFine) = xL(2) + ( iFineX2 - 1 ) * ( xR(2) - xL(2) ) / nFineX(2)
    xR_Fine(2,iFine) = xL(2) + ( iFineX2     ) * ( xR(2) - xL(2) ) / nFineX(2)
    xL_Fine(3,iFine) = xL(3) + ( iFineX3 - 1 ) * ( xR(3) - xL(3) ) / nFineX(3)
    xR_Fine(3,iFine) = xL(3) + ( iFineX3     ) * ( xR(3) - xL(3) ) / nFineX(3)
    DO k = 1, 3
      CALL CreateMesh &
             ( SubMeshX_Fine(k,iFine), nX_Fine(k) / nFineX(k), &
               nNodes, 0, xL_Fine(k,iFine), xR_Fine(k,iFine) )
    END DO

  END DO
  END DO
  END DO

  ALLOCATE( X1_0(nX_Crse(1)*nNodesX(1)) )
  ALLOCATE( X2_0(nX_Crse(2)*nNodesX(2)) )
  ALLOCATE( X3_0(nX_Crse(3)*nNodesX(3)) )

  ! Calculate sqrt(Gamma) for geometry corrections
  call ComputeGeometryX( iX_Crse_B1, iX_Crse_E1, iX_Crse_B1, iX_Crse_E1, G_Crse, &
     MeshX_Option = MeshX_Crse, &
     CoordinateSystem_Option = TRIM( CoordinateSystem_lc ) )

  call ComputeGeometryX( iX_Fine_B1, iX_Fine_E1, iX_Fine_B1, iX_Fine_E1, G_Fine, &
     MeshX_Option = MeshX_Fine, &
     CoordinateSystem_Option = TRIM( CoordinateSystem_lc ) )

  ! --- Initialize Data on Coarse Level ---

  DO iX3 = 1, nX_Crse(3)
  DO iX2 = 1, nX_Crse(2)
  DO iX1 = 1, nX_Crse(1)
  DO iNodeX = 1, nDOFX

  DO iS = 1, nSpecies
  DO iE = 1, nE
  DO iNodeE = 1, nDOFE

      iNodeX1 = NodeNumberTableX(1,iNodeX)
      iNodeX2 = NodeNumberTableX(2,iNodeX)
      iNodeX3 = NodeNumberTableX(3,iNodeX)

      X1 = NodeCoordinate( MeshX_Crse(1), iX1, iNodeX1 )
      X2 = NodeCoordinate( MeshX_Crse(2), iX2, iNodeX2 )
      X3 = NodeCoordinate( MeshX_Crse(3), iX3, iNodeX3 )

      SELECT CASE( TRIM( CoordinateSystem ) )
      CASE( 'SPHERICAL' )
        R = X1
      CASE( 'CYLINDRICAL' )
        R = SQRT( X1*X1 + X2*X2 )
      CASE( 'CARTESIAN' )
        R = SQRT( X1*X1 + X2*X2 + X3*X3 )
      CASE DEFAULT
        WRITE(*,*)
        WRITE(*,'(A5,A27,A)') &
          '', 'Invalid Coordinate System: ', TRIM( CoordinateSystem )
        STOP
      END SELECT

      iP_X1 = ( iX1 - 1 ) * nNodesX(1) + iNodeX1
      iP_X2 = ( iX2 - 1 ) * nNodesX(2) + iNodeX2
      iP_X3 = ( iX3 - 1 ) * nNodesX(3) + iNodeX3

      X1_0(iP_X1) = X1
      X2_0(iP_X2) = X2
      X3_0(iP_X3) = X3

      uPR = f_0( R )
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

        U_0   (iNodeX,iX1,iX2,iX3,iVar) = uCR(iCR)
        U_Crse(iNodeX,iX1,iX2,iX3,iVar) = uCR(iCR)

      END DO

      DO iFine = 1, nFine

        iFineX1 = MOD( (iFine-1)                            , nFineX(1) ) + 1
        iFineX2 = MOD( (iFine-1) / ( nFineX(1)             ), nFineX(2) ) + 1
        iFineX3 = MOD( (iFine-1) / ( nFineX(1) * nFineX(2) ), nFineX(3) ) + 1

        iX1_Fine = ( iX1 - 1 ) * nFineX(1) + iFineX1
        iX2_Fine = ( iX2 - 1 ) * nFineX(2) + iFineX2
        iX3_Fine = ( iX3 - 1 ) * nFineX(3) + iFineX3

        X1 = NodeCoordinate( MeshX_Fine(1), iX1_Fine, iNodeX1 )
        X2 = NodeCoordinate( MeshX_Fine(2), iX2_Fine, iNodeX2 )
        X3 = NodeCoordinate( MeshX_Fine(3), iX3_Fine, iNodeX3 )

        SELECT CASE( TRIM( CoordinateSystem ) )
        CASE( 'SPHERICAL' )
          R = X1
        CASE( 'CYLINDRICAL' )
          R = SQRT( X1*X1 + X2*X2 )
        CASE( 'CARTESIAN' )
          R = SQRT( X1*X1 + X2*X2 + X3*X3 )
        CASE DEFAULT
          WRITE(*,*)
          WRITE(*,'(A5,A27,A)') &
            '', 'Invalid Coordinate System: ', TRIM( CoordinateSystem )
          STOP
        END SELECT

        uPR = f_0( R )
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

          U_1(iNodeX,iFine,iX1,iX2,iX3,iVar) = uCR(iCR)

        END DO

      END DO

    END DO
    END DO
    END DO

  END DO
  END DO
  END DO
  END DO

  ! --- Refine ---

  Timer_Refine = 0.0_DP
  CALL TimersStart( Timer_Refine )
  CALL RefineX_TwoMoment( nX_Crse, nVar, U_Crse, U_Fine )
  CALL TimersStop( Timer_Refine )

  PRINT*, ""
  PRINT*, "After Refinement: "
  PRINT*, "  MIN/MAX/SUM U_Crse = ", MINVAL( U_Crse ), MAXVAL( U_Crse ), SUM( U_Crse )
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

      MaxError(iCR) = MAX( MaxError(iCR), &
                           ABS(   U_Fine(iNodeX,iFine,iX1,iX2,iX3,iVar) &
                                - U_1   (iNodeX,iFine,iX1,iX2,iX3,iVar) ) )

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
  WRITE(*,'(A2,A)') '', 'INFO: Refine Error'
  WRITE(*,*)
  WRITE(*,'(A4,4A12)') '', 'N', 'G1', 'G2', 'G3'
  WRITE(*,'(A4,4ES12.4E2)') '', MaxError

  CALL WriteVector( nX_Crse(1)*nNodesX(1), X1_0, 'X1_Coarse.dat' )
  CALL WriteVector( nX_Crse(2)*nNodesX(2), X2_0, 'X2_Coarse.dat' )
  CALL WriteVector( nX_Crse(3)*nNodesX(3), X3_0, 'X3_Coarse.dat' )

  CALL WriteRF( MeshX_Crse, U_Crse )
  CALL WriteVector( nDOFX*nVar*PRODUCT(nX_Crse),  &
                    RESHAPE( U_Crse(:,:,:,:,:), [nDOFX*nVar*PRODUCT(nX_Crse)] ), &
                    'U_Coarse_0.dat' )

  CALL WriteVector( nX_Crse(1)*nNodesX(1), (xL(1)+X1_0)/Two, 'X1_Fine_01.dat' )
  CALL WriteVector( nX_Crse(1)*nNodesX(1), (xR(1)+X1_0)/Two, 'X1_Fine_02.dat' )
  CALL WriteVector( nX_Crse(2)*nNodesX(2), (xL(2)+X2_0)/Two, 'X2_Fine_01.dat' )
  CALL WriteVector( nX_Crse(2)*nNodesX(2), (XR(2)+X2_0)/Two, 'X2_Fine_02.dat' )
  CALL WriteVector( nX_Crse(3)*nNodesX(3), (xL(3)+X3_0)/Two, 'X3_Fine_01.dat' )
  CALL WriteVector( nX_Crse(3)*nNodesX(3), (XR(3)+X3_0)/Two, 'X3_Fine_02.dat' )

  DO iFine = 1, nFine

    WRITE( MeshString, FMT='(i2.2)') iFine

    VectorName = 'U_Fine_' // MeshString // '.dat'

    CALL WriteRF( SubMeshX_Fine(:,iFine), U_Fine(:,iFine,:,:,:,:) )
    CALL WriteVector( nDOFX*nVar*PRODUCT(nX_Crse),  &
                      RESHAPE( U_Fine(:,iFine,:,:,:,:), [nDOFX*nVar*PRODUCT(nX_Crse)] ), &
                      TRIM( VectorName ) )

  END DO

  ! --- Coarsen ---

  DO iX3 = 1, nX_Crse(3)
  DO iX2 = 1, nX_Crse(2)
  DO iX1 = 1, nX_Crse(1)
  DO iVar = 1, nVar

    ! --- Coarsen operation needs different ordering
    DO iFine = 1, nFine
    DO iNodeX = 1, nDOFX
      U_Fine2(iNodeX,iFine,iX1,iX2,iX3,iVar) = U_Fine(iNodeX,iFine,iX1,iX2,iX3,iVar)
    END DO
    END DO

    ! --- Reset U_Crse
    DO iNodeX = 1, nDOFX
      U_Crse(iNodeX,iX1,iX2,iX3,iVar) = Zero
    END DO

  END DO
  END DO
  END DO
  END DO

  Timer_Coarsen = 0.0_DP
  CALL TimersStart( Timer_Coarsen )
  CALL CoarsenX_TwoMoment( nX_Crse, nVar, U_Fine2, U_Crse )
  CALL TimersStop( Timer_Coarsen )

  Timer_Total = Timer_Refine + Timer_Coarsen

  PRINT*, ""
  PRINT*, "After Coarsening: "
  PRINT*, "  MIN/MAX/SUM U_Crse = ", MINVAL( U_Crse ), MAXVAL( U_Crse ), SUM( U_Crse )
  PRINT*, "  MIN/MAX/SUM U_Fine = ", MINVAL( U_Fine2), MAXVAL( U_Fine2), SUM( U_Fine2)

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

      MaxError(iCR) = MAX( MaxError(iCR), &
                           ABS(   U_Crse(iNodeX,iX1,iX2,iX3,iVar) &
                                - U_0   (iNodeX,iX1,iX2,iX3,iVar) ) )

    END DO

    END DO
    END DO
    END DO

    END DO
    END DO
    END DO

  END DO

  WRITE(*,*)
  WRITE(*,'(A2,A)') '', 'INFO: Coarsen Error'
  WRITE(*,*)
  WRITE(*,'(A4,4A12)') '', 'N', 'G1', 'G2', 'G3'
  WRITE(*,'(A4,4ES12.4E2)') '', MaxError

  CALL WriteRF( MeshX_Crse, U_Crse )
  CALL WriteVector( nDOFX*nVar*PRODUCT(nX_Crse),  &
                    RESHAPE( U_Crse(:,:,:,:,:), [nDOFX*nVar*PRODUCT(nX_Crse)] ), &
                    'U_Coarse_1.dat' )

  DO k = 1, 3
    CALL DestroyMesh( MeshX_Crse(k) )
    CALL DestroyMesh( MeshX_Fine(k) )
  END DO

  DO iFine = 1, nFine
  DO k = 1, 3
    CALL DestroyMesh( SubMeshX_Fine(k,iFine) )
  END DO
  END DO

  DEALLOCATE( X1_0, X2_0, X3_0 )
  DEALLOCATE( G_Crse, G_Fine, U_0, U_Crse, U_1, U_Fine, U_Fine2 )
  DEALLOCATE( SubMeshX_Fine )

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

  FUNCTION f_0( R ) RESULT( uPR )

    REAL(DP), INTENT(in)  :: R
    REAL(DP)              :: uPR(nPR)

    REAL(DP), PARAMETER :: Sigma = 1.0e2_DP
    REAL(DP), PARAMETER :: t_0   = 1.0e0_DP

    REAL(DP) :: D_0, I_0

    D_0 = ( Sigma / t_0 )**( 1.5_DP ) &
            * EXP( - 3.0_DP * Sigma * R**2 / ( 4.0_DP *t_0 ) )

    I_0 = D_0 * R / ( 2.0_DP * t_0 )

    uPR(iPR_D ) = D_0
    uPR(iPR_I1) = I_0
    uPR(iPR_I2) = 0.0_DP
    uPR(iPR_I3) = 0.0_DP

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

    CALL InitializeMeshRefinement_TwoMoment

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


  SUBROUTINE WriteRF( MeshX_write, U_write )

    Use MeshModule, Only : &
      MeshX
    USE RadiationFieldsModule, ONLY: &
      uPR, uCR
    USE InputOutputModuleHDF, ONLY: &
      WriteFieldsHDF

    TYPE(MeshType), INTENT(in) :: MeshX_write(3)
    REAL(DP), INTENT(in) :: U_write(:,:,:,:,:)

    TYPE(MeshType) :: MeshX_save(3)
    INTEGER :: iS, iCR, iX3, iX2, iX1, iE, iNodeX, iNodeE, iNode, iVar

    MeshX_save = MeshX
    MeshX = MeshX_write

    IF (      SIZE(U_write,2) /= SIZE(uCR,3) &
         .OR. SIZE(U_write,3) /= SIZE(uCR,4) &
         .OR. SIZE(U_write,4) /= SIZE(uCR,5) ) THEN
      WRITE(*,*)
      WRITE(*,'(A5,A27)') &
        '', 'WriteRF Shape Error'
      STOP
    END IF

    DO iS = 1, nSpecies
    DO iCR = 1, nCR

    DO iX3 = 1, SIZE(uCR,5)
    DO iX2 = 1, SIZE(uCR,4)
    DO iX1 = 1, SIZE(uCR,3)
    DO iE = 1, nE

      DO iNodeX = 1, nDOFX
      DO iNodeE = 1, nDOFE

          iNode = ( iNodeX - 1 ) * nDOFE + iNodeE

          iVar =   ( iS  - 1 ) * nDOFE * nE * nCR &
                 + ( iCR - 1 ) * nDOFE * nE &
                 + ( iE  - 1 ) * nDOFE &
                 +   iNodeE

          uPR(iNode,iE,iX1,iX2,iX3,iCR,iS) = U_write(iNodeX,iX1,iX2,iX3,iVar)
          uCR(iNode,iE,iX1,iX2,iX3,iCR,iS) = U_write(iNodeX,iX1,iX2,iX3,iVar)

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

  END SUBROUTINE WriteRF

END PROGRAM MeshRefinementTest_TwoMoment
