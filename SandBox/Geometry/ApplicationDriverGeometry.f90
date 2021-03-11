PROGRAM ApplicationDriverGeometry

  USE KindModule, ONLY: &
    DP, Zero, One, Pi
  USE ProgramHeaderModule, ONLY: &
    iX_B0, iX_E0, iX_B1, iX_E1, &
    iE_B0, iE_E0, iE_B1, iE_E1, &
    iZ_B0, iZ_E0, iZ_B1, iZ_E1, nDOFX, nNodesX
  USE GeometryFieldsModule, ONLY: &
    uGF, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33
  USE WeakDerivatives, ONLY: &
    ComputeWeakDerivatives_X0, &
    ComputeWeakDerivatives_X1, &
    ComputeWeakDerivatives_X2, &
    ComputeWeakDerivatives_X3, &
    ComputeChristoffel
  USE InputOutputModuleHDF, ONLY: &
    WriteDataset1DHDF, &
    WriteDataset3DHDF, &
    CreateGroupHDF
  USE InputOutputUtilitiesModule, ONLY: &
    NodeCoordinates, Field3D
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE ReferenceElementModule, ONLY: &
    NodeNumberTable
  USE MeshModule, ONLY: &
    MeshX, CreateMesh
  USE HDF5


  IMPLICIT NONE

  CHARACTER(32) :: ProgramName
  CHARACTER(32) :: CoordinateSystem
  INTEGER       :: nSpecies, iDim
  INTEGER       :: nNodes, iX1, iX2, iX3, iNodeX, iZ2, iZ3, iZ4
  INTEGER       :: nE, bcE, nX(3), bcX(3), swX(3)
  REAL(DP)      :: xL(3), xR(3), ZoomX(3) = One
  REAL(DP)      :: eL, eR, ZoomE = One
  REAL(DP)      :: t_end
    
  REAL(DP), ALLOCATABLE :: &
      dG_dd_dX0(:,:,:,:,:,:), & 
      dG_dd_dX1(:,:,:,:,:,:), &
      dG_dd_dX2(:,:,:,:,:,:), &
      dG_dd_dX3(:,:,:,:,:,:),  &
      Gamma_udd(:,:,:,:,:,:,:)


  ProgramName = "thornado"
  CoordinateSystem = 'SPHERICAL'
      nSpecies = 1
      nNodes   = 3

      nX  = [ 100, 100, 1 ]
      xL  = [ 0.0_DP, 0.0_DP, 0.0_DP ] 
      xR  = [ 1.0_DP, Pi, 1.0_DP ] 
      bcX = [ 0, 0, 0 ]
      swX = [ 1, 1, 0]

      nE    = 1
      eL    = 0.0_DP
      eR    = 1.0_DP
      bcE   = 0
      ZoomE = 1.0_DP




  CALL InitializeDriver

    CALL ComputeWeakDerivatives_X0 &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGF, dG_dd_dX0 )
    CALL ComputeWeakDerivatives_X1 &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGF, dG_dd_dX1 ) 


    CALL ComputeWeakDerivatives_X2 &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGF, dG_dd_dX2 ) 

    CALL ComputeWeakDerivatives_X3 &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGF, dG_dd_dX3 )

    CALL ComputeChristoffel( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGF, &
                             dG_dd_dX0, dG_dd_dX1, dG_dd_dX2, dG_dd_dX3, Gamma_udd )


    CALL WriteChristoffelHDF( Gamma_udd )

CONTAINS


  SUBROUTINE InitializeDriver

    USE ProgramInitializationModule, ONLY: &
      InitializeProgram
    USE ReferenceElementModule, ONLY: &
      InitializeReferenceElement
    USE ReferenceElementModuleX, ONLY: &
      InitializeReferenceElementX
    USE ReferenceElementModuleX_Lagrange, ONLY: &
      InitializeReferenceElementX_Lagrange
    USE GeometryComputationModule, ONLY: &
      ComputeGeometryX
    USE ReferenceElementModuleE, ONLY: &
      InitializeReferenceElementE
    USE ReferenceElementModuleE_Lagrange, ONLY: &
      InitializeReferenceElementE_Lagrange
    USE ReferenceElementModuleZ, ONLY: &
      InitializeReferenceElementZ
    USE ReferenceElementModule, ONLY: &
      InitializeReferenceElement
    USE ReferenceElementModule_Lagrange, ONLY: &
      InitializeReferenceElement_Lagrange


    CALL InitializeProgram &
           ( nX_Option &
               = nX, &
             swX_Option &
               = [ 1, 1, 1 ], &
             bcX_Option &
               = bcX, &
             xL_Option &
               = xL, &
             xR_Option &
               = xR, &
             nE_Option &
               = nE, &
             swE_Option &
               = 1, &
             bcE_Option &
               = bcE, &
             eL_Option &
               = eL, &
             eR_Option &
               = eR, &
             nNodes_Option &
               = nNodes, &
             CoordinateSystem_Option &
               = TRIM( CoordinateSystem ), &
             BasicInitialization_Option &
               = .TRUE. )
    ! --- Position Space Reference Element and Geometry ---

    CALL InitializeReferenceElement
    CALL InitializeReferenceElementX
    CALL InitializeReferenceElementX_Lagrange

    CALL ComputeGeometryX &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF )




    ALLOCATE( dG_dd_dX0(1:nDOFX,0:3,0:3, iX_B0(1):iX_E0(1), &
              iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3)  ) )
    ALLOCATE( dG_dd_dX1(1:nDOFX,0:3,0:3, iX_B0(1):iX_E0(1), &
              iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3)  ) )
    ALLOCATE( dG_dd_dX2(1:nDOFX,0:3,0:3, iX_B0(1):iX_E0(1), &
              iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3)  ) )
    ALLOCATE( dG_dd_dX3(1:nDOFX,0:3,0:3, iX_B0(1):iX_E0(1), &
              iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3)  ) )
    ALLOCATE( Gamma_udd(1:nDOFX,0:3,0:3,0:3, iX_B0(1):iX_E0(1), &
              iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3)  ) )



  END SUBROUTINE InitializeDriver

  SUBROUTINE WriteChristoffelHDF( Gamma_udd )
  
    REAL(DP), INTENT(in) :: & 
      Gamma_udd &
         (1:nDOFX,0:3,0:3,0:3, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))  


    CHARACTER(6)   :: FileNumberString
    CHARACTER(256) :: FileName
    CHARACTER(256) :: GroupName
    CHARACTER(256) :: DatasetName
    INTEGER        :: mu, nu, ru
    INTEGER(HID_T) :: FILE_ID
    REAL(DP)       :: Dummy3D(2,2,2) = 0.0_DP
    CHARACTER(9),  PARAMETER :: &
      OutputDirectory = '../Output'
    CHARACTER(14), PARAMETER :: &
      GeometrySuffix  = 'Christoffel'
    CHARACTER(3) :: ChristoffelNum
    CHARACTER(1) :: muc, nuc, ruc
    INTEGER :: FileNumber = 0
    INTEGER :: HDFERR




    WRITE( FileNumberString, FMT='(i6.6)') FileNumber

    FileName &
      = OutputDirectory // '/' // &
        TRIM( ProgramName ) // '_' // &
        TRIM(GeometrySuffix) // '_' // &
        FileNumberString // '.h5'


    CALL H5OPEN_F( HDFERR )

    CALL H5FCREATE_F( TRIM( FileName ), H5F_ACC_TRUNC_F, FILE_ID, HDFERR )

    GroupName = 'Spatial Grid'

    CALL CreateGroupHDF( FileName, TRIM( GroupName ) , FILE_ID )

    DatasetName = TRIM( GroupName ) // '/X1'

    CALL WriteDataset1DHDF &
           ( NodeCoordinates(MeshX(1),nX(1),nNodesX(1)),DatasetName,FILE_ID )

    DatasetName = TRIM( GroupName ) // '/X2'

    CALL WriteDataset1DHDF &
           ( NodeCoordinates(MeshX(2),nX(2),nNodesX(2)),DatasetName, FILE_ID )

    DatasetName = TRIM( GroupName ) // '/X3'

    CALL WriteDataset1DHDF &
           ( NodeCoordinates(MeshX(3),nX(3),nNodesX(3)),DatasetName, FILE_ID )


    GroupName = 'Geometry Fields'

    CALL CreateGroupHDF( FileName, TRIM( GroupName ), FILE_ID )


    DO mu = 0, 3
    DO nu = 0, 3
    DO ru = 0, 3


     IF (mu == 0) THEN
      muc= "0"
     ELSE IF(mu==1) THEN
      muc = "1"
     ELSE IF(mu==2) THEN
      muc = "2"
     ELSE IF(mu==3) THEN
      muc = "3"
    END IF
     IF (nu == 0) THEN
      nuc= "0"
     ELSE IF(nu==1) THEN
      nuc = "1"
     ELSE IF(nu==2) THEN
      nuc = "2"
     ELSE IF(nu==3) THEN
      nuc = "3"
    END IF
     IF (ru == 0) THEN
      ruc= "0"
     ELSE IF(ru==1) THEN
      ruc = "1"
     ELSE IF(ru==2) THEN
      ruc = "2"
     ELSE IF(ru==3) THEN
      ruc = "3"
    END IF
      ChristoffelNum = muc//nuc//ruc
      DatasetName = TRIM( GroupName ) // '/' // TRIM( ChristoffelNum )
      CALL WriteDataset3DHDF &
             ( Field3D &
                 ( Gamma_udd(1:nDOFX,mu,nu,ru,1:nX(1),1:nX(2),1:nX(3)), nX, nNodesX, &
                   nDOFX, NodeNumberTableX ), &
               DatasetName, FILE_ID )

    END DO
    END DO
    END DO

  END SUBROUTINE WriteChristoffelHDF

END PROGRAM ApplicationDriverGeometry
