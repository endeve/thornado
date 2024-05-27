MODULE InitializationModule_Relativistic

  USE KindModule, ONLY: &
    DP, &
    SqrtTiny, &
    Zero, &
    One, &
    Two, &
    Pi, &
    FourPi
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    nNodesX, &
    nDOFX, &
    swX, &
    iX_B0, &
    iX_B1, &
    iX_E0, &
    iX_E1
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX_FromScaleFactors
  USE GeometryFieldsModule, ONLY: &
    nGF, &
    uGF, &
    iGF_h_1, &
    iGF_h_2, &
    iGF_h_3, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_Beta_1, &
    iGF_Beta_2, &
    iGF_Beta_3, &
    iGF_Alpha, &
    iGF_Psi
  USE FluidFieldsModule, ONLY: &
    uPF, &
    iPF_D, &
    iPF_V1, &
    iPF_V2, &
    iPF_V3, &
    iPF_E, &
    iPF_Ne, &
    uCF, &
    iCF_D, &
    iCF_S1, &
    iCF_S2, &
    iCF_S3, &
    iCF_E, &
    iCF_Ne, &
    uAF, &
    iAF_P
  USE EquationOfStateModule_IDEAL, ONLY: &
    Gamma_IDEAL, &
    ComputePressureFromPrimitive_IDEAL
  USE Euler_UtilitiesModule_Relativistic, ONLY: &
    ComputeConserved_Euler_Relativistic
  USE UnitsModule, ONLY: &
    GravitationalConstant, &
    SpeedOfLight, &
    Kilometer, &
    SolarMass, &
    Gram, &
    Centimeter, &
    Erg, &
    Second, &
    PlanckConstant, &
    AtomicMassUnit
  USE UtilitiesModule, ONLY: &
    NodeNumberX, &
    Locate, &
    Interpolate1D_Linear

  USE HDF5

  IMPLICIT NONE
  PRIVATE

  INTEGER :: HDFERR

  PUBLIC :: InitializeFields_Relativistic


CONTAINS


  SUBROUTINE InitializeFields_Relativistic

    uPF(:,:,:,:,iPF_Ne) = Zero

    WRITE(*,*)
    WRITE(*,'(A,A)') '    INFO: ', TRIM( ProgramName )

    SELECT CASE ( TRIM( ProgramName ) )

      CASE( 'ShearingDisk_Unperturbed' )

         CALL InitializeFields_Unperturbed

      CASE( 'ShearingDisk_Perturbed' )

         !CALL InitializeFields_Perturbed

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,'(A21,A)') 'Invalid ProgramName: ', ProgramName
        WRITE(*,'(A)')     'Stopping...'
        STOP

    END SELECT

  END SUBROUTINE InitializeFields_Relativistic


  SUBROUTINE InitializeFields_Unperturbed

    CHARACTER(256) :: FileName

    INTEGER(HID_T) :: FILE_ID
    INTEGER        :: iX1, iX2, iX3, iNodeX, iNodeX1, iNodeX2, iNodeX3, &
                      jNodeX, jNodeX1, iL, nX, iGF
    REAL(DP) :: X1, X2, X3
    REAL(DP), ALLOCATABLE :: PressureArr(:), DensityArr(:), V3Arr(:), &
                             AlphaArr(:), PsiArr(:), X1Arr(:)

    FileName = "/home/jbuffal/thornado_MHD_3D/Workflow/MHD/ShearingDisk/GR_norot.h5"

    ! --- Populate arrays ---

    CALL H5OPEN_F( HDFERR )

    !PRINT*, FileName
    !PRINT*, TRIM( FileName )

    CALL H5FOPEN_F( TRIM( FileName ), H5F_ACC_RDONLY_F, FILE_ID, HDFERR )

    !CALL ReadDataset1DHDF( nX, '/size', FILE_ID )

    nX = 10000

    ALLOCATE( PressureArr(nX), DensityArr(nX), V3Arr(nX), AlphaArr(nX), &
              PsiArr(nX), X1Arr(nX) )

    CALL ReadDataset1DHDF( PsiArr,      '/psi',   FILE_ID )
    CALL ReadDataset1DHDF( AlphaArr,    '/alpha', FILE_ID )
    CALL ReadDataset1DHDF( X1Arr,       '/r',     FILE_ID )
    CALL ReadDataset1DHDF( PressureArr, '/pres',  FILE_ID )
    CALL ReadDataset1DHDF( DensityArr,  '/rho',   FILE_ID )
    CALL ReadDataset1DHDF( V3Arr,       '/V3',    FILE_ID )

    X1Arr       = X1Arr       * Centimeter
    DensityArr  = DensityArr  * ( Gram / Centimeter**3 )
    PressureArr = PressureArr * ( Erg  / Centimeter**3 )
    V3Arr       = V3Arr       * ( One  / Second )

    !PRINT*, "First: ", DensityArr(1) / ( Gram / Centimeter**3 )
    !PRINT*, "Last: ", DensityArr(nX) / ( Gram / Centimeter**3 )

    ! --- Map to 3D domain ---

    DO iX3 = iX_B1(3), iX_E1(3)
    DO iX2 = iX_B1(2), iX_E1(2)
    DO iX1 = iX_B1(1), iX_E1(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)
        iNodeX2 = NodeNumberTableX(2,iNodeX)
        iNodeX3 = NodeNumberTableX(3,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
        X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )
        X3 = NodeCoordinate( MeshX(3), iX3, iNodeX3 )

        ! --- Geometry Fields ---

        !PRINT*, 'Interpolating geometry fields.'

        uGF(iNodeX,iX1,iX2,iX3,iGF_Alpha) &
          = Interpolate1D( X1Arr, AlphaArr, SIZE( X1Arr ), X1 )

        uGF(iNodeX,iX1,iX2,iX3,iGF_Psi) &
          = Interpolate1D( X1Arr, PsiArr, SIZE( X1Arr ), X1 )

        uGF(iNodeX,iX1,iX2,iX3,iGF_h_1) &
          = uGF(iNodeX,iX1,iX2,iX3,iGF_Psi)**2
        uGF(iNodeX,iX1,iX2,iX3,iGF_h_2) &
          = uGF(iNodeX,iX1,iX2,iX3,iGF_Psi)**2
        uGF(iNodeX,iX1,iX2,iX3,iGF_h_3) &
          = uGF(iNodeX,iX1,iX2,iX3,iGF_Psi)**2 * X1

        CALL ComputeGeometryX_FromScaleFactors( uGF(:,iX1,iX2,iX3,:) )

        uGF(iNodeX,iX1,iX2,iX3,iGF_Beta_1) = Zero
        uGF(iNodeX,iX1,iX2,iX3,iGF_Beta_2) = Zero
        uGF(iNodeX,iX1,iX2,iX3,iGF_Beta_3) = Zero

        ! --- Fluid Fields ---

        !PRINT*, 'Interpolating fluid fields.'
        !PRINT*, 'Density.'

        uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
          = Interpolate1D( X1Arr, DensityArr, SIZE( X1Arr ), X1 )

        !PRINT*, 'Velocity.'
        uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = Zero
        uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = Zero
        uPF(iNodeX,iX1,iX2,iX3,iPF_V3) &
          = Interpolate1D( X1Arr, V3Arr, SIZE( X1Arr ), X1 )

        !PRINT*, 'Internal energy density.'
        uPF(iNodeX,iX1,iX2,iX3,iPF_E) &
          = Interpolate1D( X1Arr, PressureArr, SIZE( X1Arr ), X1 ) &
            / ( Gamma_IDEAL - One )

      END DO

      CALL ComputePressureFromPrimitive_IDEAL &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_E), &
               uPF(:,iX1,iX2,iX3,iPF_Ne), uAF(:,iX1,iX2,iX3,iAF_P) )

      CALL ComputeConserved_Euler_Relativistic &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_V1), &
               uPF(:,iX1,iX2,iX3,iPF_V2), uPF(:,iX1,iX2,iX3,iPF_V3), &
               uPF(:,iX1,iX2,iX3,iPF_E ), uPF(:,iX1,iX2,iX3,iPF_Ne), &
               uCF(:,iX1,iX2,iX3,iCF_D ), uCF(:,iX1,iX2,iX3,iCF_S1), &
               uCF(:,iX1,iX2,iX3,iCF_S2), uCF(:,iX1,iX2,iX3,iCF_S3), &
               uCF(:,iX1,iX2,iX3,iCF_E ), uCF(:,iX1,iX2,iX3,iCF_Ne), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33), &
               uAF(:,iX1,iX2,iX3,iAF_P) )

    END DO
    END DO
    END DO

    DEALLOCATE( X1Arr, PsiArr, AlphaArr, DensityArr, V3Arr, PressureArr )

  END SUBROUTINE InitializeFields_Unperturbed


  SUBROUTINE ReadDataset1DHDF( Dataset, DatasetName, FILE_ID )

    REAL(DP),         INTENT(out) :: Dataset(:)
    CHARACTER(LEN=*), INTENT(in)  :: DatasetName
    INTEGER(HID_T),   INTENT(in)  :: FILE_ID

    INTEGER(HID_T) :: DATASET_ID
    INTEGER(HID_T) :: DATASIZE(1)

    DATASIZE = SHAPE( Dataset )

    CALL H5DOPEN_F( FILE_ID, TRIM( DatasetName ), DATASET_ID, HDFERR )

    CALL H5DREAD_F( DATASET_ID, H5T_NATIVE_DOUBLE, Dataset, DATASIZE, HDFERR )

    CALL H5DCLOSE_F( DATASET_ID, HDFERR )

  END SUBROUTINE ReadDataset1DHDF


  REAL(DP) FUNCTION Interpolate1D( x, y, n, xq )

    INTEGER,                INTENT(in) :: n
    REAL(DP), DIMENSION(n), INTENT(in) :: x, y
    REAL(DP),               INTENT(in) :: xq

    INTEGER :: i

    i = Locate( xq, x, n )

    !PRINT*, 'i: ', i

    IF( i == 0 )THEN

      ! --- Extrapolate Left ---

      Interpolate1D &
        = Interpolate1D_Linear( xq, x(1), x(2), y(1), y(2) )

      !PRINT*, 'x(1): ', x(1)
      !PRINT*, 'x(2): ', x(2)
      !PRINT*, 'y(1): ', y(1)
      !PRINT*, 'y(2): ', y(2)

    ELSE IF( i == n )THEN

      ! --- Extrapolate Right ---

      Interpolate1D &
        = Interpolate1D_Linear( xq, x(n-1), x(n), y(n-1), y(n) )

      !PRINT*, 'x(n-1): ', x(n-1)
      !PRINT*, 'x(n): ',   x(n)
      !PRINT*, 'y(n-1): ', y(n-1)
      !PRINT*, 'y(n): ',   y(n)


    ELSE

      Interpolate1D &
        = Interpolate1D_Linear( xq, x(i), x(i+1), y(i), y(i+1) )

      !PRINT*, 'x(i): ', x(i)
      !PRINT*, 'x(i+1): ', x(i+1)
      !PRINT*, 'y(i): ', y(i)
      !PRINT*, 'y(i+1): ', y(i+1)

    END IF

    RETURN

  END FUNCTION Interpolate1D


END MODULE InitializationModule_Relativistic
