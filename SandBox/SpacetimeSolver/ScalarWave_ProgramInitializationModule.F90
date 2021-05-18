MODULE ScalarWave_ProgramInitializationModule

  USE KindModule, ONLY: &
    DP
  USE UnitsModule, ONLY: &
    ActivateUnitsDisplay, &
    DescribeUnitsDisplay, &
    UnitsDisplay
  USE ProgramHeaderModule, ONLY: &
    nX, swX, xL, xR, ZoomX, &
    nE, swE, eL, eR, ZoomE, &
    nDimsX, nDimsE, nDims, &
    nDOFX,  nDOFE,  nDOF, &
    nNodesX, nNodesE, nNodes, &
    InitializeProgramHeader
  USE UtilitiesModule, ONLY: &
    InitializeWeights
  USE QuadratureModule, ONLY: &
    InitializeQuadratures
  USE PolynomialBasisModuleX_Lagrange, ONLY: &
    InitializePolynomialBasisX_Lagrange
  USE PolynomialBasisModuleX_Legendre, ONLY: &
    InitializePolynomialBasisX_Legendre
  USE PolynomialBasisMappingModule, ONLY: &
    InitializePolynomialBasisMapping
  USE MeshModule, ONLY: &
    MeshX, MeshE, &
    CreateMesh, &
    DestroyMesh
  USE ScalarFieldsModule, ONLY: &
    CreateScalarFields, &
    DestroyScalarFields
  USE DeviceModule, ONLY: &
    InitializeDevice, &
    FinalizeDevice

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  INTEGER :: mpierr

  PUBLIC :: InitializeProgram
  PUBLIC :: FinalizeProgram

CONTAINS


  SUBROUTINE InitializeProgram &
    ( ProgramName_Option, nX_Option, swX_Option, bcX_Option, xL_Option, &
      xR_Option, zoomX_Option, CoordinateSystem_Option, &
      nNodes_Option, ActivateUnits_Option )

    CHARACTER(LEN=*),       INTENT(in), OPTIONAL :: ProgramName_Option
    INTEGER,  DIMENSION(3), INTENT(in), OPTIONAL :: nX_Option
    INTEGER,  DIMENSION(3), INTENT(in), OPTIONAL :: swX_Option
    INTEGER,  DIMENSION(3), INTENT(in), OPTIONAL :: bcX_Option
    REAL(DP), DIMENSION(3), INTENT(in), OPTIONAL :: xL_Option
    REAL(DP), DIMENSION(3), INTENT(in), OPTIONAL :: xR_Option
    REAL(DP), DIMENSION(3), INTENT(in), OPTIONAL :: zoomX_Option
    INTEGER,                INTENT(in), OPTIONAL :: nNodes_Option
    CHARACTER(LEN=*),       INTENT(in), OPTIONAL :: CoordinateSystem_Option
    LOGICAL,                INTENT(in), OPTIONAL :: ActivateUnits_Option

    LOGICAL :: ActivateUnits
    INTEGER :: iDim

    CALL MPI_INIT( mpierr )

    ! --- Device ---

    CALL InitializeDevice

    CALL InitializeProgramHeader     &
           ( ProgramName_Option      &
               = ProgramName_Option, &
             nNodes_Option           &
               = nNodes_Option,      &
             nX_Option               &
               = nX_Option,          &
             swX_Option              &
               = swX_Option,         &
             bcX_Option              &
               = bcX_Option,         &
             xL_Option               &
               = xL_Option,          &
             xR_Option               &
               = xR_Option,          &
             zoomX_Option            &
               = zoomX_Option )

    ! --- Units ---

    IF( PRESENT( ActivateUnits_Option ) )THEN
      ActivateUnits = ActivateUnits_Option
    ELSE
      ActivateUnits = .FALSE.
    END IF

    IF( ActivateUnits )THEN

      CALL ActivateUnitsDisplay &
             ( CoordinateSystem_Option = CoordinateSystem_Option )
      CALL DescribeUnitsDisplay

    END IF

    WRITE(*,'(A5,A17,I1)') &
      '', 'Dimensionality = ', nDims
    WRITE(*,'(A5,A18)') &
      '', '------------------'
    WRITE(*,*)
    WRITE(*,'(A7,A9,I1,A2,A9,I1)') &
      '', 'nDimsX = ', nDimsX, ''
    WRITE(*,*)

    WRITE(*,'(A7,A9,I4.4,A2,A9,I4.4,A2,A9,I4.4,A2)') &
      '', 'nX(1) = ', nX(1), '', 'nX(2) = ', nX(2), &
      '', 'nX(3) = ', nX(3)
    WRITE(*,'(A7,A9,I4,A2,A9,I4,A2,A9,I4,A2)') &
      '', 'swX(1) = ', swX(1), '', 'swX(2) = ', swX(2), &
      '', 'swX(3) = ', swX(3)
    WRITE(*,*)

    WRITE(*,'(A5,A36)') '', 'Degrees of Freedom / Element / Field'
    WRITE(*,*)
    WRITE(*,'(A7,A9,I2.2)') &
      '', 'nNodes = ', nNodes
    WRITE(*,*)
    DO iDim = 1, 3
      WRITE(*,'(A9,A4,I1,A10,I1,A4,I2.2)') &
        '', 'i = ', iDim, ', nNodesX(', iDim, ') = ', nNodesX(iDim)
    END DO
    WRITE(*,*)
    WRITE(*,'(A7,A8,I4.4,A2,A8,I4.4,A2,A7,I4.4)') &
      '', 'nDOFX = ', nDOFX
    WRITE(*,*)

    ! --- Quadratures ---

    CALL InitializeQuadratures

    ! --- Polynomial Basis ---

    CALL InitializePolynomialBasisX_Lagrange

    CALL InitializePolynomialBasisX_Legendre


    ASSOCIATE( U => UnitsDisplay )

    ! --- Spatial Grid ---

    WRITE(*,'(A5,A)') '', 'Computational Domain'
    WRITE(*,'(A5,A)') '', '--------------------'
    WRITE(*,*)
    WRITE(*,'(A7,A)') '', 'Spatial Domain:'
    WRITE(*,'(A7,A)') '', '---------------'
    WRITE(*,*)

    iDim = 1
    WRITE(*,'(A7,A3,I1,A4,ES9.2E2,A1,A,A2,A3,I1,A4,&
              &ES9.2E2,A1,A,A2,A6,I1,A4,ES10.4E2)') &
      '',   'xL(', iDim, ') = ', xL(iDim) / U % LengthX1Unit, &
      '', TRIM( U % LengthX1Label ), &
      ', ', 'xR(', iDim, ') = ', xR(iDim) / U % LengthX1Unit, &
      '', TRIM( U % LengthX1Label ), &
      ', ', 'ZoomX(', iDim, ') = ', ZoomX(iDim)

    WRITE(*,*)

    CALL CreateMesh &
           ( MeshX(iDim), nX(iDim), nNodesX(iDim), swX(iDim), &
             xL(iDim), xR(iDim), ZoomOption = ZoomX(iDim) )

    WRITE(*,'(A9,A11,I1,A4,ES9.2E2,A1,A,A3,ES9.2E2,A1,A)') &
      '', 'MIN/MAX dx(', iDim, ') = ', &
      MINVAL( MeshX(iDim) % Width(1:nX(iDim)) ) &
        / U % LengthX1Unit, '', TRIM( U % LengthX1Label ), &
      ' / ', &
      MAXVAL( MeshX(iDim) % Width(1:nX(iDim)) ) &
        / U % LengthX1Unit, '', TRIM( U % LengthX1Label )
    WRITE(*,*)

    iDim = 2
    WRITE(*,'(A7,A3,I1,A4,ES9.2E2,A1,A,A2,A3,I1,A4,&
              &ES9.2E2,A1,A,A2,A6,I1,A4,ES10.4E2)') &
      '',   'xL(', iDim, ') = ', xL(iDim) / U % LengthX2Unit, &
      '', TRIM( U % LengthX2Label ), &
      ', ', 'xR(', iDim, ') = ', xR(iDim) / U % LengthX2Unit, &
      '', TRIM( U % LengthX2Label ), &
      ', ', 'ZoomX(', iDim, ') = ', ZoomX(iDim)

    WRITE(*,*)

    CALL CreateMesh &
           ( MeshX(iDim), nX(iDim), nNodesX(iDim), swX(iDim), &
             xL(iDim), xR(iDim), ZoomOption = ZoomX(iDim) )

    WRITE(*,'(A9,A11,I1,A4,ES9.2E2,A1,A,A3,ES9.2E2,A1,A)') &
      '', 'MIN/MAX dx(', iDim, ') = ', &
      MINVAL( MeshX(iDim) % Width(1:nX(iDim)) ) &
        / U % LengthX2Unit, '', TRIM( U % LengthX2Label ), &
      ' / ', &
      MAXVAL( MeshX(iDim) % Width(1:nX(iDim)) ) &
        / U % LengthX2Unit, '', TRIM( U % LengthX2Label )
    WRITE(*,*)

    iDim = 3
    WRITE(*,'(A7,A3,I1,A4,ES9.2E2,A1,A,A2,A3,I1,A4,&
              &ES9.2E2,A1,A,A2,A6,I1,A4,ES10.4E2)') &
      '',   'xL(', iDim, ') = ', xL(iDim) / U % LengthX3Unit, &
      '', TRIM( U % LengthX3Label ), &
      ', ', 'xR(', iDim, ') = ', xR(iDim) / U % LengthX3Unit, &
      '', TRIM( U % LengthX3Label ), &
      ', ', 'ZoomX(', iDim, ') = ', ZoomX(iDim)

    WRITE(*,*)

    CALL CreateMesh &
           ( MeshX(iDim), nX(iDim), nNodesX(iDim), swX(iDim), &
             xL(iDim), xR(iDim), ZoomOption = ZoomX(iDim) )

    WRITE(*,'(A9,A11,I1,A4,ES9.2E2,A1,A,A3,ES9.2E2,A1,A)') &
      '', 'MIN/MAX dx(', iDim, ') = ', &
      MINVAL( MeshX(iDim) % Width(1:nX(iDim)) ) &
        / U % LengthX3Unit, '', TRIM( U % LengthX3Label ), &
      ' / ', &
      MAXVAL( MeshX(iDim) % Width(1:nX(iDim)) ) &
        / U % LengthX3Unit, '', TRIM( U % LengthX3Label )
    WRITE(*,*)

    END ASSOCIATE ! U

    ! --- Geometry (Position Space) ---

    CALL CreateScalarFields ( nX, swX )

  END SUBROUTINE InitializeProgram


  SUBROUTINE FinalizeProgram

    INTEGER :: iDim

    ! --- Spatial Grid ---

    DO iDim = 1, 3

      CALL DestroyMesh( MeshX(iDim) )

    END DO

    ! --- Scalar Fields ---

    CALL DestroyScalarFields

    ! --- Device ---

    CALL FinalizeDevice

    CALL MPI_FINALIZE( mpierr )

  END SUBROUTINE FinalizeProgram


END MODULE ScalarWave_ProgramInitializationModule
