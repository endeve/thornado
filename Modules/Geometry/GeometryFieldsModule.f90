MODULE GeometryFieldsModule

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    nDOFX, nNodesX, nDOF

  IMPLICIT NONE
  PRIVATE

  CHARACTER(16), PUBLIC :: &
    CoordinateSystem = 'CARTESIAN'

  ! --- Weights for 'Phase Space' Fields ---

  REAL(DP), DIMENSION(:),         ALLOCATABLE, PUBLIC :: &
    WeightsG
  REAL(DP), DIMENSION(:,:,:,:),   ALLOCATABLE, PUBLIC :: &
    Vol
  REAL(DP), DIMENSION(:,:,:,:,:), ALLOCATABLE, PUBLIC :: &
    VolJac, &
    VolJacE

  ! --- Weights for 'Position Space' Fields ---

  REAL(DP), DIMENSION(:),       ALLOCATABLE, PUBLIC :: &
    WeightsGX, WeightsGX_X1, WeightsGX_X2, WeightsGX_X3
  REAL(DP), DIMENSION(:,:,:),   ALLOCATABLE, PUBLIC :: &
    VolX
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE, PUBLIC :: &
    VolJacX

  ! --- Spatial Geometry Fields ---

  INTEGER, PUBLIC, PARAMETER :: iGF_Phi_N    = 1  ! Newtonian Potential
  INTEGER, PUBLIC, PARAMETER :: iGF_Gm_dd_11 = 2  ! Spatial Metric Component 11
  INTEGER, PUBLIC, PARAMETER :: iGF_Gm_dd_22 = 3  ! Spatial Metric Component 22
  INTEGER, PUBLIC, PARAMETER :: iGF_Gm_dd_33 = 4  ! Spatial Metric Component 33
  INTEGER, PUBLIC, PARAMETER :: iGF_Gm_uu_11 = 5  ! Contravariant Spatial Metric Component 11
  INTEGER, PUBLIC, PARAMETER :: iGF_Alpha    = 6  ! Lapse Function
  INTEGER, PUBLIC, PARAMETER :: iGF_Beta_1   = 7  ! Shift Vector 1
  INTEGER, PUBLIC, PARAMETER :: iGF_Beta_2   = 8  ! Shift Vector 2
  INTEGER, PUBLIC, PARAMETER :: iGF_Beta_3   = 9  ! Shift Vector 3
  INTEGER, PUBLIC, PARAMETER :: iGF_CF       = 10 ! Conformal Factor
  INTEGER, PUBLIC, PARAMETER :: nGF          = 10 ! n Geometry Fields

  CHARACTER(32), DIMENSION(nGF), PUBLIC, PARAMETER :: &
    namesGF = [ 'Newtonian Potential                         ', &
                'Spatial Metric Component (11)               ', &
                'Spatial Metric Component (22)               ', &
                'Spatial Metric Component (33)               ', &
                'Contravariant Spatial Metric Component (11) ', &
                'Lapse Function                              ', &
                'Shift Vector (1)                            ', &
                'Shift Vector (2)                            ', &
                'Shift Vector (3)                            ', &
                'Conformal Factor                            ' ]

  REAL(DP), DIMENSION(:,:,:,:,:), ALLOCATABLE, PUBLIC :: uGF

  INTERFACE
    PURE REAL(DP) FUNCTION MetricFunction( X )
      USE KindModule, ONLY: DP
      REAL(DP), DIMENSION(3), INTENT(in) :: X
    END FUNCTION MetricFunction
  END INTERFACE

  PROCEDURE (MetricFunction), POINTER, PUBLIC :: a, b, c, d
  PROCEDURE (MetricFunction), POINTER, PUBLIC :: dlnadX1
  PROCEDURE (MetricFunction), POINTER, PUBLIC :: dlnbdX1
  PROCEDURE (MetricFunction), POINTER, PUBLIC :: dlncdX2

  PUBLIC :: CreateGeometryFields
  PUBLIC :: DestroyGeometryFields
  PUBLIC :: InitializeGeometryFields_CARTESIAN
  PUBLIC :: InitializeGeometryFields_SPHERICAL
  PUBLIC :: InitializeGeometryFields_CYLINDRICAL
  PUBLIC :: FinalizeGeometryFields

CONTAINS


  SUBROUTINE CreateGeometryFields( nX, swX, nE, swE )

    INTEGER, DIMENSION(3), INTENT(in) :: nX, swX
    INTEGER,               INTENT(in) :: nE, swE

    CALL CreateGeometryFieldsX( nX, swX )

    ALLOCATE( WeightsG (1:nDOF ) )

    ALLOCATE &
      ( Vol(1-swE   :nE   +swE,    1-swX(1):nX(1)+swX(1), &
            1-swX(2):nX(2)+swX(2), 1-swX(3):nX(3)+swX(3)) )

    ALLOCATE &
      ( VolJac &
          (1:nDOF, &
           1-swE   :nE   +swE,    1-swX(1):nX(1)+swX(1), &
           1-swX(2):nX(2)+swX(2), 1-swX(3):nX(3)+swX(3)) )

    ALLOCATE &
      ( VolJacE &
          (1:nDOF, &
           1-swE   :nE   +swE,    1-swX(1):nX(1)+swX(1), &
           1-swX(2):nX(2)+swX(2), 1-swX(3):nX(3)+swX(3)) )

  END SUBROUTINE CreateGeometryFields


  SUBROUTINE CreateGeometryFieldsX( nX, swX )

    INTEGER, DIMENSION(3), INTENT(in) :: nX, swX

    INTEGER :: iGF

    ALLOCATE( WeightsGX(1:nDOFX) )
    ALLOCATE( WeightsGX_X1(1:nNodesX(2)*nNodesX(3)))
    ALLOCATE( WeightsGX_X2(1:nNodesX(1)*nNodesX(3)))
    ALLOCATE( WeightsGX_X3(1:nNodesX(1)*nNodesX(2)))

    ALLOCATE &
      ( VolX(1-swX(1):nX(1)+swX(1), 1-swX(2):nX(2)+swX(2), &
             1-swX(3):nX(3)+swX(3)) )

    ALLOCATE &
      ( VolJacX(1:nDOFX, 1-swX(1):nX(1)+swX(1), &
                1-swX(2):nX(2)+swX(2), 1-swX(3):nX(3)+swX(3)) )

    WRITE(*,*)
    WRITE(*,'(A5,A15)') '', 'Geometry Fields'
    WRITE(*,*)
    DO iGF = 1, nGF
      WRITE(*,'(A5,A32)') '', TRIM( namesGF(iGF) )
    END DO

    ALLOCATE &
      ( uGF(1:nDOFX, &
            1-swX(1):nX(1)+swX(1), &
            1-swX(2):nX(2)+swX(2), &
            1-swX(3):nX(3)+swX(3), &
            1:nGF) )

    ! --- Initialize to Flat Spacetime (Cartesian) ---

    uGF(:,:,:,:,iGF_Phi_N)    = 0.0_DP
    uGF(:,:,:,:,iGF_Gm_dd_11) = 1.0_DP
    uGF(:,:,:,:,iGF_Gm_dd_22) = 1.0_DP
    uGF(:,:,:,:,iGF_Gm_dd_33) = 1.0_DP
    uGF(:,:,:,:,iGF_Alpha)    = 1.0_DP
    uGF(:,:,:,:,iGF_Beta_1)   = 0.0_DP
    uGF(:,:,:,:,iGF_Beta_2)   = 0.0_DP
    uGF(:,:,:,:,iGF_Beta_3)   = 0.0_DP
    uGF(:,:,:,:,iGF_CF)       = 1.0_DP
    uGF(:,:,:,:,iGF_Gm_uu_11) = 1.0_DP

  END SUBROUTINE CreateGeometryFieldsX


  SUBROUTINE DestroyGeometryFields

    CALL DestroyGeometryFieldsX

    DEALLOCATE( WeightsG, Vol, VolJac, VolJacE )

  END SUBROUTINE DestroyGeometryFields


  SUBROUTINE DestroyGeometryFieldsX

    DEALLOCATE( WeightsGX, WeightsGX_X1, WeightsGX_X2, WeightsGX_X3 )
    DEALLOCATE( VolX, VolJacX, uGF  )

  END SUBROUTINE DestroyGeometryFieldsX


  ! --- Coordinate System Dependent Metric Functions ---


  SUBROUTINE InitializeGeometryFields_CARTESIAN

    a => a_CARTESIAN
    b => b_CARTESIAN
    c => c_CARTESIAN
    d => SqrtDet_CARTESIAN

    dlnadX1 => dlnadX1_CARTESIAN
    dlnbdX1 => dlnbdX1_CARTESIAN
    dlncdX2 => dlncdX2_CARTESIAN

  END SUBROUTINE InitializeGeometryFields_CARTESIAN


  SUBROUTINE InitializeGeometryFields_SPHERICAL

    a => a_SPHERICAL
    b => b_SPHERICAL
    c => c_SPHERICAL
    d => SqrtDet_SPHERICAL

    dlnadX1 => dlnadX1_SPHERICAL
    dlnbdX1 => dlnbdX1_SPHERICAL
    dlncdX2 => dlncdX2_SPHERICAL

  END SUBROUTINE InitializeGeometryFields_SPHERICAL


  SUBROUTINE InitializeGeometryFields_CYLINDRICAL

    a => a_CYLINDRICAL
    b => b_CYLINDRICAL
    c => c_CYLINDRICAL
    d => SqrtDet_CYLINDRICAL

    dlnadX1 => dlnadX1_CYLINDRICAL
    dlnbdX1 => dlnbdX1_CYLINDRICAL
    dlncdX2 => dlncdX2_CYLINDRICAL

  END SUBROUTINE InitializeGeometryFields_CYLINDRICAL


  SUBROUTINE FinalizeGeometryFields

    NULLIFY( a, b, c, d )

  END SUBROUTINE FinalizeGeometryFields


  ! --- Cartesian Coordinates ---


  PURE REAL(DP) FUNCTION a_CARTESIAN( X )

    REAL(DP), DIMENSION(3), INTENT(in) :: X

    a_CARTESIAN = 1.0_DP

    RETURN
  END FUNCTION a_CARTESIAN


  PURE REAL(DP) FUNCTION b_CARTESIAN( X )

    REAL(DP), DIMENSION(3), INTENT(in) :: X

    b_CARTESIAN = 1.0_DP

    RETURN
  END FUNCTION b_CARTESIAN


  PURE REAL(DP) FUNCTION c_CARTESIAN( X )

    REAL(DP), DIMENSION(3), INTENT(in) :: X

    c_CARTESIAN = 1.0_DP

    RETURN
  END FUNCTION c_CARTESIAN


  PURE REAL(DP) FUNCTION SqrtDet_CARTESIAN( X )

    REAL(DP), DIMENSION(3), INTENT(in) :: X

    SqrtDet_CARTESIAN = 1.0_DP

    RETURN
  END FUNCTION SqrtDet_CARTESIAN


  PURE REAL(DP) FUNCTION dlnadX1_CARTESIAN( X )

    REAL(DP), DIMENSION(3), INTENT(in) :: X

    dlnadX1_CARTESIAN = 0.0_DP

    RETURN
  END FUNCTION dlnadX1_CARTESIAN


  PURE REAL(DP) FUNCTION dlnbdX1_CARTESIAN( X )

    REAL(DP), DIMENSION(3), INTENT(in) :: X

    dlnbdX1_CARTESIAN = 0.0_DP

    RETURN
  END FUNCTION dlnbdX1_CARTESIAN


  PURE REAL(DP) FUNCTION dlncdX2_CARTESIAN( X )

    REAL(DP), DIMENSION(3), INTENT(in) :: X

    dlncdX2_CARTESIAN = 0.0_DP

    RETURN
  END FUNCTION dlncdX2_CARTESIAN


  ! --- Spherical Coordinates ---


  PURE REAL(DP) FUNCTION a_SPHERICAL( X )

    REAL(DP), DIMENSION(3), INTENT(in) :: X

    a_SPHERICAL = X(1)

    RETURN
  END FUNCTION a_SPHERICAL


  PURE REAL(DP) FUNCTION b_SPHERICAL( X )

    REAL(DP), DIMENSION(3), INTENT(in) :: X

    b_SPHERICAL = X(1)

    RETURN
  END FUNCTION b_SPHERICAL


  PURE REAL(DP) FUNCTION c_SPHERICAL( X )

    REAL(DP), DIMENSION(3), INTENT(in) :: X

    c_SPHERICAL = SIN( X(2) )

    RETURN
  END FUNCTION c_SPHERICAL


  PURE REAL(DP) FUNCTION SqrtDet_SPHERICAL( X )

    REAL(DP), DIMENSION(3), INTENT(in) :: X

    SqrtDet_SPHERICAL = X(1)**2 * SIN( X(2) )

    RETURN
  END FUNCTION SqrtDet_SPHERICAL


  PURE REAL(DP) FUNCTION dlnadX1_SPHERICAL( X )

    REAL(DP), DIMENSION(3), INTENT(in) :: X

    dlnadX1_SPHERICAL = 1.0_DP / X(1)

    RETURN
  END FUNCTION dlnadX1_SPHERICAL


  PURE REAL(DP) FUNCTION dlnbdX1_SPHERICAL( X )

    REAL(DP), DIMENSION(3), INTENT(in) :: X

    dlnbdX1_SPHERICAL = 1.0_DP / X(1)

    RETURN
  END FUNCTION dlnbdX1_SPHERICAL


  PURE REAL(DP) FUNCTION dlncdX2_SPHERICAL( X )

    REAL(DP), DIMENSION(3), INTENT(in) :: X

    dlncdX2_SPHERICAL = 1.0_DP / TAN( X(2) )

    RETURN
  END FUNCTION dlncdX2_SPHERICAL


  ! --- Cylindrical Coordinates ---


  PURE REAL(DP) FUNCTION a_CYLINDRICAL( X )

    REAL(DP), DIMENSION(3), INTENT(in) :: X

    a_CYLINDRICAL = 1.0_DP

    RETURN
  END FUNCTION a_CYLINDRICAL


  PURE REAL(DP) FUNCTION b_CYLINDRICAL( X )

    REAL(DP), DIMENSION(3), INTENT(in) :: X

    b_CYLINDRICAL = X(1)

    RETURN
  END FUNCTION b_CYLINDRICAL


  PURE REAL(DP) FUNCTION c_CYLINDRICAL( X )

    REAL(DP), DIMENSION(3), INTENT(in) :: X

    c_CYLINDRICAL = 1.0_DP

    RETURN
  END FUNCTION c_CYLINDRICAL


  PURE REAL(DP) FUNCTION SqrtDet_CYLINDRICAL( X )

    REAL(DP), DIMENSION(3), INTENT(in) :: X

    SqrtDet_CYLINDRICAL = X(1)

    RETURN
  END FUNCTION SqrtDet_CYLINDRICAL


  PURE REAL(DP) FUNCTION dlnadX1_CYLINDRICAL( X )

    REAL(DP), DIMENSION(3), INTENT(in) :: X

    dlnadX1_CYLINDRICAL = 0.0_DP

    RETURN
  END FUNCTION dlnadX1_CYLINDRICAL


  PURE REAL(DP) FUNCTION dlnbdX1_CYLINDRICAL( X )

    REAL(DP), DIMENSION(3), INTENT(in) :: X

    dlnbdX1_CYLINDRICAL = 1.0_DP / X(1)

    RETURN
  END FUNCTION dlnbdX1_CYLINDRICAL


  PURE REAL(DP) FUNCTION dlncdX2_CYLINDRICAL( X )

    REAL(DP), DIMENSION(3), INTENT(in) :: X

    dlncdX2_CYLINDRICAL = 0.0_DP

    RETURN
  END FUNCTION dlncdX2_CYLINDRICAL


END MODULE GeometryFieldsModule
