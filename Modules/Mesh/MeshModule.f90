MODULE MeshModule

  USE KindModule, ONLY: &
    DP
  USE QuadratureModule, ONLY: &
    GetQuadrature

  IMPLICIT NONE
  PRIVATE

  TYPE, PUBLIC :: MeshType
    REAL(DP)                            :: Length
    REAL(DP), DIMENSION(:), ALLOCATABLE :: Center
    REAL(DP), DIMENSION(:), ALLOCATABLE :: Width
    REAL(DP), DIMENSION(:), ALLOCATABLE :: Nodes
    REAL(DP), DIMENSION(:), ALLOCATABLE :: Weights
  END type MeshType

  TYPE(MeshType), DIMENSION(3), PUBLIC :: MeshX ! Spatial Mesh
  TYPE(MeshType),               PUBLIC :: MeshE ! Spectral Mesh

  PUBLIC :: CreateMesh
  PUBLIC :: DestroyMesh
  PUBLIC :: NodeCoordinate

CONTAINS


  SUBROUTINE CreateMesh( Mesh, N, nN, SW, xL, xR, ZoomOption )

    TYPE(MeshType)                 :: Mesh
    INTEGER, INTENT(in)            :: N, nN, SW
    REAL(DP), INTENT(in)           :: xL, xR
    REAL(DP), INTENT(in), OPTIONAL :: ZoomOption

    REAL(DP), DIMENSION(nN) :: xQ, wQ

    Mesh % Length = xR - xL

    ALLOCATE( Mesh % Center(1-SW:N+SW) )
    ALLOCATE( Mesh % Width (1-SW:N+SW) )

    IF( PRESENT( ZoomOption ) .AND. ZoomOption > 1.0_DP )THEN

      CALL CreateGeometricMesh &
             ( N, SW, xL, xR, Mesh % Center, Mesh % Width, &
               Zoom = ZoomOption )

    ELSE

      CALL CreateEquidistantMesh &
             ( N, SW, xL, xR, Mesh % Center, Mesh % Width )

    END IF

    CALL GetQuadrature( nN, xQ, wQ )

    ALLOCATE( Mesh % Nodes  (1:nN) )
    Mesh % Nodes = xQ

    ALLOCATE( Mesh % Weights(1:nN) )
    Mesh % Weights = wQ

  END SUBROUTINE CreateMesh


  SUBROUTINE CreateEquidistantMesh( N, SW, xL, xR, Center, Width )

    INTEGER,                        INTENT(in)    :: N, SW
    REAL(DP),                       INTENT(in)    :: xL, xR
    REAL(DP), DIMENSION(1-SW:N+SW), INTENT(inout) :: Center, Width

    INTEGER :: i

    Width(:) = ( xR - xL ) / REAL( N )

    Center(1) = xL + 0.5_DP * Width(1)
    DO i = 2, N
      Center(i) = Center(i-1) + Width(i-1)
    END DO

    DO i = 0, 1 - SW, - 1
      Center(i) = Center(i+1) - Width(i+1)
    END DO

    DO i = N + 1, N + SW
      Center(i) = Center(i-1) + Width(i-1)
    END DO

  END SUBROUTINE CreateEquidistantMesh


  SUBROUTINE CreateGeometricMesh( N, SW, xL, xR, Center, Width, Zoom )

    INTEGER,                        INTENT(in)    :: N, SW
    REAL(DP),                       INTENT(in)    :: xL, xR, Zoom
    REAL(DP), DIMENSION(1-SW:N+SW), INTENT(inout) :: Center, Width

    INTEGER :: i

    Width (1) = ( xR - xL ) * ( Zoom - 1.0_DP ) / ( Zoom**N - 1.0_DP )
    Center(1) = xL + 0.5_DP * Width(1)
    DO i = 2, N
      Width (i) = Width(i-1) * Zoom
      Center(i) = xL + SUM( Width(1:i-1) ) + 0.5_DP * Width(i)
    END DO

    DO i = 0, 1 - SW, - 1
      Width (i) = Width(i+1) / Zoom
      Center(i) = xL - SUM( Width(i+1:1-SW) ) - 0.5_DP * Width(i)
    END DO

    DO i = N + 1, N + SW
      Width (i) = Width(i-1) * Zoom
      Center(i) = xL + SUM( Width(1:i-1) ) + 0.5_DP * Width(i)
    END DO

  END SUBROUTINE CreateGeometricMesh


  SUBROUTINE DestroyMesh( Mesh )

    TYPE(MeshType) :: Mesh

    DEALLOCATE( Mesh % Center )
    DEALLOCATE( Mesh % Width  )
    DEALLOCATE( Mesh % Nodes  )

  END SUBROUTINE DestroyMesh


  PURE REAL(DP) FUNCTION NodeCoordinate( Mesh, iC, iN )

    TYPE(MeshType), INTENT(in) :: Mesh
    INTEGER,        INTENT(in) :: iC, iN

    NodeCoordinate &
      = Mesh % Center(iC) + Mesh % Width(iC) * Mesh % Nodes(iN)

    RETURN
  END FUNCTION NodeCoordinate


END MODULE MeshModule
