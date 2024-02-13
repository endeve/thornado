MODULE MF_GeometryModule

  ! --- AMReX Modules ---

  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_geometry_module, ONLY: &
    amrex_is_all_periodic
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy
  USE amrex_parmparse_module, ONLY: &
    amrex_parmparse, &
    amrex_parmparse_build, &
    amrex_parmparse_destroy
  USE amrex_amrcore_module, ONLY: &
    amrex_geom

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    nNodesX, &
    swX, &
    bcX
  USE UtilitiesModule, ONLY: &
    NodeNumberX
  USE ReferenceElementModuleX, ONLY: &
    nDOFX_X1, &
    NodeNumberTableX
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    LX_X1_Up
  USE UnitsModule, ONLY: &
    SolarMass
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    iGF_h_1, &
    iGF_h_2, &
    iGF_h_3, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_SqrtGm, &
    iGF_Psi, &
    iGF_Beta_1, &
    iGF_Beta_2, &
    nGF, &
    CoordinateSystem
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX
  USE GravitySolutionModule_Newtonian_PointMass, ONLY: &
    ComputeGravitationalPotential
  USE LinearAlgebraModule, ONLY: &
    MatrixMatrixMultiply

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero, &
    SqrtTiny, &
    One
  USE MF_UtilitiesModule, ONLY: &
    thornado2amrex_X, &
    AllocateArray_X, &
    DeallocateArray_X
  USE InputParsingModule, ONLY: &
    nLevels, &
    UseTiling, &
    ProgramName, &
    SolveGravity_NR
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeGeometryX_MF
  PUBLIC :: ApplyBoundaryConditions_Geometry_MF
  PUBLIC :: UpdateSpatialMetric_MF

  INTERFACE ComputeGeometryX_MF
    MODULE PROCEDURE ComputeGeometryX_MF_AllLevels
    MODULE PROCEDURE ComputeGeometryX_MF_SingleLevel
  END INTERFACE ComputeGeometryX_MF

  INTERFACE ApplyBoundaryConditions_Geometry_MF
    MODULE PROCEDURE ApplyBoundaryConditions_Geometry_MF_AllLevels
    MODULE PROCEDURE ApplyBoundaryConditions_Geometry_MF_SingleLevel
  END INTERFACE ApplyBoundaryConditions_Geometry_MF

  INTERFACE UpdateSpatialMetric_MF
    MODULE PROCEDURE UpdateSpatialMetric_MF_AllLevels
    MODULE PROCEDURE UpdateSpatialMetric_MF_SingleLevel
  END INTERFACE UpdateSpatialMetric_MF

  INTEGER :: swX_GF(3)

CONTAINS


  SUBROUTINE ComputeGeometryX_MF_AllLevels( MF_uGF )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:)

    INTEGER :: iLevel

    DO iLevel = 0, nLevels-1

      CALL CreateMesh_MF( iLevel, MeshX )

      CALL ComputeGeometryX_MF_SingleLevel( MF_uGF(iLevel) )

      CALL DestroyMesh_MF( MeshX )

    END DO

  END SUBROUTINE ComputeGeometryX_MF_AllLevels


  SUBROUTINE ComputeGeometryX_MF_SingleLevel( MF_uGF )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF

    TYPE(amrex_mfiter)    :: MFI
    TYPE(amrex_box)       :: BX
    TYPE(amrex_parmparse) :: PP

    INTEGER :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)

    REAL(DP), ALLOCATABLE :: G (:,:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)

    REAL(DP) :: Mass

    Mass = Zero

    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL PP % query( 'Mass', Mass )
    CALL amrex_parmparse_destroy( PP )

    Mass = Mass * SolarMass

#if defined( THORNADO_OMP )
    !$OMP PARALLEL &
    !$OMP PRIVATE( MFI, BX, iX_B0, iX_E0, iX_B1, iX_E1, G, uGF )
#endif

    CALL amrex_mfiter_build( MFI, MF_uGF )

    DO WHILE( MFI % next() )

      uGF => MF_uGF % DataPtr( MFI )

      BX = MFI % TileBox()

      iX_B0 = BX % lo
      iX_E0 = BX % hi
      iX_B1 = iX_B0 - swX
      iX_E1 = iX_E0 + swX

      CALL AllocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
               G )

      G = Zero ! Uninitialized variables cause crash in IO in DEBUG mode

#if defined HYDRO_RELATIVISTIC

      CALL ComputeGeometryX &
             ( iX_B0, iX_E0, iX_B1, iX_E1, G, Mass_Option = Mass )

#else

      CALL ComputeGeometryX &
             ( iX_B0, iX_E0, iX_B1, iX_E1, G )

      IF( SolveGravity_NR ) &
        CALL ComputeGravitationalPotential &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G, Mass )

#endif

      CALL thornado2amrex_X &
             ( nGF, iX_B1, iX_E1, LBOUND( uGF ), iX_B1, iX_E1, uGF, G )

      CALL DeallocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
               G )

    END DO

    CALL amrex_mfiter_destroy( MFI )

#if defined( THORNADO_OMP )
    !$OMP END PARALLEL
#endif

  END SUBROUTINE ComputeGeometryX_MF_SingleLevel


  SUBROUTINE ApplyBoundaryConditions_Geometry_MF_AllLevels( MF_uGF )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:)

    INTEGER :: iLevel

    DO iLevel = 0, nLevels-1

      CALL ApplyBoundaryConditions_Geometry_MF_SingleLevel &
             ( iLevel, MF_uGF(iLevel) )

    END DO

  END SUBROUTINE ApplyBoundaryConditions_Geometry_MF_AllLevels


  SUBROUTINE ApplyBoundaryConditions_Geometry_MF_SingleLevel( iLevel, MF_uGF )

    INTEGER             , INTENT(in)    :: iLevel
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF

    swX_GF = swX

    IF( .NOT. amrex_is_all_periodic() )THEN

      CALL ApplyBoundaryConditions_Geometry_MF_X1( iLevel, MF_uGF )
      CALL ApplyBoundaryConditions_Geometry_MF_X2( iLevel, MF_uGF )

    END IF

  END SUBROUTINE ApplyBoundaryConditions_Geometry_MF_SingleLevel


  SUBROUTINE UpdateSpatialMetric_MF_AllLevels( MF_uGF, swX_Option )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:)
    INTEGER             , INTENT(in), OPTIONAL :: swX_Option(3)

    INTEGER  :: iLevel, swXX(3)

    swXX = 0
    IF( PRESENT( swX_Option ) ) &
      swXX = swX_Option

    DO iLevel = 0, nLevels-1

      CALL UpdateSpatialMetric_MF_SingleLevel &
             ( iLevel, MF_uGF(iLevel), swX_Option = swXX )

    END DO ! iLevel = 0, nLevels-1

  END SUBROUTINE UpdateSpatialMetric_MF_AllLevels


  SUBROUTINE UpdateSpatialMetric_MF_SingleLevel( iLevel, MF_uGF, swX_Option )

    INTEGER             , INTENT(in)    :: iLevel
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF
    INTEGER             , INTENT(in), OPTIONAL :: swX_Option(3)

    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)

    INTEGER  :: iNX, iX1, iX2, iX3, iNX1, iNX2, swXX(3)
    INTEGER  :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP) :: X1, X2, Psi, h1, h2, h3

    swXX = 0
    IF( PRESENT( swX_Option ) ) &
      swXX = swX_Option

    CALL amrex_mfiter_build( MFI, MF_uGF, tiling = UseTiling )

    CALL CreateMesh_MF( iLevel, MeshX )

    DO WHILE( MFI % next() )

      uGF => MF_uGF % DataPtr( MFI )

      BX = MFI % tilebox()

      iX_B0 = BX % lo
      iX_E0 = BX % hi
      iX_B1 = iX_B0 - swXX
      iX_E1 = iX_E0 + swXX

      DO iX3 = iX_B1(3), iX_E1(3)
      DO iX2 = iX_B1(2), iX_E1(2)
      DO iX1 = iX_B1(1), iX_E1(1)
      DO iNX = 1       , nDOFX

        iNX1 = NodeNumberTableX(1,iNX)
        iNX2 = NodeNumberTableX(2,iNX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNX1 )
        X2 = NodeCoordinate( MeshX(2), iX2, iNX2 )

        Psi = uGF(iX1,iX2,iX3,nDOFX*(iGF_Psi-1)+iNX)

        IF( TRIM( CoordinateSystem ) .EQ. 'SPHERICAL' )THEN

          h1 = Psi**2
          h2 = Psi**2 * ABS( X1 )
          h3 = Psi**2 * ABS( X1 * SIN( X2 ) )

        ELSE IF( TRIM( CoordinateSystem ) .EQ. 'CYLINDIRCAL' )THEN

          h1 = Psi**2
          h2 = Psi**2
          h3 = Psi**2 * ABS( X1  )

        ELSE

          h1 = Psi**2
          h2 = Psi**2
          h3 = Psi**2

        END IF

        uGF(iX1,iX2,iX3,nDOFX*(iGF_h_1-1)+iNX) = h1
        uGF(iX1,iX2,iX3,nDOFX*(iGF_h_2-1)+iNX) = h2
        uGF(iX1,iX2,iX3,nDOFX*(iGF_h_3-1)+iNX) = h3

        uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_11-1)+iNX) = MAX( h1**2, SqrtTiny )
        uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_22-1)+iNX) = MAX( h2**2, SqrtTiny )
        uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_33-1)+iNX) = MAX( h3**2, SqrtTiny )

        uGF(iX1,iX2,iX3,nDOFX*(iGF_SqrtGm-1)+iNX) &
          = MAX( h1 * h2 * h3, SqrtTiny )

      END DO
      END DO
      END DO
      END DO

    END DO ! WHILE( MFI % next() )

    CALL DestroyMesh_MF( MeshX )

    CALL amrex_mfiter_destroy( MFI )

  END SUBROUTINE UpdateSpatialMetric_MF_SingleLevel


  ! --- PRIVATE SUBROUTINES ---


  SUBROUTINE ApplyBoundaryConditions_Geometry_MF_X1 &
    ( iLevel, MF_uGF )

    INTEGER             , INTENT(in)    :: iLevel
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF

    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)

    INTEGER :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)

#if defined( THORNADO_OMP )
    !$OMP PARALLEL &
    !$OMP PRIVATE( BX, MFI, uGF, iX_B0, iX_E0, iX_B1, iX_E1 )
#endif

    CALL amrex_mfiter_build( MFI, MF_uGF, tiling = UseTiling )

    DO WHILE( MFI % next() )

      uGF => MF_uGF % DataPtr( MFI )

      BX = MFI % tilebox()

      iX_B0 = BX % lo
      iX_E0 = BX % hi
      iX_B1 = iX_B0 - swX
      iX_E1 = iX_E0 + swX

      ! --- Lower boundary ---

      IF( iX_B0(1) .EQ. amrex_geom(iLevel) % domain % lo( 1 ) )THEN

        IF( bcX(1) .EQ. 02 )THEN

          CALL ApplyBoundaryConditions_Geometry_X1_Inner_Homogeneous &
                 ( iX_B0, iX_E0, iX_B1, iX_E1, uGF )

        ELSE IF( bcX(1) .EQ. 30 .OR. bcX(1) .EQ. 03 )THEN

          CALL ApplyBoundaryConditions_Geometry_X1_Inner_Reflecting &
                 ( iX_B0, iX_E0, iX_B1, iX_E1, uGF )

        END IF

      END IF !  Lower boundary

      ! --- Higher boundary ---

      IF( iX_E0(1) .EQ. amrex_geom(iLevel) % domain % hi( 1 ) )THEN

        IF( bcX(1) .EQ. 02 )THEN

          CALL ApplyBoundaryConditions_Geometry_X1_Outer_Homogeneous &
                 ( iX_B0, iX_E0, iX_B1, iX_E1, uGF )

        ELSE IF( bcX(1) .EQ. 03 )THEN

          CALL ApplyBoundaryConditions_Geometry_X1_Outer_Reflecting &
                 ( iX_B0, iX_E0, iX_B1, iX_E1, uGF )

        ! --- Extrapolate ---

        ELSE IF( bcX(1) .EQ. 30 )THEN

          CALL ApplyBoundaryConditions_Geometry_X1_Outer_Extrapolate &
                 ( iX_B0, iX_E0, iX_B1, iX_E1, uGF )

        END IF

      END IF ! Higher boundary

    END DO

    CALL amrex_mfiter_destroy( MFI )

#if defined( THORNADO_OMP )
    !$OMP END PARALLEL
#endif

  END SUBROUTINE ApplyBoundaryConditions_Geometry_MF_X1


  SUBROUTINE ApplyBoundaryConditions_Geometry_MF_X2 &
    ( iLevel, MF_uGF )

    INTEGER             , INTENT(in)    :: iLevel
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF

    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)

    INTEGER :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)

#if defined( THORNADO_OMP )
    !$OMP PARALLEL &
    !$OMP PRIVATE( BX, MFI, uGF, iX_B0, iX_E0, iX_B1, iX_E1 )
#endif

    CALL amrex_mfiter_build( MFI, MF_uGF, tiling = UseTiling )

    DO WHILE( MFI % next() )

      uGF => MF_uGF % DataPtr( MFI )

      BX = MFI % tilebox()

      iX_B0 = BX % lo
      iX_E0 = BX % hi
      iX_B1 = iX_B0 - swX
      iX_E1 = iX_E0 + swX

      ! --- Lower boundary ---

      IF( iX_B0(2) .EQ. amrex_geom(iLevel) % domain % lo( 2 ) )THEN

        IF( bcX(2) .EQ. 02 )THEN

          CALL ApplyBoundaryConditions_Geometry_X2_Inner_Homogeneous &
                 ( iX_B0, iX_E0, iX_B1, iX_E1, uGF )

        ELSE IF( bcX(2) .EQ. 03 )THEN

          CALL ApplyBoundaryConditions_Geometry_X2_Inner_Reflecting &
                 ( iX_B0, iX_E0, iX_B1, iX_E1, uGF )

        END IF

      END IF !  Lower boundary

      ! --- Higher boundary ---

      IF( iX_E0(2) .EQ. amrex_geom(iLevel) % domain % hi( 2 ) )THEN

        IF( bcX(2) .EQ. 02 )THEN

          CALL ApplyBoundaryConditions_Geometry_X2_Outer_Homogeneous &
                 ( iX_B0, iX_E0, iX_B1, iX_E1, uGF )

        ELSE IF( bcX(2) .EQ. 03 )THEN

          CALL ApplyBoundaryConditions_Geometry_X2_Outer_Reflecting &
                 ( iX_B0, iX_E0, iX_B1, iX_E1, uGF )

        END IF

      END IF ! Higher boundary

    END DO

    CALL amrex_mfiter_destroy( MFI )

#if defined( THORNADO_OMP )
    !$OMP END PARALLEL
#endif

  END SUBROUTINE ApplyBoundaryConditions_Geometry_MF_X2


  SUBROUTINE ApplyBoundaryConditions_Geometry_X1_Inner_Homogeneous &
    ( iX_B0, iX_E0, iX_B1, iX_E1, uGF )

    INTEGER , INTENT(in)    :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout) :: uGF(iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER :: iNX, iX2, iX3, iGF, iX_B(3), iX_E(3)

    iX_B = iX_B0 - swX_GF
    iX_E = iX_E0 + swX_GF

#if defined( THORNADO_OMP )
    !$OMP PARALLEL DO COLLAPSE(4)
#endif
    DO iGF  = 1      , nGF
    DO iX3  = iX_B(3), iX_E(3)
    DO iX2  = iX_B(2), iX_E(2)
    DO iNX  = 1      , nDOFX

      uGF(iX_B0(1)-1,iX2,iX3,nDOFX*(iGF-1)+iNX) &
        = uGF(iX_B0(1),iX2,iX3,nDOFX*(iGF-1)+iNX)

    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE ApplyBoundaryConditions_Geometry_X1_Inner_Homogeneous


  SUBROUTINE ApplyBoundaryConditions_Geometry_X1_Outer_Homogeneous &
    ( iX_B0, iX_E0, iX_B1, iX_E1, uGF )

    INTEGER , INTENT(in)    :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout) :: uGF(iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER :: iNX, iX2, iX3, iGF, iX_B(3), iX_E(3)

    iX_B = iX_B0 - swX_GF
    iX_E = iX_E0 + swX_GF

#if defined( THORNADO_OMP )
    !$OMP PARALLEL DO COLLAPSE(4)
#endif
    DO iGF  = 1      , nGF
    DO iX3  = iX_B(3), iX_E(3)
    DO iX2  = iX_B(2), iX_E(2)
    DO iNX  = 1      , nDOFX

      uGF(iX_E0(1)+1,iX2,iX3,nDOFX*(iGF-1)+iNX) &
        = uGF(iX_E0(1),iX2,iX3,nDOFX*(iGF-1)+iNX)

    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE ApplyBoundaryConditions_Geometry_X1_Outer_Homogeneous


  SUBROUTINE ApplyBoundaryConditions_Geometry_X1_Inner_Reflecting &
    ( iX_B0, iX_E0, iX_B1, iX_E1, uGF )

    INTEGER , INTENT(in)    :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout) :: uGF(iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER :: iNX1, iNX2, iNX3, iX2, iX3, iGF, iX_B(3), iX_E(3)
    INTEGER :: iNX, jNX, jNX1

    iX_B = iX_B0 - swX_GF
    iX_E = iX_E0 + swX_GF

#if defined( THORNADO_OMP )
    !$OMP PARALLEL DO COLLAPSE(6) &
    !$OMP PRIVATE( iNX, jNX, jNX1 )
#endif
    DO iGF  = 1      , nGF
    DO iX3  = iX_B(3), iX_E(3)
    DO iX2  = iX_B(2), iX_E(2)
    DO iNX3 = 1      , nNodesX(3)
    DO iNX2 = 1      , nNodesX(2)
    DO iNX1 = 1      , nNodesX(1)

      jNX1 = ( nNodesX(1) - iNX1 ) + 1

      iNX = NodeNumberX( iNX1, iNX2, iNX3 )
      jNX = NodeNumberX( jNX1, iNX2, iNX3 )

      uGF(iX_B0(1)-1,iX2,iX3,nDOFX*(iGF-1)+iNX) &
        = uGF(iX_B0(1),iX2,iX3,nDOFX*(iGF-1)+jNX)

      IF( iGF .EQ. iGF_Beta_1 ) &
        uGF(iX_B0(1)-1,iX2,iX3,nDOFX*(iGF-1)+iNX) &
          = -uGF(iX_B0(1),iX2,iX3,nDOFX*(iGF-1)+jNX)

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE ApplyBoundaryConditions_Geometry_X1_Inner_Reflecting


  SUBROUTINE ApplyBoundaryConditions_Geometry_X1_Outer_Reflecting &
    ( iX_B0, iX_E0, iX_B1, iX_E1, uGF )

    INTEGER , INTENT(in)    :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout) :: uGF(iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER :: iNX1, iNX2, iNX3, iX2, iX3, iGF, iX_B(3), iX_E(3)
    INTEGER :: iNX, jNX, jNX1

    iX_B = iX_B0 - swX_GF
    iX_E = iX_E0 + swX_GF

#if defined( THORNADO_OMP )
    !$OMP PARALLEL DO COLLAPSE(6) &
    !$OMP PRIVATE( iNX, jNX, jNX1 )
#endif
    DO iGF  = 1      , nGF
    DO iX3  = iX_B(3), iX_E(3)
    DO iX2  = iX_B(2), iX_E(2)
    DO iNX3 = 1      , nNodesX(3)
    DO iNX2 = 1      , nNodesX(2)
    DO iNX1 = 1      , nNodesX(1)

      jNX1 = ( nNodesX(1) - iNX1 ) + 1

      iNX = NodeNumberX( iNX1, iNX2, iNX3 )
      jNX = NodeNumberX( jNX1, iNX2, iNX3 )

      uGF(iX_E0(1)+1,iX2,iX3,nDOFX*(iGF-1)+iNX) &
        = uGF(iX_E0(1),iX2,iX3,nDOFX*(iGF-1)+jNX)

      IF( iGF .EQ. iGF_Beta_1 ) &
        uGF(iX_E0(1)+1,iX2,iX3,nDOFX*(iGF-1)+iNX) &
          = -uGF(iX_E0(1),iX2,iX3,nDOFX*(iGF-1)+jNX)

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE ApplyBoundaryConditions_Geometry_X1_Outer_Reflecting


  SUBROUTINE ApplyBoundaryConditions_Geometry_X1_Outer_Extrapolate &
    ( iX_B0, iX_E0, iX_B1, iX_E1, uGF )

    INTEGER , INTENT(in)    :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout) :: uGF(iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    REAL(DP), ALLOCATABLE :: G_K(:,:,:,:)
    REAL(DP), ALLOCATABLE :: G_F(:,:,:,:)

    INTEGER :: iNX, iX2, iX3, iGF, nX1_X, iX_B(3), iX_E(3)

    iX_B = iX_B0 - swX_GF
    iX_E = iX_E0 + swX_GF

    nX1_X = ( iX_E(3) - iX_B(3) + 1 ) * ( iX_E(2) - iX_B(2) + 1 )

    ALLOCATE( G_K(1:nDOFX   ,iX_B(2):iX_E(2),iX_B(3):iX_E(3),1:nGF) )
    ALLOCATE( G_F(1:nDOFX_X1,iX_B(2):iX_E(2),iX_B(3):iX_E(3),1:nGF) )

#if defined( THORNADO_OMP )
    !$OMP PARALLEL DO COLLAPSE(4)
#endif
    DO iGF = 1      , nGF
    DO iX3 = iX_B(3), iX_E(3)
    DO iX2 = iX_B(2), iX_E(2)
    DO iNX = 1      , nDOFX

      G_K(iNX,iX2,iX3,iGF) = uGF(iX_E0(1),iX2,iX3,nDOFX*(iGF-1)+iNX)

    END DO
    END DO
    END DO
    END DO

    DO iGF = 1, nGF

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, One, LX_X1_Up, &
               nDOFX_X1,   G_K(1,iX_B(2),iX_B(3),iGF), &
               nDOFX, Zero,G_F(1,iX_B(2),iX_B(3),iGF), &
               nDOFX_X1 )

    END DO

#if defined( THORNADO_OMP )
    !$OMP PARALLEL DO COLLAPSE(4)
#endif
    DO iGF = 1      , nGF
    DO iX3 = iX_B(3), iX_E(3)
    DO iX2 = iX_B(2), iX_E(2)
    DO iNX = 1      , nDOFX

      uGF(iX_E0(1)+1,iX2,iX3,nDOFX*(iGF-1)+iNX) = G_F(1,iX2,iX3,iGF)

    END DO
    END DO
    END DO
    END DO

    DEALLOCATE( G_F )
    DEALLOCATE( G_K )

  END SUBROUTINE ApplyBoundaryConditions_Geometry_X1_Outer_Extrapolate


  SUBROUTINE ApplyBoundaryConditions_Geometry_X2_Inner_Homogeneous &
    ( iX_B0, iX_E0, iX_B1, iX_E1, uGF )

    INTEGER , INTENT(in)    :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout) :: uGF(iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER :: iNX, iX1, iX3, iGF, iX_B(3), iX_E(3)

    iX_B = iX_B0 - swX_GF
    iX_E = iX_E0 + swX_GF

#if defined( THORNADO_OMP )
    !$OMP PARALLEL DO COLLAPSE(4)
#endif
    DO iGF  = 1      , nGF
    DO iX3  = iX_B(3), iX_E(3)
    DO iX1  = iX_B(1), iX_E(1)
    DO iNX  = 1      , nDOFX

      uGF(iX1,iX_B0(2)-1,iX3,nDOFX*(iGF-1)+iNX) &
        = uGF(iX1,iX_B0(2),iX3,nDOFX*(iGF-1)+iNX)

    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE ApplyBoundaryConditions_Geometry_X2_Inner_Homogeneous


  SUBROUTINE ApplyBoundaryConditions_Geometry_X2_Outer_Homogeneous &
    ( iX_B0, iX_E0, iX_B1, iX_E1, uGF )

    INTEGER , INTENT(in)    :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout) :: uGF(iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER :: iNX, iX1, iX3, iGF, iX_B(3), iX_E(3)

    iX_B = iX_B0 - swX_GF
    iX_E = iX_E0 + swX_GF

#if defined( THORNADO_OMP )
    !$OMP PARALLEL DO COLLAPSE(4)
#endif
    DO iGF  = 1      , nGF
    DO iX3  = iX_B(3), iX_E(3)
    DO iX1  = iX_B(1), iX_E(1)
    DO iNX  = 1      , nDOFX

      uGF(iX1,iX_E0(2)+1,iX3,nDOFX*(iGF-1)+iNX) &
        = uGF(iX1,iX_E0(2),iX3,nDOFX*(iGF-1)+iNX)

    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE ApplyBoundaryConditions_Geometry_X2_Outer_Homogeneous


  SUBROUTINE ApplyBoundaryConditions_Geometry_X2_Inner_Reflecting &
    ( iX_B0, iX_E0, iX_B1, iX_E1, uGF )

    INTEGER , INTENT(in)    :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout) :: uGF(iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER :: iNX1, iNX2, iNX3, iX1, iX3, iGF, iX_B(3), iX_E(3)
    INTEGER :: iNX, jNX, jNX2

    iX_B = iX_B0 - swX_GF
    iX_E = iX_E0 + swX_GF

#if defined( THORNADO_OMP )
    !$OMP PARALLEL DO COLLAPSE(6) &
    !$OMP PRIVATE( iNX, jNX, jNX2 )
#endif
    DO iGF  = 1      , nGF
    DO iX3  = iX_B(3), iX_E(3)
    DO iX1  = iX_B(1), iX_E(1)
    DO iNX3 = 1      , nNodesX(3)
    DO iNX2 = 1      , nNodesX(2)
    DO iNX1 = 1      , nNodesX(1)

      jNX2 = ( nNodesX(2) - iNX2 ) + 1

      iNX = NodeNumberX( iNX1, iNX2, iNX3 )
      jNX = NodeNumberX( iNX1, jNX2, iNX3 )

      uGF(iX1,iX_B0(2)-1,iX3,nDOFX*(iGF-1)+iNX) &
        = uGF(iX1,iX_B0(1),iX3,nDOFX*(iGF-1)+jNX)

      IF( iGF .EQ. iGF_Beta_2 ) &
        uGF(iX1,iX_B0(2)-1,iX3,nDOFX*(iGF-1)+iNX) &
          = -uGF(iX1,iX_B0(2),iX3,nDOFX*(iGF-1)+jNX)

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE ApplyBoundaryConditions_Geometry_X2_Inner_Reflecting


  SUBROUTINE ApplyBoundaryConditions_Geometry_X2_Outer_Reflecting &
    ( iX_B0, iX_E0, iX_B1, iX_E1, uGF )

    INTEGER , INTENT(in)    :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout) :: uGF(iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER :: iNX1, iNX2, iNX3, iX1, iX3, iGF, iX_B(3), iX_E(3)
    INTEGER :: iNX, jNX, jNX2

    iX_B = iX_B0 - swX_GF
    iX_E = iX_E0 + swX_GF

#if defined( THORNADO_OMP )
    !$OMP PARALLEL DO COLLAPSE(6) &
    !$OMP PRIVATE( iNX, jNX, jNX2 )
#endif
    DO iGF  = 1      , nGF
    DO iX3  = iX_B(3), iX_E(3)
    DO iX1  = iX_B(1), iX_E(1)
    DO iNX3 = 1      , nNodesX(3)
    DO iNX2 = 1      , nNodesX(2)
    DO iNX1 = 1      , nNodesX(1)

      jNX2 = ( nNodesX(2) - iNX2 ) + 1

      iNX = NodeNumberX( iNX1, iNX2, iNX3 )
      jNX = NodeNumberX( iNX1, jNX2, iNX3 )

      uGF(iX1,iX_E0(2)+1,iX3,nDOFX*(iGF-1)+iNX) &
        = uGF(iX1,iX_E0(1),iX3,nDOFX*(iGF-1)+jNX)

      IF( iGF .EQ. iGF_Beta_2 ) &
        uGF(iX1,iX_E0(2)+1,iX3,nDOFX*(iGF-1)+iNX) &
          = -uGF(iX1,iX_E0(2),iX3,nDOFX*(iGF-1)+jNX)

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE ApplyBoundaryConditions_Geometry_X2_Outer_Reflecting


END MODULE MF_GeometryModule
