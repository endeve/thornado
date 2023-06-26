MODULE GeometryComputationModule

  USE KindModule, ONLY: &
    DP, Zero, Half, One, SqrtTiny
  USE ProgramHeaderModule, ONLY: &
    nDOFX
  USE ReferenceElementModuleX, ONLY: &
    NodesLX_q
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    LX_L2G
  USE MeshModule, ONLY: &
    MeshX_mod => MeshX, &
    MeshType
  USE GeometryFieldsModule, ONLY: &
    CoordinateSystem_mod => CoordinateSystem, &
    iGF_Phi_N, &
    iGF_h_1,      iGF_h_2,      iGF_h_3,      &
    iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33, &
    iGF_Beta_1,   iGF_Beta_2,   iGF_Beta_3,   &
    iGF_SqrtGm, &
    iGF_Alpha, iGF_Psi, &
    nGF
  USE LinearAlgebraModule, ONLY: &
    MatrixMatrixMultiply

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeGeometryX
  PUBLIC :: ComputeGeometryX_FromScaleFactors
  PUBLIC :: ComputeGeometryX_SpatialMetric
  PUBLIC :: LapseFunction
  PUBLIC :: ConformalFactor


CONTAINS


  SUBROUTINE ComputeGeometryX &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, &
      Mass_Option, MeshX_Option, CoordinateSystem_Option )

    INTEGER, INTENT(in)            :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout)        :: &
      G(:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    REAL(DP), INTENT(in), OPTIONAL :: &
      Mass_Option
    TYPE(MeshType), INTENT(in), OPTIONAL   :: &
      MeshX_Option(3)
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: &
      CoordinateSystem_Option

    REAL(DP) :: Mass
    TYPE(MeshType) :: MeshX(3)
    CHARACTER(24) :: CoordinateSystem

    IF( PRESENT( Mass_Option ) )THEN
      Mass = Mass_Option
    ELSE
      Mass = Zero
    END IF

    IF( PRESENT( MeshX_Option ) )THEN
      MeshX = MeshX_Option
    ELSE
      MeshX = MeshX_mod
    END IF

    IF( PRESENT(CoordinateSystem_Option) )THEN

      IF( TRIM(CoordinateSystem_Option) == 'spherical' )THEN
        CoordinateSystem = 'SPHERICAL'
      ELSE IF( TRIM(CoordinateSystem_Option) == 'cylindrical' )THEN
        CoordinateSystem = 'CYLINDRICAL'
      ELSE IF( TRIM(CoordinateSystem_Option) == 'cartesian' )THEN
        CoordinateSystem = 'CARTESIAN'
      ELSE
        PRINT*, '[ComputeGeometryX] Invalid Coordinate System: ', &
                 CoordinateSystem_Option
      END IF

    ELSE
      CoordinateSystem = CoordinateSystem_mod
    END IF

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: G, iX_B1, iX_E1 )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA &
    !$ACC COPYIN(  G, iX_B1, iX_E1 )
#endif

    SELECT CASE ( TRIM( CoordinateSystem ) )

      CASE ( 'CARTESIAN' )

        CALL ComputeGeometryX_CARTESIAN &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G, Mass, MeshX )

      CASE ( 'SPHERICAL' )

        CALL ComputeGeometryX_SPHERICAL &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G, Mass, MeshX )

      CASE ( 'CYLINDRICAL' )

        CALL ComputeGeometryX_CYLINDRICAL &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G, Mass, MeshX )

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,'(A5,A27,A)') &
          '', 'Invalid Coordinate System: ', TRIM( CoordinateSystem )
        STOP

    END SELECT

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET UPDATE FROM( G )
#elif defined( THORNADO_OACC   )
    !$ACC UPDATE HOST       ( G )
#endif

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: G, iX_B1, iX_E1 )
#elif defined( THORNADO_OACC   )
    !$ACC EXIT DATA &
    !$ACC DELETE(       G, iX_B1, iX_E1 )
#endif

  END SUBROUTINE ComputeGeometryX


  SUBROUTINE ComputeGeometryX_CARTESIAN &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, Mass, MeshX )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout) :: &
      G(:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    REAL(DP), INTENT(in)    :: &
      Mass
    TYPE(MeshType), INTENT(in) :: &
      MeshX(3)

    INTEGER :: iX1, iX2, iX3, iNodeX

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4) &
    !$OMP MAP( to:     iX_B1, iX_E1 ) &
    !$OMP MAP( tofrom: G )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC COPYIN(      iX_B1, iX_E1 ) &
    !$ACC COPY(        G )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(4)
#endif
    DO iX3    = iX_B1(3), iX_E1(3)
    DO iX2    = iX_B1(2), iX_E1(2)
    DO iX1    = iX_B1(1), iX_E1(1)
    DO iNodeX = 1       , nDOFX

      ! --- Initialize to flat spacetime ---

      G(iNodeX,iX1,iX2,iX3,iGF_Phi_N)    = Zero
      G(iNodeX,iX1,iX2,iX3,iGF_h_1)      = One
      G(iNodeX,iX1,iX2,iX3,iGF_h_2)      = One
      G(iNodeX,iX1,iX2,iX3,iGF_h_3)      = One
      G(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11) = One
      G(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22) = One
      G(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33) = One
      G(iNodeX,iX1,iX2,iX3,iGF_SqrtGm)   = One
      G(iNodeX,iX1,iX2,iX3,iGF_Alpha)    = One
      G(iNodeX,iX1,iX2,iX3,iGF_Beta_1)   = Zero
      G(iNodeX,iX1,iX2,iX3,iGF_Beta_2)   = Zero
      G(iNodeX,iX1,iX2,iX3,iGF_Beta_3)   = Zero
      G(iNodeX,iX1,iX2,iX3,iGF_Psi)      = One

      CALL ComputeGeometryX_SpatialMetric &
             ( G(iNodeX,iX1,iX2,iX3,iGF_h_1), &
               G(iNodeX,iX1,iX2,iX3,iGF_h_2), &
               G(iNodeX,iX1,iX2,iX3,iGF_h_3), &
               G(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11), &
               G(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22), &
               G(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33), &
               G(iNodeX,iX1,iX2,iX3,iGF_SqrtGm) )

    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE ComputeGeometryX_CARTESIAN


  SUBROUTINE ComputeGeometryX_SPHERICAL &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, Mass, MeshX )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout) :: &
      G(nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),nGF)
    REAL(DP), INTENT(in)    :: &
      Mass
    TYPE(MeshType), INTENT(in) :: &
      MeshX(3)

    INTEGER  :: iX1, iX2, iX3, iNodeX
    INTEGER  :: nP_X, nX(3)
    REAL(DP) :: x1L_q, x2L_q, x1G_q, x2G_q
    REAL(DP) :: h_1_L  (nDOFX,iX_B1(1):iX_E1(1), &
                              iX_B1(2):iX_E1(2), &
                              iX_B1(3):iX_E1(3))
    REAL(DP) :: h_2_L  (nDOFX,iX_B1(1):iX_E1(1), &
                              iX_B1(2):iX_E1(2), &
                              iX_B1(3):iX_E1(3))
    REAL(DP) :: h_3_L  (nDOFX,iX_B1(1):iX_E1(1), &
                              iX_B1(2):iX_E1(2), &
                              iX_B1(3):iX_E1(3))
    REAL(DP) :: Alpha_L(nDOFX,iX_B1(1):iX_E1(1), &
                              iX_B1(2):iX_E1(2), &
                              iX_B1(3):iX_E1(3))
    REAL(DP) :: Psi_L  (nDOFX,iX_B1(1):iX_E1(1), &
                              iX_B1(2):iX_E1(2), &
                              iX_B1(3):iX_E1(3))

    nX   = iX_E1 - iX_B1 + 1
    nP_X = PRODUCT( nX )

    ASSOCIATE ( CenterX1 => MeshX(1) % Center, &
                CenterX2 => MeshX(2) % Center, &
                WidthX1  => MeshX(1) % Width, &
                WidthX2  => MeshX(2) % Width )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to:    G, iX_B1, iX_E1, CenterX1, CenterX2, WidthX1, WidthX2 ) &
    !$OMP MAP( alloc: h_1_L, h_2_L, h_3_L, Alpha_L, Psi_L )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA &
    !$ACC COPYIN(     G, iX_B1, iX_E1, CenterX1, CenterX2, WidthX1, WidthX2 ) &
    !$ACC CREATE(     h_1_L, h_2_L, h_3_L, Alpha_L, Psi_L )
#endif

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4) &
    !$OMP PRIVATE( x1L_q, x2L_q, x1G_q, x2G_q )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRIVATE( x1L_q, x2L_q, x1G_q, x2G_q )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(4) &
    !$OMP PRIVATE( x1L_q, x2L_q, x1G_q, x2G_q )
#endif
    DO iX3    = iX_B1(3), iX_E1(3)
    DO iX2    = iX_B1(2), iX_E1(2)
    DO iX1    = iX_B1(1), iX_E1(1)
    DO iNodeX = 1       , nDOFX

      ! --- Compute Geometry Fields in Lobatto Points ---

      ! --- Local Coordinates (Lobatto Points) ---

      x1L_q = NodesLX_q(1,iNodeX)
      x2L_q = NodesLX_q(2,iNodeX)

      ! --- Global Coordinates ---

      x1G_q = CenterX1(iX1) + WidthX1(iX1) * x1L_q
      x2G_q = CenterX2(iX2) + WidthX2(iX2) * x2L_q

      ! --- Compute Lapse Function and Conformal Factor ---

      Alpha_L(iNodeX,iX1,iX2,iX3) &
        = LapseFunction  ( x1G_q, Mass )
      Psi_L(iNodeX,iX1,iX2,iX3) &
        = ConformalFactor( x1G_q, Mass )

      ! --- Set Geometry in Lobatto Points ---

      h_1_L(iNodeX,iX1,iX2,iX3) &
        = Psi_L(iNodeX,iX1,iX2,iX3)**2
      h_2_L(iNodeX,iX1,iX2,iX3) &
        = Psi_L(iNodeX,iX1,iX2,iX3)**2 * x1G_q
      h_3_L(iNodeX,iX1,iX2,iX3) &
        = Psi_L(iNodeX,iX1,iX2,iX3)**2 * x1G_q * SIN( x2G_q )

    END DO
    END DO
    END DO
    END DO

    ! --- Interpolate from Lobatto to Gaussian Points ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX, nP_X, nDOFX, One, LX_L2G, nDOFX, &
             h_1_L, nDOFX, Zero, &
             G(1,iX_B1(1),iX_B1(2),iX_B1(3),iGF_h_1), nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX, nP_X, nDOFX, One, LX_L2G, nDOFX, &
             h_2_L, nDOFX, Zero, &
             G(1,iX_B1(1),iX_B1(2),iX_B1(3),iGF_h_2), nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX, nP_X, nDOFX, One, LX_L2G, nDOFX, &
             h_3_L, nDOFX, Zero, &
             G(1,iX_B1(1),iX_B1(2),iX_B1(3),iGF_h_3), nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX, nP_X, nDOFX, One, LX_L2G, nDOFX, &
             Alpha_L, nDOFX, Zero, &
             G(1,iX_B1(1),iX_B1(2),iX_B1(3),iGF_Alpha), nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX, nP_X, nDOFX, One, LX_L2G, nDOFX, &
             Psi_L, nDOFX, Zero, &
             G(1,iX_B1(1),iX_B1(2),iX_B1(3),iGF_Psi), nDOFX )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4)
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(4)
#endif
    DO iX3    = iX_B1(3), iX_E1(3)
    DO iX2    = iX_B1(2), iX_E1(2)
    DO iX1    = iX_B1(1), iX_E1(1)
    DO iNodeX = 1       , nDOFX

      CALL ComputeGeometryX_SpatialMetric &
             ( G(iNodeX,iX1,iX2,iX3,iGF_h_1), &
               G(iNodeX,iX1,iX2,iX3,iGF_h_2), &
               G(iNodeX,iX1,iX2,iX3,iGF_h_3), &
               G(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11), &
               G(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22), &
               G(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33), &
               G(iNodeX,iX1,iX2,iX3,iGF_SqrtGm) )

      G(iNodeX,iX1,iX2,iX3,iGF_Beta_1) = Zero
      G(iNodeX,iX1,iX2,iX3,iGF_Beta_2) = Zero
      G(iNodeX,iX1,iX2,iX3,iGF_Beta_3) = Zero

    END DO
    END DO
    END DO
    END DO

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from: G ) &
    !$OMP MAP( release: iX_B1, iX_E1, CenterX1, CenterX2, WidthX1, WidthX2, &
    !$OMP               h_1_L, h_2_L, h_3_L, Alpha_L, Psi_L )
#elif defined( THORNADO_OACC   )
    !$ACC EXIT DATA &
    !$ACC COPYOUT( G ) &
    !$ACC DELETE( iX_B1, iX_E1, CenterX1, CenterX2, WidthX1, WidthX2, &
    !$ACC         h_1_L, h_2_L, h_3_L, Alpha_L, Psi_L )
#endif

    END ASSOCIATE

  END SUBROUTINE ComputeGeometryX_SPHERICAL


  SUBROUTINE ComputeGeometryX_CYLINDRICAL &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, Mass, MeshX )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout) :: &
      G(nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),nGF)
    REAL(DP), INTENT(in)    :: &
      Mass
    TYPE(MeshType), INTENT(in) :: &
      MeshX(3)

    INTEGER  :: iX1, iX2, iX3, iNodeX
    INTEGER  :: nP_X, nX(3)
    REAL(DP) :: x1L_q, x1G_q
    REAL(DP) :: h_1_L  (nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3))
    REAL(DP) :: h_2_L  (nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3))
    REAL(DP) :: h_3_L  (nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3))
    REAL(DP) :: Alpha_L(nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3))
    REAL(DP) :: Psi_L  (nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3))

    nX   = iX_E1 - iX_B1 + 1
    nP_X = PRODUCT( nX )

    ASSOCIATE ( CenterX1 => MeshX(1) % Center, &
                WidthX1  => MeshX(1) % Width )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: G, iX_B1, iX_E1, CenterX1, WidthX1 ) &
    !$OMP MAP( alloc: h_1_L, h_2_L, h_3_L, Alpha_L, Psi_L )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA &
    !$ACC COPYIN( G, iX_B1, iX_E1, CenterX1, WidthX1 ) &
    !$ACC CREATE( h_1_L, h_2_L, h_3_L, Alpha_L, Psi_L )
#endif

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4) &
    !$OMP PRIVATE( x1L_q, x1G_q )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRIVATE( x1L_q, x1G_q )
#elif defined( THORNADO_OMP   )
    !$OMP PARALLEL DO COLLAPSE(4) &
    !$OMP PRIVATE( x1L_q, x1G_q )
#endif
    DO iX3    = iX_B1(3), iX_E1(3)
    DO iX2    = iX_B1(2), iX_E1(2)
    DO iX1    = iX_B1(1), iX_E1(1)
    DO iNodeX = 1       , nDOFX

        ! --- Local Coordinates (Lobatto Points) ---

        x1L_q = NodesLX_q(1,iNodeX)

        ! --- Global Coordinates ---

        x1G_q = CenterX1(iX1) + WidthX1(iX1) * x1L_q

        ! --- Compute Lapse Function and Conformal Factor ---

        Alpha_L(iNodeX,iX1,iX2,iX3) &
          = One
        Psi_L(iNodeX,iX1,iX2,iX3) &
          = One

        ! --- Set Geometry in Lobatto Points ---

        h_1_L(iNodeX,iX1,iX2,iX3) &
          = Psi_L(iNodeX,iX1,iX2,iX3)**2
        h_2_L(iNodeX,iX1,iX2,iX3) &
          = Psi_L(iNodeX,iX1,iX2,iX3)**2
        h_3_L(iNodeX,iX1,iX2,iX3) &
          = Psi_L(iNodeX,iX1,iX2,iX3)**2 * x1G_q

    END DO
    END DO
    END DO
    END DO

    ! --- Interpolate from Lobatto to Gaussian Points ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX, nP_X, nDOFX, One, LX_L2G, nDOFX, &
             h_1_L, nDOFX, Zero, &
             G(1,iX_B1(1),iX_B1(2),iX_B1(3),iGF_h_1), nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX, nP_X, nDOFX, One, LX_L2G, nDOFX, &
             h_2_L, nDOFX, Zero, &
             G(1,iX_B1(1),iX_B1(2),iX_B1(3),iGF_h_2), nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX, nP_X, nDOFX, One, LX_L2G, nDOFX, &
             h_3_L, nDOFX, Zero, &
             G(1,iX_B1(1),iX_B1(2),iX_B1(3),iGF_h_3), nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX, nP_X, nDOFX, One, LX_L2G, nDOFX, &
             Alpha_L, nDOFX, Zero, &
             G(1,iX_B1(1),iX_B1(2),iX_B1(3),iGF_Alpha), nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX, nP_X, nDOFX, One, LX_L2G, nDOFX, &
             Psi_L, nDOFX, Zero, &
             G(1,iX_B1(1),iX_B1(2),iX_B1(3),iGF_Psi), nDOFX )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4)
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(4)
#endif
    DO iX3    = iX_B1(3), iX_E1(3)
    DO iX2    = iX_B1(2), iX_E1(2)
    DO iX1    = iX_B1(1), iX_E1(1)
    DO iNodeX = 1       , nDOFX

      CALL ComputeGeometryX_SpatialMetric &
             ( G(iNodeX,iX1,iX2,iX3,iGF_h_1), &
               G(iNodeX,iX1,iX2,iX3,iGF_h_2), &
               G(iNodeX,iX1,iX2,iX3,iGF_h_3), &
               G(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11), &
               G(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22), &
               G(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33), &
               G(iNodeX,iX1,iX2,iX3,iGF_SqrtGm) )

      G(iNodeX,iX1,iX2,iX3,iGF_Beta_1) = Zero
      G(iNodeX,iX1,iX2,iX3,iGF_Beta_2) = Zero
      G(iNodeX,iX1,iX2,iX3,iGF_Beta_3) = Zero

    END DO
    END DO
    END DO
    END DO

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from:    G ) &
    !$OMP MAP( release: iX_B1, iX_E1, CenterX1, WidthX1, &
    !$OMP               h_1_L, h_2_L, h_3_L, Alpha_L, Psi_L )
#elif defined( THORNADO_OACC   )
    !$ACC EXIT DATA &
    !$ACC COPYOUT(      G ) &
    !$ACC DELETE(       iX_B1, iX_E1, CenterX1, WidthX1, &
    !$ACC               h_1_L, h_2_L, h_3_L, Alpha_L, Psi_L )
#endif

    END ASSOCIATE

  END SUBROUTINE ComputeGeometryX_CYLINDRICAL


  SUBROUTINE ComputeGeometryX_SpatialMetric &
      ( h_1, h_2, h_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, SqrtGm )
#if   defined( THORNADO_OMP_OL )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC   )
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: h_1, h_2, h_3
    REAL(DP), INTENT(out) :: Gm_dd_11, Gm_dd_22, Gm_dd_33, SqrtGm

    Gm_dd_11 = MAX( h_1**2, SqrtTiny )
    Gm_dd_22 = MAX( h_2**2, SqrtTiny )
    Gm_dd_33 = MAX( h_3**2, SqrtTiny )

    SqrtGm = SQRT( Gm_dd_11 * Gm_dd_22 * Gm_dd_33 )

  END SUBROUTINE ComputeGeometryX_SpatialMetric


  SUBROUTINE ComputeGeometryX_FromScaleFactors( G )
#if   defined( THORNADO_OMP_OL )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC   )
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(inout) :: G(1:,1:)

    G(:,iGF_Gm_dd_11) = MAX( G(:,iGF_h_1)**2, SqrtTiny )
    G(:,iGF_Gm_dd_22) = MAX( G(:,iGF_h_2)**2, SqrtTiny )
    G(:,iGF_Gm_dd_33) = MAX( G(:,iGF_h_3)**2, SqrtTiny )

    G(:,iGF_SqrtGm) = G(:,iGF_h_1) * G(:,iGF_h_2) * G(:,iGF_h_3)

  END SUBROUTINE ComputeGeometryX_FromScaleFactors


  PURE REAL(DP) FUNCTION LapseFunction( R, M )
#if   defined( THORNADO_OMP_OL )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC   )
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: R, M

    ! --- Schwarzschild Metric in Isotropic Coordinates ---

    LapseFunction = ABS( ( MAX( ABS( R ), SqrtTiny ) - Half * M ) &
                       / ( MAX( ABS( R ), SqrtTiny ) + Half * M ) )

    RETURN
  END FUNCTION LapseFunction


  PURE REAL(DP) FUNCTION ConformalFactor( R, M )
#if   defined( THORNADO_OMP_OL )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC   )
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: R, M

    ! --- Schwarzschild Metric in Isotropic Coordinates ---

    ConformalFactor = One + Half * M / MAX( ABS( R ), SqrtTiny )

    RETURN
  END FUNCTION ConformalFactor


END MODULE GeometryComputationModule
