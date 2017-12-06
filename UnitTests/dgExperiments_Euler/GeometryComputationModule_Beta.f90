MODULE GeometryComputationModule_Beta

  USE KindModule, ONLY: &
    DP, Zero, One
  USE ProgramHeaderModule, ONLY: &
    nDOFX
  USE ReferenceElementModuleX, ONLY: &
    NodesLX_q
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    LX_L2G
  USE MeshModule, ONLY: &
    MeshX
  USE GeometryFieldsModule, ONLY: &
    CoordinateSystem, &
    iGF_h_1, &
    iGF_h_2, &
    iGF_h_3, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_SqrtGm, &
    nGF

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeGeometryX
  PUBLIC :: ComputeGeometryX_FromScaleFactors

CONTAINS


  SUBROUTINE ComputeGeometryX( iX_B0, iX_E0, iX_B1, iX_E1, G )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout) :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    SELECT CASE ( TRIM( CoordinateSystem ) )

      CASE ( 'CARTESIAN' )

        CALL ComputeGeometryX_CARTESIAN &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G )

      CASE ( 'SPHERICAL' )

        CALL ComputeGeometryX_SPHERICAL &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G )

      CASE ( 'CYLINDRICAL' )

        CALL ComputeGeometryX_CYLINDRICAL &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G )

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,'(A5,A27,A)') &
          '', 'Invalid Coordinate System: ', TRIM( CoordinateSystem )
        STOP

    END SELECT

    ! Set Boundary Conditions for Geometry Fields Here
    PRINT*
    PRINT*, "Boundary Conditions Not Set in ComputeGeometryX"
    PRINT*

  END SUBROUTINE ComputeGeometryX


  SUBROUTINE ComputeGeometryX_CARTESIAN( iX_B0, iX_E0, iX_B1, iX_E1, G )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout) :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER :: iX1, iX2, iX3

    DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)

          G(:,iX1,iX2,iX3,iGF_h_1) = One
          G(:,iX1,iX2,iX3,iGF_h_2) = One
          G(:,iX1,iX2,iX3,iGF_h_3) = One

          CALL ComputeGeometryX_FromScaleFactors( G(:,iX1,iX2,iX3,:) )

        END DO
      END DO
    END DO

  END SUBROUTINE ComputeGeometryX_CARTESIAN


  SUBROUTINE ComputeGeometryX_SPHERICAL( iX_B0, iX_E0, iX_B1, iX_E1, G )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout) :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    REAL(DP) :: XC(3), dX(3), xL_q(3), xG_q(3)
    REAL(DP) :: G_L(nDOFX,nGF)
    INTEGER  :: iX1, iX2, iX3, iNodeX

    DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)

        XC(2) = MeshX(2) % Center(iX2)
        dX(2) = MeshX(2) % Width (iX2)

        DO iX1 = iX_B0(1), iX_E0(1)

          XC(1) = MeshX(1) % Center(iX1)
          dX(1) = MeshX(1) % Width (iX1)

          ! --- Compute Geometry Fields in Lobatto Points ---

          DO iNodeX = 1, nDOFX

            ! --- Local Coordinates (Lobatto Points) ---

            xL_q = NodesLX_q(1:3,iNodeX)

            ! --- Global Coordinates ---

            xG_q = XC + dX * xL_q

            ! --- Set Geometry in Lobatto Points ---

            G_L(iNodeX,iGF_h_1) = One
            G_L(iNodeX,iGF_h_2) = xG_q(1)
            G_L(iNodeX,iGF_h_3) = xG_q(1) * MAX( SIN( xG_q(2) ), Zero )

          END DO

          ! --- Interpolate from Lobatto to Gaussian Points ---

          CALL DGEMV &
                 ( 'N', nDOFX, nDOFX, One, LX_L2G, nDOFX, &
                   G_L(:,iGF_h_1), 1, Zero, G(:,iX1,iX2,iX3,iGF_h_1), 1 )
          CALL DGEMV &
                 ( 'N', nDOFX, nDOFX, One, LX_L2G, nDOFX, &
                   G_L(:,iGF_h_2), 1, Zero, G(:,iX1,iX2,iX3,iGF_h_2), 1 )
          CALL DGEMV &
                 ( 'N', nDOFX, nDOFX, One, LX_L2G, nDOFX, &
                   G_L(:,iGF_h_3), 1, Zero, G(:,iX1,iX2,iX3,iGF_h_3), 1 )

          CALL ComputeGeometryX_FromScaleFactors( G(:,iX1,iX2,iX3,:) )

        END DO
      END DO
    END DO

  END SUBROUTINE ComputeGeometryX_SPHERICAL


  SUBROUTINE ComputeGeometryX_CYLINDRICAL( iX_B0, iX_E0, iX_B1, iX_E1, G )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout) :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    REAL(DP) :: XC(3), dX(3), xL_q(3), xG_q(3)
    REAL(DP) :: G_L(nDOFX,nGF)
    INTEGER  :: iX1, iX2, iX3, iNodeX

    PRINT*, "ComputeGeometryX_CYLINDRICAL"

    DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)

          XC(1) = MeshX(1) % Center(iX1)
          dX(1) = MeshX(1) % Width (iX1)

          DO iNodeX = 1, nDOFX

            ! --- Local Coordinates (Lobatto Points) ---

            xL_q = NodesLX_q(1:3,iNodeX)

            ! --- Global Coordinates ---

            xG_q = XC + dX * xL_q

            ! --- Set Geometry in Lobatto Points ---

            G_L(iNodeX,iGF_h_1) = One
            G_L(iNodeX,iGF_h_2) = One
            G_L(iNodeX,iGF_h_3) = xG_q(1)

          END DO

          ! --- Interpolate from Lobatto to Gaussian Points ---

          CALL DGEMV &
                 ( 'N', nDOFX, nDOFX, One, LX_L2G, nDOFX, &
                   G_L(:,iGF_h_1), 1, Zero, G(:,iX1,iX2,iX3,iGF_h_1), 1 )
          CALL DGEMV &
                 ( 'N', nDOFX, nDOFX, One, LX_L2G, nDOFX, &
                   G_L(:,iGF_h_2), 1, Zero, G(:,iX1,iX2,iX3,iGF_h_2), 1 )
          CALL DGEMV &
                 ( 'N', nDOFX, nDOFX, One, LX_L2G, nDOFX, &
                   G_L(:,iGF_h_3), 1, Zero, G(:,iX1,iX2,iX3,iGF_h_3), 1 )

          CALL ComputeGeometryX_FromScaleFactors( G(:,iX1,iX2,iX3,:) )

        END DO
      END DO
    END DO

  END SUBROUTINE ComputeGeometryX_CYLINDRICAL


  SUBROUTINE ComputeGeometryX_FromScaleFactors( G )

    REAL(DP), INTENT(inout) :: G(1:,1:)

    G(:,iGF_Gm_dd_11) = G(:,iGF_h_1)**2
    G(:,iGF_Gm_dd_22) = G(:,iGF_h_2)**2
    G(:,iGF_Gm_dd_33) = G(:,iGF_h_3)**2

    G(:,iGF_SqrtGm) = G(:,iGF_h_1) * G(:,iGF_h_2) * G(:,iGF_h_3)

  END SUBROUTINE ComputeGeometryX_FromScaleFactors


END MODULE GeometryComputationModule_Beta
