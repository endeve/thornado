MODULE XCFC_UtilitiesModule

  ! --- AMReX Modules ---

  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy, &
    amrex_imultifab
  USE amrex_amrcore_module, ONLY: &
    amrex_geom
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_reduce_sum

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    nNodesX
  USE ReferenceElementModuleX, ONLY: &
    WeightsX_q, &
    NodeNumberTableX
  USE UtilitiesModule, ONLY: &
    NodeNumberX
  USE LinearAlgebraModule, ONLY: &
    MatrixMatrixMultiply
  USE ReferenceElementModuleX, ONLY: &
    nDOFX_X1
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    LX_X1_Up
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
    iGF_Alpha, &
    iGF_Beta_1, &
    iGF_Beta_2, &
    iGF_Beta_3, &
    iGF_K_dd_11, &
    iGF_K_dd_12, &
    iGF_K_dd_13, &
    iGF_K_dd_22, &
    iGF_K_dd_23, &
    iGF_K_dd_33, &
    nGF

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero, &
    SqrtTiny, &
    One, &
    Two, &
    Pi
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF
  USE MaskModule, ONLY: &
    CreateFineMask, &
    DestroyFineMask, &
    IsNotLeafElement
  USE InputParsingModule, ONLY: &
    nLevels, &
    UseTiling

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MultiplyWithPsi6_MF
  PUBLIC :: UpdateConformalFactorAndMetric_MF
  PUBLIC :: UpdateGeometry_MF
  PUBLIC :: ComputeGravitationalMass_MF
  PUBLIC :: PopulateMF_uMF
  PUBLIC :: ApplyBoundaryConditions_Geometry

  INTEGER, PUBLIC :: swXX(3)

  ! --- MF: Metric Fields ---

  INTEGER, PARAMETER, PUBLIC :: iMF_Psi     = 1
  INTEGER, PARAMETER, PUBLIC :: iMF_Alpha   = 2
  INTEGER, PARAMETER, PUBLIC :: iMF_Beta_1  = 3
  INTEGER, PARAMETER, PUBLIC :: iMF_Beta_2  = 4
  INTEGER, PARAMETER, PUBLIC :: iMF_Beta_3  = 5
  INTEGER, PARAMETER, PUBLIC :: iMF_K_dd_11 = 6
  INTEGER, PARAMETER, PUBLIC :: iMF_K_dd_12 = 7
  INTEGER, PARAMETER, PUBLIC :: iMF_K_dd_13 = 8
  INTEGER, PARAMETER, PUBLIC :: iMF_K_dd_22 = 9
  INTEGER, PARAMETER, PUBLIC :: iMF_K_dd_23 = 10
  INTEGER, PARAMETER, PUBLIC :: iMF_K_dd_33 = 11
  INTEGER, PARAMETER, PUBLIC :: nMF         = 11

  ! --- GS: Gravity/Geometry Sources ---

  INTEGER, PARAMETER, PUBLIC :: iGS_E  = 1
  INTEGER, PARAMETER, PUBLIC :: iGS_S1 = 2
  INTEGER, PARAMETER, PUBLIC :: iGS_S2 = 3
  INTEGER, PARAMETER, PUBLIC :: iGS_S3 = 4
  INTEGER, PARAMETER, PUBLIC :: iGS_S  = 5
  INTEGER, PARAMETER, PUBLIC :: iGS_Mg = 6
  INTEGER, PARAMETER, PUBLIC :: nGS    = 6

CONTAINS


  SUBROUTINE MultiplyWithPsi6_MF &
    ( MF_uGF, Power, nDOFE, iE_B0, iE_E0, nS, MF_U )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:nLevels-1)
    INTEGER             , INTENT(in)    :: Power
    INTEGER             , INTENT(in)    :: nDOFE, iE_B0, iE_E0, nS
    TYPE(amrex_multifab), INTENT(inout) :: MF_U  (0:nLevels-1)

    TYPE(amrex_imultifab) :: iMF_FineMask
    TYPE(amrex_box)       :: BX
    TYPE(amrex_mfiter)    :: MFI

    REAL(DP), CONTIGUOUS, POINTER :: uGF     (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: U       (:,:,:,:)
    INTEGER , CONTIGUOUS, POINTER :: FineMask(:,:,:,:)

    INTEGER  :: iLevel, iNX, iNZ, iE, iX1, iX2, iX3, &
                iS, iFd, iComp, nDOF, nFd, nE
    INTEGER  :: iX_B0(3), iX_E0(3)

    REAL(DP) :: Psi6

    nE = iE_E0 - iE_B0 + 1

    nDOF = nDOFX * nDOFE

    nFd = MF_U(0) % nComp() / ( nDOF * nS * nE )

    DO iLevel = 0, nLevels-1

      CALL CreateFineMask( iLevel, iMF_FineMask, MF_uGF % BA, MF_uGF % DM )

#if defined( THORNADO_OMP )
      !$OMP PARALLEL &
      !$OMP PRIVATE( BX, MFI, uGF, U, FineMask, iNX, iComp, iX_B0, iX_E0, Psi6 )
#endif

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF      => MF_uGF(iLevel) % DataPtr( MFI )
        U        => MF_U  (iLevel) % DataPtr( MFI )
        FineMask => iMF_FineMask   % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi

        DO iS  = 1       , nS
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iE  = iE_B0   , iE_E0
        DO iNZ = 1       , nDOF

          IF( IsNotLeafElement( FineMask(iX1,iX2,iX3,1) ) ) CYCLE

          iNX  = MOD( ( iNZ - 1 ) / nDOFE, nDOFX ) + 1

          Psi6 = uGF(iX1,iX2,iX3,nDOFX*(iGF_Psi-1)+iNX)**6

          DO iFd = 1, nFd

            iComp = iNZ &
                      + ( iE  - 1 ) * nDOF &
                      + ( iFd - 1 ) * nDOF * nE &
                      + ( iS  - 1 ) * nDOF * nE * nFd

            U(iX1,iX2,iX3,iComp) &
              = U(iX1,iX2,iX3,iComp) * Psi6**( Power )

          END DO

        END DO
        END DO
        END DO
        END DO
        END DO
        END DO

      END DO ! WHILE( MFI % next() )

      CALL amrex_mfiter_destroy( MFI )

#if defined( THORNADO_OMP )
      !$OMP END PARALLEL
#endif

      CALL DestroyFineMask( iMF_FineMask )

    END DO ! iLevel = 0, nLevels-1

  END SUBROUTINE MultiplyWithPsi6_MF


  SUBROUTINE UpdateConformalFactorAndMetric_MF( MF_uMF, MF_uGF )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uMF(0:nLevels-1) ! Metric Fields
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:nLevels-1)

    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI

    REAL(DP), CONTIGUOUS, POINTER :: uMF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)

    INTEGER  :: iLevel, iNX, iX1, iX2, iX3, iNX1, iNX2
    INTEGER  :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP) :: X1, X2, Psi, h1, h2, h3

    DO iLevel = 0, nLevels-1

      CALL CreateMesh_MF( iLevel, MeshX )

#if defined( THORNADO_OMP )
      !$OMP PARALLEL &
      !$OMP PRIVATE( BX, MFI, uMF, uGF, iNX1, iNX2, &
      !$OMP          iX_B0, iX_E0, iX_B1, iX_E1, X1, X2, Psi, h1, h2, h3 )
#endif

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uMF => MF_uMF(iLevel) % DataPtr( MFI )
        uGF => MF_uGF(iLevel) % DataPtr( MFI )

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

          Psi = uMF(iX1,iX2,iX3,nDOFX*(iMF_Psi-1)+iNX)
          h1  = Psi**2
          h2  = Psi**2 * X1
          h3  = Psi**2 * X1 * SIN( X2 )

          uGF(iX1,iX2,iX3,nDOFX*(iGF_Psi-1)+iNX) = Psi
          uGF(iX1,iX2,iX3,nDOFX*(iGF_h_1-1)+iNX) = h1
          uGF(iX1,iX2,iX3,nDOFX*(iGF_h_2-1)+iNX) = h2
          uGF(iX1,iX2,iX3,nDOFX*(iGF_h_3-1)+iNX) = h3

          uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_11-1)+iNX) = MAX( h1**2, SqrtTiny )
          uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_22-1)+iNX) = MAX( h2**2, SqrtTiny )
          uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_33-1)+iNX) = MAX( h3**2, SqrtTiny )

          uGF(iX1,iX2,iX3,nDOFX*(iGF_SqrtGm-1)+iNX) = h1 * h2 * h3

        END DO
        END DO
        END DO
        END DO

      END DO ! WHILE( MFI % next() )

      CALL amrex_mfiter_destroy( MFI )

#if defined( THORNADO_OMP )
      !$OMP END PARALLEL
#endif

      CALL DestroyMesh_MF( MeshX )

    END DO ! iLevel = 0, nLevels-1

  END SUBROUTINE UpdateConformalFactorAndMetric_MF


  SUBROUTINE UpdateGeometry_MF( MF_uMF, MF_uGF )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uMF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:nLevels-1)

    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI

    REAL(DP), CONTIGUOUS, POINTER :: uMF (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uGF (:,:,:,:)

    INTEGER  :: iLevel, iNX, iX1, iX2, iX3
    INTEGER  :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)

    DO iLevel = 0, nLevels-1

#if defined( THORNADO_OMP )
      !$OMP PARALLEL &
      !$OMP PRIVATE( BX, MFI, uMF, uGF, iX_B0, iX_E0, iX_B1, iX_E1 )
#endif

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uMF => MF_uMF(iLevel) % DataPtr( MFI )
        uGF => MF_uGF(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = iX_B0 - swXX
        iX_E1 = iX_E0 + swXX

        DO iX3 = iX_B1(3), iX_E1(3)
        DO iX2 = iX_B1(2), iX_E1(2)
        DO iX1 = iX_B1(1), iX_E1(1)
        DO iNX = 1       , nDOFX

          uGF    (iX1,iX2,iX3,nDOFX*(iGF_Alpha-1)+iNX) &
            = uMF(iX1,iX2,iX3,nDOFX*(iMF_Alpha-1)+iNX)

          uGF    (iX1,iX2,iX3,nDOFX*(iGF_Beta_1-1)+iNX) &
            = uMF(iX1,iX2,iX3,nDOFX*(iMF_Beta_1-1)+iNX)
          uGF    (iX1,iX2,iX3,nDOFX*(iGF_Beta_2-1)+iNX) &
            = uMF(iX1,iX2,iX3,nDOFX*(iMF_Beta_2-1)+iNX)
          uGF    (iX1,iX2,iX3,nDOFX*(iGF_Beta_3-1)+iNX) &
            = uMF(iX1,iX2,iX3,nDOFX*(iMF_Beta_3-1)+iNX)

          uGF    (iX1,iX2,iX3,nDOFX*(iGF_K_dd_11-1)+iNX) &
            = uMF(iX1,iX2,iX3,nDOFX*(iMF_K_dd_11-1)+iNX)
          uGF    (iX1,iX2,iX3,nDOFX*(iGF_K_dd_12-1)+iNX) &
            = uMF(iX1,iX2,iX3,nDOFX*(iMF_K_dd_12-1)+iNX)
          uGF    (iX1,iX2,iX3,nDOFX*(iGF_K_dd_13-1)+iNX) &
            = uMF(iX1,iX2,iX3,nDOFX*(iMF_K_dd_13-1)+iNX)
          uGF    (iX1,iX2,iX3,nDOFX*(iGF_K_dd_22-1)+iNX) &
            = uMF(iX1,iX2,iX3,nDOFX*(iMF_K_dd_22-1)+iNX)
          uGF    (iX1,iX2,iX3,nDOFX*(iGF_K_dd_23-1)+iNX) &
            = uMF(iX1,iX2,iX3,nDOFX*(iMF_K_dd_23-1)+iNX)
          uGF    (iX1,iX2,iX3,nDOFX*(iGF_K_dd_33-1)+iNX) &
            = uMF(iX1,iX2,iX3,nDOFX*(iMF_K_dd_33-1)+iNX)

        END DO
        END DO
        END DO
        END DO

      END DO ! WHILE( MFI % next() )

      CALL amrex_mfiter_destroy( MFI )

#if defined( THORNADO_OMP )
      !$OMP END PARALLEL
#endif

    END DO ! iLevel = 0, nLevels-1

  END SUBROUTINE UpdateGeometry_MF


  SUBROUTINE PopulateMF_uMF( MF_uGF, MF_uMF )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uMF(0:)

    TYPE(amrex_imultifab) :: iMF_FineMask
    TYPE(amrex_box)       :: BX
    TYPE(amrex_mfiter)    :: MFI

    REAL(DP), CONTIGUOUS, POINTER :: uGF     (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uMF     (:,:,:,:)
    INTEGER , CONTIGUOUS, POINTER :: FineMask(:,:,:,:)

    INTEGER  :: iLevel, iNX, iX1, iX2, iX3
    INTEGER  :: iX_B0(3), iX_E0(3)

    DO iLevel = 0, nLevels-1

      CALL CreateFineMask( iLevel, iMF_FineMask, MF_uGF % BA, MF_uGF % DM )

#if defined( THORNADO_OMP )
      !$OMP PARALLEL &
      !$OMP PRIVATE( BX, MFI, uGF, uMF, FineMask, iX_B0, iX_E0 )
#endif

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF      => MF_uGF(iLevel) % DataPtr( MFI )
        uMF      => MF_uMF(iLevel) % DataPtr( MFI )
        FineMask => iMF_FineMask   % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi

        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1       , nDOFX

          IF( IsNotLeafElement( FineMask(iX1,iX2,iX3,1) ) ) CYCLE

          uMF    (iX1,iX2,iX3,nDOFX*(iMF_Psi-1)+iNX) &
            = uGF(iX1,iX2,iX3,nDOFX*(iGF_Psi-1)+iNX)

          uMF    (iX1,iX2,iX3,nDOFX*(iMF_Alpha-1)+iNX) &
            = uGF(iX1,iX2,iX3,nDOFX*(iGF_Alpha-1)+iNX)

          uMF    (iX1,iX2,iX3,nDOFX*(iMF_Beta_1-1)+iNX) &
            = uGF(iX1,iX2,iX3,nDOFX*(iGF_Beta_1-1)+iNX)
          uMF    (iX1,iX2,iX3,nDOFX*(iMF_Beta_2-1)+iNX) &
            = uGF(iX1,iX2,iX3,nDOFX*(iGF_Beta_2-1)+iNX)
          uMF    (iX1,iX2,iX3,nDOFX*(iMF_Beta_3-1)+iNX) &
            = uGF(iX1,iX2,iX3,nDOFX*(iGF_Beta_3-1)+iNX)

          uMF    (iX1,iX2,iX3,nDOFX*(iMF_K_dd_11-1)+iNX) &
            = uGF(iX1,iX2,iX3,nDOFX*(iGF_K_dd_11-1)+iNX)
          uMF    (iX1,iX2,iX3,nDOFX*(iMF_K_dd_12-1)+iNX) &
            = uGF(iX1,iX2,iX3,nDOFX*(iGF_K_dd_12-1)+iNX)
          uMF    (iX1,iX2,iX3,nDOFX*(iMF_K_dd_13-1)+iNX) &
            = uGF(iX1,iX2,iX3,nDOFX*(iGF_K_dd_13-1)+iNX)
          uMF    (iX1,iX2,iX3,nDOFX*(iMF_K_dd_22-1)+iNX) &
            = uGF(iX1,iX2,iX3,nDOFX*(iGF_K_dd_22-1)+iNX)
          uMF    (iX1,iX2,iX3,nDOFX*(iMF_K_dd_23-1)+iNX) &
            = uGF(iX1,iX2,iX3,nDOFX*(iGF_K_dd_23-1)+iNX)
          uMF    (iX1,iX2,iX3,nDOFX*(iMF_K_dd_33-1)+iNX) &
            = uGF(iX1,iX2,iX3,nDOFX*(iGF_K_dd_33-1)+iNX)

        END DO
        END DO
        END DO
        END DO

      END DO ! WHILE( MFI % next() )

      CALL amrex_mfiter_destroy( MFI )

#if defined( THORNADO_OMP )
      !$OMP END PARALLEL
#endif

    END DO ! iLevel = 0, nLevels-1

  END SUBROUTINE PopulateMF_uMF


  SUBROUTINE ComputeGravitationalMass_MF( MF_uGS, GravitationalMass )

    TYPE(amrex_multifab), INTENT(in)  :: MF_uGS(0:nLevels-1)
    REAL(DP)            , INTENT(out) :: GravitationalMass

    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI

    REAL(DP), CONTIGUOUS, POINTER :: uGS(:,:,:,:)

    INTEGER  :: iNX, iX1, iX2, iX3, iX_B0(3), iX_E0(3), iLevel
    REAL(DP) :: d3X, GravitationalMass_OMP

    ! --- Assuming 1D spherical symmetry ---

    GravitationalMass = Zero

    DO iLevel = 0, nLevels-1

      CALL CreateMesh_MF( iLevel, MeshX )

      GravitationalMass_OMP = Zero

#if defined( THORNADO_OMP )
      !$OMP PARALLEL &
      !$OMP PRIVATE( BX, MFI, uGS, iX_B0, iX_E0, d3X ) &
      !$OMP REDUCTION( +:GravitationalMass_OMP )
#endif

      CALL amrex_mfiter_build( MFI, MF_uGS(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGS => MF_uGS(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi

        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1       , nDOFX

          d3X = Two / Pi * MeshX(1) % Width(iX1) &
                         * MeshX(2) % Width(iX2) &
                         * MeshX(3) % Width(iX3)

          GravitationalMass_OMP &
            = GravitationalMass_OMP + d3X &
                * WeightsX_q(iNX) * uGS(iX1,iX2,iX3,nDOFX*(iGS_Mg-1)+iNX)

        END DO
        END DO
        END DO
        END DO

      END DO ! WHILE( MFI % next() )

      CALL amrex_mfiter_destroy( MFI )

#if defined( THORNADO_OMP )
      !$OMP END PARALLEL
#endif

      GravitationalMass = GravitationalMass_OMP

      CALL DestroyMesh_MF( MeshX )

    END DO ! iLevel = 0, nLevels-1

    CALL amrex_parallel_reduce_sum( GravitationalMass )

  END SUBROUTINE ComputeGravitationalMass_MF


  SUBROUTINE ApplyBoundaryConditions_Geometry( MF_uGF )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:nLevels-1)

    CALL ApplyBoundaryConditions_X1_Inner( MF_uGF )
    CALL ApplyBoundaryConditions_X1_Outer( MF_uGF )

  END SUBROUTINE ApplyBoundaryConditions_Geometry


  ! --- PRIVATE SUBROUTINES ---


  SUBROUTINE ApplyBoundaryConditions_X1_Inner( MF_uGF )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:nLevels-1)

    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)

    INTEGER :: iLevel, iX2, iX3
    INTEGER :: iX_B0(3), iX_E0(3)
    INTEGER :: iNX1, iNX2, iNX3, iNX
    INTEGER :: jNX1, jNX

    DO iLevel = 0, nLevels-1

#if defined( THORNADO_OMP )
      !$OMP PARALLEL &
      !$OMP PRIVATE( BX, MFI, uGF, iX_B0, iX_E0, iNX, jNX, jNX1 )
#endif

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi

        ! --- Inner Boundary: Reflecting ---

        IF( iX_B0(1) .EQ. amrex_geom(iLevel) % domain % lo( 1 ) )THEN

          DO iX3 = iX_B0(3), iX_E0(3)
          DO iX2 = iX_B0(2), iX_E0(2)

            DO iNX3 = 1, nNodesX(3)
            DO iNX2 = 1, nNodesX(2)
            DO iNX1 = 1, nNodesX(1)

              jNX1 = ( nNodesX(1) - iNX1 ) + 1

              iNX = NodeNumberX( iNX1, iNX2, iNX3 )
              jNX = NodeNumberX( jNX1, iNX2, iNX3 )

              uGF     (iX_B0(1)-1,iX2,iX3,nDOFX*(iGF_Alpha-1)+iNX) &
                = +uGF(iX_B0(1)  ,iX2,iX3,nDOFX*(iGF_Alpha-1)+jNX)

              uGF     (iX_B0(1)-1,iX2,iX3,nDOFX*(iGF_Psi-1)+iNX) &
                = +uGF(iX_B0(1)  ,iX2,iX3,nDOFX*(iGF_Psi-1)+jNX)

              uGF     (iX_B0(1)-1,iX2,iX3,nDOFX*(iGF_Beta_1-1)+iNX) &
                = -uGF(iX_B0(1)  ,iX2,iX3,nDOFX*(iGF_Beta_1-1)+jNX)

              uGF     (iX_B0(1)-1,iX2,iX3,nDOFX*(iGF_Beta_2-1)+iNX) &
                = +uGF(iX_B0(1)  ,iX2,iX3,nDOFX*(iGF_Beta_2-1)+jNX)

              uGF     (iX_B0(1)-1,iX2,iX3,nDOFX*(iGF_Beta_3-1)+iNX) &
                = +uGF(iX_B0(1)  ,iX2,iX3,nDOFX*(iGF_Beta_3-1)+jNX)

              uGF     (iX_B0(1)-1,iX2,iX3,nDOFX*(iGF_h_1-1)+iNX) &
                = +uGF(iX_B0(1)  ,iX2,iX3,nDOFX*(iGF_h_1-1)+jNX)

              uGF     (iX_B0(1)-1,iX2,iX3,nDOFX*(iGF_h_2-1)+iNX) &
                = +uGF(iX_B0(1)  ,iX2,iX3,nDOFX*(iGF_h_2-1)+jNX)

              uGF     (iX_B0(1)-1,iX2,iX3,nDOFX*(iGF_h_3-1)+iNX) &
                = +uGF(iX_B0(1)  ,iX2,iX3,nDOFX*(iGF_h_3-1)+jNX)

              uGF(iX_B0(1)-1,iX2,iX3,nDOFX*(iGF_Gm_dd_11-1)+iNX) &
                = MAX( uGF(iX_B0(1)-1,iX2,iX3,nDOFX*(iGF_h_1-1)+iNX)**2, &
                       SqrtTiny )

              uGF(iX_B0(1)-1,iX2,iX3,nDOFX*(iGF_Gm_dd_22-1)+iNX) &
                = MAX( uGF(iX_B0(1)-1,iX2,iX3,nDOFX*(iGF_h_2-1)+iNX)**2, &
                       SqrtTiny )

              uGF(iX_B0(1)-1,iX2,iX3,nDOFX*(iGF_Gm_dd_33-1)+iNX) &
                = MAX( uGF(iX_B0(1)-1,iX2,iX3,nDOFX*(iGF_h_3-1)+iNX)**2, &
                       SqrtTiny )

              uGF(iX_B0(1)-1,iX2,iX3,nDOFX*(iGF_SqrtGm-1)+iNX) &
                =   uGF(iX_B0(1)-1,iX2,iX3,nDOFX*(iGF_h_1-1)+iNX) &
                  * uGF(iX_B0(1)-1,iX2,iX3,nDOFX*(iGF_h_2-1)+iNX) &
                  * uGF(iX_B0(1)-1,iX2,iX3,nDOFX*(iGF_h_3-1)+iNX)

            END DO
            END DO
            END DO

          END DO
          END DO

        END IF

      END DO ! WHILE( MFI % next() )

      CALL amrex_mfiter_destroy( MFI )

#if defined( THORNADO_OMP )
      !$OMP END PARALLEL
#endif

    END DO ! iLevel = 0, nLevels-1

  END SUBROUTINE ApplyBoundaryConditions_X1_Inner


  SUBROUTINE ApplyBoundaryConditions_X1_Outer( MF_uGF )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:nLevels-1)

    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI

    REAL(DP), ALLOCATABLE :: G_K(:,:,:,:)
    REAL(DP), ALLOCATABLE :: G_F(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)

    INTEGER :: iLevel, iX2, iX3, iGF, nX1_X
    INTEGER :: iX_B0(3), iX_E0(3)
    INTEGER :: iNX

    DO iLevel = 0, nLevels-1

#if defined( THORNADO_OMP )
      !$OMP PARALLEL &
      !$OMP PRIVATE( BX, MFI, G_K, G_F, uGF, iX_B0, iX_E0, nX1_X )
#endif

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi

        ! --- Outer Boundary ---

        IF( iX_E0(1) .EQ. amrex_geom(iLevel) % domain % hi( 1 ) )THEN

          nX1_X = ( iX_E0(3) - iX_B0(3) + 1 ) * ( iX_E0(2) - iX_B0(2) + 1 )

          ALLOCATE( G_K(1:nDOFX   ,iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),1:nGF) )
          ALLOCATE( G_F(1:nDOFX_X1,iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),1:nGF) )

          DO iGF = 1       , nGF
          DO iX3 = iX_B0(3), iX_E0(3)
          DO iX2 = iX_B0(2), iX_E0(2)
          DO iNX = 1       , nDOFX

            G_K(iNX,iX2,iX3,iGF) = uGF(iX_E0(1),iX2,iX3,nDOFX*(iGF-1)+iNX)

          END DO
          END DO
          END DO
          END DO

          DO iGF = 1, nGF

            CALL MatrixMatrixMultiply &
                   ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, One, LX_X1_Up, &
                     nDOFX_X1,   G_K(1,iX_B0(2),iX_B0(3),iGF), &
                     nDOFX, Zero,G_F(1,iX_B0(2),iX_B0(3),iGF), &
                     nDOFX_X1 )

          END DO

          DO iGF = 1       , nGF
          DO iX3 = iX_B0(3), iX_E0(3)
          DO iX2 = iX_B0(2), iX_E0(2)
          DO iNX = 1       , nDOFX

            uGF(iX_E0(1)+1,iX2,iX3,nDOFX*(iGF-1)+iNX) = G_F(1,iX2,iX3,iGF)

          END DO
          END DO
          END DO
          END DO

          DEALLOCATE( G_F )
          DEALLOCATE( G_K )

        END IF

      END DO ! WHILE( MFI % next() )

      CALL amrex_mfiter_destroy( MFI )

#if defined( THORNADO_OMP )
      !$OMP END PARALLEL
#endif

    END DO ! iLevel = 0, nLevels-1

  END SUBROUTINE ApplyBoundaryConditions_X1_Outer


END MODULE XCFC_UtilitiesModule
