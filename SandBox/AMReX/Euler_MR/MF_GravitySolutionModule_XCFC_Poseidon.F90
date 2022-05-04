MODULE MF_GravitySolutionModule_XCFC_Poseidon

  ! --- AMReX Modules ---

  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nDOFX
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
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
    iGF_Psi

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    SqrtTiny
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF
  USE InputParsingModule, ONLY: &
    nLevels, &
    UseTiling

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: UpdateConformalFactorAndMetric_MF

CONTAINS


  SUBROUTINE UpdateConformalFactorAndMetric_MF( MF_uMF, MF_uGF )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uMF(0:nLevels-1) ! Metric Fields
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:nLevels-1)

    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI

    REAL(DP), CONTIGUOUS, POINTER :: uMF (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uGF (:,:,:,:)

    INTEGER  :: iLevel, iNX, iX1, iX2, iX3, iNX1, iNX2
    INTEGER  :: iX_B0(3), iX_E0(3)
    INTEGER  :: iLo_G(4), iLo_M(4)
    REAL(DP) :: X1, X2, Psi, h1, h2, h3

    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      CALL CreateMesh_MF( iLevel, MeshX )

      DO WHILE( MFI % next() )

        uMF => MF_uMF(iLevel) % DataPtr( MFI )
        uGF => MF_uGF(iLevel) % DataPtr( MFI )

        iLo_M = LBOUND( uMF )
        iLo_G = LBOUND( uGF )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi

        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)

          DO iNX = 1, nDOFX

            iNX1 = NodeNumberTableX(1,iNX)
            iNX2 = NodeNumberTableX(2,iNX)

            X1 = NodeCoordinate( MeshX(1), iX1, iNX1 )
            X2 = NodeCoordinate( MeshX(2), iX2, iNX2 )

            Psi = uMF(iX1,iX2,iX3,nDOFX*(1-1)+iNX) ! 1st component of uMF is Psi
            h1 = Psi**2
            h2 = Psi**2 * X1
            h3 = Psi**2 * X1 * SIN( X2 )

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

      END DO

      CALL DestroyMesh_MF( MeshX )

      CALL amrex_mfiter_destroy( MFI )

    END DO

  END SUBROUTINE UpdateConformalFactorAndMetric_MF


END MODULE MF_GravitySolutionModule_XCFC_Poseidon
