MODULE MF_TwoMoment_XCFC_UtilitiesModule

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
    nDOFX, &
    nDOFZ, &
    iE_E0, &
    iE_B0, &
    nDOFE, &
    nE
  USE MeshModule, ONLY: &
    MeshX, &
    MeshE
  USE GeometryFieldsModule, ONLY: &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_SqrtGm, &
    iGF_Alpha, &
    iGF_Beta_1, &
    iGF_Beta_2, &
    iGF_Beta_3, &
    iGF_Psi
  USE GeometryFieldsModuleE, ONLY: &
    iGE_Ep3, &
    uGE
  USE FluidFieldsModule, ONLY: &
    iCF_D, &
    iCF_S1, &
    iCF_S2, &
    iCF_S3, &
    iCF_E, &
    iCF_Ne, &
    iPF_D, &
    iPF_V1, &
    iPF_V2, &
    iPF_V3, &
    iPF_E, &
    iPF_Ne, &
    nPF
  USE XCFC_UtilitiesModule, ONLY: &
    iGS_E, &
    iGS_S1, &
    iGS_S2, &
    iGS_S3, &
    iGS_S, &
    iGS_Mg
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    nCR,    &
    iCR_N,  &
    iCR_G1, &
    iCR_G2, &
    iCR_G3
  USE Euler_UtilitiesModule_Relativistic, ONLY: &
    ComputePrimitive_Euler_Relativistic
  USE EquationOfStateModule, ONLY: &
    ComputePressureFromPrimitive
  USE Euler_ErrorModule, ONLY: &
    DescribeError_Euler
  USE ReferenceElementModuleE, ONLY: &
    WeightsE

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero, &
    One, &
    Two, &
    Three, &
    FourPi
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF
  USE InputParsingModule, ONLY: &
    nLevels, &
    UseTiling

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeConformalFactorSourcesAndMg_XCFC_TwoMoment_MF
  PUBLIC :: ComputePressureTensorTrace_XCFC_TwoMoment_MF

CONTAINS


  SUBROUTINE ComputeConformalFactorSourcesAndMg_XCFC_TwoMoment_MF &
    ( MF_uGF, MF_uCF, MF_uCR, MF_uGS )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in)    :: MF_uCF(0:nLevels-1) ! Psi^6 * U
    TYPE(amrex_multifab), INTENT(in)    :: MF_uCR(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGS(0:nLevels-1)

    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCR(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uGS(:,:,:,:)

    INTEGER  :: iLevel, iNX, iX1, iX2, iX3, ErrorExists
    INTEGER  :: iX_B0(3), iX_E0(3)
    REAL(DP) :: Psi6
    REAL(DP) :: uPF(nPF), LorentzFactor, BetaDotV, Enthalpy, Pressure

    INTEGER, ALLOCATABLE :: ITERATION(:,:,:,:)
    INTEGER, ALLOCATABLE :: iErr     (:,:,:,:)

    REAL(DP) :: E, S_i(3), E_int, S_i_int(3)
    REAL(DP) :: N, G_d_1, G_d_2, G_d_3, vG
    REAL(DP) :: V_u_1, V_u_2, V_u_3
    REAL(DP) :: V_d_1, V_d_2, V_d_3
    REAL(DP) :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    INTEGER  :: iD_N, iD_G1, iD_G2, iD_G3, iE, iN_E, iS, iN_Z

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

    ASSOCIATE &
      ( dZ1 => MeshE  % Width )

    E       = Zero
    E_int   = Zero
    S_i     = Zero
    S_i_int = Zero

    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCF => MF_uCF(iLevel) % DataPtr( MFI )
        uCR => MF_uCR(iLevel) % DataPtr( MFI )
        uGS => MF_uGS(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi

        ErrorExists = 0

        ALLOCATE( ITERATION(1:nDOFX,iX_B0(1):iX_E0(1), &
                                    iX_B0(2):iX_E0(2), &
                                    iX_B0(3):iX_E0(3)) )
        ALLOCATE( iErr     (1:nDOFX,iX_B0(1):iX_E0(1), &
                                    iX_B0(2):iX_E0(2), &
                                    iX_B0(3):iX_E0(3)) )

        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1       , nDOFX

          uGS       (iX1,iX2,iX3,nDOFX*(iGS_E-1)+iNX) &
            =  ( uCF(iX1,iX2,iX3,nDOFX*(iCF_E-1)+iNX) &
               + uCF(iX1,iX2,iX3,nDOFX*(iCF_D-1)+iNX) )

          uGS    (iX1,iX2,iX3,nDOFX*(iGS_S1-1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_S1-1)+iNX)

          uGS    (iX1,iX2,iX3,nDOFX*(iGS_S2-1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_S2-1)+iNX)

          uGS    (iX1,iX2,iX3,nDOFX*(iGS_S3-1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_S3-1)+iNX)

          ! Assume Psi^(iStage) ~ Psi^(iStage+1) for Poseidon BCs

          Psi6 = uGF(iX1,iX2,iX3,nDOFX*(iGF_Psi-1)+iNX)**6

          ITERATION(iNX,iX1,iX2,iX3) = 0
          iErr     (iNX,iX1,iX2,iX3) = 0

          uCF    (iX1,iX2,iX3,nDOFX*(iCF_D -1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_D -1)+iNX) / Psi6
          uCF    (iX1,iX2,iX3,nDOFX*(iCF_S1-1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_S1-1)+iNX) / Psi6
          uCF    (iX1,iX2,iX3,nDOFX*(iCF_S2-1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_S2-1)+iNX) / Psi6
          uCF    (iX1,iX2,iX3,nDOFX*(iCF_S3-1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_S3-1)+iNX) / Psi6
          uCF    (iX1,iX2,iX3,nDOFX*(iCF_E -1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_E -1)+iNX) / Psi6
          uCF    (iX1,iX2,iX3,nDOFX*(iCF_Ne-1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_Ne-1)+iNX) / Psi6

          CALL ComputePrimitive_Euler_Relativistic &
                 ( uCF(iX1,iX2,iX3,nDOFX*(iCF_D -1)+iNX), &
                   uCF(iX1,iX2,iX3,nDOFX*(iCF_S1-1)+iNX), &
                   uCF(iX1,iX2,iX3,nDOFX*(iCF_S2-1)+iNX), &
                   uCF(iX1,iX2,iX3,nDOFX*(iCF_S3-1)+iNX), &
                   uCF(iX1,iX2,iX3,nDOFX*(iCF_E -1)+iNX), &
                   uCF(iX1,iX2,iX3,nDOFX*(iCF_Ne-1)+iNX), &
                   uPF(iPF_D ), &
                   uPF(iPF_V1), &
                   uPF(iPF_V2), &
                   uPF(iPF_V3), &
                   uPF(iPF_E ), &
                   uPF(iPF_Ne), &
                   uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_11-1)+iNX), &
                   uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_22-1)+iNX), &
                   uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_33-1)+iNX), &
                   ITERATION_Option = ITERATION(iNX,iX1,iX2,iX3), &
                   iErr_Option      = iErr     (iNX,iX1,iX2,iX3) )

          uCF    (iX1,iX2,iX3,nDOFX*(iCF_D -1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_D -1)+iNX) * Psi6
          uCF    (iX1,iX2,iX3,nDOFX*(iCF_S1-1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_S1-1)+iNX) * Psi6
          uCF    (iX1,iX2,iX3,nDOFX*(iCF_S2-1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_S2-1)+iNX) * Psi6
          uCF    (iX1,iX2,iX3,nDOFX*(iCF_S3-1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_S3-1)+iNX) * Psi6
          uCF    (iX1,iX2,iX3,nDOFX*(iCF_E -1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_E -1)+iNX) * Psi6
          uCF    (iX1,iX2,iX3,nDOFX*(iCF_Ne-1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_Ne-1)+iNX) * Psi6

          ErrorExists = ErrorExists + iErr(iNX,iX1,iX2,iX3)

          CALL ComputePressureFromPrimitive &
                 ( uPF(iPF_D), uPF(iPF_E), uPF(iPF_Ne), Pressure )

          LorentzFactor &
            = One / SQRT( One                              &
                - ( uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_11-1)+iNX) &
                      * uPF(iPF_V1)**2 &
                  + uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_22-1)+iNX) &
                      * uPF(iPF_V2)**2 &
                  + uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_33-1)+iNX) &
                      * uPF(iPF_V3)**2 ) )

          BetaDotV =   uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_11-1)+iNX) &
                         * uGF(iX1,iX2,iX3,nDOFX*(iGF_Beta_1-1)+iNX) &
                         * uPF(iPF_V1) &
                     + uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_22-1)+iNX) &
                         * uGF(iX1,iX2,iX3,nDOFX*(iGF_Beta_2-1)+iNX) &
                         * uPF(iPF_V2) &
                     + uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_33-1)+iNX) &
                         * uGF(iX1,iX2,iX3,nDOFX*(iGF_Beta_3-1)+iNX) &
                         * uPF(iPF_V3)

          Enthalpy = uPF(iPF_D) + uPF(iPF_E) + Pressure

          uGS(iX1,iX2,iX3,nDOFX*(iGS_Mg-1)+iNX) &
            = ( Enthalpy * ( Two * LorentzFactor**2 &
                  * ( One - BetaDotV &
                              / uGF(iX1,iX2,iX3,nDOFX*(iGF_Alpha-1)+iNX) ) &
                      - One ) &
                + Two * Pressure ) &
               * uGF(iX1,iX2,iX3,nDOFX*(iGF_Alpha -1)+iNX) &
               * uGF(iX1,iX2,iX3,nDOFX*(iGF_SqrtGm-1)+iNX)

          DO iS   = 1, nSpecies
          DO iE   = 1, nE
          DO iN_E = 1, nDOFE

            iN_Z = (iNX-1) * nDOFE + iN_E

            iD_N = ( iS - 1 ) * nCR * ( iE_E0 - iE_B0 + 1 ) * nDOFZ &
                 + ( iCR_N - 1 ) * ( iE_E0 - iE_B0 + 1 ) * nDOFZ &
                 + ( iE - 1 ) * nDOFZ + iN_Z
            iD_G1 = ( iS - 1 ) * nCR * ( iE_E0 - iE_B0 + 1 ) * nDOFZ &
                 + ( iCR_G1 - 1 ) * ( iE_E0 - iE_B0 + 1 ) * nDOFZ &
                 + ( iE - 1 ) * nDOFZ + iN_Z
            iD_G2 = ( iS - 1 ) * nCR * ( iE_E0 - iE_B0 + 1 ) * nDOFZ &
                 + ( iCR_G2 - 1 ) * ( iE_E0 - iE_B0 + 1 ) * nDOFZ &
                 + ( iE - 1 ) * nDOFZ + iN_Z
            iD_G3 = ( iS - 1 ) * nCR * ( iE_E0 - iE_B0 + 1 ) * nDOFZ &
                 + ( iCR_G3 - 1 ) * ( iE_E0 - iE_B0 + 1 ) * nDOFZ &
                 + ( iE - 1 ) * nDOFZ + iN_Z

            N     = uCR(iX1,iX2,iX3,iD_N)
            G_d_1 = uCR(iX1,iX2,iX3,iD_G1)
            G_d_2 = uCR(iX1,iX2,iX3,iD_G2)
            G_d_3 = uCR(iX1,iX2,iX3,iD_G3)

            V_u_1 = uPF (iPF_V1)
            V_u_2 = uPF (iPF_V2)
            V_u_3 = uPF (iPF_V3)

            Gm_dd_11 = uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_11-1)+iNX)
            Gm_dd_22 = uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_22-1)+iNX)
            Gm_dd_33 = uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_33-1)+iNX)

            V_d_1 = V_u_1 * Gm_dd_11
            V_d_2 = V_u_2 * Gm_dd_22
            V_d_3 = V_u_3 * Gm_dd_33

            vG = V_u_1 * G_d_1 + V_u_2 * G_d_2 + V_u_3 * G_d_3

            E_int      = LorentzFactor * N + vG
            S_i_int(1) = LorentzFactor * V_d_1 * N + G_d_1
            S_i_int(2) = LorentzFactor * V_d_2 * N + G_d_2
            S_i_int(3) = LorentzFactor * V_d_3 * N + G_d_3

            E = E &
              + FourPi * dZ1(iE) * WeightsE(iN_E) &
              * uGE(iN_E,iE,iGE_Ep3) * E_int
            S_i(1) = S_i(1) &
                   + FourPi * dZ1(iE) * WeightsE(iN_E) &
                   * uGE(iN_E,iE,iGE_Ep3) * S_i_int(1)
            S_i(2) = S_i(2) &
                   + FourPi * dZ1(iE) * WeightsE(iN_E) &
                   * uGE(iN_E,iE,iGE_Ep3) * S_i_int(2)
            S_i(3) = S_i(3) &
                   + FourPi * dZ1(iE) * WeightsE(iN_E) &
                   * uGE(iN_E,iE,iGE_Ep3) * S_i_int(3)

          END DO
          END DO
          END DO

          uGS(iX1,iX2,iX3,nDOFX*(iGS_E-1)+iNX) &
            = uGS(iX1,iX2,iX3,nDOFX*(iGS_E-1)+iNX) + E

          uGS(iX1,iX2,iX3,nDOFX*(iGS_S1-1)+iNX) &
            = uGS(iX1,iX2,iX3,nDOFX*(iGS_S1-1)+iNX) + S_i(1)

          uGS(iX1,iX2,iX3,nDOFX*(iGS_S2-1)+iNX) &
            = uGS(iX1,iX2,iX3,nDOFX*(iGS_S2-1)+iNX) + S_i(2)

          uGS(iX1,iX2,iX3,nDOFX*(iGS_S3-1)+iNX) &
            = uGS(iX1,iX2,iX3,nDOFX*(iGS_S3-1)+iNX) + S_i(3)

          E   = Zero
          S_i = Zero

        END DO
        END DO
        END DO
        END DO

        IF( ErrorExists .GT. 0 )THEN

          WRITE(*,*) 'ERROR'
          WRITE(*,*) '-----'
          WRITE(*,*) &
            '    MODULE: MF_TwoMoment_XCFC_UtilitiesModule'
          WRITE(*,*) &
            'SUBROUTINE: ComputeConformalFactorSourcesAndMg_XCFC_TwoMoment_MF'

          CALL CreateMesh_MF( iLevel, MeshX )

          DO iX3 = iX_B0(3), iX_E0(3)
          DO iX2 = iX_B0(2), iX_E0(2)
          DO iX1 = iX_B0(1), iX_E0(1)
          DO iNX = 1       , nDOFX

            Psi6 = uGF(iX1,iX2,iX3,nDOFX*(iGF_Psi-1)+iNX)**6

            CALL DescribeError_Euler &
                   ( iErr(iNX,iX1,iX2,iX3), &
                     Int_Option &
                       = [ ITERATION(iNX,iX1,iX2,iX3), 99999999, &
                           iX_B0(1), iX_B0(2), iX_B0(3), &
                           iX_E0(1), iX_E0(2), iX_E0(3), &
                           iNX, iX1, iX2, iX3 ], &
                     Real_Option &
                       = [ MeshX(1) % Center(iX1), &
                           MeshX(2) % Center(iX2), &
                           MeshX(3) % Center(iX3), &
                           MeshX(1) % Width (iX1), &
                           MeshX(2) % Width (iX2), &
                           MeshX(3) % Width (iX3), &
                           uCF(iX1,iX2,iX3,nDOFX*(iCF_D -1)+iNX) / Psi6, &
                           uCF(iX1,iX2,iX3,nDOFX*(iCF_S1-1)+iNX) / Psi6, &
                           uCF(iX1,iX2,iX3,nDOFX*(iCF_S2-1)+iNX) / Psi6, &
                           uCF(iX1,iX2,iX3,nDOFX*(iCF_S3-1)+iNX) / Psi6, &
                           uCF(iX1,iX2,iX3,nDOFX*(iCF_E -1)+iNX) / Psi6, &
                           uCF(iX1,iX2,iX3,nDOFX*(iCF_Ne-1)+iNX) / Psi6, &
                           uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_11-1)+iNX), &
                           uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_22-1)+iNX), &
                           uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_33-1)+iNX) ], &
                     Char_Option = [ 'NA' ], &
                     Message_Option &
                       = 'Calling from ComputeConformalFactorSourcesAndMg_XCFC_TwoMoment_MF' )

          END DO
          END DO
          END DO
          END DO

          CALL DestroyMesh_MF( MeshX )

        END IF

        DEALLOCATE( ITERATION )
        DEALLOCATE( iErr      )

      END DO ! WHILE( MFI % next() )

      CALL amrex_mfiter_destroy( MFI )

    END DO ! iLevel = 0, nLevels-1

    END ASSOCIATE

#endif

  END SUBROUTINE ComputeConformalFactorSourcesAndMg_XCFC_TwoMoment_MF


  SUBROUTINE ComputePressureTensorTrace_XCFC_TwoMoment_MF &
    ( MF_uGF, MF_uCF, MF_uCR, MF_uGS )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in)    :: MF_uCF(0:nLevels-1) ! Psi^6 * U
    TYPE(amrex_multifab), INTENT(in)    :: MF_uCR(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGS(0:nLevels-1)

    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCR(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uGS(:,:,:,:)

    INTEGER :: iLevel, iNX, iX1, iX2, iX3
    INTEGER :: iX_B0(3), iX_E0(3)
    INTEGER :: ErrorExists
    INTEGER, ALLOCATABLE :: ITERATION(:,:,:,:)
    INTEGER, ALLOCATABLE :: iErr     (:,:,:,:)

    REAL(DP) :: uPF(nPF), Pressure, Psi6

    REAL(DP) :: S, S_int, N, G_d_1, G_d_2, G_d_3, vG
    REAL(DP) :: LorentzFactor, V_u_1, V_u_2, V_u_3
    INTEGER  :: iD_N, iD_G1, iD_G2, iD_G3, iE, iN_E, iS, iN_Z

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

    ASSOCIATE &
      ( dZ1 => MeshE  % Width )

    S = Zero

    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCF => MF_uCF(iLevel) % DataPtr( MFI )
        uCR => MF_uCR(iLevel) % DataPtr( MFI )
        uGS => MF_uGS(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi

        ALLOCATE( ITERATION(1:nDOFX,iX_B0(1):iX_E0(1), &
                                    iX_B0(2):iX_E0(2), &
                                    iX_B0(3):iX_E0(3)) )
        ALLOCATE( iErr     (1:nDOFX,iX_B0(1):iX_E0(1), &
                                    iX_B0(2):iX_E0(2), &
                                    iX_B0(3):iX_E0(3)) )

        ErrorExists = 0

        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1       , nDOFX

          ITERATION(iNX,iX1,iX2,iX3) = 0
          iErr     (iNX,iX1,iX2,iX3) = 0

          Psi6 = uGF(iX1,iX2,iX3,nDOFX*(iGF_Psi-1)+iNX)**6

          ! --- Compute trace of stress tensor ---

          uCF    (iX1,iX2,iX3,nDOFX*(iCF_D -1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_D -1)+iNX) / Psi6
          uCF    (iX1,iX2,iX3,nDOFX*(iCF_S1-1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_S1-1)+iNX) / Psi6
          uCF    (iX1,iX2,iX3,nDOFX*(iCF_S2-1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_S2-1)+iNX) / Psi6
          uCF    (iX1,iX2,iX3,nDOFX*(iCF_S3-1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_S3-1)+iNX) / Psi6
          uCF    (iX1,iX2,iX3,nDOFX*(iCF_E -1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_E -1)+iNX) / Psi6
          uCF    (iX1,iX2,iX3,nDOFX*(iCF_Ne-1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_Ne-1)+iNX) / Psi6

          CALL ComputePrimitive_Euler_Relativistic &
                 ( uCF(iX1,iX2,iX3,nDOFX*(iCF_D -1)+iNX), &
                   uCF(iX1,iX2,iX3,nDOFX*(iCF_S1-1)+iNX), &
                   uCF(iX1,iX2,iX3,nDOFX*(iCF_S2-1)+iNX), &
                   uCF(iX1,iX2,iX3,nDOFX*(iCF_S3-1)+iNX), &
                   uCF(iX1,iX2,iX3,nDOFX*(iCF_E -1)+iNX), &
                   uCF(iX1,iX2,iX3,nDOFX*(iCF_Ne-1)+iNX), &
                   uPF(iPF_D ), &
                   uPF(iPF_V1), &
                   uPF(iPF_V2), &
                   uPF(iPF_V3), &
                   uPF(iPF_E ), &
                   uPF(iPF_Ne), &
                   uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_11-1)+iNX), &
                   uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_22-1)+iNX), &
                   uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_33-1)+iNX), &
                   ITERATION_Option = ITERATION(iNX,iX1,iX2,iX3), &
                   iErr_Option      = iErr     (iNX,iX1,iX2,iX3) )

          uCF    (iX1,iX2,iX3,nDOFX*(iCF_D -1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_D -1)+iNX) * Psi6
          uCF    (iX1,iX2,iX3,nDOFX*(iCF_S1-1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_S1-1)+iNX) * Psi6
          uCF    (iX1,iX2,iX3,nDOFX*(iCF_S2-1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_S2-1)+iNX) * Psi6
          uCF    (iX1,iX2,iX3,nDOFX*(iCF_S3-1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_S3-1)+iNX) * Psi6
          uCF    (iX1,iX2,iX3,nDOFX*(iCF_E -1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_E -1)+iNX) * Psi6
          uCF    (iX1,iX2,iX3,nDOFX*(iCF_Ne-1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_Ne-1)+iNX) * Psi6

          ErrorExists = ErrorExists + iErr(iNX,iX1,iX2,iX3)

          CALL ComputePressureFromPrimitive &
                 ( uPF(iPF_D), uPF(iPF_E), uPF(iPF_Ne), Pressure )

          uGS(iX1,iX2,iX3,nDOFX*(iGS_S-1)+iNX) &
            =   (  uCF(iX1,iX2,iX3,nDOFX*(iCF_S1-1)+iNX) / Psi6 * uPF(iPF_V1) &
                 + uCF(iX1,iX2,iX3,nDOFX*(iCF_S2-1)+iNX) / Psi6 * uPF(iPF_V2) &
                 + uCF(iX1,iX2,iX3,nDOFX*(iCF_S3-1)+iNX) / Psi6 * uPF(iPF_V3) &
                 + Three * Pressure ) * Psi6

          LorentzFactor &
            = One / SQRT( One                              &
                - ( uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_11-1)+iNX) &
                      * uPF(iPF_V1)**2 &
                  + uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_22-1)+iNX) &
                      * uPF(iPF_V2)**2 &
                  + uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_33-1)+iNX) &
                      * uPF(iPF_V3)**2 ) )

          DO iS   = 1, nSpecies
          DO iE   = 1, nE
          DO iN_E = 1, nDOFE

            iN_Z = ( iNX - 1 ) * nDOFE + iN_E

            iD_N  = ( iS - 1 ) * nCR * ( iE_E0 - iE_B0 + 1 ) * nDOFZ &
                      + ( iCR_N - 1 ) * ( iE_E0 - iE_B0 + 1 ) * nDOFZ &
                      + ( iE - 1 ) * nDOFZ + iN_Z
            iD_G1 = ( iS - 1 ) * nCR * ( iE_E0 - iE_B0 + 1 ) * nDOFZ &
                      + ( iCR_G1 - 1 ) * ( iE_E0 - iE_B0 + 1 ) * nDOFZ &
                      + ( iE - 1 ) * nDOFZ + iN_Z
            iD_G2 = ( iS - 1 ) * nCR * ( iE_E0 - iE_B0 + 1 ) * nDOFZ &
                      + ( iCR_G2 - 1 ) * ( iE_E0 - iE_B0 + 1 ) * nDOFZ &
                      + ( iE - 1 ) * nDOFZ + iN_Z
            iD_G3 = ( iS - 1 ) * nCR * ( iE_E0 - iE_B0 + 1 ) * nDOFZ &
                      + ( iCR_G3 - 1 ) * ( iE_E0 - iE_B0 + 1 ) * nDOFZ &
                      + ( iE - 1 ) * nDOFZ + iN_Z

            N     = uCR(iX1,iX2,iX3,iD_N )
            G_d_1 = uCR(iX1,iX2,iX3,iD_G1)
            G_d_2 = uCR(iX1,iX2,iX3,iD_G2)
            G_d_3 = uCR(iX1,iX2,iX3,iD_G3)

            V_u_1 = uPF (iPF_V1)
            V_u_2 = uPF (iPF_V2)
            V_u_3 = uPF (iPF_V3)

            vG = V_u_1 * G_d_1 + V_u_2 * G_d_2 + V_u_3 * G_d_3

            S_int = LorentzFactor * N + vG

            S &
              = S + FourPi * dZ1(iE) * WeightsE(iN_E) &
                      * uGE(iN_E,iE,iGE_Ep3) * S_int

          END DO
          END DO
          END DO

          uGS(iX1,iX2,iX3,nDOFX*(iGS_S-1)+iNX) &
            = uGS(iX1,iX2,iX3,nDOFX*(iGS_S-1)+iNX) + S

          S = Zero

        END DO
        END DO
        END DO
        END DO

        IF( ErrorExists .NE. 0 )THEN

          WRITE(*,*) &
            'ERROR'
          WRITE(*,*) &
            '-----'
          WRITE(*,*) &
            '    MODULE: Poseidon_UtilitiesModule'
          WRITE(*,*) &
            'SUBROUTINE: ComputePressureTensorTrace_XCFC_TwoMoment_MF'

          CALL CreateMesh_MF( iLevel, MeshX )

          DO iX3 = iX_B0(3), iX_E0(3)
          DO iX2 = iX_B0(2), iX_E0(2)
          DO iX1 = iX_B0(1), iX_E0(1)
          DO iNX = 1       , nDOFX

            Psi6 = uGF(iX1,iX2,iX3,nDOFX*(iGF_Psi-1)+iNX)**6

            CALL DescribeError_Euler &
                   ( iErr(iNX,iX1,iX2,iX3), &
                     Int_Option &
                       = [ ITERATION(iNX,iX1,iX2,iX3), 99999999, &
                           iX_B0(1), iX_B0(2), iX_B0(3), &
                           iX_E0(1), iX_E0(2), iX_E0(3), &
                           iNX, iX1, iX2, iX3 ], &
                     Real_Option &
                       = [ MeshX(1) % Center(iX1), &
                           MeshX(2) % Center(iX2), &
                           MeshX(3) % Center(iX3), &
                           MeshX(1) % Width (iX1), &
                           MeshX(2) % Width (iX2), &
                           MeshX(3) % Width (iX3), &
                           uCF(iX1,iX2,iX3,nDOFX*(iCF_D -1)+iNX) / Psi6, &
                           uCF(iX1,iX2,iX3,nDOFX*(iCF_S1-1)+iNX) / Psi6, &
                           uCF(iX1,iX2,iX3,nDOFX*(iCF_S2-1)+iNX) / Psi6, &
                           uCF(iX1,iX2,iX3,nDOFX*(iCF_S3-1)+iNX) / Psi6, &
                           uCF(iX1,iX2,iX3,nDOFX*(iCF_E -1)+iNX) / Psi6, &
                           uCF(iX1,iX2,iX3,nDOFX*(iCF_Ne-1)+iNX) / Psi6, &
                           uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_11-1)+iNX), &
                           uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_22-1)+iNX), &
                           uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_33-1)+iNX) ], &
                     Char_Option = [ 'NA' ], &
                     Message_Option &
                       = 'Calling from ComputePressureTensorTrace_XCFC_TwoMoment_MF' )

          END DO
          END DO
          END DO
          END DO

          CALL DestroyMesh_MF( MeshX )

        END IF

        DEALLOCATE( ITERATION )
        DEALLOCATE( iErr      )

      END DO ! WHILE( MFI % next() )

      CALL amrex_mfiter_destroy( MFI )

    END DO ! iLevel = 0, nLevels-1

    END ASSOCIATE

#endif

  END SUBROUTINE ComputePressureTensorTrace_XCFC_TwoMoment_MF




END MODULE MF_TwoMoment_XCFC_UtilitiesModule
