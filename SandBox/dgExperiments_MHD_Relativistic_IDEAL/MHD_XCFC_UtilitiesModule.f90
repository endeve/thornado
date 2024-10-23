MODULE MHD_XCFC_UtilitiesModule

  USE KindModule, ONLY: &
    DP, &
    Zero, &
    One, &
    Two, &
    Three, &
    FourPi
  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    nDimsX
  USE ReferenceElementModuleX, ONLY: &
    WeightsX_q, &
    NodesX1
  USE MeshModule, ONLY: &
    MeshX
  USE GeometryFieldsModule, ONLY: &
    iGF_Phi_N, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_SqrtGm, &
    iGF_Alpha, &
    iGF_Beta_1, &
    iGF_Beta_2, &
    iGF_Beta_3, &
    iGF_Psi, &
    nGF
  USE MagnetofluidFieldsModule, ONLY: &
    iCM_D, &
    iCM_S1, &
    iCM_S2, &
    iCM_S3, &
    iCM_E, &
    iCM_Ne, &
    iCM_B1, &
    iCM_B2, &
    iCM_B3, &
    iCM_Chi, &
    nCM, &
    iPM_D, &
    iPM_V1, &
    iPM_V2, &
    iPM_V3, &
    iPM_E, &
    iPM_Ne, &
    iPM_B1, &
    iPM_B2, &
    iPM_B3, &
    iPM_Chi, &
    nPM, &
    iAM_P
  USE MHD_UtilitiesModule_Relativistic, ONLY: &
    ComputeConserved_MHD_Relativistic, &
    ComputePrimitive_MHD_Relativistic
  USE EquationOfStateModule, ONLY: &
    ComputePressureFromPrimitive
  USE MHD_ErrorModule, ONLY: &
    DescribeError_MHD
  USE XCFC_UtilitiesModule, ONLY: &
    iGS_E, &
    iGS_S1, &
    iGS_S2, &
    iGS_S3, &
    iGS_S, &
    iGS_Mg, &
    nGS, &
    nMF, &
    swX_GS, &
    MultiplyWithPsi6, &
    UpdateConformalFactorAndMetric_XCFC, &
    UpdateLapseShiftCurvature_XCFC, &
    ApplyBoundaryConditions_Geometry_XCFC
  USE GravitySolutionModule_XCFC, ONLY: &
    ComputeConformalFactor_XCFC, &
    ComputeLapseShiftCurvature_XCFC
  USE TimersModule_MHD, ONLY: &
    TimersStart_MHD, &
    TimersStop_MHD, &
    Timer_GS_ComputeSourceTerms

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeConformalFactorSourcesAndMg_XCFC_MHD
  PUBLIC :: ComputePressureTensorTrace_XCFC_MHD
  PUBLIC :: InitializeMetric_MHD
  PUBLIC :: ComputeNewtonianPotential_SphericalSymmetry

  ! https://amrex-codes.github.io/amrex/docs_html/Basics.html#fine-mask
  INTEGER, PARAMETER :: iLeaf    = 0
  INTEGER, PARAMETER :: iNotLeaf = 1

CONTAINS


  SUBROUTINE ComputeConformalFactorSourcesAndMg_XCFC_MHD &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, GS, EvolveOnlyMagnetic, Mask_Option )

    INTEGER,  INTENT(in)    :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: G (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: U (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:) ! psi^6*U
    REAL(DP), INTENT(inout) :: GS(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)
    LOGICAL,  INTENT(in)    :: EvolveOnlyMagnetic
    INTEGER , INTENT(in), OPTIONAL :: &
      Mask_Option(iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    REAL(DP) :: uGF(nGF), uCM(nCM), uPM(nPM), Psi6, Pressure, &
                LorentzFactor, EnthalpyDensity, BetaDotV
    INTEGER :: iNX, iX1, iX2, iX3, iGF, iCM

    INTEGER :: ITERATION(1:nDOFX,iX_B0(1):iX_E0(1), &
                                 iX_B0(2):iX_E0(2), &
                                 iX_B0(3):iX_E0(3))
    INTEGER :: iErr     (1:nDOFX,iX_B0(1):iX_E0(1), &
                                 iX_B0(2):iX_E0(2), &
                                 iX_B0(3):iX_E0(3))
    INTEGER :: Mask             (iX_B1(1):iX_E1(1), &
                                 iX_B1(2):iX_E1(2), &
                                 iX_B1(3):iX_E1(3),1)

    CALL TimersStart_MHD( Timer_GS_ComputeSourceTerms )

    IF( PRESENT( Mask_Option ) )THEN

      Mask = Mask_Option

    ELSE

      Mask = iLeaf

    END IF

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      ITERATION(:,iX1,iX2,iX3) = 0
      iErr     (:,iX1,iX2,iX3) = 0

      IF( IsNotLeafElement( Mask(iX1,iX2,iX3,1) ) ) CYCLE

      DO iNX = 1, nDOFX

        GS(iNX,iX1,iX2,iX3,iGS_E) &
          = U(iNX,iX1,iX2,iX3,iCM_E) + U(iNX,iX1,iX2,iX3,iCM_D)

        GS(iNX,iX1,iX2,iX3,iGS_S1) = U(iNX,iX1,iX2,iX3,iCM_S1)
        GS(iNX,iX1,iX2,iX3,iGS_S2) = U(iNX,iX1,iX2,iX3,iCM_S2)
        GS(iNX,iX1,iX2,iX3,iGS_S3) = U(iNX,iX1,iX2,iX3,iCM_S3)

        ! --- Compute gravitational mass integrand
        !     (Eq. (12.45) from Rezzolla & Zanotti, Relativistic Hyrodynamics.
        !     valid only for stationary, axisymmetric spacetimes) ---

        DO iGF = 1, nGF

          uGF(iGF) = G(iNX,iX1,iX2,iX3,iGF)

        END DO

        Psi6 = uGF(iGF_Psi)**6

        DO iCM = 1, nCM

          uCM(iCM) = U(iNX,iX1,iX2,iX3,iCM) / Psi6

        END DO

        CALL ComputePrimitive_MHD_Relativistic &
               ( uCM(iCM_D  ), &
                 uCM(iCM_S1 ), &
                 uCM(iCM_S2 ), &
                 uCM(iCM_S3 ), &
                 uCM(iCM_E  ), &
                 uCM(iCM_Ne ), &
                 uCM(iCM_B1 ), &
                 uCM(iCM_B2 ), &
                 uCM(iCM_B3 ), &
                 uCM(iCM_Chi), &
                 uPM(iPM_D  ), &
                 uPM(iPM_V1 ), &
                 uPM(iPM_V2 ), &
                 uPM(iPM_V3 ), &
                 uPM(iPM_E  ), &
                 uPM(iPM_Ne ), &
                 uPM(iPM_B1 ), &
                 uPM(iPM_B2 ), &
                 uPM(iPM_B3 ), &
                 uPM(iPM_Chi), &
                 uGF(iGF_Gm_dd_11), &
                 uGF(iGF_Gm_dd_22), &
                 uGF(iGF_Gm_dd_33), &
                 uGF(iGF_Alpha   ), &
                 uGF(iGF_Beta_1  ), &
                 uGF(iGF_Beta_2  ), &
                 uGF(iGF_Beta_3  ), &
                 EvolveOnlyMagnetic )

        DO iCM = 1, nCM

          U(iNX,iX1,iX2,iX3,iCM) = uCM(iCM) * Psi6

        END DO

         CALL ComputePressureFromPrimitive &
                ( uPM(iPM_D), uPM(iPM_E), uPM(iPM_Ne), Pressure )

         LorentzFactor &
           = One / SQRT( One                              &
               - ( uGF(iGF_Gm_dd_11) * uPM(iPM_V1)**2 &
                 + uGF(iGF_Gm_dd_22) * uPM(iPM_V2)**2 &
                 + uGF(iGF_Gm_dd_33) * uPM(iPM_V3)**2 ) )

         BetaDotV =   uGF(iGF_Gm_dd_11) * uGF(iGF_Beta_1) * uPM(iPM_V1) &
                    + uGF(iGF_Gm_dd_22) * uGF(iGF_Beta_2) * uPM(iPM_V2) &
                    + uGF(iGF_Gm_dd_33) * uGF(iGF_Beta_3) * uPM(iPM_V3)

         EnthalpyDensity = uPM(iPM_D) + uPM(iPM_E) + Pressure

         GS(iNX,iX1,iX2,iX3,iGS_Mg) &
           = ( Two * EnthalpyDensity * LorentzFactor**2 &
                 * ( One - BetaDotV / uGF(iGF_Alpha) ) &
                     - EnthalpyDensity + Two * Pressure ) &
               * uGF(iGF_Alpha) * uGF(iGF_SqrtGm)

      END DO

    END DO
    END DO
    END DO

    IF( ANY( iErr .NE. 0 ) )THEN

      WRITE(*,*) 'ERROR'
      WRITE(*,*) '-----'
      WRITE(*,*) '    MODULE: MHD_XCFC_UtilitiesModule'
      WRITE(*,*) 'SUBROUTINE: ComputeConformalFactorSourcesAndMg_XCFC_MHD'

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)

        IF( IsNotLeafElement( Mask(iX1,iX2,iX3,1) ) ) CYCLE

        DO iNX = 1, nDOFX

          Psi6 = G(iNX,iX1,iX2,iX3,iGF_Psi)**6

          CALL DescribeError_MHD &
            ( iErr(iNX,iX1,iX2,iX3), &
              Int_Option = [ ITERATION(iNX,iX1,iX2,iX3), 99999999, &
                             iX_B0(1), iX_B0(2), iX_B0(3), &
                             iX_E0(1), iX_E0(2), iX_E0(3), &
                             iNX, iX1, iX2, iX3 ], &
              Real_Option = [ MeshX(1) % Center(iX1), &
                              MeshX(2) % Center(iX2), &
                              MeshX(3) % Center(iX3), &
                              MeshX(1) % Width (iX1), &
                              MeshX(2) % Width (iX2), &
                              MeshX(3) % Width (iX3), &
                              U(iNX,iX1,iX2,iX3,iCM_D ) / Psi6, &
                              U(iNX,iX1,iX2,iX3,iCM_S1) / Psi6, &
                              U(iNX,iX1,iX2,iX3,iCM_S2) / Psi6, &
                              U(iNX,iX1,iX2,iX3,iCM_S3) / Psi6, &
                              U(iNX,iX1,iX2,iX3,iCM_E ) / Psi6, &
                              U(iNX,iX1,iX2,iX3,iCM_Ne) / Psi6, &
                              G(iNX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                              G(iNX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                              G(iNX,iX1,iX2,iX3,iGF_Gm_dd_33) ], &
              Char_Option = [ 'NA' ], &
              Message_Option &
                = 'Calling from ComputeConformalFactorSourcesAndMg_XCFC_MHD' )

        END DO

      END DO
      END DO
      END DO

    END IF

    CALL TimersStop_MHD( Timer_GS_ComputeSourceTerms )

  END SUBROUTINE ComputeConformalFactorSourcesAndMg_XCFC_MHD


  SUBROUTINE ComputePressureTensorTrace_XCFC_MHD &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, GS, EvolveOnlyMagnetic, Mask_Option )

    INTEGER,  INTENT(in)    :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: G (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: U (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:) ! psi^6*U
    REAL(DP), INTENT(inout) :: GS(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)
    LOGICAL,  INTENT(in)    :: EvolveOnlyMagnetic
    INTEGER , INTENT(in), OPTIONAL :: &
      Mask_Option(iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    REAL(DP) :: uGF(nGF), uCM(nCM), uPM(nPM), Psi6, Pressure
    INTEGER  :: iNX, iX1, iX2, iX3, iGF, iCM

    INTEGER :: ITERATION(1:nDOFX,iX_B0(1):iX_E0(1), &
                                 iX_B0(2):iX_E0(2), &
                                 iX_B0(3):iX_E0(3))
    INTEGER :: iErr     (1:nDOFX,iX_B0(1):iX_E0(1), &
                                 iX_B0(2):iX_E0(2), &
                                 iX_B0(3):iX_E0(3))
    INTEGER :: Mask             (iX_B1(1):iX_E1(1), &
                                 iX_B1(2):iX_E1(2), &
                                 iX_B1(3):iX_E1(3),1)

    CALL TimersStart_MHD( Timer_GS_ComputeSourceTerms )

    IF( PRESENT( Mask_Option ) )THEN

      Mask = Mask_Option

    ELSE

      Mask = iLeaf

    END IF

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      ITERATION(:,iX1,iX2,iX3) = 0
      iErr     (:,iX1,iX2,iX3) = 0

      IF( IsNotLeafElement( Mask(iX1,iX2,iX3,1) ) ) CYCLE

      DO iNX = 1, nDOFX

        DO iGF = 1, nGF

          uGF(iGF) = G(iNX,iX1,iX2,iX3,iGF)

        END DO

        Psi6 = uGF(iGF_Psi)**6

        DO iCM = 1, nCM

          uCM(iCM) = U(iNX,iX1,iX2,iX3,iCM) / Psi6

        END DO

        ! --- Compute trace of stress tensor ---

        CALL ComputePrimitive_MHD_Relativistic &
               ( uCM(iCM_D  ), &
                 uCM(iCM_S1 ), &
                 uCM(iCM_S2 ), &
                 uCM(iCM_S3 ), &
                 uCM(iCM_E  ), &
                 uCM(iCM_Ne ), &
                 uCM(iCM_B1 ), &
                 uCM(iCM_B2 ), &
                 uCM(iCM_B3 ), &
                 uCM(iCM_Chi), &
                 uPM(iPM_D  ), &
                 uPM(iPM_V1 ), &
                 uPM(iPM_V2 ), &
                 uPM(iPM_V3 ), &
                 uPM(iPM_E  ), &
                 uPM(iPM_Ne ), &
                 uPM(iPM_B1 ), &
                 uPM(iPM_B2 ), &
                 uPM(iPM_B3 ), &
                 uPM(iPM_Chi), &
                 uGF(iGF_Gm_dd_11), &
                 uGF(iGF_Gm_dd_22), &
                 uGF(iGF_Gm_dd_33), &
                 uGF(iGF_Alpha   ), &
                 uGF(iGF_Beta_1  ), &
                 uGF(iGF_Beta_2  ), &
                 uGF(iGF_Beta_3  ), &
                 EvolveOnlyMagnetic )

         DO iCM = 1, nCM

          U(iNX,iX1,iX2,iX3,iCM) = uCM(iCM) * Psi6

        END DO

        CALL ComputePressureFromPrimitive &
               ( uPM(iPM_D), uPM(iPM_E), uPM(iPM_Ne), Pressure )

        GS(iNX,iX1,iX2,iX3,iGS_S) &
          = Psi6 * ( uCM(iCM_S1) * uPM(iPM_V1) &
                   + uCM(iCM_S2) * uPM(iPM_V2) &
                   + uCM(iCM_S3) * uPM(iPM_V3) &
                   + Three * Pressure )

      END DO

    END DO
    END DO
    END DO

    IF( ANY( iErr .NE. 0 ) )THEN

      WRITE(*,*) 'ERROR'
      WRITE(*,*) '-----'
      WRITE(*,*) '    MODULE: MHD_XCFC_UtilitiesModule'
      WRITE(*,*) 'SUBROUTINE: ComputePressureTensorTrace_XCFC_MHD'

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)

        IF( IsNotLeafElement( Mask(iX1,iX2,iX3,1) ) ) CYCLE

        DO iNX = 1, nDOFX

          IF( iErr(iNX,iX1,iX2,iX3) .NE. 0 )THEN

            WRITE(*,'(2x,A,4I5.4)') 'iNX, iX1, iX2, iX3 = ', iNX, iX1, iX2, iX3
            CALL DescribeError_MHD( iErr(iNX,iX1,iX2,iX3) )

          END IF

        END DO

      END DO
      END DO
      END DO

    END IF

    CALL TimersStop_MHD( Timer_GS_ComputeSourceTerms )

  END SUBROUTINE ComputePressureTensorTrace_XCFC_MHD


  SUBROUTINE InitializeMetric_MHD &
    ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCM, uPM, uAM, EvolveOnlyMagnetic )

    INTEGER , INTENT(in)    :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout) :: uGF(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
                               uCM(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
                               uPM(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
                               uAM(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    LOGICAL,  INTENT(in)    :: EvolveOnlyMagnetic


    INTEGER             :: iNX, iX1, iX2, iX3
    INTEGER             :: ITER
    INTEGER , PARAMETER :: MAX_ITER  = 10
    REAL(DP), PARAMETER :: TOLERANCE = 1.0e-13_DP
    REAL(DP)            :: dAlpha, dPsi
    LOGICAL             :: CONVERGED

    REAL(DP) :: uGS(nDOFX,iX_B0(1):iX_E0(1), &
                          iX_B0(2):iX_E0(2), &
                          iX_B0(3):iX_E0(3),nGS)
    REAL(DP) :: uMF(nDOFX,iX_B0(1)-swX_GS(1):iX_E0(1)+swX_GS(1), &
                          iX_B0(2)-swX_GS(2):iX_E0(2)+swX_GS(2), &
                          iX_B0(3)-swX_GS(3):iX_E0(3)+swX_GS(3),nMF)

    REAL(DP) :: dAl1(nDOFX,iX_B0(1):iX_E0(1), &
                           iX_B0(2):iX_E0(2), &
                           iX_B0(3):iX_E0(3))
    REAL(DP) :: dCF1(nDOFX,iX_B0(1):iX_E0(1), &
                           iX_B0(2):iX_E0(2), &
                           iX_B0(3):iX_E0(3))
    REAL(DP) :: dAl2(nDOFX,iX_B0(1):iX_E0(1), &
                           iX_B0(2):iX_E0(2), &
                           iX_B0(3):iX_E0(3))
    REAL(DP) :: dCF2(nDOFX,iX_B0(1):iX_E0(1), &
                           iX_B0(2):iX_E0(2), &
                           iX_B0(3):iX_E0(3))

    CONVERGED = .FALSE.
    ITER      = 0

    DO WHILE( .NOT. CONVERGED )

      ITER = ITER + 1

      dAl1 = uGF(:,iX_B0(1):iX_E0(1), &
                   iX_B0(2):iX_E0(2), &
                   iX_B0(3):iX_E0(3),iGF_Alpha)
      dCF1 = uGF(:,iX_B0(1):iX_E0(1), &
                   iX_B0(2):iX_E0(2), &
                   iX_B0(3):iX_E0(3),iGF_Psi  )

      CALL MultiplyWithPsi6( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCM, +1 )

      CALL ComputeConformalFactor &
             ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCM, uMF, uGS, EvolveOnlyMagnetic )

      CALL ComputeLapseShiftCurvature &
             ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCM, uMF, uGS, EvolveOnlyMagnetic )

      CALL MultiplyWithPsi6( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCM, -1 )

      dAl2 = uGF(:,iX_B0(1):iX_E0(1), &
                   iX_B0(2):iX_E0(2), &
                   iX_B0(3):iX_E0(3),iGF_Alpha)
      dCF2 = uGF(:,iX_B0(1):iX_E0(1), &
                   iX_B0(2):iX_E0(2), &
                   iX_B0(3):iX_E0(3),iGF_Psi  )

      dAlpha = MAXVAL( ABS( dAl2 - dAl1 ) )
      dPsi   = MAXVAL( ABS( dCF2 - dCF1 ) )

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iNX = 1       , nDOFX

        CALL ComputeConserved_MHD_Relativistic &
               ( uPM(iNX,iX1,iX2,iX3,iPM_D       ), &
                 uPM(iNX,iX1,iX2,iX3,iPM_V1      ), &
                 uPM(iNX,iX1,iX2,iX3,iPM_V2      ), &
                 uPM(iNX,iX1,iX2,iX3,iPM_V3      ), &
                 uPM(iNX,iX1,iX2,iX3,iPM_E       ), &
                 uPM(iNX,iX1,iX2,iX3,iPM_Ne      ), &
                 uPM(iNX,iX1,iX2,iX3,iPM_B1      ), &
                 uPM(iNX,iX1,iX2,iX3,iPM_B2      ), &
                 uPM(iNX,iX1,iX2,iX3,iPM_B3      ), &
                 uPM(iNX,iX1,iX2,iX3,iPM_Chi     ), &
                 uCM(iNX,iX1,iX2,iX3,iCM_D       ), &
                 uCM(iNX,iX1,iX2,iX3,iCM_S1      ), &
                 uCM(iNX,iX1,iX2,iX3,iCM_S2      ), &
                 uCM(iNX,iX1,iX2,iX3,iCM_S3      ), &
                 uCM(iNX,iX1,iX2,iX3,iCM_E       ), &
                 uCM(iNX,iX1,iX2,iX3,iCM_Ne      ), &
                 uCM(iNX,iX1,iX2,iX3,iCM_B1      ), &
                 uCM(iNX,iX1,iX2,iX3,iCM_B2      ), &
                 uCM(iNX,iX1,iX2,iX3,iCM_B3      ), &
                 uCM(iNX,iX1,iX2,iX3,iCM_Chi     ), &
                 uGF(iNX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 uGF(iNX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 uGF(iNX,iX1,iX2,iX3,iGF_Gm_dd_33), &
                 uGF(iNX,iX1,iX2,iX3,iGF_Alpha   ), &
                 uGF(iNX,iX1,iX2,iX3,iGF_Beta_1  ), &
                 uGF(iNX,iX1,iX2,iX3,iGF_Beta_2  ), &
                 uGF(iNX,iX1,iX2,iX3,iGF_Beta_3  ), &
                 uAM(iNX,iX1,iX2,iX3,iAM_P       ), &
                 EvolveOnlyMagnetic )

      END DO
      END DO
      END DO
      END DO

      IF( MAX( dAlpha, dPsi ) .LT. TOLERANCE ) CONVERGED = .TRUE.

      IF( ITER .EQ. MAX_ITER )THEN

        WRITE(*,*) 'Could not initialize fields. Exiting...'
        STOP

      END IF

    END DO

  END SUBROUTINE InitializeMetric_MHD


  SUBROUTINE ComputeNewtonianPotential_SphericalSymmetry &
    ( iX_B0, iX_E0, iX_B1, iX_E1, P, G )

    INTEGER,  INTENT(in)    :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: P(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER  :: iX1
    REAL(DP) :: X1C, dX, X1q(nDOFX), dM, dPhi
    REAL(DP) :: EnclosedMass(iX_B0(1):iX_E0(1))

    IF( .NOT. nDimsX .EQ. 1 ) RETURN

    ! --- Compute enclosed mass ---

    dM = Zero

    DO iX1 = iX_B0(1), iX_E0(1)

      X1C = MeshX(1) % Center(iX1)
      dX  = MeshX(1) % Width (iX1)

      X1q = X1C + NodesX1 * dX

      dM &
        = dM &
            + FourPi * dX &
                 * SUM( WeightsX_q * X1q**2 * P(:,iX1,1,1,iPM_D) )

      EnclosedMass(iX1) = dM

    END DO

    ! --- Compute Newtonian gravitational potential ---

    G(:,iX_E0(1),1,1,iGF_Phi_N) &
      = -EnclosedMass(iX_E0(1)) / MeshX(1) % Center(iX_E0(1))

    dPhi = Zero

    DO iX1 = iX_E0(1)-1, iX_B0(1), -1

      X1C = MeshX(1) % Center(iX1)
      dX  = MeshX(1) % Width (iX1)

      dPhi = dPhi - EnclosedMass(iX1) / X1C**2 * dX

      G(:,iX1,1,1,iGF_Phi_N) = G(:,iX_E0(1),1,1,iGF_Phi_N) + dPhi

    END DO

  END SUBROUTINE ComputeNewtonianPotential_SphericalSymmetry


  ! --- PRIVATE SUBROUTINES ---


  LOGICAL FUNCTION IsNotLeafElement( Element )

    INTEGER, INTENT(in) :: Element

    IF( Element .EQ. iNotLeaf )THEN
      IsNotLeafElement = .TRUE.
    ELSE
      IsNotLeafElement = .FALSE.
    END IF

    RETURN
  END FUNCTION IsNotLeafElement


  SUBROUTINE ComputeConformalFactor &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, Ustar, M, GS, EvolveOnlyMagnetic )

    INTEGER , INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout) :: &
      G    (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
      Ustar(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
      M    (1:,iX_B0(1)-swX_GS(1):, &
               iX_B0(2)-swX_GS(2):, &
               iX_B0(3)-swX_GS(3):,1:), &
      GS   (1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)
    LOGICAL,  INTENT(in)    :: EvolveOnlyMagnetic

    CALL ComputeConformalFactorSourcesAndMg_XCFC_MHD &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, Ustar, GS, EvolveOnlyMagnetic )

    CALL ComputeConformalFactor_XCFC &
           ( iX_B0, iX_E0, iX_B1, iX_E1, GS, M )

    CALL UpdateConformalFactorAndMetric_XCFC &
           ( iX_B0, iX_E0, iX_B1, iX_E1, M, G )

  END SUBROUTINE ComputeConformalFactor


  SUBROUTINE ComputeLapseShiftCurvature &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, Ustar, M, GS, EvolveOnlyMagnetic )

    INTEGER , INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout) :: &
      G    (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
      Ustar(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
      M    (1:,iX_B0(1)-swX_GS(1):, &
               iX_B0(2)-swX_GS(2):, &
               iX_B0(3)-swX_GS(3):,1:), &
      GS   (1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)
    LOGICAL, INTENT(in)     :: EvolveOnlyMagnetic

    CALL ComputePressureTensorTrace_XCFC_MHD &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, Ustar, GS, EvolveOnlyMagnetic )

    CALL ComputeLapseShiftCurvature_XCFC &
           ( iX_B0, iX_E0, iX_B1, iX_E1, GS, M )

    CALL UpdateLapseShiftCurvature_XCFC &
           ( iX_B0, iX_E0, iX_B1, iX_E1, M, G )

    CALL ApplyBoundaryConditions_Geometry_XCFC &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G )

  END SUBROUTINE ComputeLapseShiftCurvature


END MODULE MHD_XCFC_UtilitiesModule
