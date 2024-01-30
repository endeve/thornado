MODULE Euler_XCFC_UtilitiesModule

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
  USE FluidFieldsModule, ONLY: &
    iCF_D, &
    iCF_S1, &
    iCF_S2, &
    iCF_S3, &
    iCF_E, &
    iCF_Ne, &
    nCF, &
    iPF_D, &
    iPF_V1, &
    iPF_V2, &
    iPF_V3, &
    iPF_E, &
    iPF_Ne, &
    nPF, &
    iAF_P
  USE Euler_UtilitiesModule_Relativistic, ONLY: &
    ComputeConserved_Euler_Relativistic, &
    ComputePrimitive_Euler_Relativistic
  USE EquationOfStateModule, ONLY: &
    ComputePressureFromPrimitive
  USE Euler_ErrorModule, ONLY: &
    DescribeError_Euler
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
  USE TimersModule_Euler, ONLY: &
    TimersStart_Euler, &
    TimersStop_Euler, &
    Timer_GS_ComputeSourceTerms

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeConformalFactorSourcesAndMg_XCFC_Euler
  PUBLIC :: ComputePressureTensorTrace_XCFC_Euler
  PUBLIC :: InitializeMetric_Euler
  PUBLIC :: ComputeNewtonianPotential_SphericalSymmetry

  ! https://amrex-codes.github.io/amrex/docs_html/Basics.html#fine-mask
  INTEGER, PARAMETER :: iLeaf    = 0
  INTEGER, PARAMETER :: iNotLeaf = 1

CONTAINS


  SUBROUTINE ComputeConformalFactorSourcesAndMg_XCFC_Euler &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, GS, Mask_Option )

    INTEGER,  INTENT(in)    :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: G (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: U (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:) ! psi^6*U
    REAL(DP), INTENT(inout) :: GS(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)
    INTEGER , INTENT(in), OPTIONAL :: &
      Mask_Option(iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    REAL(DP) :: uGF(nGF), uCF(nCF), uPF(nPF), Psi6, Pressure, &
                LorentzFactor, Enthalpy, BetaDotV
    INTEGER :: iNX, iX1, iX2, iX3, iGF, iCF

    INTEGER :: ITERATION(1:nDOFX,iX_B0(1):iX_E0(1), &
                                 iX_B0(2):iX_E0(2), &
                                 iX_B0(3):iX_E0(3))
    INTEGER :: iErr     (1:nDOFX,iX_B0(1):iX_E0(1), &
                                 iX_B0(2):iX_E0(2), &
                                 iX_B0(3):iX_E0(3))
    INTEGER :: Mask             (iX_B1(1):iX_E1(1), &
                                 iX_B1(2):iX_E1(2), &
                                 iX_B1(3):iX_E1(3),1)

    CALL TimersStart_Euler( Timer_GS_ComputeSourceTerms )

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
          = U(iNX,iX1,iX2,iX3,iCF_E) + U(iNX,iX1,iX2,iX3,iCF_D)

        GS(iNX,iX1,iX2,iX3,iGS_S1) = U(iNX,iX1,iX2,iX3,iCF_S1)
        GS(iNX,iX1,iX2,iX3,iGS_S2) = U(iNX,iX1,iX2,iX3,iCF_S2)
        GS(iNX,iX1,iX2,iX3,iGS_S3) = U(iNX,iX1,iX2,iX3,iCF_S3)

        ! --- Compute gravitational mass ---

        DO iGF = 1, nGF

          uGF(iGF) = G(iNX,iX1,iX2,iX3,iGF)

        END DO

        Psi6 = uGF(iGF_Psi)**6

        DO iCF = 1, nCF

          uCF(iCF) = U(iNX,iX1,iX2,iX3,iCF) / Psi6

        END DO

        CALL ComputePrimitive_Euler_Relativistic &
               ( uCF(iCF_D ), &
                 uCF(iCF_S1), &
                 uCF(iCF_S2), &
                 uCF(iCF_S3), &
                 uCF(iCF_E ), &
                 uCF(iCF_Ne), &
                 uPF(iPF_D ), &
                 uPF(iPF_V1), &
                 uPF(iPF_V2), &
                 uPF(iPF_V3), &
                 uPF(iPF_E ), &
                 uPF(iPF_Ne), &
                 uGF(iGF_Gm_dd_11), &
                 uGF(iGF_Gm_dd_22), &
                 uGF(iGF_Gm_dd_33), &
                 ITERATION_Option = ITERATION(iNX,iX1,iX2,iX3), &
                 iErr_Option      = iErr     (iNX,iX1,iX2,iX3) )

        DO iCF = 1, nCF

          U(iNX,iX1,iX2,iX3,iCF) = uCF(iCF) * Psi6

        END DO

         CALL ComputePressureFromPrimitive &
                ( uPF(iPF_D), uPF(iPF_E), uPF(iPF_Ne), Pressure )

         LorentzFactor &
           = One / SQRT( One                              &
               - ( uGF(iGF_Gm_dd_11) * uPF(iPF_V1)**2 &
                 + uGF(iGF_Gm_dd_22) * uPF(iPF_V2)**2 &
                 + uGF(iGF_Gm_dd_33) * uPF(iPF_V3)**2 ) )

         BetaDotV =   uGF(iGF_Gm_dd_11) * uGF(iGF_Beta_1) * uPF(iPF_V1) &
                    + uGF(iGF_Gm_dd_22) * uGF(iGF_Beta_2) * uPF(iPF_V2) &
                    + uGF(iGF_Gm_dd_33) * uGF(iGF_Beta_3) * uPF(iPF_V3)

         Enthalpy = uPF(iPF_D) + uPF(iPF_E) + Pressure

         GS(iNX,iX1,iX2,iX3,iGS_Mg) &
           = ( Enthalpy * ( Two * LorentzFactor**2             &
                 * ( One - BetaDotV / uGF(iGF_Alpha) ) - One ) &
                 + Two * Pressure ) * uGF(iGF_Alpha) * uGF(iGF_SqrtGm)

      END DO

    END DO
    END DO
    END DO

    IF( ANY( iErr .NE. 0 ) )THEN

      WRITE(*,*) 'ERROR'
      WRITE(*,*) '-----'
      WRITE(*,*) '    MODULE: Euler_XCFC_UtilitiesModule'
      WRITE(*,*) 'SUBROUTINE: ComputeConformalFactorSourcesAndMg_XCFC_Euler'

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)

        IF( IsNotLeafElement( Mask(iX1,iX2,iX3,1) ) ) CYCLE

        DO iNX = 1, nDOFX

          Psi6 = G(iNX,iX1,iX2,iX3,iGF_Psi)**6

          CALL DescribeError_Euler &
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
                              U(iNX,iX1,iX2,iX3,iCF_D ) / Psi6, &
                              U(iNX,iX1,iX2,iX3,iCF_S1) / Psi6, &
                              U(iNX,iX1,iX2,iX3,iCF_S2) / Psi6, &
                              U(iNX,iX1,iX2,iX3,iCF_S3) / Psi6, &
                              U(iNX,iX1,iX2,iX3,iCF_E ) / Psi6, &
                              U(iNX,iX1,iX2,iX3,iCF_Ne) / Psi6, &
                              G(iNX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                              G(iNX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                              G(iNX,iX1,iX2,iX3,iGF_Gm_dd_33) ], &
              Char_Option = [ 'NA' ], &
              Message_Option &
                = 'Calling from ComputeConformalFactorSourcesAndMg_XCFC_Euler' )

        END DO

      END DO
      END DO
      END DO

    END IF

    CALL TimersStop_Euler( Timer_GS_ComputeSourceTerms )

  END SUBROUTINE ComputeConformalFactorSourcesAndMg_XCFC_Euler


  SUBROUTINE ComputePressureTensorTrace_XCFC_Euler &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, GS, Mask_Option )

    INTEGER,  INTENT(in)    :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: G (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: U (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:) ! psi^6*U
    REAL(DP), INTENT(inout) :: GS(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)
    INTEGER , INTENT(in), OPTIONAL :: &
      Mask_Option(iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    REAL(DP) :: uGF(nGF), uCF(nCF), uPF(nPF), Psi6, Pressure
    INTEGER  :: iNX, iX1, iX2, iX3, iGF, iCF

    INTEGER :: ITERATION(1:nDOFX,iX_B0(1):iX_E0(1), &
                                 iX_B0(2):iX_E0(2), &
                                 iX_B0(3):iX_E0(3))
    INTEGER :: iErr     (1:nDOFX,iX_B0(1):iX_E0(1), &
                                 iX_B0(2):iX_E0(2), &
                                 iX_B0(3):iX_E0(3))
    INTEGER :: Mask             (iX_B1(1):iX_E1(1), &
                                 iX_B1(2):iX_E1(2), &
                                 iX_B1(3):iX_E1(3),1)

    CALL TimersStart_Euler( Timer_GS_ComputeSourceTerms )

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

        DO iCF = 1, nCF

          uCF(iCF) = U(iNX,iX1,iX2,iX3,iCF) / Psi6

        END DO

        ! --- Compute trace of stress tensor ---

        CALL ComputePrimitive_Euler_Relativistic &
               ( uCF(iCF_D ), &
                 uCF(iCF_S1), &
                 uCF(iCF_S2), &
                 uCF(iCF_S3), &
                 uCF(iCF_E ), &
                 uCF(iCF_Ne), &
                 uPF(iPF_D ), &
                 uPF(iPF_V1), &
                 uPF(iPF_V2), &
                 uPF(iPF_V3), &
                 uPF(iPF_E ), &
                 uPF(iPF_Ne), &
                 uGF(iGF_Gm_dd_11), &
                 uGF(iGF_Gm_dd_22), &
                 uGF(iGF_Gm_dd_33), &
                 ITERATION_Option = ITERATION(iNX,iX1,iX2,iX3), &
                 iErr_Option      = iErr     (iNX,iX1,iX2,iX3) )

        DO iCF = 1, nCF

          U(iNX,iX1,iX2,iX3,iCF) = uCF(iCF) * Psi6

        END DO

        CALL ComputePressureFromPrimitive &
               ( uPF(iPF_D), uPF(iPF_E), uPF(iPF_Ne), Pressure )

        GS(iNX,iX1,iX2,iX3,iGS_S) &
          = Psi6 * ( uCF(iCF_S1) * uPF(iPF_V1) &
                   + uCF(iCF_S2) * uPF(iPF_V2) &
                   + uCF(iCF_S3) * uPF(iPF_V3) &
                   + Three * Pressure )

      END DO

    END DO
    END DO
    END DO

    IF( ANY( iErr .NE. 0 ) )THEN

      WRITE(*,*) 'ERROR'
      WRITE(*,*) '-----'
      WRITE(*,*) '    MODULE: Euler_XCFC_UtilitiesModule'
      WRITE(*,*) 'SUBROUTINE: ComputePressureTensorTrace_XCFC_Euler'

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)

        IF( IsNotLeafElement( Mask(iX1,iX2,iX3,1) ) ) CYCLE

        DO iNX = 1, nDOFX

          IF( iErr(iNX,iX1,iX2,iX3) .NE. 0 )THEN

            WRITE(*,'(2x,A,4I5.4)') 'iNX, iX1, iX2, iX3 = ', iNX, iX1, iX2, iX3
            CALL DescribeError_Euler( iErr(iNX,iX1,iX2,iX3) )

          END IF

        END DO

      END DO
      END DO
      END DO

    END IF

    CALL TimersStop_Euler( Timer_GS_ComputeSourceTerms )

  END SUBROUTINE ComputePressureTensorTrace_XCFC_Euler


  SUBROUTINE InitializeMetric_Euler &
    ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uPF, uAF )

    INTEGER , INTENT(in)    :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout) :: uGF(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
                               uCF(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
                               uPF(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
                               uAF(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)


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

      CALL MultiplyWithPsi6( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, +1 )

      CALL ComputeConformalFactor &
             ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uMF, uGS )

      CALL ComputeLapseShiftCurvature &
             ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uMF, uGS )

      CALL MultiplyWithPsi6( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, -1 )

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

        CALL ComputeConserved_Euler_Relativistic &
               ( uPF(iNX,iX1,iX2,iX3,iPF_D       ), &
                 uPF(iNX,iX1,iX2,iX3,iPF_V1      ), &
                 uPF(iNX,iX1,iX2,iX3,iPF_V2      ), &
                 uPF(iNX,iX1,iX2,iX3,iPF_V3      ), &
                 uPF(iNX,iX1,iX2,iX3,iPF_E       ), &
                 uPF(iNX,iX1,iX2,iX3,iPF_Ne      ), &
                 uCF(iNX,iX1,iX2,iX3,iCF_D       ), &
                 uCF(iNX,iX1,iX2,iX3,iCF_S1      ), &
                 uCF(iNX,iX1,iX2,iX3,iCF_S2      ), &
                 uCF(iNX,iX1,iX2,iX3,iCF_S3      ), &
                 uCF(iNX,iX1,iX2,iX3,iCF_E       ), &
                 uCF(iNX,iX1,iX2,iX3,iCF_Ne      ), &
                 uGF(iNX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 uGF(iNX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 uGF(iNX,iX1,iX2,iX3,iGF_Gm_dd_33), &
                 uAF(iNX,iX1,iX2,iX3,iAF_P       ) )

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

  END SUBROUTINE InitializeMetric_Euler


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
                 * SUM( WeightsX_q * X1q**2 * P(:,iX1,1,1,iPF_D) )

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
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, Ustar, M, GS )

    INTEGER , INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout) :: &
      G    (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
      Ustar(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
      M    (1:,iX_B0(1)-swX_GS(1):, &
               iX_B0(2)-swX_GS(2):, &
               iX_B0(3)-swX_GS(3):,1:), &
      GS   (1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)

    CALL ComputeConformalFactorSourcesAndMg_XCFC_Euler &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, Ustar, GS )

    CALL ComputeConformalFactor_XCFC &
           ( iX_B0, iX_E0, iX_B1, iX_E1, GS, M )

    CALL UpdateConformalFactorAndMetric_XCFC &
           ( iX_B0, iX_E0, iX_B1, iX_E1, M, G )

  END SUBROUTINE ComputeConformalFactor


  SUBROUTINE ComputeLapseShiftCurvature &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, Ustar, M, GS )

    INTEGER , INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout) :: &
      G    (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
      Ustar(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
      M    (1:,iX_B0(1)-swX_GS(1):, &
               iX_B0(2)-swX_GS(2):, &
               iX_B0(3)-swX_GS(3):,1:), &
      GS   (1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)

    CALL ComputePressureTensorTrace_XCFC_Euler &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, Ustar, GS )

    CALL ComputeLapseShiftCurvature_XCFC &
           ( iX_B0, iX_E0, iX_B1, iX_E1, GS, M )

    CALL UpdateLapseShiftCurvature_XCFC &
           ( iX_B0, iX_E0, iX_B1, iX_E1, M, G )

    CALL ApplyBoundaryConditions_Geometry_XCFC &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G )

  END SUBROUTINE ComputeLapseShiftCurvature


END MODULE Euler_XCFC_UtilitiesModule
