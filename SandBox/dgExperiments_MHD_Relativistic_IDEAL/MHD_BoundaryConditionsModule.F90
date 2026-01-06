MODULE MHD_BoundaryConditionsModule

  USE KindModule, ONLY: &
    DP, &
    Zero, &
    One
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    bcX, &
    swX, &
    nDOFX, &
    nNodesX
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE InputOutputUtilitiesModule, ONLY: &
    FromField3D
  USE UtilitiesModule, ONLY: &
    NodeNumberX, &
    Locate, &
    Interpolate1D_Linear
  USE MagnetofluidFieldsModule, ONLY: &
    nCM, &
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
    nPM, &
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
    nAM, &
    iAM_P, &
    nDM, &
    iDM_IC_D, &
    iDM_IC_S1, &
    iDM_IC_S2, &
    iDM_IC_S3, &
    iDM_IC_E, &
    iDM_IC_Ne, &
    iDM_IC_B1, &
    iDM_IC_B2, &
    iDM_IC_B3, &
    iDM_IC_Chi
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX_FromScaleFactors
  USE GeometryFieldsModule, ONLY: &
    nGF, &
    iGF_h_1, &
    iGF_h_2, &
    iGF_h_3, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_Alpha, &
    iGF_Beta_1, &
    iGF_Beta_2, &
    iGF_Beta_3, &
    iGF_Psi, &
    iGF_SqrtGm
  USE UnitsModule, ONLY: &
    Gram, &
    Centimeter, &
    Erg, &
    Second, &
    Gauss
  USE EquationOfStateModule_IDEAL, ONLY: &
    Gamma_IDEAL, &
    ComputePressureFromPrimitive_IDEAL
  USE MHD_UtilitiesModule_Relativistic, ONLY: &
    ComputeConserved_MHD_Relativistic, &
    ComputePrimitive_MHD_Relativistic
  USE TimersModule_MHD, ONLY: &
    TimersStart_MHD, &
    TimersStop_MHD, &
    Timer_MHD_BoundaryConditions, &
    Timer_MHD_BC_ApplyBC, &
    Timer_MHD_BC_CopyIn, &
    Timer_MHD_BC_CopyOut
  USE QuadratureModule, ONLY: &
    GetQuadrature

  USE HDF5

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ApplyBoundaryConditions_MHD
  PUBLIC :: ApplyInnerBC_MHD
  PUBLIC :: ApplyOuterBC_MHD

  INTEGER, PARAMETER, PUBLIC :: iApplyBC_MHD_Both  = 0
  INTEGER, PARAMETER, PUBLIC :: iApplyBC_MHD_Inner = 1
  INTEGER, PARAMETER, PUBLIC :: iApplyBC_MHD_Outer = 2
  INTEGER, PARAMETER, PUBLIC :: iApplyBC_MHD_None  = 3

  REAL(DP), PUBLIC :: ExpD
  REAL(DP), PUBLIC :: ExpE

  CHARACTER(9),  PARAMETER :: &
    OutputDirectory    = '../Output'
  CHARACTER(18), PARAMETER :: &
    MagnetofluidSuffix = 'MagnetofluidFields'
  INTEGER :: FileNumber = 0
  INTEGER :: HDFERR

CONTAINS


  LOGICAL FUNCTION ApplyInnerBC_MHD( iApplyBC )

    INTEGER, INTENT(in) :: iApplyBC

    ApplyInnerBC_MHD = .FALSE.
    IF( iApplyBC .EQ. iApplyBC_MHD_Inner .OR. &
        iApplyBC .EQ. iApplyBC_MHD_Both ) &
    ApplyInnerBC_MHD = .TRUE.

  END FUNCTION ApplyInnerBC_MHD


  LOGICAL FUNCTION ApplyOuterBC_MHD( iApplyBC )

    INTEGER, INTENT(in) :: iApplyBC

    ApplyOuterBC_MHD = .FALSE.
    IF( iApplyBC .EQ. iApplyBC_MHD_Outer .OR. &
        iApplyBC .EQ. iApplyBC_MHD_Both ) &
    ApplyOuterBC_MHD = .TRUE.

  END FUNCTION ApplyOuterBC_MHD


  SUBROUTINE ApplyBoundaryConditions_MHD &
    ( t, iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, iApplyBC_Option )

    REAL(DP), INTENT(in)           :: t
    INTEGER,  INTENT(in)           :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)           :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(in)           :: &
      D(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout)        :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    INTEGER,  INTENT(in), OPTIONAL :: &
      iApplyBC_Option(3)

    INTEGER :: iApplyBC(3)

    CALL TimersStart_MHD( Timer_MHD_BoundaryConditions)

    iApplyBC = iApplyBC_MHD_Both

    IF( PRESENT( iApplyBC_Option ) ) &
      iApplyBC = iApplyBC_Option

    CALL TimersStart_MHD( Timer_MHD_BC_CopyIn )

    CALL TimersStop_MHD( Timer_MHD_BC_CopyIn )

    CALL TimersStart_MHD( Timer_MHD_BC_ApplyBC )

    CALL ApplyBC_MHD_X1( t, iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, iApplyBC(1) )

    CALL ApplyBC_MHD_X2( t, iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, iApplyBC(2) )

    CALL ApplyBC_MHD_X3( t, iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, iApplyBC(3) )

    CALL TimersStop_MHD( Timer_MHD_BC_ApplyBC )

    CALL TimersStart_MHD( Timer_MHD_BC_CopyOut )

    CALL TimersStop_MHD( Timer_MHD_BC_CopyOut )

    CALL TimersStop_MHD( Timer_MHD_BoundaryConditions )

  END SUBROUTINE ApplyBoundaryConditions_MHD


  SUBROUTINE ApplyBC_MHD_X1( t, iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, iApplyBC )

    REAL(DP), INTENT(in)    :: t
    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), &
      iApplyBC
    REAL(DP), INTENT(in)    :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(in)    :: &
      D(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    CHARACTER(6)   :: FileNumberString
    CHARACTER(256) :: FileName
    CHARACTER(256) :: GroupName
    CHARACTER(256) :: DatasetName

    INTEGER(HID_T) :: FILE_ID

    REAL(DP), ALLOCATABLE :: PressureArr(:), DensityArr(:), V3Arr(:), &
                             AlphaArr(:), PsiArr(:), X1Arr(:)

    REAL(DP) :: P(1:nDOFX,iX_B1(1):iX_E1(1), &
                          iX_B1(2):iX_E1(2), &
                          iX_B1(3):iX_E1(3),1:nPM)

    REAL(DP) :: A(1:nDOFX,iX_B1(1):iX_E1(1), &
                          iX_B1(2):iX_E1(2), &
                          iX_B1(3):iX_E1(3),1:nAM)

    REAL(DP) :: CD_I(1:nDOFX,iX_B1(1):iX_E1(1), &
                             iX_B1(2):iX_E1(2), &
                             iX_B1(3):iX_E1(3))

    REAL(DP) :: CS1_I(1:nDOFX,iX_B1(1):iX_E1(1), &
                              iX_B1(2):iX_E1(2), &
                              iX_B1(3):iX_E1(3))

    REAL(DP) :: CS2_I(1:nDOFX,iX_B1(1):iX_E1(1), &
                              iX_B1(2):iX_E1(2), &
                              iX_B1(3):iX_E1(3))

    REAL(DP) :: CS3_I(1:nDOFX,iX_B1(1):iX_E1(1), &
                              iX_B1(2):iX_E1(2), &
                              iX_B1(3):iX_E1(3))

    REAL(DP) :: CE_I(1:nDOFX,iX_B1(1):iX_E1(1), &
                             iX_B1(2):iX_E1(2), &
                             iX_B1(3):iX_E1(3))

    REAL(DP) :: CNe_I(1:nDOFX,iX_B1(1):iX_E1(1), &
                              iX_B1(2):iX_E1(2), &
                              iX_B1(3):iX_E1(3))

    REAL(DP) :: CB1_I(1:nDOFX,iX_B1(1):iX_E1(1), &
                              iX_B1(2):iX_E1(2), &
                              iX_B1(3):iX_E1(3))

    REAL(DP) :: CB2_I(1:nDOFX,iX_B1(1):iX_E1(1), &
                              iX_B1(2):iX_E1(2), &
                              iX_B1(3):iX_E1(3))

    REAL(DP) :: CB3_I(1:nDOFX,iX_B1(1):iX_E1(1), &
                              iX_B1(2):iX_E1(2), &
                              iX_B1(3):iX_E1(3))

    REAL(DP) :: CChi_I(1:nDOFX,iX_B1(1):iX_E1(1), &
                               iX_B1(2):iX_E1(2), &
                               iX_B1(3):iX_E1(3))

    REAL(DP) :: Dataset3D(iX_B1(1):iX_E1(1)*nNodesX(1), &
                          iX_B1(2):iX_E1(2)*nNodesX(2), &
                          iX_B1(3):iX_E1(3)*nNodesX(3))

    INTEGER  :: nX
    INTEGER  :: iCM, iX1, iX2, iX3
    INTEGER  :: iNX, iNX_L, iNX_R, iNX_0
    INTEGER  :: iNX1, iNX1_L, iNX1_R, &
                iNX2, iNX2_L, iNX2_R, &
                iNX3, iNX3_L, iNX3_R, &
                jNX, jNX1
    REAL(DP) :: X1, X2, X3, &
                X1_L, X1_R, X1_L_C, X1_R_C, X1_I, X1_O
    REAL(DP) :: dX1_L, dX1_R
    REAL(DP) :: V1, V2, V3, VSq, W
    REAL(DP) :: CB1, CB2, CB3, VdotB
    REAL(DP) :: D_0, E_0, R_0, R_q
    REAL(DP), DIMENSION(nNodesX(1)) :: xQ1_L, wQ1_L, xQ1_R, wQ1_R
    REAL(DP), DIMENSION(nNodesX(2)) :: xQ2_L, wQ2_L, xQ2_R, wQ2_R
    REAL(DP), DIMENSION(nNodesX(3)) :: xQ3_L, wQ3_L, xQ3_R, wQ3_R

    SELECT CASE ( bcX(1) )

    CASE ( 0 ) ! No Boundary Condition

    CASE ( 1 ) ! Periodic

#ifndef THORNADO_USE_AMREX

      ! --- Inner Boundary ---

      IF( ApplyInnerBC_MHD( iApplyBC ) )THEN

        DO iCM = 1, nCM
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = 1, swX(1)
        DO iNX = 1, nDOFX

          U(iNX,iX_B0(1)-iX1,iX2,iX3,iCM) &
            = U(iNX,iX_E0(1)-(iX1-1),iX2,iX3,iCM)

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

      ! --- Outer Boundary ---

      IF( ApplyOuterBC_MHD( iApplyBC ) )THEN

        DO iCM = 1, nCM
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = 1, swX(1)
        DO iNX = 1, nDOFX

          U(iNX,iX_E0(1)+iX1,iX2,iX3,iCM) &
            = U(iNX,iX_B0(1)+(iX1-1),iX2,iX3,iCM)

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

#endif

    CASE ( 2 ) ! Homogeneous

      ! --- Inner Boundary ---

      IF( ApplyInnerBC_MHD( iApplyBC ) )THEN

        DO iCM = 1, nCM
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = 1, swX(1)
        DO iNX = 1, nDOFX

          U(iNX,iX_B0(1)-iX1,iX2,iX3,iCM) &
            = U(iNX,iX_B0(1),iX2,iX3,iCM)

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

      ! --- Outer Boundary ---

      IF( ApplyOuterBC_MHD( iApplyBC ) )THEN

        DO iCM = 1, nCM
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = 1, swX(1)
        DO iNX = 1, nDOFX

          U(iNX,iX_E0(1)+iX1,iX2,iX3,iCM) &
            = U(iNX,iX_E0(1),iX2,iX3,iCM)

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

    CASE ( 3 ) ! Reflecting

      ! --- Inner Boundary ---

      IF( ApplyInnerBC_MHD( iApplyBC ) )THEN

        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = 1, swX(1)

          DO iNX3 = 1, nNodesX(3)
          DO iNX2 = 1, nNodesX(2)
          DO iNX1 = 1, nNodesX(1)

            jNX1 = ( nNodesX(1) - iNX1 ) + 1

            iNX = NodeNumberX( iNX1, iNX2, iNX3 )
            jNX = NodeNumberX( jNX1, iNX2, iNX3 )

            U(iNX,iX_B0(1)-iX1,iX2,iX3,iCM_D) &
              = + U(jNX,iX_B0(1),iX2,iX3,iCM_D)
            U(iNX,iX_B0(1)-iX1,iX2,iX3,iCM_S1) &
              = - U(jNX,iX_B0(1),iX2,iX3,iCM_S1)
            U(iNX,iX_B0(1)-iX1,iX2,iX3,iCM_S2) &
              = + U(jNX,iX_B0(1),iX2,iX3,iCM_S2)
            U(iNX,iX_B0(1)-iX1,iX2,iX3,iCM_S3) &
              = + U(jNX,iX_B0(1),iX2,iX3,iCM_S3)
            U(iNX,iX_B0(1)-iX1,iX2,iX3,iCM_E) &
              = + U(jNX,iX_B0(1),iX2,iX3,iCM_E)
            U(iNX,iX_B0(1)-iX1,iX2,iX3,iCM_Ne) &
              = + U(jNX,iX_B0(1),iX2,iX3,iCM_Ne)

          END DO
          END DO
          END DO

        END DO
        END DO
        END DO

      END IF

      ! --- Outer Boundary ---

      IF( ApplyOuterBC_MHD( iApplyBC ) )THEN

        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = 1, swX(1)

          DO iNX3 = 1, nNodesX(3)
          DO iNX2 = 1, nNodesX(2)
          DO iNX1 = 1, nNodesX(1)

            jNX1 = ( nNodesX(1) - iNX1 ) + 1

            iNX = NodeNumberX( iNX1, iNX2, iNX3 )
            jNX = NodeNumberX( jNX1, iNX2, iNX3 )

            U(iNX,iX_E0(1)+iX1,iX2,iX3,iCM_D) &
              = + U(jNX,iX_E0(1),iX2,iX3,iCM_D)
            U(iNX,iX_E0(1)+iX1,iX2,iX3,iCM_S1) &
              = - U(jNX,iX_E0(1),iX2,iX3,iCM_S1)
            U(iNX,iX_E0(1)+iX1,iX2,iX3,iCM_S2) &
              = + U(jNX,iX_E0(1),iX2,iX3,iCM_S2)
            U(iNX,iX_E0(1)+iX1,iX2,iX3,iCM_S3) &
              = + U(jNX,iX_E0(1),iX2,iX3,iCM_S3)
            U(iNX,iX_E0(1)+iX1,iX2,iX3,iCM_E) &
              = + U(jNX,iX_E0(1),iX2,iX3,iCM_E)
            U(iNX,iX_E0(1)+iX1,iX2,iX3,iCM_Ne) &
              = + U(jNX,iX_E0(1),iX2,iX3,iCM_Ne)

          END DO
          END DO
          END DO

        END DO
        END DO
        END DO

      END IF

    CASE ( 30 ) ! Reflecting (Inner), Zero (Outer)

      IF( ApplyInnerBC_MHD( iApplyBC ) )THEN

        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = 1, swX(1)

            DO iNX3 = 1, nNodesX(3)
            DO iNX2 = 1, nNodesX(2)
            DO iNX1 = 1, nNodesX(1)

              jNX1 = ( nNodesX(1) - iNX1 ) + 1

              iNX = NodeNumberX( iNX1, iNX2, iNX3 )
              jNX = NodeNumberX( jNX1, iNX2, iNX3 )

              U(iNX,iX_B0(1)-iX1,iX2,iX3,iCM_D) &
                = + U(jNX,iX_B0(1),iX2,iX3,iCM_D)
              U(iNX,iX_B0(1)-iX1,iX2,iX3,iCM_S1) &
                = - U(jNX,iX_B0(1),iX2,iX3,iCM_S1)
              U(iNX,iX_B0(1)-iX1,iX2,iX3,iCM_S2) &
                = + U(jNX,iX_B0(1),iX2,iX3,iCM_S2)
              U(iNX,iX_B0(1)-iX1,iX2,iX3,iCM_S3) &
                = + U(jNX,iX_B0(1),iX2,iX3,iCM_S3)
              U(iNX,iX_B0(1)-iX1,iX2,iX3,iCM_E) &
                = + U(jNX,iX_B0(1),iX2,iX3,iCM_E)
              U(iNX,iX_B0(1)-iX1,iX2,iX3,iCM_Ne) &
                = + U(jNX,iX_B0(1),iX2,iX3,iCM_Ne)

          END DO
          END DO
          END DO

        END DO
        END DO
        END DO

      END IF

    CASE ( 31 ) ! Reflecting (Inner), Homogeneous (Outer)

      ! --- Inner Boundary ---

      IF( ApplyInnerBC_MHD( iApplyBC ) )THEN

        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = 1, swX(1)

          DO iNX3 = 1, nNodesX(3)
          DO iNX2 = 1, nNodesX(2)
          DO iNX1 = 1, nNodesX(1)

            jNX1 = ( nNodesX(1) - iNX1 ) + 1

            iNX = NodeNumberX( iNX1, iNX2, iNX3 )
            jNX = NodeNumberX( jNX1, iNX2, iNX3 )

            U(iNX,iX_B0(1)-iX1,iX2,iX3,iCM_D) &
              = + U(jNX,iX_B0(1),iX2,iX3,iCM_D)
            U(iNX,iX_B0(1)-iX1,iX2,iX3,iCM_S1) &
              = - U(jNX,iX_B0(1),iX2,iX3,iCM_S1)
            U(iNX,iX_B0(1)-iX1,iX2,iX3,iCM_S2) &
              = + U(jNX,iX_B0(1),iX2,iX3,iCM_S2)
            U(iNX,iX_B0(1)-iX1,iX2,iX3,iCM_S3) &
              = + U(jNX,iX_B0(1),iX2,iX3,iCM_S3)
            U(iNX,iX_B0(1)-iX1,iX2,iX3,iCM_E) &
              = + U(jNX,iX_B0(1),iX2,iX3,iCM_E)
            U(iNX,iX_B0(1)-iX1,iX2,iX3,iCM_Ne) &
              = + U(jNX,iX_B0(1),iX2,iX3,iCM_Ne)

          END DO
          END DO
          END DO

        END DO
        END DO
        END DO

      END IF

     ! --- Outer Boundary ---

     IF( ApplyOuterBC_MHD( iApplyBC ) )THEN

        DO iCM = 1, nCM
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = 1, swX(1)
        DO iNX = 1, nDOFX

            U(iNX,iX_E0(1)+iX1,iX2,iX3,iCM) &
              = U(iNX,iX_E0(1),iX2,iX3,iCM)

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

    CASE ( 23 ) ! Homogeneous (Inner), Reflecting (Outer)

     ! --- Inner Boundary ---

     IF( ApplyInnerBC_MHD( iApplyBC ) )THEN

        DO iCM = 1, nCM
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = 1, swX(1)
        DO iNX = 1, nDOFX

            U(iNX,iX_B0(1)-iX1,iX2,iX3,iCM) &
              = U(iNX,iX_B0(1),iX2,iX3,iCM)

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

      ! --- Outer Boundary ---

      IF( ApplyOuterBC_MHD( iApplyBC ) )THEN

        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = 1, swX(1)

          DO iNX3 = 1, nNodesX(3)
          DO iNX2 = 1, nNodesX(2)
          DO iNX1 = 1, nNodesX(1)

            jNX1 = ( nNodesX(1) - iNX1 ) + 1

            iNX = NodeNumberX( iNX1, iNX2, iNX3 )
            jNX = NodeNumberX( jNX1, iNX2, iNX3 )

            U(iNX,iX_E0(1)+iX1,iX2,iX3,iCM_D) &
              = + U(jNX,iX_E0(1),iX2,iX3,iCM_D)
            U(iNX,iX_E0(1)+iX1,iX2,iX3,iCM_S1) &
              = - U(jNX,iX_E0(1),iX2,iX3,iCM_S1)
            U(iNX,iX_E0(1)+iX1,iX2,iX3,iCM_S2) &
              = + U(jNX,iX_E0(1),iX2,iX3,iCM_S2)
            U(iNX,iX_E0(1)+iX1,iX2,iX3,iCM_S3) &
              = + U(jNX,iX_E0(1),iX2,iX3,iCM_S3)
            U(iNX,iX_E0(1)+iX1,iX2,iX3,iCM_E) &
              = + U(jNX,iX_E0(1),iX2,iX3,iCM_E)
            U(iNX,iX_E0(1)+iX1,iX2,iX3,iCM_Ne) &
              = + U(jNX,iX_E0(1),iX2,iX3,iCM_Ne)

          END DO
          END DO
          END DO

        END DO
        END DO
        END DO

      END IF

    CASE ( 11 ) ! Custom BCs for Accretion Problem

      ! --- Inner Boundary ---

      IF( ApplyInnerBC_MHD( iApplyBC ) )THEN

        ASSOCIATE( X1_C  => MeshX(1) % Center, &
                   dX1   => MeshX(1) % Width,  &
                   eta_q => MeshX(1) % Nodes )

        R_0 = X1_C(1) + dX1(1) * eta_q(1)

        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = 1, swX(1)

          DO iNX3 = 1, nNodesX(3)
          DO iNX2 = 1, nNodesX(2)
          DO iNX1 = 1, nNodesX(1)

            iNX   = NodeNumberX( iNX1, iNX2, iNX3 )
            iNX_0 = NodeNumberX( 1,       iNX2, iNX3 )

            D_0 = U(iNX_0,1,iX2,iX3,iCM_D)
            E_0 = U(iNX_0,1,iX2,iX3,iCM_E)

            R_q = NodeCoordinate &
                    ( X1_C( iX_B0(1) - iX1 ), dX1( iX_B0(1) - iX1 ), &
                      eta_q( iNX1 ) )

            U(iNX,iX_B0(1)-iX1,iX2,iX3,iCM_D) &
              = D_0 * ( R_0 / R_q )**3

            U(iNX,iX_B0(1)-iX1,iX2,iX3,iCM_E) &
              = E_0 * ( R_0 / R_q )**4

          END DO
          END DO
          END DO

        END DO
        END DO
        END DO

        END ASSOCIATE

      END IF

    CASE ( 100 ) ! Custom BCs for Relativistic SAS Problem

      ! --- Inner Boundary ---

      IF( ApplyInnerBC_MHD( iApplyBC ) )THEN

        ASSOCIATE( X1_C  => MeshX(1) % Center, &
                   dX1   => MeshX(1) % Width,  &
                   eta_q => MeshX(1) % Nodes )

        R_0 = X1_C(1) + dX1(1) * eta_q(1)

        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = 1, swX(1)

          DO iNX3 = 1, nNodesX(3)
          DO iNX2 = 1, nNodesX(2)
          DO iNX1 = 1, nNodesX(1)

            iNX   = NodeNumberX( iNX1, iNX2, iNX3 )
            iNX_0 = NodeNumberX( 1,       iNX2, iNX3 )

            D_0 = U(iNX_0,1,iX2,iX3,iCM_D)
            E_0 = U(iNX_0,1,iX2,iX3,iCM_E)

            R_q = NodeCoordinate &
                    ( X1_C( iX_B0(1) - iX1 ), dX1( iX_B0(1) - iX1 ), &
                      eta_q( iNX1 ) )

            U(iNX,iX_B0(1)-iX1,iX2,iX3,iCM_D) &
              = D_0 * ( R_0 / R_q )**( ExpD )

            U(iNX,iX_B0(1)-iX1,iX2,iX3,iCM_E) &
              = E_0 * ( R_0 / R_q )**( ExpE )

          END DO
          END DO
          END DO

        END DO
        END DO
        END DO

        END ASSOCIATE

      END IF

    CASE ( 12 ) ! No Boundary Condition (Inner), Homogeneous (Outer)

     ! --- Outer Boundary ---

     IF( ApplyOuterBC_MHD( iApplyBC ) )THEN

        DO iCM = 1, nCM
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = 1, swX(1)
        DO iNX = 1, nDOFX

            U(iNX,iX_E0(1)+iX1,iX2,iX3,iCM) &
              = U(iNX,iX_E0(1),iX2,iX3,iCM)

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

    CASE( 41 ) ! BCs for constant rotation, 'shearing disk'-like problem

      ASSOCIATE( X1_C => MeshX(1) % Center, &
                 dX1  => MeshX(1) % Width )

      IF( ApplyInnerBC_MHD( iApplyBC ) )THEN

        DO iCM   = 1, nCM
        DO iX3   = iX_B0(3), iX_E0(3)
        DO iX2   = iX_B0(2), iX_E0(2)
        DO iX1   = 1, swX(1)
        DO iNX_L = 1, nDOFX

          ! Only BCs for B^{3}, others are initial values.

          IF( ( iCM .EQ. 7 ) .OR. ( iCM .EQ. 9 ) )THEN

            ! Get necessary cell and node values for the ghost cell (L)
            ! and the first compute cell (R).

            iNX1_L = NodeNumberTableX(1,iNX_L)
            X1_L   = NodeCoordinate( MeshX(1), iX_B0(1)-iX1, iNX1_L )

            X1_R_C = X1_C(iX_B0(1))
            dX1_R  = dX1 (iX_B0(1))

            ! Get the quadrature weights and points for the first compute
            ! cell (R).

            CALL GetQuadrature( nNodesX(1), xQ1_R, wQ1_R, 'Gaussian' )
            CALL GetQuadrature( nNodesX(2), xQ2_R, wQ2_R, 'Gaussian' )
            CALL GetQuadrature( nNodesX(3), xQ3_R, wQ3_R, 'Gaussian' )

            U(iNX_L,iX_B0(1)-iX1,iX2,iX3,iCM) = Zero

            DO iNX_R = 1, nDOFX

              ! Get node numbers for the first compute cell (R).

              iNX1_R = NodeNumberTableX(1,iNX_R)
              iNX2_R = NodeNumberTableX(2,iNX_R)
              iNX3_R = NodeNumberTableX(3,iNX_R)

              ! Set each ghost cell (L) node to the cell average of
              ! r * B^{3} in the first compute cell (R).

              U(iNX_L,iX_B0(1)-iX1,iX2,iX3,iCM) &
                = U(iNX_L,iX_B0(1)-iX1,iX2,iX3,iCM) &
                  + ( One / X1_L ) * ( One / X1_R_C ) &
                    * wQ1_R(iNX1_R) * wQ2_R(iNX2_R) * wQ3_R(iNX3_R) &
                    * ( X1_R_C + dX1_R * xQ1_R(iNX1_R) )**2 &
                    * U(iNX_R,iX_B0(1),iX2,iX3,iCM)

            END DO

          END IF

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

      ! --- Outer Boundary ---

      IF( ApplyOuterBC_MHD( iApplyBC ) )THEN

        DO iCM   = 1, nCM
        DO iX3   = iX_B0(3), iX_E0(3)
        DO iX2   = iX_B0(2), iX_E0(2)
        DO iX1   = 1, swX(1)
        DO iNX_R = 1, nDOFX

           ! Only BCs for B^{3}, others are initial values.

          IF( ( iCM .EQ. 7 ) .OR. ( iCM .EQ. 9 ) )THEN

            ! Get necessary cell and node values for the ghost cell (R)
            ! and the first compute cell (L).

            iNX1_R = NodeNumberTableX(1,iNX_R)
            X1_R   = NodeCoordinate( MeshX(1), iX_E0(1)+iX1, iNX1_R )

            X1_L_C = X1_C(iX_E0(1))
            dX1_L  = dX1 (iX_E0(1))

            ! Get the quadrature weights and points for the last compute
            ! cell (L).

            CALL GetQuadrature( nNodesX(1), xQ1_L, wQ1_L, 'Gaussian' )
            CALL GetQuadrature( nNodesX(2), xQ2_L, wQ2_L, 'Gaussian' )
            CALL GetQuadrature( nNodesX(3), xQ3_L, wQ3_L, 'Gaussian' )

            U(iNX_R,iX_E0(1)+iX1,iX2,iX3,iCM) = Zero

            DO iNX_L = 1, nDOFX

              ! Get node numbers for the last compute cell (L).

              iNX1_L = NodeNumberTableX(1,iNX_L)
              iNX2_L = NodeNumberTableX(2,iNX_L)
              iNX3_L = NodeNumberTableX(3,iNX_L)

              ! Set each ghost cell (R) node to the cell average of
              ! r * B^{3} in the last compute cell (L).

              U(iNX_R,iX_E0(1)+iX1,iX2,iX3,iCM) &
                = U(iNX_R,iX_E0(1)+iX1,iX2,iX3,iCM) &
                  + ( One / X1_R ) * ( One / X1_L_C ) &
                    * wQ1_L(iNX1_L) * wQ2_L(iNX2_L) * wQ3_L(iNX3_L) &
                    * ( X1_L_C + dX1_L * xQ1_L(iNX1_L) )**2 &
                    * U(iNX_L,iX_E0(1),iX2,iX3,iCM)

            END DO

          END IF

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

      END ASSOCIATE

    CASE( 42 ) ! Periodic in all fields except cleaning field = 0.

      ! --- Inner Boundary ---

      IF( ApplyInnerBC_MHD( iApplyBC ) )THEN

        DO iCM = 1, nCM
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = 1, swX(1)
        DO iNX = 1, nDOFX

          IF( iCM .NE. iCM_Chi )THEN

            U(iNX,iX_B0(1)-iX1,iX2,iX3,iCM) &
              = U(iNX,iX_E0(1)-(iX1-1),iX2,iX3,iCM)

          ELSE

            U(iNX,iX_B0(1)-iX1,iX2,iX3,iCM) &
              = Zero

          END IF

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

      ! --- Outer Boundary ---

      IF( ApplyOuterBC_MHD( iApplyBC ) )THEN

        DO iCM = 1, nCM
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = 1, swX(1)
        DO iNX = 1, nDOFX

          IF( iCM .NE. iCM_Chi )THEN

            U(iNX,iX_E0(1)+iX1,iX2,iX3,iCM) &
              = U(iNX,iX_B0(1)+(iX1-1),iX2,iX3,iCM)

          ELSE

            U(iNX,iX_E0(1)+iX1,iX2,iX3,iCM) &
              = Zero

          END IF

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

    CASE( 44 ) ! Custom BCs for relativistic shearing disk.

      IF( t .NE. Zero )THEN

        ! --- Inner Boundary --

        IF( ApplyInnerBC_MHD( iApplyBC ) )THEN

          DO iX3 = iX_B0(3), iX_E0(3)
          DO iX2 = iX_B0(2), iX_E0(2)
          DO iX1 = 1, swX(1)
          DO iNX = 1, nDOFX

            !PRINT*, 'Inner Boundary Location:     ', iNX, iX_B0(1) - iX1,     iX2, iX3
            !PRINT*, 'Last Compute Cell Location:  ', iNX, iX_E0(1) - (iX1-1), iX2, iX3
            !PRINT*

            IF( .TRUE. )THEN

              !PRINT*, 'Radial velocities too low. Using initial density for inner boundary.'

              U(iNX,iX_B0(1)-iX1,iX2,iX3,1) &
                = D(iNX,iX_B0(1)-iX1,iX2,iX3,iDM_IC_D)           
 
              U(iNX,iX_B0(1)-iX1,iX2,iX3,2) &
                = D(iNX,iX_B0(1)-iX1,iX2,iX3,iDM_IC_S1)           
 
              U(iNX,iX_B0(1)-iX1,iX2,iX3,3) &
                = D(iNX,iX_B0(1)-iX1,iX2,iX3,iDM_IC_S2)           
 
              U(iNX,iX_B0(1)-iX1,iX2,iX3,4) &
                = D(iNX,iX_B0(1)-iX1,iX2,iX3,iDM_IC_S3)           
 
              U(iNX,iX_B0(1)-iX1,iX2,iX3,5) &
                = D(iNX,iX_B0(1)-iX1,iX2,iX3,iDM_IC_E)           
 
            ELSE

              U(iNX,iX_B0(1)-iX1,iX2,iX3,1) &
                = U(iNX,iX_E0(1)-(iX1-1),iX2,iX3,1) &
                  - D(iNX,iX_E0(1)-(iX1-1),iX2,iX3,iDM_IC_D) &
                  + D(iNX,iX_B0(1)-iX1,iX2,iX3,iDM_IC_D)
 
              U(iNX,iX_B0(1)-iX1,iX2,iX3,2) &
                = U(iNX,iX_E0(1)-(iX1-1),iX2,iX3,2) &
                  - D(iNX,iX_E0(1)-(iX1-1),iX2,iX3,iDM_IC_S1) &
                  + D(iNX,iX_B0(1)-iX1,iX2,iX3,iDM_IC_S1)
  
              U(iNX,iX_B0(1)-iX1,iX2,iX3,3) &
                = U(iNX,iX_E0(1)-(iX1-1),iX2,iX3,3) &
                  - D(iNX,iX_E0(1)-(iX1-1),iX2,iX3,iDM_IC_S2) &
                  + D(iNX,iX_B0(1)-iX1,iX2,iX3,iDM_IC_S2)
  
              U(iNX,iX_B0(1)-iX1,iX2,iX3,4) &
                = U(iNX,iX_E0(1)-(iX1-1),iX2,iX3,4) &
                  - D(iNX,iX_E0(1)-(iX1-1),iX2,iX3,iDM_IC_S3) &
                  + D(iNX,iX_B0(1)-iX1,iX2,iX3,iDM_IC_S3)
  
              U(iNX,iX_B0(1)-iX1,iX2,iX3,5) &
                = U(iNX,iX_E0(1)-(iX1-1),iX2,iX3,5) &
                  - D(iNX,iX_E0(1)-(iX1-1),iX2,iX3,iDM_IC_E) &
                  + D(iNX,iX_B0(1)-iX1,iX2,iX3,iDM_IC_E)
  
            END IF

            U(iNX,iX_B0(1)-iX1,iX2,iX3,7) &
              = U(iNX,iX_E0(1)-(iX1-1),iX2,iX3,7)

            U(iNX,iX_B0(1)-iX1,iX2,iX3,8) &
              = U(iNX,iX_E0(1)-(iX1-1),iX2,iX3,8)

            U(iNX,iX_B0(1)-iX1,iX2,iX3,9) &
              = U(iNX,iX_E0(1)-(iX1-1),iX2,iX3,9)

            U(iNX,iX_B0(1)-iX1,iX2,iX3,10) &
              = U(iNX,iX_E0(1)-(iX1-1),iX2,iX3,10)

          END DO
          END DO
          END DO
          END DO

        END IF

        ! --- Outer Boundary --

        IF( ApplyOuterBC_MHD( iApplyBC ) )THEN

          DO iX3 = iX_B0(3), iX_E0(3)
          DO iX2 = iX_B0(2), iX_E0(2)
          DO iX1 = 1, swX(1)
          DO iNX = 1, nDOFX

            CALL ComputePrimitive_MHD_Relativistic &
                   ( U(iNX,iX_E0(1)+iX1,iX2,iX3,1 ), &
                     U(iNX,iX_E0(1)+iX1,iX2,iX3,2 ), &
                     U(iNX,iX_E0(1)+iX1,iX2,iX3,3 ), &
                     U(iNX,iX_E0(1)+iX1,iX2,iX3,4 ), &
                     U(iNX,iX_E0(1)+iX1,iX2,iX3,5 ), &
                     U(iNX,iX_E0(1)+iX1,iX2,iX3,6 ), &
                     U(iNX,iX_E0(1)+iX1,iX2,iX3,7 ), &
                     U(iNX,iX_E0(1)+iX1,iX2,iX3,8 ), &
                     U(iNX,iX_E0(1)+iX1,iX2,iX3,9 ), &
                     U(iNX,iX_E0(1)+iX1,iX2,iX3,10), &
                     P(iNX,iX_E0(1)+iX1,iX2,iX3,1 ), &
                     P(iNX,iX_E0(1)+iX1,iX2,iX3,2 ), &
                     P(iNX,iX_E0(1)+iX1,iX2,iX3,3 ), &
                     P(iNX,iX_E0(1)+iX1,iX2,iX3,4 ), &
                     P(iNX,iX_E0(1)+iX1,iX2,iX3,5 ), &
                     P(iNX,iX_E0(1)+iX1,iX2,iX3,6 ), &
                     P(iNX,iX_E0(1)+iX1,iX2,iX3,7 ), &
                     P(iNX,iX_E0(1)+iX1,iX2,iX3,8 ), &
                     P(iNX,iX_E0(1)+iX1,iX2,iX3,9 ), &
                     P(iNX,iX_E0(1)+iX1,iX2,iX3,10), &
                     G(iNX,iX_E0(1)+iX1,iX2,iX3,iGF_Gm_dd_11), &
                     G(iNX,iX_E0(1)+iX1,iX2,iX3,iGF_Gm_dd_22), &
                     G(iNX,iX_E0(1)+iX1,iX2,iX3,iGF_Gm_dd_33), &
                     G(iNX,iX_E0(1)+iX1,iX2,iX3,iGF_Alpha   ), &
                     G(iNX,iX_E0(1)+iX1,iX2,iX3,iGF_Beta_1  ), &
                     G(iNX,iX_E0(1)+iX1,iX2,iX3,iGF_Beta_2  ), &
                     G(iNX,iX_E0(1)+iX1,iX2,iX3,iGF_Beta_3  ), &
                     .FALSE. )

            CALL ComputePrimitive_MHD_Relativistic &
                   ( U(iNX,iX_B0(1)+(iX1-1),iX2,iX3,1 ), &
                     U(iNX,iX_B0(1)+(iX1-1),iX2,iX3,2 ), &
                     U(iNX,iX_B0(1)+(iX1-1),iX2,iX3,3 ), &
                     U(iNX,iX_B0(1)+(iX1-1),iX2,iX3,4 ), &
                     U(iNX,iX_B0(1)+(iX1-1),iX2,iX3,5 ), &
                     U(iNX,iX_B0(1)+(iX1-1),iX2,iX3,6 ), &
                     U(iNX,iX_B0(1)+(iX1-1),iX2,iX3,7 ), &
                     U(iNX,iX_B0(1)+(iX1-1),iX2,iX3,8 ), &
                     U(iNX,iX_B0(1)+(iX1-1),iX2,iX3,9 ), &
                     U(iNX,iX_B0(1)+(iX1-1),iX2,iX3,10), &
                     P(iNX,iX_B0(1)+(iX1-1),iX2,iX3,1 ), &
                     P(iNX,iX_B0(1)+(iX1-1),iX2,iX3,2 ), &
                     P(iNX,iX_B0(1)+(iX1-1),iX2,iX3,3 ), &
                     P(iNX,iX_B0(1)+(iX1-1),iX2,iX3,4 ), &
                     P(iNX,iX_B0(1)+(iX1-1),iX2,iX3,5 ), &
                     P(iNX,iX_B0(1)+(iX1-1),iX2,iX3,6 ), &
                     P(iNX,iX_B0(1)+(iX1-1),iX2,iX3,7 ), &
                     P(iNX,iX_B0(1)+(iX1-1),iX2,iX3,8 ), &
                     P(iNX,iX_B0(1)+(iX1-1),iX2,iX3,9 ), &
                     P(iNX,iX_B0(1)+(iX1-1),iX2,iX3,10), &
                     G(iNX,iX_B0(1)+(iX1-1),iX2,iX3,iGF_Gm_dd_11), &
                     G(iNX,iX_B0(1)+(iX1-1),iX2,iX3,iGF_Gm_dd_22), &
                     G(iNX,iX_B0(1)+(iX1-1),iX2,iX3,iGF_Gm_dd_33), &
                     G(iNX,iX_B0(1)+(iX1-1),iX2,iX3,iGF_Alpha   ), &
                     G(iNX,iX_B0(1)+(iX1-1),iX2,iX3,iGF_Beta_1  ), &
                     G(iNX,iX_B0(1)+(iX1-1),iX2,iX3,iGF_Beta_2  ), &
                     G(iNX,iX_B0(1)+(iX1-1),iX2,iX3,iGF_Beta_3  ), &
                     .FALSE. )

            !PRINT*, 'Outer Boundary Location:     ', iNX, iX_E0(1) + iX1,     iX2, iX3
            !PRINT*, 'First Compute Cell Location: ', iNX, iX_B0(1) + (iX1-1), iX2, iX3
            !PRINT*

            IF( .TRUE. )THEN

              !PRINT*, 'Radial velocities too low. Using initial density for outer boundary.'

              U(iNX,iX_E0(1)+iX1,iX2,iX3,1) &
                = D(iNX,iX_E0(1)+iX1,iX2,iX3,iDM_IC_D)

              U(iNX,iX_E0(1)+iX1,iX2,iX3,2) &
                = D(iNX,iX_E0(1)+iX1,iX2,iX3,iDM_IC_S1)

              U(iNX,iX_E0(1)+iX1,iX2,iX3,3) &
                = D(iNX,iX_E0(1)+iX1,iX2,iX3,iDM_IC_S2)

              U(iNX,iX_E0(1)+iX1,iX2,iX3,4) &
                = D(iNX,iX_E0(1)+iX1,iX2,iX3,iDM_IC_S3)
 
              U(iNX,iX_E0(1)+iX1,iX2,iX3,5) &
                = D(iNX,iX_E0(1)+iX1,iX2,iX3,iDM_IC_E)
 
            ELSE

              U(iNX,iX_E0(1)+iX1,iX2,iX3,1) &
                = U(iNX,iX_B0(1)+(iX1-1),iX2,iX3,1) &
                  - D(iNX,iX_B0(1)+(iX1-1),iX2,iX3,iDM_IC_D) &
                  + D(iNX,iX_E0(1)+iX1,iX2,iX3,iDM_IC_D)
 
              U(iNX,iX_E0(1)+iX1,iX2,iX3,2) &
                = U(iNX,iX_B0(1)+(iX1-1),iX2,iX3,2) &
                  - D(iNX,iX_B0(1)+(iX1-1),iX2,iX3,iDM_IC_S1) &
                  + D(iNX,iX_E0(1)+iX1,iX2,iX3,iDM_IC_S1)
  
              U(iNX,iX_E0(1)+iX1,iX2,iX3,3) &
                = U(iNX,iX_B0(1)+(iX1-1),iX2,iX3,3) &
                  - D(iNX,iX_B0(1)+(iX1-1),iX2,iX3,iDM_IC_S2) &
                  + D(iNX,iX_E0(1)+iX1,iX2,iX3,iDM_IC_S2)
  
              U(iNX,iX_E0(1)+iX1,iX2,iX3,4) &
                = U(iNX,iX_B0(1)+(iX1-1),iX2,iX3,4) &
                  - D(iNX,iX_B0(1)+(iX1-1),iX2,iX3,iDM_IC_S3) &
                  + D(iNX,iX_E0(1)+iX1,iX2,iX3,iDM_IC_S3)
  
              U(iNX,iX_E0(1)+iX1,iX2,iX3,5) &
                = U(iNX,iX_B0(1)+(iX1-1),iX2,iX3,5) &
                  - D(iNX,iX_B0(1)+(iX1-1),iX2,iX3,iDM_IC_E) &
                  + D(iNX,iX_E0(1)+iX1,iX2,iX3,iDM_IC_E)
 
            END IF
 
              U(iNX,iX_E0(1)+iX1,iX2,iX3,7) &
                = U(iNX,iX_B0(1)+(iX1-1),iX2,iX3,7)
  
              U(iNX,iX_E0(1)+iX1,iX2,iX3,8) &
                = U(iNX,iX_B0(1)+(iX1-1),iX2,iX3,8)
  
              U(iNX,iX_E0(1)+iX1,iX2,iX3,9) &
                = U(iNX,iX_B0(1)+(iX1-1),iX2,iX3,9)
  
              U(iNX,iX_E0(1)+iX1,iX2,iX3,10) &
                = U(iNX,iX_B0(1)+(iX1-1),iX2,iX3,10)
  
          END DO
          END DO
          END DO
          END DO

        END IF

      ELSE

        RETURN

      END IF

    CASE DEFAULT

    END SELECT

  END SUBROUTINE ApplyBC_MHD_X1


  SUBROUTINE ApplyBC_MHD_X2( t, iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, iApplyBC )

    REAL(DP), INTENT(in)    :: t
    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), &
      iApplyBC
    REAL(DP), INTENT(in)    :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(in)    :: &
      D(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:) 
    REAL(DP), INTENT(inout) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER :: iCM, iX1, iX2, iX3
    INTEGER :: iNX, iNX1, iNX2, iNX3, jNX, jNX2

    SELECT CASE ( bcX(2) )

    CASE ( 0 ) ! No Boundary Condition

    CASE ( 1 ) ! Periodic

#ifndef THORNADO_USE_AMREX

      ! --- Inner Boundary ---

      IF( ApplyInnerBC_MHD( iApplyBC ) )THEN

        DO iCM = 1, nCM
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = 1, swX(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1, nDOFX

          U(iNX,iX1,iX_B0(2)-iX2,iX3,iCM) &
            = U(iNX,iX1,iX_E0(2)-(iX2-1),iX3,iCM)

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

      ! --- Outer Boundary ---

      IF( ApplyOuterBC_MHD( iApplyBC ) )THEN

        DO iCM = 1, nCM
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = 1, swX(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1, nDOFX

          U(iNX,iX1,iX_E0(2)+iX2,iX3,iCM) &
            = U(iNX,iX1,iX_B0(2)+(iX2-1),iX3,iCM)

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

#endif

    CASE ( 2 ) ! Homogeneous

      ! --- Inner Boundary ---

      IF( ApplyInnerBC_MHD( iApplyBC ) )THEN

        DO iCM = 1, nCM
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = 1, swX(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1, nDOFX

          U(iNX,iX1,iX_B0(2)-iX2,iX3,iCM) &
            = U(iNX,iX1,iX_B0(2),iX3,iCM)

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

      ! --- Outer Boundary ---

      IF( ApplyOuterBC_MHD( iApplyBC ) )THEN

        DO iCM = 1, nCM
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = 1, swX(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1, nDOFX

          U(iNX,iX1,iX_E0(2)+iX2,iX3,iCM) &
            = U(iNX,iX1,iX_E0(2),iX3,iCM)

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

    CASE ( 3 ) ! Reflecting

      ! --- Inner Boundary ---

      IF( ApplyInnerBC_MHD( iApplyBC ) )THEN

        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = 1, swX(2)
        DO iX1 = iX_B0(1), iX_E0(1)

          DO iNX3 = 1, nNodesX(3)
          DO iNX2 = 1, nNodesX(2)
          DO iNX1 = 1, nNodesX(1)

            jNX2 = ( nNodesX(2) - iNX2 ) + 1

            iNX = NodeNumberX( iNX1, iNX2, iNX3 )
            jNX = NodeNumberX( iNX1, jNX2, iNX3 )

            U(iNX,iX1,iX_B0(2)-iX2,iX3,iCM_D) &
              = + U(jNX,iX1,iX_B0(2),iX3,iCM_D)
            U(iNX,iX1,iX_B0(2)-iX2,iX3,iCM_S1) &
              = + U(jNX,iX1,iX_B0(2),iX3,iCM_S1)
            U(iNX,iX1,iX_B0(2)-iX2,iX3,iCM_S2) &
              = - U(jNX,iX1,iX_B0(2),iX3,iCM_S2)
            U(iNX,iX1,iX_B0(2)-iX2,iX3,iCM_S3) &
              = + U(jNX,iX1,iX_B0(2),iX3,iCM_S3)
            U(iNX,iX1,iX_B0(2)-iX2,iX3,iCM_E) &
              = + U(jNX,iX1,iX_B0(2),iX3,iCM_E)
            U(iNX,iX1,iX_B0(2)-iX2,iX3,iCM_Ne) &
              = + U(jNX,iX1,iX_B0(2),iX3,iCM_Ne)

          END DO
          END DO
          END DO

        END DO
        END DO
        END DO

      END IF

      ! --- Outer Boundary ---

      IF( ApplyOuterBC_MHD( iApplyBC ) )THEN

        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = 1, swX(2)
        DO iX1 = iX_B0(1), iX_E0(1)

          DO iNX3 = 1, nNodesX(3)
          DO iNX2 = 1, nNodesX(2)
          DO iNX1 = 1, nNodesX(1)

            jNX2 = ( nNodesX(2) - iNX2 ) + 1

            iNX = NodeNumberX( iNX1, iNX2, iNX3 )
            jNX = NodeNumberX( iNX1, jNX2, iNX3 )

            U(iNX,iX1,iX_E0(2)+iX2,iX3,iCM_D) &
              = + U(jNX,iX1,iX_E0(2),iX3,iCM_D)
            U(iNX,iX1,iX_E0(2)+iX2,iX3,iCM_S1) &
              = + U(jNX,iX1,iX_E0(2),iX3,iCM_S1)
            U(iNX,iX1,iX_E0(2)+iX2,iX3,iCM_S2) &
              = - U(jNX,iX1,iX_E0(2),iX3,iCM_S2)
            U(iNX,iX1,iX_E0(2)+iX2,iX3,iCM_S3) &
              = + U(jNX,iX1,iX_E0(2),iX3,iCM_S3)
            U(iNX,iX1,iX_E0(2)+iX2,iX3,iCM_E) &
              = + U(jNX,iX1,iX_E0(2),iX3,iCM_E)
            U(iNX,iX1,iX_E0(2)+iX2,iX3,iCM_Ne) &
              = + U(jNX,iX1,iX_E0(2),iX3,iCM_Ne)

          END DO
          END DO
          END DO

        END DO
        END DO
        END DO

      END IF

    CASE ( 30 ) ! Reflecting (Inner), Zero (Outer)

      IF( ApplyInnerBC_MHD( iApplyBC ) )THEN

        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = 1, swX(2)
        DO iX1 = iX_B0(1), iX_E0(1)

          DO iNX3 = 1, nNodesX(3)
          DO iNX2 = 1, nNodesX(2)
          DO iNX1 = 1, nNodesX(1)

            jNX2 = ( nNodesX(2) - iNX2 ) + 1

            iNX = NodeNumberX( iNX1, iNX2, iNX3 )
            jNX = NodeNumberX( iNX1, jNX2, iNX3 )

            U(iNX,iX1,iX_B0(2)-iX2,iX3,iCM_D) &
              = + U(jNX,iX1,iX_B0(2),iX3,iCM_D)
            U(iNX,iX1,iX_B0(2)-iX2,iX3,iCM_S1) &
              = + U(jNX,iX1,iX_B0(2),iX3,iCM_S1)
            U(iNX,iX1,iX_B0(2)-iX2,iX3,iCM_S2) &
              = - U(jNX,iX1,iX_B0(2),iX3,iCM_S2)
            U(iNX,iX1,iX_B0(2)-iX2,iX3,iCM_S3) &
              = + U(jNX,iX1,iX_B0(2),iX3,iCM_S3)
            U(iNX,iX1,iX_B0(2)-iX2,iX3,iCM_E) &
              = + U(jNX,iX1,iX_B0(2),iX3,iCM_E)
            U(iNX,iX1,iX_B0(2)-iX2,iX3,iCM_Ne) &
              = + U(jNX,iX1,iX_B0(2),iX3,iCM_Ne)

          END DO
          END DO
          END DO

        END DO
        END DO
        END DO

      END IF

    CASE ( 31 ) ! Reflecting (Inner), Homogeneous (Outer)

      ! --- Inner Boundary ---

      IF( ApplyInnerBC_MHD( iApplyBC ) )THEN

        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = 1, swX(2)
        DO iX1 = iX_B0(1), iX_E0(1)

          DO iNX3 = 1, nNodesX(3)
          DO iNX2 = 1, nNodesX(2)
          DO iNX1 = 1, nNodesX(1)

            jNX2 = ( nNodesX(2) - iNX2 ) + 1

            iNX = NodeNumberX( iNX1, iNX2, iNX3 )
            jNX = NodeNumberX( iNX1, jNX2, iNX3 )

            U(iNX,iX1,iX_B0(2)-iX2,iX3,iCM_D) &
              = + U(jNX,iX1,iX_B0(2),iX3,iCM_D)
            U(iNX,iX1,iX_B0(2)-iX2,iX3,iCM_S1) &
              = + U(jNX,iX1,iX_B0(2),iX3,iCM_S1)
            U(iNX,iX1,iX_B0(2)-iX2,iX3,iCM_S2) &
              = - U(jNX,iX1,iX_B0(2),iX3,iCM_S2)
            U(iNX,iX1,iX_B0(2)-iX2,iX3,iCM_S3) &
              = + U(jNX,iX1,iX_B0(2),iX3,iCM_S3)
            U(iNX,iX1,iX_B0(2)-iX2,iX3,iCM_E) &
              = + U(jNX,iX1,iX_B0(2),iX3,iCM_E)
            U(iNX,iX1,iX_B0(2)-iX2,iX3,iCM_Ne) &
              = + U(jNX,iX1,iX_B0(2),iX3,iCM_Ne)

          END DO
          END DO
          END DO

        END DO
        END DO
        END DO

      END IF

      ! --- Outer Boundary ---

      IF( ApplyOuterBC_MHD( iApplyBC ) )THEN

        DO iCM = 1, nCM
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = 1, swX(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1, nDOFX

          U(iNX,iX1,iX_E0(2)+iX2,iX3,iCM) &
            = U(iNX,iX1,iX_E0(2),iX3,iCM)

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

    CASE ( 12 ) ! No Boundary Condition (Inner), Homogeneous (Outer)

      ! --- Outer Boundary ---

      IF( ApplyOuterBC_MHD( iApplyBC ) )THEN

        DO iCM = 1, nCM
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = 1, swX(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1, nDOFX

          U(iNX,iX1,iX_E0(2)+iX2,iX3,iCM) &
            = U(iNX,iX1,iX_E0(2),iX3,iCM)

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

    CASE( 42 ) ! Periodic in all fields except cleaning field = 0.

      ! --- Inner Boundary ---

      IF( ApplyInnerBC_MHD( iApplyBC ) )THEN

        DO iCM = 1, nCM
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = 1, swX(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1, nDOFX

          IF( iCM .NE. iCM_Chi )THEN

            U(iNX,iX1,iX_B0(2)-iX2,iX3,iCM) &
              = U(iNX,iX1,iX_E0(2)-(iX2-1),iX3,iCM)

          ELSE

            U(iNX,iX1,iX_B0(2)-iX2,iX3,iCM) &
              = Zero

          END IF

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

      ! --- Outer Boundary ---

      IF( ApplyOuterBC_MHD( iApplyBC ) )THEN

        DO iCM = 1, nCM
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = 1, swX(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1, nDOFX

          IF( iCM .NE. iCM_Chi )THEN

            U(iNX,iX1,iX_E0(2)+iX2,iX3,iCM) &
              = U(iNX,iX1,iX_B0(2)+(iX2-1),iX3,iCM)

          ELSE

            U(iNX,iX1,iX_E0(2)+iX2,iX3,iCM) &
              = Zero

          END IF

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

    CASE DEFAULT

    END SELECT

  END SUBROUTINE ApplyBC_MHD_X2


  SUBROUTINE ApplyBC_MHD_X3( t, iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, iApplyBC )

    REAL(DP), INTENT(in)    :: t
    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), &
      iApplyBC
    REAL(DP), INTENT(in)    :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(in)    :: &
      D(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER :: iCM, iX1, iX2, iX3
    INTEGER :: iNX

    SELECT CASE ( bcX(3) )

    CASE ( 0 ) ! No Boundary Condition

    CASE ( 1 ) ! Periodic

#ifndef THORNADO_USE_AMREX

      ! --- Inner Boundary ---

      IF( ApplyInnerBC_MHD( iApplyBC ) )THEN

        DO iCM = 1, nCM
        DO iX3 = 1, swX(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1, nDOFX

          U(iNX,iX1,iX2,iX_B0(3)-iX3,iCM) &
            = U(iNX,iX1,iX2,iX_E0(3)-(iX3-1),iCM)

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

      ! --- Outer Boundary ---

      IF( ApplyOuterBC_MHD( iApplyBC ) )THEN

        DO iCM = 1, nCM
        DO iX3 = 1, swX(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1, nDOFX

            U(iNX,iX1,iX2,iX_E0(3)+iX3,iCM) &
              = U(iNX,iX1,iX2,iX_B0(3)+(iX3-1),iCM)

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

#endif

    CASE ( 2 ) ! Homogeneous

      ! --- Inner Boundary ---

      IF( ApplyInnerBC_MHD( iApplyBC ) )THEN

        DO iCM = 1, nCM
        DO iX3 = 1, swX(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1, nDOFX

          U(iNX,iX1,iX2,iX_B0(3)-iX3,iCM) &
            = U(iNX,iX1,iX2,iX_B0(3),iCM)


        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

      ! --- Outer Boundary ---

      IF( ApplyOuterBC_MHD( iApplyBC ) )THEN

        DO iCM = 1, nCM
        DO iX3 = 1, swX(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1, nDOFX

          U(iNX,iX1,iX2,iX_E0(3)+iX3,iCM) &
            = U(iNX,iX1,iX2,iX_E0(3),iCM)

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

    CASE ( 12 ) ! No Boundary Condition (Inner), Homogeneous (Outer)

      ! --- Outer Boundary ---

      IF( ApplyOuterBC_MHD( iApplyBC ) )THEN

        DO iCM = 1, nCM
        DO iX3 = 1, swX(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1, nDOFX

          U(iNX,iX1,iX2,iX_E0(3)+iX3,iCM) &
            = U(iNX,iX1,iX2,iX_E0(3),iCM)

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

    CASE DEFAULT

    END SELECT

  END SUBROUTINE ApplyBC_MHD_X3


END MODULE MHD_BoundaryConditionsModule
