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
    iAM_P
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
    iGF_Psi
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
    ComputeConserved_MHD_Relativistic
  USE TimersModule_MHD, ONLY: &
    TimersStart_MHD, &
    TimersStop_MHD, &
    Timer_MHD_BoundaryConditions, &
    Timer_MHD_BC_ApplyBC, &
    Timer_MHD_BC_CopyIn, &
    Timer_MHD_BC_CopyOut

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
    ( t, iX_B0, iX_E0, iX_B1, iX_E1, U, iApplyBC_Option )

    REAL(DP), INTENT(in)           :: t
    INTEGER,  INTENT(in)           :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
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

    CALL ApplyBC_MHD_X1( t, iX_B0, iX_E0, iX_B1, iX_E1, U, iApplyBC(1) )

    CALL ApplyBC_MHD_X2( t, iX_B0, iX_E0, iX_B1, iX_E1, U, iApplyBC(2) )

    CALL ApplyBC_MHD_X3( t, iX_B0, iX_E0, iX_B1, iX_E1, U, iApplyBC(3) )

    CALL TimersStop_MHD( Timer_MHD_BC_ApplyBC )

    CALL TimersStart_MHD( Timer_MHD_BC_CopyOut )

    CALL TimersStop_MHD( Timer_MHD_BC_CopyOut )

    CALL TimersStop_MHD( Timer_MHD_BoundaryConditions )

  END SUBROUTINE ApplyBoundaryConditions_MHD


  SUBROUTINE ApplyBC_MHD_X1( t, iX_B0, iX_E0, iX_B1, iX_E1, U, iApplyBC )

    REAL(DP), INTENT(in)    :: t
    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), &
      iApplyBC
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

    REAL(DP) :: G(1:nDOFX,iX_B1(1):iX_E1(1), &
                          iX_B1(2):iX_E1(2), &
                          iX_B1(3):iX_E1(3),1:nGF)

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
    INTEGER  :: iNX, iNX_0
    INTEGER  :: iNX1, iNX2, iNX3, jNX, jNX1
    REAL(DP) :: X1, X2, X3
    REAL(DP) :: V1, V2, V3, VSq, W
    REAL(DP) :: CB1, CB2, CB3, VdotB
    REAL(DP) :: D_0, E_0, R_0, R_q

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

    CASE( 44 ) ! Custom BCs for relativistic shearing disk.

      IF( t .NE. Zero )THEN

        IF( ApplyOuterBC_MHD( iApplyBC ) .OR. ApplyInnerBC_MHD( iApplyBC ) )THEN

          FileName = "/home/jbuffal/thornado_MHD_3D/Workflow/MHD/ShearingDisk/GR_LR_diffrot.h5"

          ! --- Populate arrays ---

          CALL H5OPEN_F( HDFERR )

          CALL H5FOPEN_F( TRIM( FileName ), H5F_ACC_RDONLY_F, FILE_ID, HDFERR )

          nX = 10000

          ALLOCATE( PressureArr(nX), DensityArr(nX), V3Arr(nX), AlphaArr(nX), &
                    PsiArr(nX), X1Arr(nX) )

          CALL ReadDataset1DHDF( PsiArr,      '/psi',   FILE_ID )
          CALL ReadDataset1DHDF( AlphaArr,    '/alpha', FILE_ID )
          CALL ReadDataset1DHDF( X1Arr,       '/r',     FILE_ID )
          CALL ReadDataset1DHDF( PressureArr, '/pres',  FILE_ID )
          CALL ReadDataset1DHDF( DensityArr,  '/rho',   FILE_ID )
          CALL ReadDataset1DHDF( V3Arr,       '/V3',    FILE_ID )

          X1Arr       = X1Arr       * Centimeter
          DensityArr  = DensityArr  * ( Gram / Centimeter**3 )
          PressureArr = PressureArr * ( Erg  / Centimeter**3 )
          V3Arr       = V3Arr       * ( One  / Second )

          ! --- Map to inner boundary. ---

          IF( ApplyInnerBC_MHD( iApplyBC ) )THEN

            DO iX3 = iX_B0(3), iX_E0(3)
            DO iX2 = iX_B0(2), iX_E0(2)

            DO iX1 = 1, swX(1)

              DO iNX = 1, nDOFX

                iNX1 = NodeNumberTableX(1,iNX)
                iNX2 = NodeNumberTableX(2,iNX)
                iNX3 = NodeNumberTableX(3,iNX)

                X1 = NodeCoordinate( MeshX(1), iX_B0(1)-iX1, iNX1 )
                X2 = NodeCoordinate( MeshX(2), iX2, iNX2 )
                X3 = NodeCoordinate( MeshX(3), iX3, iNX3 )

                ! --- Geometry Fields ---

                G(iNX,iX_B0(1)-iX1,iX2,iX3,iGF_Alpha) &
                  = Interpolate1D( X1Arr, AlphaArr, SIZE( X1Arr ), X1 )

                G(iNX,iX_B0(1)-iX1,iX2,iX3,iGF_Psi) &
                  = Interpolate1D( X1Arr, PsiArr, SIZE( X1Arr ), X1 )

                G(iNX,iX_B0(1)-iX1,iX2,iX3,iGF_Beta_1) = Zero
                G(iNX,iX_B0(1)-iX1,iX2,iX3,iGF_Beta_2) = Zero
                G(iNX,iX_B0(1)-iX1,iX2,iX3,iGF_Beta_3) = Zero

                G(iNX,iX_B0(1)-iX1,iX2,iX3,iGF_h_1) &
                  = G(iNX,iX_B0(1)-iX1,iX2,iX3,iGF_Psi)**2
                G(iNX,iX_B0(1)-iX1,iX2,iX3,iGF_h_2) &
                  = G(iNX,iX_B0(1)-iX1,iX2,iX3,iGF_Psi)**2
                G(iNX,iX_B0(1)-iX1,iX2,iX3,iGF_h_3) &
                  = G(iNX,iX_B0(1)-iX1,iX2,iX3,iGF_Psi)**2 * X1

              END DO

              CALL ComputeGeometryX_FromScaleFactors( G(:,iX_B0(1)-iX1,iX2,iX3,:) )

              DO iNX = 1, nDOFX

                iNX1 = NodeNumberTableX(1,iNX)
                iNX2 = NodeNumberTableX(2,iNX)
                iNX3 = NodeNumberTableX(3,iNX)

                X1 = NodeCoordinate( MeshX(1), iX_B0(1)-iX1, iNX1 )
                X2 = NodeCoordinate( MeshX(2), iX2, iNX2 )
                X3 = NodeCoordinate( MeshX(3), iX3, iNX3 )

                ! --- Fluid Fields ---

                P(iNX,iX_B0(1)-iX1,iX2,iX3,iPM_D) &
                  = Interpolate1D( X1Arr, DensityArr, SIZE( X1Arr ), X1 )

                V1 = Zero
                V2 = Zero
                V3 = Interpolate1D( X1Arr, V3Arr, SIZE( X1Arr ), X1 )

                VSq = G(iNX,iX_B0(1)-iX1,iX2,iX3,iGF_Gm_dd_11) * V1**2 &
                      + G(iNX,iX_B0(1)-iX1,iX2,iX3,iGF_Gm_dd_22) * V2**2 &
                      + G(iNX,iX_B0(1)-iX1,iX2,iX3,iGF_Gm_dd_33) * V3**2

                W = One / SQRT( One - VSq )

                P(iNX,iX_B0(1)-iX1,iX2,iX3,iPM_V1) = V1
                P(iNX,iX_B0(1)-iX1,iX2,iX3,iPM_V2) = V2
                P(iNX,iX_B0(1)-iX1,iX2,iX3,iPM_V3) = V3

                P(iNX,iX_B0(1)-iX1,iX2,iX3,iPM_E) &
                  = Interpolate1D( X1Arr, PressureArr, SIZE( X1Arr ), X1 ) &
                    / ( Gamma_IDEAL - One )

                P(iNX,iX_B0(1)-iX1,iX2,iX3,iPM_Ne) = Zero

                CB1 = Zero
                CB2 = 2.0 * 1.0d13 * Gauss
                CB3 = Zero

                VdotB = G(iNX,iX_B0(1)-iX1,iX2,iX3,iGF_Gm_dd_11) * V1 * CB1 &
                        + G(iNX,iX_B0(1)-iX1,iX2,iX3,iGF_Gm_dd_22) * V2 * CB2 &
                        + G(iNX,iX_B0(1)-iX1,iX2,iX3,iGF_Gm_dd_33) * V3 * CB3

                P(iNX,iX_B0(1)-iX1,iX2,iX3,iPM_B1) = W * VdotB * V1 + CB1 / W
                P(iNX,iX_B0(1)-iX1,iX2,iX3,iPM_B2) = W * VdotB * V2 + CB2 / W
                P(iNX,iX_B0(1)-iX1,iX2,iX3,iPM_B3) = W * VdotB * V3 + CB3 / W

                P(iNX,iX_B0(1)-iX1,iX2,iX3,iPM_Chi) = Zero

              END DO

              CALL ComputePressureFromPrimitive_IDEAL &
                     ( P(:,iX_B0(1)-iX1,iX2,iX3,iPM_D ), P(:,iX_B0(1)-iX1,iX2,iX3,iPM_E ), &
                       P(:,iX_B0(1)-iX1,iX2,iX3,iPM_Ne), A(:,iX_B0(1)-iX1,iX2,iX3,iAM_P) )

              CALL ComputeConserved_MHD_Relativistic &
                     ( P(:,iX_B0(1)-iX1,iX2,iX3,iPM_D ), P(:,iX_B0(1)-iX1,iX2,iX3,iPM_V1),  &
                       P(:,iX_B0(1)-iX1,iX2,iX3,iPM_V2), P(:,iX_B0(1)-iX1,iX2,iX3,iPM_V3),  &
                       P(:,iX_B0(1)-iX1,iX2,iX3,iPM_E ), P(:,iX_B0(1)-iX1,iX2,iX3,iPM_Ne),  &
                       P(:,iX_B0(1)-iX1,iX2,iX3,iPM_B1), P(:,iX_B0(1)-iX1,iX2,iX3,iPM_B2),  &
                       P(:,iX_B0(1)-iX1,iX2,iX3,iPM_B3), P(:,iX_B0(1)-iX1,iX2,iX3,iPM_Chi), &
                       CD_I (:,iX_B0(1)-iX1,iX2,iX3),    CS1_I(:,iX_B0(1)-iX1,iX2,iX3),     &
                       CS2_I(:,iX_B0(1)-iX1,iX2,iX3),    CS3_I(:,iX_B0(1)-iX1,iX2,iX3),     &
                       CE_I (:,iX_B0(1)-iX1,iX2,iX3),    CNe_I(:,iX_B0(1)-iX1,iX2,iX3),     &
                       CB1_I(:,iX_B0(1)-iX1,iX2,iX3),    CB2_I(:,iX_B0(1)-iX1,iX2,iX3),     &
                       CB3_I(:,iX_B0(1)-iX1,iX2,iX3),    CChi_I(:,iX_B0(1)-iX1,iX2,iX3),    &
                       G(:,iX_B0(1)-iX1,iX2,iX3,iGF_Gm_dd_11), &
                       G(:,iX_B0(1)-iX1,iX2,iX3,iGF_Gm_dd_22), &
                       G(:,iX_B0(1)-iX1,iX2,iX3,iGF_Gm_dd_33), &
                       G(:,iX_B0(1)-iX1,iX2,iX3,iGF_Alpha   ), &
                       G(:,iX_B0(1)-iX1,iX2,iX3,iGF_Beta_1  ), &
                       G(:,iX_B0(1)-iX1,iX2,iX3,iGF_Beta_2  ), &
                       G(:,iX_B0(1)-iX1,iX2,iX3,iGF_Beta_3  ), &
                       A(:,iX_B0(1)-iX1,iX2,iX3,iAM_P), &
                       .FALSE. )

            END DO

            DO iX1 = 1, swX(1)

              DO iNX = 1, nDOFX

                iNX1 = NodeNumberTableX(1,iNX)
                iNX2 = NodeNumberTableX(2,iNX)
                iNX3 = NodeNumberTableX(3,iNX)

                X1 = NodeCoordinate( MeshX(1), iX_E0(1)-(iX1-1), iNX1 )
                X2 = NodeCoordinate( MeshX(2), iX2, iNX2 )
                X3 = NodeCoordinate( MeshX(3), iX3, iNX3 )

                ! --- Geometry Fields ---

                G(iNX,iX_E0(1)-(iX1-1),iX2,iX3,iGF_Alpha) &
                  = Interpolate1D( X1Arr, AlphaArr, SIZE( X1Arr ), X1 )

                G(iNX,iX_E0(1)-(iX1-1),iX2,iX3,iGF_Psi) &
                  = Interpolate1D( X1Arr, PsiArr, SIZE( X1Arr ), X1 )

                G(iNX,iX_E0(1)-(iX1-1),iX2,iX3,iGF_Beta_1) = Zero
                G(iNX,iX_E0(1)-(iX1-1),iX2,iX3,iGF_Beta_2) = Zero
                G(iNX,iX_E0(1)-(iX1-1),iX2,iX3,iGF_Beta_3) = Zero

                G(iNX,iX_E0(1)-(iX1-1),iX2,iX3,iGF_h_1) &
                  = G(iNX,iX_E0(1)-(iX1-1),iX2,iX3,iGF_Psi)**2
                G(iNX,iX_E0(1)-(iX1-1),iX2,iX3,iGF_h_2) &
                  = G(iNX,iX_E0(1)-(iX1-1),iX2,iX3,iGF_Psi)**2
                G(iNX,iX_E0(1)-(iX1-1),iX2,iX3,iGF_h_3) &
                  = G(iNX,iX_E0(1)-(iX1-1),iX2,iX3,iGF_Psi)**2 * X1

              END DO

              CALL ComputeGeometryX_FromScaleFactors( G(:,iX_E0(1)-(iX1-1),iX2,iX3,:) )

              DO iNX = 1, nDOFX

                iNX1 = NodeNumberTableX(1,iNX)
                iNX2 = NodeNumberTableX(2,iNX)
                iNX3 = NodeNumberTableX(3,iNX)

                X1 = NodeCoordinate( MeshX(1), iX_E0(1)-(iX1-1), iNX1 )
                X2 = NodeCoordinate( MeshX(2), iX2, iNX2 )
                X3 = NodeCoordinate( MeshX(3), iX3, iNX3 )

                ! --- Fluid Fields ---

                P(iNX,iX_E0(1)-(iX1-1),iX2,iX3,iPM_D) &
                  = Interpolate1D( X1Arr, DensityArr, SIZE( X1Arr ), X1 )

                V1 = Zero
                V2 = Zero
                V3 = Interpolate1D( X1Arr, V3Arr, SIZE( X1Arr ), X1 )

                VSq = G(iNX,iX_E0(1)-(iX1-1),iX2,iX3,iGF_Gm_dd_11) * V1**2 &
                      + G(iNX,iX_E0(1)-(iX1-1),iX2,iX3,iGF_Gm_dd_22) * V2**2 &
                      + G(iNX,iX_E0(1)-(iX1-1),iX2,iX3,iGF_Gm_dd_33) * V3**2

                W = One / SQRT( One - VSq )

                P(iNX,iX_E0(1)-(iX1-1),iX2,iX3,iPM_V1) = V1
                P(iNX,iX_E0(1)-(iX1-1),iX2,iX3,iPM_V2) = V2
                P(iNX,iX_E0(1)-(iX1-1),iX2,iX3,iPM_V3) = V3

                P(iNX,iX_E0(1)-(iX1-1),iX2,iX3,iPM_E) &
                  = Interpolate1D( X1Arr, PressureArr, SIZE( X1Arr ), X1 ) &
                    / ( Gamma_IDEAL - One )

                P(iNX,iX_E0(1)-(iX1-1),iX2,iX3,iPM_Ne) = Zero

                CB1 = Zero
                CB2 = 2.0 * 1.0d13 * Gauss
                CB3 = Zero

                VdotB = G(iNX,iX_E0(1)-(iX1-1),iX2,iX3,iGF_Gm_dd_11) * V1 * CB1 &
                        + G(iNX,iX_E0(1)-(iX1-1),iX2,iX3,iGF_Gm_dd_22) * V2 * CB2 &
                        + G(iNX,iX_E0(1)-(iX1-1),iX2,iX3,iGF_Gm_dd_33) * V3 * CB3

                P(iNX,iX_E0(1)-(iX1-1),iX2,iX3,iPM_B1) = W * VdotB * V1 + CB1 / W
                P(iNX,iX_E0(1)-(iX1-1),iX2,iX3,iPM_B2) = W * VdotB * V2 + CB2 / W
                P(iNX,iX_E0(1)-(iX1-1),iX2,iX3,iPM_B3) = W * VdotB * V3 + CB3 / W

                P(iNX,iX_E0(1)-(iX1-1),iX2,iX3,iPM_Chi) = Zero

              END DO

              CALL ComputePressureFromPrimitive_IDEAL &
                     ( P(:,iX_E0(1)-(iX1-1),iX2,iX3,iPM_D ), P(:,iX_E0(1)-(iX1-1),iX2,iX3,iPM_E ), &
                       P(:,iX_E0(1)-(iX1-1),iX2,iX3,iPM_Ne), A(:,iX_E0(1)-(iX1-1),iX2,iX3,iAM_P) )

              CALL ComputeConserved_MHD_Relativistic &
                     ( P(:,iX_E0(1)-(iX1-1),iX2,iX3,iPM_D ), P(:,iX_E0(1)-(iX1-1),iX2,iX3,iPM_V1),  &
                       P(:,iX_E0(1)-(iX1-1),iX2,iX3,iPM_V2), P(:,iX_E0(1)-(iX1-1),iX2,iX3,iPM_V3),  &
                       P(:,iX_E0(1)-(iX1-1),iX2,iX3,iPM_E ), P(:,iX_E0(1)-(iX1-1),iX2,iX3,iPM_Ne),  &
                       P(:,iX_E0(1)-(iX1-1),iX2,iX3,iPM_B1), P(:,iX_E0(1)-(iX1-1),iX2,iX3,iPM_B2),  &
                       P(:,iX_E0(1)-(iX1-1),iX2,iX3,iPM_B3), P(:,iX_E0(1)-(iX1-1),iX2,iX3,iPM_Chi), &
                       CD_I (:,iX_E0(1)-(iX1-1),iX2,iX3),    CS1_I(:,iX_E0(1)-(iX1-1),iX2,iX3),     &
                       CS2_I(:,iX_E0(1)-(iX1-1),iX2,iX3),    CS3_I(:,iX_E0(1)-(iX1-1),iX2,iX3),     &
                       CE_I (:,iX_E0(1)-(iX1-1),iX2,iX3),    CNe_I(:,iX_E0(1)-(iX1-1),iX2,iX3),     &
                       CB1_I(:,iX_E0(1)-(iX1-1),iX2,iX3),    CB2_I(:,iX_E0(1)-(iX1-1),iX2,iX3),     &
                       CB3_I(:,iX_E0(1)-(iX1-1),iX2,iX3),    CChi_I(:,iX_E0(1)-(iX1-1),iX2,iX3),    &
                       G(:,iX_E0(1)-(iX1-1),iX2,iX3,iGF_Gm_dd_11), &
                       G(:,iX_E0(1)-(iX1-1),iX2,iX3,iGF_Gm_dd_22), &
                       G(:,iX_E0(1)-(iX1-1),iX2,iX3,iGF_Gm_dd_33), &
                       G(:,iX_E0(1)-(iX1-1),iX2,iX3,iGF_Alpha   ), &
                       G(:,iX_E0(1)-(iX1-1),iX2,iX3,iGF_Beta_1  ), &
                       G(:,iX_E0(1)-(iX1-1),iX2,iX3,iGF_Beta_2  ), &
                       G(:,iX_E0(1)-(iX1-1),iX2,iX3,iGF_Beta_3  ), &
                       A(:,iX_E0(1)-(iX1-1),iX2,iX3,iAM_P), &
                       .FALSE. )

            END DO

            END DO
            END DO

          END IF

          IF ( ApplyOuterBC_MHD( iApplyBC ) )THEN

          ! --- Map to outer boundary. ---

            DO iX3 = iX_B0(3), iX_E0(3)
            DO iX2 = iX_B0(2), iX_E0(2)

            DO iX1 = 1, swX(1)

              DO iNX = 1, nDOFX

                iNX1 = NodeNumberTableX(1,iNX)
                iNX2 = NodeNumberTableX(2,iNX)
                iNX3 = NodeNumberTableX(3,iNX)

                X1 = NodeCoordinate( MeshX(1), iX_E0(1)+iX1, iNX1 )
                X2 = NodeCoordinate( MeshX(2), iX2, iNX2 )
                X3 = NodeCoordinate( MeshX(3), iX3, iNX3 )

                ! --- Geometry Fields ---

                G(iNX,iX_E0(1)+iX1,iX2,iX3,iGF_Alpha) &
                  = Interpolate1D( X1Arr, AlphaArr, SIZE( X1Arr ), X1 )

                G(iNX,iX_E0(1)+iX1,iX2,iX3,iGF_Psi) &
                  = Interpolate1D( X1Arr, PsiArr, SIZE( X1Arr ), X1 )

                G(iNX,iX_E0(1)+iX1,iX2,iX3,iGF_Beta_1) = Zero
                G(iNX,iX_E0(1)+iX1,iX2,iX3,iGF_Beta_2) = Zero
                G(iNX,iX_E0(1)+iX1,iX2,iX3,iGF_Beta_3) = Zero

                G(iNX,iX_E0(1)+iX1,iX2,iX3,iGF_h_1) &
                  = G(iNX,iX_E0(1)+iX1,iX2,iX3,iGF_Psi)**2
                G(iNX,iX_E0(1)+iX1,iX2,iX3,iGF_h_2) &
                  = G(iNX,iX_E0(1)+iX1,iX2,iX3,iGF_Psi)**2
                G(iNX,iX_E0(1)+iX1,iX2,iX3,iGF_h_3) &
                  = G(iNX,iX_E0(1)+iX1,iX2,iX3,iGF_Psi)**2 * X1

              END DO

              CALL ComputeGeometryX_FromScaleFactors( G(:,iX_E0(1)+iX1,iX2,iX3,:) )

              DO iNX = 1, nDOFX

                iNX1 = NodeNumberTableX(1,iNX)
                iNX2 = NodeNumberTableX(2,iNX)
                iNX3 = NodeNumberTableX(3,iNX)

                X1 = NodeCoordinate( MeshX(1), iX_E0(1)+iX1, iNX1 )
                X2 = NodeCoordinate( MeshX(2), iX2, iNX2 )
                X3 = NodeCoordinate( MeshX(3), iX3, iNX3 )

                ! --- Fluid Fields ---

                P(iNX,iX_E0(1)+iX1,iX2,iX3,iPM_D) &
                  = Interpolate1D( X1Arr, DensityArr, SIZE( X1Arr ), X1 )

                V1 = Zero
                V2 = Zero
                V3 = Interpolate1D( X1Arr, V3Arr, SIZE( X1Arr ), X1 )

                VSq = G(iNX,iX_E0(1)+iX1,iX2,iX3,iGF_Gm_dd_11) * V1**2 &
                      + G(iNX,iX_E0(1)+iX1,iX2,iX3,iGF_Gm_dd_22) * V2**2 &
                      + G(iNX,iX_E0(1)+iX1,iX2,iX3,iGF_Gm_dd_33) * V3**2

                W = One / SQRT( One - VSq )

                P(iNX,iX_E0(1)+iX1,iX2,iX3,iPM_V1) = V1
                P(iNX,iX_E0(1)+iX1,iX2,iX3,iPM_V2) = V2
                P(iNX,iX_E0(1)+iX1,iX2,iX3,iPM_V3) = V3

                P(iNX,iX_E0(1)+iX1,iX2,iX3,iPM_E) &
                  = Interpolate1D( X1Arr, PressureArr, SIZE( X1Arr ), X1 ) &
                    / ( Gamma_IDEAL - One )

                P(iNX,iX_E0(1)+iX1,iX2,iX3,iPM_Ne) = Zero

                CB1 = Zero
                CB2 = 2.0 * 1.0d13 * Gauss
                CB3 = Zero

                VdotB = G(iNX,iX_E0(1)+iX1,iX2,iX3,iGF_Gm_dd_11) * V1 * CB1 &
                        + G(iNX,iX_E0(1)+iX1,iX2,iX3,iGF_Gm_dd_22) * V2 * CB2 &
                        + G(iNX,iX_E0(1)+iX1,iX2,iX3,iGF_Gm_dd_33) * V3 * CB3

                P(iNX,iX_E0(1)+iX1,iX2,iX3,iPM_B1) = W * VdotB * V1 + CB1 / W
                P(iNX,iX_E0(1)+iX1,iX2,iX3,iPM_B2) = W * VdotB * V2 + CB2 / W
                P(iNX,iX_E0(1)+iX1,iX2,iX3,iPM_B3) = W * VdotB * V3 + CB3 / W

                P(iNX,iX_E0(1)+iX1,iX2,iX3,iPM_Chi) = Zero

              END DO

              CALL ComputePressureFromPrimitive_IDEAL &
                     ( P(:,iX_E0(1)+iX1,iX2,iX3,iPM_D ), P(:,iX_E0(1)+iX1,iX2,iX3,iPM_E ), &
                       P(:,iX_E0(1)+iX1,iX2,iX3,iPM_Ne), A(:,iX_E0(1)+iX1,iX2,iX3,iAM_P) )

              CALL ComputeConserved_MHD_Relativistic &
                     ( P(:,iX_E0(1)+iX1,iX2,iX3,iPM_D ), P(:,iX_E0(1)+iX1,iX2,iX3,iPM_V1),  &
                       P(:,iX_E0(1)+iX1,iX2,iX3,iPM_V2), P(:,iX_E0(1)+iX1,iX2,iX3,iPM_V3),  &
                       P(:,iX_E0(1)+iX1,iX2,iX3,iPM_E ), P(:,iX_E0(1)+iX1,iX2,iX3,iPM_Ne),  &
                       P(:,iX_E0(1)+iX1,iX2,iX3,iPM_B1), P(:,iX_E0(1)+iX1,iX2,iX3,iPM_B2),  &
                       P(:,iX_E0(1)+iX1,iX2,iX3,iPM_B3), P(:,iX_E0(1)+iX1,iX2,iX3,iPM_Chi), &
                       CD_I (:,iX_E0(1)+iX1,iX2,iX3),    CS1_I(:,iX_E0(1)+iX1,iX2,iX3),     &
                       CS2_I(:,iX_E0(1)+iX1,iX2,iX3),    CS3_I(:,iX_E0(1)+iX1,iX2,iX3),     &
                       CE_I (:,iX_E0(1)+iX1,iX2,iX3),    CNe_I(:,iX_E0(1)+iX1,iX2,iX3),     &
                       CB1_I(:,iX_E0(1)+iX1,iX2,iX3),    CB2_I(:,iX_E0(1)+iX1,iX2,iX3),     &
                       CB3_I(:,iX_E0(1)+iX1,iX2,iX3),    CChi_I(:,iX_E0(1)+iX1,iX2,iX3),    &
                       G(:,iX_E0(1)+iX1,iX2,iX3,iGF_Gm_dd_11), &
                       G(:,iX_E0(1)+iX1,iX2,iX3,iGF_Gm_dd_22), &
                       G(:,iX_E0(1)+iX1,iX2,iX3,iGF_Gm_dd_33), &
                       G(:,iX_E0(1)+iX1,iX2,iX3,iGF_Alpha   ), &
                       G(:,iX_E0(1)+iX1,iX2,iX3,iGF_Beta_1  ), &
                       G(:,iX_E0(1)+iX1,iX2,iX3,iGF_Beta_2  ), &
                       G(:,iX_E0(1)+iX1,iX2,iX3,iGF_Beta_3  ), &
                       A(:,iX_E0(1)+iX1,iX2,iX3,iAM_P), &
                       .FALSE. )

            END DO

            DO iX1 = 1, swX(1)

              DO iNX = 1, nDOFX

                iNX1 = NodeNumberTableX(1,iNX)
                iNX2 = NodeNumberTableX(2,iNX)
                iNX3 = NodeNumberTableX(3,iNX)

                X1 = NodeCoordinate( MeshX(1), iX_B0(1)+(iX1-1), iNX1 )
                X2 = NodeCoordinate( MeshX(2), iX2, iNX2 )
                X3 = NodeCoordinate( MeshX(3), iX3, iNX3 )

                ! --- Geometry Fields ---

                G(iNX,iX_B0(1)+(iX1-1),iX2,iX3,iGF_Alpha) &
                  = Interpolate1D( X1Arr, AlphaArr, SIZE( X1Arr ), X1 )

                G(iNX,iX_B0(1)+(iX1-1),iX2,iX3,iGF_Psi) &
                  = Interpolate1D( X1Arr, PsiArr, SIZE( X1Arr ), X1 )

                G(iNX,iX_B0(1)+(iX1-1),iX2,iX3,iGF_Beta_1) = Zero
                G(iNX,iX_B0(1)+(iX1-1),iX2,iX3,iGF_Beta_2) = Zero
                G(iNX,iX_B0(1)+(iX1-1),iX2,iX3,iGF_Beta_3) = Zero

                G(iNX,iX_B0(1)+(iX1-1),iX2,iX3,iGF_h_1) &
                  = G(iNX,iX_B0(1)+(iX1-1),iX2,iX3,iGF_Psi)**2
                G(iNX,iX_B0(1)+(iX1-1),iX2,iX3,iGF_h_2) &
                  = G(iNX,iX_B0(1)+(iX1-1),iX2,iX3,iGF_Psi)**2
                G(iNX,iX_B0(1)+(iX1-1),iX2,iX3,iGF_h_3) &
                  = G(iNX,iX_B0(1)+(iX1-1),iX2,iX3,iGF_Psi)**2 * X1

              END DO

              CALL ComputeGeometryX_FromScaleFactors( G(:,iX_B0(1)+(iX1-1),iX2,iX3,:) )

              DO iNX = 1, nDOFX

                iNX1 = NodeNumberTableX(1,iNX)
                iNX2 = NodeNumberTableX(2,iNX)
                iNX3 = NodeNumberTableX(3,iNX)

                X1 = NodeCoordinate( MeshX(1), iX_B0(1)+(iX1-1), iNX1 )
                X2 = NodeCoordinate( MeshX(2), iX2, iNX2 )
                X3 = NodeCoordinate( MeshX(3), iX3, iNX3 )

                ! --- Fluid Fields ---

                P(iNX,iX_B0(1)+(iX1-1),iX2,iX3,iPM_D) &
                  = Interpolate1D( X1Arr, DensityArr, SIZE( X1Arr ), X1 )

                V1 = Zero
                V2 = Zero
                V3 = Interpolate1D( X1Arr, V3Arr, SIZE( X1Arr ), X1 )

                VSq = G(iNX,iX_B0(1)+(iX1-1),iX2,iX3,iGF_Gm_dd_11) * V1**2 &
                      + G(iNX,iX_B0(1)+(iX1-1),iX2,iX3,iGF_Gm_dd_22) * V2**2 &
                      + G(iNX,iX_B0(1)+(iX1-1),iX2,iX3,iGF_Gm_dd_33) * V3**2

                W = One / SQRT( One - VSq )

                P(iNX,iX_B0(1)+(iX1-1),iX2,iX3,iPM_V1) = V1
                P(iNX,iX_B0(1)+(iX1-1),iX2,iX3,iPM_V2) = V2
                P(iNX,iX_B0(1)+(iX1-1),iX2,iX3,iPM_V3) = V3

                P(iNX,iX_B0(1)+(iX1-1),iX2,iX3,iPM_E) &
                  = Interpolate1D( X1Arr, PressureArr, SIZE( X1Arr ), X1 ) &
                    / ( Gamma_IDEAL - One )

                P(iNX,iX_B0(1)+(iX1-1),iX2,iX3,iPM_Ne) = Zero

                CB1 = Zero
                CB2 = 2.0 * 1.0d13 * Gauss
                CB3 = Zero

                VdotB = G(iNX,iX_B0(1)+(iX1-1),iX2,iX3,iGF_Gm_dd_11) * V1 * CB1 &
                        + G(iNX,iX_B0(1)+(iX1-1),iX2,iX3,iGF_Gm_dd_22) * V2 * CB2 &
                        + G(iNX,iX_B0(1)+(iX1-1),iX2,iX3,iGF_Gm_dd_33) * V3 * CB3

                P(iNX,iX_B0(1)+(iX1-1),iX2,iX3,iPM_B1) = W * VdotB * V1 + CB1 / W
                P(iNX,iX_B0(1)+(iX1-1),iX2,iX3,iPM_B2) = W * VdotB * V2 + CB2 / W
                P(iNX,iX_B0(1)+(iX1-1),iX2,iX3,iPM_B3) = W * VdotB * V3 + CB3 / W

                P(iNX,iX_B0(1)+(iX1-1),iX2,iX3,iPM_Chi) = Zero

              END DO

              CALL ComputePressureFromPrimitive_IDEAL &
                     ( P(:,iX_B0(1)+(iX1-1),iX2,iX3,iPM_D ), P(:,iX_B0(1)+(iX1-1),iX2,iX3,iPM_E ), &
                       P(:,iX_B0(1)+(iX1-1),iX2,iX3,iPM_Ne), A(:,iX_B0(1)+(iX1-1),iX2,iX3,iAM_P) )

              CALL ComputeConserved_MHD_Relativistic &
                     ( P(:,iX_B0(1)+(iX1-1),iX2,iX3,iPM_D ), P(:,iX_B0(1)+(iX1-1),iX2,iX3,iPM_V1),  &
                       P(:,iX_B0(1)+(iX1-1),iX2,iX3,iPM_V2), P(:,iX_B0(1)+(iX1-1),iX2,iX3,iPM_V3),  &
                       P(:,iX_B0(1)+(iX1-1),iX2,iX3,iPM_E ), P(:,iX_B0(1)+(iX1-1),iX2,iX3,iPM_Ne),  &
                       P(:,iX_B0(1)+(iX1-1),iX2,iX3,iPM_B1), P(:,iX_B0(1)+(iX1-1),iX2,iX3,iPM_B2),  &
                       P(:,iX_B0(1)+(iX1-1),iX2,iX3,iPM_B3), P(:,iX_B0(1)+(iX1-1),iX2,iX3,iPM_Chi), &
                       CD_I (:,iX_B0(1)+(iX1-1),iX2,iX3),    CS1_I(:,iX_B0(1)+(iX1-1),iX2,iX3),     &
                       CS2_I(:,iX_B0(1)+(iX1-1),iX2,iX3),    CS3_I(:,iX_B0(1)+(iX1-1),iX2,iX3),     &
                       CE_I (:,iX_B0(1)+(iX1-1),iX2,iX3),    CNe_I(:,iX_B0(1)+(iX1-1),iX2,iX3),     &
                       CB1_I(:,iX_B0(1)+(iX1-1),iX2,iX3),    CB2_I(:,iX_B0(1)+(iX1-1),iX2,iX3),     &
                       CB3_I(:,iX_B0(1)+(iX1-1),iX2,iX3),    CChi_I(:,iX_B0(1)+(iX1-1),iX2,iX3),    &
                       G(:,iX_B0(1)+(iX1-1),iX2,iX3,iGF_Gm_dd_11), &
                       G(:,iX_B0(1)+(iX1-1),iX2,iX3,iGF_Gm_dd_22), &
                       G(:,iX_B0(1)+(iX1-1),iX2,iX3,iGF_Gm_dd_33), &
                       G(:,iX_B0(1)+(iX1-1),iX2,iX3,iGF_Alpha   ), &
                       G(:,iX_B0(1)+(iX1-1),iX2,iX3,iGF_Beta_1  ), &
                       G(:,iX_B0(1)+(iX1-1),iX2,iX3,iGF_Beta_2  ), &
                       G(:,iX_B0(1)+(iX1-1),iX2,iX3,iGF_Beta_3  ), &
                       A(:,iX_B0(1)+(iX1-1),iX2,iX3,iAM_P), &
                       .FALSE. )

            END DO

            END DO
            END DO

          END IF

          DEALLOCATE( X1Arr, PsiArr, AlphaArr, DensityArr, V3Arr, PressureArr )

        END IF

        ! --- Inner Boundary --

        IF( ApplyInnerBC_MHD( iApplyBC ) )THEN

          DO iX3 = iX_B0(3), iX_E0(3)
          DO iX2 = iX_B0(2), iX_E0(2)
          DO iX1 = 1, swX(1)
          DO iNX = 1, nDOFX

            U(iNX,iX_B0(1)-iX1,iX2,iX3,1) &
              = U(iNX,iX_E0(1)-(iX1-1),iX2,iX3,1) &
                - CD_I(iNX,iX_E0(1)-(iX1-1),iX2,iX3) &
                + CD_I(iNX,iX_B0(1)-iX1,iX2,iX3)

            U(iNX,iX_B0(1)-iX1,iX2,iX3,2) &
              = U(iNX,iX_E0(1)-(iX1-1),iX2,iX3,2) &
                - CS1_I(iNX,iX_E0(1)-(iX1-1),iX2,iX3) &
                + CS1_I(iNX,iX_B0(1)-iX1,iX2,iX3)

            U(iNX,iX_B0(1)-iX1,iX2,iX3,3) &
              = U(iNX,iX_E0(1)-(iX1-1),iX2,iX3,3) &
                - CS2_I(iNX,iX_E0(1)-(iX1-1),iX2,iX3) &
                + CS2_I(iNX,iX_B0(1)-iX1,iX2,iX3)

            U(iNX,iX_B0(1)-iX1,iX2,iX3,4) &
              = U(iNX,iX_E0(1)-(iX1-1),iX2,iX3,4) &
                - CS3_I(iNX,iX_E0(1)-(iX1-1),iX2,iX3) &
                + CS3_I(iNX,iX_B0(1)-iX1,iX2,iX3)

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

            U(iNX,iX_E0(1)+iX1,iX2,iX3,1) &
              = U(iNX,iX_B0(1)+(iX1-1),iX2,iX3,1) &
                - CD_I(iNX,iX_B0(1)+(iX1-1),iX2,iX3) &
                + CD_I(iNX,iX_E0(1)+iX1,iX2,iX3)

            U(iNX,iX_E0(1)+iX1,iX2,iX3,2) &
              = U(iNX,iX_B0(1)+(iX1-1),iX2,iX3,2) &
                - CS1_I(iNX,iX_B0(1)+(iX1-1),iX2,iX3) &
                + CS1_I(iNX,iX_E0(1)+iX1,iX2,iX3)

            U(iNX,iX_E0(1)+iX1,iX2,iX3,3) &
              = U(iNX,iX_B0(1)+(iX1-1),iX2,iX3,3) &
                - CS2_I(iNX,iX_B0(1)+(iX1-1),iX2,iX3) &
                + CS2_I(iNX,iX_E0(1)+iX1,iX2,iX3)

            U(iNX,iX_E0(1)+iX1,iX2,iX3,4) &
              = U(iNX,iX_B0(1)+(iX1-1),iX2,iX3,4) &
                - CS3_I(iNX,iX_B0(1)+(iX1-1),iX2,iX3) &
                + CS3_I(iNX,iX_E0(1)+iX1,iX2,iX3)

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


  SUBROUTINE ApplyBC_MHD_X2( t, iX_B0, iX_E0, iX_B1, iX_E1, U, iApplyBC )

    REAL(DP), INTENT(in)    :: t
    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), &
      iApplyBC
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

    CASE DEFAULT

    END SELECT

  END SUBROUTINE ApplyBC_MHD_X2


  SUBROUTINE ApplyBC_MHD_X3( t, iX_B0, iX_E0, iX_B1, iX_E1, U, iApplyBC )

    REAL(DP), INTENT(in)    :: t
    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), &
      iApplyBC
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


  SUBROUTINE ReadDataset3DHDF( Dataset, DatasetName, FILE_ID )

    REAL(DP),         INTENT(out) :: Dataset(:,:,:)
    CHARACTER(LEN=*), INTENT(in)  :: DatasetName
    INTEGER(HID_T),   INTENT(in)  :: FILE_ID

    INTEGER(HID_T) :: DATASET_ID
    INTEGER(HID_T) :: DATASIZE(3)

    DATASIZE = SHAPE( Dataset )

    CALL H5DOPEN_F( FILE_ID, TRIM( DatasetName ), DATASET_ID, HDFERR )

    CALL H5DREAD_F( DATASET_ID, H5T_NATIVE_DOUBLE, Dataset, DATASIZE, HDFERR )

    CALL H5DCLOSE_F( DATASET_ID, HDFERR )

  END SUBROUTINE ReadDataset3DHDF


  SUBROUTINE ReadDataset1DHDF( Dataset, DatasetName, FILE_ID )

    REAL(DP),         INTENT(out) :: Dataset(:)
    CHARACTER(LEN=*), INTENT(in)  :: DatasetName
    INTEGER(HID_T),   INTENT(in)  :: FILE_ID

    INTEGER(HID_T) :: DATASET_ID
    INTEGER(HID_T) :: DATASIZE(1)

    DATASIZE = SHAPE( Dataset )

    CALL H5DOPEN_F( FILE_ID, TRIM( DatasetName ), DATASET_ID, HDFERR )

    CALL H5DREAD_F( DATASET_ID, H5T_NATIVE_DOUBLE, Dataset, DATASIZE, HDFERR )

    CALL H5DCLOSE_F( DATASET_ID, HDFERR )

  END SUBROUTINE ReadDataset1DHDF


  REAL(DP) FUNCTION Interpolate1D( x, y, n, xq )

    INTEGER,                INTENT(in) :: n
    REAL(DP), DIMENSION(n), INTENT(in) :: x, y
    REAL(DP),               INTENT(in) :: xq

    INTEGER :: i

    i = Locate( xq, x, n )

    !PRINT*, 'i: ', i

    IF( i == 0 )THEN

      ! --- Extrapolate Left ---

      Interpolate1D &
        = Interpolate1D_Linear( xq, x(1), x(2), y(1), y(2) )

      !PRINT*, 'x(1): ', x(1)
      !PRINT*, 'x(2): ', x(2)
      !PRINT*, 'y(1): ', y(1)
      !PRINT*, 'y(2): ', y(2)

    ELSE IF( i == n )THEN

      ! --- Extrapolate Right ---

      Interpolate1D &
        = Interpolate1D_Linear( xq, x(n-1), x(n), y(n-1), y(n) )

      !PRINT*, 'x(n-1): ', x(n-1)
      !PRINT*, 'x(n): ',   x(n)
      !PRINT*, 'y(n-1): ', y(n-1)
      !PRINT*, 'y(n): ',   y(n)


    ELSE

      Interpolate1D &
        = Interpolate1D_Linear( xq, x(i), x(i+1), y(i), y(i+1) )

      !PRINT*, 'x(i): ', x(i)
      !PRINT*, 'x(i+1): ', x(i+1)
      !PRINT*, 'y(i): ', y(i)
      !PRINT*, 'y(i+1): ', y(i+1)

    END IF

    RETURN

  END FUNCTION Interpolate1D


END MODULE MHD_BoundaryConditionsModule
