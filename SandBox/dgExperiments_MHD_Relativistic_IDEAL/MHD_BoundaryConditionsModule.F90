MODULE MHD_BoundaryConditionsModule

  USE KindModule, ONLY: &
    DP, &
    Zero
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    bcX, &
    swX, &
    nDOFX, &
    nNodesX, &
    nX
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE InputOutputUtilitiesModule, ONLY: &
    FromField3D
  USE UtilitiesModule, ONLY: &
    NodeNumberX
  USE MagnetofluidFieldsModule, ONLY: &
    nCM, &
    namesCM, &
    unitsCM, &
    iCM_D, &
    iCM_S1, &
    iCM_S2, &
    iCM_S3, &
    iCM_E, &
    iCM_Ne, &
    iCM_B1, &
    iCM_B2, &
    iCM_B3, &
    iCM_Chi

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

    iApplyBC = iApplyBC_MHD_Both

    IF( PRESENT( iApplyBC_Option ) ) &
      iApplyBC = iApplyBC_Option

    CALL ApplyBC_MHD_X1( t, iX_B0, iX_E0, iX_B1, iX_E1, U, iApplyBC(1) )

    CALL ApplyBC_MHD_X2( t, iX_B0, iX_E0, iX_B1, iX_E1, U, iApplyBC(2) )

    CALL ApplyBC_MHD_X3( t, iX_B0, iX_E0, iX_B1, iX_E1, U, iApplyBC(3) )

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

    REAL(DP) :: CD_I(1:nDOFX,1:nX(1), &
                             1:nX(2), &
                             1:nX(3))

    REAL(DP) :: S1_I(1:nDOFX,1:nX(1), &
                             1:nX(2), &
                             1:nX(3))

    REAL(DP) :: S2_I(1:nDOFX,1:nX(1), &
                             1:nX(2), &
                             1:nX(3))

    REAL(DP) :: S3_I(1:nDOFX,1:nX(1), &
                             1:nX(2), &
                             1:nX(3))

    REAL(DP) :: Dataset3D(nX(1)*nNodesX(1), &
                          nX(2)*nNodesX(2), &
                          nX(3)*nNodesX(3))

    INTEGER  :: iCM, iX1, iX2, iX3
    INTEGER  :: iNX, iNX_0
    INTEGER  :: iNX1, iNX2, iNX3, jNX, jNX1
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

          WRITE( FileNumberString, FMT='(i6.6)') FileNumber

          FileName &
            = OutputDirectory // '/' // &
              TRIM( ProgramName ) // '_' // &
              MagnetofluidSuffix // '_' // &
              FileNumberString // '.h5'

          CALL H5OPEN_F( HDFERR )

          CALL H5FOPEN_F( TRIM( FileName ), H5F_ACC_RDONLY_F, FILE_ID, HDFERR )

          GroupName = 'Magnetofluid Fields/Conserved'

          DatasetName = TRIM( GroupName ) // '/' // TRIM( namesCM(1) )

          CALL ReadDataset3DHDF( Dataset3D, DatasetName, FILE_ID )

          CD_I(1:nDOFX,1:nX(1),1:nX(2),1:nX(3)) &
            = FromField3D( Dataset3D, nX, nNodesX, nDOFX, NodeNumberTableX ) &
                * unitsCM(1)

          DatasetName = TRIM( GroupName ) // '/' // TRIM( namesCM(2) )

          CALL ReadDataset3DHDF( Dataset3D, DatasetName, FILE_ID )

          S1_I(1:nDOFX,1:nX(1),1:nX(2),1:nX(3)) &
            = FromField3D( Dataset3D, nX, nNodesX, nDOFX, NodeNumberTableX ) &
                * unitsCM(2)

          DatasetName = TRIM( GroupName ) // '/' // TRIM( namesCM(3) )

          CALL ReadDataset3DHDF( Dataset3D, DatasetName, FILE_ID )

          S2_I(1:nDOFX,1:nX(1),1:nX(2),1:nX(3)) &
            = FromField3D( Dataset3D, nX, nNodesX, nDOFX, NodeNumberTableX ) &
                * unitsCM(3)

          DatasetName = TRIM( GroupName ) // '/' // TRIM( namesCM(4) )

          CALL ReadDataset3DHDF( Dataset3D, DatasetName, FILE_ID )

          S3_I(1:nDOFX,1:nX(1),1:nX(2),1:nX(3)) &
            = FromField3D( Dataset3D, nX, nNodesX, nDOFX, NodeNumberTableX ) &
                * unitsCM(4)

        END IF

        ! --- Inner Boundary --

        IF( ApplyOuterBC_MHD( iApplyBC ) )THEN

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
                - S1_I(iNX,iX_E0(1)-(iX1-1),iX2,iX3) &
                + S1_I(iNX,iX_B0(1)-iX1,iX2,iX3)

            U(iNX,iX_B0(1)-iX1,iX2,iX3,3) &
              = U(iNX,iX_E0(1)-(iX1-1),iX2,iX3,3) &
                - S2_I(iNX,iX_E0(1)-(iX1-1),iX2,iX3) &
                + S2_I(iNX,iX_B0(1)-iX1,iX2,iX3)

            U(iNX,iX_B0(1)-iX1,iX2,iX3,4) &
              = U(iNX,iX_E0(1)-(iX1-1),iX2,iX3,4) &
                - S3_I(iNX,iX_E0(1)-(iX1-1),iX2,iX3) &
                + S3_I(iNX,iX_B0(1)-iX1,iX2,iX3)

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
                - S1_I(iNX,iX_B0(1)+(iX1-1),iX2,iX3) &
                + S1_I(iNX,iX_E0(1)+iX1,iX2,iX3)

            U(iNX,iX_E0(1)+iX1,iX2,iX3,3) &
              = U(iNX,iX_B0(1)+(iX1-1),iX2,iX3,3) &
                - S2_I(iNX,iX_B0(1)+(iX1-1),iX2,iX3) &
                + S2_I(iNX,iX_E0(1)+iX1,iX2,iX3)

            U(iNX,iX_E0(1)+iX1,iX2,iX3,4) &
              = U(iNX,iX_B0(1)+(iX1-1),iX2,iX3,4) &
                - S3_I(iNX,iX_B0(1)+(iX1-1),iX2,iX3) &
                + S3_I(iNX,iX_E0(1)+iX1,iX2,iX3)

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


END MODULE MHD_BoundaryConditionsModule
