MODULE InitializationModule

  USE KindModule, ONLY: &
    DP
  USE UnitsModule, ONLY: &
    GravitationalConstant
  USE ProgramHeaderModule, ONLY: &
    iX_B0, iX_B1, iX_E0, iX_E1, &
    nDOFX
  USE UtilitiesModule, ONLY: &
    Locate
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    uGF, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33
  USE FluidFieldsModule, ONLY: &
    uPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne, &
    uCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne
  USE EulerEquationsUtilitiesModule_Beta, ONLY: &
    ComputeConserved

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFields

CONTAINS


  SUBROUTINE InitializeFields &
    ( FileName, Gamma, CollapseTime, CentralDensity, CentralPressure )

    CHARACTER(*), INTENT(in) :: FileName
    REAL(DP),     INTENT(in) :: Gamma
    REAL(DP),     INTENT(in) :: CollapseTime
    REAL(DP),     INTENT(in) :: CentralDensity
    REAL(DP),     INTENT(in) :: CentralPressure

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX, iNodeX1, i_Lo, i_Hi
    INTEGER, PARAMETER :: N = 203
    REAL(DP) :: Kappa, C_X, C_D, C_V, R_q, X_q
    REAL(DP) :: w_Lo, w_Hi
    REAL(DP), ALLOCATABLE :: X(:), D(:), V(:)

    Kappa = CentralPressure / CentralDensity**Gamma

    C_X = SQRT( Kappa ) &
            * GravitationalConstant**( 0.5_DP * ( 1.0_DP - Gamma ) ) &
            * CollapseTime**( 2.0_DP - Gamma )

    C_D = GravitationalConstant**( - 1.0_DP ) &
            * CollapseTime**( - 2.0_DP )

    C_V = SQRT( Kappa ) &
            * GravitationalConstant**( 0.5_DP * ( 1.0_DP - Gamma ) ) &
            * CollapseTime**( 1.0_DP - Gamma )

    ! --- Get Array Size Here ---

    ALLOCATE( X(N), D(N), V(N) )

    CALL ReadYahilProfile( FileName, N, X, D, V )

    DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E1(1)

          DO iNodeX = 1, nDOFX

            iNodeX1 = NodeNumberTableX(1,iNodeX)

            R_q = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

            X_q = R_q / C_X

            i_Lo = Locate( X_q, X, N )
            i_Hi = i_Lo + 1

            w_Lo = ( X(i_Hi) - X_q ) / ( X(i_Hi) - X(i_Lo) )
            w_Hi = ( X_q - X(i_Lo) ) / ( X(i_Hi) - X(i_Lo) )

            uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
              = ( D(i_Lo) * w_Lo + D(i_Hi) * w_Hi ) * C_D

          END DO

        END DO
      END DO
    END DO

    DEALLOCATE( X, D, V )

  END SUBROUTINE InitializeFields


  SUBROUTINE ReadYahilProfile( FILE_NAME, N, X, D, V )

    !===========================================
    !             ReadYahilProfile
    !===========================================
    !  Purpose: This subroutine reads a file and returns the 
    !           position, density, and velocity
    !  Author: Kristopher Andrew
    !  Date: 6/29/18
    !===========================================

    CHARACTER(LEN=*), INTENT(in)    :: FILE_NAME
    INTEGER,          INTENT(in)    :: N
    REAL(DP),         INTENT(inout) :: X(N), D(N), V(N)

    INTEGER :: i

    !--- Open the file.
    OPEN(UNIT = 1, FILE = FILE_NAME, STATUS="OLD", ACTION="READ")

    !--- Read all the files. Discard the header
    READ(1, *)
    DO i = 1, N
      READ(1, *) X(i), D(i), V(i)
    END DO

    !--- Close the file
    CLOSE(UNIT = 1)

  END SUBROUTINE ReadYahilProfile


END MODULE InitializationModule
