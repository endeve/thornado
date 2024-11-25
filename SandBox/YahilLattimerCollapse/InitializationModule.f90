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
  USE Euler_UtilitiesModule_NonRelativistic, ONLY: &
    ComputeConserved_Euler_NonRelativistic

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

    INTEGER  :: iX1, iX2, iX3, N
    INTEGER  :: iNodeX, iNodeX1, i_Lo, i_Hi
    REAL(DP) :: Kappa, C_X, C_D, C_V, R_q, X_q
    REAL(DP) :: w_Lo, w_Hi
    REAL(DP), ALLOCATABLE :: X(:), D(:), V(:)

    Kappa = CentralPressure / CentralDensity**( Gamma )

    C_X = SQRT( Kappa ) &
            * GravitationalConstant**( 0.5_DP * ( 1.0_DP - Gamma ) ) &
            * CollapseTime**( 2.0_DP - Gamma )

    C_D = GravitationalConstant**( - 1.0_DP ) &
            * CollapseTime**( - 2.0_DP )

    C_V = SQRT( Kappa ) &
            * GravitationalConstant**( 0.5_DP * ( 1.0_DP - Gamma ) ) &
            * CollapseTime**( 1.0_DP - Gamma )

    CALL length( FileName, N )

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

            uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
              = ( V(i_Lo) * w_Lo + V(i_Hi) * w_Hi ) * C_V

            uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP

            uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP

            uPF(iNodeX,iX1,iX2,iX3,iPF_E) &
              = Kappa * uPF(iNodeX,iX1,iX2,iX3,iPF_D) ** Gamma &
                  / ( Gamma - 1.0_DP )

          END DO

          CALL ComputeConserved_Euler_NonRelativistic &
               ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_V1), &
                 uPF(:,iX1,iX2,iX3,iPF_V2), uPF(:,iX1,iX2,iX3,iPF_V3), &
                 uPF(:,iX1,iX2,iX3,iPF_E ), uPF(:,iX1,iX2,iX3,iPF_Ne), &
                 uCF(:,iX1,iX2,iX3,iCF_D ), uCF(:,iX1,iX2,iX3,iCF_S1), &
                 uCF(:,iX1,iX2,iX3,iCF_S2), uCF(:,iX1,iX2,iX3,iCF_S3), &
                 uCF(:,iX1,iX2,iX3,iCF_E ), uCF(:,iX1,iX2,iX3,iCF_Ne), &
                 uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33) )

        END DO
      END DO
    END DO

    DEALLOCATE( X, D, V )

  END SUBROUTINE InitializeFields
  
  subroutine length( file, N )

    !===========================================
    !                    Length
    !===========================================
    !  Purpose: Given a file name, it will return the 
    !           number of line in the file
    !  Author: Kristopher Andrew
    !  Date: 7/3/18                              
    !===========================================

    character(len=*), intent(in)  :: file
    integer,          intent(out) :: N

    integer :: stat

    N = 0

    OPEN(UNIT=1,FILE=file,FORM="FORMATTED",STATUS="OLD",ACTION="READ")

    READ(1, *) ! --- Ignore First Line
    DO
       READ(1,*,IOSTAT = stat)
       IF(stat /= 0)THEN
          CLOSE(UNIT=1)
          RETURN
       END IF
       N = N+1
    END DO
  
  END SUBROUTINE length


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
