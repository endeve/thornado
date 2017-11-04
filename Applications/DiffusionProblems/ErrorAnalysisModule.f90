MODULE ErrorAnalysisModule

  USE KindModule, ONLY: &
    DP, Zero, Third, Half, One, Three, &
    Pi, TwoPi, FourPi
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    nDOF, nE, nX, nNodes
  USE ReferenceElementModule, ONLY: &
    NodeNumberTable1D3D
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE RadiationFieldsModule, ONLY: &
    WeightsR, nSpecies, &
    uCR, nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3

  IMPLICIT NONE
  PRIVATE

  REAL(DP), DIMENSION(:),             ALLOCATABLE :: One_Error
  REAL(DP), DIMENSION(:),             ALLOCATABLE :: Inf_Error
  REAL(DP), DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE :: uCR_Exact

  PUBLIC :: InitializeErrorAnalysis
  PUBLIC :: FinalizeErrorAnalysis

CONTAINS


  SUBROUTINE InitializeErrorAnalysis( Time )

    REAL(DP), INTENT(in) :: Time

    ALLOCATE( One_Error(nCR) )
    ALLOCATE( Inf_Error(nCR) )
    ALLOCATE( uCR_Exact(nDOF,nE,nX(1),nX(2),nX(3),nCR,nSpecies) )

    SELECT CASE ( TRIM( ProgramName ) )

      CASE( 'SineWaveDiffusion1D' )

        CALL ComputeExactSolution_SineWaveDiffusion1D( Time )

      CASE( 'SineWaveDamping1D' )

        CALL ComputeExactSolution_SineWaveDamping1D( Time )

    END SELECT

  END SUBROUTINE InitializeErrorAnalysis


  SUBROUTINE FinalizeErrorAnalysis

    INTEGER :: iE, iX1, iX2, iX3, iCR, iS
    REAL(DP), DIMENSION(nDOF) :: Error

    One_Error = 0.0_DP
    Inf_Error = 0.0_DP

    DO iS = 1, nSpecies
      DO iCR = 1, nCR
        DO iX3 = 1, nX(3)
          DO iX2 = 1, nX(2)
            DO iX1 = 1, nX(1)
              DO iE = 1, nE

                Error(:) &
                  = ABS( uCR(:,iE,iX1,iX2,iX3,iCR,iS) &
                         - uCR_Exact(:,iE,iX1,iX2,iX3,iCR,iS) )

                One_Error(iCR) &
                  = One_Error(iCR) &
                      + DOT_PRODUCT( WeightsR, Error )

                Inf_Error(iCR) &
                  = MAX( Inf_Error(iCR), MAXVAL( Error ) )

              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

    One_Error &
      = One_Error / REAL( nDOF * nE * PRODUCT( nX ) * nSpecies )

    WRITE(*,*)
    WRITE(*,'(A4,A)') '', '--- Error Analysis ---'
    WRITE(*,*)
    WRITE(*,'(A4,A10,3I5.4)') '', 'nX = ', nX
    WRITE(*,'(A4,A10,1I5.4)') '', 'nE = ', nE
    WRITE(*,'(A4,A10,1I5.4)') '', 'nNodes = ', nNodes
    WRITE(*,*)
    WRITE(*,'(A4,A12,A4,A12)') '', 'L_1', '', 'L_Inf'
    WRITE(*,*)
    WRITE(*,'(A4,A4,ES12.6E2,A4,ES12.6E2)') &
      '', 'N: ',  One_Error(iCR_N ), '', Inf_Error(iCR_N )
    WRITE(*,'(A4,A4,ES12.6E2,A4,ES12.6E2)') &
      '', 'G1: ', One_Error(iCR_G1), '', Inf_Error(iCR_G1)
    WRITE(*,'(A4,A4,ES12.6E2,A4,ES12.6E2)') &
      '', 'G2: ', One_Error(iCR_G2), '', Inf_Error(iCR_G2)
    WRITE(*,'(A4,A4,ES12.6E2,A4,ES12.6E2)') &
      '', 'G3: ', One_Error(iCR_G3), '', Inf_Error(iCR_G3)
    WRITE(*,*)

    DEALLOCATE( One_Error, Inf_Error, uCR_Exact )

  END SUBROUTINE FinalizeErrorAnalysis


  SUBROUTINE ComputeExactSolution_SineWaveDiffusion1D( Time )

    REAL(DP), INTENT(in) :: Time

    INTEGER  :: iS, iX1, iX2, iX3, iE
    INTEGER  :: iNodeX1, iNode
    REAL(DP) :: Kappa, X1

    Kappa = 1.0d+2

    DO iS = 1, nSpecies

      DO iX3 = 1, nX(3) 
        DO iX2 = 1, nX(2)
          DO iX1 = 1, nX(1)
            DO iE = 1, nE

              DO iNode = 1, nDOF

                iNodeX1 = NodeNumberTable1D3D(2,iNode)
                X1      = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

                uCR_Exact(iNode,iE,iX1,iX2,iX3,iCR_N,iS) &
                  = Three * SQRT( FourPi ) &
                    * ( One + EXP( - Pi**2 * Time / ( 27.0_DP * Kappa ) ) &
                              * SIN( Third * Pi * X1 ) )

                uCR_Exact(iNode,iE,iX1,iX2,iX3,iCR_G1,iS) &
                  = - ( Pi * SQRT( FourPi ) / ( Three * Kappa ) ) &
                    * EXP( - Pi**2 * Time / ( 27.0_DP * Kappa ) ) &
                    * COS( Third * Pi * X1 )

                uCR_Exact(iNode,iE,iX1,iX2,iX3,iCR_G2,iS) &
                  = Zero

                uCR_Exact(iNode,iE,iX1,iX2,iX3,iCR_G3,iS) &
                  = Zero

              END DO

            END DO
          END DO
        END DO
      END DO

    END DO

  END SUBROUTINE ComputeExactSolution_SineWaveDiffusion1D


  SUBROUTINE ComputeExactSolution_SineWaveDamping1D( Time )

    REAL(DP), INTENT(in) :: Time

    INTEGER  :: iS, iX1, iX2, iX3, iE
    INTEGER  :: iNodeX1, iNode
    REAL(DP) :: Kappa, X1

    Kappa = 1.0d+1

    DO iS = 1, nSpecies

      DO iX3 = 1, nX(3) 
        DO iX2 = 1, nX(2)
          DO iX1 = 1, nX(1)
            DO iE = 1, nE

              DO iNode = 1, nDOF

                iNodeX1 = NodeNumberTable1D3D(2,iNode)
                X1      = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

                uCR_Exact(iNode,iE,iX1,iX2,iX3,iCR_N,iS) &
                  = ( One + Half * SIN( TwoPi * ( X1 - Time ) ) ) &
                      * EXP( - Kappa * Time )

                uCR_Exact(iNode,iE,iX1,iX2,iX3,iCR_G1,iS) &
                  = uCR_Exact(iNode,iE,iX1,iX2,iX3,iCR_N,iS)

                uCR_Exact(iNode,iE,iX1,iX2,iX3,iCR_G2,iS) &
                  = Zero

                uCR_Exact(iNode,iE,iX1,iX2,iX3,iCR_G3,iS) &
                  = Zero

              END DO

            END DO
          END DO
        END DO
      END DO

    END DO

  END SUBROUTINE ComputeExactSolution_SineWaveDamping1D


END MODULE ErrorAnalysisModule
