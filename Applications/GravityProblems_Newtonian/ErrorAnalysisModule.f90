MODULE ErrorAnalysisModule

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    nDOFX, nX, nNodes
  USE FluidFieldsModule, ONLY: &
    WeightsF, &
    uPF, nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E

  IMPLICIT NONE
  PRIVATE

  REAL(DP), DIMENSION(:),         ALLOCATABLE :: One_Error
  REAL(DP), DIMENSION(:),         ALLOCATABLE :: Inf_Error
  REAL(DP), DIMENSION(:,:,:,:,:), ALLOCATABLE :: uPF_0

  PUBLIC :: InitializeErrorAnalysis
  PUBLIC :: FinalizeErrorAnalysis

CONTAINS


  SUBROUTINE InitializeErrorAnalysis

    ALLOCATE( One_Error(nPF) )
    ALLOCATE( Inf_Error(nPF) )
    ALLOCATE( uPF_0(nDOFX,nX(1),nX(2),nX(3),nPF) )

    uPF_0 = uPF(1:nDOFX,1:nX(1),1:nX(2),1:nX(3),1:nPF)

  END SUBROUTINE InitializeErrorAnalysis


  SUBROUTINE FinalizeErrorAnalysis

    INTEGER :: iX1, iX2, iX3, iPF
    REAL(DP), DIMENSION(nDOFX) :: Error

    One_Error = 0.0_DP
    Inf_Error = 0.0_DP

    DO iPF = 1, nPF
      DO iX3 = 1, nX(3)
        DO iX2 = 1, nX(2)
          DO iX1 = 1, nX(1)

            Error(:) &
              = ABS( uPF(:,iX1,iX2,iX3,iPF) &
                     - uPF_0(:,iX1,iX2,iX3,iPF) )

            One_Error(iPF) &
              = One_Error(iPF) &
                  + DOT_PRODUCT( WeightsF, Error )

            Inf_Error(iPF) &
              = MAX( Inf_Error(iPF), MAXVAL( Error ) )

          END DO
        END DO
      END DO
    END DO

    One_Error &
      = One_Error / REAL( nDOFX * PRODUCT( nX ) )

    WRITE(*,*)
    WRITE(*,'(A4,A)') '', '--- Error Analysis ---'
    WRITE(*,*)
    WRITE(*,'(A4,A10,3I5.4)') '', 'nX = ', nX
    WRITE(*,'(A4,A10,1I5.4)') '', 'nNodes = ', nNodes
    WRITE(*,*)
    WRITE(*,'(A4,A12,A4,A12)') '', 'L_1', '', 'L_Inf'
    WRITE(*,*)
    WRITE(*,'(A4,A4,ES12.6E2,A4,ES12.6E2)') &
      '', 'D: ',  One_Error(iPF_D ), '', Inf_Error(iPF_D )
    WRITE(*,'(A4,A4,ES12.6E2,A4,ES12.6E2)') &
      '', 'V1: ', One_Error(iPF_V1), '', Inf_Error(iPF_V1)
    WRITE(*,'(A4,A4,ES12.6E2,A4,ES12.6E2)') &
      '', 'V2: ', One_Error(iPF_V2), '', Inf_Error(iPF_V2)
    WRITE(*,'(A4,A4,ES12.6E2,A4,ES12.6E2)') &
      '', 'V3: ', One_Error(iPF_V3), '', Inf_Error(iPF_V3)
    WRITE(*,'(A4,A4,ES12.6E2,A4,ES12.6E2)') &
      '', 'E: ',  One_Error(iPF_E ), '', Inf_Error(iPF_E )
    WRITE(*,*)

    DEALLOCATE( One_Error, Inf_Error, uPF_0 )

  END SUBROUTINE FinalizeErrorAnalysis

END MODULE ErrorAnalysisModule
