MODULE ReferenceElementModuleX

  USE KindModule, ONLY: &
    DP, &
    One, &
    Half
  USE QuadratureModule, ONLY: &
    GetQuadrature
  USE ProgramHeaderModule, ONLY: &
    nNodesX, nDOFX, nDimsX

  IMPLICIT NONE
  PRIVATE

  INTEGER,               PUBLIC :: nDOFX_X1
  INTEGER,               PUBLIC :: nDOFX_X2
  INTEGER,               PUBLIC :: nDOFX_X3
  INTEGER,  ALLOCATABLE, PUBLIC :: NodeNumberTableX(:,:)
  INTEGER,  ALLOCATABLE, PUBLIC :: NodeNumberTableX_X1(:,:)
  INTEGER,  ALLOCATABLE, PUBLIC :: NodeNumberTableX_X2(:,:)
  INTEGER,  ALLOCATABLE, PUBLIC :: NodeNumberTableX_X3(:,:)
  INTEGER,  ALLOCATABLE, PUBLIC :: NodeNumberTableX3D(:,:,:)
  REAL(DP), ALLOCATABLE, PUBLIC :: NodesX1(:), WeightsX1(:)
  REAL(DP), ALLOCATABLE, PUBLIC :: NodesX2(:), WeightsX2(:)
  REAL(DP), ALLOCATABLE, PUBLIC :: NodesX3(:), WeightsX3(:)
  REAL(DP), ALLOCATABLE, PUBLIC :: WeightsX_X1(:)
  REAL(DP), ALLOCATABLE, PUBLIC :: WeightsX_X2(:)
  REAL(DP), ALLOCATABLE, PUBLIC :: WeightsX_X3(:)
  REAL(DP), ALLOCATABLE, PUBLIC :: WeightsX_q(:)
  REAL(DP), ALLOCATABLE, PUBLIC :: NodesX_q(:,:)
  REAL(DP), ALLOCATABLE, PUBLIC :: NodesLX1(:), WeightsLX1(:)
  REAL(DP), ALLOCATABLE, PUBLIC :: NodesLX2(:), WeightsLX2(:)
  REAL(DP), ALLOCATABLE, PUBLIC :: NodesLX3(:), WeightsLX3(:)
  REAL(DP), ALLOCATABLE, PUBLIC :: NodesLX_q(:,:)
  REAL(DP), ALLOCATABLE, PUBLIC :: CoarseToFine_X1(:)
  REAL(DP), ALLOCATABLE, PUBLIC :: CoarseToFine_X2(:)
  REAL(DP), ALLOCATABLE, PUBLIC :: CoarseToFine_X3(:)

  PUBLIC :: InitializeReferenceElementX
  PUBLIC :: FinalizeReferenceElementX

CONTAINS


  SUBROUTINE InitializeReferenceElementX

    INTEGER :: iNodeX1, iNodeX2, iNodeX3, iNodeX

    REAL(DP) :: xi_X1, xi_X2, xi_X3

    nDOFX_X1 = nNodesX(2) * nNodesX(3)
    nDOFX_X2 = nNodesX(1) * nNodesX(3)
    nDOFX_X3 = nNodesX(1) * nNodesX(2)

    ALLOCATE &
      ( NodeNumberTableX(3,nDOFX) )
    ALLOCATE &
      ( NodeNumberTableX3D(nNodesX(1),nNodesX(2),nNodesX(3)) )

    iNodeX = 0
    DO iNodeX3 = 1, nNodesX(3)
    DO iNodeX2 = 1, nNodesX(2)
    DO iNodeX1 = 1, nNodesX(1)

      iNodeX = iNodeX + 1

      NodeNumberTableX(1:3,iNodeX) &
        = [ iNodeX1, iNodeX2, iNodeX3 ]

      NodeNumberTableX3D(iNodeX1,iNodeX2,iNodeX3) &
        = iNodeX

    END DO
    END DO
    END DO

    ALLOCATE( NodeNumberTableX_X1(2,nDOFX_X1) )

    iNodeX = 0
    DO iNodeX3 = 1, nNodesX(3)
    DO iNodeX2 = 1, nNodesX(2)

      iNodeX = iNodeX + 1

      NodeNumberTableX_X1(1:2,iNodeX) &
        = [ iNodeX2, iNodeX3 ]

    END DO
    END DO

    ALLOCATE( NodeNumberTableX_X2(2,nDOFX_X2) )

    iNodeX = 0
    DO iNodeX3 = 1, nNodesX(3)
    DO iNodeX1 = 1, nNodesX(1)

      iNodeX = iNodeX + 1

      NodeNumberTableX_X2(1:2,iNodeX) &
        = [ iNodeX1, iNodeX3 ]

    END DO
    END DO

    ALLOCATE( NodeNumberTableX_X3(2,nDOFX_X3) )

    iNodeX = 0
    DO iNodeX2 = 1, nNodesX(2)
    DO iNodeX1 = 1, nNodesX(1)

      iNodeX = iNodeX + 1

      NodeNumberTableX_X3(1:2,iNodeX) &
        = [ iNodeX1, iNodeX2 ]

    END DO
    END DO

    ! --- Gaussian Quadrature Points and Weights ---

    ALLOCATE( NodesX1(nNodesX(1)), WeightsX1(nNodesX(1)) )
    ALLOCATE( NodesX2(nNodesX(2)), WeightsX2(nNodesX(2)) )
    ALLOCATE( NodesX3(nNodesX(3)), WeightsX3(nNodesX(3)) )

    CALL GetQuadrature( nNodesX(1), NodesX1, WeightsX1 )
    CALL GetQuadrature( nNodesX(2), NodesX2, WeightsX2 )
    CALL GetQuadrature( nNodesX(3), NodesX3, WeightsX3 )

    ALLOCATE( WeightsX_X1(nDOFX_X1) )

    DO iNodeX = 1, nDOFX_X1

      iNodeX2 = NodeNumberTableX_X1(1,iNodeX)
      iNodeX3 = NodeNumberTableX_X1(2,iNodeX)

      WeightsX_X1(iNodeX) = WeightsX2(iNodeX2) * WeightsX3(iNodeX3)

    END DO

    ALLOCATE( WeightsX_X2(nDOFX_X2) )

    DO iNodeX = 1, nDOFX_X2

      iNodeX1 = NodeNumberTableX_X2(1,iNodeX)
      iNodeX3 = NodeNumberTableX_X2(2,iNodeX)

      WeightsX_X2(iNodeX) = WeightsX1(iNodeX1) * WeightsX3(iNodeX3)

    END DO

    ALLOCATE( WeightsX_X3(nDOFX_X3) )

    DO iNodeX = 1, nDOFX_X3

      iNodeX1 = NodeNumberTableX_X3(1,iNodeX)
      iNodeX2 = NodeNumberTableX_X3(2,iNodeX)

      WeightsX_X3(iNodeX) = WeightsX1(iNodeX1) * WeightsX2(iNodeX2)

    END DO

    ALLOCATE( WeightsX_q(nDOFX) )
    ALLOCATE( NodesX_q(3,nDOFX) )

    DO iNodeX = 1, nDOFX

      iNodeX1 = NodeNumberTableX(1,iNodeX)
      iNodeX2 = NodeNumberTableX(2,iNodeX)
      iNodeX3 = NodeNumberTableX(3,iNodeX)

      WeightsX_q(iNodeX) &
        = WeightsX1(iNodeX1) * WeightsX2(iNodeX2) * WeightsX3(iNodeX3)

      NodesX_q(1:3,iNodeX) &
        = [ NodesX1(iNodeX1), NodesX2(iNodeX2), NodesX3(iNodeX3) ]

    END DO

    ! --- Lobatto Quadrature Points and Weights ---

    ALLOCATE( NodesLX1(nNodesX(1)), WeightsLX1(nNodesX(1)) )
    ALLOCATE( NodesLX2(nNodesX(2)), WeightsLX2(nNodesX(2)) )
    ALLOCATE( NodesLX3(nNodesX(3)), WeightsLX3(nNodesX(3)) )

    CALL GetQuadrature( nNodesX(1), NodesLX1, WeightsLX1, 'Lobatto' )
    CALL GetQuadrature( nNodesX(2), NodesLX2, WeightsLX2, 'Lobatto' )
    CALL GetQuadrature( nNodesX(3), NodesLX3, WeightsLX3, 'Lobatto' )

    ALLOCATE( NodesLX_q(3,nDOFX) )

    DO iNodeX = 1, nDOFX

      iNodeX1 = NodeNumberTableX(1,iNodeX)
      iNodeX2 = NodeNumberTableX(2,iNodeX)
      iNodeX3 = NodeNumberTableX(3,iNodeX)

      NodesLX_q(1:3,iNodeX) &
        = [ NodesLX1(iNodeX1), NodesLX2(iNodeX2), NodesLX3(iNodeX3) ]

    END DO

    ! --- Coarse element points interpolated to fine element faces ---

    ALLOCATE( CoarseToFine_X1(2**(nDimsX)*nDOFX_X1) )
    ALLOCATE( CoarseToFine_X2(2**(nDimsX)*nDOFX_X2) )
    ALLOCATE( CoarseToFine_X3(2**(nDimsX)*nDOFX_X3) )

    CALL ComputeCoarseToFine_X1( CoarseToFine_X1 )
    CALL ComputeCoarseToFine_X2( CoarseToFine_X2 )
    CALL ComputeCoarseToFine_X3( CoarseToFine_X3 )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: WeightsX_q, WeightsX_X1, WeightsX_X2, WeightsX_X3 )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( WeightsX_q, WeightsX_X1, WeightsX_X2, WeightsX_X3 )
#endif

  END SUBROUTINE InitializeReferenceElementX


  SUBROUTINE FinalizeReferenceElementX

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: WeightsX_q, WeightsX_X1, WeightsX_X2, WeightsX_X3 )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( WeightsX_q, WeightsX_X1, WeightsX_X2, WeightsX_X3 )
#endif

    DEALLOCATE( CoarseToFine_X3 )
    DEALLOCATE( CoarseToFine_X2 )
    DEALLOCATE( CoarseToFine_X1 )
    DEALLOCATE( NodeNumberTableX )
    DEALLOCATE( NodeNumberTableX_X1 )
    DEALLOCATE( NodeNumberTableX_X2 )
    DEALLOCATE( NodeNumberTableX_X3 )
    DEALLOCATE( NodeNumberTableX3D )
    DEALLOCATE( NodesX1, WeightsX1 )
    DEALLOCATE( NodesX2, WeightsX2 )
    DEALLOCATE( NodesX3, WeightsX3 )
    DEALLOCATE( WeightsX_X1 )
    DEALLOCATE( WeightsX_X2 )
    DEALLOCATE( WeightsX_X3 )
    DEALLOCATE( WeightsX_q )
    DEALLOCATE( NodesX_q )

  END SUBROUTINE FinalizeReferenceElementX


  ! --- PRIVATE SUBROUTINES ---


  REAL(DP) FUNCTION EtaOfXi_1D( k, xi ) RESULT( eta )

    INTEGER,  INTENT(in) :: k
    REAL(DP), INTENT(in) :: xi

    eta = Half * ( xi + Half * ( -One )**( k ) )

    RETURN
  END FUNCTION EtaOfXi_1D


  SUBROUTINE ComputeCoarseToFine_X1( CoarseToFine_X1 )

    REAL(DP), INTENT(out) :: CoarseToFine_X1(:)

    INTEGER :: iNodeX_X1, iNodeX2, iNodeX3

    DO iNodeX3 = 1, nNodesX(3)
    DO iNodeX2 = 1, nNodesX(2)

      iNodeX_X1 = iNodeX2 + iNodeX3 - 1

      IF( nDimsX .EQ. 1 )THEN

        CoarseToFine_X1(iNodeX_X1)          = EtaOfXi_1D( 1, -Half )
        CoarseToFine_X1(iNodeX_X1+nDOFX_X1) = EtaOfXi_1D( 2, +Half )

      ELSE IF( nDimsX .EQ. 2 )THEN

        CoarseToFine_X1(iNodeX_X1+0*nDOFX_X1) &
          =   EtaOfXi_1D( 1, -Half ) &
            * EtaOfXi_1D( 1, NodesX2(iNodeX2) )

        CoarseToFine_X1(iNodeX_X1+1*nDOFX_X1) &
          =   EtaOfXi_1D( 2, +Half ) &
            * EtaOfXi_1D( 1, NodesX2(iNodeX2) )

        CoarseToFine_X1(iNodeX_X1+2*nDOFX_X1) &
          =   EtaOfXi_1D( 1, -Half ) &
            * EtaOfXi_1D( 2, NodesX2(iNodeX2) )

        CoarseToFine_X1(iNodeX_X1+3*nDOFX_X1) &
          =   EtaOfXi_1D( 2, +Half ) &
            * EtaOfXi_1D( 2, NodesX2(iNodeX2) )

      ELSE

        CoarseToFine_X1(iNodeX_X1+0*nDOFX_X1) &
          =   EtaOfXi_1D( 1, -Half ) &
            * EtaOfXi_1D( 1, NodesX2(iNodeX2) ) &
            * EtaOfXi_1D( 1, NodesX3(iNodeX3) )

        CoarseToFine_X1(iNodeX_X1+1*nDOFX_X1) &
          =   EtaOfXi_1D( 2, +Half ) &
            * EtaOfXi_1D( 1, NodesX2(iNodeX2) ) &
            * EtaOfXi_1D( 1, NodesX3(iNodeX3) )

        CoarseToFine_X1(iNodeX_X1+2*nDOFX_X1) &
          =   EtaOfXi_1D( 1, -Half ) &
            * EtaOfXi_1D( 2, NodesX2(iNodeX2) ) &
            * EtaOfXi_1D( 1, NodesX3(iNodeX3) )

        CoarseToFine_X1(iNodeX_X1+3*nDOFX_X1) &
          =   EtaOfXi_1D( 2, +Half ) &
            * EtaOfXi_1D( 2, NodesX2(iNodeX2) ) &
            * EtaOfXi_1D( 1, NodesX3(iNodeX3) )

        CoarseToFine_X1(iNodeX_X1+4*nDOFX_X1) &
          =   EtaOfXi_1D( 1, -Half ) &
            * EtaOfXi_1D( 1, NodesX2(iNodeX2) ) &
            * EtaOfXi_1D( 2, NodesX3(iNodeX3) )

        CoarseToFine_X1(iNodeX_X1+5*nDOFX_X1) &
          =   EtaOfXi_1D( 2, +Half ) &
            * EtaOfXi_1D( 1, NodesX2(iNodeX2) ) &
            * EtaOfXi_1D( 2, NodesX3(iNodeX3) )

        CoarseToFine_X1(iNodeX_X1+6*nDOFX_X1) &
          =   EtaOfXi_1D( 1, -Half ) &
            * EtaOfXi_1D( 2, NodesX2(iNodeX2) ) &
            * EtaOfXi_1D( 2, NodesX3(iNodeX3) )

        CoarseToFine_X1(iNodeX_X1+7*nDOFX_X1) &
          =   EtaOfXi_1D( 2, +Half ) &
            * EtaOfXi_1D( 2, NodesX2(iNodeX2) ) &
            * EtaOfXi_1D( 2, NodesX3(iNodeX3) )

      END IF

    END DO
    END DO

  END SUBROUTINE ComputeCoarseToFine_X1


  SUBROUTINE ComputeCoarseToFine_X2( CoarseToFine_X2 )

    REAL(DP), INTENT(out) :: CoarseToFine_X2(:)

    INTEGER :: iNodeX_X2, iNodeX1, iNodeX3, iNodeX

    DO iNodeX_X2 = 1, nDOFX_X2

      IF( nDimsX .EQ. 1 )THEN

        CoarseToFine_X2(iNodeX_X2) = One

      ELSE IF( nDimsX .EQ. 2 )THEN

        DO iNodeX1 = 1, nNodesX(1)

          CoarseToFine_X2(iNodeX1+0*nDOFX_X2) &
            =   EtaOfXi_1D( 1, NodesX1(iNodeX1) ) &
              * EtaOfXi_1D( 1, -Half )

          CoarseToFine_X2(iNodeX1+1*nDOFX_X2) &
            =   EtaOfXi_1D( 2, NodesX1(iNodeX1) ) &
              * EtaOfXi_1D( 1, -Half )

          CoarseToFine_X2(iNodeX1+2*nDOFX_X2) &
            =   EtaOfXi_1D( 1, NodesX1(iNodeX1) ) &
              * EtaOfXi_1D( 2, +Half )

          CoarseToFine_X2(iNodeX1+3*nDOFX_X2) &
            =   EtaOfXi_1D( 2, NodesX1(iNodeX1) ) &
              * EtaOfXi_1D( 2, +Half )

        END DO

      ELSE

        iNodeX = 0

        DO iNodeX3 = 1, nNodesX(3)
        DO iNodeX1 = 1, nNodesX(1)

          iNodeX = iNodeX + 1

          CoarseToFine_X2(iNodeX+0*nDOFX_X2) &
            =   EtaOfXi_1D( 1, NodesX1(iNodeX1) ) &
              * EtaOfXi_1D( 1, -Half ) &
              * EtaOfXi_1D( 1, NodesX3(iNodeX3) )

          CoarseToFine_X2(iNodeX+1*nDOFX_X2) &
            =   EtaOfXi_1D( 2, NodesX1(iNodeX1) ) &
              * EtaOfXi_1D( 1, -Half ) &
              * EtaOfXi_1D( 1, NodesX3(iNodeX3) )

          CoarseToFine_X2(iNodeX+2*nDOFX_X2) &
            =   EtaOfXi_1D( 1, NodesX1(iNodeX1) ) &
              * EtaOfXi_1D( 2, +Half ) &
              * EtaOfXi_1D( 1, NodesX3(iNodeX3) )

          CoarseToFine_X2(iNodeX+3*nDOFX_X2) &
            =   EtaOfXi_1D( 2, NodesX1(iNodeX1) ) &
              * EtaOfXi_1D( 2, +Half ) &
              * EtaOfXi_1D( 1, NodesX3(iNodeX3) )

          CoarseToFine_X2(iNodeX+4*nDOFX_X2) &
            =   EtaOfXi_1D( 1, NodesX1(iNodeX1) ) &
              * EtaOfXi_1D( 1, -Half ) &
              * EtaOfXi_1D( 2, NodesX3(iNodeX3) )

          CoarseToFine_X2(iNodeX+5*nDOFX_X2) &
            =   EtaOfXi_1D( 2, NodesX1(iNodeX1) ) &
              * EtaOfXi_1D( 1, -Half ) &
              * EtaOfXi_1D( 2, NodesX3(iNodeX3) )

          CoarseToFine_X2(iNodeX+6*nDOFX_X2) &
            =   EtaOfXi_1D( 1, NodesX1(iNodeX1) ) &
              * EtaOfXi_1D( 2, +Half ) &
              * EtaOfXi_1D( 2, NodesX3(iNodeX3) )

          CoarseToFine_X2(iNodeX+7*nDOFX_X2) &
            =   EtaOfXi_1D( 2, NodesX1(iNodeX1) ) &
              * EtaOfXi_1D( 2, +Half ) &
              * EtaOfXi_1D( 2, NodesX3(iNodeX3) )

        END DO
        END DO

      END IF

    END DO

  END SUBROUTINE ComputeCoarseToFine_X2


  SUBROUTINE ComputeCoarseToFine_X3( CoarseToFine_X3 )

    REAL(DP), INTENT(out) :: CoarseToFine_X3(:)

    INTEGER :: iNodeX_X3, iNodeX1, iNodeX2, iNodeX

    DO iNodeX_X3 = 1, nDOFX_X3

      IF( nDimsX .EQ. 1 )THEN

        CoarseToFine_X3 = One

      ELSE IF( nDimsX .EQ. 2 )THEN

        CoarseToFine_X3 = One

      ELSE

        iNodeX = 0

        DO iNodeX2 = 1, nNodesX(2)
        DO iNodeX1 = 1, nNodesX(1)

          iNodeX = iNodeX + 1

          CoarseToFine_X3(iNodeX+0*nDOFX_X3) &
            =   EtaOfXi_1D( 1, NodesX1(iNodeX1) ) &
              * EtaOfXi_1D( 1, NodesX2(iNodeX2) ) &
              * EtaOfXi_1D( 1, +Half )

          CoarseToFine_X3(iNodeX+1*nDOFX_X3) &
            =   EtaOfXi_1D( 2, NodesX1(iNodeX1) ) &
              * EtaOfXi_1D( 1, NodesX2(iNodeX2) ) &
              * EtaOfXi_1D( 1, +Half )

          CoarseToFine_X3(iNodeX+2*nDOFX_X3) &
            =   EtaOfXi_1D( 1, NodesX1(iNodeX1) ) &
              * EtaOfXi_1D( 2, NodesX2(iNodeX2) ) &
              * EtaOfXi_1D( 1, +Half )

          CoarseToFine_X3(iNodeX+3*nDOFX_X3) &
            =   EtaOfXi_1D( 2, NodesX1(iNodeX1) ) &
              * EtaOfXi_1D( 2, NodesX2(iNodeX2) ) &
              * EtaOfXi_1D( 1, +Half )

          CoarseToFine_X3(iNodeX+4*nDOFX_X3) &
            =   EtaOfXi_1D( 1, NodesX1(iNodeX1) ) &
              * EtaOfXi_1D( 1, NodesX2(iNodeX2) ) &
              * EtaOfXi_1D( 2, -Half )

          CoarseToFine_X3(iNodeX+5*nDOFX_X3) &
            =   EtaOfXi_1D( 2, NodesX1(iNodeX1) ) &
              * EtaOfXi_1D( 1, NodesX2(iNodeX2) ) &
              * EtaOfXi_1D( 2, -Half )

          CoarseToFine_X3(iNodeX+6*nDOFX_X3) &
            =   EtaOfXi_1D( 1, NodesX1(iNodeX1) ) &
              * EtaOfXi_1D( 2, NodesX2(iNodeX2) ) &
              * EtaOfXi_1D( 2, -Half )

          CoarseToFine_X3(iNodeX+7*nDOFX_X3) &
            =   EtaOfXi_1D( 2, NodesX1(iNodeX1) ) &
              * EtaOfXi_1D( 2, NodesX2(iNodeX2) ) &
              * EtaOfXi_1D( 2, -Half )

        END DO
        END DO

      END IF

    END DO

  END SUBROUTINE ComputeCoarseToFine_X3


END MODULE ReferenceElementModuleX
