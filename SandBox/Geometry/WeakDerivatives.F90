MODULE WeakDerivatives

  USE KindModule, ONLY: &
    DP, Zero, Half, One, Two, &
    SqrtTiny
  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    nDOFE, &
    nX,    &
    nDOFZ
  USE LinearAlgebraModule, ONLY: &
    MatrixMatrixMultiply
  USE ReferenceElementModuleX, ONLY: &
    nDOFX_X1, nDOFX_X2, nDOFX_X3, &
    WeightsX_q, &
    WeightsX_X1, &
    WeightsX_X2, &
    WeightsX_X3
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    dLXdX1_q, LX_X1_Dn, LX_X1_Up, &
    dLXdX2_q, LX_X2_Dn, LX_X2_Up, &
    dLXdX3_q, LX_X3_Dn, LX_X3_Up
  USE ReferenceElementModule, ONLY: &
    nDOF_X1, nDOF_X2, nDOF_X3, &
    Weights_q, &
    Weights_X1, &
    Weights_X2, &
    Weights_X3
  USE ReferenceElementModule_Lagrange, ONLY: &
    L_E_Dn,  L_E_Up, &
    dLdE_q, &
    L_X1_Dn, L_X1_Up, &
    dLdX1_q, &
    L_X2_Dn, L_X2_Up, &
    dLdX2_q, &
    L_X3_Dn, L_X3_Up, &
    dLdX3_q
  USE MeshModule, ONLY: &
    MeshE, &
    NodeCoordinate, &
    MeshX
  USE GeometryFieldsModule, ONLY: &
    nGF, &
    iGF_h_1, &
    iGF_h_2, &
    iGF_h_3, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_SqrtGm, &
    iGF_Alpha, &
    iGF_Beta_1, &
    iGF_Beta_2, &
    iGF_Beta_3
  USE FluidFieldsModule, ONLY: &
    nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne
  USE TwoMoment_ClosureModule, ONLY: &
    FluxFactor, &
    EddingtonFactor
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE Euler_UtilitiesModule_Relativistic, ONLY: &
    ComputePrimitive_Euler_Relativistic
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PARAMETER :: UpwindFlux = .FALSE.


  PUBLIC :: ComputeWeakDerivatives_X0
  PUBLIC :: ComputeWeakDerivatives_X1
  PUBLIC :: ComputeWeakDerivatives_X2
  PUBLIC :: ComputeWeakDerivatives_X3
  PUBLIC :: ComputeChristoffel

CONTAINS



  FUNCTION FaceVelocity_X1 &
    ( V1_L, V2_L, V3_L, V1_R, V2_R, V3_R )

    REAL(DP), INTENT(in) :: V1_L, V2_L, V3_L
    REAL(DP), INTENT(in) :: V1_R, V2_R, V3_R
    REAL(DP)             :: FaceVelocity_X1(1:3)

    ! --- Average Left and Right States ---

    FaceVelocity_X1(1) = Half * ( V1_L + V1_R )
    FaceVelocity_X1(2) = Half * ( V2_L + V2_R )
    FaceVelocity_X1(3) = Half * ( V3_L + V3_R )

    RETURN
  END FUNCTION FaceVelocity_X1

  FUNCTION FaceVelocity_X2 &
    ( V1_L, V2_L, V3_L, V1_R, V2_R, V3_R )

    REAL(DP), INTENT(in) :: V1_L, V2_L, V3_L
    REAL(DP), INTENT(in) :: V1_R, V2_R, V3_R
    REAL(DP)             :: FaceVelocity_X2(1:3)

    ! --- Average Left and Right States ---

    FaceVelocity_X2(1) = Half * ( V1_L + V1_R )
    FaceVelocity_X2(2) = Half * ( V2_L + V2_R )
    FaceVelocity_X2(3) = Half * ( V3_L + V3_R )

    RETURN
  END FUNCTION FaceVelocity_X2

  FUNCTION FaceVelocity_X3 &
    ( V1_L, V2_L, V3_L, V1_R, V2_R, V3_R )

    REAL(DP), INTENT(in) :: V1_L, V2_L, V3_L
    REAL(DP), INTENT(in) :: V1_R, V2_R, V3_R
    REAL(DP)             :: FaceVelocity_X3(1:3)

    ! --- Average Left and Right States ---

    FaceVelocity_X3(1) = Half * ( V1_L + V1_R )
    FaceVelocity_X3(2) = Half * ( V2_L + V2_R )
    FaceVelocity_X3(3) = Half * ( V3_L + V3_R )

    RETURN
  END FUNCTION FaceVelocity_X3

  


  SUBROUTINE ComputeWeakDerivatives_X0 &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, dG_dd_dX0, Verbose_Option  )

    INTEGER,  INTENT(in)  :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)  :: &
      GX (1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(out) :: &
      dG_dd_dX0 &
         (1:nDOFX,0:3,0:3, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)) 
    LOGICAL,          INTENT(in), OPTIONAL :: Verbose_Option

    dG_dd_dX0 = Zero
  END SUBROUTINE ComputeWeakDerivatives_X0 


  SUBROUTINE ComputeWeakDerivatives_X1 &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, dG_dd_dX1  )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)  :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)  :: &
      GX (1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(out) :: &
      dG_dd_dX1 &
         (1:nDOFX,0:3,0:3, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)) 

    INTEGER  :: nK(4), nK_X1(4), nX, nX_X1
    INTEGER  :: iNodeX, iNodeX1
    INTEGER  :: i, iZ2, iZ3, iZ4, iCF, iGF
    REAL(DP) :: B_u_X1(3), B_d_X1(3), A_X1, B_u_K(3), B_d_K(3), A_K
    REAL(DP) :: &
      G_munu_K(nDOFX   ,7,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
           iZ_B1(2):iZ_E1(2)), &
      G_munu_F(nDOFX_X1,7,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
           iZ_B0(2):iZ_E0(2)+1), &
      H_munu_K(nDOFX   ,7,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
           iZ_B0(2):iZ_E0(2)), &
      H_munu_F(nDOFX_X1,7,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
           iZ_B0(2):iZ_E0(2)+1)
    REAL(DP) :: & 
      dG_dd_dX1_Temp &
         (1:nDOFX,7, &
          iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),iZ_B0(2):iZ_E0(2)) 
    CHARACTER(len=40) :: name1, name2
    CHARACTER(len=1):: nds
    CHARACTER(len=2)::nxn1
    CHARACTER(len=3)::nxn2
    LOGICAL :: Verbose

    IF( iZ_E0(2) .EQ. iZ_B0(2) )THEN
      dG_dd_dX1 = Zero
      RETURN
    END IF

    PRINT*, "      ComputeWeakDerivatives_X1"
    

    nK    = iZ_E0 - iZ_B0 + 1 ! Number of Elements per Phase Space Dimension
    nK_X1 = nK + [0,1,0,0]    ! Number of X1 Faces per Phase Space Dimension
    nX    = PRODUCT( nK   (2:4) ) ! Number of Elements in Position Space
    nX_X1 = PRODUCT( nK_X1(2:4) ) ! Number of X1 Faces in Position Space
    
    ! --- Permute Geometry Fields ---


    DO iZ2 = iZ_B1(2), iZ_E1(2)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)

      DO iNodeX = 1, nDOFX





        A_K = GX(iNodeX,iZ2,iZ3,iZ4,iGF_Alpha)

        B_u_K(1) =  GX(iNodeX,iZ2,iZ3,iZ4,iGF_Beta_1)
        B_u_K(2) =  GX(iNodeX,iZ2,iZ3,iZ4,iGF_Beta_2)
        B_u_K(3) =  GX(iNodeX,iZ2,iZ3,iZ4,iGF_Beta_3)

        B_d_K(1) =  GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_11) * B_u_K(1)
        B_d_K(2) =  GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_22) * B_u_K(2)
        B_d_K(3) =  GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_33) * B_u_K(3)

        G_munu_K(iNodeX,1,iZ3,iZ4,iZ2) = -A_K**2 + B_d_K(1) * B_u_K(1) + B_d_K(2) * B_u_K(2) + B_d_K(3) * B_u_K(3)

        G_munu_K(iNodeX,2,iZ3,iZ4,iZ2) = B_d_K(1)

        G_munu_K(iNodeX,3,iZ3,iZ4,iZ2) = B_d_K(2) 
        
        G_munu_K(iNodeX,4,iZ3,iZ4,iZ2) = B_d_K(3) 
        
        G_munu_K(iNodeX,5,iZ3,iZ4,iZ2) = GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_11)

        G_munu_K(iNodeX,6,iZ3,iZ4,iZ2) = GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_22)

        G_munu_K(iNodeX,7,iZ3,iZ4,iZ2) = GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_33)
      END DO

    END DO
    END DO
    END DO
    ! --- Interpolate Geometry Fields ---







      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X1, 7*nX_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
               G_munu_K(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)-1), nDOFX, Zero, &
               G_munu_F(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)), nDOFX_X1 )

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X1, 7*nX_X1, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
               G_munu_K(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)), nDOFX, Half, &
               G_munu_F(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)), nDOFX_X1 )
    ! --- Compute Metric Components ---

    ! --- Permute Fluid Fields ---


    ! --- Interpolate Fluid Fields ---

  
    DO iZ2 = iZ_B0(2), iZ_E0(2) + 1
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)

      DO iNodeX = 1, nDOFX_X1




        H_munu_F(iNodeX,1,iZ3,iZ4,iZ2) = WeightsX_X1(iNodeX) * G_munu_F(iNodeX,1,iZ3,iZ4,iZ2)
        
        H_munu_F(iNodeX,2,iZ3,iZ4,iZ2) = WeightsX_X1(iNodeX) * G_munu_F(iNodeX,2,iZ3,iZ4,iZ2)
        
        H_munu_F(iNodeX,3,iZ3,iZ4,iZ2) = WeightsX_X1(iNodeX) * G_munu_F(iNodeX,3,iZ3,iZ4,iZ2)
        
        H_munu_F(iNodeX,4,iZ3,iZ4,iZ2) = WeightsX_X1(iNodeX) * G_munu_F(iNodeX,4,iZ3,iZ4,iZ2)

        H_munu_F(iNodeX,5,iZ3,iZ4,iZ2) = WeightsX_X1(iNodeX) * G_munu_F(iNodeX,5,iZ3,iZ4,iZ2)

        H_munu_F(iNodeX,6,iZ3,iZ4,iZ2) = WeightsX_X1(iNodeX) * G_munu_F(iNodeX,6,iZ3,iZ4,iZ2)

        H_munu_F(iNodeX,7,iZ3,iZ4,iZ2) = WeightsX_X1(iNodeX) * G_munu_F(iNodeX,7,iZ3,iZ4,iZ2)
     END DO

    END DO
    END DO
    END DO


    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 7*nX, nDOFX_X1, - One, LX_X1_Dn, nDOFX_X1, &
             H_munu_F(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)  ), nDOFX_X1, Zero, &
             dG_dd_dX1_Temp, nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 7*nX, nDOFX_X1, + One, LX_X1_Up, nDOFX_X1, &
             H_munu_F(1,1,iZ_B0(3),iZ_B0(4),iZ_B0(2)+1), nDOFX_X1, One,  &
             dG_dd_dX1_Temp, nDOFX )

    ! --- Volume Term ---
    ! -------------------

    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)

      DO iNodeX = 1, nDOFX



        H_munu_K(iNodeX,1,iZ3,iZ4,iZ2) = WeightsX_q(iNodeX) * G_munu_K(iNodeX,1,iZ3,iZ4,iZ2)

        H_munu_K(iNodeX,2,iZ3,iZ4,iZ2) = WeightsX_q(iNodeX) * G_munu_K(iNodeX,2,iZ3,iZ4,iZ2)

        H_munu_K(iNodeX,3,iZ3,iZ4,iZ2) = WeightsX_q(iNodeX) * G_munu_K(iNodeX,3,iZ3,iZ4,iZ2)

        H_munu_K(iNodeX,4,iZ3,iZ4,iZ2) = WeightsX_q(iNodeX) * G_munu_K(iNodeX,4,iZ3,iZ4,iZ2)

        H_munu_K(iNodeX,5,iZ3,iZ4,iZ2) = WeightsX_q(iNodeX) * G_munu_K(iNodeX,5,iZ3,iZ4,iZ2)

        H_munu_K(iNodeX,6,iZ3,iZ4,iZ2) = WeightsX_q(iNodeX) * G_munu_K(iNodeX,6,iZ3,iZ4,iZ2)

        H_munu_K(iNodeX,7,iZ3,iZ4,iZ2) = WeightsX_q(iNodeX) * G_munu_K(iNodeX,7,iZ3,iZ4,iZ2)


      END DO

    END DO
    END DO
    END DO


    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 7*nX, nDOFX, - One, dLXdX1_q, nDOFX, &
             H_munu_K, nDOFX, One, dG_dd_dX1_Temp, nDOFX )

    ASSOCIATE( dZ2 => MeshX(1) % Width )

    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)

      DO iNodeX = 1, nDOFX


        dG_dd_dX1_Temp(iNodeX,1,iZ3,iZ4,iZ2) &
         = dG_dd_dX1_Temp(iNodeX,1,iZ3,iZ4,iZ2) &
             / ( WeightsX_q(iNodeX) * dZ2(iZ2) )

        dG_dd_dX1_Temp(iNodeX,2,iZ3,iZ4,iZ2) &
         = dG_dd_dX1_Temp(iNodeX,2,iZ3,iZ4,iZ2) &
             / ( WeightsX_q(iNodeX) * dZ2(iZ2) )

        dG_dd_dX1_Temp(iNodeX,3,iZ3,iZ4,iZ2) &
         = dG_dd_dX1_Temp(iNodeX,3,iZ3,iZ4,iZ2) &
             / ( WeightsX_q(iNodeX) * dZ2(iZ2) )

        dG_dd_dX1_Temp(iNodeX,4,iZ3,iZ4,iZ2) &
         = dG_dd_dX1_Temp(iNodeX,4,iZ3,iZ4,iZ2) &
             / ( WeightsX_q(iNodeX) * dZ2(iZ2) )


        dG_dd_dX1_Temp(iNodeX,5,iZ3,iZ4,iZ2) &
         = dG_dd_dX1_Temp(iNodeX,5,iZ3,iZ4,iZ2) &
             / ( WeightsX_q(iNodeX) * dZ2(iZ2) )

        dG_dd_dX1_Temp(iNodeX,6,iZ3,iZ4,iZ2) &
         = dG_dd_dX1_Temp(iNodeX,6,iZ3,iZ4,iZ2) &
             / ( WeightsX_q(iNodeX) * dZ2(iZ2) )

        dG_dd_dX1_Temp(iNodeX,7,iZ3,iZ4,iZ2) &
         = dG_dd_dX1_Temp(iNodeX,7,iZ3,iZ4,iZ2) &
             / ( WeightsX_q(iNodeX) * dZ2(iZ2) )



      END DO

    END DO
    END DO
    END DO
    END ASSOCIATE

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX



        dG_dd_dX1(iNodeX,0,0,iZ2,iZ3,iZ4) = dG_dd_dX1_Temp(iNodeX,1,iZ3,iZ4,iZ2) 

        dG_dd_dX1(iNodeX,0,1,iZ2,iZ3,iZ4) = dG_dd_dX1_Temp(iNodeX,2,iZ3,iZ4,iZ2) 

        dG_dd_dX1(iNodeX,0,2,iZ2,iZ3,iZ4) = dG_dd_dX1_Temp(iNodeX,3,iZ3,iZ4,iZ2) 

        dG_dd_dX1(iNodeX,0,3,iZ2,iZ3,iZ4) = dG_dd_dX1_Temp(iNodeX,4,iZ3,iZ4,iZ2) 

        dG_dd_dX1(iNodeX,1,0,iZ2,iZ3,iZ4) = dG_dd_dX1_Temp(iNodeX,2,iZ3,iZ4,iZ2) 

        dG_dd_dX1(iNodeX,2,0,iZ2,iZ3,iZ4) = dG_dd_dX1_Temp(iNodeX,3,iZ3,iZ4,iZ2) 

        dG_dd_dX1(iNodeX,3,0,iZ2,iZ3,iZ4) = dG_dd_dX1_Temp(iNodeX,4,iZ3,iZ4,iZ2) 

        dG_dd_dX1(iNodeX,1,1,iZ2,iZ3,iZ4) = dG_dd_dX1_Temp(iNodeX,5,iZ3,iZ4,iZ2) 

        dG_dd_dX1(iNodeX,2,2,iZ2,iZ3,iZ4) = dG_dd_dX1_Temp(iNodeX,6,iZ3,iZ4,iZ2) 

        dG_dd_dX1(iNodeX,3,3,iZ2,iZ3,iZ4) = dG_dd_dX1_Temp(iNodeX,7,iZ3,iZ4,iZ2) 

        dG_dd_dX1(iNodeX,1,2,iZ2,iZ3,iZ4) = 0.0_DP

        dG_dd_dX1(iNodeX,2,1,iZ2,iZ3,iZ4) = 0.0_DP 

        dG_dd_dX1(iNodeX,1,3,iZ2,iZ3,iZ4) = 0.0_DP

        dG_dd_dX1(iNodeX,3,1,iZ2,iZ3,iZ4) = 0.0_DP

        dG_dd_dX1(iNodeX,2,3,iZ2,iZ3,iZ4) = 0.0_DP

        dG_dd_dX1(iNodeX,3,2,iZ2,iZ3,iZ4) = 0.0_DP






!print*,iZ2, dG_dd_dX1(iNodeX,2,2,iZ2,iZ3,iZ4)

      END DO

    END DO
    END DO
    END DO
!STOP
  END SUBROUTINE ComputeWeakDerivatives_X1

  SUBROUTINE ComputeWeakDerivatives_X2 &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, dG_dd_dX2, Verbose_Option  )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)  :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)  :: &
      GX (1:nDOFX, &
          iZ_B1(2):iZ_E1(2), &
          iZ_B1(3):iZ_E1(3), &
          iZ_B1(4):iZ_E1(4), &
          1:nGF)
    REAL(DP), INTENT(out) :: &
      dG_dd_dX2 &
         (1:nDOFX,0:3,0:3, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)) 
    LOGICAL,          INTENT(in), OPTIONAL :: Verbose_Option

    INTEGER  :: iNodeX
    INTEGER  :: i, iZ2, iZ3, iZ4, iGF, iCF
    INTEGER  :: nK(4), nK_X2(4), nX, nX_X2
    REAL(DP) :: uPF_L(nPF), uPF_R(nPF), uPF_K(nPF)
    REAL(DP) :: &
      GX_K(nDOFX,nGF, &
           iZ_B0(2):iZ_E0(2), &
           iZ_B0(4):iZ_E0(4), &
           iZ_B1(3):iZ_E1(3)), &
      G_munu_K(nDOFX   ,7,iZ_B0(2):iZ_E0(2),iZ_B0(4):iZ_E0(4), &
           iZ_B1(3):iZ_E1(3)), &
      G_munu_F(nDOFX_X2,7,iZ_B0(2):iZ_E0(2),iZ_B0(4):iZ_E0(4), &
           iZ_B0(3):iZ_E0(3)+1), &
      H_munu_K(nDOFX   ,7,iZ_B0(2):iZ_E0(2),iZ_B0(4):iZ_E0(4), &
           iZ_B0(3):iZ_E0(3)), &
      H_munu_F(nDOFX_X2,7,iZ_B0(2):iZ_E0(2),iZ_B0(4):iZ_E0(4), &
           iZ_B0(3):iZ_E0(3)+1)
    REAL(DP) :: &
      GX_F(nDOFX_X2,nGF, &
           iZ_B0(2):iZ_E0(2), &
           iZ_B0(4):iZ_E0(4), &
           iZ_B0(3):iZ_E1(3))
    REAL(DP) :: B_u_X2(3), B_d_X2(3), A_X2, B_u_K(3), B_d_K(3), A_K
    REAL(DP) :: & 
      dG_dd_dX2_Temp &
         (1:nDOFX,7, &
          iZ_B0(2):iZ_E0(2),iZ_B0(4):iZ_E0(4),iZ_B0(3):iZ_E0(3)) 
    LOGICAL :: Verbose

    IF( iZ_E0(3) .EQ. iZ_B0(3) )THEN
      dG_dd_dX2 = Zero
      RETURN
    END IF

    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    IF (Verbose) THEN
    PRINT*, "      ComputeWeakDerivatives_X2"
    END IF

    nK    = iZ_E0 - iZ_B0 + 1 ! Number of Elements per Phase Space Dimension
    nK_X2 = nK + [0,0,1,0]    ! Number of X2 Faces per Phase Space Dimension
    nX    = PRODUCT( nK   (2:4) ) ! Number of Elements in Position Space
    nX_X2 = PRODUCT( nK_X2(2:4) ) ! Number of X2 Faces in Position Space

    ! --- Permute Geometry Fields ---


    DO iZ3 = iZ_B1(3), iZ_E1(3)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX



        A_K = GX(iNodeX,iZ2,iZ3,iZ4,iGF_Alpha)

        B_u_K(1) =  GX(iNodeX,iZ2,iZ3,iZ4,iGF_Beta_1)
        B_u_K(2) =  GX(iNodeX,iZ2,iZ3,iZ4,iGF_Beta_2)
        B_u_K(3) =  GX(iNodeX,iZ2,iZ3,iZ4,iGF_Beta_3)

        B_d_K(1) =  GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_11) * B_u_K(1)
        B_d_K(2) =  GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_22) * B_u_K(2)
        B_d_K(3) =  GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_33) * B_u_K(3)

        G_munu_K(iNodeX,1,iZ2,iZ4,iZ3) = -A_K**2 + B_d_K(1) * B_u_K(1) + B_d_K(2) * B_u_K(2) + B_d_K(3) * B_u_K(3)

        G_munu_K(iNodeX,2,iZ2,iZ4,iZ3) = B_d_K(1)

        G_munu_K(iNodeX,3,iZ2,iZ4,iZ3) = B_d_K(2) 
        
        G_munu_K(iNodeX,4,iZ2,iZ4,iZ3) = B_d_K(3) 
        
        G_munu_K(iNodeX,5,iZ2,iZ4,iZ3) = GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_11)

        G_munu_K(iNodeX,6,iZ2,iZ4,iZ3) = GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_22)

        G_munu_K(iNodeX,7,iZ2,iZ4,iZ3) = GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_33)

      END DO

    END DO
    END DO
    END DO


    !---------------------
    ! --- Surface Term ---
    !---------------------


    ! --- Interpolate Geometry Fields on Shared Face ---

    ! --- Face States (Average of Left and Right States) ---











      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X2, 7*nX_X2, nDOFX, One,  LX_X2_Up, nDOFX_X2, &
               G_munu_K(1,1,iZ_B0(2),iZ_B0(4),iZ_B0(3)-1), nDOFX, Zero, &
               G_munu_F(1,1,iZ_B0(2),iZ_B0(4),iZ_B0(3)), nDOFX_X2 )

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X2, 7*nX_X2, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
               G_munu_K(1,1,iZ_B0(2),iZ_B0(4),iZ_B0(3)), nDOFX, Half, &
               G_munu_F(1,1,iZ_B0(2),iZ_B0(4),iZ_B0(3)), nDOFX_X2 )



    ! --- Compute Metric Components from Scale Factors ---



    ! --- Permute Fluid Fields ---




    ! --- Interpolate Fluid Fields ---


    ! --- Interpolate Left State ---



    DO iZ3 = iZ_B0(3), iZ_E1(3)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX_X2



        H_munu_F(iNodeX,1,iZ2,iZ4,iZ3) = WeightsX_X2(iNodeX) * G_munu_F(iNodeX,1,iZ2,iZ4,iZ3)
        
        H_munu_F(iNodeX,2,iZ2,iZ4,iZ3) = WeightsX_X2(iNodeX) * G_munu_F(iNodeX,2,iZ2,iZ4,iZ3)
        
        H_munu_F(iNodeX,3,iZ2,iZ4,iZ3) = WeightsX_X2(iNodeX) * G_munu_F(iNodeX,3,iZ2,iZ4,iZ3)
        
        H_munu_F(iNodeX,4,iZ2,iZ4,iZ3) = WeightsX_X2(iNodeX) * G_munu_F(iNodeX,4,iZ2,iZ4,iZ3)

        H_munu_F(iNodeX,5,iZ2,iZ4,iZ3) = WeightsX_X2(iNodeX) * G_munu_F(iNodeX,5,iZ2,iZ4,iZ3)

        H_munu_F(iNodeX,6,iZ2,iZ4,iZ3) = WeightsX_X2(iNodeX) * G_munu_F(iNodeX,6,iZ2,iZ4,iZ3)

        H_munu_F(iNodeX,7,iZ2,iZ4,iZ3) = WeightsX_X2(iNodeX) * G_munu_F(iNodeX,7,iZ2,iZ4,iZ3)


      END DO

    END DO
    END DO
    END DO








    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 7*nX, nDOFX_X2, - One, LX_X2_Dn, nDOFX_X2, &
             H_munu_F(1,1,iZ_B0(2),iZ_B0(4),iZ_B0(3)  ), nDOFX_X2, Zero, &
             dG_dd_dX2_Temp, nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 7*nX, nDOFX_X2, + One, LX_X2_Up, nDOFX_X2, &
             H_munu_F(1,1,iZ_B0(2),iZ_B0(4),iZ_B0(3)+1), nDOFX_X2, One,  &
             dG_dd_dX2_Temp, nDOFX )










    ! -------------------
    ! --- Volume Term ---
    ! -------------------

    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX



        H_munu_K(iNodeX,1,iZ2,iZ4,iZ3) = WeightsX_q(iNodeX) * G_munu_K(iNodeX,1,iZ2,iZ4,iZ3)

        H_munu_K(iNodeX,2,iZ2,iZ4,iZ3) = WeightsX_q(iNodeX) * G_munu_K(iNodeX,2,iZ2,iZ4,iZ3)

        H_munu_K(iNodeX,3,iZ2,iZ4,iZ3) = WeightsX_q(iNodeX) * G_munu_K(iNodeX,3,iZ2,iZ4,iZ3)

        H_munu_K(iNodeX,4,iZ2,iZ4,iZ3) = WeightsX_q(iNodeX) * G_munu_K(iNodeX,4,iZ2,iZ4,iZ3)

        H_munu_K(iNodeX,5,iZ2,iZ4,iZ3) = WeightsX_q(iNodeX) * G_munu_K(iNodeX,5,iZ2,iZ4,iZ3)

        H_munu_K(iNodeX,6,iZ2,iZ4,iZ3) = WeightsX_q(iNodeX) * G_munu_K(iNodeX,6,iZ2,iZ4,iZ3)

        H_munu_K(iNodeX,7,iZ2,iZ4,iZ3) = WeightsX_q(iNodeX) * G_munu_K(iNodeX,7,iZ2,iZ4,iZ3)


      END DO

    END DO
    END DO
    END DO










    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 7*nX, nDOFX, - One, dLXdX2_q, nDOFX, &
             H_munu_K, nDOFX, One, dG_dd_dX2_Temp, nDOFX )




    ASSOCIATE( dZ3 => MeshX(2) % Width )

    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX








        dG_dd_dX2_Temp(iNodeX,1,iZ2,iZ4,iZ3) &
         = dG_dd_dX2_Temp(iNodeX,1,iZ2,iZ4,iZ3) &
             / ( WeightsX_q(iNodeX) * dZ3(iZ3) )

        dG_dd_dX2_Temp(iNodeX,2,iZ2,iZ4,iZ3) &
         = dG_dd_dX2_Temp(iNodeX,2,iZ2,iZ4,iZ3) &
             / ( WeightsX_q(iNodeX) * dZ3(iZ3) )

        dG_dd_dX2_Temp(iNodeX,3,iZ2,iZ4,iZ3) &
         = dG_dd_dX2_Temp(iNodeX,3,iZ2,iZ4,iZ3) &
             / ( WeightsX_q(iNodeX) * dZ3(iZ3) )

        dG_dd_dX2_Temp(iNodeX,4,iZ2,iZ4,iZ3) &
         = dG_dd_dX2_Temp(iNodeX,4,iZ2,iZ4,iZ3) &
             / ( WeightsX_q(iNodeX) * dZ3(iZ3) )


        dG_dd_dX2_Temp(iNodeX,5,iZ2,iZ4,iZ3) &
         = dG_dd_dX2_Temp(iNodeX,5,iZ2,iZ4,iZ3) &
             / ( WeightsX_q(iNodeX) * dZ3(iZ3) )

        dG_dd_dX2_Temp(iNodeX,6,iZ2,iZ4,iZ3) &
         = dG_dd_dX2_Temp(iNodeX,6,iZ2,iZ4,iZ3) &
             / ( WeightsX_q(iNodeX) * dZ3(iZ3) )

        dG_dd_dX2_Temp(iNodeX,7,iZ2,iZ4,iZ3) &
         = dG_dd_dX2_Temp(iNodeX,7,iZ2,iZ4,iZ3) &
             / ( WeightsX_q(iNodeX) * dZ3(iZ3) )




      END DO

    END DO
    END DO
    END DO


    END ASSOCIATE

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX




        dG_dd_dX2(iNodeX,0,0,iZ2,iZ3,iZ4) = dG_dd_dX2_Temp(iNodeX,1,iZ2,iZ4,iZ3) 

        dG_dd_dX2(iNodeX,0,1,iZ2,iZ3,iZ4) = dG_dd_dX2_Temp(iNodeX,2,iZ2,iZ4,iZ3) 

        dG_dd_dX2(iNodeX,0,2,iZ2,iZ3,iZ4) = dG_dd_dX2_Temp(iNodeX,3,iZ2,iZ4,iZ3) 

        dG_dd_dX2(iNodeX,0,3,iZ2,iZ3,iZ4) = dG_dd_dX2_Temp(iNodeX,4,iZ2,iZ4,iZ3) 

        dG_dd_dX2(iNodeX,1,0,iZ2,iZ3,iZ4) = dG_dd_dX2_Temp(iNodeX,2,iZ2,iZ4,iZ3) 

        dG_dd_dX2(iNodeX,2,0,iZ2,iZ3,iZ4) = dG_dd_dX2_Temp(iNodeX,3,iZ2,iZ4,iZ3) 

        dG_dd_dX2(iNodeX,3,0,iZ2,iZ3,iZ4) = dG_dd_dX2_Temp(iNodeX,4,iZ2,iZ4,iZ3) 

        dG_dd_dX2(iNodeX,1,1,iZ2,iZ3,iZ4) = dG_dd_dX2_Temp(iNodeX,5,iZ2,iZ4,iZ3) 

        dG_dd_dX2(iNodeX,2,2,iZ2,iZ3,iZ4) = dG_dd_dX2_Temp(iNodeX,6,iZ2,iZ4,iZ3) 

        dG_dd_dX2(iNodeX,3,3,iZ2,iZ3,iZ4) = dG_dd_dX2_Temp(iNodeX,7,iZ2,iZ4,iZ3) 

        dG_dd_dX2(iNodeX,1,2,iZ2,iZ3,iZ4) = 0.0_DP

        dG_dd_dX2(iNodeX,2,1,iZ2,iZ3,iZ4) = 0.0_DP 

        dG_dd_dX2(iNodeX,1,3,iZ2,iZ3,iZ4) = 0.0_DP

        dG_dd_dX2(iNodeX,3,1,iZ2,iZ3,iZ4) = 0.0_DP

        dG_dd_dX2(iNodeX,2,3,iZ2,iZ3,iZ4) = 0.0_DP

        dG_dd_dX2(iNodeX,3,2,iZ2,iZ3,iZ4) = 0.0_DP

      END DO

    END DO
    END DO
    END DO

  END SUBROUTINE ComputeWeakDerivatives_X2


  SUBROUTINE ComputeWeakDerivatives_X3 &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, dG_dd_dX3, Verbose_Option  )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)  :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)  :: &
      GX (1:nDOFX, &
          iZ_B1(2):iZ_E1(2), &
          iZ_B1(3):iZ_E1(3), &
          iZ_B1(4):iZ_E1(4), &
          1:nGF)
    REAL(DP), INTENT(out) :: &
      dG_dd_dX3 &
         (1:nDOFX,0:3,0:3, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)) 
    LOGICAL,          INTENT(in), OPTIONAL :: Verbose_Option

    INTEGER  :: iNodeX
    INTEGER  :: i, iZ2, iZ3, iZ4, iGF, iCF
    INTEGER  :: nK(4), nK_X3(4), nX, nX_X3
    REAL(DP) :: uPF_L(nPF), uPF_R(nPF), uPF_K(nPF)
    REAL(DP) :: B_u_X3(3), B_d_X3(3), A_X3, B_u_K(3), B_d_K(3), A_K
    REAL(DP) :: &
      G_munu_K(nDOFX   ,7,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3), &
           iZ_B1(4):iZ_E1(4)), &
      G_munu_F(nDOFX_X3,7,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3), &
           iZ_B0(4):iZ_E0(4)+1), &
      H_munu_K(nDOFX   ,7,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3), &
           iZ_B0(4):iZ_E0(4)), &
      H_munu_F(nDOFX_X3,7,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3), &
           iZ_B0(4):iZ_E0(4)+1)
    REAL(DP) :: & 
      dG_dd_dX3_Temp &
         (1:nDOFX,7, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))

    LOGICAL :: Verbose

    IF( iZ_E0(4) .EQ. iZ_B0(4) )THEN
      dG_dd_dX3 = Zero
      RETURN
    END IF

    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    IF (Verbose) THEN
    PRINT*, "      ComputeWeakDerivatives_X3"
    END IF

    nK    = iZ_E0 - iZ_B0 + 1 ! Number of Elements per Phase Space Dimension
    nK_X3 = nK + [0,0,0,1]    ! Number of X2 Faces per Phase Space Dimension
    nX    = PRODUCT( nK   (2:4) ) ! Number of Elements in Position Space
    nX_X3 = PRODUCT( nK_X3(2:4) ) ! Number of X2 Faces in Position Space

    ! --- Permute Geometry Fields ---


    DO iZ4 = iZ_B1(4), iZ_E1(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX


        A_K = GX(iNodeX,iZ2,iZ3,iZ4,iGF_Alpha)

        B_u_K(1) =  GX(iNodeX,iZ2,iZ3,iZ4,iGF_Beta_1)
        B_u_K(2) =  GX(iNodeX,iZ2,iZ3,iZ4,iGF_Beta_2)
        B_u_K(3) =  GX(iNodeX,iZ2,iZ3,iZ4,iGF_Beta_3)

        B_d_K(1) =  GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_11) * B_u_K(1)
        B_d_K(2) =  GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_22) * B_u_K(2)
        B_d_K(3) =  GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_33) * B_u_K(3)

        G_munu_K(iNodeX,1,iZ2,iZ3,iZ4) = -A_K**2 + B_d_K(1) * B_u_K(1) + B_d_K(2) * B_u_K(2) + B_d_K(3) * B_u_K(3)

        G_munu_K(iNodeX,2,iZ2,iZ3,iZ4) = B_d_K(1)

        G_munu_K(iNodeX,3,iZ2,iZ3,iZ4) = B_d_K(2) 
        
        G_munu_K(iNodeX,4,iZ2,iZ3,iZ4) = B_d_K(3) 
        
        G_munu_K(iNodeX,5,iZ2,iZ3,iZ4) = GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_11)

        G_munu_K(iNodeX,6,iZ2,iZ3,iZ4) = GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_22)

        G_munu_K(iNodeX,7,iZ2,iZ3,iZ4) = GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_33)


      END DO

    END DO
    END DO
    END DO


    !---------------------
    ! --- Surface Term ---
    !---------------------


    ! --- Interpolate Geometry Fields on Shared Face ---

    ! --- Face States (Average of Left and Right States) ---









      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X3, 7*nX_X3, nDOFX, One,  LX_X3_Up, nDOFX_X3, &
               G_munu_K(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4)-1), nDOFX, Zero, &
               G_munu_F(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4)), nDOFX_X3 )

      CALL MatrixMatrixMultiply &
             ( 'N', 'N', nDOFX_X3, 7*nX_X3, nDOFX, Half, LX_X3_Dn, nDOFX_X3, &
               G_munu_K(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4)), nDOFX, Half, &
               G_munu_F(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4)), nDOFX_X3 )



    ! --- Compute Metric Components from Scale Factors ---



    ! --- Permute Fluid Fields ---



    ! --- Interpolate Fluid Fields ---

    ! --- Interpolate Left State ---


    ! --- Interpolate Right State ---


    ! --- Compute Face Velocity Components ---

    DO iZ4 = iZ_B0(4), iZ_E1(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX_X3

        ! --- Left States ---




        H_munu_F(iNodeX,1,iZ2,iZ3,iZ4) = WeightsX_X3(iNodeX) * G_munu_F(iNodeX,1,iZ2,iZ3,iZ4)
        
        H_munu_F(iNodeX,2,iZ2,iZ3,iZ4) = WeightsX_X3(iNodeX) * G_munu_F(iNodeX,2,iZ2,iZ3,iZ4)
        
        H_munu_F(iNodeX,3,iZ2,iZ3,iZ4) = WeightsX_X3(iNodeX) * G_munu_F(iNodeX,3,iZ2,iZ3,iZ4)
        
        H_munu_F(iNodeX,4,iZ2,iZ3,iZ4) = WeightsX_X3(iNodeX) * G_munu_F(iNodeX,4,iZ2,iZ3,iZ4)

        H_munu_F(iNodeX,5,iZ2,iZ3,iZ4) = WeightsX_X3(iNodeX) * G_munu_F(iNodeX,5,iZ2,iZ3,iZ4)

        H_munu_F(iNodeX,6,iZ2,iZ3,iZ4) = WeightsX_X3(iNodeX) * G_munu_F(iNodeX,6,iZ2,iZ3,iZ4)

        H_munu_F(iNodeX,7,iZ2,iZ3,iZ4) = WeightsX_X3(iNodeX) * G_munu_F(iNodeX,7,iZ2,iZ3,iZ4)


      END DO

    END DO
    END DO
    END DO




    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 7*nX, nDOFX_X3, - One, LX_X3_Dn, nDOFX_X3, &
             H_munu_F(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4)  ), nDOFX_X3, Zero, &
             dG_dd_dX3_Temp, nDOFX )

    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 7*nX, nDOFX_X3, + One, LX_X3_Up, nDOFX_X3, &
             H_munu_F(1,1,iZ_B0(2),iZ_B0(3),iZ_B0(4)+1), nDOFX_X3, One,  &
             dG_dd_dX3_Temp, nDOFX )


    ! -------------------
    ! --- Volume Term ---
    ! -------------------

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX




        H_munu_K(iNodeX,1,iZ2,iZ3,iZ4) = WeightsX_q(iNodeX) * G_munu_K(iNodeX,1,iZ2,iZ3,iZ4)

        H_munu_K(iNodeX,2,iZ2,iZ3,iZ4) = WeightsX_q(iNodeX) * G_munu_K(iNodeX,2,iZ2,iZ3,iZ4)

        H_munu_K(iNodeX,3,iZ2,iZ3,iZ4) = WeightsX_q(iNodeX) * G_munu_K(iNodeX,3,iZ2,iZ3,iZ4)

        H_munu_K(iNodeX,4,iZ2,iZ3,iZ4) = WeightsX_q(iNodeX) * G_munu_K(iNodeX,4,iZ2,iZ3,iZ4)

        H_munu_K(iNodeX,5,iZ2,iZ3,iZ4) = WeightsX_q(iNodeX) * G_munu_K(iNodeX,5,iZ2,iZ3,iZ4)

        H_munu_K(iNodeX,6,iZ2,iZ3,iZ4) = WeightsX_q(iNodeX) * G_munu_K(iNodeX,6,iZ2,iZ3,iZ4)

        H_munu_K(iNodeX,7,iZ2,iZ3,iZ4) = WeightsX_q(iNodeX) * G_munu_K(iNodeX,7,iZ2,iZ3,iZ4)




      END DO

    END DO
    END DO
    END DO






    CALL MatrixMatrixMultiply &
           ( 'T', 'N', nDOFX, 7*nX, nDOFX, - One, dLXdX3_q, nDOFX, &
             H_munu_K, nDOFX, One, dG_dd_dX3_Temp, nDOFX )



    ASSOCIATE( dZ4 => MeshX(3) % Width )

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX





        dG_dd_dX3_Temp(iNodeX,1,iZ2,iZ3,iZ4) &
         = dG_dd_dX3_Temp(iNodeX,1,iZ2,iZ3,iZ4) &
             / ( WeightsX_q(iNodeX) * dZ4(iZ4) )

        dG_dd_dX3_Temp(iNodeX,2,iZ2,iZ3,iZ4) &
         = dG_dd_dX3_Temp(iNodeX,2,iZ2,iZ3,iZ4) &
             / ( WeightsX_q(iNodeX) * dZ4(iZ4) )

        dG_dd_dX3_Temp(iNodeX,3,iZ2,iZ3,iZ4) &
         = dG_dd_dX3_Temp(iNodeX,3,iZ2,iZ3,iZ4) &
             / ( WeightsX_q(iNodeX) * dZ4(iZ4) )

        dG_dd_dX3_Temp(iNodeX,4,iZ2,iZ3,iZ4) &
         = dG_dd_dX3_Temp(iNodeX,4,iZ2,iZ3,iZ4) &
             / ( WeightsX_q(iNodeX) * dZ4(iZ4) )


        dG_dd_dX3_Temp(iNodeX,5,iZ2,iZ3,iZ4) &
         = dG_dd_dX3_Temp(iNodeX,5,iZ2,iZ3,iZ4) &
             / ( WeightsX_q(iNodeX) * dZ4(iZ4) )

        dG_dd_dX3_Temp(iNodeX,6,iZ2,iZ3,iZ4) &
         = dG_dd_dX3_Temp(iNodeX,6,iZ2,iZ3,iZ4) &
             / ( WeightsX_q(iNodeX) * dZ4(iZ4) )

        dG_dd_dX3_Temp(iNodeX,7,iZ2,iZ3,iZ4) &
         = dG_dd_dX3_Temp(iNodeX,7,iZ2,iZ3,iZ4) &
             / ( WeightsX_q(iNodeX) * dZ4(iZ4) )




      END DO

    END DO
    END DO
    END DO


    END ASSOCIATE

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX

        dG_dd_dX3(iNodeX,0,0,iZ2,iZ3,iZ4) = dG_dd_dX3_Temp(iNodeX,1,iZ2,iZ3,iZ4) 

        dG_dd_dX3(iNodeX,0,1,iZ2,iZ3,iZ4) = dG_dd_dX3_Temp(iNodeX,2,iZ2,iZ3,iZ4) 

        dG_dd_dX3(iNodeX,0,2,iZ2,iZ3,iZ4) = dG_dd_dX3_Temp(iNodeX,3,iZ2,iZ3,iZ4) 

        dG_dd_dX3(iNodeX,0,3,iZ2,iZ3,iZ4) = dG_dd_dX3_Temp(iNodeX,4,iZ2,iZ3,iZ4) 

        dG_dd_dX3(iNodeX,1,0,iZ2,iZ3,iZ4) = dG_dd_dX3_Temp(iNodeX,2,iZ2,iZ3,iZ4) 

        dG_dd_dX3(iNodeX,2,0,iZ2,iZ3,iZ4) = dG_dd_dX3_Temp(iNodeX,3,iZ2,iZ3,iZ4) 

        dG_dd_dX3(iNodeX,3,0,iZ2,iZ3,iZ4) = dG_dd_dX3_Temp(iNodeX,4,iZ2,iZ3,iZ4) 

        dG_dd_dX3(iNodeX,1,1,iZ2,iZ3,iZ4) = dG_dd_dX3_Temp(iNodeX,5,iZ2,iZ3,iZ4) 

        dG_dd_dX3(iNodeX,2,2,iZ2,iZ3,iZ4) = dG_dd_dX3_Temp(iNodeX,6,iZ2,iZ3,iZ4) 

        dG_dd_dX3(iNodeX,3,3,iZ2,iZ3,iZ4) = dG_dd_dX3_Temp(iNodeX,7,iZ2,iZ3,iZ4) 

        dG_dd_dX3(iNodeX,1,2,iZ2,iZ3,iZ4) = 0.0_DP

        dG_dd_dX3(iNodeX,2,3,iZ2,iZ3,iZ4) = 0.0_DP 

        dG_dd_dX3(iNodeX,2,1,iZ2,iZ3,iZ4) = 0.0_DP

        dG_dd_dX3(iNodeX,3,2,iZ2,iZ3,iZ4) = 0.0_DP

        dG_dd_dX3(iNodeX,1,3,iZ2,iZ3,iZ4) = 0.0_DP

        dG_dd_dX3(iNodeX,3,1,iZ2,iZ3,iZ4) = 0.0_DP

      END DO

    END DO
    END DO
    END DO


  END SUBROUTINE ComputeWeakDerivatives_X3

  SUBROUTINE ComputeAlpha( V_u_1, V_u_2, V_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                           U_u, dU_d_dX0, dU_d_dX1, dU_d_dX2, dU_d_dX3, &
                           dU_u_dX0, dU_u_dX1, dU_u_dX2, dU_u_dX3, C, Alpha )

    REAL(DP), INTENT(in)  :: V_u_1, V_u_2, V_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), INTENT(in)  :: U_u(4), dU_d_dX0(4), dU_d_dX1(4), dU_d_dX2(4), dU_d_dX3(4)
    REAL(DP), INTENT(in)  :: dU_u_dX0(4), dU_u_dX1(4), dU_u_dX2(4), dU_u_dX3(4), C     
    REAL(DP), INTENT(out) :: Alpha


    REAL(DP) :: dU_d_dX(4,4), dU_u_dX(4,4), W, Vsq, Alpha_Eig, Alpha_A, A(4,4), Lambda(4), WORK(11), l_mu(4)
    INTEGER  :: mu, nu, INFO, B  
 
    Vsq =  V_u_1**2 * Gm_dd_11 + V_u_2**2 * Gm_dd_22 + V_u_3**2 * Gm_dd_33 
 
    W = 1.0_DP / SQRT( 1.0_DP - Vsq )

    dU_d_dX(1,1) = dU_d_dX0(1)
    dU_d_dX(1,2) = dU_d_dX0(2)
    dU_d_dX(1,3) = dU_d_dX0(3)
    dU_d_dX(1,4) = dU_d_dX0(4)

    dU_d_dX(2,1) = dU_d_dX1(1)
    dU_d_dX(2,2) = dU_d_dX1(2)
    dU_d_dX(2,3) = dU_d_dX1(3)
    dU_d_dX(2,4) = dU_d_dX1(4)

    dU_d_dX(3,1) = dU_d_dX2(1)
    dU_d_dX(3,2) = dU_d_dX2(2)
    dU_d_dX(3,3) = dU_d_dX2(3)
    dU_d_dX(3,4) = dU_d_dX2(4)

    dU_d_dX(4,1) = dU_d_dX3(1)
    dU_d_dX(4,2) = dU_d_dX3(2)
    dU_d_dX(4,3) = dU_d_dX3(3)
    dU_d_dX(4,4) = dU_d_dX3(4)

    dU_u_dX(1,1) = dU_u_dX0(1)
    dU_u_dX(1,2) = dU_u_dX0(2)
    dU_u_dX(1,3) = dU_u_dX0(3)
    dU_u_dX(1,4) = dU_u_dX0(4)

    dU_u_dX(2,1) = dU_u_dX1(1)
    dU_u_dX(2,2) = dU_u_dX1(2)
    dU_u_dX(2,3) = dU_u_dX1(3)
    dU_u_dX(2,4) = dU_u_dX1(4)

    dU_u_dX(3,1) = dU_u_dX2(1)
    dU_u_dX(3,2) = dU_u_dX2(2)
    dU_u_dX(3,3) = dU_u_dX2(3)
    dU_u_dX(3,4) = dU_u_dX2(4)

    dU_u_dX(4,1) = dU_u_dX3(1)
    dU_u_dX(4,2) = dU_u_dX3(2)
    dU_u_dX(4,3) = dU_u_dX3(3)
    dU_u_dX(4,4) = dU_u_dX3(4)

    B = 0

IF (B==0) THEN

    A(1,1) =  dU_d_dX(1,1) + dU_d_dX(1,1)
    A(1,2) =  dU_d_dX(1,2) + dU_d_dX(2,1)
    A(1,3) =  dU_d_dX(1,3) + dU_d_dX(3,1)
    A(1,4) =  dU_d_dX(1,4) + dU_d_dX(4,1)

    A(2,1) =  dU_d_dX(2,1) + dU_d_dX(1,2)
    A(2,2) =  dU_d_dX(2,2) + dU_d_dX(2,2)
    A(2,3) =  dU_d_dX(2,3) + dU_d_dX(3,2)
    A(2,4) =  dU_d_dX(2,4) + dU_d_dX(4,2)

    A(3,1) =  dU_d_dX(3,1) + dU_d_dX(1,3)
    A(3,2) =  dU_d_dX(3,2) + dU_d_dX(2,3)
    A(3,3) =  dU_d_dX(3,3) + dU_d_dX(3,3)
    A(3,4) =  dU_d_dX(3,4) + dU_d_dX(4,3)

    A(4,1) =  dU_d_dX(4,1) + dU_d_dX(1,4)
    A(4,2) =  dU_d_dX(4,2) + dU_d_dX(2,4)
    A(4,3) =  dU_d_dX(4,3) + dU_d_dX(3,4)
    A(4,4) =  dU_d_dX(4,4) + dU_d_dX(4,4)

    A = 0.5_DP * A
  
    CALL DSYEV( 'N', 'U', 4, A, 4, Lambda, WORK, 11, INFO )

    Alpha_Eig = MAXVAL( ABS( Lambda ) ) / W

    Alpha_A = 0.0_DP

    DO mu = 1, 4 
    DO nu = 1, 4

      Alpha_A = Alpha_A + U_u(nu) * dU_d_dX(nu,mu) * U_u(nu) * dU_u_dX(nu,mu)

    END DO
    END DO

    Alpha_A = SQRT(ABS(Alpha_A)) / W

    Alpha = C * ( Alpha_Eig + Alpha_A )

END IF


IF (B==1) THEN

    l_mu(1) = W * V_u_1 + W * V_u_2 + W * V_u_3  
    IF (Vsq .NE. 0.0_DP) THEN

      l_mu(2) = ( 1.0_DP + ( W - 1.0_DP ) * V_u_1**2 / Vsq ) &
              + ( W - 1.0_DP ) * V_u_1 * V_u_2 / Vsq &
              + ( W - 1.0_DP ) * V_u_1 * V_u_3 / Vsq 

      l_mu(3) = ( 1.0_DP + ( W - 1.0_DP ) * V_u_2**2 / Vsq ) &
              + ( W - 1.0_DP ) * V_u_1 * V_u_2 / Vsq &
              + ( W - 1.0_DP ) * V_u_2 * V_u_3 / Vsq 

      l_mu(4) = ( 1.0_DP + ( W - 1.0_DP ) * V_u_3**2 / Vsq ) &
              + ( W - 1.0_DP ) * V_u_1 * V_u_3 / Vsq &
              + ( W - 1.0_DP ) * V_u_2 * V_u_3 / Vsq 
    ELSE

      l_mu(2) = 0.0_DP
      l_mu(3) = 0.0_DP
      l_mu(4) = 0.0_DP

    END IF

    Alpha = 0.0_DP
 


    DO mu = 1, 4 
    DO nu = 1, 4

      ! Alpha = Alpha + ( l_mu(nu) * U_u(mu) + l_mu(mu) * l_mu(nu) ) * dU_dX(mu,nu)
       Alpha = Alpha + 0.5_DP * ( U_u(mu) + l_mu(mu) ) * ( dU_d_dX(mu,nu) + dU_d_dX(nu,mu) ) &
             * ( U_u(nu) + l_mu(nu))
    END DO
    END DO
 
    Alpha = C * ABS( Alpha ) / W

END IF




  END SUBROUTINE ComputeAlpha

  SUBROUTINE ComputeChristoffel( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, &
                                 dG_dd_dX0, dG_dd_dX1, dG_dd_dX2, dG_dd_dX3, Gamma_udd )
   

    INTEGER, INTENT(in) :: iZ_B0(4), iZ_E0(4)
    INTEGER, INTENT(in) :: iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)  :: &
      GX (1:nDOFX, &
          iZ_B1(2):iZ_E1(2), &
          iZ_B1(3):iZ_E1(3), &
          iZ_B1(4):iZ_E1(4), &
          1:nGF)
    REAL(DP), INTENT(in) :: & 
      dG_dd_dX0 &
         (1:nDOFX,0:3,0:3, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)), & 
      dG_dd_dX1 &
         (1:nDOFX,0:3,0:3, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)), & 
      dG_dd_dX2 &
         (1:nDOFX,0:3,0:3, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)), & 
      dG_dd_dX3 &
         (1:nDOFX,0:3,0:3, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)) 

    REAL(DP), INTENT(out) :: & 
      Gamma_udd &
         (1:nDOFX,0:3,0:3,0:3, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))



    REAL(DP) :: &
      G_uu_munu (1:nDOFX,0:3,0:3, &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4) )
     REAL(DP) :: dG_munu_dXrho &
         (1:nDOFX,0:3,0:3,0:3, &
          iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
     INTEGER :: iZ2, iZ3, iZ4, iNodeX, mu, nu, rho, sig, iNodeX1, iNodeX2
     REAL(DP) :: G_dd_11, G_dd_22, G_dd_33, B(3), A, X1, X2


    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX


        iNodeX1 = NodeNumberTableX(1,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iZ2, iNodeX1 )

        iNodeX2 = NodeNumberTableX(2,iNodeX)

        X2 = NodeCoordinate( MeshX(2), iZ3, iNodeX2 )

        dG_munu_dXrho(iNodeX,0,:,:,iZ2,iZ3,iZ4) = dG_dd_dX0(iNodeX,:,:,iZ2,iZ3,iZ4)
        dG_munu_dXrho(iNodeX,1,:,:,iZ2,iZ3,iZ4) = dG_dd_dX1(iNodeX,:,:,iZ2,iZ3,iZ4)
        dG_munu_dXrho(iNodeX,2,:,:,iZ2,iZ3,iZ4) = dG_dd_dX2(iNodeX,:,:,iZ2,iZ3,iZ4)
        dG_munu_dXrho(iNodeX,3,:,:,iZ2,iZ3,iZ4) = dG_dd_dX3(iNodeX,:,:,iZ2,iZ3,iZ4)




        G_dd_11 = GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_11)
        G_dd_22 = GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_22)
        G_dd_33 = GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_33)


 
    !    print*, iZ2,100.0_DP*ABS(2.0_DP*X1 - dG_munu_dXrho(iNodeX,1,2,2,iZ2,iZ3,iZ4)) / (2.0_DP*X1),&
     !           2.0_DP * X1, dG_munu_dXrho(iNodeX,1,2,2,iZ2,iZ3,iZ4)
       ! print*, iZ2, 100.0_DP*ABS(2.0_DP*X1*G_dd_33/G_dd_22 - dG_munu_dXrho(iNodeX,1,3,3,iZ2,iZ3,iZ4)) / (2.0*X1*G_dd_33/G_dd_22), (2.0_DP * X1*G_dd_33/G_dd_22), dG_munu_dXrho(iNodeX,1,3,3,iZ2,iZ3,iZ4)

        A = GX(iNodeX,iZ2,iZ3,iZ4,iGF_Alpha)
        B(1) = GX(iNodeX,iZ2,iZ3,iZ4,iGF_Beta_1)
        B(2) = GX(iNodeX,iZ2,iZ3,iZ4,iGF_Beta_2)
        B(3) = GX(iNodeX,iZ2,iZ3,iZ4,iGF_Beta_3)

        G_uu_munu(iNodeX,:,:,iZ2,iZ3,iZ4) = 0.0_DP
          
        G_uu_munu(iNodeX,0,0,iZ2,iZ3,iZ4) = -1.0_DP / A**2
        G_uu_munu(iNodeX,1,0,iZ2,iZ3,iZ4) = B(1) / A**2
        G_uu_munu(iNodeX,2,0,iZ2,iZ3,iZ4) = B(2) / A**2
        G_uu_munu(iNodeX,3,0,iZ2,iZ3,iZ4) = B(3) / A**2
        G_uu_munu(iNodeX,0,1,iZ2,iZ3,iZ4) = B(1) / A**2
        G_uu_munu(iNodeX,0,2,iZ2,iZ3,iZ4) = B(2) / A**2
        G_uu_munu(iNodeX,0,3,iZ2,iZ3,iZ4) = B(3) / A**2
        G_uu_munu(iNodeX,1,1,iZ2,iZ3,iZ4) = 1.0_DP / G_dd_11 - B(1) * B(1) / A**2
        G_uu_munu(iNodeX,2,2,iZ2,iZ3,iZ4) = 1.0_DP / G_dd_22 - B(2) * B(2) / A**2
        G_uu_munu(iNodeX,3,3,iZ2,iZ3,iZ4) = 1.0_DP / G_dd_33 - B(3) * B(3) / A**2
      END DO

    END DO
    END DO
    END DO


    Gamma_udd = 0.0_DP

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX


        DO rho = 0,3
        DO mu = 0,3
        DO nu = 0,3
        DO sig = 0,3

        Gamma_udd(iNodeX,rho,mu,nu,iZ2,iZ3,iZ4) = Gamma_udd(iNodeX,rho,mu,nu,iZ2,iZ3,iZ4) &
                                                + 0.5_DP * G_uu_munu(iNodeX,sig,rho,iZ2,iZ3,iZ4) &
                                                * ( dG_munu_dXrho(iNodeX,nu,mu,sig,iZ2,iZ3,iZ4) &
                                                + dG_munu_dXrho(iNodeX,mu,nu,sig,iZ2,iZ3,iZ4) &
                                                - dG_munu_dXrho(iNodeX,sig,mu,nu,iZ2,iZ3,iZ4) )


        END DO
        END DO
        END DO
        END DO

      END DO

    END DO
    END DO
    END DO

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX


       ! DO rho = 0,3
       ! DO mu = 0,3
       ! DO nu = 0,3


        !G_dd_22 = GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_22)
        !G_dd_22 = 1.0_DP/SQRT(G_dd_22)
        !print*,iZ2, 100.0_DP*ABS(G_dd_22 - Gamma_udd(iNodeX,2,1,2,iZ2,iZ3,iZ4))/G_dd_22, G_dd_22, Gamma_udd(iNodeX,2,1,2,iZ2,iZ3,iZ4)
           

        !END DO
        !END DO
        !END DO

      END DO

    END DO
    END DO
    END DO
  END SUBROUTINE ComputeChristoffel




END MODULE WeakDerivatives
