MODULE FluidRadiationCouplingSolutionModule_Penalization

  USE KindModule, ONLY: &
    DP, FourPi
  USE UnitsModule, ONLY: &
    Millisecond, &
    BoltzmannConstant, &
    PlanckConstant, &
    SpeedOfLight, &
    Kilometer, &
    Gram, &
    Centimeter, &
    MeV
  USE ProgramHeaderModule, ONLY: &
    nX, nNodesX, nDOFX, &
    nE, nNodesE, nDOFE, &
    nDOF
  USE UtilitiesModule, ONLY: &
    WriteVector, &
    WriteMatrix
  USE FluidFieldsModule, ONLY: &
    nPF, uPF, iPF_D, &
    nAF, uAF, iAF_T, iAF_Ye, iAF_Me, iAF_Mp, iAF_Mn
  USE RadiationFieldsModule, ONLY: &
    nCR, iCR_N, nSpecies, &
    nPR, uPR, iPR_D, iPR_I1, iPR_I2, iPR_I3
  USE MomentEquationsUtilitiesModule, ONLY: &
    ComputePrimitiveMoments
  USE FluidRadiationCouplingUtilitiesModule, ONLY: &
    InitializeNodes, &
    InitializeNodesX, &
    InitializeWeights, &
    MapForward_FluidField, &
    MapBackward_FluidField, &
    MapForward_RadiationField, &
    MapBackward_RadiationField, &
    FermiDirac
  USE OpacityModule, ONLY: &
    ComputeScatteringOpacity_NES
  USE MeshModule, ONLY: &
    MeshE 
   
  IMPLICIT NONE
  PRIVATE

  REAL(DP), DIMENSION(:,:,:,:),       ALLOCATABLE, PUBLIC :: absLambda
  REAL(DP), DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE, PUBLIC :: rhsCR_C

  INTEGER                                 :: nNodesX_G
  INTEGER                                 :: nNodesE_G
  REAL(DP), DIMENSION(:),     ALLOCATABLE :: E_N, W2_N, W3_N
  REAL(DP), DIMENSION(:),     ALLOCATABLE :: D_N, T_N, Y_N
  REAL(DP), DIMENSION(:),     ALLOCATABLE :: Me_N, Mnu_N
  REAL(DP), DIMENSION(:),     ALLOCATABLE :: absLambda_N
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE :: X_N
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE :: FD
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE :: J_N, H1_N, H2_N, H3_N
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE :: RHS_J
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: R0_In, R0_Out

  PUBLIC :: InitializeFluidRadiationCoupling
  PUBLIC :: FinalizeFluidRadiationCoupling
  PUBLIC :: ComputeRHS_C
 
CONTAINS


  SUBROUTINE InitializeFluidRadiationCoupling

    ALLOCATE( absLambda(nDOFX,nX(1),nX(2),nX(3)) )
    ALLOCATE( rhsCR_C(nDOF,nE,nX(1),nX(2),nX(3),nCR,nSpecies) )

    nNodesX_G = PRODUCT(nX) * nDOFX
    nNodesE_G =         nE  * nDOFE

    ALLOCATE( E_N(nNodesE_G) )
    CALL InitializeNodes( E_N )
 
    ALLOCATE( W2_N(nNodesE_G), W3_N(nNodesE_G) )
    CALL InitializeWeights( W2_N, W3_N )

    ALLOCATE( X_N(nNodesX_G,3) )
    CALL InitializeNodesX( X_N )

    ALLOCATE( D_N        (nNodesX_G) )
    ALLOCATE( T_N        (nNodesX_G) )
    ALLOCATE( Y_N        (nNodesX_G) )
    ALLOCATE( Me_N       (nNodesX_G) )
    ALLOCATE( Mnu_N      (nNodesX_G) )
    ALLOCATE( absLambda_N(nNodesX_G) )

    ALLOCATE( FD   (nNodesE_G, nNodesX_G) )
    ALLOCATE( J_N  (nNodesE_G, nNodesX_G) )
    ALLOCATE( H1_N (nNodesE_G, nNodesX_G) )
    ALLOCATE( H2_N (nNodesE_G, nNodesX_G) )
    ALLOCATE( H3_N (nNodesE_G, nNodesX_G) )
    ALLOCATE( RHS_J(nNodesE_G, nNodesX_G) )

    ALLOCATE &
      ( R0_In (nNodesE_G,nNodesE_G,nNodesX_G), &
        R0_Out(nNodesE_G,nNodesE_G,nNodesX_G) )

  END SUBROUTINE InitializeFluidRadiationCoupling


  SUBROUTINE FinalizeFluidRadiationCoupling

    DEALLOCATE( absLambda )
    DEALLOCATE( rhsCR_C )

    DEALLOCATE( E_N, W2_N, W3_N )
    DEALLOCATE( X_N )
    DEALLOCATE( D_N, T_N, Y_N )
    DEALLOCATE( Me_N, Mnu_N )
    DEALLOCATE( absLambda_N )
    DEALLOCATE( FD )
    DEALLOCATE( J_N, H1_N, H2_N, H3_N )
    DEALLOCATE( RHS_J )
    DEALLOCATE( R0_In, R0_Out )

  END SUBROUTINE FinalizeFluidRadiationCoupling


  SUBROUTINE ComputeRHS_C

    CALL ComputePrimitiveMoments &
           ( iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ] )

    CALL MapForward_FluidField &
           ( uPF(1:nDOFX,1:nX(1),1:nX(2),1:nX(3),iPF_D),  D_N(1:nNodesX_G) )
    CALL MapForward_FluidField &
           ( uAF(1:nDOFX,1:nX(1),1:nX(2),1:nX(3),iAF_T),  T_N(1:nNodesX_G) )
    CALL MapForward_FluidField &
           ( uAF(1:nDOFX,1:nX(1),1:nX(2),1:nX(3),iAF_Ye), Y_N(1:nNodesX_G) )
    CALL MapForward_FluidField &
           ( uAF(1:nDOFX,1:nX(1),1:nX(2),1:nX(3),iAF_Me), Me_N(1:nNodesX_G) )
    CALL MapForward_FluidField &
           ( uAF  (1:nDOFX,1:nX(1),1:nX(2),1:nX(3),iAF_Me) &
             + uAF(1:nDOFX,1:nX(1),1:nX(2),1:nX(3),iAF_Mp) &
             - uAF(1:nDOFX,1:nX(1),1:nX(2),1:nX(3),iAF_Mn), &
             Mnu_N(1:nNodesX_G) )

    CALL MapForward_RadiationField &
           ( uPR(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),iPR_D, 1), &
             J_N(1:nNodesE_G,1:nNodesX_G) )
    CALL MapForward_RadiationField &
           ( uPR(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),iPR_I1,1), &
             H1_N(1:nNodesE_G,1:nNodesX_G) )
    CALL MapForward_RadiationField &
           ( uPR(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),iPR_I2,1), &
             H2_N(1:nNodesE_G,1:nNodesX_G) )
    CALL MapForward_RadiationField &
           ( uPR(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),iPR_I3,1), &
             H3_N(1:nNodesE_G,1:nNodesX_G) )

    CALL SetRates

    CALL SetEquilibrium

    CALL SetAbsLambda

    CALL SetRHS

    CALL MapBackward_FluidField &
           ( absLambda(1:nDOFX,1:nX(1),1:nX(2),1:nX(3)), &
             absLambda_N(1:nNodesX_G) ) 

    CALL MapBackward_RadiationField &
           ( rhsCR_C(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),iCR_N,1), &
             RHS_J(1:nNodesE_G,1:nNodesX_G) )

  END SUBROUTINE ComputeRHS_C


  SUBROUTINE SetRates

    ASSOCIATE &
      ( Eta_N => Me_N / ( BoltzmannConstant * T_N ) )

    CALL ComputeScatteringOpacity_NES &
           ( E_N, T_N, Eta_N, R0_In, R0_Out )

    END ASSOCIATE

  END SUBROUTINE SetRates


  SUBROUTINE SetEquilibrium

    INTEGER :: iX

    ASSOCIATE( kT => BoltzmannConstant * T_N )

    DO iX = 1, nNodesX_G

      FD(:,iX) = FermiDirac( E_N, Mnu_N(iX), kT(iX) )

    END DO

    END ASSOCIATE

  END SUBROUTINE SetEquilibrium


  SUBROUTINE SetAbsLambda

    INTEGER :: iX, i

    INTEGER :: INFO, WORKL
    REAL(DP), DIMENSION(:), ALLOCATABLE:: WORK
    REAL(DP), DIMENSION(1,nNodesE_G) :: doummy
  
    REAL(DP), DIMENSION(nNodesE_G) :: dV, EigendC, EigendC_R, EigendC_I
    REAL(DP), DIMENSION(nNodesE_G,nNodesE_G) :: diag_dV, dC, diag_FD, &
                                                R_In_H, R_Out_H

    diag_dV = 0.0_DP
    DO i = 1, nNodesE_G
      diag_dV(i,i) = W2_N(i)
    END DO

    DO iX = 1, nNodesX_G

      diag_FD = 0.0_DP
      DO i = 1, nNodesE_G
        diag_FD(i,i) = FD(i,iX)
      END DO

      R_In_H  = MATMUL( R0_In(:,:,iX),  diag_dV )
      R_Out_H = MATMUL( R0_Out(:,:,iX), diag_dV )   

      dC = L_FUN( FD(:,iX), R_In_H, R_Out_H ) &
             -  MATMUL( diag_FD, R_In_H - R_Out_H )

      WORKL = 4 * nNodesE_G
      ALLOCATE( WORK( WORKL ) )

      CALL dgeev("N","N",nNodesE_G,dC,nNodesE_G,EigendC_R,EigendC_I,&
                 doummy, 1, doummy, 1,WORK, WORKL,INFO)
      DEALLOCATE( WORK )

      EigendC = EigendC_R

      absLambda_N(iX) = MAXVAL( ABS( EigendC ) )

    END DO
 
  END SUBROUTINE SetAbsLambda


  SUBROUTINE SetRHS

    INTEGER :: iX, iE

    DO iX = 1, nNodesX_G
      DO iE = 1, nNodesE_G

        RHS_J(iE,iX) &
          = (FourPi-J_N(iE,iX))*SUM( W2_N(:)*R0_In(:,iE,iX)*J_N(:,iX) ) &
            - J_N(iE,iX)*SUM( W2_N(:)*R0_Out(:,iE,iX)*(FourPi-J_N(:,iX)) )

      END DO
    END DO

  END SUBROUTINE SetRHS


  PURE FUNCTION L_FUN( N, R_In_H, R_Out_H )

    REAL(DP), DIMENSION(nNodesE_G),       INTENT(in) :: N
    REAL(DP), DIMENSION(nNodesE_G,nNodesE_G), INTENT(in) :: R_In_H, R_Out_H

    REAL(DP), DIMENSION(nNodesE_G,nNodesE_G) :: L_FUN
    REAL(DP), DIMENSION(nNodesE_G,nNodesE_G) :: diag_K, Rinout_H

    INTEGER :: i

    Rinout_H = R_In_H - R_Out_H

    diag_K = 0.0_DP
    DO i = 1, nNodesE_G
        diag_K(i,i) = SUM( R_Out_H(i,:) ) + SUM( Rinout_H(i,:) * N(:) )
    END DO

    L_FUN = R_In_H - diag_K

    RETURN

  END FUNCTION L_FUN


END MODULE FluidRadiationCouplingSolutionModule_Penalization
