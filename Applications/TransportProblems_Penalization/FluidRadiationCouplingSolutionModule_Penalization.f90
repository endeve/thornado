MODULE FluidRadiationCouplingSolutionModule_Penalization

  USE KindModule, ONLY: &
    DP, FourPi
  USE UnitsModule, ONLY: &
    Millisecond, &
    BoltzmannConstant, &
    PlanckConstant, &
    SpeedOfLight
  USE ProgramHeaderModule, ONLY: &
    nX, nNodesX, nDOFX, &
    nE, nNodesE, nDOFE, &
    nDOF
  USE FluidFieldsModule, ONLY: &
    nPF, uPF, iPF_D, &
    nAF, uAF, iAF_T, iAF_Ye, iAF_Me, &
    iAF_E, iAF_Mp, iAF_Mn
  USE RadiationFieldsModule, ONLY: &
    nPR, uPR, iPR_D
  USE FluidRadiationCouplingUtilitiesModule, ONLY: &
    InitializeNodes, &
    InitializeNodesX, &
    InitializeWeights, &
    MapForward_FluidField, &
    MapBackward_FluidField, &
    MapForward_RadiationField, &
    MapBackward_RadiationField, &
    InitializeFluidFields, &
    InitializeRadiationFields, &
    FinalizeFluidFields, &
    FinalizeRadiationFields
  USE OpacityModule, ONLY: &
    ComputeScatteringOpacity_NES
  USE EquationOfStateModule, ONLY: &
    ComputeTemperatureFromSpecificInternalEnergy
    
  IMPLICIT NONE
  PRIVATE

  INTEGER,  PUBLIC                                :: nNodesX_G
  INTEGER,  PUBLIC                                :: nNodesE_G
  REAL(DP), PUBLIC, DIMENSION(:,:),   ALLOCATABLE :: Neq
  REAL(DP), PUBLIC, DIMENSION(:,:),   ALLOCATABLE :: uPF_N, uAF_N
  REAL(DP), PUBLIC, DIMENSION(:,:,:), ALLOCATABLE :: uPR_N
  REAL(DP), PUBLIC, DIMENSION(:,:,:), ALLOCATABLE :: R0_In, R0_Out

  INTEGER,  PARAMETER :: iOld = 1
  INTEGER,  PARAMETER :: iNew = 2
  REAL(DP), DIMENSION(:),     ALLOCATABLE :: E_N, W2_N, W3_N
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE :: X_N
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: UVEC
!  REAL(DP), DIMENSION(:),     ALLOCATABLE :: D_N, T_N, Y_N
!  REAL(DP), DIMENSION(:,:),   ALLOCATABLE :: J_N, H1_N, H2_N, H3_N

  PUBLIC :: CoupleFluidRadiation
  PUBLIC :: InitializeFluidRadiationCoupling
  PUBLIC :: FinalizeFluidRadiationCoupling
  PUBLIC :: SetRates
  PUBLIC :: RHSLAMP
 
CONTAINS


  SUBROUTINE CoupleFluidRadiation( dt )

    REAL(DP), INTENT(in) :: dt
   
    CALL UpdateFluidRadiationFields_N( dt )

  END SUBROUTINE CoupleFluidRadiation


  SUBROUTINE InitializeFluidRadiationCoupling

    nNodesX_G = PRODUCT(nX) * nDOFX
    nNodesE_G =         nE  * nDOFE

    PRINT*, 'In InitializeFluidRadiationCoupling, after allocate'

    ALLOCATE( E_N(nNodesE_G) )
    CALL InitializeNodes( E_N )

    ALLOCATE( W2_N(nNodesE_G), W3_N(nNodesE_G) )
    CALL InitializeWeights( W2_N, W3_N )

    ALLOCATE( X_N(nNodesX_G,3) )
    CALL InitializeNodesX( X_N )

    ALLOCATE( Neq(nNodesE_G, nNodesX_G) )
    PRINT*, 'Neq max is: ', MAXVAL( Neq )
    Neq = 0.0

    ALLOCATE( uPF_N(nPF, nNodesX_G) )
    ALLOCATE( uAF_N(nAF, nNodesX_G) )
    CALL InitializeFluidFields( uPF_N, uAF_N )

    ALLOCATE( uPR_N(nNodesE_G, nPR, nNodesX_G) )
    CALL InitializeRadiationFields( uPR_N )

!    ALLOCATE( D_N(nNodesX_G) )
!    CALL MapForward_FluidField &
!           ( uPF(1:nDOFX,1:nX(1),1:nX(2),1:nX(3),iPF_D), &
!             D_N(1:nNodesX_G) )
!
!    ALLOCATE( T_N(nNodesX_G) )
!    CALL MapForward_FluidField &
!            ( uAF(1:nDOFX,1:nX(1),1:nX(2),1:nX(3),iAF_T), &
!              T_N(1:nNodesX_G) )
!
!    ALLOCATE( Y_N(nNodesX_G) )
!    CALL MapForward_FluidField &
!            ( uAF(1:nDOFX,1:nX(1),1:nX(2),1:nX(3),iAF_Ye), &
!              Y_N(1:nNodesX_G) )
!
!    ALLOCATE( J_N(nNodesE_G, nNodesX_G) )
!    CALL MapForward_RadiationField &
!           ( uPR(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),iPR_D,1), &
!             J_N(1:nNodesE_G,1:nNodesX_G) )
!
    ALLOCATE &
      ( R0_In (nNodesE_G,nNodesE_G,nNodesX_G), &
        R0_Out(nNodesE_G,nNodesE_G,nNodesX_G) )

    PRINT*, 'R0_In shape:', SHAPE( R0_In )
    PRINT*, 'R0_In max is: ', MAXVAL( R0_In )
    PRINT*, 'Neq shape:', SHAPE( Neq )
    PRINT*, 'Neq max is: ', MAXVAL( Neq )
    PRINT*, ' '
  
  END SUBROUTINE InitializeFluidRadiationCoupling


  SUBROUTINE FinalizeFluidRadiationCoupling

    CALL FinalizeFluidFields( uPF_N, uAF_N )

    CALL FinalizeRadiationFields( uPR_N )

    DEALLOCATE( E_N, W2_N, W3_N )
    DEALLOCATE( X_N )
!    DEALLOCATE( X_N, D_N, T_N, Y_N )
!    DEALLOCATE( J_N )
    DEALLOCATE( Neq )
    DEALLOCATE( uPF_N, uAF_N, uPR_N )
    DEALLOCATE( R0_In, R0_Out )

  END SUBROUTINE FinalizeFluidRadiationCoupling


  SUBROUTINE UpdateFluidRadiationFields_N( dt )

    REAL(DP), INTENT(in) :: dt

    INTEGER                        :: iX
    REAL(DP), DIMENSION(nNodesE_G) :: RHS

    ASSOCIATE( hc3 => ( PlanckConstant * SpeedOfLight )**3 )

    CALL SetRates

    DO iX = 1, nNodesX_G

      CALL ComputeRHS( dt, RHS, iX )

      ! --- Update J_N ---
 
!      J_N(:,iX) = J_N(:,iX) + RHS(:)

      uPR_N(:,iPR_D,iX) = uPR_N(:,iPR_D,iX) + RHS(:)
     
      ! --- Update Specific Internal Energy ---

      uAF_N(iAF_E,iX) &
        = uAF_N(iAF_E,iX) &
            - dt * SUM( W3_N(:) * RHS(:) ) / hc3 / uPF_N(iPF_D,iX)

    END DO

     ! --- Update Temperature ---

     CALL ComputeTemperatureFromSpecificInternalEnergy &
            ( uPF_N(iPF_D, 1:nNodesX_G), &
              uAF_N(iAF_E, 1:nNodesX_G), &
              uAF_N(iAF_Ye,1:nNodesX_G), &
              uAF_N(iAF_T, 1:nNodesX_G) )

    END ASSOCIATE ! hc3

  END SUBROUTINE UpdateFluidRadiationFields_N


  SUBROUTINE ComputeRHS( dt, RHS, iX )

    REAL(DP), INTENT(in) :: dt
    INTEGER,  INTENT(in) :: iX
    REAL(DP), DIMENSION(nNodesE_G), INTENT(out) :: RHS

    REAL(DP)  :: LAMP
    REAL(DP), DIMENSION(nNodesE_G,nNodesE_G) :: dC

    LAMP = RHSLAMP &
           ( nNodesE_G, R0_In(:,:,iX), R0_Out(:,:,iX), Neq(:,iX) )


    RHS = ( dt / ( 1.0_DP + LAMP * dt ) ) &
          * RHSVEC &
            ( nNodesE_G, W2_N(:), uPR_N(:,iPR_D,iX), &
              R0_In(:,:,iX), R0_Out(:,:,iX) )
  
    IF( MAXVAL( ABS( RHS ) ) > 1.0 ) THEN
       PRINT*, "ERROR in ComputeRHS: beyond boundary!"
       STOP
    END IF
 
 
  END SUBROUTINE ComputeRHS


  SUBROUTINE SetRates

    INTEGER :: i

    ASSOCIATE &
      ( kT => BoltzmannConstant * uAF_N(iAF_T,1:nNodesX_G) )

    ASSOCIATE &
      ( T_N   => uAF_N(iAF_T, 1:nNodesX_G), &
        Eta_N => uAF_N(iAF_Me,1:nNodesX_G) / kT, &
        Mnu   => uAF_N(iAF_Me,1:nNodesX_G) &
                   + uAF_N(iAF_Mp,1:nNodesX_G) &
                   - uAF_N(iAF_Mn,1:nNodesX_G) )

    CALL ComputeScatteringOpacity_NES &
           ( E_N, T_N, Eta_N, R0_In, R0_Out )

    DO i = 1, nNodesX_G

      CALL ComputeFermiDistribution &
             ( nNodesE_G, kT(i), Mnu(i), E_N, Neq(:,i) )
    
    END DO     

    END ASSOCIATE ! T_N, etc.
    END ASSOCIATE ! kT

  END SUBROUTINE SetRates


  SUBROUTINE ComputeFermiDistribution( nNodesE, kT, Mnu, E, Neq )

    INTEGER,                  INTENT(in) :: nNodesE
    REAL(DP),                 INTENT(in) :: kT, Mnu
    REAL(DP), DIMENSION(nNodesE), INTENT(in) :: E 
    REAL(DP), DIMENSION(nNodesE), INTENT(out) :: Neq

    INTEGER :: i

    DO i = 1,nNodesE

      Neq(i) = MAX( 1.0d-100, &
                    FourPi / ( EXP( (E(i)-Mnu)/kT ) + 1.0_DP ) )

    END DO

  END SUBROUTINE ComputeFermiDistribution
 
 
  PURE FUNCTION RHSVEC( nNode, W, D, R_In, R_Out )

    INTEGER,                  INTENT(in) :: nNode
    REAL(DP), DIMENSION(nNode),   INTENT(in) :: W, D
    REAL(DP), DIMENSION(nNode,nNode), INTENT(in) :: R_In, R_Out
    REAL(DP), DIMENSION(nNode)               :: RHSVEC

    INTEGER :: i

    DO i = 1, nNode

      RHSVEC(i) &
        ! --- In-Scattering Term ---
        = ( FourPi - D(i) ) * SUM( W(:) * R_In(:,i) * D(:) ) &
        ! --- Out-Scattering Term ---
          - D(i) * SUM( W(:) * R_Out(:,i) * ( FourPi - D(:) ) )

    END DO

    RETURN

  END FUNCTION RHSVEC
  

  PURE FUNCTION RHSLAMP( nNodesE, R_In_H, R_Out_H, Neq )

    INTEGER,                      INTENT(in) :: nNodesE
    REAL(DP), DIMENSION(nNodesE), INTENT(in) :: Neq
    REAL(DP), DIMENSION(nNodesE,nNodesE), INTENT(in) :: R_In_H, R_Out_H

    REAL(DP) :: RHSLAMP
    REAL(DP), DIMENSION(nNodesE,nNodesE) :: LL, MID, dC, diag_N

    INTEGER :: i,j

    diag_N = 0.0_DP

    DO i = 1, nNodesE

      diag_N(i,i) = Neq(i)
  
    END DO

    LL  = L_FUN( nNodesE, Neq, R_In_H, R_Out_H )

    MID = MatrixMultiMatrix( nNodesE, diag_N, LL )
 
    dC  = LL - MID

    RHSLAMP = MAXVAL( ABS( EIGEN( nNodesE, dC ) ) )
   
    RETURN

  END FUNCTION RHSLAMP

  
  PURE FUNCTION L_FUN( N_g, N, R_In, R_Out )
  
    INTEGER,                      INTENT(in) :: N_g
    REAL(DP), DIMENSION(N_g),     INTENT(in) :: N
    REAL(DP), DIMENSION(N_g,N_g), INTENT(in) :: R_In, R_Out
    
    REAL(DP), DIMENSION(N_g,N_g) :: L_FUN
    REAL(DP), DIMENSION(N_g,N_g) :: diag_K, Rinout

    INTEGER :: i

    Rinout = R_In - R_Out

    diag_K = 0.0_DP

    DO i = 1, N_g

        diag_K(i,i) = SUM( R_Out(i,:) ) + SUM( Rinout(i,:) * N(:) )

    END DO
   
    L_FUN = R_In - diag_K

    RETURN

  END FUNCTION L_FUN

  
  PURE FUNCTION MatrixMultiMatrix( N, MatA, MatB )
 
   INTEGER,                  INTENT(in) :: N
   REAL(DP), DIMENSION(N,N), INTENT(in) :: MatA, MatB

   REAL(DP), DIMENSION(N,N) :: MatrixMultiMatrix

   INTEGER :: i, j, k

   DO j = 1, N
     DO i = 1, N
   
     MatrixMultiMatrix(i,j) = SUM( MatA(i,:) * MatB(:,j) )

     END DO
   END DO   
  
   RETURN
 
  END FUNCTION MatrixMultiMatrix


  PURE FUNCTION EIGEN( N, Mat )

    INTEGER,                  INTENT(in) :: N
    REAL(DP), DIMENSION(N,N), INTENT(in) :: Mat

    REAL(DP), DIMENSION(N) :: EIGEN

    INTEGER :: i
 
    DO i = 1, N

      EIGEN(i) = Mat(i,i)

    END DO
    
    RETURN

  END FUNCTION EIGEN


END MODULE FluidRadiationCouplingSolutionModule_Penalization
