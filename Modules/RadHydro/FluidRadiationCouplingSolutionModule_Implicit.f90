MODULE FluidRadiationCouplingSolutionModule_Implicit

  USE KindModule, ONLY: &
    DP, Pi
  USE UnitsModule, ONLY: &
    BoltzmannConstant, &
    Centimeter, &
    MeV
  USE ProgramHeaderModule, ONLY: &
    nX, nNodesX, nDOFX, &
    nE, nNodesE, nDOFE
  USE UtilitiesModule, ONLY: &
    WriteVector, &
    WriteMatrix
  USE FluidFieldsModule, ONLY: &
    iPF_D, iPF_E, iPF_Ne, nPF, &
    iAF_T, iAF_E, iAF_Ye, iAF_Me, iAF_Mp, iAF_Mn, nAF
  USE RadiationFieldsModule, ONLY: &
    iPR_D, nPR
  USE FluidRadiationCouplingUtilitiesModule, ONLY: &
    InitializeNodes, &
    InitializeWeights, &
    InitializeFluidFields, &
    FinalizeFluidFields, &
    InitializeRadiationFields, &
    FinalizeRadiationFields
  USE EquationOfStateModule, ONLY: &
    ComputeThermodynamicStates_Auxiliary, &
    ComputeChemicalPotentials
  USE OpacityModule, ONLY: &
    ComputeAbsorptionCoefficients

  IMPLICIT NONE
  PRIVATE

  INTEGER :: nNodesX_G, nNodesE_G, iFRC_Ne, iFRC_E
  INTEGER, PARAMETER :: iOld = 0, iNew = 1
  REAL(DP), DIMENSION(:),     ALLOCATABLE :: E_N, W2_N, W3_N
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE :: uPF_N, uAF_N
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE :: Chi, f_FD
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: uPR_N

  REAL(DP), DIMENSION(:,:),   ALLOCATABLE :: FVEC_FRC, dU_FRC
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: U_FRC, FJAC_FRC

  PUBLIC :: CoupleFluidRadiation_Implicit_EmissionAbsorption

CONTAINS


  SUBROUTINE CoupleFluidRadiation_Implicit_EmissionAbsorption &
               ( dt, iX_Begin, iX_End )

    REAL(DP),              INTENT(in) :: dt
    INTEGER, DIMENSION(3), INTENT(in) :: iX_Begin, iX_End

    CALL InitializeFluidRadiationCoupling

    CALL CoupleFluidRadiation_EmissionAbsorption( dt )

    CALL FinalizeFluidRadiationCoupling

  END SUBROUTINE CoupleFluidRadiation_Implicit_EmissionAbsorption


  SUBROUTINE InitializeFluidRadiationCoupling

    nNodesX_G = PRODUCT(nX) * nDOFX
    nNodesE_G =         nE  * nDOFE

    iFRC_Ne = nNodesE_G + 1
    iFRC_E  = nNodesE_G + 2

    ALLOCATE( E_N(nNodesE_G) )
    CALL InitializeNodes( E_N )

!!$    CALL WriteVector( SIZE( E_N ), E_N / MeV, 'E.dat' )

    ALLOCATE( W2_N(nNodesE_G), W3_N(nNodesE_G) )
    CALL InitializeWeights( W2_N, W3_N )

!!$    CALL WriteVector( SIZE( W2_N ), W2_N / MeV**3, 'W2.dat' )
!!$    CALL WriteVector( SIZE( W3_N ), W3_N / MeV**4, 'W3.dat' )

    ALLOCATE( uPF_N(nPF, nNodesX_G), uAF_N(nAF, nNodesX_G) )
    CALL InitializeFluidFields( uPF_N, uAF_N )

    ALLOCATE( uPR_N(nNodesE_G, nPR, nNodesX_G) )
    CALL InitializeRadiationFields( uPR_N )

    ALLOCATE( Chi(nNodesE_G, nNodesX_G), f_FD(nNodesE_G, nNodesX_G) )

    ALLOCATE( U_FRC   (nNodesE_G+2, nNodesX_G, 0:1) )
    ALLOCATE( FVEC_FRC(nNodesE_G+2, nNodesX_G) )
    ALLOCATE( dU_FRC  (nNodesE_G+2, nNodesX_G) )
    ALLOCATE( FJAC_FRC(nNodesE_G+2, nNodesE_G+2, nNodesX_G) )

  END SUBROUTINE InitializeFluidRadiationCoupling


  SUBROUTINE FinalizeFluidRadiationCoupling

    CALL FinalizeFluidFields( uPF_N, uAF_N )

    CALL FinalizeRadiationFields( uPR_N )

    DEALLOCATE( E_N, W2_N, W3_N, uPF_N, uAF_N, uPR_N, Chi, f_FD )

    DEALLOCATE( U_FRC, FVEC_FRC, dU_FRC, FJAC_FRC )

  END SUBROUTINE FinalizeFluidRadiationCoupling


  SUBROUTINE CoupleFluidRadiation_EmissionAbsorption( dt )

    REAL(DP), INTENT(in) :: dt

!!$    CALL WriteVector( SIZE( uPR_N(:,1,1) ), uPR_N(:,1,1), 'N_0.dat' )

    CALL SetStates_FRC &
           ( uPR_N(:,iPR_D,:), uPF_N(iPF_Ne,:), uPF_N(iPF_E,:), iOld )

! --- Start iteration loop here

    CALL ComputeThermodynamicStates_Auxiliary &
           ( uPF_N(iPF_D,:), uPF_N(iPF_E,:), uPF_N(iPF_Ne,:), &
             uAF_N(iAF_T,:), uAF_N(iAF_E,:), uAF_N(iAF_Ye,:) )

    CALL SetStates_FRC &
           ( uPR_N(:,iPR_D,:), uPF_N(iPF_Ne,:), uPF_N(iPF_E,:), iNew )

    CALL SetRates_EmissionAbsorption ! --- Computes Absoption Coefficients

!!$    CALL WriteVector &
!!$           ( SIZE( Chi(:,1) ), Chi(:,1) / (1.0_DP/Centimeter), 'Chi.dat' )

    CALL SetEquilibrium_EmissionAbsorption ! --- Computes FD Distribution

!!$    CALL WriteVector( SIZE( f_FD(:,1) ), f_FD(:,1), 'f_FD.dat' )

    CALL SetFVEC_EmissionAbsorption( dt ) ! --- Equations Vector

!!$    CALL WriteVector( SIZE( FVEC_FRC(:,1) ), FVEC_FRC(:,1), 'FVEC.dat' )

    CALL SetFJAC_EmissionAbsorption( dt ) ! --- Jacobian Matrix

!!$    CALL WriteMatrix &
!!$           ( SIZE( FJAC_FRC, DIM = 1 ), SIZE( FJAC_FRC, DIM = 2 ), &
!!$             FJAC_FRC, 'FJAC.dat' )

    CALL SolveLinearSystems_FRC

!!$    CALL WriteVector( SIZE( dU_FRC(:,1) ), dU_FRC(:,1), 'dU.dat' )

    U_FRC(:,:,iNew) = U_FRC(:,:,iNew) + dU_FRC(:,:)

    CALL GetStates_FRC &
           ( uPR_N(:,iPR_D,:), uPF_N(iPF_Ne,:), uPF_N(iPF_E,:), iNew )

! --- End Iteration Loop Here

!!$    CALL WriteVector( SIZE( uPR_N(:,1,1) ), uPR_N(:,1,1), 'N_1.dat' )

  END SUBROUTINE CoupleFluidRadiation_EmissionAbsorption


  SUBROUTINE SetRates_EmissionAbsorption

    INTEGER :: iX

    DO iX = 1, nNodesX_G

      CALL ComputeAbsorptionCoefficients &
             ( E_N, [ uPF_N(iPF_D,iX) ], [ uAF_N(iAF_T,iX) ], &
               [ uAF_N(iAF_Ye,iX) ], Chi(:,iX) )

    END DO

  END SUBROUTINE SetRates_EmissionAbsorption


  SUBROUTINE SetEquilibrium_EmissionAbsorption

    INTEGER  :: iX
    REAL(DP) :: Mnu

    CALL ComputeChemicalPotentials &
           ( uPF_N(iPF_D, :), uAF_N(iAF_T, :), uAF_N(iAF_Ye,:), &
             uAF_N(iAF_Me,:), uAF_N(iAF_Mp,:), uAF_N(iAF_Mn,:) )

    DO iX = 1, nNodesX_G

      Mnu = uAF_N(iAF_Me,iX) + uAF_N(iAF_Mp,iX) - uAF_N(iAF_Mn,iX)

      f_FD(:,iX) &
        = FermiDirac( E_N, Mnu, BoltzmannConstant * uAF_N(iAF_T,iX) )

    END DO

  END SUBROUTINE SetEquilibrium_EmissionAbsorption


  SUBROUTINE SetFVEC_EmissionAbsorption( dt )

    REAL(DP), INTENT(in) :: dt

    INTEGER :: iX, iE

    DO iX = 1, nNodesX_G

      DO iE = 1, nNodesE_G

        FVEC_FRC(iE,iX) &
          = ( 1.0_DP + dt * Chi(iE,iX) ) * U_FRC(iE,iX,iNew) &
              - ( U_FRC(iE,iX,iOld) &
                    + dt * Chi(iE,iX) * 4.0_DP * Pi * f_FD(iE,iX) )

      END DO

      FVEC_FRC(iFRC_Ne,iX) &
        = ( U_FRC(iFRC_Ne,iX,iNew) - U_FRC(iFRC_Ne,iX,iOld) ) !&
!            - dt * DOT_PRODUCT &
!                   ( W2_N, Chi(:,iX) * ( 4.0_DP * Pi * f_FD(:,iX) &
!                                         - U_FRC(1:nNodesE_G,iX,1) ) )

      FVEC_FRC(iFRC_E,iX) &
        = ( U_FRC(iFRC_E,iX,iNew) - U_FRC(iFRC_E,iX,iOld) ) !&
!            - dt * DOT_PRODUCT &
!                   ( W3_N, Chi(:,iX) * ( 4.0_DP * Pi * f_FD(:,iX) &
!                                         - U_FRC(1:nNodesE_G,iX,1) ) )

    END DO

  END SUBROUTINE SetFVEC_EmissionAbsorption


  SUBROUTINE SetFJAC_EmissionAbsorption( dt )

    REAL(DP), INTENT(in) :: dt

    INTEGER :: iX, iE, jE

    FJAC_FRC = 0.0_DP

    DO iX = 1, nNodesX_G

      DO iE = 1, nNodesE_G

        FJAC_FRC(iE,iE,iX) &
          = 1.0_DP + dt * Chi(iE,iX)

      END DO

      FJAC_FRC(iFRC_Ne,iFRC_Ne,iX) = 1.0_DP

      FJAC_FRC(iFRC_E, iFRC_E, iX) = 1.0_DP

    END DO

  END SUBROUTINE SetFJAC_EmissionAbsorption


  SUBROUTINE SolveLinearSystems_FRC

    INTEGER :: &
      iX, LDA, INFO, iE
    INTEGER,  DIMENSION(nNodesE_G+2) :: &
      IPIV
    REAL(DP), DIMENSION(nNodesE_G+2) :: &
      b
    REAL(DP), DIMENSION(nNodesE_G+2,nNodesE_G+2) :: &
      A

    LDA = nNodesE_G + 2

    DO iX = 1, nNodesX_G

      b(:)   = - FVEC_FRC(:,iX)
      A(:,:) =   FJAC_FRC(:,:,iX)

      CALL DGESV( LDA, 1, A, LDA, IPIV, b, LDA, INFO )

      dU_FRC(:,iX) = b(:)

    END DO

  END SUBROUTINE SolveLinearSystems_FRC


  PURE ELEMENTAL REAL(DP) FUNCTION FermiDirac( E, Mu, kT )

    REAL(DP), INTENT(in) :: E, Mu, kT

    FermiDirac = 1.0_DP / ( EXP( ( E - Mu ) / kT ) + 1.0_DP )

    RETURN
  END FUNCTION FermiDirac


  SUBROUTINE SetStates_FRC( N, Ne, E, iState )

    REAL(DP), DIMENSION(nNodesE_G,nNodesX_G), INTENT(in) :: N
    REAL(DP), DIMENSION          (nNodesX_G), INTENT(in) :: Ne, E
    INTEGER,                                  INTENT(in) :: iState

    INTEGER :: iX

    DO iX = 1, nNodesX_G

      U_FRC(:,iX,iState) = State_FRC( N(:,iX), Ne(iX), E(iX) )

    END DO

  END SUBROUTINE SetStates_FRC


  PURE FUNCTION State_FRC( N, Ne, E )

    REAL(DP), DIMENSION(nNodesE_G+2)           :: State_FRC
    REAL(DP), DIMENSION(nNodesE_G), INTENT(in) :: N
    REAL(DP),                       INTENT(in) :: Ne, E

    State_FRC = [ N(:), Ne, E ]

    RETURN
  END FUNCTION State_FRC


  SUBROUTINE GetStates_FRC( N, Ne, E, iState )

    REAL(DP), DIMENSION(nNodesE_G,nNodesX_G), INTENT(out) :: N
    REAL(DP), DIMENSION          (nNodesX_G), INTENT(out) :: Ne, E
    INTEGER,                                  INTENT(in)  :: iState

    INTEGER :: iX, iE

    DO iX = 1, nNodesX_G
      DO iE = 1, nNodesE_G
        N(iE,iX) = U_FRC(iE,iX,iState)
      END DO
      Ne(iX) = U_FRC(iFRC_Ne,iX,iState)
      E (iX) = U_FRC(iFRC_E, iX,iState)
    END DO

  END SUBROUTINE GetStates_FRC


END MODULE FluidRadiationCouplingSolutionModule_Implicit
