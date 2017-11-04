MODULE FluidRadiationCouplingSolutionModule_NES

  USE KindModule, ONLY: &
    DP, FourPi
  USE UnitsModule, ONLY: &
    BoltzmannConstant, &
    PlanckConstant, &
    SpeedOfLight
  USE ProgramHeaderModule, ONLY: &
    nX, nNodesX, nDOFX, &
    nE, nNodesE, nDOFE
  USE UtilitiesModule, ONLY: &
    WriteVector, &
    WriteMatrix
  USE FluidFieldsModule, ONLY: &
    uCF, nCF, &
    uPF, nPF, iPF_D, iPF_E, iPF_Ne, &
    uAF, nAF, iAF_T, iAF_E, iAF_Ye, iAF_Me
  USE EulerEquationsUtilitiesModule, ONLY: &
    ComputeConserved, &
    ComputePrimitive
  USE RadiationFieldsModule, ONLY: &
    nPR, iPR_D, iPR_I1, iPR_I2, iPR_I3
  USE MomentEquationsUtilitiesModule, ONLY: &
    ComputeConservedMoments, &
    ComputePrimitiveMoments
  USE FluidRadiationCouplingUtilitiesModule, ONLY: &
    InitializeNodes, &
    InitializeWeights, &
    InitializeFluidFields, &
    InitializeRadiationFields, &
    FinalizeFluidFields, &
    FinalizeRadiationFields, &
    ENORM
  USE EquationOfStateModule, ONLY: &
    ComputeTemperatureFromSpecificInternalEnergy, &
    ComputeThermodynamicStates_Primitive
  USE OpacityModule, ONLY: &
    ComputeScatteringOpacity_NES

  IMPLICIT NONE
  PRIVATE

  LOGICAL :: EvolveFluid
  INTEGER :: nNodesX_G, nNodesE_G
  INTEGER,  PARAMETER :: iOld = 1
  INTEGER,  PARAMETER :: iNew = 2
  REAL(DP), DIMENSION(:),     ALLOCATABLE :: E_N, W2_N, W3_N
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE :: uPF_N, uAF_N
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE :: FVEC, dUVEC
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: uPR_N
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: R0_In, R0_Out
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: UVEC, FJAC

  PUBLIC :: CoupleFluidRadiation_NES

CONTAINS


  SUBROUTINE CoupleFluidRadiation_NES &
               ( dt, iX_B0, iX_E0, iX_B1, iX_E1, U_F, dU_F, U_R, dU_R, &
                 EvolveFluid_Option )

    REAL(DP), INTENT(in)  :: &
      dt
    INTEGER,  INTENT(in)  :: &
      iX_B0(3), iX_B1(3), iX_E0(3), iX_E1(3)
    REAL(DP), INTENT(in)  :: &
      U_F (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(out) :: &
      dU_F(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)
    REAL(DP), INTENT(in)  :: &
      U_R (1:,1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:,1:)
    REAL(DP), INTENT(out) :: &
      dU_R(1:,1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:,1:)
    LOGICAL,  INTENT(in), OPTIONAL :: &
      EvolveFluid_Option

    INTEGER :: iX1, iX2, iX3

    EvolveFluid = .TRUE.
    IF( PRESENT( EvolveFluid_Option ) ) &
      EvolveFluid = EvolveFluid_Option

    CALL ComputePrimitiveMoments &
           ( iX_Begin = iX_B0, iX_End = iX_E0 )

    DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)

          CALL ComputePrimitive &
                 ( uCF(:,iX1,iX2,iX3,1:nCF), uPF(:,iX1,iX2,iX3,1:nPF) )

        END DO
      END DO
    END DO


    CALL InitializeFluidRadiationCoupling

    CALL CoupleFluidRadiation( dt )

    CALL FinalizeFluidRadiationCoupling

    DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)

          CALL ComputeTemperatureFromSpecificInternalEnergy &
                 ( uPF(1:nDOFX,iX1,iX2,iX3,iPF_D ), &
                   uAF(1:nDOFX,iX1,iX2,iX3,iAF_E ), &
                   uAF(1:nDOFX,iX1,iX2,iX3,iAF_Ye), &
                   uAF(1:nDOFX,iX1,iX2,iX3,iAF_T ) )

          CALL ComputeThermodynamicStates_Primitive &
                 ( uPF(1:nDOFX,iX1,iX2,iX3,iPF_D ), &
                   uAF(1:nDOFX,iX1,iX2,iX3,iAF_T ), &
                   uAF(1:nDOFX,iX1,iX2,iX3,iAF_Ye), &
                   uPF(1:nDOFX,iX1,iX2,iX3,iPF_E ), &
                   uAF(1:nDOFX,iX1,iX2,iX3,iAF_E ), &
                   uPF(1:nDOFX,iX1,iX2,iX3,iPF_Ne) )

          CALL ComputeConserved &
                 ( uPF(:,iX1,iX2,iX3,1:nPF), uCF(:,iX1,iX2,iX3,1:nCF) )

        END DO
      END DO
    END DO

    CALL ComputeConservedMoments &
           ( iX_Begin = iX_B0, iX_End = iX_E0 )

  END SUBROUTINE CoupleFluidRadiation_NES


  SUBROUTINE InitializeFluidRadiationCoupling

    nNodesX_G = PRODUCT(nX) * nDOFX
    nNodesE_G =         nE  * nDOFE

    ALLOCATE( E_N(nNodesE_G) )
    CALL InitializeNodes( E_N )

    ALLOCATE( W2_N(nNodesE_G), W3_N(nNodesE_G) )
    CALL InitializeWeights( W2_N, W3_N )

    ALLOCATE( uPF_N(nPF, nNodesX_G) )
    ALLOCATE( uAF_N(nAF, nNodesX_G) )
    CALL InitializeFluidFields( uPF_N, uAF_N )

    ALLOCATE( uPR_N(nNodesE_G, nPR, nNodesX_G) )
    CALL InitializeRadiationFields( uPR_N )

    ALLOCATE &
      ( R0_In (nNodesE_G,nNodesE_G,nNodesX_G), &
        R0_Out(nNodesE_G,nNodesE_G,nNodesX_G) )

    ALLOCATE(  FVEC(1:nNodesE_G,1:nNodesX_G) )
    ALLOCATE( dUVEC(1:nNodesE_G,1:nNodesX_G) )
    ALLOCATE(  UVEC(1:nNodesE_G,1:nNodesX_G,iOld:iNew) )
    ALLOCATE(  FJAC(1:nNodesE_G,1:nNodesE_G,1:nNodesX_G) )

  END SUBROUTINE InitializeFluidRadiationCoupling


  SUBROUTINE FinalizeFluidRadiationCoupling

    CALL FinalizeFluidFields( uPF_N, uAF_N )

    CALL FinalizeRadiationFields( uPR_N )

    DEALLOCATE( E_N, W2_N, W3_N )
    DEALLOCATE( uPF_N, uAF_N, uPR_N )
    DEALLOCATE( R0_In, R0_Out )
    DEALLOCATE( FVEC, dUVEC, UVEC, FJAC )

  END SUBROUTINE FinalizeFluidRadiationCoupling


  SUBROUTINE CoupleFluidRadiation( dt )

    REAL(DP), INTENT(in) :: dt

    LOGICAL, DIMENSION(1:nNodesX_G) :: Converged
    INTEGER                         :: iter, iX, iX_Max, iE
    REAL(DP)                        :: Norm, MaxNorm
    REAL(DP),             PARAMETER :: NewtonTol = 1.0d-10

    CALL SetUVEC( uPR_N(1:nNodesE_G,iPR_D,1:nNodesX_G), iOld )

    CALL SetUVEC( uPR_N(1:nNodesE_G,iPR_D,1:nNodesX_G), iNew )

    ! --- Compute Scattering Opacities ---

    CALL SetRates

    Converged = .FALSE.
    iter      = 0

    DO WHILE ( .NOT. ALL( Converged ) )

      iter = iter + 1

      ! --- Set Equation Vectors ---

      CALL SetFVEC( dt, Converged )

      ! --- Set Jacobian Matrix ---

      CALL SetFJAC( dt, Converged )

      ! --- Invert for Correction ---

      CALL SolveLinearSystems( Converged )

      MaxNorm = 0.0_DP
      DO iX = 1, nNodesX_G

        UVEC(:,iX,iNew) = UVEC(:,iX,iNew) + dUVEC(:,iX)

        Norm = ENORM( dUVEC(:,iX) / ( UVEC(:,iX,iNew) + EPSILON(1.0_DP) ) )
        IF( Norm <= NewtonTol ) Converged(iX) = .TRUE.

        IF( Norm >= MaxNorm )THEN

          MaxNorm = Norm
          iX_Max  = iX

        END IF

      END DO

      IF( iter > 9 )THEN

        WRITE(*,*)
        WRITE(*,'(A8,A)') ' ', 'CoupleFluidRadiation_NES'
        WRITE(*,*)
        WRITE(*,'(A10,A12,I6.6,A2,A11,ES10.4E2,A2,A9,I6.6)') &
          ' ', 'Iteration = ', iter, &
          ', ', '||dU/U|| = ', MaxNorm, &
          ', ', 'iX_MAX = ', iX_MAX
        WRITE(*,*)
        WRITE(*,'(A12,I4.4,A19,I4.4)') &
          '', COUNT( Converged .EQV. .TRUE. ), &
          ' converged, out of ', SIZE( Converged )
        WRITE(*,*)

        WRITE(*,*)
        WRITE(*,'(A12,A22,ES10.4E2,A3,ES10.4E2)') &
          '', 'MIN/MAX    ||dU/U|| = ', &
          MINVAL( ABS( dUVEC(:,iX_Max) &
                       / (UVEC(:,iX_Max,iNew)+EPSILON(1.0_DP)) ) ),' / ', &
          MAXVAL( ABS( dUVEC(:,iX_Max) &
                       / (UVEC(:,iX_Max,iNew)+EPSILON(1.0_DP)) ) )
        WRITE(*,'(A12,A22,ES10.4E2,A3,ES10.4E2)') &
          '', 'MIN/MAX ||FVEC(U)|| = ', &
          MINVAL( ABS( FVEC(:,iX_Max) ) ), ' / ', &
          MAXVAL( ABS( FVEC(:,iX_Max) ) )
        WRITE(*,*)

      END IF

    END DO

    CALL UpdateFluidFields( dt )

    CALL GetUVEC( uPR_N(1:nNodesE_G,iPR_D,1:nNodesX_G), iNew )

  END SUBROUTINE CoupleFluidRadiation


  SUBROUTINE SetRates

    ASSOCIATE &
      ( kT => BoltzmannConstant * uAF_N(iAF_T,1:nNodesX_G) )

    ASSOCIATE &
      ( T_N   => uAF_N(iAF_T, 1:nNodesX_G), &
        Eta_N => uAF_N(iAF_Me,1:nNodesX_G) / kT )

    CALL ComputeScatteringOpacity_NES &
           ( E_N, T_N, Eta_N, R0_In, R0_Out )

!!$    CALL WriteVector &
!!$           ( SIZE( E_N ), E_N, 'E.dat' )
!!$    CALL WriteMatrix &
!!$           ( SIZE( E_N ), SIZE( E_N ), R0_In (:,:,1), 'R_In.dat' )
!!$    CALL WriteMatrix &
!!$           ( SIZE( E_N ), SIZE( E_N ), R0_Out(:,:,1), 'R_Ot.dat' )

    END ASSOCIATE ! T_N, etc.
    END ASSOCIATE ! kT

  END SUBROUTINE SetRates


  SUBROUTINE SetUVEC( D, iState )

    REAL(DP), DIMENSION(nNodesE_G,nNodesX_G), INTENT(in) :: D
    INTEGER,                                  INTENT(in) :: iState

    INTEGER :: iX, iE

    DO iX = 1, nNodesX_G
      DO iE = 1, nNodesE_G

        UVEC(iE,iX,iState) = D(iE,iX)

      END DO
    END DO

  END SUBROUTINE SetUVEC


  SUBROUTINE GetUVEC( D, iState )

    REAL(DP), DIMENSION(nNodesE_G,nNodesX_G), INTENT(out) :: D
    INTEGER,                                  INTENT(in)  :: iState

    INTEGER :: iX, iE

    DO iX = 1, nNodesX_G
      DO iE = 1, nNodesE_G

        D(iE,iX) = UVEC(iE,iX,iState)

      END DO
    END DO

  END SUBROUTINE GetUVEC


  SUBROUTINE SetFVEC( dt, Converged )

    REAL(DP),                      INTENT(in) :: dt
    LOGICAL, DIMENSION(nNodesX_G), INTENT(in) :: Converged

    INTEGER :: iX

    DO iX = 1, nNodesX_G

      IF( Converged(iX) ) CYCLE

      FVEC(:,iX) &
        = ( UVEC(:,iX,iNew) - UVEC(:,iX,iOld) ) &
            - dt * RHSVEC( nNodesE_G, W2_N, UVEC(:,iX,iNew), &
                           R0_In(:,:,iX), R0_Out(:,:,iX) )

    END DO

  END SUBROUTINE SetFVEC


  SUBROUTINE SetFJAC( dt, Converged )

    REAL(DP),                      INTENT(in) :: dt
    LOGICAL, DIMENSION(nNodesX_G), INTENT(in) :: Converged

    INTEGER :: iX, iE

    DO iX = 1, nNodesX_G

      FJAC(:,:,iX) = 0.0_DP

      IF( Converged(iX) ) CYCLE

      ! --- Diagonal Elements ---

      DO iE = 1, nNodesE_G

        FJAC(iE,iE,iX) &
          = 1.0_DP &
            + dt * SUM( W2_N(:) &
                        * ( R0_In(:,iE,iX) * UVEC(:,iX,iNew) &
                            + R0_Out(:,iE,iX) * (FourPi-UVEC(:,iX,iNew)) ) )

      END DO

      ! --- Off-Diagonal Elements ---

      DO iE = 1, nNodesE_G

        FJAC(:,iE,iX) &
          = FJAC(:,iE,iX) &
             - dt * W2_N(iE) &
               * ( (FourPi-UVEC(:,iX,iNew)) * R0_Out(:,iE,iX) &
                   + UVEC(:,iX,iNew) * R0_Out(:,iE,iX) )

      END DO

    END DO

  END SUBROUTINE SetFJAC


  SUBROUTINE SolveLinearSystems( Converged )

    LOGICAL, DIMENSION(nNodesX_G), INTENT(in) :: Converged

    INTEGER                        :: iX, INFO
    INTEGER,  DIMENSION(nNodesE_G) :: IPIV

    DO iX = 1, nNodesX_G

      dUVEC(:,iX) = 0.0_DP

      IF( Converged(iX) ) CYCLE

      dUVEC(:,iX) = - FVEC(:,iX)

      CALL DGESV &
             ( nNodesE_G, 1, FJAC(:,:,iX), nNodesE_G, &
               IPIV, dUVEC(:,iX), nNodesE_G, INFO )

    END DO

  END SUBROUTINE SolveLinearSystems


  SUBROUTINE UpdateFluidFields( dt )

    REAL(DP), INTENT(in) :: dt

    INTEGER                        :: iX
    REAL(DP), DIMENSION(nNodesE_G) :: RHS

    ASSOCIATE( hc3 => ( PlanckConstant * SpeedOfLight )**3 )

    DO iX = 1, nNodesX_G

      RHS &
        = RHSVEC &
            ( nNodesE_G, W2_N(:), UVEC(:,iX,iNew), &
              R0_In(:,:,iX), R0_Out(:,:,iX) )

      uAF_N(iAF_E,iX) &
        = uAF_N(iAF_E,iX) &
            - dt * SUM( W3_N(:) * RHS(:) ) / hc3 / uPF_N(iPF_D,iX)

    END DO

    END ASSOCIATE ! hc3

  END SUBROUTINE UpdateFluidFields


  PURE FUNCTION RHSVEC( N, W, D, R_In, R_Out )

    INTEGER,                  INTENT(in) :: N
    REAL(DP), DIMENSION(N),   INTENT(in) :: W, D
    REAL(DP), DIMENSION(N,N), INTENT(in) :: R_In, R_Out
    REAL(DP), DIMENSION(N)               :: RHSVEC

    INTEGER :: i

    DO i = 1, N

      RHSVEC(i) &
        ! --- In-Scattering Term ---
        = ( FourPi - D(i) ) * SUM( W(:) * R_In(:,i) * D(:) ) &
        ! --- Out-Scattering Term ---
          - D(i) * SUM( W(:) * R_Out(:,i) * ( FourPi - D(:) ) )

    END DO

    RETURN
  END FUNCTION RHSVEC


END MODULE FluidRadiationCouplingSolutionModule_NES
