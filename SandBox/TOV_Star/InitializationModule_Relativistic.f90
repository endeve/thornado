MODULE InitializationModule_Relativistic

  USE KindModule, ONLY: &
    DP, &
    SqrtTiny, &
    Zero, &
    One, &
    Two, &
    Pi, &
    FourPi
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    nNodesX, &
    nDOFX, &
    swX, &
    iX_B0, &
    iX_B1, &
    iX_E0, &
    iX_E1
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE GravitySolutionModule_XCFC_Poseidon, ONLY: &
    SolveGravity_XCFC_Poseidon
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX_FromScaleFactors
  USE GeometryFieldsModule, ONLY: &
    nGF, &
    uGF, &
    iGF_h_1, &
    iGF_h_2, &
    iGF_h_3, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_Beta_1, &
    iGF_Beta_2, &
    iGF_Beta_3, &
    iGF_Alpha, &
    iGF_Psi
  USE FluidFieldsModule, ONLY: &
    uPF, &
    iPF_D, &
    iPF_V1, &
    iPF_V2, &
    iPF_V3, &
    iPF_E, &
    iPF_Ne, &
    uCF, &
    iCF_D, &
    iCF_S1, &
    iCF_S2, &
    iCF_S3, &
    iCF_E, &
    iCF_Ne, &
    uAF, &
    iAF_P
  USE EquationOfStateModule_IDEAL, ONLY: &
    Gamma_IDEAL, &
    ComputePressureFromPrimitive_IDEAL
  USE Euler_UtilitiesModule_Relativistic, ONLY: &
    ComputeConserved_Euler_Relativistic
  USE UnitsModule, ONLY: &
    GravitationalConstant, &
    SpeedOfLight, &
    Kilometer, &
    SolarMass, &
    Gram, &
    Centimeter, &
    Erg, &
    Second, &
    PlanckConstant, &
    AtomicMassUnit
  USE UtilitiesModule, ONLY: &
    NodeNumberX, &
    Locate, &
    Interpolate1D_Linear
  USE Poseidon_UtilitiesModule, ONLY: &
    ComputeSourceTerms_Poseidon

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFields_Relativistic

  REAL(DP), PARAMETER :: &
    PolytropicConstant_TOV &
!!$      = 1.0_DP / 20.0_DP * ( 3.0_DP / Pi )**( 2.0_DP / 3.0_DP ) &
!!$          * ( PlanckConstant / ( Erg * Second ) )**2 &
!!$          / ( AtomicMassUnit / Gram )**( 8.0_DP / 3.0_DP ) &
!!$          * ( Erg / Centimeter**3 ) &
!!$          / ( Gram / Centimeter**3 )**( 5.0_DP/ 3.0_DP )
      = 1.455e5_DP * Erg * Centimeter**3 / Gram**2


CONTAINS


  SUBROUTINE InitializeFields_Relativistic

    uPF(:,:,:,:,iPF_Ne) = Zero

    WRITE(*,*)
    WRITE(*,'(A,A)') '    INFO: ', TRIM( ProgramName )

    SELECT CASE ( TRIM( ProgramName ) )

      CASE( 'StaticTOV' )

         CALL InitializeFields_StaticTOV

      CASE( 'DynamicTOV' )

         CALL InitializeFields_DynamicTOV

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,'(A21,A)') 'Invalid ProgramName: ', ProgramName
        WRITE(*,'(A)')     'Stopping...'
        STOP

    END SELECT

  END SUBROUTINE InitializeFields_Relativistic


  SUBROUTINE InitializeFields_StaticTOV

    REAL(DP), PARAMETER :: &
!!$      CentralDensity = 3.301e14_DP * ( Gram / Centimeter**3 ), &
      CentralDensity = 7.906e14_DP * ( Gram / Centimeter**3 ), &
      dX1            = 1.0e-4_DP * Kilometer, &
      TolF           = 1.0e-15_DP

    INTEGER  :: iX1, iX2, iX3, iNodeX, iNodeX1, iNodeX2, iNodeX3, &
                jNodeX, jNodeX1, iL, ITER, nX, iGF
    REAL(DP) :: X1, X2
    REAL(DP) :: Pressure , E1 , E2, Psi, Alpha, Phi
    REAL(DP) :: PressureN, E1N, E2N
    REAL(DP) :: CentralPressure, Psi0, Alpha0
    REAL(DP) :: GravitationalMass, Radius, dAlpha, dPsi, Alpha_A, Psi_A, dF

    REAL(DP), ALLOCATABLE :: PressureArr(:), DensityArr(:), &
                             AlphaArr(:), PsiArr(:), X1Arr(:)

    INTEGER, PARAMETER :: nMaxIter = 1000

    REAL(DP), ALLOCATABLE :: SourceTerms_Poseidon(:,:,:,:,:)
    LOGICAL               :: CONVERGED

    CentralPressure = PolytropicConstant_TOV * CentralDensity**( Gamma_IDEAL )
    Psi0            = 1.4_DP ! --- Initial guess ---
    Alpha0          = 0.8_DP ! --- Initial guess ---

    WRITE(*,*)
    WRITE(*,'(6x,A,ES10.3E3,A)' ) &
      'Polytropic Constant = ', PolytropicConstant_TOV &
                                  / ( Erg / Centimeter**3 &
                                  / ( Gram / Centimeter**3 )** &
                                    ( Gamma_IDEAL &
                                    ) ), &
      ' [ erg / cm^3 / ( g / cm^3 )^( Gamma ) ]'
    WRITE(*,'(6x,A,ES10.3E3,A)')  &
      'Central Density     = ', CentralDensity &
                                  / ( Gram / Centimeter**3 ), &
      ' [ g / cm^3 ]'
    WRITE(*,'(6x,A,ES10.3E3,A)')  &
      'Central Pressure    = ', CentralPressure &
                                  / ( Erg / Centimeter**3 ), &
      ' [ erg / cm^3 ]'
    WRITE(*,'(6x,A,ES10.3E3,A)')  &
      'dX                  = ', dX1 &
                                  / ( Kilometer ), &
      ' [ km ]'
    WRITE(*,'(6x,A,ES10.3E3)')    &
      'TolF                = ', TolF
    WRITE(*,*)

    ! --- Find geometry fields at center by iteratively integrating outward ---

    dF   = 1.1_DP * TolF
    ITER = 0
    DO WHILE( dF .GT. TolF .AND. ITER .LT. nMaxIter )

      ITER = ITER + 1

      IF( MOD( ITER, 100 ) .EQ. 0 ) PRINT*, 'Iteration ', ITER

      CALL IntegrateOutwards &
             ( dX1, CentralPressure, Psi0, Alpha0, &
               GravitationalMass, Radius, dAlpha, dPsi, Alpha_A, Psi_A, nX )

      ! --- Update guess for central values ---

      Alpha0 = Alpha0 + dAlpha
      Psi0   = Psi0   + dPsi

      dF = MAX( ABS( dAlpha / Alpha_A ), ABS( dPsi / Psi_A ) )

    END DO

    WRITE(*,'(6x,A,I4.4)') &
      'nIterations         = ', ITER
    WRITE(*,'(6x,A,ES13.6E3)'  )  &
      'dF                  = ', dF
    WRITE(*,'(6x,A,ES13.6E3)'  )  &
      'Alpha0              = ', Alpha0
    WRITE(*,'(6x,A,ES13.6E3)'  )  &
      'Psi0                = ', Psi0
    WRITE(*,'(6x,A,ES13.6E3,A)')  &
      'Radius              = ', Radius / Kilometer, ' km'
    WRITE(*,'(6x,A,ES13.6E3,A)')  &
      'Gravitational Mass  = ', GravitationalMass / SolarMass, ' Msun'

    ! --- Populate arrays ---

    ALLOCATE( PressureArr(nX), DensityArr(nX), AlphaArr(nX), &
              PsiArr(nX), X1Arr(nX) )

    ! --- Central values ---

    X1                = SqrtTiny * Kilometer
    Pressure          = CentralPressure
    E1                = Zero
    E2                = Zero
    Psi               = Psi0
    Phi               = Alpha0 * Psi0
    GravitationalMass = E1 / SpeedOfLight**2

    DO iX1 = 1, nX

      X1Arr      (iX1) = X1
      PressureArr(iX1) = Pressure
      DensityArr (iX1) = ( Pressure &
                             / PolytropicConstant_TOV )**( One / Gamma_IDEAL )

      ! --- Explicit steps ---

      PressureN = Pressure + dX1 * dpdr  ( Pressure, Phi, Psi, E1, E2, X1 )
      E1N       = E1       + dX1 * dE1dr ( Pressure, Psi, X1 )
      E2N       = E2       + dX1 * dE2dr ( Pressure, Phi, Psi, X1 )

      ! --- Implicit steps ---

      X1  = X1  + dX1
      Psi = Psi + dX1 * dPsidr( PressureN, E1N, X1 )
      Phi = Phi + dX1 * dPhidr( PressureN, E2N, X1 )

      PsiArr  (iX1) = Psi
      AlphaArr(iX1) = Phi / Psi

      Pressure = PressureN
      E1       = E1N
      E2       = E2N

    END DO

    ! --- Map to 3D domain ---

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E1(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)
        iNodeX2 = NodeNumberTableX(2,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
        X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

        iL = Locate( X1, X1Arr, nX )

        ! --- Geometry Fields ---

        uGF(iNodeX,iX1,iX2,iX3,iGF_Alpha) &
          = Interpolate1D_Linear( X1, X1Arr(iL), X1Arr(iL+1), &
                                  AlphaArr(iL), AlphaArr(iL+1) )

        uGF(iNodeX,iX1,iX2,iX3,iGF_Psi) &
          = Interpolate1D_Linear( X1, X1Arr(iL), X1Arr(iL+1), &
                                  PsiArr(iL), PsiArr(iL+1) )

        uGF(iNodeX,iX1,iX2,iX3,iGF_h_1) &
          = uGF(iNodeX,iX1,iX2,iX3,iGF_Psi)**2
        uGF(iNodeX,iX1,iX2,iX3,iGF_h_2) &
          = uGF(iNodeX,iX1,iX2,iX3,iGF_Psi)**2 * X1
        uGF(iNodeX,iX1,iX2,iX3,iGF_h_3) &
          = uGF(iNodeX,iX1,iX2,iX3,iGF_Psi)**2 * X1 * SIN( X2 )

        CALL ComputeGeometryX_FromScaleFactors( uGF(:,iX1,iX2,iX3,:) )

        uGF(iNodeX,iX1,iX2,iX3,iGF_Beta_1) = Zero
        uGF(iNodeX,iX1,iX2,iX3,iGF_Beta_2) = Zero
        uGF(iNodeX,iX1,iX2,iX3,iGF_Beta_3) = Zero

        ! --- Fluid Fields ---

        uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
          = Interpolate1D_Linear( X1, X1Arr(iL), X1Arr(iL+1), &
                                  DensityArr(iL), DensityArr(iL+1) )

        uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = Zero
        uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = Zero
        uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = Zero

        uPF(iNodeX,iX1,iX2,iX3,iPF_E) &
          = Interpolate1D_Linear &
              ( X1, X1Arr(iL), X1Arr(iL+1), &
                PressureArr(iL), PressureArr(iL+1) ) / ( Gamma_IDEAL - One )

      END DO

      CALL ComputePressureFromPrimitive_IDEAL &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_E), &
               uPF(:,iX1,iX2,iX3,iPF_Ne), uAF(:,iX1,iX2,iX3,iAF_P) )

      CALL ComputeConserved_Euler_Relativistic &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_V1), &
               uPF(:,iX1,iX2,iX3,iPF_V2), uPF(:,iX1,iX2,iX3,iPF_V3), &
               uPF(:,iX1,iX2,iX3,iPF_E ), uPF(:,iX1,iX2,iX3,iPF_Ne), &
               uCF(:,iX1,iX2,iX3,iCF_D ), uCF(:,iX1,iX2,iX3,iCF_S1), &
               uCF(:,iX1,iX2,iX3,iCF_S2), uCF(:,iX1,iX2,iX3,iCF_S3), &
               uCF(:,iX1,iX2,iX3,iCF_E ), uCF(:,iX1,iX2,iX3,iCF_Ne), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33), &
               uAF(:,iX1,iX2,iX3,iAF_P) )

    END DO
    END DO
    END DO

    ! --- Apply reflecting boundary conditions to geometry fields (X1) ---

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = 1, swX(1)

        DO iNodeX3 = 1, nNodesX(3)
        DO iNodeX2 = 1, nNodesX(2)
        DO iNodeX1 = 1, nNodesX(1)

          jNodeX1 = ( nNodesX(1) - iNodeX1 ) + 1

          iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )
          jNodeX = NodeNumberX( jNodeX1, iNodeX2, iNodeX3 )

          DO iGF = 1, nGF

            uGF(iNodeX,iX_B0(1)-iX1,iX2,iX3,iGF) &
              = uGF(jNodeX,iX_B0(1),iX2,iX3,iGF)

          END DO

        END DO
        END DO
        END DO

    END DO
    END DO
    END DO

    DEALLOCATE( X1Arr, PsiArr, AlphaArr, DensityArr, PressureArr )

  END SUBROUTINE InitializeFields_StaticTOV


  SUBROUTINE InitializeFields_DynamicTOV

    REAL(DP), PARAMETER :: &
!!$      CentralDensity = 3.301e14_DP * ( Gram / Centimeter**3 ), &
      CentralDensity = 7.906e14_DP * ( Gram / Centimeter**3 ), &
      dX1            = 1.0e-4_DP * Kilometer, &
      TolF           = 1.0e-15_DP

    INTEGER  :: iX1, iX2, iX3, iNodeX, iNodeX1, iNodeX2, iNodeX3, &
                jNodeX, jNodeX1, iL, ITER, nX, iGF
    REAL(DP) :: X1, X2
    REAL(DP) :: Pressure , E1 , E2, Psi, Alpha, Phi
    REAL(DP) :: PressureN, E1N, E2N
    REAL(DP) :: CentralPressure, Psi0, Alpha0
    REAL(DP) :: GravitationalMass, Radius, dAlpha, dPsi, Alpha_A, Psi_A, dF

    REAL(DP), ALLOCATABLE :: PressureArr(:), DensityArr(:), &
                             AlphaArr(:), PsiArr(:), X1Arr(:)

    INTEGER, PARAMETER :: nMaxIter = 1000

    REAL(DP), ALLOCATABLE :: SourceTerms_Poseidon(:,:,:,:,:)
    LOGICAL               :: CONVERGED

    CentralPressure = PolytropicConstant_TOV * CentralDensity**( Gamma_IDEAL )
    Psi0            = 1.4_DP ! --- Initial guess ---
    Alpha0          = 0.8_DP ! --- Initial guess ---

    WRITE(*,*)
    WRITE(*,'(6x,A,ES10.3E3,A)' ) &
      'Polytropic Constant = ', PolytropicConstant_TOV &
                                  / ( Erg / Centimeter**3 &
                                  / ( Gram / Centimeter**3 )** &
                                    ( Gamma_IDEAL &
                                    ) ), &
      ' [ erg / cm^3 / ( g / cm^3 )^( Gamma ) ]'
    WRITE(*,'(6x,A,ES10.3E3,A)')  &
      'Central Density     = ', CentralDensity &
                                  / ( Gram / Centimeter**3 ), &
      ' [ g / cm^3 ]'
    WRITE(*,'(6x,A,ES10.3E3,A)')  &
      'Central Pressure    = ', CentralPressure &
                                  / ( Erg / Centimeter**3 ), &
      ' [ erg / cm^3 ]'
    WRITE(*,'(6x,A,ES10.3E3,A)')  &
      'dX                  = ', dX1 &
                                  / ( Kilometer ), &
      ' [ km ]'
    WRITE(*,'(6x,A,ES10.3E3)')    &
      'TolF                = ', TolF
    WRITE(*,*)

    ! --- Find geometry fields at center by iteratively integrating outward ---

    dF   = 1.1_DP * TolF
    ITER = 0
    DO WHILE( dF .GT. TolF .AND. ITER .LT. nMaxIter )

      ITER = ITER + 1

      IF( MOD( ITER, 100 ) .EQ. 0 ) PRINT*, 'Iteration ', ITER

      CALL IntegrateOutwards &
             ( dX1, CentralPressure, Psi0, Alpha0, &
               GravitationalMass, Radius, dAlpha, dPsi, Alpha_A, Psi_A, nX )

      ! --- Update guess for central values ---

      Alpha0 = Alpha0 + dAlpha
      Psi0   = Psi0   + dPsi

      dF = MAX( ABS( dAlpha / Alpha_A ), ABS( dPsi / Psi_A ) )

    END DO

    WRITE(*,'(6x,A,I4.4)') &
      'nIterations         = ', ITER
    WRITE(*,'(6x,A,ES13.6E3)'  )  &
      'dF                  = ', dF
    WRITE(*,'(6x,A,ES13.6E3)'  )  &
      'Alpha0              = ', Alpha0
    WRITE(*,'(6x,A,ES13.6E3)'  )  &
      'Psi0                = ', Psi0
    WRITE(*,'(6x,A,ES13.6E3,A)')  &
      'Radius              = ', Radius / Kilometer, ' km'
    WRITE(*,'(6x,A,ES13.6E3,A)')  &
      'Gravitational Mass  = ', GravitationalMass / SolarMass, ' Msun'

    ! --- Populate arrays ---

    ALLOCATE( PressureArr(nX), DensityArr(nX), AlphaArr(nX), &
              PsiArr(nX), X1Arr(nX) )

    ! --- Central values ---

    X1                = SqrtTiny * Kilometer
    Pressure          = CentralPressure
    E1                = Zero
    E2                = Zero
    Psi               = Psi0
    Phi               = Alpha0 * Psi0
    GravitationalMass = E1 / SpeedOfLight**2

    DO iX1 = 1, nX

      X1Arr      (iX1) = X1
      PressureArr(iX1) = Pressure
      DensityArr (iX1) = ( Pressure &
                             / PolytropicConstant_TOV )**( One / Gamma_IDEAL )

      ! --- Explicit steps ---

      PressureN = Pressure + dX1 * dpdr  ( Pressure, Phi, Psi, E1, E2, X1 )
      E1N       = E1       + dX1 * dE1dr ( Pressure, Psi, X1 )
      E2N       = E2       + dX1 * dE2dr ( Pressure, Phi, Psi, X1 )

      ! --- Implicit steps ---

      X1  = X1  + dX1
      Psi = Psi + dX1 * dPsidr( PressureN, E1N, X1 )
      Phi = Phi + dX1 * dPhidr( PressureN, E2N, X1 )

      PsiArr  (iX1) = Psi
      AlphaArr(iX1) = Phi / Psi

      Pressure = PressureN
      E1       = E1N
      E2       = E2N

    END DO

    ! --- Map to 3D domain ---

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E1(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)
        iNodeX2 = NodeNumberTableX(2,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
        X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

        iL = Locate( X1, X1Arr, nX )

        ! --- Geometry Fields ---

        uGF(iNodeX,iX1,iX2,iX3,iGF_Alpha) &
          = Interpolate1D_Linear( X1, X1Arr(iL), X1Arr(iL+1), &
                                  AlphaArr(iL), AlphaArr(iL+1) )

        uGF(iNodeX,iX1,iX2,iX3,iGF_Psi) &
          = Interpolate1D_Linear( X1, X1Arr(iL), X1Arr(iL+1), &
                                  PsiArr(iL), PsiArr(iL+1) )

        uGF(iNodeX,iX1,iX2,iX3,iGF_h_1) &
          = uGF(iNodeX,iX1,iX2,iX3,iGF_Psi)**2
        uGF(iNodeX,iX1,iX2,iX3,iGF_h_2) &
          = uGF(iNodeX,iX1,iX2,iX3,iGF_Psi)**2 * X1
        uGF(iNodeX,iX1,iX2,iX3,iGF_h_3) &
          = uGF(iNodeX,iX1,iX2,iX3,iGF_Psi)**2 * X1 * SIN( X2 )

        CALL ComputeGeometryX_FromScaleFactors( uGF(:,iX1,iX2,iX3,:) )

        uGF(iNodeX,iX1,iX2,iX3,iGF_Beta_1) = Zero
        uGF(iNodeX,iX1,iX2,iX3,iGF_Beta_2) = Zero
        uGF(iNodeX,iX1,iX2,iX3,iGF_Beta_3) = Zero

        ! --- Fluid Fields ---

        uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
          = Interpolate1D_Linear( X1, X1Arr(iL), X1Arr(iL+1), &
                                  DensityArr(iL), DensityArr(iL+1) )

        uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = Zero
        uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = Zero
        uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = Zero

        uPF(iNodeX,iX1,iX2,iX3,iPF_E) &
          = Interpolate1D_Linear &
              ( X1, X1Arr(iL), X1Arr(iL+1), &
                PressureArr(iL), PressureArr(iL+1) ) / ( Gamma_IDEAL - One )

      END DO

      CALL ComputePressureFromPrimitive_IDEAL &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_E), &
               uPF(:,iX1,iX2,iX3,iPF_Ne), uAF(:,iX1,iX2,iX3,iAF_P) )

      CALL ComputeConserved_Euler_Relativistic &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_V1), &
               uPF(:,iX1,iX2,iX3,iPF_V2), uPF(:,iX1,iX2,iX3,iPF_V3), &
               uPF(:,iX1,iX2,iX3,iPF_E ), uPF(:,iX1,iX2,iX3,iPF_Ne), &
               uCF(:,iX1,iX2,iX3,iCF_D ), uCF(:,iX1,iX2,iX3,iCF_S1), &
               uCF(:,iX1,iX2,iX3,iCF_S2), uCF(:,iX1,iX2,iX3,iCF_S3), &
               uCF(:,iX1,iX2,iX3,iCF_E ), uCF(:,iX1,iX2,iX3,iCF_Ne), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33), &
               uAF(:,iX1,iX2,iX3,iAF_P) )

    END DO
    END DO
    END DO

    ! --- Apply reflecting boundary conditions to geometry fields (X1) ---

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = 1, swX(1)

        DO iNodeX3 = 1, nNodesX(3)
        DO iNodeX2 = 1, nNodesX(2)
        DO iNodeX1 = 1, nNodesX(1)

          jNodeX1 = ( nNodesX(1) - iNodeX1 ) + 1

          iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )
          jNodeX = NodeNumberX( jNodeX1, iNodeX2, iNodeX3 )

          DO iGF = 1, nGF

            uGF(iNodeX,iX_B0(1)-iX1,iX2,iX3,iGF) &
              = uGF(jNodeX,iX_B0(1),iX2,iX3,iGF)

          END DO

        END DO
        END DO
        END DO

    END DO
    END DO
    END DO

    DEALLOCATE( X1Arr, PsiArr, AlphaArr, DensityArr, PressureArr )

    ! --- Iterate to incorporate gravity in initial conditions ---

    ALLOCATE( SourceTerms_Poseidon(1:nDOFX,iX_B0(1):iX_E0(1), &
                                           iX_B0(2):iX_E0(2), &
                                           iX_B0(3):iX_E0(3),6) )

    CONVERGED = .FALSE.
    ITER = 0

    DO WHILE( .NOT. CONVERGED )

      ITER = ITER + 1

      CALL ComputeSourceTerms_Poseidon &
             ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, SourceTerms_Poseidon )

      dAlpha = MINVAL( uGF(:,:,:,:,iGF_Alpha) )
      dPsi   = MAXVAL( uGF(:,:,:,:,iGF_Psi  ) )

      CALL SolveGravity_XCFC_Poseidon &
             ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, SourceTerms_Poseidon )

      dAlpha = ABS( dAlpha - MINVAL( uGF(:,:,:,:,iGF_Alpha) ) ) &
                 / MINVAL( uGF(:,:,:,:,iGF_Alpha) )
      dPsi   = ABS( dPsi   - MAXVAL( uGF(:,:,:,:,iGF_Psi)   ) ) &
                 / MAXVAL( uGF(:,:,:,:,iGF_Psi)   )

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E1(1)

        CALL ComputeConserved_Euler_Relativistic &
               ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_V1), &
                 uPF(:,iX1,iX2,iX3,iPF_V2), uPF(:,iX1,iX2,iX3,iPF_V3), &
                 uPF(:,iX1,iX2,iX3,iPF_E ), uPF(:,iX1,iX2,iX3,iPF_Ne), &
                 uCF(:,iX1,iX2,iX3,iCF_D ), uCF(:,iX1,iX2,iX3,iCF_S1), &
                 uCF(:,iX1,iX2,iX3,iCF_S2), uCF(:,iX1,iX2,iX3,iCF_S3), &
                 uCF(:,iX1,iX2,iX3,iCF_E ), uCF(:,iX1,iX2,iX3,iCF_Ne), &
                 uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33), &
                 uAF(:,iX1,iX2,iX3,iAF_P) )

      END DO
      END DO
      END DO

      IF( MAX( dAlpha, dPsi ) .LT. 1.0e-13_DP ) CONVERGED = .TRUE.

      IF( ITER .EQ. 50 )THEN

        WRITE(*,*) 'Could not initialize fields. Exiting...'
        STOP

      END IF

    END DO

    DEALLOCATE( SourceTerms_Poseidon )

  END SUBROUTINE InitializeFields_DynamicTOV


  ! --- Auxiliary utilities for TOV problem ---

  SUBROUTINE IntegrateOutwards &
    ( dX1, CentralPressure, Psi0, Alpha0, &
      GravitationalMass, Radius, dAlpha, dPsi, Alpha_A, Psi_A, nX )

    REAL(DP), INTENT(in)  :: dX1, CentralPressure, Psi0, Alpha0
    REAL(DP), INTENT(out) :: GravitationalMass, Radius, &
                             dAlpha, dPsi, Alpha_A, Psi_A
    INTEGER,  INTENT(out) :: nX

    REAL(DP) :: Pressure , E1 , E2 , Psi , Phi, X1
    REAL(DP) :: PressureN, E1N, E2N
    REAL(DP) :: Alpha

    ! --- Set inner boundary values ---

    X1                = SqrtTiny * Kilometer
    Pressure          = CentralPressure
    E1                = Zero
    E2                = Zero
    Psi               = Psi0
    Phi               = Alpha0 * Psi0
    GravitationalMass = Zero

    nX = 1

    DO WHILE( Pressure .GT. 1.0e-8_DP * CentralPressure .AND. X1 .LT. 1.01e1_DP * Kilometer )

      ! --- Explicit steps ---

      PressureN = Pressure + dX1 * dpdr ( Pressure, Phi, Psi, E1, E2, X1 )
      E1N       = E1       + dX1 * dE1dr( Pressure, Psi, X1 )
      E2N       = E2       + dX1 * dE2dr( Pressure, Phi, Psi, X1 )

      ! --- Implicit steps ---

      X1  = X1  + dX1
      Psi = Psi + dX1 * dPsidr( PressureN, E1N, X1 )
      Phi = Phi + dX1 * dPhidr( PressureN, E2N, X1 )

      Pressure = PressureN
      E1       = E1N
      E2       = E2N

      GravitationalMass &
        = GravitationalMass &
            + FourPi * X1**2 * ( Enthalpy( Pressure ) + Two * Pressure ) &
                * Phi * Psi**5 * dX1

      nX = nX + 1

    END DO

    Radius = X1

    Alpha   &
      = Phi / Psi

    Alpha_A &
      =  ( One - GravitationalMass / ( Two * SpeedOfLight**2 * Radius ) ) &
       / ( One + GravitationalMass / ( Two * SpeedOfLight**2 * Radius ) )

    Psi     &
      = Psi

    Psi_A   &
      = One + GravitationalMass / ( Two * SpeedOfLight**2 * Radius )

    dAlpha = Alpha_A - Alpha
    dPsi   = Psi_A - Psi

  END SUBROUTINE IntegrateOutwards


  REAL(DP) FUNCTION dpdr( Pressure, Phi, Psi, E1, E2, X1  )

    REAL(DP), INTENT(in) :: Pressure, Phi, Psi, E1, E2, X1

    REAL(DP) :: Lapse

    Lapse = Phi / Psi

    dpdr = - Enthalpy( Pressure ) &
             * ( dPhidr( Pressure, E2, X1 ) &
                   - Lapse * dPsidr( Pressure, E1, X1 ) ) / Phi

    RETURN
  END FUNCTION dpdr


  REAL(DP) FUNCTION Enthalpy( Pressure )

    REAL(DP), INTENT(in) :: Pressure

    Enthalpy = ( Pressure / PolytropicConstant_TOV )**( One / Gamma_IDEAL ) &
                 * SpeedOfLight**2 + Pressure / ( Gamma_IDEAL - One ) + Pressure

    RETURN
  END FUNCTION Enthalpy


  REAL(DP) FUNCTION dPhidr( Pressure, E2, X1 )

    REAL(DP), INTENT(in) :: Pressure, E2, X1

    dPhidr = GravitationalConstant / ( Two * SpeedOfLight**4 ) * E2 / X1**2

    RETURN
  END FUNCTION dPhidr


  REAL(DP) FUNCTION dPsidr( Pressure, E1, X1 )

    REAL(DP), INTENT(in) :: Pressure, E1, X1

    dPsidr = -GravitationalConstant / ( Two * SpeedOfLight**4 ) * E1 / X1**2

    RETURN
  END FUNCTION dPsidr


  REAL(DP) FUNCTION dE1dr( Pressure, Psi, X1 )

    REAL(DP), INTENT(in) :: Pressure, Psi, X1

    dE1dr = FourPi * X1**2 * f_E1( Pressure ) * Psi**5

    RETURN
  END FUNCTION dE1dr


  REAL(DP) FUNCTION dE2dr( Pressure, Phi, Psi, X1 )

    REAL(DP), INTENT(in) :: Pressure, Phi, Psi, X1

    dE2dr = FourPi * X1**2 * f_E2( Pressure ) * Phi * Psi**4

    RETURN
  END FUNCTION dE2dr


  REAL(DP) FUNCTION f_E1( Pressure )

    REAL(DP), INTENT(in) :: Pressure

    f_E1 = Enthalpy( Pressure ) - Pressure

    RETURN
  END FUNCTION f_E1


  REAL(DP) FUNCTION f_E2( Pressure )

    REAL(DP), INTENT(in) :: Pressure

    f_E2 = f_E1( Pressure ) + 6.0_DP * Pressure

    RETURN
  END FUNCTION f_E2


  REAL(DP) FUNCTION LorentzFactor( Psi, V )

    REAL(DP), INTENT(in) :: Psi, V

    LorentzFactor = One / SQRT( One - Psi**4 * ( V / SpeedOfLight )**2 )

    RETURN
  END FUNCTION LorentzFactor

  ! --- End of auxiliary utilities for TOV problem ---

END MODULE InitializationModule_Relativistic
