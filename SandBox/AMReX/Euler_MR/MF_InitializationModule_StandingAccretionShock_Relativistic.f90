MODULE MF_InitializationModule_StandingAccretionShock_Relativistic

  ! --- AMReX Modules ---

  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_amrcore_module, ONLY: &
    amrex_geom
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor, &
    amrex_parallel_communicator, amrex_parallel_myproc, &
    amrex_parallel_reduce_min
  USE amrex_parmparse_module, ONLY: &
    amrex_parmparse, &
    amrex_parmparse_build, &
    amrex_parmparse_destroy

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    nNodesX, &
    swX, &
    nDimsX
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX, &
    WeightsX1
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    nGF, &
    iGF_Alpha, &
    iGF_Psi, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33
  USE GeometryComputationModule, ONLY: &
    LapseFunction, &
    ConformalFactor
  USE FluidFieldsModule, ONLY: &
    nCF, &
    iCF_D, &
    iCF_S1, &
    iCF_S2, &
    iCF_S3, &
    iCF_E, &
    iPF_Ne, &
    nPF, &
    iPF_D, &
    iPF_V1, &
    iPF_V2, &
    iPF_V3, &
    iPF_E, &
    iCF_Ne, &
    nAF, &
    iAF_P
  USE Euler_BoundaryConditionsModule, ONLY: &
    ExpD, &
    ExpE
  USE Euler_UtilitiesModule, ONLY: &
    ComputeConserved_Euler
  USE EquationOfStateModule, ONLY: &
    ComputePressureFromPrimitive
  USE UnitsModule, ONLY: &
    Kilometer, &
    Second, &
    SolarMass, &
    Gram, &
    Centimeter, &
    Erg, &
    SpeedOfLight, &
    GravitationalConstant, &
    Millisecond
  USE UtilitiesModule, ONLY: &
    NodeNumberX, &
    Locate
  USE Euler_ErrorModule, ONLY: &
    DescribeError_Euler

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero, &
    SqrtTiny, &
    Half, &
    One, &
    Two, &
    Three, &
    Four, &
    Five, &
    Pi, &
    TwoPi, &
    FourPi
  USE InputParsingModule, ONLY: &
    nLevels, &
    xL, &
    xR, &
    Gamma_IDEAL, &
    UseTiling, &
    t_end, &
    StepNo, &
    iOS_CPP
!!$  USE MF_AccretionShockUtilitiesModule, ONLY: &
!!$    FileName_Nodal1DIC_SAS
  USE InputOutputModuleAMReX, ONLY: &
    WriteFieldsAMReX_PlotFile
  USE MF_FieldsModule, ONLY: &
    MF_uPF, &
    MF_uAF, &
    MF_uDF
  USE MF_Euler_UtilitiesModule, ONLY: &
    ComputeFromConserved_Euler_MF

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFields_Euler_StandingAccretionShock_Relativistic_MF

CONTAINS


  SUBROUTINE InitializeFields_Euler_StandingAccretionShock_Relativistic_MF &
    ( iLevel, MF_uGF, MF_uCF )

    INTEGER             , INTENT(in)    :: iLevel
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF

    ! --- thornado ---

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNX, iNX1, iNX2
    REAL(DP) :: X1, X2
    REAL(DP) :: uGF_K(nDOFX,nGF)
    REAL(DP) :: uCF_K(nDOFX,nCF)
    REAL(DP) :: uPF_K(nDOFX,nPF)
    REAL(DP) :: uAF_K(nDOFX,nAF)

    ! --- AMReX ---

    INTEGER                       :: lo_G(4), hi_G(4), nX(3)
    INTEGER                       :: lo_F(4), hi_F(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    TYPE(amrex_parmparse)         :: PP

    ! --- Problem-dependent Parameters ---

    REAL(DP) :: MassPNS, ShockRadius, AccretionRate, PolytropicConstant
    LOGICAL  :: ApplyPerturbation
    INTEGER  :: PerturbationOrder
    REAL(DP) :: PerturbationAmplitude
    REAL(DP) :: rPerturbationInner
    REAL(DP) :: rPerturbationOuter

    INTEGER  :: iX1_1, iX1_2, iNX1_1, iNX1_2, indC
    INTEGER  :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), iX_B(3), iX_E(3)
    REAL(DP) :: X1_1, X1_2, D_1, D_2, V_1, V_2, P_1, P_2, uPert
    REAL(DP) :: D0, V0, P0
    REAL(DP) :: K_1, K_2, Mdot, AdvectionTime, rCenter, sigma
    REAL(DP) :: Alpha
    REAL(DP) :: Psi
    REAL(DP), ALLOCATABLE :: D    (:,:)
    REAL(DP), ALLOCATABLE :: V    (:,:)
    REAL(DP), ALLOCATABLE :: P    (:,:)
    REAL(DP), ALLOCATABLE :: Field(:,:)
    LOGICAL               :: InitializeFromFile, ResetEndTime
    INTEGER, PARAMETER    :: nX_LeastSquares = 5
    CHARACTER(LEN=:), ALLOCATABLE :: PerturbationType

    ApplyPerturbation     = .FALSE.
    PerturbationOrder     = 0
    PerturbationAmplitude = Zero
    rPerturbationInner    = Zero
    rPerturbationOuter    = Zero
    InitializeFromFile    = .FALSE.
    ResetEndTime          = .FALSE.
    rCenter               = Zero
    sigma                 = SqrtTiny
    PerturbationType      = 'StepFunction'
    CALL amrex_parmparse_build( PP, 'SAS' )
      CALL PP % get  ( 'Mass', &
                        MassPNS )
      CALL PP % get  ( 'AccretionRate', &
                        AccretionRate )
      CALL PP % get  ( 'ShockRadius', &
                        ShockRadius )
      CALL PP % query( 'ApplyPerturbation', &
                        ApplyPerturbation )
      CALL PP % query( 'PerturbationOrder', &
                        PerturbationOrder )
      CALL PP % query( 'PerturbationAmplitude', &
                        PerturbationAmplitude )
      CALL PP % query( 'rPerturbationInner', &
                        rPerturbationInner )
      CALL PP % query( 'rPerturbationOuter', &
                        rPerturbationOuter )
      CALL PP % query( 'rCenter', &
                       rCenter )
      CALL PP % query( 'sigma', &
                       sigma )
      CALL PP % query( 'InitializeFromFile', &
                        InitializeFromFile )
      CALL PP % query( 'ResetEndTime', &
                        ResetEndTime )
      CALL PP % query( 'PerturbationType', &
                        PerturbationType )
    CALL amrex_parmparse_destroy( PP )

    MassPNS            = MassPNS            * SolarMass
    AccretionRate      = AccretionRate      * ( SolarMass / Second )
    ShockRadius        = ShockRadius        * Kilometer
    rPerturbationInner = rPerturbationInner * Kilometer
    rPerturbationOuter = rPerturbationOuter * Kilometer
    rCenter            = rCenter            * ShockRadius
    sigma              = ( ShockRadius - rCenter ) / sigma
    PolytropicConstant = 2.0e14_DP &
                           * ( Erg / Centimeter**3 &
                           / ( Gram / Centimeter**3 )**( Gamma_IDEAL ) )

    Mdot = AccretionRate
    K_1  = PolytropicConstant

    IF( iLevel .EQ. 0 .AND. amrex_parallel_ioprocessor() )THEN

      WRITE(*,*)

      WRITE(*,'(6x,A,L)') &
        'InitializeFromFile:              ', &
        InitializeFromFile

      WRITE(*,'(6x,A,ES9.2E3,A)') &
        'Shock radius:                    ', &
        ShockRadius / Kilometer, &
        ' km'

      WRITE(*,'(6x,A,ES9.2E3,A)') &
        'PNS Mass:                        ', &
        MassPNS / SolarMass, &
        ' Msun'

      WRITE(*,'(6x,A,ES9.2E3,A)') &
        'Accretion Rate:                  ', &
        AccretionRate / ( SolarMass / Second ), &
        ' Msun/s'

      WRITE(*,'(6x,A,ES9.2E3,A)') &
        'Polytropic Constant (pre-shock): ', &
        K_1 / ( ( Erg / Centimeter**3 ) &
          / ( Gram / Centimeter**3 )**( Gamma_IDEAL ) ), &
        ' erg/cm^3 / (g/cm^3)^( Gamma )'

      WRITE(*,*)

      WRITE(*,'(6x,A,L)') &
        'Apply Perturbation:           ', &
        ApplyPerturbation

      WRITE(*,'(6x,A,A)') &
        'Perturbation Type:            ', &
        TRIM( PerturbationType )

      WRITE(*,'(6x,A,I1)') &
        'Perturbation order:           ', &
        PerturbationOrder

      WRITE(*,'(6x,A,ES9.2E3)') &
        'Perturbation amplitude:       ', &
         PerturbationAmplitude

      WRITE(*,'(6x,A,ES9.2E3,A)') &
        'Inner radius of perturbation: ', &
        rPerturbationInner / Kilometer, ' km'

      WRITE(*,'(6x,A,ES9.2E3,A)') &
        'Outer radius of perturbation: ', &
        rPerturbationOuter / Kilometer, ' km'

      WRITE(*,'(6x,A,ES9.2E3,A)') &
        'Center of perturbation:       ', &
        rCenter / Kilometer, ' km'

      WRITE(*,'(6x,A,ES9.2E3,A)') &
        'Width of perturbation:        ', &
        sigma / Kilometer, ' km'

      WRITE(*,'(6x,A,L)') &
        'Reset end-time:               ', &
        ResetEndTime

    END IF

    uGF_K = Zero
    uCF_K = Zero
    uPF_K = Zero
    uAF_K = Zero

    iX_B1 = amrex_geom(iLevel) % domain % lo - swX
    iX_E1 = amrex_geom(iLevel) % domain % hi + swX

    nX = amrex_geom(iLevel) % domain % hi &
           - amrex_geom(iLevel) % domain % lo + 1

    ALLOCATE( D    (1:nNodesX(1),iX_B1(1):iX_E1(1)) )
    ALLOCATE( V    (1:nNodesX(1),iX_B1(1):iX_E1(1)) )
    ALLOCATE( P    (1:nNodesX(1),iX_B1(1):iX_E1(1)) )
    ALLOCATE( Field(1:nNodesX(1),iX_B1(1):iX_E1(1)) )

    IF( InitializeFromFile )THEN

!!$      CALL ReadFluidFieldsFromFile( iX_B1, iX_E1, D, V, P )

    ELSE

      ! --- Quantities with _1 are pre-shock, those with _2 are post-shock ---

      CALL LocateFirstUnShockedElement &
             ( iX_B1, iX_E1, ShockRadius, &
               iX1_1, iX1_2, iNX1_1, iNX1_2, X1_1, X1_2 )

      ! --- Pre-shock Fields ---

      X1 = NodeCoordinate( MeshX(1), iX_E1(1), nNodesX(1) )

      ! --- Use Newtonian values as initial guesses ---

      V0 = -SQRT( Two * GravitationalConstant * MassPNS / X1 )
      D0 = -Mdot / ( FourPi * X1**2 * V0 )
      P0 = K_1 * D0**( Gamma_IDEAL )

      DO iX1 = iX_E1(1), iX1_1, -1

        DO iNX1 = nNodesX(1), 1, -1

          X1 = NodeCoordinate( MeshX(1), iX1, iNX1 )

          IF( X1 .LE. ShockRadius ) CYCLE

          Alpha = LapseFunction  ( X1, MassPNS )
          Psi   = ConformalFactor( X1, MassPNS )

          CALL NewtonRaphson_SAS &
                 ( X1, MassPNS, Mdot, K_1, &
                   Alpha, Psi, D0, V0, P0, &
                   D(iNX1,iX1), V(iNX1,iX1), P(iNX1,iX1) )

          D0 = D(iNX1,iX1)
          V0 = V(iNX1,iX1)
          P0 = P(iNX1,iX1)

        END DO

      END DO

      ! --- Apply Jump Conditions ---

      D_1 = D(iNX1_1,iX1_1)
      V_1 = V(iNX1_1,iX1_1)
      P_1 = P(iNX1_1,iX1_1)

      CALL ApplyJumpConditions_SAS &
             ( MassPNS, Mdot, ShockRadius, &
               LapseFunction  ( X1_1, MassPNS ), &
               ConformalFactor( X1_1, MassPNS ), D_1, V_1, P_1, &
               LapseFunction  ( X1_2, MassPNS ), &
               ConformalFactor( X1_2, MassPNS ), D_2, V_2, P_2 )

      K_2 = P_2 / D_2**( Gamma_IDEAL )

      IF( amrex_parallel_ioprocessor() )THEN

        WRITE(*,*)
        WRITE(*,'(6x,A)') 'Jump Conditions'
        WRITE(*,'(6x,A)') '---------------'
        WRITE(*,*)
        WRITE(*,'(8x,A)') 'Pre-shock:'
        WRITE(*,'(10x,A,I4.4)')       'iX1      = ', iX1_1
        WRITE(*,'(10x,A,I2.2)')       'iNX1     = ', iNX1_1
        WRITE(*,'(10x,A,ES13.6E3,A)') 'X1       = ', X1_1 / Kilometer, '  km'
        WRITE(*,'(10x,A,ES13.6E3,A)') 'Density  = ', &
          D_1 / ( Gram / Centimeter**3 ), '  g/cm^3'
        WRITE(*,'(10x,A,ES14.6E3,A)') 'Velocity = ', &
          V_1 / ( Kilometer / Second ), ' km/s'
        WRITE(*,'(10x,A,ES13.6E3,A)') 'Pressure = ', &
          P_1 / ( Erg / Centimeter**3 ), '  erg/cm^3'
        WRITE(*,*)
        WRITE(*,'(8x,A)') 'Post-shock:'
        WRITE(*,'(10x,A,I4.4)')       'iX1      = ', iX1_2
        WRITE(*,'(10x,A,I2.2)')       'iNX1     = ', iNX1_2
        WRITE(*,'(10x,A,ES13.6E3,A)') 'X1       = ', X1_2 / Kilometer, '  km'
        WRITE(*,'(10x,A,ES13.6E3,A)') 'Density  = ', &
          D_2 / ( Gram / Centimeter**3 ), '  g/cm^3'
        WRITE(*,'(10x,A,ES14.6E3,A)') 'Velocity = ', &
          V_2 / ( Kilometer / Second ), ' km/s'
        WRITE(*,'(10x,A,ES13.6E3,A)') 'Pressure = ', &
          P_2 / ( Erg / Centimeter**3 ), '  erg/cm^3'
        WRITE(*,*)

      END IF

      ! --- Post-shock Fields ---

      AdvectionTime = Zero

      D0 = D_2
      V0 = V_2
      P0 = P_2

      DO iX1 = iX1_2, iX_B1(1), -1

        DO iNX1 = nNodesX(1), 1, -1

          X1 = NodeCoordinate( MeshX(1), iX1, iNX1 )

          IF( X1 .GT. ShockRadius ) CYCLE

          Alpha = LapseFunction  ( X1, MassPNS )
          Psi   = ConformalFactor( X1, MassPNS )

          CALL NewtonRaphson_SAS &
                 ( X1, MassPNS, Mdot, K_2, &
                   Alpha, Psi, D0, V0, P0, &
                   D(iNX1,iX1), V(iNX1,iX1), P(iNX1,iX1) )

          D0 = D(iNX1,iX1)
          V0 = V(iNX1,iX1)
          P0 = P(iNX1,iX1)

          AdvectionTime &
            = AdvectionTime &
                + WeightsX1(iNX1) * MeshX(1) % Width(iX1) / ABS( V(iNX1,iX1) )

        END DO

      END DO

      IF( ResetEndTime )THEN
        t_end = 50.0_DP * AdvectionTime
        WRITE(*,'(6x,A,ES9.2E3,A)') &
          't_end: ', &
          t_end / Millisecond, ' ms'
      END IF

    END IF

    indC = Locate( rCenter, MeshX(1) % Center, nX(1) )

    ! --- Map to 3D domain ---

    CALL amrex_mfiter_build( MFI, MF_uGF, tiling = UseTiling )

    DO WHILE( MFI % next() )

      uGF => MF_uGF % DataPtr( MFI )
      uCF => MF_uCF % DataPtr( MFI )

      BX = MFI % tilebox()

      lo_G = LBOUND( uGF )
      hi_G = UBOUND( uGF )

      lo_F = LBOUND( uCF )
      hi_F = UBOUND( uCF )

      iX_B0 = BX % lo
      iX_E0 = BX % hi
      iX_B1 = BX % lo - swX
      iX_E1 = BX % hi + swX

      iX_B(1) = iX_B0(1)
      iX_E(1) = iX_E0(1)

      IF( BX % lo(1) .EQ. amrex_geom(iLevel) % domain % lo(1) ) &
        iX_B(1) = iX_B1(1)
      IF( BX % hi(1) .EQ. amrex_geom(iLevel) % domain % hi(1) ) &
        iX_E(1) = iX_E1(1)

      ! IF( TRIM( PerturbedField ) .EQ. 'Pressure' )THEN
        Field = P
      ! END IF

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B (1), iX_E (1)

        uGF_K &
          = RESHAPE( uGF(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )

        DO iNX = 1, nDOFX

          iNX1 = NodeNumberTableX(1,iNX)
          iNX2 = NodeNumberTableX(2,iNX)

          uPF_K(iNX,iPF_D ) = D(iNX1,iX1)
          uPF_K(iNX,iPF_V1) = V(iNX1,iX1)
          uPF_K(iNX,iPF_V2) = Zero
          uPF_K(iNX,iPF_V3) = Zero
          uAF_K(iNX,iAF_P ) = P(iNX1,iX1)
          uPF_K(iNX,iPF_E ) = P(iNX1,iX1) / ( Gamma_IDEAL - One )

          IF( ApplyPerturbation )THEN

            X1 = NodeCoordinate( MeshX(1), iX1, iNX1 )
            X2 = NodeCoordinate( MeshX(2), iX2, iNX2 )

            uPert = Field(iNX1,iX1)

            IF( TRIM( PerturbationType ) .EQ. 'Gaussian' )THEN

              IF( X1 .LT. ShockRadius )THEN

                IF( nDimsX .EQ. 1 .OR. PerturbationOrder .EQ. 0 )THEN

                  uPert &
                    = Field(iNX1,iX1) &
                        + Field(1,indC) * PerturbationAmplitude &
                        * EXP( -( X1 - rCenter )**2 / ( Two * sigma**2 ) )

                ELSE

                  uPert &
                    = Field(iNX1,iX1) &
                        + Field(1,indC) * PerturbationAmplitude * COS( X2 ) &
                        * EXP( -( X1 - rCenter )**2 / ( Two * sigma**2 ) )

                END IF ! nDimsX .EQ. 1 .OR. PerturbationOrder .EQ. 0

              END IF ! X1 .LT. ShockRadius

            END IF ! Gaussian

            IF( TRIM( PerturbationType ) .EQ. 'StepFunction' )THEN

              IF( X1 .GE. rPerturbationInner &
                    .AND. X1 .LE. rPerturbationOuter )THEN

                IF( PerturbationOrder .EQ. 0 ) &
                  uPert &
                    = Field(iNX1,iX1) &
                        * ( One + PerturbationAmplitude )

                IF( PerturbationOrder .EQ. 1 ) &
                  uPert &
                    = Field(iNX1,iX1) &
                        * ( One + PerturbationAmplitude * COS( X2 ) )

              END IF ! rInner <= X1 <= rOuter

            END IF ! StepFunction

            ! IF( TRIM( PerturbedField ) .EQ. 'Pressure' )THEN
              uPF_K(iNX,iPF_E ) = uPert / ( Gamma_IDEAL - One )
            ! END IF
             uPF_K(iNX,iPF_D) = uPert

          END IF ! Apply perturbation

        END DO !iNX

        CALL ComputePressureFromPrimitive &
               ( uPF_K(:,iPF_D), uPF_K(:,iPF_E), uPF_K(:,iPF_Ne), &
                 uAF_K(:,iAF_P) )

        CALL ComputeConserved_Euler &
               ( uPF_K(:,iPF_D ), uPF_K(:,iPF_V1), uPF_K(:,iPF_V2), &
                 uPF_K(:,iPF_V3), uPF_K(:,iPF_E ), uPF_K(:,iPF_Ne), &
                 uCF_K(:,iCF_D ), uCF_K(:,iCF_S1), uCF_K(:,iCF_S2), &
                 uCF_K(:,iCF_S3), uCF_K(:,iCF_E ), uCF_K(:,iCF_Ne), &
                 uGF_K(:,iGF_Gm_dd_11), &
                 uGF_K(:,iGF_Gm_dd_22), &
                 uGF_K(:,iGF_Gm_dd_33), &
                 uAF_K(:,iAF_P) )

        uCF(iX1,iX2,iX3,lo_F(4):hi_F(4)) &
          = RESHAPE( uCF_K, [ hi_F(4) - lo_F(4) + 1 ] )

      END DO
      END DO
      END DO

    END DO

    DEALLOCATE( Field   )
    DEALLOCATE( P )
    DEALLOCATE( V )
    DEALLOCATE( D )

    IF( .NOT. InitializeFromFile ) &
      CALL ComputeExtrapolationExponents &
             ( MF_uCF, iLevel, &
               amrex_geom(iLevel) % domain % lo - swX, &
               amrex_geom(iLevel) % domain % hi + swX, &
               nX_LeastSquares )

    IF( amrex_parallel_ioprocessor() )THEN

      IF( .NOT. InitializeFromFile ) &
        WRITE(*,'(6x,A,ES13.6E3,A)') &
          'Advection time:  ', &
          AdvectionTime / Millisecond, ' ms'

      WRITE(*,'(6x,A,I2.2)') &
        'nX_LeastSquares: ', &
        nX_LeastSquares

    END IF

  END SUBROUTINE InitializeFields_Euler_StandingAccretionShock_Relativistic_MF


  ! --- Auxiliary utilities for SAS problem ---


  SUBROUTINE ComputeExtrapolationExponents &
    ( MF_uCF, iLevel, iLo, iHi, nX_LeastSquares )

    TYPE(amrex_multifab), INTENT(in) :: MF_uCF
    INTEGER             , INTENT(in) :: iLevel, iLo(3), iHi(3)
    INTEGER             , INTENT(in) :: nX_LeastSquares

    REAL(DP) :: lnR   (nNodesX(1),iLo(1):iHi(1))
    REAL(DP) :: lnD   (nNodesX(1),iLo(1):iHi(1))
    REAL(DP) :: lnE   (nNodesX(1),iLo(1):iHi(1))
    REAL(DP) :: lnR_LS(nNodesX(1),nX_LeastSquares)
    REAL(DP) :: lnD_LS(nNodesX(1),nX_LeastSquares)
    REAL(DP) :: lnE_LS(nNodesX(1),nX_LeastSquares)

    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    INTEGER  :: iX1, iNX1
    REAL(DP) :: n

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX
    INTEGER            :: iX_B0(3), iX_B1(3)

    LOGICAL :: FirstTime

    lnR    = HUGE( One )
    lnD    = HUGE( One )
    lnE    = HUGE( One )
    lnR_LS = HUGE( One )
    lnD_LS = HUGE( One )
    lnE_LS = HUGE( One )

    ! --- Make local copies of X1, D, and tau ---

    CALL amrex_mfiter_build( MFI, MF_uCF, tiling = UseTiling )

    FirstTime = .TRUE.

    DO WHILE( MFI % next() )

      uCF => MF_uCF % DataPtr( MFI )
      BX = MFI % tilebox()

      iX_B0 = BX % lo
      iX_B1 = iX_B0 - swX

      IF( iX_B1(1) .EQ. iLo(1) .AND. FirstTime )THEN

        DO iX1 = iLo(1), iLo(1)+nX_LeastSquares-1

          DO iNX1 = 1, nNodesX(1)

            lnR(iNX1,iX1) &
              = LOG( NodeCoordinate( MeshX(1), iX1, iNX1 ) )

            lnD(iNX1,iX1) &
              = LOG( uCF(iX1,iX_B0(2),iX_B0(3),nDOFX*(iCF_D-1)+iNX1) )

            lnE(iNX1,iX1) &
              = LOG( uCF(iX1,iX_B0(2),iX_B0(3),nDOFX*(iCF_E-1)+iNX1) )

          END DO ! iNX1 = 1, nNodesX(1)

        END DO ! iX1 = iLo(1), iLo91)+nX_LeastSquares-1

        FirstTime = .FALSE.

      END IF

    END DO

    CALL amrex_mfiter_destroy( MFI )

    lnR_LS = lnR(:,iLo(1):iLo(1)+nX_LeastSquares-1)
    lnD_LS = lnD(:,iLo(1):iLo(1)+nX_LeastSquares-1)
    lnE_LS = lnE(:,iLo(1):iLo(1)+nX_LeastSquares-1)


    DO iX1 = 1, nX_LeastSquares

      CALL amrex_parallel_reduce_min( lnR_LS(:,iX1), nNodesX(1) )
      CALL amrex_parallel_reduce_min( lnD_LS(:,iX1), nNodesX(1) )
      CALL amrex_parallel_reduce_min( lnE_LS(:,iX1), nNodesX(1) )

    END DO

    n = DBLE( nNodesX(1) ) * DBLE( nX_LeastSquares )

    ! --- Expression for exponents from:
    !     https://mathworld.wolfram.com/LeastSquaresFittingPowerLaw.html ---

    ExpD = -( n * SUM( lnR_LS * lnD_LS ) - SUM( lnR_LS ) * SUM( lnD_LS ) ) &
             / ( n * SUM( lnR_LS**2 ) - SUM( lnR_LS )**2 )

    ExpE = -( n * SUM( lnR_LS * lnE_LS ) - SUM( lnR_LS ) * SUM( lnE_LS ) ) &
             / ( n * SUM( lnR_LS**2 ) - SUM( lnR_LS )**2 )

  END SUBROUTINE ComputeExtrapolationExponents


!!$  SUBROUTINE ReadFluidFieldsFromFile( iX_B1, iX_E1, D, V, P )
!!$
!!$    INTEGER,  INTENT(in)  :: iX_B1(3), iX_E1(3)
!!$    REAL(DP), INTENT(out) :: D(1:,iX_B1(1):), V(1:,iX_B1(1):), P(1:,iX_B1(1):)
!!$
!!$    CHARACTER(LEN=16) :: FMT
!!$    INTEGER           :: iX1
!!$
!!$    OPEN( UNIT = 101, FILE = TRIM( FileName_Nodal1DIC_SAS ) )
!!$
!!$    READ(101,*) ExpD
!!$    READ(101,*) ExpE
!!$    READ(101,*) FMT
!!$
!!$    DO iX1 = iX_B1(1), iX_E1(1)
!!$
!!$      READ(101,TRIM(FMT)) D(:,iX1)
!!$      READ(101,TRIM(FMT)) V(:,iX1)
!!$      READ(101,TRIM(FMT)) P(:,iX1)
!!$
!!$    END DO
!!$
!!$    CLOSE( 101 )
!!$
!!$  END SUBROUTINE ReadFluidFieldsFromFile


  SUBROUTINE NewtonRaphson_SAS &
    ( X1, MassPNS, Mdot, K, Alpha, Psi, D0, V0, P0, D, V, P )

    REAL(DP), INTENT(in)  :: X1, MassPNS, Mdot, K, &
                             Alpha, Psi, D0, V0, P0
    REAL(DP), INTENT(out) :: D ,V ,P

    REAL(DP) :: W
    REAL(DP) :: Jac(3,3), invJac(3,3)
    REAL(DP) :: f(3), uO(3), uN(3), du(3)

    LOGICAL             :: CONVERGED
    INTEGER             :: ITER
    REAL(DP), PARAMETER :: Tolu = 1.0e-16_DP
    REAL(DP), PARAMETER :: Tolf = 1.0e-16_DP
    REAL(DP), PARAMETER :: BernoulliConstant = SpeedOfLight**2
    INTEGER,  PARAMETER :: MAX_ITER = 4 - INT( LOG( Tolu ) /  LOG( Two ) )

    REAL(DP) :: a1, b1, b2, c1

    a1 = FourPi * Alpha * Psi**6 * X1**2 / Mdot * D0 * V0
    b1 = Alpha * SpeedOfLight**2 / BernoulliConstant
    b2 = Alpha * Gamma_IDEAL / ( Gamma_IDEAL - One ) &
           * P0 /( D0 * BernoulliConstant )
    c1 = P0 / ( K * D0**( Gamma_IDEAL ) )

    uO(1) = One
    uO(2) = One
    uO(3) = One

    CONVERGED = .FALSE.
    ITER      = 0
    DO WHILE( .NOT. CONVERGED .AND. ITER .LT. MAX_ITER )

      ITER = ITER + 1

      W = LorentzFactor( Psi, V0 * uO(2) )

      f(1) = a1 * uO(1) * W * uO(2) + One
      f(2) = b1 * W + b2 / uO(1) * W * uO(3) - One
      f(3) = c1 * uO(3) - uO(1)**( Gamma_IDEAL )

      Jac(1,1) = a1 * W * uO(2)
      Jac(1,2) = a1 * uO(1) * W**3
      Jac(1,3) = Zero
      Jac(2,1) = -b2 / uO(1)**2 * W * uO(3)
      Jac(2,2) = ( b1 + b2 / uO(1) * uO(3) ) &
                   * W**3 * Psi**4 * V0**2 / SpeedOfLight**2 * uO(2)
      Jac(2,3) = b2 / uO(1) * W
      Jac(3,1) = -Gamma_IDEAL * uO(1)**( Gamma_IDEAL - One )
      Jac(3,2) = Zero
      Jac(3,3) = c1

      InvJac = Inv3x3( Jac )

      uN = uO - MATMUL( InvJac, f )

      du = uN - uO

      IF( MAXVAL( ABS( du / uO ) ) .LT. Tolu ) CONVERGED = .TRUE.

      uO = uN

    END DO

    D = uN(1) * D0
    V = uN(2) * V0
    P = uN(3) * P0

  END SUBROUTINE NewtonRaphson_SAS


  SUBROUTINE ApplyJumpConditions_SAS &
    ( MassPNS, AccretionRate, ShockRadius, &
      Alpha_1, Psi_1, D_1, V_1, P_1, Alpha_2, Psi_2, D_2, V_2, P_2 )

    REAL(DP), INTENT(in)  :: MassPNS, AccretionRate, ShockRadius, &
                             Alpha_1, Psi_1, D_1, V_1, P_1, Alpha_2, Psi_2
    REAL(DP), INTENT(out) :: D_2, V_2, P_2

    REAL(DP) :: D1, D2, D3

    REAL(DP) :: W_1, h_1, W_2
    REAL(DP) :: Jac(3,3), invJac(3,3)
    REAL(DP) :: f(3), uO(3), uN(3), du(3)

    LOGICAL             :: CONVERGED
    INTEGER             :: ITER
    REAL(DP), PARAMETER :: Tolu = 1.0e-16_DP
    REAL(DP), PARAMETER :: Tolf = 1.0e-16_DP
    INTEGER,  PARAMETER :: MAX_ITER = 4 - INT( LOG( Tolu ) /  LOG( Two ) )

    REAL(DP) :: a1, b1, b2, b3, c1, c2
    REAL(DP) :: D20, V20, P20

    W_1 = LorentzFactor( Psi_1, V_1 )
    h_1 = SpeedOfLight**2 + Gamma_IDEAL / ( Gamma_IDEAL - One ) * P_1 / D_1

    D1 = D_1 * W_1 * V_1
    D2 = D_1 * h_1 / SpeedOfLight**2 * W_1**2 * V_1**2 + Psi_1**( -4 ) * P_1
    D3 = D_1 * h_1 / Alpha_1 * W_1**2 * V_1

    D20 &
      = ( Gamma_IDEAL + One ) / ( Gamma_IDEAL - One ) &
          * ( AccretionRate / FourPi ) &
          * ( Two * GravitationalConstant * MassPNS )**( -Half ) &
          * ShockRadius**( -Three / Two )

    V20 &
      = - ( Gamma_IDEAL - One ) / ( Gamma_IDEAL + One ) &
          * ( Two * GravitationalConstant * MassPNS )**( Half ) &
          * ShockRadius**( -Half )

    P20 &
      = Two / ( Gamma_IDEAL + One ) &
          * ( AccretionRate / FourPi ) &
          * ( Two * GravitationalConstant * MassPNS )**( Half ) &
          * ShockRadius**( -Five / Two )

    a1 = D20 * V20 / D1
    b1 = D20 * V20**2 / D2
    b2 = SpeedOfLight**( -2 ) * Gamma_IDEAL / ( Gamma_IDEAL - One ) &
           * P20 * V20**2 / D2
    b3 = Psi_2**( -4 ) * P20 / D2
    c1 = D20 * SpeedOfLight**2 * V20 * Alpha_2**( -1 ) / D3
    c2 = Gamma_IDEAL / ( Gamma_IDEAL - One ) * P20 * V20 * Alpha_2**( -1 ) / D3

    uO(1) = One
    uO(2) = One
    uO(3) = One

    CONVERGED = .FALSE.
    ITER      = 0
    DO WHILE( .NOT. CONVERGED .AND. ITER .LT. MAX_ITER )

      ITER = ITER + 1

      W_2 = LorentzFactor( Psi_2, V20 * uO(2) )

      f(1) = a1 * uO(1) * W_2 * uO(2) - One
      f(2) = b1 * uO(1) * W_2**2 * uO(2)**2 + b2 * W_2**2 * uO(2)**2 * uO(3) &
               + b3 * uO(3) - One
      f(3) = c1 * uO(1) * W_2**2 * uO(2) + c2 * W_2**2 * uO(2) * uO(3) - One

      Jac(1,1) = a1 * W_2 * uO(2)
      Jac(1,2) = a1 * uO(1) * W_2**3
      Jac(1,3) = Zero
      Jac(2,1) = b1 * W_2**2 * uO(2)**2
      Jac(2,2) = Two * ( b1 * uO(1) + b2 * uO(3) ) * W_2**4 * uO(2)
      Jac(2,3) = b2 * W_2**2 * uO(2)**2 + b3
      Jac(3,1) = c1 * W_2**2 * uO(2)
      Jac(3,2) = ( c1 * uO(1) + c2 * uO(3) ) * ( Two * W_2**2 - One ) * W_2**2
      Jac(3,3) = c2 * W_2**2 * uO(2)

      InvJac = Inv3x3( Jac )

      uN = uO - MATMUL( InvJac, f )

      du = uN - uO

      IF( MAXVAL( ABS( du / uO ) ) .LT. Tolu ) CONVERGED = .TRUE.

      uO = uN

    END DO

    D_2 = uN(1) * D20
    V_2 = uN(2) * V20
    P_2 = uN(3) * P20

  END SUBROUTINE ApplyJumpConditions_SAS


  ! --- From: http://fortranwiki.org/fortran/show/Matrix+inversion ---
  FUNCTION Inv3x3( A ) RESULT( invA )

    ! --- Performs a direct calculation of the inverse of a 3×3 matrix ---

    REAL(DP), INTENT(in) :: A   (3,3)
    REAL(DP)             :: invA(3,3)
    REAL(DP)             :: InvDet

    ! --- Calculate the inverse of the determinant of the matrix ---

    InvDet = One / ( A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)     &
                       - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1) &
                       + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1) )

    ! --- Calculate the inverse of the matrix ---

    invA(1,1) = +InvDet * ( A(2,2)*A(3,3) - A(2,3)*A(3,2) )
    invA(2,1) = -InvDet * ( A(2,1)*A(3,3) - A(2,3)*A(3,1) )
    invA(3,1) = +InvDet * ( A(2,1)*A(3,2) - A(2,2)*A(3,1) )
    invA(1,2) = -InvDet * ( A(1,2)*A(3,3) - A(1,3)*A(3,2) )
    invA(2,2) = +InvDet * ( A(1,1)*A(3,3) - A(1,3)*A(3,1) )
    invA(3,2) = -InvDet * ( A(1,1)*A(3,2) - A(1,2)*A(3,1) )
    invA(1,3) = +InvDet * ( A(1,2)*A(2,3) - A(1,3)*A(2,2) )
    invA(2,3) = -InvDet * ( A(1,1)*A(2,3) - A(1,3)*A(2,1) )
    invA(3,3) = +InvDet * ( A(1,1)*A(2,2) - A(1,2)*A(2,1) )

    RETURN
  END FUNCTION Inv3x3


  REAL(DP) FUNCTION LorentzFactor( Psi, V )

    REAL(DP), INTENT(in) :: Psi, V

    LorentzFactor = One / SQRT( One - Psi**4 * ( V / SpeedOfLight )**2 )

    RETURN
  END FUNCTION LorentzFactor


  SUBROUTINE LocateFirstUnShockedElement &
    ( iX_B1, iX_E1, ShockRadius, &
      iX1_1, iX1_2, iNX1_1, iNX1_2, X1_1, X1_2 )

    INTEGER,        INTENT(in)  :: iX_B1(3), iX_E1(3)
    REAL(DP),       INTENT(in)  :: ShockRadius
    INTEGER,        INTENT(out) :: iX1_1, iX1_2, iNX1_1, iNX1_2
    REAL(DP),       INTENT(out) :: X1_1, X1_2

    REAL(DP) :: X1, dX1
    INTEGER  :: iX1, iNX1
    LOGICAL  :: FirstPreShockElement = .FALSE.

    X1 = Zero

    DO iX1 = iX_B1(1), iX_E1(1)

      DO iNX1 = 1, nNodesX(1)

        dX1 = NodeCoordinate( MeshX(1), iX1, iNX1 ) - X1
        X1  = NodeCoordinate( MeshX(1), iX1, iNX1 )

        IF( X1 .LE. ShockRadius ) CYCLE

        IF( X1 .GT. ShockRadius .AND. .NOT. FirstPreShockElement )THEN

          iX1_1  = iX1
          iNX1_1 = iNX1
          X1_1   = X1
          X1_2   = X1 - dX1

          IF( iNX1_1 .EQ. 1 )THEN

            iX1_2  = iX1_1 - 1
            iNX1_2 = nNodesX(1)

          ELSE

            iX1_2  = iX1_1
            iNX1_2 = iNX1_1 - 1

          END IF

          FirstPreShockElement = .TRUE.

        END IF

      END DO

    END DO

  END SUBROUTINE LocateFirstUnShockedElement


END MODULE MF_InitializationModule_StandingAccretionShock_Relativistic
