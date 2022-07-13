MODULE InitializationModule_Relativistic

  USE KindModule, ONLY: &
    DP,       &
    Zero,     &
    Half,     &
    One,      &
    Two,      &
    Three,    &
    Four,     &
    Pi,       &
    TwoPi,    &
    FourPi
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    nNodesX,     &
    nDimsX,      &
    nDOFX,       &
    iX_B0,       &
    iX_B1,       &
    iX_E0,       &
    iX_E1
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX, &
    WeightsX_q
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    uGF,          &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33
  USE FluidFieldsModule, ONLY: &
    nPF,    &
    uPF,    &
    iPF_D,  &
    iPF_V1, &
    iPF_V2, &
    iPF_V3, &
    iPF_E,  &
    iPF_Ne, &
    uCF,    &
    iCF_D,  &
    iCF_S1, &
    iCF_S2, &
    iCF_S3, &
    iCF_E,  &
    iCF_Ne, &
    uAF,    &
    iAF_P
  USE EquationOfStateModule_IDEAL, ONLY: &
    Gamma_IDEAL, &
    ComputePressureFromPrimitive_IDEAL
  USE Euler_UtilitiesModule_Relativistic, ONLY: &
    ComputeConserved_Euler_Relativistic
  USE UnitsModule, ONLY: &
    SpeedOfLight
  USE QuadratureModule, ONLY: &
    GetQuadrature
  USE PolynomialBasisModule_Lagrange, ONLY: &
    LagrangeP

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFields_Relativistic


CONTAINS


  SUBROUTINE InitializeFields_Relativistic &
               ( AdvectionProfile_Option, &
                 RiemannProblemName_Option, &
                 nDetCells_Option, Eblast_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: AdvectionProfile_Option
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: RiemannProblemName_Option
    INTEGER,          INTENT(in), OPTIONAL :: nDetCells_Option
    REAL(DP),         INTENT(in), OPTIONAL :: Eblast_Option

    CHARACTER(LEN=64) :: AdvectionProfile = 'SineWave'
    CHARACTER(LEN=64) :: RiemannProblemName = 'Sod'

    ! --- Sedov-Taylor Blast Wave (Defaults) ---
    INTEGER  :: nDetCells = 1
    REAL(DP) :: Eblast    = 1.0e-3_DP

    uPF(:,:,:,:,iPF_Ne) = Zero

    IF( PRESENT( AdvectionProfile_Option ) ) &
      AdvectionProfile = TRIM( AdvectionProfile_Option )

    IF( PRESENT( RiemannProblemName_Option ) ) &
      RiemannProblemName = TRIM( RiemannProblemName_Option )

    IF( PRESENT( nDetCells_Option ) ) &
      nDetCells = nDetCells_Option
    IF( PRESENT( Eblast_Option ) ) &
      Eblast = Eblast_Option

    WRITE(*,*)
    WRITE(*,'(A,A)') '    INFO: ', TRIM( ProgramName )

    SELECT CASE ( TRIM( ProgramName ) )

      CASE( 'SlopeLimiterTest' )

        CALL InitializeFields_SlopeLimiterTest

      CASE( 'Advection' )

        CALL InitializeFields_Advection &
               ( TRIM( AdvectionProfile ) )

      CASE( 'Advection2D' )

        CALL InitializeFields_Advection2D &
               ( TRIM( AdvectionProfile ) )

      CASE( 'RiemannProblem' )

        CALL InitializeFields_RiemannProblem &
               ( TRIM( RiemannProblemName ), &
                 nDetCells_Option = nDetCells, &
                 Eblast_Option    = Eblast )

      CASE( 'RiemannProblem2D' )

        CALL InitializeFields_RiemannProblem2D &
               ( TRIM( RiemannProblemName ) )

      CASE( 'RiemannProblemSpherical' )

        CALL InitializeFields_RiemannProblemSpherical &
               ( TRIM( RiemannProblemName ) )

      CASE( 'SedovTaylorBlastWave' )

        CALL InitializeFields_SedovTaylorBlastWave &
               ( nDetCells, Eblast )

      CASE( 'KelvinHelmholtzInstability' )

         CALL InitializeFields_KelvinHelmholtzInstability

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,'(A21,A)') 'Invalid ProgramName: ', ProgramName
        WRITE(*,'(A)')     'Stopping...'
        STOP

    END SELECT

#if defined(THORNADO_OMP_OL)
  !$OMP TARGET UPDATE TO( uCF )
#elif defined(THORNADO_OACC)
  !$ACC UPDATE DEVICE   ( uCF )
#endif

  END SUBROUTINE InitializeFields_Relativistic


  SUBROUTINE InitializeFields_SlopeLimiterTest

    INTEGER       :: iX1, iX2, iX3, iPF
    INTEGER       :: iNodeX, iNodeX1, iNodeX2, iNodeX3
    REAL(DP)      :: X, X1, X2
    REAL(DP)      :: a0, a1, a2, a3, a4
    REAL(DP)      :: X0, theta
    CHARACTER(32) :: Problem

    INTEGER  :: M, M1, M2, M3, nDOFQ
    INTEGER  :: qNodeX, qNodeX1, qNodeX2, qNodeX3
    REAL(DP) :: etaG(nNodesX(1)), wG(nNodesX(1))
    REAL(DP), ALLOCATABLE :: etaQ_X1(:), wQ_X1(:)
    REAL(DP), ALLOCATABLE :: etaQ_X2(:), wQ_X2(:)
    REAL(DP), ALLOCATABLE :: etaQ_X3(:), wQ_X3(:)
    REAL(DP), ALLOCATABLE :: etaQ(:),    wQ(:)
    REAL(DP), ALLOCATABLE :: u0(:,:)
    INTEGER,  ALLOCATABLE :: NodeNumberTableQ(:,:)

    REAL(DP), ALLOCATABLE :: InterpolationMatrix(:,:)

    M = 5
    nDOFQ = M**nDimsX

    M1 = M
    M2 = 1
    M3 = 1
    IF( nDimsX .GT. 1 ) M2 = M
    IF( nDimsX .GT. 2 ) M3 = M

    nDOFQ = M1 * M2 * M3

    ALLOCATE( etaQ_X1(M1), wQ_X1(M1) )
    ALLOCATE( etaQ_X2(M2), wQ_X2(M2) )
    ALLOCATE( etaQ_X3(M3), wQ_X3(M3) )
    ALLOCATE( etaQ(nDOFQ), wQ(nDOFQ) )
    ALLOCATE( u0(nDOFQ,nPF) )
    ALLOCATE( InterpolationMatrix(nDOFX,nDOFQ) )
    ALLOCATE( NodeNumberTableQ(3,nDOFQ) )

    CALL GetQuadrature( M1        , etaQ_X1, wQ_X1 )
    CALL GetQuadrature( M2        , etaQ_X2, wQ_X2 )
    CALL GetQuadrature( M3        , etaQ_X3, wQ_X3 )
    CALL GetQuadrature( nNodesX(1), etaG, wG )

    qNodeX = 0
    DO qNodeX3 = 1, M3
    DO qNodeX2 = 1, M2
    DO qNodeX1 = 1, M1

      qNodeX = qNodeX + 1

      NodeNumberTableQ(1:3,qNodeX) &
        = [ qNodeX1, qNodeX2, qNodeX3 ]

      etaQ(qNodeX) = etaQ_X1(qNodeX1) * etaQ_X2(qNodeX2) * etaQ_X3(qNodeX3)
      wQ  (qNodeX) = wQ_X1  (qNodeX1) * wQ_X2  (qNodeX2) * wQ_X3  (qNodeX3)

    END DO
    END DO
    END DO

    DO iNodeX = 1, nDOFX

      iNodeX1 = NodeNumberTableX(1,iNodeX)
      iNodeX2 = NodeNumberTableX(2,iNodeX)
      iNodeX3 = NodeNumberTableX(3,iNodeX)

      DO qNodeX = 1, nDOFQ

        qNodeX1 = NodeNumberTableQ(1,qNodeX)
        qNodeX2 = NodeNumberTableQ(2,qNodeX)
        qNodeX3 = NodeNumberTableQ(3,qNodeX)

        InterpolationMatrix(iNodeX,qNodeX) &
          = wQ(qNodeX) &
              * LagrangeP( etaQ_X1(qNodeX1), iNodeX1, etaG, nNodesX(1) ) &
              * LagrangeP( etaQ_X2(qNodeX2), iNodeX2, etaG, nNodesX(2) ) &
              * LagrangeP( etaQ_X3(qNodeX3), iNodeX3, etaG, nNodesX(3) )

      END DO

    END DO

    Problem = 'SmoothTANH'

    ! --- Coefficients for polynomial problem ---
    a0 = 1.0_DP
    a1 = 0.1_DP
    a2 = 4.0_DP
    a3 = -0.4_DP
    a4 = -1.0_DP

    ! --- Scale length for TANH problem ---
    X0 = 0.01_DP

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO qNodeX = 1, nDOFQ

        qNodeX1 = NodeNumberTableQ(1,qNodeX)
        qNodeX2 = NodeNumberTableQ(2,qNodeX)

        X1 = etaQ_X1(qNodeX1) * MeshX(1) % Width(iX1) + MeshX(1) % Center(iX1)
        X2 = etaQ_X2(qNodeX2) * MeshX(2) % Width(iX2) + MeshX(2) % Center(iX2)

        X = X1

        u0(qNodeX,iPF_V1) = Zero
        u0(qNodeX,iPF_V2) = Zero
        u0(qNodeX,iPF_V3) = Zero
        u0(qNodeX,iPF_E ) = One / ( Gamma_IDEAL - One )
        u0(qNodeX,iPF_Ne) = Zero

        SELECT CASE( Problem )

          CASE( 'Polynomial' )

            u0(qNodeX,iPF_D) = a0*X**0 + a1*X**1 + a2*X**2 + a3*X**3 + a4*X**4

          CASE( 'SmoothTANH' )

            u0(qNodeX,iPF_D) = 2.0d0 + TANH( X / X0 )

          CASE( 'SinCosTANH' )

            theta = Half * ( One - TANH( X / X0 ) )

            u0(qNodeX,iPF_D) = 1.0d0 + theta * COS( TwoPi * X ) &
                                 + ( One - theta ) * SIN( TwoPi * X )

          CASE( 'IsolatedContact' )

            IF( X .LT. Zero )THEN

               u0(qNodeX,iPF_D) = 0.1_DP

           ELSE

               u0(qNodeX,iPF_D) = 0.5_DP

           END IF

        END SELECT

      END DO

      DO iPF = 1, nPF

        uPF(:,iX1,iX2,iX3,iPF) &
          = MATMUL( InterpolationMatrix, u0(:,iPF) ) / WeightsX_q

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

    DEALLOCATE( NodeNumberTableQ )
    DEALLOCATE( InterpolationMatrix )
    DEALLOCATE( u0 )
    DEALLOCATE( etaQ   , wQ )
    DEALLOCATE( etaQ_X3, wQ_X3 )
    DEALLOCATE( etaQ_X2, wQ_X2 )
    DEALLOCATE( etaQ_X1, wQ_X1 )

  END SUBROUTINE InitializeFields_SlopeLimiterTest


  SUBROUTINE InitializeFields_Advection( AdvectionProfile )

    CHARACTER(LEN=*), INTENT(in) :: AdvectionProfile

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX, iNodeX1
    REAL(DP) :: X1

    WRITE(*,*)
    WRITE(*,'(A4,A,A)') &
      '', 'Advection Profile: ', TRIM( AdvectionProfile )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

        SELECT CASE( TRIM( AdvectionProfile ) )

          CASE( 'SineWave' )

            uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = One + 0.1_DP * SIN( TwoPi * X1 )
            uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.1_DP
            uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
            uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
            uAF(iNodeX,iX1,iX2,iX3,iAF_P ) = 1.0_DP
            uPF(iNodeX,iX1,iX2,iX3,iPF_E )  &
              = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - One )

          CASE( 'TopHat' )

            IF( X1 .GT. 0.45 .AND. X1 .LT. 0.55 )THEN

              uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 2.0_DP

            ELSE

              uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 1.0_DP

            END IF

              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.1_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
              uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = 1.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_E)  &
                = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - One )

          CASE DEFAULT

            WRITE(*,*)
            WRITE(*,'(A,A)') &
              'Invalid choice for AdvectionProfile: ', AdvectionProfile
            WRITE(*,'(A)') 'Valid choices:'
            WRITE(*,'(A)') '  SineWave'
            WRITE(*,'(A)') '  TopHat'
            WRITE(*,*)
            WRITE(*,'(A)') 'Stopping...'
            STOP

        END SELECT

      END DO

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

  END SUBROUTINE InitializeFields_Advection


  SUBROUTINE InitializeFields_Advection2D( AdvectionProfile )

    CHARACTER(LEN=*), INTENT(in) :: AdvectionProfile

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX, iNodeX1, iNodeX2
    REAL(DP) :: X1, X2

    WRITE(*,*)
    WRITE(*,'(A4,A,A)') &
      '', 'Advection Profile: ', TRIM( AdvectionProfile )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)
        iNodeX2 = NodeNumberTableX(2,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
        X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

        SELECT CASE( TRIM( AdvectionProfile ) )

          CASE( 'SineWaveX1' )

            uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = One + 0.1_DP * SIN( TwoPi * X1 )
            uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.1_DP
            uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
            uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
            uAF(iNodeX,iX1,iX2,iX3,iAF_P ) = 1.0_DP
            uPF(iNodeX,iX1,iX2,iX3,iPF_E )  &
              = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - One )

          CASE( 'SineWaveX2' )

            uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = One + 0.1_DP * SIN( TwoPi * X2 )
            uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
            uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.1_DP
            uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
            uAF(iNodeX,iX1,iX2,iX3,iAF_P ) = 1.0_DP
            uPF(iNodeX,iX1,iX2,iX3,iPF_E )  &
              = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - One )

          CASE( 'SineWaveX1X2' )

            uPF(iNodeX,iX1,iX2,iX3,iPF_D)  &
              = One + 0.1_DP * SIN( SQRT( Two ) * TwoPi * ( X1 + X2 ) )
            uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.1_DP / SQRT( Two )
            uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.1_DP / SQRT( Two )
            uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
            uAF(iNodeX,iX1,iX2,iX3,iAF_P ) = 1.0_DP
            uPF(iNodeX,iX1,iX2,iX3,iPF_E )  &
              = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - One )

          CASE DEFAULT

            WRITE(*,*)
            WRITE(*,'(A,A)') &
              'Invalid choice for AdvectionProfile: ', AdvectionProfile
            WRITE(*,'(A)') 'Valid choices:'
            WRITE(*,'(A)') '  SineWaveX1'
            WRITE(*,'(A)') '  SineWaveX2'
            WRITE(*,'(A)') '  SineWaveX1X2'
            WRITE(*,*)
            WRITE(*,'(A)') 'Stopping...'
            STOP

        END SELECT

      END DO

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

  END SUBROUTINE InitializeFields_Advection2D


  SUBROUTINE InitializeFields_RiemannProblem &
               ( RiemannProblemName, &
                 nDetCells_Option, Eblast_Option )

    CHARACTER(LEN=*), INTENT(in)           :: RiemannProblemName
    INTEGER,          INTENT(in), OPTIONAL :: nDetCells_Option
    REAL(DP),         INTENT(in), OPTIONAL :: Eblast_Option

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX, iNodeX1
    REAL(DP) :: X1, XD, Vs

    REAL(DP) :: LeftState(nPF), RightState(nPF)

    WRITE(*,*)
    WRITE(*,'(A4,A,A)') &
      '', 'Riemann Problem Name: ', TRIM( RiemannProblemName )
    WRITE(*,*)

    SELECT CASE( TRIM( RiemannProblemName ) )

      CASE( 'Sod' )

        XD = 0.5_DP

        LeftState(iPF_D ) = 1.0_DP
        LeftState(iPF_V1) = 0.0_DP
        LeftState(iPF_V2) = 0.0_DP
        LeftState(iPF_V3) = 0.0_DP
        LeftState(iPF_E ) = 1.0_DP / ( Gamma_IDEAL - One )

        RightState(iPF_D ) = 0.125_DP
        RightState(iPF_V1) = 0.0_DP
        RightState(iPF_V2) = 0.0_DP
        RightState(iPF_V3) = 0.0_DP
        RightState(iPF_E ) = 0.1_DP / ( Gamma_IDEAL - One )

      CASE( 'IsolatedShock' )

        Vs = 0.01_DP
        XD = Half

        RightState(iPF_D)  = 1.0_DP
        RightState(iPF_V1) = -0.9_DP
        RightState(iPF_V2) = 0.0_DP
        RightState(iPF_V3) = 0.0_DP
        RightState(iPF_E)  = 1.0_DP / ( Gamma_IDEAL - One )

        CALL ComputeLeftState &
               ( Vs,                 &
                 RightState(iPF_D ), &
                 RightState(iPF_V1), &
                 RightState(iPF_E ) * ( Gamma_IDEAL - One ), &
                 LeftState (iPF_D ), &
                 LeftState (iPF_V1), &
                 LeftState (iPF_E ) )

        LeftState(iPF_V2) = 0.0_DP
        LeftState(iPF_V3) = 0.0_DP

      CASE( 'IsolatedContact' )

        Vs = 0.01_DP
        XD = 0.5_DP

        LeftState(iPF_D ) = 5.9718209694880811e0_DP
        LeftState(iPF_V1) = Vs
        LeftState(iPF_V2) = 0.0_DP
        LeftState(iPF_V3) = 0.0_DP
        LeftState(iPF_E ) = 1.0_DP / ( Gamma_IDEAL - One )

        RightState(iPF_D ) = 1.0_DP
        RightState(iPF_V1) = Vs
        RightState(iPF_V2) = 0.0_DP
        RightState(iPF_V3) = 0.0_DP
        RightState(iPF_E ) = 1.0_DP / ( Gamma_IDEAL - One )

      CASE( 'MBProblem1' )

        XD = 0.5_DP

        LeftState(iPF_D ) = 1.0_DP
        LeftState(iPF_V1) = 0.9_DP
        LeftState(iPF_V2) = 0.0_DP
        LeftState(iPF_V3) = 0.0_DP
        LeftState(iPF_E ) = 1.0_DP / ( Gamma_IDEAL - One )

        RightState(iPF_D ) = 1.0_DP
        RightState(iPF_V1) = 0.0_DP
        RightState(iPF_V2) = 0.0_DP
        RightState(iPF_V3) = 0.0_DP
        RightState(iPF_E ) = 10.0_DP / ( Gamma_IDEAL - One )

      CASE( 'MBProblem4' )

        XD = 0.5_DP

        LeftState(iPF_D ) = 1.0_DP
        LeftState(iPF_V1) = 0.0_DP
        LeftState(iPF_V2) = 0.0_DP
        LeftState(iPF_V3) = 0.0_DP
        LeftState(iPF_E ) = 1.0e3_DP / ( Gamma_IDEAL - One )

        RightState(iPF_D ) = 1.0_DP
        RightState(iPF_V1) = 0.0_DP
        RightState(iPF_V2) = 0.0_DP
        RightState(iPF_V3) = 0.0_DP
        RightState(iPF_E ) = 1.0e-2_DP / ( Gamma_IDEAL - One )

      CASE( 'PerturbedShockTube' )

        XD = 0.5_DP

        LeftState(iPF_D ) = 5.0_DP
        LeftState(iPF_V1) = 0.0_DP
        LeftState(iPF_V2) = 0.0_DP
        LeftState(iPF_V3) = 0.0_DP
        LeftState(iPF_E ) = 50.0_DP / ( Gamma_IDEAL - One )

        RightState(iPF_D ) = 0.0_DP ! --- Dummy ---
        RightState(iPF_V1) = 0.0_DP
        RightState(iPF_V2) = 0.0_DP
        RightState(iPF_V3) = 0.0_DP
        RightState(iPF_E ) = 5.0_DP / ( Gamma_IDEAL - One )

      CASE( 'ShockReflection' )

        XD = 1.0_DP

        LeftState(iPF_D ) = 1.0_DP
        LeftState(iPF_V1) = 0.99999_DP
        LeftState(iPF_V2) = 0.0_DP
        LeftState(iPF_V3) = 0.0_DP
        LeftState(iPF_E ) = 0.01_DP / ( Gamma_IDEAL - One )

        ! --- All of these are dummies ---
        RightState(iPF_D ) = 0.0_DP
        RightState(iPF_V1) = 0.0_DP
        RightState(iPF_V2) = 0.0_DP
        RightState(iPF_V3) = 0.0_DP
        RightState(iPF_E ) = 0.0_DP

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,'(A,A)') &
          'Invalid choice for RiemannProblemName: ', RiemannProblemName
        WRITE(*,'(A)') 'Valid choices:'
        WRITE(*,'(A)') &
          "  'Sod' - &
          Sod's shock tube"
        WRITE(*,'(A)') &
          "  'MBProblem1' - &
          Mignone & Bodo (2005) MNRAS, 364, 126, Problem 1"
        WRITE(*,'(A)') &
          "  'MBProblem4' - &
          Mignone & Bodo (2005) MNRAS, 364, 126, Problem 4"
        WRITE(*,'(A)') &
          "  'PerturbedShockTube' - &
          Del Zanna & Bucciantini (2002) AA, 390, 1177, &
          Sinusoidal density perturbation"
        WRITE(*,'(A)') &
          "  'ShockReflection' - &
          Del Zanna & Bucciantini (2002) AA, 390, 1177, &
          Planar shock reflection"
        WRITE(*,'(A)') 'Stopping...'
        STOP

    END SELECT

    IF( TRIM( RiemannProblemName ) .EQ. 'IsolatedShock' )THEN

      WRITE(*,'(6x,A,ES14.6E3)') 'Shock Velocity = ', Vs
      WRITE(*,*)

    END IF

    WRITE(*,'(6x,A,F8.6)') 'Gamma_IDEAL = ', Gamma_IDEAL
    WRITE(*,*)
    WRITE(*,'(6x,A,F8.6)') 'XD = ', XD
    WRITE(*,*)
    WRITE(*,'(6x,A)') 'Right State:'
    WRITE(*,*)
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_D  = ', RightState(iPF_D )
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_V1 = ', RightState(iPF_V1)
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_V2 = ', RightState(iPF_V2)
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_V3 = ', RightState(iPF_V3)
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_E  = ', RightState(iPF_E )
    WRITE(*,*)
    WRITE(*,'(6x,A)') 'Left State:'
    WRITE(*,*)
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_D  = ', LeftState(iPF_D )
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_V1 = ', LeftState(iPF_V1)
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_V2 = ', LeftState(iPF_V2)
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_V3 = ', LeftState(iPF_V3)
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_E  = ', LeftState(iPF_E )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

        IF( X1 .LE. XD )THEN

          uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = LeftState(iPF_D )
          uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = LeftState(iPF_V1)
          uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = LeftState(iPF_V2)
          uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = LeftState(iPF_V3)
          uPF(iNodeX,iX1,iX2,iX3,iPF_E ) = LeftState(iPF_E )

        ELSE

          uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = RightState(iPF_D )
          uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = RightState(iPF_V1)
          uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = RightState(iPF_V2)
          uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = RightState(iPF_V3)
          uPF(iNodeX,iX1,iX2,iX3,iPF_E ) = RightState(iPF_E )

          IF( TRIM( RiemannProblemName ) .EQ. 'PerturbedShockTube' ) &
            uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
              = 2.0_DP + 0.3_DP * SIN( 50.0_DP * X1 )

        END IF

      END DO

      CALL ComputePressureFromPrimitive_IDEAL &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_E ), &
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

  END SUBROUTINE InitializeFields_RiemannProblem


  SUBROUTINE InitializeFields_RiemannProblem2D( RiemannProblemName )

    CHARACTER(LEN=*), INTENT(in) :: RiemannProblemName

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX, iNodeX1, iNodeX2
    REAL(DP) :: X1, X2, X1D, X2D, Vs

    REAL(DP) :: NE(nPF), NW(nPF), SW(nPF), SE(nPF)

    WRITE(*,*)
    WRITE(*,'(A4,A,A)') &
      '', '2D Riemann Problem Name: ', TRIM( RiemannProblemName )
    WRITE(*,*)

    SELECT CASE( TRIM( RiemannProblemName ) )

      CASE( 'DzB2002' )

        X1D = 0.5_DP
        X2D = 0.5_DP

        NE(iPF_D ) = 0.1_DP
        NE(iPF_V1) = 0.0_DP
        NE(iPF_V2) = 0.0_DP
        NE(iPF_V3) = 0.0_DP
        NE(iPF_E ) = 0.01_DP / ( Gamma_IDEAL - One )

        NW(iPF_D ) = 0.1_DP
        NW(iPF_V1) = 0.99_DP
        NW(iPF_V2) = 0.0_DP
        NW(iPF_V3) = 0.0_DP
        NW(iPF_E ) = 1.0_DP / ( Gamma_IDEAL - One )

        SW(iPF_D ) = 0.5_DP
        SW(iPF_V1) = 0.0_DP
        SW(iPF_V2) = 0.0_DP
        SW(iPF_V3) = 0.0_DP
        SW(iPF_E ) = 1.0_DP / ( Gamma_IDEAL - One )

        SE(iPF_D ) = 0.1_DP
        SE(iPF_V1) = 0.0_DP
        SE(iPF_V2) = 0.99_DP
        SE(iPF_V3) = 0.0_DP
        SE(iPF_E ) = 1.0_DP / ( Gamma_IDEAL - One )

      CASE( 'IsolatedShock' )

        Vs = 0.01_DP

        X1D = 0.5_DP
        X2D = 0.5_DP

        NE(iPF_D ) = 1.0_DP
        NE(iPF_V1) = -0.9_DP
        NE(iPF_V2) = 0.0_DP
        NE(iPF_V3) = 0.0_DP
        NE(iPF_E ) = 1.0_DP / ( Gamma_IDEAL - One )

        CALL ComputeLeftState &
               ( Vs, &
                 NE(iPF_D ), &
                 NE(iPF_V1), &
                 NE(iPF_E ) * ( Gamma_IDEAL - One ), &
                 NW(iPF_D ), &
                 NW(iPF_V1), &
                 NW(iPF_E ) )

        NW(iPF_V2) = 0.0_DP
        NW(iPF_V3) = 0.0_DP

        SE(iPF_D ) = 1.0_DP
        SE(iPF_V1) = -0.9_DP
        SE(iPF_V2) = 0.0_DP
        SE(iPF_V3) = 0.0_DP
        SE(iPF_E ) = 1.0_DP / ( Gamma_IDEAL - One )

        CALL ComputeLeftState &
               ( Vs, &
                 SE(iPF_D ), &
                 SE(iPF_V1), &
                 SE(iPF_E ) * ( Gamma_IDEAL - One ), &
                 SW(iPF_D ), &
                 SW(iPF_V1), &
                 SW(iPF_E ) )

        SW(iPF_V2) = 0.0_DP
        SW(iPF_V3) = 0.0_DP

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,'(A,A)') &
          'Invalid choice for RiemannProblemName: ', TRIM( RiemannProblemName )
        WRITE(*,'(A)') 'Valid choices:'
        WRITE(*,'(A)') &
          "  'DzB2002' - &
          Blast wave from Del-Zanna & Bucciantini (2002)"
        WRITE(*,'(A)') &
          "  'IsolatedShock'"
        WRITE(*,'(A)') 'Stopping...'
        STOP

    END SELECT

    IF( TRIM( RiemannProblemName ) .EQ. 'IsolatedShock' )THEN

      WRITE(*,'(6x,A,ES14.6E3)') 'Shock Velocity = ', Vs
      WRITE(*,*)

    END IF

    WRITE(*,'(6x,A,F8.6)') 'Gamma_IDEAL = ', Gamma_IDEAL
    WRITE(*,*)
    WRITE(*,'(6x,A,F8.6)') 'X1D = ', X1D
    WRITE(*,'(6x,A,F8.6)') 'X2D = ', X2D
    WRITE(*,*)
    WRITE(*,'(6x,A)') 'NE:'
    WRITE(*,*)
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_D  = ', NE(iPF_D )
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_V1 = ', NE(iPF_V1)
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_V2 = ', NE(iPF_V2)
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_V3 = ', NE(iPF_V3)
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_E  = ', NE(iPF_E )
    WRITE(*,*)
    WRITE(*,'(6x,A)') 'NW:'
    WRITE(*,*)
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_D  = ', NW(iPF_D )
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_V1 = ', NW(iPF_V1)
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_V2 = ', NW(iPF_V2)
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_V3 = ', NW(iPF_V3)
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_E  = ', NW(iPF_E )
    WRITE(*,*)
    WRITE(*,'(6x,A)') 'SE:'
    WRITE(*,*)
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_D  = ', SE(iPF_D )
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_V1 = ', SE(iPF_V1)
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_V2 = ', SE(iPF_V2)
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_V3 = ', SE(iPF_V3)
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_E  = ', SE(iPF_E )
    WRITE(*,*)
    WRITE(*,'(6x,A)') 'SW:'
    WRITE(*,*)
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_D  = ', SW(iPF_D )
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_V1 = ', SW(iPF_V1)
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_V2 = ', SW(iPF_V2)
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_V3 = ', SW(iPF_V3)
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_E  = ', SW(iPF_E )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)
        iNodeX2 = NodeNumberTableX(2,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
        X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

        ! --- NE ---
        IF     ( X1 .GT. X1D .AND. X2 .GT. X2D )THEN

          uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = NE(iPF_D )
          uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = NE(iPF_V1)
          uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = NE(iPF_V2)
          uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = NE(iPF_V3)
          uPF(iNodeX,iX1,iX2,iX3,iPF_E ) = NE(iPF_E )

        ! --- NW ---
        ELSE IF( X1 .LE. X1D .AND. X2 .GT. X2D )THEN

          uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = NW(iPF_D )
          uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = NW(iPF_V1)
          uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = NW(iPF_V2)
          uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = NW(iPF_V3)
          uPF(iNodeX,iX1,iX2,iX3,iPF_E ) = NW(iPF_E )

        ! --- SW ---
        ELSE IF( X1 .LE. X1D .AND. X2 .LE. X2D )THEN

          uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = SW(iPF_D )
          uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = SW(iPF_V1)
          uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = SW(iPF_V2)
          uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = SW(iPF_V3)
          uPF(iNodeX,iX1,iX2,iX3,iPF_E ) = SW(iPF_E )

        ! --- SE ---
        ELSE

          uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = SE(iPF_D )
          uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = SE(iPF_V1)
          uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = SE(iPF_V2)
          uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = SE(iPF_V3)
          uPF(iNodeX,iX1,iX2,iX3,iPF_E ) = SE(iPF_E )

        END IF

      END DO

      CALL ComputePressureFromPrimitive_IDEAL &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_E ), &
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


  END SUBROUTINE InitializeFields_RiemannProblem2D


  SUBROUTINE InitializeFields_RiemannProblemSpherical( RiemannProblemName )

    CHARACTER(LEN=*), INTENT(in) :: RiemannProblemName

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX, iNodeX1
    REAL(DP) :: X1, XD

    REAL(DP) :: LeftState(nPF), RightState(nPF)

    WRITE(*,*)
    WRITE(*,'(A4,A,A)') &
      '', 'Riemann Problem Name: ', &
        TRIM( RiemannProblemName )
    WRITE(*,*)

    SELECT CASE( TRIM( RiemannProblemName ) )

      CASE( 'SphericalSod' )

        XD = 1.0_DP

        LeftState(iPF_D ) = 1.0_DP
        LeftState(iPF_V1) = 0.0_DP
        LeftState(iPF_V2) = 0.0_DP
        LeftState(iPF_V3) = 0.0_DP
        LeftState(iPF_E ) = 1.0_DP / ( Gamma_IDEAL - One )

        RightState(iPF_D ) = 0.125_DP
        RightState(iPF_V1) = 0.0_DP
        RightState(iPF_V2) = 0.0_DP
        RightState(iPF_V3) = 0.0_DP
        RightState(iPF_E ) = 0.1_DP / ( Gamma_IDEAL - One )

      CASE( 'ConstantSpherical' )

        XD = 1.0_DP

        LeftState(iPF_D ) = 1.0_DP
        LeftState(iPF_V1) = 0.0_DP
        LeftState(iPF_V2) = 0.0_DP
        LeftState(iPF_V3) = 0.0_DP
        LeftState(iPF_E ) = 1.0_DP / ( Gamma_IDEAL - One )

        RightState(iPF_D ) = LeftState(iPF_D )
        RightState(iPF_V1) = LeftState(iPF_V1)
        RightState(iPF_V2) = LeftState(iPF_V2)
        RightState(iPF_V3) = LeftState(iPF_V3)
        RightState(iPF_E ) = LeftState(iPF_E )

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,'(A,A)') &
          'Invalid choice for RiemannProblemName: ', RiemannProblemName
        WRITE(*,'(A)') 'Valid choices:'
        WRITE(*,'(A)') &
          "  'SphericalSod' - &
          Spherical Sod's shock tube"
        WRITE(*,'(A)') &
          "  'ConstantSpherical' - &
          All fields constant"
        WRITE(*,'(A)') 'Stopping...'
        STOP

    END SELECT

    WRITE(*,'(6x,A,F8.6)') 'Gamma_IDEAL = ', Gamma_IDEAL
    WRITE(*,*)
    WRITE(*,'(6x,A,F8.6)') 'XD = ', XD
    WRITE(*,*)
    WRITE(*,'(6x,A)') 'Right State:'
    WRITE(*,*)
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_D  = ', RightState(iPF_D )
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_V1 = ', RightState(iPF_V1)
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_V2 = ', RightState(iPF_V2)
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_V3 = ', RightState(iPF_V3)
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_E  = ', RightState(iPF_E )
    WRITE(*,*)
    WRITE(*,'(6x,A)') 'Left State:'
    WRITE(*,*)
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_D  = ', LeftState(iPF_D )
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_V1 = ', LeftState(iPF_V1)
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_V2 = ', LeftState(iPF_V2)
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_V3 = ', LeftState(iPF_V3)
    WRITE(*,'(8x,A,ES14.6E3)') 'PF_E  = ', LeftState(iPF_E )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

        SELECT CASE ( TRIM( RiemannProblemName ) )

          CASE( 'SphericalSod' )

            IF( X1 .LE. XD )THEN

              uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = LeftState(iPF_D )
              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = LeftState(iPF_V1)
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = LeftState(iPF_V2)
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = LeftState(iPF_V3)
              uPF(iNodeX,iX1,iX2,iX3,iPF_E)  = LeftState(iPF_E )

            ELSE

              uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = RightState(iPF_D )
              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = RightState(iPF_V1)
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = RightState(iPF_V2)
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = RightState(iPF_V3)
              uPF(iNodeX,iX1,iX2,iX3,iPF_E ) = RightState(iPF_E )

            END IF

          CASE( 'ConstantSpherical' )

            uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = LeftState(iPF_D )
            uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = LeftState(iPF_V1)
            uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = LeftState(iPF_V2)
            uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = LeftState(iPF_V3)
            uPF(iNodeX,iX1,iX2,iX3,iPF_E)  = LeftState(iPF_E )

        END SELECT

      END DO

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

  END SUBROUTINE InitializeFields_RiemannProblemSpherical


  SUBROUTINE InitializeFields_SedovTaylorBlastWave( nDetCells, Eblast )

    INTEGER,  INTENT(in) :: nDetCells
    REAL(DP), INTENT(in) :: Eblast

    INTEGER  :: iX1, iX2, iX3, iNodeX1, iNodeX
    REAL(DP) :: X1, X_D

    X_D = DBLE( nDetCells ) * MeshX(1) % Width(1)
    WRITE(*,*)
    WRITE(*,'(A,I4.4)')      '     nDetCells:              ', nDetCells
    WRITE(*,'(A,ES23.16E3)') '     Initial blast radius:   ', X_D
    WRITE(*,'(A,ES23.16E3)') '     Blast energy:           ', Eblast
    WRITE(*,'(A,ES23.16E3)') '     Initial blast pressure: ', &
                                     ( Gamma_IDEAL - One ) &
                                       * Eblast / ( FourPi / Three * X_D**3 )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

        IF( X1 <= X_D)THEN

          uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 1.0_DP
          uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
          uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
          uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
          uPF(iNodeX,iX1,iX2,iX3,iPF_E)  &
            = Eblast / ( FourPi / Three * X_D**3 )
          uAF(iNodeX,iX1,iX2,iX3,iAF_P)  &
            = ( Gamma_IDEAL - One ) * uPF(iNodeX,iX1,iX2,iX3,iPF_E)

        ELSE

          uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 1.0_DP
          uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
          uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
          uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
          uPF(iNodeX,iX1,iX2,iX3,iPF_E)  &
            = 1.0d-5
          uAF(iNodeX,iX1,iX2,iX3,iAF_P)  &
            = ( Gamma_IDEAL - One ) * uPF(iNodeX,iX1,iX2,iX3,iPF_E)

        END IF

      END DO

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

  END SUBROUTINE InitializeFields_SedovTaylorBlastWave


  ! --- Relativistic Kelvin-Helmholtz instability a la
  !     Beckwith & Stone (2011), ApjS, 193, 6 (typo in Eq. (63)) ---
  SUBROUTINE InitializeFields_KelvinHelmholtzInstability

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX, iNodeX1, iNodeX2
    REAL(DP) :: X1, X2
    REAL(DP) :: rho0, rho1
    REAL(DP) :: Vshear, a, X2_Offset, sigma, A0, Vz

    INTEGER, ALLOCATABLE :: Seed(:)
    INTEGER              :: SeedNum = 1

    CALL RANDOM_SEED( SIZE = SeedNum )
    ALLOCATE( Seed(SeedNum) )
    Seed(1) = 1
    CALL RANDOM_SEED( PUT = Seed )

    rho0 = 0.505_DP
    rho1 = 0.495_DP

    Vshear    = 0.5_DP
    a         = 0.01_DP
    X2_Offset = 0.5_DP
    sigma     = 0.1_DP

    A0 = 0.1_DP

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)
        iNodeX2 = NodeNumberTableX(2,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
        X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

        ! --- Top ---

        IF( X2 .GT. Zero )THEN

          uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
            = rho0 + rho1 * TANH( ( X2 - X2_Offset ) / a )
          uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
            = Vshear      * TANH( ( X2 - X2_Offset ) / a )

          ! --- This is where the typo is. The following expression is
          !     taken from Radice & Rezzolla, 2012, AA, 547, A26, Eq. (48) ---

          uPF(iNodeX,iX1,iX2,iX3,iPF_V2) &
            = A0 * Vshear * SIN( TwoPi * X1 ) &
                * EXP( -( ( X2 - X2_Offset ) / sigma )**2 )

        ! --- Bottom ---

        ELSE

          uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
            = rho0 - rho1 * TANH( ( X2 + X2_Offset ) / a )
          uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
            = -Vshear     * TANH( ( X2 + X2_Offset ) / a )
          uPF(iNodeX,iX1,iX2,iX3,iPF_V2) &
            = -A0 * Vshear * SIN( TwoPi * X1 ) &
                * EXP( -( ( X2 + X2_Offset ) / sigma )**2 )

        END IF

        IF( nDimsX .EQ. 2 )THEN

          uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP

        ELSE

          CALL RANDOM_NUMBER( Vz )
          uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.01_DP * Vz

        END IF

        uAF(iNodeX,iX1,iX2,iX3,iAF_P) = One
        uPF(iNodeX,iX1,iX2,iX3,iPF_E) &
          = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - One )

      END DO

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

    DEALLOCATE( Seed )

  END SUBROUTINE InitializeFields_KelvinHelmholtzInstability


  ! --- Auxiliary functions/subroutines for computing left state ---


  SUBROUTINE ComputeLeftState( Vs, DR, VR, PR, DL, VL, PL )

    REAL(DP), INTENT(in)  :: Vs, DR, VR, PR
    REAL(DP), INTENT(out) ::     DL, VL, PL

    CALL ApplyJumpConditions_LeftState( Vs, DR, VR, PR, DL, VL, PL )

    ! --- Return energy-density instead of pressure ---
    PL = PL / ( Gamma_IDEAL - One )

  END SUBROUTINE ComputeLeftState


  SUBROUTINE ApplyJumpConditions_LeftState( Vs, DR, VR, PR, DL, VL, PL )

    REAL(DP), INTENT(in)  :: Vs, DR, VR, PR
    REAL(DP), INTENT(out) ::     DL, VL, PL

    REAL(DP), PARAMETER :: EPS = 1.0e-15_DP

    REAL(DP), PARAMETER :: ToldV = EPS
    REAL(DP), PARAMETER :: TolF  = EPS
    INTEGER,  PARAMETER :: nMaxIter = 1000

    INTEGER :: ITERATION
    REAL(DP) :: D, V, P, F
    REAL(DP) :: Vmin, Vmax, Fmin, Fmax, VV, FF

    IF( VR .LT. Zero )THEN

      Vmin = VR   + EPS
      Vmax = +One - EPS

    ELSE

      Vmin = -One + EPS
      Vmax = VR   - EPS

    END IF

    D = Density ( Vs, DR, VR, Vmin )
    P = Pressure( Vs, DR, VR, PR, D, Vmin )
    Fmin = PostShockVelocity( Vs, DR, VR, PR, D, Vmin, P )

    D = Density( Vs, DR, VR, Vmax )
    P = Pressure( Vs, DR, VR, PR, D, Vmax )
    Fmax = PostShockVelocity( Vs, DR, VR, PR, D, Vmax, P )

    IF( .NOT. Fmin * Fmax .LT. Zero )THEN

      WRITE(*,*) 'Root not bracketed. Stopping...'
      WRITE(*,*) 'Fmin = ', Fmin
      WRITE(*,*) 'Fmax = ', Fmax
      STOP

    END IF

    IF( Fmin .GT. Zero )THEN

      VV = Vmax
      FF = Fmax

      Vmax = Vmin
      Vmin = VV

      Fmax = Fmin
      Fmin = FF

    END IF

    ITERATION = 0
    DO WHILE( ITERATION .LT. nMaxIter )

      ITERATION = ITERATION + 1

      V = ( Vmin + Vmax ) / Two

      D = Density ( Vs, DR, VR, V )
      P = Pressure( Vs, DR, VR, PR, D, V )

      F = PostShockVelocity( Vs, DR, VR, PR, D, V, P )

      IF( ABS( V - Vmin ) / MAX( ABS( Vmax ), ABS( Vmin ) ) .LT. ToldV ) EXIT

      IF( F .GT. Zero )THEN

        Vmax = V
        Fmax = F

     ELSE

        Vmin = V
        Fmin = F

     END IF

    END DO

!!$    WRITE(*,*) 'Converged at iteration ', ITERATION
!!$    WRITE(*,*) '|F|:  ' , ABS( F )
!!$    WRITE(*,*) 'dV/V: ', ABS( V - Vmax ) / ABS( Vmax )

    VL = V
    DL = Density ( Vs, DR, VR, VL )
    PL = Pressure( Vs, DR, VR, PR, DL, VL )

  END SUBROUTINE ApplyJumpConditions_LeftState


  REAL(DP) FUNCTION Density( Vs, DR, VR, VL )

    REAL(DP), INTENT(in) :: Vs, DR, VR, VL

    REAL(DP) :: WR, WL

    WR = LorentzFactor( One, VR )
    WL = LorentzFactor( One, VL )

    Density = DR * ( WR * ( VR - Vs ) ) / ( WL * ( VL - Vs ) )

    RETURN
  END FUNCTION Density


  REAL(DP) FUNCTION Pressure( Vs, DR, VR, PR, DL, VL )

    REAL(DP), INTENT(in) :: Vs, DR, VR, PR, DL, VL

    REAL(DP) :: WR, WL, tau

    WR = LorentzFactor( One, VR )
    WL = LorentzFactor( One, VL )

    tau = Gamma_IDEAL / ( Gamma_IDEAL - One )

    Pressure = ( PR * ( One + tau * WR**2 * VR * ( VR - Vs ) ) &
                 - DL * WL**2 * VL**2 + DR * WR**2 * VR**2 &
                 + Vs * ( DL * WL**2 * VL - DR * WR**2 * VR ) ) &
               / ( One + tau * WL**2 * VL * ( VL - Vs ) )

    RETURN
  END FUNCTION Pressure


  REAL(DP) FUNCTION PostShockVelocity( Vs, DR, VR, PR, DL, VL, PL )

    REAL(DP), INTENT(in) :: Vs, DR, VR, PR, DL, VL, PL

    REAL(DP) :: WR, WL, tau

    WR = LorentzFactor( One, VR )
    WL = LorentzFactor( One, VL )

    tau = Gamma_IDEAL / ( Gamma_IDEAL - One )

    PostShockVelocity &
      = ( DL + tau * PL ) * WL**2 * ( VL - Vs ) &
          - ( DR + tau * PR ) * WR**2 * ( VR - Vs ) + Vs * ( PL - PR )

    RETURN
  END FUNCTION PostShockVelocity


  REAL(DP) FUNCTION LorentzFactor( Psi, V )

    REAL(DP), INTENT(in) :: Psi, V

    LorentzFactor = One / SQRT( One - Psi**4 * ( V / SpeedOfLight )**2 )

    RETURN
  END FUNCTION LorentzFactor

  ! --- End of auxiliary functions/subroutines for computing left state ---


END MODULE InitializationModule_Relativistic
