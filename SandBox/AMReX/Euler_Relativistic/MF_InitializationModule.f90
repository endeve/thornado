MODULE MF_InitializationModule

  ! --- AMReX Modules ---

  USE amrex_fort_module,      ONLY: &
    amrex_real
  USE amrex_box_module,       ONLY: &
    amrex_box
  USE amrex_multifab_module,  ONLY: &
    amrex_multifab, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy
  USE amrex_parallel_module,  ONLY: &
    amrex_parallel_ioprocessor
  USE amrex_parmparse_module, ONLY: &
    amrex_parmparse, &
    amrex_parmparse_build, &
    amrex_parmparse_destroy

  ! --- thornado Modules ---

  USE KindModule,              ONLY: &
    DP, Half, One, Pi, TwoPi, FourPi
  USE ProgramHeaderModule,     ONLY: &
    nDOFX, nX, nNodesX, swX
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE MeshModule,              ONLY: &
    MeshType,                        &
    CreateMesh,                      &
    DestroyMesh,                     &
    NodeCoordinate
  USE GeometryFieldsModule,    ONLY: &
    nGF, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33
  USE FluidFieldsModule,       ONLY: &
    nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iPF_Ne, &
    nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iCF_Ne, &
    nAF, iAF_P
  USE Euler_UtilitiesModule,   ONLY: &
    Euler_ComputeConserved
  USE EquationOfStateModule,   ONLY: &
    ComputePressureFromPrimitive

  ! --- Local Modules ---
  USE MyAmrModule, ONLY: &
    nLevels, xL, xR, Gamma_IDEAL

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MF_InitializeFields


CONTAINS


  SUBROUTINE MF_InitializeFields( ProgramName, MF_uGF, MF_uCF )

    CHARACTER(LEN=*),     INTENT(in   ) :: ProgramName
    TYPE(amrex_multifab), INTENT(in   ) :: MF_uGF(0:nLevels)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels)

    IF( amrex_parallel_ioprocessor() )THEN
      WRITE(*,*)
      WRITE(*,'(A4,A,A)') '', 'Initializing: ', TRIM( ProgramName )
    END IF

    SELECT CASE ( TRIM( ProgramName ) )

      CASE( 'Sod_Relativistic' )

        CALL InitializeFields_Sod_Relativistic( MF_uGF, MF_uCF )

      CASE( 'KelvinHelmholtz_Relativistic' )

        CALL InitializeFields_KelvinHelmholtz_Relativistic( MF_uGF, MF_uCF )

      CASE( 'RiemannProblem_2D' )

        CALL InitializeFields_RiemannProblem_2D( MF_uGF, MF_uCF )

      CASE DEFAULT

        IF( amrex_parallel_ioprocessor() )THEN
          WRITE(*,*)
          WRITE(*,'(A4,A,A)') '', 'Unknown Program: ', TRIM( ProgramName )
        END IF

    END SELECT

  END SUBROUTINE MF_InitializeFields


  SUBROUTINE InitializeFields_Sod_Relativistic( MF_uGF, MF_uCF )

    TYPE(amrex_multifab), INTENT(in   ) :: MF_uGF(0:nLevels)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels)

    INTEGER            :: iDim, iLevel
    INTEGER            :: iX1, iX2, iX3
    INTEGER            :: iNodeX
    INTEGER            :: lo_G(4), hi_G(4)
    INTEGER            :: lo_F(4), hi_F(4)
    REAL(amrex_real)   :: X1
    REAL(amrex_real)   :: uGF_K(nDOFX,nGF)
    REAL(amrex_real)   :: uCF_K(nDOFX,nCF)
    REAL(amrex_real)   :: uPF_K(nDOFX,nPF)
    REAL(amrex_real)   :: uAF_K(nDOFX,nAF)
    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI
    TYPE(MeshType)     :: MeshX(3)
    REAL(amrex_real), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(amrex_real), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    uGF_K = 0.0d0
    uPF_K = 0.0d0
    uCF_K = 0.0d0

    DO iDim = 1, 3

      CALL CreateMesh &
             ( MeshX(iDim), nX(iDim), nNodesX(iDim), 0, &
               xL(iDim), xR(iDim) )

    END DO

    DO iLevel = 0, nLevels

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = .TRUE. )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCF => MF_uCF(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        lo_G = LBOUND( uGF )
        hi_G = UBOUND( uGF )

        lo_F = LBOUND( uCF )
        hi_F = UBOUND( uCF )

        DO iX3 = BX % lo(3), BX % hi(3)
        DO iX2 = BX % lo(2), BX % hi(2)
        DO iX1 = BX % lo(1), BX % hi(1)

          uGF_K &
            = RESHAPE( uGF(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )

          X1 = MeshX(1) % Center(iX1)

          DO iNodeX = 1, nDOFX

            IF( X1 .LE. 0.5_amrex_real ) THEN

              uPF_K(iNodeX,iPF_D)  = 1.0_amrex_real
              uPF_K(iNodeX,iPF_V1) = 0.0_amrex_real
              uPF_K(iNodeX,iPF_V2) = 0.0_amrex_real
              uPF_K(iNodeX,iPF_V3) = 0.0_amrex_real
              uPF_K(iNodeX,iPF_E)  = 1.0_amrex_real &
                                       / (Gamma_IDEAL - 1.0_amrex_real)

            ELSE

              uPF_K(iNodeX,iPF_D)  = 0.125_amrex_real
              uPF_K(iNodeX,iPF_V1) = 0.0_amrex_real
              uPF_K(iNodeX,iPF_V2) = 0.0_amrex_real
              uPF_K(iNodeX,iPF_V3) = 0.0_amrex_real
              uPF_K(iNodeX,iPF_E)  = 0.1_amrex_real &
                                       / (Gamma_IDEAL - 1.0_amrex_real)

            END IF

          END DO

          CALL ComputePressureFromPrimitive &
                 ( uPF_K(:,iPF_D), uPF_K(:,iPF_E), uPF_K(:,iPF_Ne), &
                   uAF_K(:,iAF_P) )

          CALL Euler_ComputeConserved &
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

      CALL amrex_mfiter_destroy( MFI )

    END DO

    DO iDim = 1, 3

      CALL DestroyMesh( MeshX(iDim) )

    END DO

  END SUBROUTINE InitializeFields_Sod_Relativistic


  ! --- Relativistic 2D Kelvin-Helmholtz instability a la
  !     Beckwith & Stone (2011), ApjS, 193, 6 (typo in Eq. (63)) ---
  SUBROUTINE InitializeFields_KelvinHelmholtz_Relativistic( MF_uGF, MF_uCF )

    TYPE(amrex_multifab), INTENT(in   ) :: MF_uGF(0:nLevels)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels)

    INTEGER          :: iX1, iX2, iX3
    INTEGER          :: iNodeX, iNodeX1, iNodeX2
    REAL(amrex_real) :: X1, X2
    REAL(amrex_real) :: rho0, rho1
    REAL(amrex_real) :: Vshear, a, X2_Offset, sigma, A0

    INTEGER            :: iDim, iLevel
    INTEGER            :: lo_G(4), hi_G(4)
    INTEGER            :: lo_F(4), hi_F(4)
    REAL(amrex_real)   :: uGF_K(nDOFX,nGF)
    REAL(amrex_real)   :: uCF_K(nDOFX,nCF)
    REAL(amrex_real)   :: uPF_K(nDOFX,nPF)
    REAL(amrex_real)   :: uAF_K(nDOFX,nAF)
    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI
    TYPE(MeshType)     :: MeshX(3)
    REAL(amrex_real), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(amrex_real), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    rho0 = 0.505d0
    rho1 = 0.495d0

    Vshear    = 0.5d0
    a         = 0.01d0
    X2_Offset = 0.5d0
    sigma     = 0.1d0

    A0 = 0.1d0

    uGF_K = 0.0d0
    uPF_K = 0.0d0
    uCF_K = 0.0d0

    DO iDim = 1, 3

      CALL CreateMesh &
             ( MeshX(iDim), nX(iDim), nNodesX(iDim), 0, &
               xL(iDim), xR(iDim) )

    END DO

    DO iLevel = 0, nLevels

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = .TRUE. )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCF => MF_uCF(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        lo_G = LBOUND( uGF )
        hi_G = UBOUND( uGF )

        lo_F = LBOUND( uCF )
        hi_F = UBOUND( uCF )

        DO iX3 = BX % lo(3), BX % hi(3)
        DO iX2 = BX % lo(2), BX % hi(2)
        DO iX1 = BX % lo(1), BX % hi(1)

          uGF_K &
            = RESHAPE( uGF(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )

          DO iNodeX = 1, nDOFX

            iNodeX1 = NodeNumberTableX(1,iNodeX)
            iNodeX2 = NodeNumberTableX(2,iNodeX)

            X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
            X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

            ! --- Top ---
            IF( X2 .GT. 0.0d0 )THEN
              uPF_K(iNodeX,iPF_D) &
                = rho0 + rho1 * TANH( ( X2 - X2_Offset ) / a )
              uPF_K(iNodeX,iPF_V1) &
                = Vshear      * TANH( ( X2 - X2_Offset ) / a )

              ! --- This is where the typo is. The following expression is
              !     taken from Radice & Rezzolla, 2012, AA, 547, A26, Eq. (48)
              uPF_K(iNodeX,iPF_V2) &
                = A0 * Vshear * SIN( 2.0d0 * Pi * X1 ) &
                    * EXP( -( ( X2 - X2_Offset ) / sigma )**2 )

            ! --- Bottom ---
            ELSE
              uPF_K(iNodeX,iPF_D) &
                = rho0 - rho1 * TANH( ( X2 + X2_Offset ) / a )
              uPF_K(iNodeX,iPF_V1) &
                = -Vshear     * TANH( ( X2 + X2_Offset ) / a )
              uPF_K(iNodeX,iPF_V2) &
                = -A0 * Vshear * SIN( 2.0d0 * Pi * X1 ) &
                    * EXP( -( ( X2 + X2_Offset ) / sigma )**2 )

             END IF

            uPF_K(iNodeX,iPF_V3) = 0.0d0
            uPF_K(iNodeX,iPF_E)  = 1.0d0 / ( Gamma_IDEAL - One )

          END DO

          CALL ComputePressureFromPrimitive &
                 ( uPF_K(:,iPF_D), uPF_K(:,iPF_E), uPF_K(:,iPF_Ne), &
                   uAF_K(:,iAF_P) )

          CALL Euler_ComputeConserved &
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

      CALL amrex_mfiter_destroy( MFI )

    END DO

    DO iDim = 1, 3

      CALL DestroyMesh( MeshX(iDim) )

    END DO

  END SUBROUTINE InitializeFields_KelvinHelmholtz_Relativistic


  ! --- Relativistic 2D Riemann problem from
  !     Del Zanna & Bucciantini, (2002), A&A, 390, 1177 ---
  SUBROUTINE InitializeFields_RiemannProblem_2D( MF_uGF, MF_uCF )

    TYPE(amrex_multifab), INTENT(in   ) :: MF_uGF(0:nLevels)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels)

    INTEGER          :: iX1, iX2, iX3
    INTEGER          :: iNodeX, iNodeX1, iNodeX2
    REAL(amrex_real) :: X1, X2

    INTEGER            :: iDim, iLevel
    INTEGER            :: lo_G(4), hi_G(4)
    INTEGER            :: lo_F(4), hi_F(4)
    REAL(amrex_real)   :: uGF_K(nDOFX,nGF)
    REAL(amrex_real)   :: uCF_K(nDOFX,nCF)
    REAL(amrex_real)   :: uPF_K(nDOFX,nPF)
    REAL(amrex_real)   :: uAF_K(nDOFX,nAF)
    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI
    TYPE(MeshType)     :: MeshX(3)
    REAL(amrex_real), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(amrex_real), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    uGF_K = 0.0d0
    uPF_K = 0.0d0
    uCF_K = 0.0d0

    DO iDim = 1, 3

      CALL CreateMesh &
             ( MeshX(iDim), nX(iDim), nNodesX(iDim), 0, &
               xL(iDim), xR(iDim) )

    END DO

    DO iLevel = 0, nLevels

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = .TRUE. )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCF => MF_uCF(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        lo_G = LBOUND( uGF )
        hi_G = UBOUND( uGF )

        lo_F = LBOUND( uCF )
        hi_F = UBOUND( uCF )

        DO iX3 = BX % lo(3), BX % hi(3)
        DO iX2 = BX % lo(2), BX % hi(2)
        DO iX1 = BX % lo(1), BX % hi(1)

          uGF_K &
            = RESHAPE( uGF(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )

          DO iNodeX = 1, nDOFX

            iNodeX1 = NodeNumberTableX(1,iNodeX)
            iNodeX2 = NodeNumberTableX(2,iNodeX)

            X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
            X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

            ! --- NE ---
            IF     ( X1 .GT. 0.5_amrex_real .AND. X2 .GT. 0.5_amrex_real )THEN
              uPF_K(iNodeX,iPF_D ) = 0.1_amrex_real
              uPF_K(iNodeX,iPF_V1) = 0.0_amrex_real
              uPF_K(iNodeX,iPF_V2) = 0.0_amrex_real
              uPF_K(iNodeX,iPF_E ) = 0.01_amrex_real &
                                       / ( Gamma_IDEAL - 1.0_amrex_real )
            ! --- NW ---
            ELSE IF( X1 .LE. 0.5_amrex_real .AND. X2 .GT. 0.5_amrex_real )THEN
              uPF_K(iNodeX,iPF_D ) = 0.1_amrex_real
              uPF_K(iNodeX,iPF_V1) = 0.99_amrex_real
              uPF_K(iNodeX,iPF_V2) = 0.0_amrex_real
              uPF_K(iNodeX,iPF_E ) = 1.0_amrex_real &
                                       / ( Gamma_IDEAL - 1.0_amrex_real )
            ! --- SW ---
            ELSE IF( X1 .LE. 0.5_amrex_real .AND. X2 .LE. 0.5_amrex_real )THEN
              uPF_K(iNodeX,iPF_D ) = 0.5_amrex_real
              uPF_K(iNodeX,iPF_V1) = 0.0_amrex_real
              uPF_K(iNodeX,iPF_V2) = 0.0_amrex_real
              uPF_K(iNodeX,iPF_E ) = 1.0_amrex_real &
                                       / ( Gamma_IDEAL - 1.0_amrex_real )

            ! --- SE ---
            ELSE
              uPF_K(iNodeX,iPF_D ) = 0.1_amrex_real
              uPF_K(iNodeX,iPF_V1) = 0.0_amrex_real
              uPF_K(iNodeX,iPF_V2) = 0.99_amrex_real
              uPF_K(iNodeX,iPF_E ) = 1.0_amrex_real &
                                       / ( Gamma_IDEAL - 1.0_amrex_real )
            END IF

          END DO

          uPF_K(:,iPF_V3) = 0.0_amrex_real

          CALL ComputePressureFromPrimitive &
                 ( uPF_K(:,iPF_D), uPF_K(:,iPF_E), uPF_K(:,iPF_Ne), &
                   uAF_K(:,iAF_P) )

          CALL Euler_ComputeConserved &
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

      CALL amrex_mfiter_destroy( MFI )

    END DO

    DO iDim = 1, 3

      CALL DestroyMesh( MeshX(iDim) )

    END DO

  END SUBROUTINE InitializeFields_RiemannProblem_2D


END MODULE MF_InitializationModule
