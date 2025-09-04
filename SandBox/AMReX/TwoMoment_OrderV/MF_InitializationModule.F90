MODULE MF_InitializationModule

  ! --- AMReX Modules ---

  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor
  USE amrex_parmparse_module, ONLY: &
    amrex_parmparse, &
    amrex_parmparse_build, &
    amrex_parmparse_destroy
  USE amrex_geometry_module, ONLY: &
    amrex_geometry
  USE amrex_amrcore_module, ONLY: &
    amrex_geom

  ! --- thornado Modules ---

  USE EquationOfStateModule, ONLY: &
    ComputePressureFromPrimitive

  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    nDOFX, &
    nX, &
    nE, &
    swX, &
    nNodesX, &
    nDOFE, &
    nDOFZ, &
    iZ_B0, &
    iZ_E0
  USE RadiationFieldsModule, ONLY: &
    iCR_N, &
    uPR, &
    uCR, &
    iCR_G1, &
    iCR_G2, &
    iCR_G3, &
    nCR, &
    iPR_D, &
    iPR_I1, &
    iPR_I2, &
    iPR_I3, &
    nPR, &
    nCR, &
    nSpecies
  USE FluidFieldsModule, ONLY: &
    uPF, &
    nCF, &
    iCF_D, &
    iCF_S1, &
    iCF_S2, &
    iCF_S3, &
    iCF_E, &
    iCF_Ne, &
    nPF, &
    iPF_D, &
    iPF_V1, &
    iPF_V2, &
    iPF_V3, &
    iPF_E, &
    iPF_Ne, &
    nAF, &
    iAF_P, &
    iAF_T, &
    iAF_Ye, &
    iAF_S, &
    iAF_E, &
    iAF_Gm, &
    iAF_Me, &
    iAF_Mp, &
    iAF_Mn, &
    iAF_Xp, &
    iAF_Xn, &
    iAF_Xa, &
    iAF_Xh
  USE Euler_UtilitiesModule_NonRelativistic, ONLY: &
    ComputeConserved_Euler_NonRelativistic
  USE GeometryFieldsModule, ONLY: &
    uGF, &
    nGF, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_h_1, &
    iGF_h_2, &
    iGF_h_3, &
    iGF_Phi_N, &
    iGF_Psi, &
    iGF_SqrtGm
  USE TwoMoment_OpacityModule, ONLY: &
   uOP, &
   iOP_Sigma
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE MeshModule, ONLY: &
    MeshType, &
    MeshX, &
    MeshE, &
    NodeCoordinate
  USE ReferenceElementModule, ONLY: &
    NodeNumberTable
  USE TwoMoment_UtilitiesModule, ONLY: &
    ComputeConserved_TwoMoment
  USE UnitsModule, ONLY: &
    UnitsDisplay, &
    SolarMass, &
    SpeedOfLight, &
    Erg, &
    Gram, &
    Centimeter, &
    MeV, &
    BoltzmannConstant, &
    GravitationalConstant

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero, &
    One, &
    Pi, &
    Half, &
    TwoPi, &
    Three
  USE MF_ErrorModule, ONLY: &
    DescribeError_MF
  USE InputParsingModule, ONLY: &
    nLevels, &
    UseTiling

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFields_MF

CONTAINS


  SUBROUTINE InitializeFields_MF &
    ( iLevel, MF_uGF, MF_uCR, MF_uCF )

    INTEGER,              INTENT(in)    :: iLevel
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR, MF_uCF

    IF( iLevel .EQ. 0 .AND. amrex_parallel_ioprocessor() )THEN

      WRITE(*,*)
      WRITE(*,'(4x,A,A)') 'INFO: Initial Conditions'
      WRITE(*,'(4x,A,A)') '------------------------'
      WRITE(*,*)

    END IF

    SELECT CASE( TRIM( ProgramName ) )

      CASE( 'SineWaveStreaming' )

        CALL InitializeFields_SineWaveStreaming &
               ( iLevel, MF_uGF, MF_uCR, MF_uCF )

      CASE( 'GaussianDiffusion' )

        CALL InitializeFields_GaussianDiffusion &
               ( iLevel, MF_uGF, MF_uCR, MF_uCF )

      CASE( 'TransparentShock' )

        CALL InitializeFields_TransparentShock &
               ( iLevel, MF_uGF, MF_uCR, MF_uCF )

      CASE( 'StreamingDopplerShift' )

        CALL InitializeFields_StreamingDopplerShift &
               ( iLevel, MF_uGF, MF_uCR, MF_uCF )

      CASE( 'TransparentVortex' )

        CALL InitializeFields_TransparentVortex &
               ( iLevel, MF_uGF, MF_uCR, MF_uCF )

      CASE DEFAULT

        CALL DescribeError_MF &
               ( 301, Message_Option &
                        = 'Invalid ProgramName: ' // TRIM( ProgramName ) )

    END SELECT

  END SUBROUTINE InitializeFields_MF


  ! --- PRIVATE SUBROUTINES ---


  SUBROUTINE InitializeFields_SineWaveStreaming &
    ( iLevel, MF_uGF, MF_uCR, MF_uCF )

    INTEGER             , INTENT(in)    :: iLevel
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF

    ! --- thornado ---

    INTEGER        :: iDim, iE
    INTEGER        :: iX1, iX2, iX3, iZ1, iZ2, iZ3, iZ4, iS, iNodeZ, iSpecies
    INTEGER        :: iNodeX, iNodeX1, iNodeX2, iNodeX3, &
                      iNodeZ2, iNodeZ3, iNodeZ4, iNodeE
    REAL(DP)       :: X1, X2, X3, X_2D, L
    REAL(DP)       :: uCR_K(nDOFZ,nE,nCR,nSpecies)
    REAL(DP)       :: uPR_K(nDOFZ,nE,nPR,nSpecies)
    REAL(DP)       :: uGF_K(nDOFX,nGF)
    REAL(DP)       :: uPF_K(nDOFX,nPF)
    REAL(DP)       :: uCF_K(nDOFX,nCF)
    REAL(DP)       :: uAF_K(nDOFX,nAF)

    ! --- AMReX ---

    INTEGER                       :: lo_C(4), hi_C(4)
    INTEGER                       :: lo_G(4), hi_G(4)
    INTEGER                       :: lo_F(4), hi_F(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCR(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    ! --- Problem-dependent parameters ---

    TYPE(amrex_parmparse) :: PP
    CHARACTER(:), ALLOCATABLE :: Direction
    REAL(DP)    , ALLOCATABLE :: V_0(:)

    Direction = 'X'
    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL PP % query ( 'Direction', &
                         Direction )
      CALL PP % getarr( 'V_0', &
                         V_0 )
    CALL amrex_parmparse_destroy( PP )

    uCR_K = Zero
    uPF_K = Zero
    uCF_K = Zero
    uGF_K = Zero
    uAF_K = Zero

    CALL amrex_mfiter_build( MFI, MF_uCR, tiling = UseTiling )

    DO WHILE( MFI % next() )

      uGF => MF_uGF % DataPtr( MFI )
      uCR => MF_uCR % DataPtr( MFI )
      uCF => MF_uCF % DataPtr( MFI )

      BX = MFI % tilebox()

      lo_G = LBOUND( uGF )
      hi_G = UBOUND( uGF )

      lo_C = LBOUND( uCR )
      hi_C = UBOUND( uCR )

      lo_F = LBOUND( uCF )
      hi_F = UBOUND( uCF )

      DO iX3 = BX % lo(3), BX % hi(3)
      DO iX2 = BX % lo(2), BX % hi(2)
      DO iX1 = BX % lo(1), BX % hi(1)

        uGF_K &
          = RESHAPE( uGF(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )

        ! --- Fluid Fields ---

        DO iNodeX = 1, nDOFX

          uPF_K(iNodeX,iPF_D ) = 1.0_DP
          uPF_K(iNodeX,iPF_V1) = V_0(1)
          uPF_K(iNodeX,iPF_V2) = V_0(2)
          uPF_K(iNodeX,iPF_V3) = V_0(3)
          uPF_K(iNodeX,iPF_E ) = 0.1_DP
          uPF_K(iNodeX,iPF_Ne) = 0.0_DP

        END DO

        CALL ComputeConserved_Euler_NonRelativistic &
              (  uPF_K(:,iPF_D ), &
                 uPF_K(:,iPF_V1), &
                 uPF_K(:,iPF_V2), &
                 uPF_K(:,iPF_V3), &
                 uPF_K(:,iPF_E ), &
                 uPF_K(:,iPF_Ne), &
                 uCF_K(:,iCF_D ), &
                 uCF_K(:,iCF_S1), &
                 uCF_K(:,iCF_S2), &
                 uCF_K(:,iCF_S3), &
                 uCF_K(:,iCF_E ), &
                 uCF_K(:,iCF_Ne), &
                 uGF_K(:,iGF_Gm_dd_11), &
                 uGF_K(:,iGF_Gm_dd_22), &
                 uGF_K(:,iGF_Gm_dd_33))

        uCF(iX1,iX2,iX3,lo_F(4):hi_F(4)) &
          = RESHAPE( uCF_K, [ hi_F(4) - lo_F(4) + 1 ] )

        DO iNodeZ = 1, nDOFZ

          DO iS  = 1, nSpecies
          DO iZ1 = 1, nE

            iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1

            iNodeZ2 = NodeNumberTable(2,iNodeZ)
            iNodeZ3 = NodeNumberTable(3,iNodeZ)
            iNodeZ4 = NodeNumberTable(4,iNodeZ)

            X1 = NodeCoordinate( MeshX(1), iX1, iNodeZ2 )
            X2 = NodeCoordinate( MeshX(2), iX2, iNodeZ3 )
            X3 = NodeCoordinate( MeshX(3), iX3, iNodeZ4 )

            SELECT CASE( TRIM( Direction ) )
              CASE( 'X' )

                uPR_K( iNodeZ, iZ1, iPR_D, iS ) &
                  = 0.50_DP + 0.49_DP * SIN( TwoPi * X1 )

                uPR_K( iNodeZ, iZ1, iPR_I1, iS ) &
                  = uPR_K( iNodeZ, iZ1, iPR_D, iS )

                uPR_K( iNodeZ, iZ1, iPR_I2, iS ) &
                  = 0.0_DP

                uPR_K( iNodeZ, iZ1, iPR_I3, iS ) &
                  = 0.0_DP

              CASE( 'Y' )

                uPR_K( iNodeZ, iZ1, iPR_D, iS ) &
                  = 0.50_DP + 0.49_DP * SIN( TwoPi * X2 )

                uPR_K( iNodeZ, iZ1, iPR_I1, iS ) &
                  = 0.0_DP

                uPR_K( iNodeZ, iZ1, iPR_I2, iS ) &
                  = uPR_K( iNodeZ, iZ1, iPR_D, iS )

                uPR_K( iNodeZ, iZ1, iPR_I3, iS ) &
                  = 0.0_DP

              CASE( 'Z' )

                uPR_K( iNodeZ, iZ1, iPR_D, iS ) &
                  = 0.50_DP + 0.49_DP * SIN( TwoPi * X3 )

                uPR_K( iNodeZ, iZ1, iPR_I1, iS ) &
                  = 0.0_DP

                uPR_K( iNodeZ, iZ1, iPR_I2, iS ) &
                  = 0.0_DP

                uPR_K( iNodeZ, iZ1, iPR_I3, iS ) &
                  = uPR_K( iNodeZ, iZ1, iPR_D, iS )


              CASE( 'XY' )

                X_2D = SQRT( 2.0_DP ) * X1 &
                         + SQRT( 2.0_DP ) * X2

                L  = SQRT( 2.0_DP )

                uPR_K( iNodeZ, iZ1, iPR_D, iS ) &
                  =0.50_DP + 0.49_DP * SIN( TwoPi * X_2D / L )

                uPR_K( iNodeZ, iZ1, iPR_I1, iS ) &
                  = SQRT( 2.0_DP ) / 2.0_DP * uPR_K( iNodeZ, iZ1, iPR_D, iS )

                uPR_K( iNodeZ, iZ1, iPR_I2, iS ) &
                  = SQRT( 2.0_DP ) / 2.0_DP * uPR_K( iNodeZ, iZ1, iPR_D, iS )

                uPR_K( iNodeZ, iZ1, iPR_I3, iS ) &
                  = 0.0_DP

              CASE DEFAULT

                WRITE(*,*)
                WRITE(*,'(A8,A)')    &
                  '', 'InitializeFields_SineWaveStreaming'
                WRITE(*,'(A8,A,A2)') &
                  '', 'Invalid Direction: ', TRIM( Direction )
                WRITE(*,*)
                STOP

          END SELECT

          CALL ComputeConserved_TwoMoment &
                 ( uPR_K(iNodeZ,iZ1,iPR_D ,iS), &
                   uPR_K(iNodeZ,iZ1,iPR_I1,iS), &
                   uPR_K(iNodeZ,iZ1,iPR_I2,iS), &
                   uPR_K(iNodeZ,iZ1,iPR_I3,iS), &
                   uCR_K(iNodeZ,iZ1,iCR_N ,iS), &
                   uCR_K(iNodeZ,iZ1,iCR_G1,iS), &
                   uCR_K(iNodeZ,iZ1,iCR_G2,iS), &
                   uCR_K(iNodeZ,iZ1,iCR_G3,iS), &
                   uPF_K(iNodeX,iPF_V1), &
                   uPF_K(iNodeX,iPF_V2), &
                   uPF_K(iNodeX,iPF_V3), &
                   uGF_K(iNodeX,iGF_Gm_dd_11), &
                   uGF_K(iNodeX,iGF_Gm_dd_22), &
                   uGF_K(iNodeX,iGF_Gm_dd_33) )


          END DO 
          END DO 

        END DO 

        uCR(iX1,iX2,iX3,lo_C(4):hi_C(4)) &
          = RESHAPE( uCR_K, [ hi_C(4) - lo_C(4) + 1 ] )

      END DO ! iX1 = BX % lo(1), BX % hi(1)
      END DO ! iX2 = BX % lo(2), BX % hi(2)
      END DO ! iX3 = BX % lo(3), BX % hi(3)

    END DO ! WHILE( MFI % next() )

    CALL amrex_mfiter_destroy( MFI )

  END SUBROUTINE InitializeFields_SineWaveStreaming


  SUBROUTINE InitializeFields_GaussianDiffusion &
    ( iLevel, MF_uGF, MF_uCR, MF_uCF )

    INTEGER,              INTENT(in   ) :: iLevel
    TYPE(amrex_multifab), INTENT(in   ) :: MF_uGF
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF


    ! --- thornado ---
    INTEGER        :: iDim
    INTEGER        :: iX1, iX2, iX3, iZ1, iZ2, iZ3, iZ4, iS, iNodeZ, iSpecies
    INTEGER        :: iNodeX, iNodeX1, iNodeX2, iNodeX3, iNodeZ2, iNodeZ3, iNodeE
    REAL(DP)       :: X1, X2, X3, t_0, D_min, D_0, X1_0, X2_0
    REAL(DP)       :: uCR_K( nDOFZ, nE, nCR, nSpecies )
    REAL(DP)       :: uPR_K( nDOFZ, nE, nPR, nSpecies )
    REAL(DP)       :: uGF_K( nDOFX, nGF )
    REAL(DP)       :: uPF_K( nDOFX, nPF )
    REAL(DP)       :: uCF_K( nDOFX, nCF )
    REAL(DP)       :: uAF_K( nDOFX, nAF )

    ! --- AMReX ---
    INTEGER                       :: lo_C(4), hi_C(4)
    INTEGER                       :: lo_G(4), hi_G(4)
    INTEGER                       :: lo_F(4), hi_F(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCR(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    REAL(DP)                      :: Three, EN, Sigma
    TYPE(amrex_parmparse) :: PP
    REAL(DP)    , ALLOCATABLE :: V_0(:)

    uCR_K = Zero
    uPF_K = Zero
    uCF_K = Zero
    uGF_K = Zero
    uAF_K = Zero

    D_min = 1.0d-06
    t_0   = 5.0_DP
    X1_0  = One
    X2_0  = One

    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL PP % getarr( 'V_0', &
                         V_0 )
    CALL amrex_parmparse_destroy( PP )

    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL PP % get( 'Sigma', &
                         Sigma )                         
    CALL amrex_parmparse_destroy( PP )
    Three = 3.0_DP


      CALL amrex_mfiter_build( MFI, MF_uCR, tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF % DataPtr( MFI )
        uCR => MF_uCR % DataPtr( MFI )
        uCF => MF_uCF % DataPtr( MFI )

        BX = MFI % tilebox()

        lo_G = LBOUND( uGF )
        hi_G = UBOUND( uGF )

        lo_C = LBOUND( uCR )
        hi_C = UBOUND( uCR )

        lo_F = LBOUND( uCF )
        hi_F = UBOUND( uCF )

        DO iX3 = BX % lo(3), BX % hi(3)
        DO iX2 = BX % lo(2), BX % hi(2)
        DO iX1 = BX % lo(1), BX % hi(1)

          uGF_K &
            = RESHAPE( uGF(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )

          DO iNodeX = 1, nDOFX

            uPF_K(iNodeX,iPF_D ) = 1.0_DP
            uPF_K(iNodeX,iPF_V1) = V_0(1)
            uPF_K(iNodeX,iPF_V2) = V_0(2)
            uPF_K(iNodeX,iPF_V3) = V_0(3)
            uPF_K(iNodeX,iPF_E ) = 0.1_DP
            uPF_K(iNodeX,iPF_Ne) = 0.0_DP

          END DO

        CALL ComputeConserved_Euler_NonRelativistic &
              (  uPF_K(:,iPF_D ), &
                 uPF_K(:,iPF_V1), &
                 uPF_K(:,iPF_V2), &
                 uPF_K(:,iPF_V3), &
                 uPF_K(:,iPF_E ), &
                 uPF_K(:,iPF_Ne), &
                 uCF_K(:,iCF_D ), &
                 uCF_K(:,iCF_S1), &
                 uCF_K(:,iCF_S2), &
                 uCF_K(:,iCF_S3), &
                 uCF_K(:,iCF_E ), &
                 uCF_K(:,iCF_Ne), &
                 uGF_K(:,iGF_Gm_dd_11), &
                 uGF_K(:,iGF_Gm_dd_22), &
                 uGF_K(:,iGF_Gm_dd_33))

        uCF(iX1,iX2,iX3,lo_F(4):hi_F(4)) &
          = RESHAPE( uCF_K, [ hi_F(4) - lo_F(4) + 1 ] )

          DO iNodeZ = 1, nDOFZ

            DO iS = 1, nSpecies
            DO iZ1 = 1, nE

              iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1

              iNodeZ2 = NodeNumberTable(2,iNodeZ)
              iNodeZ3 = NodeNumberTable(3,iNodeZ)

              X1 = NodeCoordinate( MeshX(1), iX1, iNodeZ2 )

              D_0 = One / ( Three * uOP(iNodeZ,iZ1,iX1,iX2,iX3,iOP_Sigma,iS) )

              EN = EXP( - ( (X1-X1_0)*(X1-X1_0) ) / ( 4.0_DP * t_0 * D_0 ) )

              uPR_K( iNodeZ, iZ1, iPR_D, iS ) &
                =  EN

              uPR_K( iNodeZ, iZ1, iPR_I1, iS ) &
                = Half * (X1-X1_0) * uPR_K( iNodeZ, iZ1, iPR_D, iS )/t_0

              uPR_K( iNodeZ, iZ1, iPR_I2, iS ) &
                = 0.0_DP
              uPR_K( iNodeZ, iZ1, iPR_I3, iS ) &
                = 0.0_DP

              IF( uPR_K( iNodeZ, iZ1, iPR_D, iS ) < D_min ) THEN

                uPR_K( iNodeZ, iZ1, iPR_D, iS ) = D_min
                uPR_K( iNodeZ, iZ1, iPR_I1, iS ) = 0.0d-00
                uPR_K( iNodeZ, iZ1, iPR_I2, iS ) = 0.0d-00
                uPR_K( iNodeZ, iZ1, iPR_I2, iS ) = 0.0d-00

              END IF

          CALL ComputeConserved_TwoMoment &
                 ( uPR_K(iNodeZ,iZ1,iPR_D ,iS), &
                   uPR_K(iNodeZ,iZ1,iPR_I1,iS), &
                   uPR_K(iNodeZ,iZ1,iPR_I2,iS), &
                   uPR_K(iNodeZ,iZ1,iPR_I3,iS), &
                   uCR_K(iNodeZ,iZ1,iCR_N ,iS), &
                   uCR_K(iNodeZ,iZ1,iCR_G1,iS), &
                   uCR_K(iNodeZ,iZ1,iCR_G2,iS), &
                   uCR_K(iNodeZ,iZ1,iCR_G3,iS), &
                   uPF_K(iNodeX,iPF_V1), &
                   uPF_K(iNodeX,iPF_V2), &
                   uPF_K(iNodeX,iPF_V3), &
                   uGF_K(iNodeX,iGF_Gm_dd_11), &
                   uGF_K(iNodeX,iGF_Gm_dd_22), &
                   uGF_K(iNodeX,iGF_Gm_dd_33) )

          END DO 
          END DO 

        END DO 

        uCR(iX1,iX2,iX3,lo_C(4):hi_C(4)) &
          = RESHAPE( uCR_K, [ hi_C(4) - lo_C(4) + 1 ] )

        END DO
        END DO
        END DO

      END DO

      CALL amrex_mfiter_destroy( MFI )


  END SUBROUTINE InitializeFields_GaussianDiffusion




  SUBROUTINE InitializeFields_TransparentShock &
    ( iLevel, MF_uGF, MF_uCR, MF_uCF )

    INTEGER             , INTENT(in)    :: iLevel
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF

    ! --- thornado ---

    INTEGER        :: iDim
    INTEGER        :: iX1, iX2, iX3, iZ1, iZ2, iZ3, iZ4, iS, iNodeZ, iSpecies
    INTEGER        :: iNodeX, iNodeX1, iNodeX2, iNodeX3, &
                      iNodeZ2, iNodeZ3, iNodeZ4, iNodeE
    REAL(DP)       :: X1, X2, X3
    REAL(DP)       :: uCR_K(nDOFZ,nE,nCR,nSpecies)
    REAL(DP)       :: uPR_K(nDOFZ,nE,nPR,nSpecies)
    REAL(DP)       :: uGF_K(nDOFX,nGF)
    REAL(DP)       :: uPF_K(nDOFX,nPF)
    REAL(DP)       :: uCF_K(nDOFX,nCF)
    REAL(DP)       :: uAF_K(nDOFX,nAF)

    ! --- AMReX ---

    INTEGER                       :: lo_C(4), hi_C(4), loX(3), hiX(3), loZ(4), hiZ(4)
    INTEGER                       :: lo_G(4), hi_G(4)
    INTEGER                       :: lo_F(4), hi_F(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCR(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    REAL(DP)                      :: ShockWidth, E

    ! --- Problem-dependent parameters ---

    TYPE(amrex_parmparse) :: PP
    REAL(DP)    , ALLOCATABLE :: V_0(:)

    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL PP % getarr( 'V_0', &
                         V_0 )
    CALL amrex_parmparse_destroy( PP )

    ShockWidth = 0.01_DP

    uCR_K = Zero
    uPF_K = Zero
    uCF_K = Zero
    uGF_K = Zero
    uAF_K = Zero

    CALL amrex_mfiter_build( MFI, MF_uCR, tiling = UseTiling )

    DO WHILE( MFI % next() )

      uGF => MF_uGF % DataPtr( MFI )
      uCR => MF_uCR % DataPtr( MFI )
      uCF => MF_uCF % DataPtr( MFI )

      BX = MFI % tilebox()

      lo_G = LBOUND( uGF )
      hi_G = UBOUND( uGF )

      lo_C = LBOUND( uCR )
      hi_C = UBOUND( uCR )

      lo_F = LBOUND( uCF )
      hi_F = UBOUND( uCF )

      loX = BX % lo
      hiX = BX % hi

      DO iX3 = loX(3), hiX(3)
      DO iX2 = loX(2), hiX(2)
      DO iX1 = loX(1), hiX(1)

        uGF_K &
          = RESHAPE( uGF(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )

        DO iNodeX = 1, nDOFX

          X1 = NodeCoordinate( MeshX(1), iX1, iNodeX )

          uPF_K(iNodeX,iPF_D ) = 1.0_DP
          uPF_K(iNodeX,iPF_V2) = 0.0_DP
          uPF_K(iNodeX,iPF_V3) = 0.0_DP
          uPF_K(iNodeX,iPF_E ) = 0.1_DP
          uPF_K(iNodeX,iPF_Ne) = 0.0_DP

          uPF_K(iNodeX,iPF_V1) &
            = Half * V_0(1) * ( One + TANH( (X1-1.0_DP)/ShockWidth ) )
        
        END DO

        CALL ComputeConserved_Euler_NonRelativistic &
              (  uPF_K(:,iPF_D ), &
                 uPF_K(:,iPF_V1), &
                 uPF_K(:,iPF_V2), &
                 uPF_K(:,iPF_V3), &
                 uPF_K(:,iPF_E ), &
                 uPF_K(:,iPF_Ne), &
                 uCF_K(:,iCF_D ), &
                 uCF_K(:,iCF_S1), &
                 uCF_K(:,iCF_S2), &
                 uCF_K(:,iCF_S3), &
                 uCF_K(:,iCF_E ), &
                 uCF_K(:,iCF_Ne), &
                 uGF_K(:,iGF_Gm_dd_11), &
                 uGF_K(:,iGF_Gm_dd_22), &
                 uGF_K(:,iGF_Gm_dd_33))

        uCF(iX1,iX2,iX3,lo_F(4):hi_F(4)) &
          = RESHAPE( uCF_K, [ hi_F(4) - lo_F(4) + 1 ] )

        DO iNodeZ = 1, nDOFZ

          DO iS  = 1, nSpecies
          DO iZ1 = 1, nE

            iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1

            iNodeZ2 = NodeNumberTable(2,iNodeZ)
            iNodeZ3 = NodeNumberTable(3,iNodeZ)
            iNodeZ4 = NodeNumberTable(4,iNodeZ)

            X1 = NodeCoordinate( MeshX(1), iX1, iNodeZ2 )

            iNodeE = MOD( (iNodeZ-1)        , nDOFE ) + 1

            E = NodeCoordinate( MeshE, iZ1, iNodeE )

            uPR_K( iNodeZ, iZ1, iPR_D, iS ) &
               = 1.0d-8
            uPR_K( iNodeZ, iZ1, iPR_I1, iS ) &
               = 0.0_DP
            uPR_K( iNodeZ, iZ1, iPR_I2, iS ) &
              = 0.0_DP
            uPR_K( iNodeZ, iZ1, iPR_I3, iS ) &
              = 0.0_DP

            CALL ComputeConserved_TwoMoment &
                   ( uPR_K(iNodeZ,iZ1,iPR_D ,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I1,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I2,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I3,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_N ,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G1,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G2,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G3,iS), &
                     uPF_K(iNodeX,iPF_V1), &
                     uPF_K(iNodeX,iPF_V2), &
                     uPF_K(iNodeX,iPF_V3), &
                     uGF_K(iNodeX,iGF_Gm_dd_11), &
                     uGF_K(iNodeX,iGF_Gm_dd_22), &
                     uGF_K(iNodeX,iGF_Gm_dd_33) )

          END DO ! iZ1 = 1, nE
          END DO ! iS  = 1, nSpecies

        END DO ! iNodeZ = 1, nDOFZ

        uCR(iX1,iX2,iX3,lo_C(4):hi_C(4)) &
          = RESHAPE( uCR_K, [ hi_C(4) - lo_C(4) + 1 ] )

      END DO ! iX1 = BX % lo(1), BX % hi(1)
      END DO ! iX2 = BX % lo(2), BX % hi(2)
      END DO ! iX3 = BX % lo(3), BX % hi(3)

!!!! Applying boundary conditions

    DO iX3 = loX(3), hiX(3)
    DO iX2 = loX(2), hiX(2)
    DO iX1 = loX(1)-sWX(1), loX(1)-1

      DO iNodeX = 1, nDOFX

        uPF_K(iNodeX,iPF_D ) = One
        uPF_K(iNodeX,iPF_V1) = Zero
        uPF_K(iNodeX,iPF_V2) = Zero
        uPF_K(iNodeX,iPF_V3) = Zero
        uPF_K(iNodeX,iPF_E ) = 1.0d-1
        uPF_K(iNodeX,iPF_Ne) = Zero

        CALL ComputeConserved_Euler_NonRelativistic &
              (  uPF_K(:,iPF_D ), &
                 uPF_K(:,iPF_V1), &
                 uPF_K(:,iPF_V2), &
                 uPF_K(:,iPF_V3), &
                 uPF_K(:,iPF_E ), &
                 uPF_K(:,iPF_Ne), &
                 uCF_K(:,iCF_D ), &
                 uCF_K(:,iCF_S1), &
                 uCF_K(:,iCF_S2), &
                 uCF_K(:,iCF_S3), &
                 uCF_K(:,iCF_E ), &
                 uCF_K(:,iCF_Ne), &
                 uGF_K(:,iGF_Gm_dd_11), &
                 uGF_K(:,iGF_Gm_dd_22), &
                 uGF_K(:,iGF_Gm_dd_33))

      END DO

      uCF(iX1,iX2,iX3,lo_F(4):hi_F(4)) &
        = RESHAPE( uCF_K, [ hi_F(4) - lo_F(4) + 1 ] )

    END DO
    END DO
    END DO
    
    ! --- BC for Radiation Fields ---

    loZ(2:4) = loX
    hiZ(2:4) = hiX
 
    loZ(1) = iZ_B0(1)
    hiZ(1) = iZ_E0(1)

    DO iZ4 = loZ(4), hiZ(4)
    DO iZ3 = loZ(3), hiZ(3)
    DO iZ2 = loZ(2)-sWX(1), loZ(2)-1

      DO iNodeZ = 1, nDOFZ

        DO iS  = 1, nSpecies

          DO iZ1 = loZ(1), hiZ(1)

            iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1

            iNodeE = MOD( (iNodeZ-1)        , nDOFE ) + 1

            E = NodeCoordinate( MeshE, iZ1, iNodeE )

            uPR_K( iNodeZ, iZ1, iPR_D, iS ) &
              = One / ( EXP( E / Three - Three ) + One )
            uPR_K( iNodeZ, iZ1, iPR_I1, iS ) &
              = 0.999_DP * uPR_K( iNodeZ, iZ1, iPR_D, iS )
            uPR_K( iNodeZ, iZ1, iPR_I2, iS ) &
              = 0.0_DP
            uPR_K( iNodeZ, iZ1, iPR_I3, iS ) &
              = 0.0_DP

            CALL ComputeConserved_TwoMoment &
                   ( uPR_K(iNodeZ,iZ1,iPR_D ,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I1,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I2,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I3,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_N ,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G1,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G2,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G3,iS), &
                     uPF_K(iNodeX,iPF_V1), &
                     uPF_K(iNodeX,iPF_V2), &
                     uPF_K(iNodeX,iPF_V3), &
                     uGF_K(iNodeX,iGF_Gm_dd_11), &
                     uGF_K(iNodeX,iGF_Gm_dd_22), &
                     uGF_K(iNodeX,iGF_Gm_dd_33) )
          END DO
        END DO
      END DO

      uCR(iZ2,iZ3,iZ4,lo_C(4):hi_C(4)) &
        = RESHAPE( uCR_K, [ hi_C(4) - lo_C(4) + 1 ] )

    END DO
    END DO
    END DO
    
!!!!!!

    END DO ! WHILE( MFI % next() )

    CALL amrex_mfiter_destroy( MFI )

  END SUBROUTINE InitializeFields_TransparentShock




  SUBROUTINE InitializeFields_StreamingDopplerShift &
    ( iLevel, MF_uGF, MF_uCR, MF_uCF )

    INTEGER             , INTENT(in)    :: iLevel
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF


    ! --- thornado ---

    REAL(DP), PARAMETER :: X_0 = 2.0_DP
    REAL(DP), PARAMETER :: X_1 = 3.5_DP
    REAL(DP), PARAMETER :: X_2 = 6.5_DP
    REAL(DP), PARAMETER :: X_3 = 8.0_DP
    REAL(DP), PARAMETER :: L_X = 6.0_DP

    INTEGER        :: iDim
    INTEGER        :: iX1, iX2, iX3, iZ1, iZ2, iZ3, iZ4, iS, iNodeZ, iSpecies
    INTEGER        :: iNodeX, iNodeX1, iNodeX2, iNodeX3, &
                      iNodeZ2, iNodeZ3, iNodeZ4, iNodeE
    REAL(DP)       :: X1, X2, X3
    REAL(DP)       :: uCR_K(nDOFZ,nE,nCR,nSpecies)
    REAL(DP)       :: uPR_K(nDOFZ,nE,nPR,nSpecies)
    REAL(DP)       :: uGF_K(nDOFX,nGF)
    REAL(DP)       :: uPF_K(nDOFX,nPF)
    REAL(DP)       :: uCF_K(nDOFX,nCF)
    REAL(DP)       :: uAF_K(nDOFX,nAF)
    REAL(DP)       :: E

    ! --- AMReX ---

    INTEGER                       :: lo_C(4), hi_C(4), loX(3), hiX(3), loZ(4), hiZ(4)
    INTEGER                       :: lo_G(4), hi_G(4)
    INTEGER                       :: lo_F(4), hi_F(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCR(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    ! --- Problem-dependent parameters ---

    TYPE(amrex_parmparse) :: PP
    REAL(DP)    , ALLOCATABLE :: V_0(:)

    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL PP % getarr( 'V_0', &
                         V_0 )
    CALL amrex_parmparse_destroy( PP )


    uCR_K = Zero
    uPF_K = Zero
    uCF_K = Zero
    uGF_K = Zero
    uAF_K = Zero

    CALL amrex_mfiter_build( MFI, MF_uCR, tiling = UseTiling )

    DO WHILE( MFI % next() )

      uGF => MF_uGF % DataPtr( MFI )
      uCR => MF_uCR % DataPtr( MFI )
      uCF => MF_uCF % DataPtr( MFI )

      BX = MFI % tilebox()

      lo_G = LBOUND( uGF )
      hi_G = UBOUND( uGF )

      lo_C = LBOUND( uCR )
      hi_C = UBOUND( uCR )

      lo_F = LBOUND( uCF )
      hi_F = UBOUND( uCF )

      loX = BX % lo
      hiX = BX % hi

      DO iX3 = loX(3), hiX(3)
      DO iX2 = loX(2), hiX(2)
      DO iX1 = loX(1), hiX(1)

        uGF_K &
          = RESHAPE( uGF(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )

        DO iNodeX = 1, nDOFX

          X1 = NodeCoordinate( MeshX(1), iX1, iNodeX )

          uPF_K(iNodeX,iPF_D ) = 1.0_DP
          uPF_K(iNodeX,iPF_V2) = 0.0_DP
          uPF_K(iNodeX,iPF_V3) = 0.0_DP
          uPF_K(iNodeX,iPF_E ) = 0.1_DP
          uPF_K(iNodeX,iPF_Ne) = 0.0_DP

          IF( X1 .LT. X_0 )THEN
            uPF_K(iNodeX,iPF_V1) = 0.0_DP
          ELSEIF( X1 .GE. X_0 .AND. X1 .LT. X_1 )THEN
            uPF_K(iNodeX,iPF_V1) = V_0(1) * SIN( TwoPi * ( X1 - X_0 ) / L_X )**2
          ELSEIF( X1 .GE. X_1 .AND. X1 .LT. X_2 )THEN
            uPF_K(iNodeX,iPF_V1) = V_0(1)
          ELSEIF( X1 .GE. X_2 .AND. X1 .LT. X_3 )THEN
            uPF_K(iNodeX,iPF_V1)= V_0(1) * SIN( TwoPi * ( X1 - X_0 ) / L_X )**2
          ELSE
            uPF_K(iNodeX,iPF_V1) = 0.0_DP
          END IF
        
        CALL ComputeConserved_Euler_NonRelativistic &
              (  uPF_K(:,iPF_D ), &
                 uPF_K(:,iPF_V1), &
                 uPF_K(:,iPF_V2), &
                 uPF_K(:,iPF_V3), &
                 uPF_K(:,iPF_E ), &
                 uPF_K(:,iPF_Ne), &
                 uCF_K(:,iCF_D ), &
                 uCF_K(:,iCF_S1), &
                 uCF_K(:,iCF_S2), &
                 uCF_K(:,iCF_S3), &
                 uCF_K(:,iCF_E ), &
                 uCF_K(:,iCF_Ne), &
                 uGF_K(:,iGF_Gm_dd_11), &
                 uGF_K(:,iGF_Gm_dd_22), &
                 uGF_K(:,iGF_Gm_dd_33))

        uCF(iX1,iX2,iX3,lo_F(4):hi_F(4)) &
          = RESHAPE( uCF_K, [ hi_F(4) - lo_F(4) + 1 ] )
        END DO

        DO iNodeZ = 1, nDOFZ

          DO iS  = 1, nSpecies
          DO iZ1 = 1, nE

            iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1

            iNodeZ2 = NodeNumberTable(2,iNodeZ)
            iNodeZ3 = NodeNumberTable(3,iNodeZ)
            iNodeZ4 = NodeNumberTable(4,iNodeZ)

            X1 = NodeCoordinate( MeshX(1), iX1, iNodeZ2 )

            iNodeE = MOD( (iNodeZ-1)        , nDOFE ) + 1

            E = NodeCoordinate( MeshE, iZ1, iNodeE )

            uPR_K( iNodeZ, iZ1, iPR_D, iS ) &
               = 1.0d-40
            uPR_K( iNodeZ, iZ1, iPR_I1, iS ) &
               = 0.0_DP
            uPR_K( iNodeZ, iZ1, iPR_I2, iS ) &
              = 0.0_DP
            uPR_K( iNodeZ, iZ1, iPR_I3, iS ) &
              = 0.0_DP

            CALL ComputeConserved_TwoMoment &
                   ( uPR_K(iNodeZ,iZ1,iPR_D ,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I1,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I2,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I3,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_N ,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G1,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G2,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G3,iS), &
                     uPF_K(iNodeX,iPF_V1), &
                     uPF_K(iNodeX,iPF_V2), &
                     uPF_K(iNodeX,iPF_V3), &
                     uGF_K(iNodeX,iGF_Gm_dd_11), &
                     uGF_K(iNodeX,iGF_Gm_dd_22), &
                     uGF_K(iNodeX,iGF_Gm_dd_33) )

          END DO ! iZ1 = 1, nE
          END DO ! iS  = 1, nSpecies

        END DO ! iNodeZ = 1, nDOFZ

        uCR(iX1,iX2,iX3,lo_C(4):hi_C(4)) &
          = RESHAPE( uCR_K, [ hi_C(4) - lo_C(4) + 1 ] )

      END DO ! iX1 = BX % lo(1), BX % hi(1)
      END DO ! iX2 = BX % lo(2), BX % hi(2)
      END DO ! iX3 = BX % lo(3), BX % hi(3)

!!!! Applying boundary conditions

    DO iX3 = loX(3), hiX(3)
    DO iX2 = loX(2), hiX(2)
    DO iX1 = loX(1)-sWX(1), loX(1)-1

      DO iNodeX = 1, nDOFX

        uPF_K(iNodeX,iPF_D ) = One
        uPF_K(iNodeX,iPF_V1) = Zero
        uPF_K(iNodeX,iPF_V2) = Zero
        uPF_K(iNodeX,iPF_V3) = Zero
        uPF_K(iNodeX,iPF_E ) = 1.0d-1
        uPF_K(iNodeX,iPF_Ne) = Zero

        CALL ComputeConserved_Euler_NonRelativistic &
              (  uPF_K(:,iPF_D ), &
                 uPF_K(:,iPF_V1), &
                 uPF_K(:,iPF_V2), &
                 uPF_K(:,iPF_V3), &
                 uPF_K(:,iPF_E ), &
                 uPF_K(:,iPF_Ne), &
                 uCF_K(:,iCF_D ), &
                 uCF_K(:,iCF_S1), &
                 uCF_K(:,iCF_S2), &
                 uCF_K(:,iCF_S3), &
                 uCF_K(:,iCF_E ), &
                 uCF_K(:,iCF_Ne), &
                 uGF_K(:,iGF_Gm_dd_11), &
                 uGF_K(:,iGF_Gm_dd_22), &
                 uGF_K(:,iGF_Gm_dd_33))

      END DO

      uCF(iX1,iX2,iX3,lo_F(4):hi_F(4)) &
        = RESHAPE( uCF_K, [ hi_F(4) - lo_F(4) + 1 ] )

    END DO
    END DO
    END DO
    
    ! --- BC for Radiation Fields ---

    loZ(2:4) = loX
    hiZ(2:4) = hiX
 
    loZ(1) = iZ_B0(1)
    hiZ(1) = iZ_E0(1)

    DO iZ4 = loZ(4), hiZ(4)
    DO iZ3 = loZ(3), hiZ(3)
    DO iZ2 = loZ(2)-sWX(1), loZ(2)-1

      DO iNodeZ = 1, nDOFZ

        DO iS  = 1, nSpecies

          DO iZ1 = loZ(1), hiZ(1)
          
            iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1

            iNodeE = MOD( (iNodeZ-1)        , nDOFE ) + 1

            E = NodeCoordinate( MeshE, iZ1, iNodeE )

            uPR_K( iNodeZ, iZ1, iPR_D, iS ) &
              = One / ( EXP( E / Three - Three ) + One )
            ! Set up for Fermi-Dirac
            uPR_K( iNodeZ, iZ1, iPR_I1, iS ) &
              = 0.999_DP * uPR_K( iNodeZ, iZ1, iPR_D, iS )
            uPR_K( iNodeZ, iZ1, iPR_I2, iS ) &
              = 0.0_DP
            uPR_K( iNodeZ, iZ1, iPR_I3, iS ) &
              = 0.0_DP

            CALL ComputeConserved_TwoMoment &
                   ( uPR_K(iNodeZ,iZ1,iPR_D ,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I1,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I2,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I3,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_N ,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G1,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G2,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G3,iS), &
                     uPF_K(iNodeX,iPF_V1), &
                     uPF_K(iNodeX,iPF_V2), &
                     uPF_K(iNodeX,iPF_V3), &
                     uGF_K(iNodeX,iGF_Gm_dd_11), &
                     uGF_K(iNodeX,iGF_Gm_dd_22), &
                     uGF_K(iNodeX,iGF_Gm_dd_33) )
          END DO
        END DO
      END DO

      uCR(iZ2,iZ3,iZ4,lo_C(4):hi_C(4)) &
        = RESHAPE( uCR_K, [ hi_C(4) - lo_C(4) + 1 ] )

    END DO
    END DO
    END DO
    
!!!!!!

    END DO ! WHILE( MFI % next() )

    CALL amrex_mfiter_destroy( MFI )

  END SUBROUTINE InitializeFields_StreamingDopplerShift



  SUBROUTINE InitializeFields_TransparentVortex &
    ( iLevel, MF_uGF, MF_uCR, MF_uCF )

    INTEGER             , INTENT(in)    :: iLevel
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF


    ! --- thornado ---


    INTEGER        :: iDim
    INTEGER        :: iX1, iX2, iX3, iZ1, iZ2, iZ3, iZ4, iS, iNodeZ, iSpecies
    INTEGER        :: iNodeX, iNodeX1, iNodeX2, iNodeX3, &
                      iNodeZ2, iNodeZ3, iNodeZ4, iNodeE
    REAL(DP)       :: X1, X2, X3
    REAL(DP)       :: uCR_K(nDOFZ,nE,nCR,nSpecies)
    REAL(DP)       :: uPR_K(nDOFZ,nE,nPR,nSpecies)
    REAL(DP)       :: uGF_K(nDOFX,nGF)
    REAL(DP)       :: uPF_K(nDOFX,nPF)
    REAL(DP)       :: uCF_K(nDOFX,nCF)
    REAL(DP)       :: uAF_K(nDOFX,nAF)
    REAL(DP)       :: Beta, R
    REAL(DP)       :: E, Mu_0
    ! --- AMReX ---

    INTEGER                       :: lo_C(4), hi_C(4), loX(3), hiX(3), loZ(4), hiZ(4)
    INTEGER                       :: lo_G(4), hi_G(4)
    INTEGER                       :: lo_F(4), hi_F(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCR(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    ! --- Problem-dependent parameters ---

    TYPE(amrex_parmparse) :: PP
    REAL(DP)    , ALLOCATABLE :: V_0(:)

    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL PP % getarr( 'V_0', &
                         V_0 )
    CALL amrex_parmparse_destroy( PP )

    Beta = SQRT( V_0(1)**2 + V_0(2)**2 + V_0(3)**2 )

    uCR_K = Zero
    uPF_K = Zero
    uCF_K = Zero
    uGF_K = Zero
    uAF_K = Zero

    CALL amrex_mfiter_build( MFI, MF_uCR, tiling = UseTiling )

    DO WHILE( MFI % next() )

      uGF => MF_uGF % DataPtr( MFI )
      uCR => MF_uCR % DataPtr( MFI )
      uCF => MF_uCF % DataPtr( MFI )

      BX = MFI % tilebox()

      lo_G = LBOUND( uGF )
      hi_G = UBOUND( uGF )

      lo_C = LBOUND( uCR )
      hi_C = UBOUND( uCR )

      lo_F = LBOUND( uCF )
      hi_F = UBOUND( uCF )

      loX = BX % lo
      hiX = BX % hi

      DO iX3 = loX(3), hiX(3)
      DO iX2 = loX(2), hiX(2)
      DO iX1 = loX(1), hiX(1)

        uGF_K &
          = RESHAPE( uGF(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )

        DO iNodeX = 1, nDOFX

          iNodeX1 = NodeNumberTableX(1,iNodeX)
          iNodeX2 = NodeNumberTableX(2,iNodeX)

          X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
          X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )
          R  = SQRT( X1**2 + X2**2 )

          uPF_K(iNodeX,iPF_D ) = 1.0_DP
          uPF_K(iNodeX,iPF_V1) &
            = - X2 * Beta * EXP( Half * ( One - R**2 ) )
          uPF_K(iNodeX,iPF_V2) &
            = + X1 * Beta * EXP( Half * ( One - R**2 ) )
          uPF_K(iNodeX,iPF_V3) = 0.0_DP
          uPF_K(iNodeX,iPF_E ) = 0.1_DP
          uPF_K(iNodeX,iPF_Ne) = 0.0_DP

        CALL ComputeConserved_Euler_NonRelativistic &
              (  uPF_K(:,iPF_D ), &
                 uPF_K(:,iPF_V1), &
                 uPF_K(:,iPF_V2), &
                 uPF_K(:,iPF_V3), &
                 uPF_K(:,iPF_E ), &
                 uPF_K(:,iPF_Ne), &
                 uCF_K(:,iCF_D ), &
                 uCF_K(:,iCF_S1), &
                 uCF_K(:,iCF_S2), &
                 uCF_K(:,iCF_S3), &
                 uCF_K(:,iCF_E ), &
                 uCF_K(:,iCF_Ne), &
                 uGF_K(:,iGF_Gm_dd_11), &
                 uGF_K(:,iGF_Gm_dd_22), &
                 uGF_K(:,iGF_Gm_dd_33))

        END DO

          uCF(iX1,iX2,iX3,lo_F(4):hi_F(4)) &
              = RESHAPE( uCF_K, [ hi_F(4) - lo_F(4) + 1 ] )


        DO iNodeZ = 1, nDOFZ

          DO iS  = 1, nSpecies
          DO iZ1 = 1, nE

            iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1

            iNodeZ2 = NodeNumberTable(2,iNodeZ)
            iNodeZ3 = NodeNumberTable(3,iNodeZ)
            iNodeZ4 = NodeNumberTable(4,iNodeZ)

            X1 = NodeCoordinate( MeshX(1), iX1, iNodeZ2 )

            iNodeE = MOD( (iNodeZ-1)        , nDOFE ) + 1

            E = NodeCoordinate( MeshE, iZ1, iNodeE )

            uPR_K( iNodeZ, iZ1, iPR_D, iS ) &
               = 1.0d-8
            uPR_K( iNodeZ, iZ1, iPR_I1, iS ) &
               = 0.0_DP
            uPR_K( iNodeZ, iZ1, iPR_I2, iS ) &
              = 0.0_DP
            uPR_K( iNodeZ, iZ1, iPR_I3, iS ) &
              = 0.0_DP

            CALL ComputeConserved_TwoMoment &
                   ( uPR_K(iNodeZ,iZ1,iPR_D ,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I1,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I2,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I3,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_N ,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G1,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G2,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G3,iS), &
                     uPF_K(iNodeX,iPF_V1), &
                     uPF_K(iNodeX,iPF_V2), &
                     uPF_K(iNodeX,iPF_V3), &
                     uGF_K(iNodeX,iGF_Gm_dd_11), &
                     uGF_K(iNodeX,iGF_Gm_dd_22), &
                     uGF_K(iNodeX,iGF_Gm_dd_33) )

          END DO ! iZ1 = 1, nE
          END DO ! iS  = 1, nSpecies

        END DO ! iNodeZ = 1, nDOFZ

        uCR(iX1,iX2,iX3,lo_C(4):hi_C(4)) &
          = RESHAPE( uCR_K, [ hi_C(4) - lo_C(4) + 1 ] )

      END DO ! iX1 = BX % lo(1), BX % hi(1)
      END DO ! iX2 = BX % lo(2), BX % hi(2)
      END DO ! iX3 = BX % lo(3), BX % hi(3)


!!!! Applying boundary conditions

! Loop for X direction
    DO iX3 = loX(3), hiX(3)
    DO iX2 = loX(2), hiX(2)
    DO iX1 = loX(1)-sWX(1), loX(1)-1

      DO iNodeX = 1, nDOFX

        uPF_K(iNodeX,iPF_D ) = One
        uPF_K(iNodeX,iPF_V1) = Zero
        uPF_K(iNodeX,iPF_V2) = Zero
        uPF_K(iNodeX,iPF_V3) = Zero
        uPF_K(iNodeX,iPF_E ) = 1.0d-1
        uPF_K(iNodeX,iPF_Ne) = Zero

        CALL ComputeConserved_Euler_NonRelativistic &
              (  uPF_K(:,iPF_D ), &
                 uPF_K(:,iPF_V1), &
                 uPF_K(:,iPF_V2), &
                 uPF_K(:,iPF_V3), &
                 uPF_K(:,iPF_E ), &
                 uPF_K(:,iPF_Ne), &
                 uCF_K(:,iCF_D ), &
                 uCF_K(:,iCF_S1), &
                 uCF_K(:,iCF_S2), &
                 uCF_K(:,iCF_S3), &
                 uCF_K(:,iCF_E ), &
                 uCF_K(:,iCF_Ne), &
                 uGF_K(:,iGF_Gm_dd_11), &
                 uGF_K(:,iGF_Gm_dd_22), &
                 uGF_K(:,iGF_Gm_dd_33))

      END DO

      uCF(iX1,iX2,iX3,lo_F(4):hi_F(4)) &
        = RESHAPE( uCF_K, [ hi_F(4) - lo_F(4) + 1 ] )

    END DO
    END DO
    END DO
!Loop for y direction
    DO iX3 = loX(3), hiX(3)
    DO iX2 = loX(2)-sWX(2), loX(2)-1
    DO iX1 = loX(1), hiX(1)

      DO iNodeX = 1, nDOFX

        uPF_K(iNodeX,iPF_D ) = One
        uPF_K(iNodeX,iPF_V1) = Zero
        uPF_K(iNodeX,iPF_V2) = Zero
        uPF_K(iNodeX,iPF_V3) = Zero
        uPF_K(iNodeX,iPF_E ) = 1.0d-1
        uPF_K(iNodeX,iPF_Ne) = Zero

        CALL ComputeConserved_Euler_NonRelativistic &
              (  uPF_K(:,iPF_D ), &
                 uPF_K(:,iPF_V1), &
                 uPF_K(:,iPF_V2), &
                 uPF_K(:,iPF_V3), &
                 uPF_K(:,iPF_E ), &
                 uPF_K(:,iPF_Ne), &
                 uCF_K(:,iCF_D ), &
                 uCF_K(:,iCF_S1), &
                 uCF_K(:,iCF_S2), &
                 uCF_K(:,iCF_S3), &
                 uCF_K(:,iCF_E ), &
                 uCF_K(:,iCF_Ne), &
                 uGF_K(:,iGF_Gm_dd_11), &
                 uGF_K(:,iGF_Gm_dd_22), &
                 uGF_K(:,iGF_Gm_dd_33))

      END DO

      uCF(iX1,iX2,iX3,lo_F(4):hi_F(4)) &
        = RESHAPE( uCF_K, [ hi_F(4) - lo_F(4) + 1 ] )

    END DO
    END DO
    END DO
    
    
    ! --- BC for Radiation Fields ---


!! Loop Radiation for X
    loZ(2:4) = loX
    hiZ(2:4) = hiX
 
    loZ(1) = iZ_B0(1)
    hiZ(1) = iZ_E0(1)

    DO iZ4 = loZ(4), hiZ(4)
    DO iZ3 = loZ(3), hiZ(3)
    DO iZ2 = loZ(2)-sWX(1), loZ(2)-1

      DO iNodeZ = 1, nDOFZ

        DO iS  = 1, nSpecies

          DO iZ1 = loZ(1), hiZ(1)
          
            iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1

            iNodeE = MOD( (iNodeZ-1)        , nDOFE ) + 1

            E = NodeCoordinate( MeshE, iZ1, iNodeE )

            Mu_0 = 0.9_DP

            uPR_K( iNodeZ, iZ1, iPR_D, iS ) &
              = 0.5_DP * ( One - Mu_0 ) / ( EXP( E / Three - Three ) + One )
            uPR_K( iNodeZ, iZ1, iPR_I1, iS ) &
              = 0.5_DP * ( One + Mu_0 ) * uPR_K( iNodeZ, iZ1, iPR_D, iS )
            uPR_K( iNodeZ, iZ1, iPR_I2, iS ) &
              = 0.0_DP
            uPR_K( iNodeZ, iZ1, iPR_I3, iS ) &
              = 0.0_DP

            CALL ComputeConserved_TwoMoment &
                   ( uPR_K(iNodeZ,iZ1,iPR_D ,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I1,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I2,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I3,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_N ,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G1,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G2,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G3,iS), &
                     uPF_K(iNodeX,iPF_V1), &
                     uPF_K(iNodeX,iPF_V2), &
                     uPF_K(iNodeX,iPF_V3), &
                     uGF_K(iNodeX,iGF_Gm_dd_11), &
                     uGF_K(iNodeX,iGF_Gm_dd_22), &
                     uGF_K(iNodeX,iGF_Gm_dd_33) )
          END DO
        END DO
      END DO

      uCR(iZ2,iZ3,iZ4,lo_C(4):hi_C(4)) &
        = RESHAPE( uCR_K, [ hi_C(4) - lo_C(4) + 1 ] )

    END DO
    END DO
    END DO
    
!!! Loop for y direction Radiation

    DO iZ4 = loZ(4), hiZ(4)
    DO iZ3 = loZ(3)-sWX(2), loZ(3)-1
    DO iZ2 = loZ(2), hiZ(2)

      DO iNodeZ = 1, nDOFZ

        DO iS  = 1, nSpecies

          DO iZ1 = loZ(1), hiZ(1)
          
            iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1

            iNodeE = MOD( (iNodeZ-1)        , nDOFE ) + 1

            E = NodeCoordinate( MeshE, iZ1, iNodeE )

            Mu_0 = 0.9_DP

            uPR_K( iNodeZ, iZ1, iPR_D, iS ) &
              = 0.5_DP * ( One - Mu_0 ) / ( EXP( E / Three - Three ) + One )
            uPR_K( iNodeZ, iZ1, iPR_I1, iS ) &
              = 0.5_DP * ( One + Mu_0 ) * uPR_K( iNodeZ, iZ1, iPR_D, iS )
            uPR_K( iNodeZ, iZ1, iPR_I2, iS ) &
              = 0.0_DP
            uPR_K( iNodeZ, iZ1, iPR_I3, iS ) &
              = 0.0_DP

            CALL ComputeConserved_TwoMoment &
                   ( uPR_K(iNodeZ,iZ1,iPR_D ,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I1,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I2,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I3,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_N ,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G1,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G2,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G3,iS), &
                     uPF_K(iNodeX,iPF_V1), &
                     uPF_K(iNodeX,iPF_V2), &
                     uPF_K(iNodeX,iPF_V3), &
                     uGF_K(iNodeX,iGF_Gm_dd_11), &
                     uGF_K(iNodeX,iGF_Gm_dd_22), &
                     uGF_K(iNodeX,iGF_Gm_dd_33) )
          END DO
        END DO
      END DO

      uCR(iZ2,iZ3,iZ4,lo_C(4):hi_C(4)) &
        = RESHAPE( uCR_K, [ hi_C(4) - lo_C(4) + 1 ] )

    END DO
    END DO
    END DO


!!!!!!

    END DO ! WHILE( MFI % next() )

    CALL amrex_mfiter_destroy( MFI )


  END SUBROUTINE InitializeFields_TransparentVortex



END MODULE MF_InitializationModule
