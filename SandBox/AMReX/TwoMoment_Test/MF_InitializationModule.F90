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

  USE Euler_UtilitiesModule_Relativistic, ONLY: &
    ComputeConserved_Euler_Relativistic
  USE EquationOfStateModule, ONLY: &
    ComputePressureFromPrimitive
  USE EquationOfStateModule_TABLE, ONLY: &
    ComputeThermodynamicStates_Primitive_TABLE, &
    ApplyEquationOfState_TABLE
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    nX, &
    swX, &
    nNodesX, &
    nDOFE, &
    nDOFZ
  USE RadiationFieldsModule, ONLY: &
    nPR, &
    iCR_N, &
    iCR_G1, &
    iCR_G2, &
    iCR_G3, &
    nCR, &
    iPR_D, &
    iPR_I1, &
    iPR_I2, &
    iPR_I3
  USE FluidFieldsModule, ONLY: &
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
  USE GeometryFieldsModule, ONLY: &
    nGF, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_Alpha, &
    iGF_Beta_1, &
    iGF_Beta_2, &
    iGF_Beta_3, &
    iGF_h_1, &
    iGF_h_2, &
    iGF_h_3, &
    iGF_Phi_N, &
    iGF_Psi, &
    iGF_SqrtGm
  USE TwoMoment_OpacityModule, ONLY: &
   uOP, &
   iOP_Sigma
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
    TwoPi
  USE MF_ErrorModule, ONLY: &
    DescribeError_MF
  USE InputParsingModule, ONLY: &
    ProgramName, &
    nLevels, &
    xL, &
    xR, &
    nE, &
    nSpecies, &
    Mass, &
    UseTiling, &
    R0

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
      CASE( 'SineWaveDiffusion' )

        CALL InitializeFields_SineWaveDiffusion &
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
      CASE( 'HomogeneousSphere1D' )

        CALL InitializeFields_HomogeneousSphere1D &
               ( iLevel, MF_uGF, MF_uCR, MF_uCF )

      CASE( 'HomogeneousSphereGR' )

        CALL InitializeFields_HomogeneousSphereGR &
               ( iLevel, MF_uGF, MF_uCR, MF_uCF )
      CASE( 'ShadowCasting' )

        CALL InitializeFields_ShadowCasting &
               ( iLevel, MF_uGF, MF_uCR, MF_uCF )

      CASE( 'RadiatingSphere' )

        CALL InitializeFields_RadiatingSphere &
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

    INTEGER                       :: lo_C(4), hi_C(4)
    INTEGER                       :: lo_G(4), hi_G(4)
    INTEGER                       :: lo_F(4), hi_F(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCR(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    REAL(DP)                      :: VSq, W

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

    VSq = V_0(1) * V_0(1) + V_0(2) * V_0(2) + V_0(3) * V_0(3)
    W = One / SQRT( One - VSq )

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

        CALL ComputePressureFromPrimitive &
                 ( uPF_K(:,iPF_D), uPF_K(:,iPF_E), uPF_K(:,iPF_Ne), &
                   uAF_K(:,iAF_P) )

        CALL ComputeConserved_Euler_Relativistic &
               ( uPF_K(:,iPF_D ), &
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
                 uGF_K(:,iGF_Gm_dd_33), &
                 uAF_K(:,iAF_P)      )

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
                  = W * uPR_K( iNodeZ, iZ1, iPR_D, iS )

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
                  = W * uPR_K( iNodeZ, iZ1, iPR_D, iS )

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
                  = W * uPR_K( iNodeZ, iZ1, iPR_D, iS )

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
                   uGF_K(iNodeX,iGF_Gm_dd_33), &
                   0.0_DP, 0.0_DP, 0.0_DP, & ! off-diagonal components
                   uGF_K(iNodeX,iGF_Alpha) , &
                   uGF_K(iNodeX,iGF_Beta_1), &
                   uGF_K(iNodeX,iGF_Beta_2), &
                   uGF_K(iNodeX,iGF_Beta_3) )

          END DO ! iZ1 = 1, nE
          END DO ! iS  = 1, nSpecies

        END DO ! iNodeZ = 1, nDOFZ

        uCR(iX1,iX2,iX3,lo_C(4):hi_C(4)) &
          = RESHAPE( uCR_K, [ hi_C(4) - lo_C(4) + 1 ] )

      END DO ! iX1 = BX % lo(1), BX % hi(1)
      END DO ! iX2 = BX % lo(2), BX % hi(2)
      END DO ! iX3 = BX % lo(3), BX % hi(3)

    END DO ! WHILE( MFI % next() )

    CALL amrex_mfiter_destroy( MFI )

  END SUBROUTINE InitializeFields_SineWaveStreaming


  SUBROUTINE InitializeFields_SineWaveDiffusion &
    ( iLevel, MF_uGF, MF_uCR, MF_uCF )

    INTEGER,              INTENT(in   ) :: iLevel
    TYPE(amrex_multifab), INTENT(in   ) :: MF_uGF
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF


    ! --- thornado ---
    INTEGER        :: iDim
    INTEGER        :: iX1, iX2, iX3, iZ1, iZ2, iZ3, iZ4, iS, iNodeZ, iSpecies
    INTEGER        :: iNodeX, iNodeX1, iNodeX2, iNodeX3, iNodeZ2, iNodeE
    REAL(DP)       :: X1, X2, X3
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
    REAL(DP)                      :: W, Third
    TYPE(amrex_parmparse) :: PP
    REAL(DP)    , ALLOCATABLE :: V_0(:)

    uCR_K = Zero
    uPF_K = Zero
    uCF_K = Zero
    uGF_K = Zero
    uAF_K = Zero

    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL PP % getarr( 'V_0', &
                         V_0 )
    CALL amrex_parmparse_destroy( PP )
    Third = 1.0_DP / 3.0_DP

    W = 1.0_DP - ( V_0(1)*V_0(1) + V_0(2)*V_0(2) + V_0(3)*V_0(3) )
    W = 1.0_DP / SQRT( W )



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
        CALL ComputePressureFromPrimitive &
                 ( uPF_K(:,iPF_D), uPF_K(:,iPF_E), uPF_K(:,iPF_Ne), &
                   uAF_K(:,iAF_P) )

        CALL ComputeConserved_Euler_Relativistic &
               ( uPF_K(:,iPF_D ), &
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
                 uGF_K(:,iGF_Gm_dd_33), &
                 uAF_K(:,iAF_P)      )

          uCF(iX1,iX2,iX3,lo_F(4):hi_F(4)) &
            = RESHAPE( uCF_K, [ hi_F(4) - lo_F(4) + 1 ] )

          DO iNodeZ = 1, nDOFZ

            DO iS = 1, nSpecies
            DO iZ1 = 1, nE

              iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1

              iNodeZ2 = NodeNumberTable(2,iNodeZ)

              X1 = NodeCoordinate( MeshX(1), iX1, iNodeZ2 )

              uPR_K( iNodeZ, iZ1, iPR_D, iS ) &
                =  0.49_DP * SIN( Third * Pi * X1 ) + 0.5_DP

              uPR_K( iNodeZ, iZ1, iPR_I1, iS ) &
                = - W * ( 0.49_DP * Pi / ( 9.0_DP * uOP(iNodeZ,iZ1,iX1,iX2,iX3,iOP_Sigma,iS) ) ) &
                  * COS( Third * Pi * X1 )
              uPR_K( iNodeZ, iZ1, iPR_I2, iS ) &
                = 0.0_DP

              uPR_K( iNodeZ, iZ1, iPR_I3, iS ) &
                = 0.0_DP


            CALL ComputeConserved_TwoMoment &
                   ( uPR_K(iNodeZ,iZ1,iPR_D,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I1,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I2,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I3,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_N,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G1,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G2,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G3,iS), &
                     uPF_K(iNodeX,iPF_V1),        &
                     uPF_K(iNodeX,iPF_V2),        &
                     uPF_K(iNodeX,iPF_V3),        &
                     uGF_K(iNodeX,iGF_Gm_dd_11), &
                     uGF_K(iNodeX,iGF_Gm_dd_22), &
                     uGF_K(iNodeX,iGF_Gm_dd_33), &
                     0.0_DP, 0.0_DP, 0.0_DP,     &

                     1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP )


            END DO
            END DO
!Reshape here instead of up top look at Hydro example
          END DO

            uCR(iX1,iX2,iX3,lo_C(4):hi_C(4)) &
              = RESHAPE( uCR_K, [ hi_C(4) - lo_C(4) + 1 ] )


        END DO
        END DO
        END DO

      END DO

      CALL amrex_mfiter_destroy( MFI )


  END SUBROUTINE InitializeFields_SineWaveDiffusion

  SUBROUTINE InitializeFields_GaussianDiffusion &
    ( iLevel, MF_uGF, MF_uCR, MF_uCF )

    INTEGER,              INTENT(in   ) :: iLevel
    TYPE(amrex_multifab), INTENT(in   ) :: MF_uGF
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF


    ! --- thornado ---
    INTEGER        :: iDim
    INTEGER        :: iX1, iX2, iX3, iZ1, iZ2, iZ3, iZ4, iS, iNodeZ, iSpecies
    INTEGER        :: iNodeX, iNodeX1, iNodeX2, iNodeX3, iNodeZ2, iNodeE
    REAL(DP)       :: X1, X2, X3
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
    REAL(DP)                      :: W, Third, EN
    TYPE(amrex_parmparse) :: PP
    REAL(DP)    , ALLOCATABLE :: V_0(:)

    uCR_K = Zero
    uPF_K = Zero
    uCF_K = Zero
    uGF_K = Zero
    uAF_K = Zero

    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL PP % getarr( 'V_0', &
                         V_0 )
    CALL amrex_parmparse_destroy( PP )
    Third = 1.0_DP / 3.0_DP

    W = 1.0_DP - ( V_0(1)*V_0(1) + V_0(2)*V_0(2) + V_0(3)*V_0(3) )
    W = 1.0_DP / SQRT( W )



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
        CALL ComputePressureFromPrimitive &
                 ( uPF_K(:,iPF_D), uPF_K(:,iPF_E), uPF_K(:,iPF_Ne), &
                   uAF_K(:,iAF_P) )

        CALL ComputeConserved_Euler_Relativistic &
               ( uPF_K(:,iPF_D ), &
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
                 uGF_K(:,iGF_Gm_dd_33), &
                 uAF_K(:,iAF_P)      )

          uCF(iX1,iX2,iX3,lo_F(4):hi_F(4)) &
            = RESHAPE( uCF_K, [ hi_F(4) - lo_F(4) + 1 ] )

          DO iNodeZ = 1, nDOFZ

            DO iS = 1, nSpecies
            DO iZ1 = 1, nE

              iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1

              iNodeZ2 = NodeNumberTable(2,iNodeZ)

              X1 = NodeCoordinate( MeshX(1), iX1, iNodeZ2 )

              EN = exp( -9.0_DP * X1**2 ) 
 
              uPR_K( iNodeZ, iZ1, iPR_D, iS ) &
                =  EN / W 

              uPR_K( iNodeZ, iZ1, iPR_I1, iS ) &
                = 0.0_DP

              uPR_K( iNodeZ, iZ1, iPR_I2, iS ) &
                = 0.0_DP

              uPR_K( iNodeZ, iZ1, iPR_I3, iS ) &
                = 0.0_DP


            CALL ComputeConserved_TwoMoment &
                   ( uPR_K(iNodeZ,iZ1,iPR_D,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I1,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I2,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I3,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_N,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G1,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G2,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G3,iS), &
                     uPF_K(iNodeX,iPF_V1),        &
                     uPF_K(iNodeX,iPF_V2),        &
                     uPF_K(iNodeX,iPF_V3),        &
                     uGF_K(iNodeX,iGF_Gm_dd_11), &
                     uGF_K(iNodeX,iGF_Gm_dd_22), &
                     uGF_K(iNodeX,iGF_Gm_dd_33), &
                     0.0_DP, 0.0_DP, 0.0_DP,     &

                     1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP )


            END DO
            END DO
!Reshape here instead of up top look at Hydro example
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

    INTEGER                       :: lo_C(4), hi_C(4)
    INTEGER                       :: lo_G(4), hi_G(4)
    INTEGER                       :: lo_F(4), hi_F(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCR(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    REAL(DP)                      :: VSq, W, ShockWidth, E, S

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

      DO iX3 = BX % lo(3), BX % hi(3)
      DO iX2 = BX % lo(2), BX % hi(2)
      DO iX1 = BX % lo(1), BX % hi(1)

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

        CALL ComputePressureFromPrimitive &
                 ( uPF_K(:,iPF_D), uPF_K(:,iPF_E), uPF_K(:,iPF_Ne), &
                   uAF_K(:,iAF_P) )

        CALL ComputeConserved_Euler_Relativistic &
               ( uPF_K(:,iPF_D ), &
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
                 uGF_K(:,iGF_Gm_dd_33), &
                 uAF_K(:,iAF_P)      )

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

            W = 1.0_DP - (uPF_K(iNodeX,iPF_V1)**2 +  uPF_K(iNodeX,iPF_V2)**2 + uPF_K(iNodeX,iPF_V3)**2 )

            W = 1.0_DP / SQRT( W )

            IF(iX1 .EQ. 1 .AND. iNodeE .EQ. 1) THEN
            print*, E
            END IF

            uPR_K( iNodeZ, iZ1, iPR_D, iS ) &
               = 1.0d-45
              != 1.0_DP / ( EXP( E / 3.0_DP - 3.0_DP ) + 1.0_DP )
            uPR_K( iNodeZ, iZ1, iPR_I1, iS ) &
               = 0.0_DP
              != 0.99_DP * uPR_K( iNodeZ, iZ1, iPR_D, iS )
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
                   uGF_K(iNodeX,iGF_Gm_dd_33), &
                   0.0_DP, 0.0_DP, 0.0_DP, & ! off-diagonal components
                   uGF_K(iNodeX,iGF_Alpha) , &
                   uGF_K(iNodeX,iGF_Beta_1), &
                   uGF_K(iNodeX,iGF_Beta_2), &
                   uGF_K(iNodeX,iGF_Beta_3) )

          END DO ! iZ1 = 1, nE
          END DO ! iS  = 1, nSpecies

        END DO ! iNodeZ = 1, nDOFZ

        uCR(iX1,iX2,iX3,lo_C(4):hi_C(4)) &
          = RESHAPE( uCR_K, [ hi_C(4) - lo_C(4) + 1 ] )

      END DO ! iX1 = BX % lo(1), BX % hi(1)
      END DO ! iX2 = BX % lo(2), BX % hi(2)
      END DO ! iX3 = BX % lo(3), BX % hi(3)

    END DO ! WHILE( MFI % next() )

    CALL amrex_mfiter_destroy( MFI )
  END SUBROUTINE InitializeFields_TransparentShock

  SUBROUTINE InitializeFields_StreamingDopplerShift &
    ( iLevel, MF_uGF, MF_uCR, MF_uCF )


    INTEGER,              INTENT(in)    :: iLevel
    TYPE(amrex_multifab), INTENT(in   ) :: MF_uGF
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF


    REAL(DP), PARAMETER :: X_0 = 2.0_DP
    REAL(DP), PARAMETER :: X_1 = 3.5_DP
    REAL(DP), PARAMETER :: X_2 = 6.5_DP
    REAL(DP), PARAMETER :: X_3 = 8.0_DP
    REAL(DP), PARAMETER :: L_X = 6.0_DP



    ! --- thornado ---
    INTEGER        :: iDim
    INTEGER        :: iX1, iX2, iX3, iZ1, iZ2, iZ3, iZ4, iS, iNodeZ, iSpecies
    INTEGER        :: iNodeX, iNodeX1, iNodeX2, iNodeX3, iNodeZ2, iNodeE
    REAL(DP)       :: X1, X2, X3
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
    INTEGER                       :: iX_B(3), iX_E(3)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCR(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    REAL(DP)                      :: W, Third, E, Mu_0, S
    CHARACTER(len=40) :: name1
    CHARACTER(len=1):: nds
    CHARACTER(len=2)::nxn1
    CHARACTER(len=3)::nxn2

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

    Third = 1.0_DP / 3.0_DP

    W = 1.0_DP - ( V_0(1)*V_0(1) + V_0(2)*V_0(2) + V_0(3)*V_0(3) )
    W = 1.0_DP / SQRT( W )

    Mu_0 = 0.8_DP




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

        iX_B = BX % lo
        iX_E = BX % hi



        DO iX3 = iX_B(3), iX_E(3)
        DO iX2 = iX_B(2), iX_E(2)
        DO iX1 = iX_B(1), iX_E(1)

          uGF_K &
            = RESHAPE( uGF(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )

          DO iNodeX = 1, nDOFX



            iNodeX1 = NodeNumberTableX(1,iNodeX)
            iNodeX2 = NodeNumberTableX(2,iNodeX)
            iNodeX3 = NodeNumberTableX(3,iNodeX)

            X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
            X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )
            X3 = NodeCoordinate( MeshX(3), iX3, iNodeX3 )

            uPF_K(iNodeX,iPF_D ) = 1.0_DP

            SELECT CASE( TRIM( Direction ) )

            CASE( 'X' )

              IF( X1 .LT. X_0 )THEN
                uPF_K(iNodeX,iPF_V1) &
                  = 0.0_DP
              ELSEIF( X1 .GE. X_0 .AND. X1 .LT. X_1 )THEN
                uPF_K(iNodeX,iPF_V1) &
                  = V_0(1) * SIN( TwoPi * ( X1 - X_0 ) / L_X )**2
              ELSEIF( X1 .GE. X_1 .AND. X1 .LT. X_2 )THEN
               uPF_K(iNodeX,iPF_V1) &
                  = V_0(1)
              ELSEIF( X1 .GE. X_2 .AND. X1 .LT. X_3 )THEN
                uPF_K(iNodeX,iPF_V1) &
                  = V_0(1) * SIN( TwoPi * ( X1 - X_0 ) / L_X )**2
              ELSE
                uPF_K(iNodeX,iPF_V1) &
                 = 0.0_DP
              END IF

              uPF_K(iNodeX,iPF_V2) = V_0(2)
              uPF_K(iNodeX,iPF_V3) = V_0(3)


            CASE( 'Y' )

              uPF_K(iNodeX,iPF_V1) = V_0(1)
              IF( X2 .LT. X_0 )THEN
                uPF_K(iNodeX,iPF_V2) &
                  = 0.0_DP
              ELSEIF( X2 .GE. X_0 .AND. X2 .LT. X_1 )THEN
                uPF_K(iNodeX,iPF_V2) &
                  = V_0(2) * SIN( TwoPi * ( X2 - X_0 ) / L_X )**2
              ELSEIF( X2 .GE. X_1 .AND. X2 .LT. X_2 )THEN
                uPF_K(iNodeX,iPF_V2) &
                  = V_0(2)
              ELSEIF( X2 .GE. X_2 .AND. X2 .LT. X_3 )THEN
                uPF_K(iNodeX,iPF_V2) &
                  = V_0(2) * SIN( TwoPi * ( X2 - X_0 ) / L_X )**2
              ELSE
                uPF_K(iNodeX,iPF_V2) &
                  = 0.0_DP
              END IF
              uPF_K(iNodeX,iPF_V3) = V_0(3)

            CASE( 'Z' )

              uPF_K(iNodeX,iPF_V1) = V_0(1)
              uPF_K(iNodeX,iPF_V2) = V_0(2)
              IF( X3 .LT. X_0 )THEN
                uPF_K(iNodeX,iPF_V3) &
                  = 0.0_DP
              ELSEIF( X3 .GE. X_0 .AND. X3 .LT. X_1 )THEN
                uPF_K(iNodeX,iPF_V3) &
                  = V_0(3) * SIN( TwoPi * ( X3 - X_0 ) / L_X )**2
              ELSEIF( X3 .GE. X_1 .AND. X3 .LT. X_2 )THEN
                uPF_K(iNodeX,iPF_V3) &
                  = V_0(3)
              ELSEIF( X3 .GE. X_2 .AND. X3 .LT. X_3 )THEN
                uPF_K(iNodeX,iPF_V3) &
                  = V_0(3) * SIN( TwoPi * ( X3 - X_0 ) / L_X )**2
              ELSE
                uPF_K(iNodeX,iPF_V3) &
                  = 0.0_DP
              END IF

            CASE DEFAULT

              WRITE(*,*)
              WRITE(*,'(A8,A)')    '', 'InitializeFields_StreamingDopplerShift'
              WRITE(*,'(A8,A,A2)') '', 'Invalid Direction: ', TRIM( Direction )
              WRITE(*,*)
              STOP

            END SELECT


            uPF_K(iNodeX,iPF_E ) = 0.1_DP
            uPF_K(iNodeX,iPF_Ne) = 0.0_DP

          END DO

        CALL ComputePressureFromPrimitive &
                 ( uPF_K(:,iPF_D), uPF_K(:,iPF_E), uPF_K(:,iPF_Ne), &
                   uAF_K(:,iAF_P) )

        CALL ComputeConserved_Euler_Relativistic &
               ( uPF_K(:,iPF_D ), &
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
                 uGF_K(:,iGF_Gm_dd_33), &
                 uAF_K(:,iAF_P)      )
          uCF(iX1,iX2,iX3,lo_F(4):hi_F(4)) &
            = RESHAPE( uCF_K, [ hi_F(4) - lo_F(4) + 1 ] )

          DO iNodeZ = 1, nDOFZ

            DO iS = 1, nSpecies
            DO iZ1 = 1, nE

              iNodeE = MOD( (iNodeZ-1)        , nDOFE ) + 1

              iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1

              iNodeZ2 = NodeNumberTable(2,iNodeZ)



              W = 1.0_DP - (uPF_K(iNodeX,iPF_V1)**2 +  uPF_K(iNodeX,iPF_V2)**2 + uPF_K(iNodeX,iPF_V3)**2 )

              W = 1.0_DP / SQRT( W )
              X1 = NodeCoordinate( MeshX(1), iX1, iNodeZ2 )

              E = NodeCoordinate( MeshE, iZ1, iNodeE )

              IF(     TRIM( Direction ) .EQ. 'X' )THEN

                uPR_K( iNodeZ, iZ1, iPR_D, iS ) &
                  ! = 1.0d-40
                  = 1.0_DP / ( EXP( E / 3.0_DP - 3.0_DP ) + 1.0_DP )
                uPR_K( iNodeZ, iZ1, iPR_I1, iS ) &
                  ! = 0.0_DP
                  = 0.99_DP * W * uPR_K( iNodeZ, iZ1, iPR_D, iS )
                uPR_K( iNodeZ, iZ1, iPR_I2, iS ) &
                  = 0.0_DP
                uPR_K( iNodeZ, iZ1, iPR_I3, iS ) &
                  = 0.0_DP
              ELSEIF( TRIM( Direction ) .EQ. 'Y' )THEN

                uPR_K( iNodeZ, iZ1, iPR_D, iS ) &
                  = 1.0_DP / ( EXP( E / 3.0_DP - 3.0_DP ) + 1.0_DP )
                uPR_K( iNodeZ, iZ1, iPR_I1, iS ) &
                  = 0.0_DP
                uPR_K( iNodeZ, iZ1, iPR_I2, iS ) &
                  = 0.99_DP * W * uPR_K( iNodeZ, iZ1, iPR_D, iS )
                uPR_K( iNodeZ, iZ1, iPR_I3, iS ) &
                  = 0.0_DP
              ELSEIF( TRIM( Direction ) .EQ. 'Z' )THEN

                uPR_K( iNodeZ, iZ1, iPR_D, iS ) &
                  = 1.0_DP / ( EXP( E / 3.0_DP - 3.0_DP ) + 1.0_DP )
                uPR_K( iNodeZ, iZ1, iPR_I1, iS ) &
                  = 0.0_DP
                uPR_K( iNodeZ, iZ1, iPR_I2, iS ) &
                  = 0.0_DP
                uPR_K( iNodeZ, iZ1, iPR_I3, iS ) &
                  = 0.99_DP * W * uPR_K( iNodeZ, iZ1, iPR_D, iS )

              END IF

              CALL ComputeConserved_TwoMoment &
                   ( uPR_K(iNodeZ,iZ1,iPR_D,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I1,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I2,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I3,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_N,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G1,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G2,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G3,iS), &
                     uPF_K(iNodeX,iPF_V1),        &
                     uPF_K(iNodeX,iPF_V2),        &
                     uPF_K(iNodeX,iPF_V3),        &
                     uGF_K(iNodeX,iGF_Gm_dd_11), &
                     uGF_K(iNodeX,iGF_Gm_dd_22), &
                     uGF_K(iNodeX,iGF_Gm_dd_33), &
                     0.0_DP, 0.0_DP, 0.0_DP,     &
                     1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP )
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


  END SUBROUTINE InitializeFields_StreamingDopplerShift



  SUBROUTINE InitializeFields_HomogeneousSphere1D &
    ( iLevel, MF_uGF, MF_uCR, MF_uCF )

    INTEGER,              INTENT(in   ) :: iLevel
    TYPE(amrex_multifab), INTENT(in   ) :: MF_uGF
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF


    ! --- thornado ---
    INTEGER        :: iDim
    INTEGER        :: iX1, iX2, iX3, iZ1, iZ2, iZ3, iZ4, iS, iNodeZ, iSpecies
    INTEGER        :: iNodeX, iNodeX1, iNodeX2, iNodeX3, iNodeZ2, iNodeE
    REAL(DP)       :: X1, X2, X3, V_0(3)
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
    REAL(DP)                      :: W

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

          DO iNodeX = 1, nDOFX

            uPF_K(iNodeX,iPF_D ) = 1.0_DP
            uPF_K(iNodeX,iPF_V1) = 0.0_DP
            uPF_K(iNodeX,iPF_V2) = 0.0_DP
            uPF_K(iNodeX,iPF_V3) = 0.0_DP
            uPF_K(iNodeX,iPF_E ) = 0.1_DP
            uPF_K(iNodeX,iPF_Ne) = 0.0_DP

          END DO
        CALL ComputePressureFromPrimitive &
                 ( uPF_K(:,iPF_D), uPF_K(:,iPF_E), uPF_K(:,iPF_Ne), &
                   uAF_K(:,iAF_P) )

        CALL ComputeConserved_Euler_Relativistic &
               ( uPF_K(:,iPF_D ), &
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
                 uGF_K(:,iGF_Gm_dd_33), &
                 uAF_K(:,iAF_P)      )


          uCF(iX1,iX2,iX3,lo_F(4):hi_F(4)) &
            = RESHAPE( uCF_K, [ hi_F(4) - lo_F(4) + 1 ] )

          DO iNodeZ = 1, nDOFZ

            DO iS = 1, nSpecies
            DO iZ1 = 1, nE

              iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1

              iNodeZ2 = NodeNumberTable(2,iNodeZ)

              uPR_K( iNodeZ, iZ1, iPR_D, iS ) &
                = 10d-8
              uPR_K( iNodeZ, iZ1, iPR_I1, iS ) &
                 = 0.0_DP
              uPR_K( iNodeZ, iZ1, iPR_I2, iS ) &
                = 0.0_DP

              uPR_K( iNodeZ, iZ1, iPR_I3, iS ) &
                = 0.0_DP


            CALL ComputeConserved_TwoMoment &
                   ( uPR_K(iNodeZ,iZ1,iPR_D,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I1,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I2,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I3,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_N,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G1,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G2,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G3,iS), &
                     uPF_K(iNodeX,iPF_V1),        &
                     uPF_K(iNodeX,iPF_V2),        &
                     uPF_K(iNodeX,iPF_V3),        &
                     uGF_K(iNodeX,iGF_Gm_dd_11), &
                     uGF_K(iNodeX,iGF_Gm_dd_22), &
                     uGF_K(iNodeX,iGF_Gm_dd_33), &
                     0.0_DP, 0.0_DP, 0.0_DP,     &
                     uGF_K(iNodeX,iGF_Alpha), &
                     uGF_K(iNodeX,iGF_Beta_1), &
                     uGF_K(iNodeX,iGF_Beta_2), &
                     uGF_K(iNodeX,iGF_Beta_3) )


            END DO
            END DO
!Reshape here instead of up top look at Hydro example
          END DO

            uCR(iX1,iX2,iX3,lo_C(4):hi_C(4)) &
              = RESHAPE( uCR_K, [ hi_C(4) - lo_C(4) + 1 ] )


        END DO
        END DO
        END DO

      END DO

      CALL amrex_mfiter_destroy( MFI )


  END SUBROUTINE InitializeFields_HomogeneousSphere1D


  SUBROUTINE InitializeFields_HomogeneousSphereGR &
    ( iLevel, MF_uGF, MF_uCR, MF_uCF )

    INTEGER,              INTENT(in   ) :: iLevel
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF


    ! --- thornado ---
    INTEGER        :: iDim
    INTEGER        :: iX1, iX2, iX3, iZ1, iZ2, iZ3, iZ4, iS, iNodeZ, iSpecies
    INTEGER        :: iNodeX, iNodeX1, iNodeX2, iNodeX3, iNodeZ2, iNodeE
    REAL(DP)       :: X1, X2, X3, V_0(3)
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
    REAL(DP)                      :: W, R, Theta, E

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
        DO iX1 = BX % lo(1) - swX(1), BX % hi(1) + swX(1)

          uGF_K &
            = RESHAPE( uGF(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )

          DO iNodeX = 1, nDOFX

            iNodeX1 = NodeNumberTableX(1,iNodeX)

            iNodeX2 = NodeNumberTableX(2,iNodeX)

            R = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
            theta = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

            uPF_K(iNodeX,iPF_D ) = 1.0d12 * Gram / Centimeter**3
            uPF_K(iNodeX,iPF_V1) = 0.0_DP
            uPF_K(iNodeX,iPF_V2) = 0.0_DP
            uPF_K(iNodeX,iPF_V3) = 0.0_DP
            uPF_K(iNodeX,iPF_E ) = 1.00d25 * Erg / Centimeter**3
            uPF_K(iNodeX,iPF_Ne) = 0.0_DP

            CALL ComputeAlphaPsi( Mass, R, R0, theta, uGF_K(iNodeX,:) )
          END DO
        CALL ComputePressureFromPrimitive &
                 ( uPF_K(:,iPF_D), uPF_K(:,iPF_E), uPF_K(:,iPF_Ne), &
                   uAF_K(:,iAF_P) )

        CALL ComputeConserved_Euler_Relativistic &
               ( uPF_K(:,iPF_D ), &
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
                 uGF_K(:,iGF_Gm_dd_33), &
                 uAF_K(:,iAF_P)      )

          uCF(iX1,iX2,iX3,lo_F(4):hi_F(4)) &
            = RESHAPE( uCF_K, [ hi_F(4) - lo_F(4) + 1 ] )

          uGF(iX1,iX2,iX3,lo_G(4):hi_G(4)) &
            = RESHAPE( uGF_K, [ hi_G(4) - lo_G(4) + 1 ] )

          DO iNodeZ = 1, nDOFZ

            DO iS = 1, nSpecies
            DO iZ1 = 1, nE

              iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1

              iNodeZ2 = NodeNumberTable(2,iNodeZ)
              iNodeE = MOD( (iNodeZ-1)        , nDOFE ) + 1
              E = NodeCoordinate( MeshE, iZ1, iNodeE )

              IF (iX1 .EQ. 1) THEN
              print*, E / MeV
              END IF

              uPR_K( iNodeZ, iZ1, iPR_D, iS ) &
                = 1.0d-8
              uPR_K( iNodeZ, iZ1, iPR_I1, iS ) &
                 = 0.0_DP
              uPR_K( iNodeZ, iZ1, iPR_I2, iS ) &
                = 0.0_DP

              uPR_K( iNodeZ, iZ1, iPR_I3, iS ) &
                = 0.0_DP


            CALL ComputeConserved_TwoMoment &
                   ( uPR_K(iNodeZ,iZ1,iPR_D,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I1,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I2,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I3,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_N,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G1,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G2,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G3,iS), &
                     uPF_K(iNodeX,iPF_V1),        &
                     uPF_K(iNodeX,iPF_V2),        &
                     uPF_K(iNodeX,iPF_V3),        &
                     uGF_K(iNodeX,iGF_Gm_dd_11), &
                     uGF_K(iNodeX,iGF_Gm_dd_22), &
                     uGF_K(iNodeX,iGF_Gm_dd_33), &
                     0.0_DP, 0.0_DP, 0.0_DP,     &
                     uGF_K(iNodeX,iGF_Alpha), &
                     uGF_K(iNodeX,iGF_Beta_1), &
                     uGF_K(iNodeX,iGF_Beta_2), &
                     uGF_K(iNodeX,iGF_Beta_3) )


            END DO
            END DO
!Reshape here instead of up top look at Hydro example
          END DO

            uCR(iX1,iX2,iX3,lo_C(4):hi_C(4)) &
              = RESHAPE( uCR_K, [ hi_C(4) - lo_C(4) + 1 ] )


        END DO
        END DO
        END DO

      END DO

      CALL amrex_mfiter_destroy( MFI )


  END SUBROUTINE InitializeFields_HomogeneousSphereGR

  SUBROUTINE InitializeFields_ShadowCasting &
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

    INTEGER                       :: lo_C(4), hi_C(4)
    INTEGER                       :: lo_G(4), hi_G(4)
    INTEGER                       :: lo_F(4), hi_F(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCR(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    REAL(DP)                      :: VSq, W




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

        DO iNodeX = 1, nDOFX

          uPF_K(iNodeX,iPF_D ) = 1.0_DP
          uPF_K(iNodeX,iPF_V1) = 0.0_DP
          uPF_K(iNodeX,iPF_V2) = 0.0_DP
          uPF_K(iNodeX,iPF_V3) = 0.0_DP
          uPF_K(iNodeX,iPF_E ) = 0.1_DP
          uPF_K(iNodeX,iPF_Ne) = 0.0_DP

        END DO

        CALL ComputePressureFromPrimitive &
                 ( uPF_K(:,iPF_D), uPF_K(:,iPF_E), uPF_K(:,iPF_Ne), &
                   uAF_K(:,iAF_P) )

        CALL ComputeConserved_Euler_Relativistic &
               ( uPF_K(:,iPF_D ), &
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
                 uGF_K(:,iGF_Gm_dd_33), &
                 uAF_K(:,iAF_P)      )

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


            uPR_K( iNodeZ, iZ1, iPR_D, iS ) &
              = 10d-10
   
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
                     uGF_K(iNodeX,iGF_Gm_dd_33), &
                     0.0_DP, 0.0_DP, 0.0_DP, & ! off-diagonal components
                     uGF_K(iNodeX,iGF_Alpha) , &
                     uGF_K(iNodeX,iGF_Beta_1), &
                     uGF_K(iNodeX,iGF_Beta_2), &
                     uGF_K(iNodeX,iGF_Beta_3) )

          END DO ! iZ1 = 1, nE
          END DO ! iS  = 1, nSpecies

        END DO ! iNodeZ = 1, nDOFZ

        uCR(iX1,iX2,iX3,lo_C(4):hi_C(4)) &
          = RESHAPE( uCR_K, [ hi_C(4) - lo_C(4) + 1 ] )

      END DO ! iX1 = BX % lo(1), BX % hi(1)
      END DO ! iX2 = BX % lo(2), BX % hi(2)
      END DO ! iX3 = BX % lo(3), BX % hi(3)

    END DO ! WHILE( MFI % next() )

    CALL amrex_mfiter_destroy( MFI )

  END SUBROUTINE InitializeFields_ShadowCasting

  SUBROUTINE InitializeFields_RadiatingSphere &
    ( iLevel, MF_uGF, MF_uCR, MF_uCF )

    INTEGER,              INTENT(in   ) :: iLevel
    TYPE(amrex_multifab), INTENT(in   ) :: MF_uGF
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF


    ! --- thornado ---
    INTEGER        :: iDim
    INTEGER        :: iX1, iX2, iX3, iZ1, iZ2, iZ3, iZ4, iS, iNodeZ, iSpecies
    INTEGER        :: iNodeX, iNodeX1, iNodeX2, iNodeX3, iNodeZ2, iNodeE
    REAL(DP)       :: X1, X2, X3, V_0(3)
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

    REAL(DP) :: V_Max = 0.2_DP

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

          DO iNodeX = 1, nDOFX

            uPF_K(iNodeX,iPF_D ) = 1.0_DP
            uPF_K(iNodeX,iPF_V2) = 0.0_DP
            uPF_K(iNodeX,iPF_V3) = 0.0_DP
            uPF_K(iNodeX,iPF_E ) = 0.1_DP
            uPF_K(iNodeX,iPF_Ne) = 0.0_DP

            iNodeX1 = NodeNumberTableX(1,iNodeX)

            X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

            IF( X1 <= 135.0_DP )THEN

              uPF_K(iNodeX,iPF_V1) &
                = 0.0_DP

            ELSEIF( X1 > 135.0_DP .AND. X1 <= 150.0_DP )THEN

              uPF_K(iNodeX,iPF_V1) &
                = - V_Max * ( X1 - 135.0_DP ) / 15.0_DP

            ELSEIF( X1 > 150.0_DP )THEN

              uPF_K(iNodeX,iPF_V1) &
                = - V_Max * ( 150.0_DP / X1 )**2

            END IF

          END DO
        CALL ComputePressureFromPrimitive &
                 ( uPF_K(:,iPF_D), uPF_K(:,iPF_E), uPF_K(:,iPF_Ne), &
                   uAF_K(:,iAF_P) )

        CALL ComputeConserved_Euler_Relativistic &
               ( uPF_K(:,iPF_D ), &
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
                 uGF_K(:,iGF_Gm_dd_33), &
                 uAF_K(:,iAF_P)      )


          uCF(iX1,iX2,iX3,lo_F(4):hi_F(4)) &
            = RESHAPE( uCF_K, [ hi_F(4) - lo_F(4) + 1 ] )

          DO iNodeZ = 1, nDOFZ

            DO iS = 1, nSpecies
            DO iZ1 = 1, nE

              iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1

              iNodeZ2 = NodeNumberTable(2,iNodeZ)

              uPR_K( iNodeZ, iZ1, iPR_D, iS ) &
                = 10d-40
              uPR_K( iNodeZ, iZ1, iPR_I1, iS ) &
                 = 0.0_DP
              uPR_K( iNodeZ, iZ1, iPR_I2, iS ) &
                = 0.0_DP

              uPR_K( iNodeZ, iZ1, iPR_I3, iS ) &
                = 0.0_DP


            CALL ComputeConserved_TwoMoment &
                   ( uPR_K(iNodeZ,iZ1,iPR_D,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I1,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I2,iS), &
                     uPR_K(iNodeZ,iZ1,iPR_I3,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_N,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G1,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G2,iS), &
                     uCR_K(iNodeZ,iZ1,iCR_G3,iS), &
                     uPF_K(iNodeX,iPF_V1),        &
                     uPF_K(iNodeX,iPF_V2),        &
                     uPF_K(iNodeX,iPF_V3),        &
                     uGF_K(iNodeX,iGF_Gm_dd_11), &
                     uGF_K(iNodeX,iGF_Gm_dd_22), &
                     uGF_K(iNodeX,iGF_Gm_dd_33), &
                     0.0_DP, 0.0_DP, 0.0_DP,     &
                     uGF_K(iNodeX,iGF_Alpha), &
                     uGF_K(iNodeX,iGF_Beta_1), &
                     uGF_K(iNodeX,iGF_Beta_2), &
                     uGF_K(iNodeX,iGF_Beta_3) )


            END DO
            END DO
!Reshape here instead of up top look at Hydro example
          END DO

            uCR(iX1,iX2,iX3,lo_C(4):hi_C(4)) &
              = RESHAPE( uCR_K, [ hi_C(4) - lo_C(4) + 1 ] )


        END DO
        END DO
        END DO

      END DO
      CALL amrex_mfiter_destroy( MFI )


  END SUBROUTINE InitializeFields_RadiatingSphere

  SUBROUTINE ComputeAlphaPsi( M, R, R0, theta, G )

    REAL(DP), INTENT(in) :: M, R, R0, theta
    REAL(DP), INTENT(inout) :: G(nGF)


    REAL(DP) :: Phi, Psi, a, b

    a = GravitationalConstant * M  / ( 2.0_DP * ( R0 )**3 )
    b = GravitationalConstant * M
      !Need to add constants here
    IF ( R .LE. R0 ) THEN

      G(iGF_Phi_N) =  a * ( ( R )**2 - 3.0_DP * ( R0 ) **2 )

    ELSE

      G(iGF_Phi_N) = - b / ( R )

    END IF

    G(iGF_Phi_N) = G(iGF_Phi_N) / ( SpeedOfLight )**2



    G(iGF_Alpha) = 1.0_DP + G(iGF_Phi_N)

    G(iGF_Psi) = 1.0_DP - G(iGF_Phi_N) / 2.0_DP


    G(iGF_h_1) = G(iGF_Psi)**2
    G(iGF_h_2) = G(iGF_Psi)**2 * R
    G(iGF_h_3) = G(iGF_Psi)**2 * R * sin(theta)


    G(iGF_Gm_dd_11) = G(iGF_h_1)**2
    G(iGF_Gm_dd_22) = G(iGF_h_2)**2
    G(iGF_Gm_dd_33) = G(iGF_h_3)**2

    G(iGF_SqrtGm) = SQRT( G(iGF_Gm_dd_11) * G(iGF_Gm_dd_22) * G(iGF_Gm_dd_33) )


  END SUBROUTINE ComputeAlphaPsi


END MODULE MF_InitializationModule
