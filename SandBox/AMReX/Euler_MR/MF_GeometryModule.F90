MODULE MF_GeometryModule

  ! --- AMReX Modules ---

  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_geometry_module, ONLY: &
    amrex_is_all_periodic
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy
  USE amrex_parmparse_module, ONLY: &
    amrex_parmparse, &
    amrex_parmparse_build, &
    amrex_parmparse_destroy
  USE amrex_amrcore_module, ONLY: &
    amrex_geom

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nDOFX
  USE ReferenceElementModuleX, ONLY: &
    nDOFX_X1
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    LX_X1_Dn, &
    LX_X1_Up
  USE GeometryFieldsModule, ONLY: &
    nGF
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX
  USE LinearAlgebraModule, ONLY: &
    MatrixMatrixMultiply

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero, &
    One
  USE MF_UtilitiesModule, ONLY: &
    thornado2amrex_X
  USE InputParsingModule, ONLY: &
    nLevels, &
    nMaxLevels, &
    swX, &
    UseTiling
  USE MF_Euler_TimersModule, ONLY: &
    TimersStart_AMReX_Euler, &
    TimersStop_AMReX_Euler, &
    Timer_AMReX_Euler_Allocate

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeGeometryX_MF
  PUBLIC :: ApplyBoundaryConditions_Geometry_MF

CONTAINS


  SUBROUTINE ComputeGeometryX_MF( MF_uGF )

    TYPE(amrex_multifab), INTENT(in) :: MF_uGF

    TYPE(amrex_mfiter)    :: MFI
    TYPE(amrex_box)       :: BX
    TYPE(amrex_parmparse) :: PP

    INTEGER :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)

    REAL(DP), ALLOCATABLE :: G (:,:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)

    INTEGER :: iLo_G(3), iHi_G(3)

    REAL(DP) :: Mass

    Mass = Zero
    CALL amrex_parmparse_build( PP, 'thornado' )

      CALL PP % query( 'Mass', Mass )

    CALL amrex_parmparse_destroy( PP )

    CALL amrex_mfiter_build( MFI, MF_uGF )

    DO WHILE( MFI % next() )

      uGF => MF_uGF % DataPtr( MFI )

      BX = MFI % TileBox()

      iX_B0 = BX % lo
      iX_E0 = BX % hi
      iX_B1 = iX_B0 - swX
      iX_E1 = iX_E0 + swX

      iLo_G = iX_B1
      iHi_G = iX_E1

      CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_Allocate )

      ALLOCATE( G(1:nDOFX,iX_B1(1):iX_E1(1), &
                          iX_B1(2):iX_E1(2), &
                          iX_B1(3):iX_E1(3), &
                  1:nGF) )

      CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_Allocate )

#if defined HYDRO_RELATIVISTIC

      CALL ComputeGeometryX &
             ( iX_B0, iX_E0, iX_B1, iX_E1, G, Mass_Option = Mass )

#else

      CALL ComputeGeometryX &
             ( iX_B0, iX_E0, iX_B1, iX_E1, G )

#endif

      CALL thornado2amrex_X &
             ( nGF, iX_B1, iX_E1, LBOUND( uGF ), iX_B1, iX_E1, uGF, G )

      CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_Allocate )

      DEALLOCATE( G )

      CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_Allocate )

    END DO

    CALL amrex_mfiter_destroy( MFI )

  END SUBROUTINE ComputeGeometryX_MF


  SUBROUTINE ApplyBoundaryConditions_Geometry_MF( MF_uGF )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:nMaxLevels-1)

    IF( .NOT. amrex_is_all_periodic() )THEN

      CALL ApplyBoundaryConditions_Geometry_MF_X1( MF_uGF )

    END IF

  END SUBROUTINE ApplyBoundaryConditions_Geometry_MF


  SUBROUTINE ApplyBoundaryConditions_Geometry_MF_X1( MF_uGF )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:nMaxLevels-1)

    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)

    INTEGER  :: iLevel, iNX, iX2, iX3, iGF, nX1_X
    INTEGER  :: iX_B0(3), iX_E0(3)

    REAL(DP), ALLOCATABLE :: G_K(:,:,:,:)
    REAL(DP), ALLOCATABLE :: G_F(:,:,:,:)

    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi

        nX1_X = ( iX_E0(3) - iX_B0(3) + 1 ) * ( iX_E0(2) - iX_B0(2) + 1 )

        ALLOCATE( G_K(1:nDOFX   ,iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),1:nGF) )
        ALLOCATE( G_F(1:nDOFX_X1,iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),1:nGF) )

        ! --- Lower boundary ---

        IF( iX_B0(1) .EQ. amrex_geom(iLevel) % domain % lo( 1 ) )THEN

          DO iGF = 1       , nGF
          DO iX3 = iX_B0(3), iX_E0(3)
          DO iX2 = iX_B0(2), iX_E0(2)
          DO iNX = 1       , nDOFX

            G_K(iNX,iX2,iX3,iGF) = uGF(iX_B0(1),iX2,iX3,nDOFX*(iGF-1)+iNX)

          END DO
          END DO
          END DO
          END DO

          DO iGF = 1, nGF

            CALL MatrixMatrixMultiply &
                   ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, One, LX_X1_Dn, &
                     nDOFX_X1,   G_K(1,iX_B0(2),iX_B0(3),iGF), &
                     nDOFX, Zero,G_F(1,iX_B0(2),iX_B0(3),iGF), &
                     nDOFX_X1 )

          END DO

          DO iGF = 1       , nGF
          DO iX3 = iX_B0(3), iX_E0(3)
          DO iX2 = iX_B0(2), iX_E0(2)
          DO iNX = 1       , nDOFX

            uGF(iX_B0(1)-1,iX2,iX3,nDOFX*(iGF-1)+iNX) = G_F(1,iX2,iX3,iGF)

          END DO
          END DO
          END DO
          END DO

        END IF

        ! --- Upper boundary ---

        IF( iX_E0(1) .EQ. amrex_geom(iLevel) % domain % hi( 1 ) )THEN

          IF( .NOT. ALLOCATED( G_K ) ) &
            ALLOCATE( G_K(1:nDOFX   ,iX_B0(2):iX_E0(2), &
                                     iX_B0(3):iX_E0(3),1:nGF) )
          IF( .NOT. ALLOCATED( G_F ) ) &
            ALLOCATE( G_F(1:nDOFX_X1,iX_B0(2):iX_E0(2), &
                                     iX_B0(3):iX_E0(3),1:nGF) )

          DO iGF = 1       , nGF
          DO iX3 = iX_B0(3), iX_E0(3)
          DO iX2 = iX_B0(2), iX_E0(2)
          DO iNX = 1       , nDOFX

            G_K(iNX,iX2,iX3,iGF) = uGF(iX_E0(1),iX2,iX3,nDOFX*(iGF-1)+iNX)

          END DO
          END DO
          END DO
          END DO

          DO iGF = 1, nGF

            CALL MatrixMatrixMultiply &
                   ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, One, LX_X1_Up, &
                     nDOFX_X1,   G_K(1,iX_B0(2),iX_B0(3),iGF), &
                     nDOFX, Zero,G_F(1,iX_B0(2),iX_B0(3),iGF), &
                     nDOFX_X1 )

          END DO

          DO iGF = 1       , nGF
          DO iX3 = iX_B0(3), iX_E0(3)
          DO iX2 = iX_B0(2), iX_E0(2)
          DO iNX = 1       , nDOFX

            uGF(iX_E0(1)+1,iX2,iX3,nDOFX*(iGF-1)+iNX) = G_F(1,iX2,iX3,iGF)

          END DO
          END DO
          END DO
          END DO

        END IF

        DEALLOCATE( G_F )
        DEALLOCATE( G_K )

      END DO

      CALL amrex_mfiter_destroy( MFI )

    END DO

  END SUBROUTINE ApplyBoundaryConditions_Geometry_MF_X1


END MODULE MF_GeometryModule
