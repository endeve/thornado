MODULE MF_InitializationModule

  ! --- AMReX Modules ---

  USE amrex_base_module

  ! --- thornado Modules ---

  USE KindModule,              ONLY: &
    DP, Half, One, Pi, TwoPi
  USE ProgramHeaderModule,     ONLY: &
    nDOFX, nX, nNodesX
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
    nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iCF_Ne
  USE Euler_UtilitiesModule,   ONLY: &
    ComputeConserved

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MF_InitializeFields

CONTAINS


  SUBROUTINE MF_InitializeFields( ProgramName, MF_uGF, MF_uCF )

    CHARACTER(LEN=*),     INTENT(in   ) :: ProgramName
    TYPE(amrex_multifab), INTENT(in   ) :: MF_uGF
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF

    IF( amrex_parallel_ioprocessor() )THEN
      WRITE(*,*)
      WRITE(*,'(A4,A,A)') '', 'Initializing: ', TRIM( ProgramName )
    END IF

    SELECT CASE ( TRIM( ProgramName ) )

      CASE ( 'IsentropicVortex' )

        CALL InitializeFields_IsentropicVortex( MF_uGF, MF_uCF )

      CASE DEFAULT

        IF( amrex_parallel_ioprocessor() )THEN
          WRITE(*,*)
          WRITE(*,'(A4,A,A)') '', 'Unknown Program: ', TRIM( ProgramName )
        END IF

    END SELECT

  END SUBROUTINE MF_InitializeFields


  SUBROUTINE InitializeFields_IsentropicVortex( MF_uGF, MF_uCF )

    TYPE(amrex_multifab), INTENT(in   ) :: MF_uGF
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF

    REAL(amrex_real), PARAMETER :: Beta  = 5.0_amrex_real
    REAL(amrex_real), PARAMETER :: Gamma = 1.4_amrex_real

    INTEGER            :: iDim
    INTEGER            :: iX1, iX2, iX3
    INTEGER            :: iNodeX
    INTEGER            :: iNodeX1, iNodeX2
    INTEGER            :: lo_G(4), hi_G(4)
    INTEGER            :: lo_F(4), hi_F(4)
    REAL(amrex_real)   :: X1, X2, R
    REAL(amrex_real)   :: uGF_K(nDOFX,nGF)
    REAL(amrex_real)   :: uPF_K(nDOFX,nPF)
    REAL(amrex_real)   :: uCF_K(nDOFX,nCF)
    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI_F
    TYPE(amrex_mfiter) :: MFI_G
    TYPE(MeshType)     :: MeshX(3)
    REAL(amrex_real), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    REAL(amrex_real), CONTIGUOUS, POINTER :: uGF(:,:,:,:)

    DO iDim = 1, 3

      CALL CreateMesh &
             ( MeshX(iDim), nX(iDim), nNodesX(iDim), 0, &
               amrex_problo(iDim), amrex_probhi(iDim) )

    END DO

    CALL amrex_mfiter_build( MFI_G, MF_uGF, tiling = .TRUE. )
    CALL amrex_mfiter_build( MFI_F, MF_uCF, tiling = .TRUE. )

    DO WHILE( MFI_G % next() .AND. MFI_F % next() )

      uGF => MF_uGF % DataPtr( MFI_G )
      uCF => MF_uCF % DataPtr( MFI_F )

      BX = MFI_F % tilebox()

      lo_F = LBOUND( uCF )
      hi_F = UBOUND( uCF )

      lo_G = LBOUND( uGF )
      hi_G = UBOUND( uGF )

      DO iX3 = BX % lo(3), BX % hi(3)
      DO iX2 = BX % lo(2), BX % hi(2)
      DO iX1 = BX % lo(1), BX % hi(1)

        uGF_K(1:nDOFX,1:nGF) &
          = RESHAPE( uGF(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )

        DO iNodeX = 1, nDOFX

          iNodeX1 = NodeNumberTableX(1,iNodeX)
          iNodeX2 = NodeNumberTableX(2,iNodeX)

          X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
          X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

          R = SQRT( X1**2 + X2**2 )

          uPF_K(iNodeX,iPF_D) &
            = ( One - ( Gamma - One ) * Beta**2 &
                      / ( 8.0_DP * Gamma * Pi**2 ) * EXP( One - R**2 ) &
              )**( One / ( Gamma - One ) )

          uPF_K(iNodeX,iPF_V1) &
            = One - X2 * ( Beta / TwoPi ) * EXP( Half * ( One - R**2 ) )

          uPF_K(iNodeX,iPF_V2) &
            = One + X1 * ( Beta / TwoPi ) * EXP( Half * ( One - R**2 ) )

          uPF_K(iNodeX,iPF_V3) &
            = 0.0_DP

          uPF_K(iNodeX,iPF_E) &
            = uPF_K(iNodeX,iPF_D)**Gamma / ( Gamma - One )

        END DO

        CALL ComputeConserved &
               ( uPF_K(:,iPF_D ), uPF_K(:,iPF_V1), uPF_K(:,iPF_V2), &
                 uPF_K(:,iPF_V3), uPF_K(:,iPF_E ), uPF_K(:,iPF_Ne), &
                 uCF_K(:,iCF_D ), uCF_K(:,iCF_S1), uCF_K(:,iCF_S2), &
                 uCF_K(:,iCF_S3), uCF_K(:,iCF_E ), uCF_K(:,iCF_Ne), &
                 uGF_K(:,iGF_Gm_dd_11), &
                 uGF_K(:,iGF_Gm_dd_22), &
                 uGF_K(:,iGF_Gm_dd_33) )

        uCF(iX1,iX2,iX3,lo_F(4):hi_F(4)) &
          = RESHAPE( uCF_K(1:nDOFX,1:nCF), [ hi_F(4) - lo_F(4) + 1 ] )

      END DO
      END DO
      END DO

    END DO

    CALL amrex_mfiter_destroy( MFI_G )
    CALL amrex_mfiter_destroy( MFI_F )

    DO iDim = 1, 3

      CALL DestroyMesh( MeshX(iDim) )

    END DO

  END SUBROUTINE InitializeFields_IsentropicVortex


END MODULE MF_InitializationModule
