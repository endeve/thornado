MODULE MF_Euler_UtilitiesModule

  ! --- AMReX Modules ---

  USE amrex_base_module

  ! --- thornado Modules ---

  USE ProgramHeaderModule,   ONLY: &
    nDOFX
  USE GeometryFieldsModule,  ONLY: &
    nGF, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33
  USE FluidFieldsModule,     ONLY: &
    nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iPF_Ne, &
    nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iCF_Ne, &
    nAF, iAF_P, iAF_Cs
  USE Euler_UtilitiesModule, ONLY: &
    ComputePrimitive, &
    ComputeTimeStep
  USE EquationOfStateModule, ONLY: &
    ComputePressureFromPrimitive, &
    ComputeSoundSpeedFromPrimitive

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MF_ComputeFromConserved
  PUBLIC :: MF_ComputeTimeStep

CONTAINS


  SUBROUTINE MF_ComputeFromConserved( MF_uGF, MF_uCF, MF_uPF, MF_uAF )

    TYPE(amrex_multifab), INTENT(in   ) :: MF_uGF
    TYPE(amrex_multifab), INTENT(in   ) :: MF_uCF
    TYPE(amrex_multifab), INTENT(inout) :: MF_uPF
    TYPE(amrex_multifab), INTENT(inout) :: MF_uAF

    INTEGER            :: iX1, iX2, iX3
    INTEGER            :: lo_G(4), hi_G(4)
    INTEGER            :: lo_C(4), hi_C(4)
    INTEGER            :: lo_P(4), hi_P(4)
    INTEGER            :: lo_A(4), hi_A(4)
    REAL(amrex_real)   :: uGF_K(nDOFX,nGF)
    REAL(amrex_real)   :: uCF_K(nDOFX,nCF)
    REAL(amrex_real)   :: uPF_K(nDOFX,nPF)
    REAL(amrex_real)   :: uAF_K(nDOFX,nAF)
    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI
    REAL(amrex_real), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(amrex_real), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    REAL(amrex_real), CONTIGUOUS, POINTER :: uPF(:,:,:,:)
    REAL(amrex_real), CONTIGUOUS, POINTER :: uAF(:,:,:,:)

    IF( amrex_parallel_ioprocessor() )THEN
      WRITE(*,*)
      WRITE(*,'(A4,A)') '', 'MF_ComputeFromConserved'
    END IF

    uGF_K = 0.0d0
    uCF_K = 0.0d0
    uPF_K = 0.0d0
    uAF_K = 0.0d0

    CALL amrex_mfiter_build( MFI, MF_uGF, tiling = .TRUE. )

    DO WHILE( MFI % next() )

      uGF => MF_uGF % DataPtr( MFI )
      uCF => MF_uCF % DataPtr( MFI )
      uPF => MF_uPF % DataPtr( MFI )
      uAF => MF_uAF % DataPtr( MFI )

      BX = MFI % tilebox()

      lo_G = LBOUND( uGF ); hi_G = UBOUND( uGF )
      lo_C = LBOUND( uCF ); hi_C = UBOUND( uCF )
      lo_P = LBOUND( uPF ); hi_P = UBOUND( uPF )
      lo_A = LBOUND( uAF ); hi_A = UBOUND( uAF )

      DO iX3 = BX % lo(3), BX % hi(3)
      DO iX2 = BX % lo(2), BX % hi(2)
      DO iX1 = BX % lo(1), BX % hi(1)

        uGF_K &
          = RESHAPE( uGF(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )

        uCF_K &
          = RESHAPE( uCF(iX1,iX2,iX3,lo_C(4):hi_C(4)), [ nDOFX, nCF ] )

        ! --- Primitive Fluid ---

        CALL ComputePrimitive &
               ( uCF_K(:,iCF_D ), uCF_K(:,iCF_S1), uCF_K(:,iCF_S2), &
                 uCF_K(:,iCF_S3), uCF_K(:,iCF_E ), uCF_K(:,iCF_Ne), &
                 uPF_K(:,iPF_D ), uPF_K(:,iPF_V1), uPF_K(:,iPF_V2), &
                 uPF_K(:,iPF_V3), uPF_K(:,iPF_E ), uPF_K(:,iPF_Ne), &
                 uGF_K(:,iGF_Gm_dd_11), &
                 uGF_K(:,iGF_Gm_dd_22), &
                 uGF_K(:,iGF_Gm_dd_33) )

        uPF(iX1,iX2,iX3,:) &
          = RESHAPE( uPF_K, [hi_P(4)-lo_P(4)+1] )

        ! --- Auxiliary Fluid ---

        CALL ComputePressureFromPrimitive &
               ( uPF_K(:,iPF_D ), uPF_K(:,iPF_E), uPF_K(:,iPF_Ne), &
                 uAF_K(:,iAF_P) )

        CALL ComputeSoundSpeedFromPrimitive &
               ( uPF_K(:,iPF_D ), uPF_K(:,iPF_E), uPF_K(:,iPF_Ne), &
                 uAF_K(:,iAF_Cs) )

        uAF(iX1,iX2,iX3,:) &
          = RESHAPE( uAF_K, [hi_A(4)-lo_A(4)+1] )

      END DO
      END DO
      END DO

    END DO

    CALL amrex_mfiter_destroy( MFI )

  END SUBROUTINE MF_ComputeFromConserved


  SUBROUTINE MF_ComputeTimeStep( nLevels, MF_uGF, MF_uCF, CFL, TimeStep )

    INTEGER,              INTENT(in)  :: &
      nLevels
    TYPE(amrex_multifab), INTENT(in)  :: &
      MF_uGF(0:nlevels), &
      MF_uCF(0:nLevels)
    REAL(amrex_real),     INTENT(in)  :: &
      CFL
    REAL(amrex_real),     INTENT(out) :: &
      TimeStep(0:nLevels)

    INTEGER            :: iLevel
    INTEGER            :: iX1, iX2, iX3, iX_B(3), iX_E(3)
    INTEGER            :: lo_G(4), hi_G(4)
    INTEGER            :: lo_C(4), hi_C(4)
    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI
    REAL(amrex_real), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(amrex_real), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    REAL(amrex_real), ALLOCATABLE         :: G(:,:,:,:,:)
    REAL(amrex_real), ALLOCATABLE         :: U(:,:,:,:,:)

    DO iLevel = 0, nLevels

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = .TRUE. )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCF => MF_uCF(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        lo_G = LBOUND( uGF ); hi_G = UBOUND( uGF )
        lo_C = LBOUND( uCF ); hi_C = UBOUND( uCF )

        iX_B = BX % lo
        iX_E = BX % hi

        ALLOCATE( G(1:nDOFX,iX_B(1):iX_E(1), &
                            iX_B(2):iX_E(2), &
                            iX_B(3):iX_E(3),1:nGF) )
        ALLOCATE( U(1:nDOFX,iX_B(1):iX_E(1), &
                            iX_B(2):iX_E(2), &
                            iX_B(3):iX_E(3),1:nCF) )

        DO iX3 = BX % lo(3), BX % hi(3)
        DO iX2 = BX % lo(2), BX % hi(2)
        DO iX1 = BX % lo(1), BX % hi(1)

           G(1:nDOFX,iX1,iX2,iX3,1:nGF) &
             = RESHAPE( uGF(iX1,iX2,iX3,lo_G(4):hi_G(4)), [ nDOFX, nGF ] )
           U(1:nDOFX,iX1,iX2,iX3,1:nCF) &
             = RESHAPE( uCF(iX1,iX2,iX3,lo_C(4):hi_C(4)), [ nDOFX, nCF ] )

        END DO
        END DO
        END DO


        CALL ComputeTimeStep &
               ( iX_B, iX_E, &
                 G(:,iX_B(1):iX_E(1),iX_B(2):iX_E(2),iX_B(3):iX_E(3),:), &
                 U(:,iX_B(1):iX_E(1),iX_B(2):iX_E(2),iX_B(3):iX_E(3),:), &
                 CFL, TimeStep(iLevel) )

        DEALLOCATE( G )
        DEALLOCATE( U )

      END DO

    END DO

  END SUBROUTINE MF_ComputeTimeStep


END MODULE MF_Euler_UtilitiesModule
