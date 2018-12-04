MODULE InputOutputModuleAMReX

  ! --- AMReX Modules ---

  USE amrex_base_module

  ! --- thornado Modules ---

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    nDOFX
  USE ReferenceElementModuleX, ONLY: &
    WeightsX_q
  USE GeometryFieldsModule, ONLY: &
    nGF

  IMPLICIT NONE
  PRIVATE

  INTEGER :: PlotFileNumber = 0

  PUBLIC :: WriteFieldsAMReX_PlotFile

CONTAINS


  SUBROUTINE WriteFieldsAMReX_PlotFile &
    ( Time, GEOM, MF_uGF_Option, MF_uCF_Option, MF_uPF_Option, MF_uAF_Option )

    REAL(DP),             INTENT(in)           :: Time
    TYPE(amrex_geometry), INTENT(in)           :: GEOM
    TYPE(amrex_multifab), INTENT(in), OPTIONAL :: MF_uGF_Option
    TYPE(amrex_multifab), INTENT(in), OPTIONAL :: MF_uCF_Option
    TYPE(amrex_multifab), INTENT(in), OPTIONAL :: MF_uPF_Option
    TYPE(amrex_multifab), INTENT(in), OPTIONAL :: MF_uAF_Option

    LOGICAL :: WriteGF
    LOGICAL :: WriteFF_C
    INTEGER :: MF_nComp
    TYPE(amrex_multifab) :: MF_PF

    ! --- Only for debugging ---
    LOGICAL :: DEBUG = .FALSE.
    TYPE(amrex_mfiter) :: MFI
    REAL(amrex_real), CONTIGUOUS, POINTER :: u(:,:,:,:)
    INTEGER :: iX1, iX2, iX3, iComp


    IF( PRESENT( MF_uGF_Option ) )THEN
      WriteGF  = .TRUE.
    ELSE
      WriteGF  = .FALSE.
    END IF

    IF( amrex_parallel_ioprocessor() )THEN

      WRITE(*,*)
      WRITE(*,'(A4,A18,I9.8)') '', 'Writing PlotFile: ', PlotFileNumber

    END IF

    IF( WriteGF )THEN

      CALL amrex_multifab_build &
             ( MF_PF, MF_uGF_Option % BA, MF_uGF_Option % DM, nGF, 0 )

      CALL MF_ComputeElementAverage( nGF, MF_uGF_Option, MF_PF )

      IF( DEBUG )THEN

        CALL amrex_mfiter_build( MFI, MF_uGF_Option, tiling = .TRUE. )

        u => MF_PF % DataPtr( MFI )
        iX1 = 1
        iX2 = 1
        iX3 = 1
        iComp = 8
        WRITE(*,*) 'iX1, iX2, iX3: ', iX1, iX2, iX3
        WRITE(*,*) 'iComp: ', iComp
        WRITE(*,*) 'Cell-average: ', u(iX1,iX2,iX3,iComp)

      END IF

      CALL amrex_multifab_destroy( MF_PF )

    END IF

  END SUBROUTINE WriteFieldsAMReX_PlotFile


  SUBROUTINE MF_ComputeElementAverage( nComp, MF, MF_A )

    INTEGER,              INTENT(in   ) :: nComp
    TYPE(amrex_multifab), INTENT(in   ) :: MF
    TYPE(amrex_multifab), INTENT(inout) :: MF_A

    INTEGER                               :: iX1, iX2, iX3, iComp
    INTEGER                               :: lo(4), hi(4)
    REAL(amrex_real)                      :: u_K(nDOFX,nComp)
    TYPE(amrex_box)                       :: BX
    TYPE(amrex_mfiter)                    :: MFI
    REAL(amrex_real), CONTIGUOUS, POINTER :: u(:,:,:,:), u_A(:,:,:,:)

    LOGICAL :: DEBUG = .FALSE.

    CALL amrex_mfiter_build( MFI, MF, tiling = .TRUE. )

    DO WHILE( MFI % next() )

      u   => MF   % DataPtr( MFI )
      IF( DEBUG ) &
        WRITE(*,*) 'SHAPE( u   ): ', SHAPE( u )

      u_A => MF_A % DataPtr( MFI )
      IF( DEBUG ) &
        WRITE(*,*) 'SHAPE( u_A ): ', SHAPE( u_A )

      BX = MFI % tilebox()

      lo = LBOUND( u ); hi = UBOUND( u )

      DO iX3 = BX % lo(3), BX % hi(3)
        DO iX2 = BX % lo(2), BX % hi(2)
          DO iX1 = BX % lo(1), BX % hi(1)

            IF( DEBUG ) &
              WRITE(*,*) 'iX1, iX2, iX3 = ', iX1, iX2, iX3

            u_K(1:nDOFX,1:nComp) &
              = RESHAPE( u(iX1,iX2,iX3,lo(4):hi(4)), [ nDOFX, nComp ] )

            ! --- Compute cell-average ---
            DO iComp = 1, nComp

              IF( DEBUG )THEN
                WRITE(*,*) 'iComp: ', iComp
                WRITE(*,*) 'Nodal values:', u_K(1:nDOFX,iComp)
              END IF

              u_A(iX1,iX2,iX3,iComp) = SUM( WeightsX_q * u_K(:,iComp) )

              IF( DEBUG )THEN
                WRITE(*,*) 'Cell-average:', u_A(iX1,iX2,iX3,iComp)
                WRITE(*,*)
              END IF

            END DO

          END DO
        END DO
      END DO

    END DO

    CALL amrex_mfiter_destroy( MFI )

  END SUBROUTINE MF_ComputeElementAverage


END MODULE InputOutputModuleAMReX
