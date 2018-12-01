MODULE InitializationModule

  ! --- AMReX Modules ---

  USE amrex_base_module

  ! --- thornado Modules ---

  USE ProgramHeaderModule,     ONLY: &
    nDOFX, nX, nNodesX
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE MeshModule,              ONLY: &
    MeshType,                        &
    CreateMesh,                      &
    DestroyMesh,                     &
    NodeCoordinate
  USE FluidFieldsModule,       ONLY: &
    nCF

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFields

CONTAINS


  SUBROUTINE InitializeFields( ProgramName, GEOM, MF_uCF )

    CHARACTER(LEN=*),     INTENT(in   ) :: ProgramName
    TYPE(amrex_geometry), INTENT(in   ) :: GEOM
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF

    IF( amrex_parallel_ioprocessor() )THEN
      WRITE(*,*)
      WRITE(*,'(A4,A,A)') '', 'Initializing: ', TRIM( ProgramName )
    END IF

    SELECT CASE ( TRIM( ProgramName ) )

      CASE ( 'IsentropicVortex' )

        CALL InitializeFields_IsentropicVortex( ProgramName, GEOM, MF_uCF )

      CASE DEFAULT

        IF( amrex_parallel_ioprocessor() )THEN
          WRITE(*,*)
          WRITE(*,'(A4,A,A)') '', 'Unknown Program: ', TRIM( ProgramName )
        END IF

    END SELECT

  END SUBROUTINE InitializeFields


  SUBROUTINE InitializeFields_IsentropicVortex( ProgramName, GEOM, MF_uCF )

    CHARACTER(LEN=*),     INTENT(in   ) :: ProgramName
    TYPE(amrex_geometry), INTENT(in   ) :: GEOM
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF

    INTEGER            :: iDim
    INTEGER            :: iX1, iX2, iX3
    INTEGER            :: iNodeX
    INTEGER            :: iNodeX1, iNodeX2, iNodeX3
    INTEGER            :: lo(4), hi(4)
    REAL(amrex_real)   :: X1, X2, X3
    REAL(amrex_real)   :: uCF_K(nDOFX,nCF)
    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI
    TYPE(MeshType)     :: MeshX(3)
    REAL(amrex_real), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    DO iDim = 1, 3

      CALL CreateMesh &
             ( MeshX(iDim), nX(iDim), nNodesX(iDim), 0, &
               amrex_problo(iDim), amrex_probhi(iDim) )

    END DO

    CALL amrex_mfiter_build( MFI, MF_uCF, tiling = .TRUE. )

    DO WHILE( MFI % next() )

      uCF => MF_uCF % DataPtr( MFI )

      BX = MFI % tilebox()
      lo = LBOUND( uCF )
      hi = UBOUND( uCF )

      PRINT*, "SHAPE( uCF ) = ", SHAPE( uCF )

      PRINT*, "lo = ", lo, BX % lo
      PRINT*, "hi = ", hi, BX % hi

      PRINT*, "amrex_problo = ", amrex_problo
      PRINT*, "amrex_probhi = ", amrex_probhi

      DO iX3 = BX % lo(3), BX % hi(3)
      DO iX2 = BX % lo(2), BX % hi(2)
      DO iX1 = BX % lo(1), BX % hi(1)

        DO iNodeX = 1, nDOFX

          iNodeX1 = NodeNumberTableX(1,iNodeX)
          iNodeX2 = NodeNumberTableX(2,iNodeX)
          iNodeX3 = NodeNumberTableX(3,iNodeX)

          X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
          X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )
          X3 = NodeCoordinate( MeshX(3), iX3, iNodeX3 )



        END DO

        uCF_K(1:nDOFX,1:nCF) &
          = RESHAPE( uCF(iX1,iX2,iX3,lo(4):hi(4)), [ nDOFX, nCF ] )

        uCF(iX1,iX2,iX3,lo(4):hi(4)) &
          = RESHAPE( uCF_K(1:nDOFX,1:nCF), [ hi(4) - lo(4) + 1 ] )

      END DO
      END DO
      END DO

    END DO

    CALL amrex_mfiter_destroy( MFI )

    DO iDim = 1, 3

      CALL DestroyMesh( MeshX(iDim) )

    END DO

  END SUBROUTINE InitializeFields_IsentropicVortex


END MODULE InitializationModule
