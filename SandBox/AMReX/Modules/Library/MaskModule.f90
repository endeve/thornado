MODULE MaskModule

  USE ISO_C_BINDING

  ! --- AMReX Modules ---

  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_boxarray_module, ONLY: &
    amrex_boxarray
  USE amrex_distromap_module, ONLY: &
    amrex_distromap
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_imultifab, &
    amrex_imultifab_build, &
    amrex_imultifab_destroy, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy
  USE amrex_amrcore_module, ONLY: &
    amrex_geom

  ! --- Local Modules ---

  USE InputParsingModule, ONLY: &
    nLevels, &
    swX, &
    UseTiling

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: CreateFineMask
  PUBLIC :: DestroyFineMask
  PUBLIC :: IsLeafElement

  INTEGER :: iLeaf  = 0
  INTEGER :: iTrunk = 1

  INTERFACE

    SUBROUTINE amrex_fi_createfinemask_thornado &
      ( iMF_Mask, CrseBA, CrseDM, FineBA, iCoarse, iFine, swXX, geom ) BIND(c)

        IMPORT
        IMPLICIT NONE

        TYPE(c_ptr)           :: iMF_Mask
        TYPE(c_ptr)   , VALUE :: CrseBA
        TYPE(c_ptr)   , VALUE :: CrseDM
        TYPE(c_ptr)   , VALUE :: FineBA
        INTEGER(c_int), VALUE :: iCoarse
        INTEGER(c_int), VALUE :: iFine
        INTEGER(c_int), VALUE :: swXX
        TYPE(c_ptr)   , VALUE :: geom

    END SUBROUTINE amrex_fi_createfinemask_thornado


    SUBROUTINE amrex_fi_destroyfinemask_thornado( iMF_Mask ) BIND(c)

        IMPORT
        IMPLICIT NONE

        TYPE(c_ptr) :: iMF_Mask

    END SUBROUTINE amrex_fi_destroyfinemask_thornado

  END INTERFACE

CONTAINS


  SUBROUTINE CreateFineMask( iLevel, iMF_Mask, BA, DM )

    INTEGER              , INTENT(in)    :: iLevel
    TYPE(amrex_imultifab), INTENT(inout) :: iMF_Mask
    TYPE(amrex_boxarray) , INTENT(in)    :: BA(0:)
    TYPE(amrex_distromap), INTENT(in)    :: DM(0:)

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX
    INTEGER, CONTIGUOUS, POINTER :: Mask(:,:,:,:)

    INTEGER :: iX1, iX2, iX3
    INTEGER :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)

    IF( nLevels .GT. 1 .AND. iLevel .LT. nLevels-1 )THEN

      ! This uses swX(1) because amrex uses IntVect(int), which creates
      ! a vector of dimension `nDimsX` and fills each value with `int`
      iMF_Mask % owner = .TRUE.
      iMF_Mask % nc    = 1
      CALL amrex_fi_createfinemask_thornado &
             ( iMF_Mask % p, BA(iLevel) % p, DM(iLevel) % p, BA(iLevel+1) % p, &
               iLeaf, iTrunk, swX(1), amrex_geom(iLevel) % p )

      ! --- Fix physical boundary ghost cells ---

      CALL amrex_mfiter_build( MFI, iMF_Mask, tiling = UseTiling )

      DO WHILE( MFI % next() )

        Mask => iMF_Mask % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = iX_B0 - swX
        iX_E1 = iX_E0 + swX

        ! --- X1 ---

        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)

          Mask(iX_B1(1),iX2,iX3,1) = Mask(iX_B0(1),iX2,iX3,1)
          Mask(iX_E1(1),iX2,iX3,1) = Mask(iX_E0(1),iX2,iX3,1)

        END DO
        END DO

        ! --- X2 ---

        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX1 = iX_B0(1), iX_E0(1)

          Mask(iX1,iX_B1(2),iX3,1) = Mask(iX1,iX_B0(2),iX3,1)
          Mask(iX1,iX_E1(2),iX3,1) = Mask(iX1,iX_E0(2),iX3,1)

        END DO
        END DO

        ! --- X3 ---

        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)

          Mask(iX1,iX2,iX_B1(3),1) = Mask(iX1,iX2,iX_B0(3),1)
          Mask(iX1,iX2,iX_E1(3),1) = Mask(iX1,iX2,iX_E0(3),1)

        END DO
        END DO

      END DO ! WHILE( MFI % next() )

      CALL amrex_mfiter_destroy( MFI )

    ELSE

      CALL amrex_imultifab_build( iMF_Mask, BA(iLevel), DM(iLevel), 1, swX )
      CALL iMF_Mask % SetVal( iLeaf )

    END IF

  END SUBROUTINE CreateFineMask


  SUBROUTINE DestroyFineMask( iLevel, iMF_Mask )

    INTEGER              , INTENT(in)    :: iLevel
    TYPE(amrex_imultifab), INTENT(inout) :: iMF_Mask

    CALL amrex_fi_destroyfinemask_thornado( iMF_Mask % p )
    CALL amrex_imultifab_destroy          ( iMF_Mask )

  END SUBROUTINE DestroyFineMask


  LOGICAL FUNCTION IsLeafElement( iX )

    INTEGER, INTENT(in) :: iX

    IF( iX .EQ. iLeaf )THEN
      IsLeafElement = .TRUE.
    ELSE
      IsLeafElement = .FALSE.
    END IF

    RETURN
  END FUNCTION IsLeafElement


END MODULE MaskModule
