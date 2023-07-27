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
  USE amrex_geometry_module, ONLY: &
    amrex_is_periodic

  ! --- Local Modules ---

  USE InputParsingModule, ONLY: &
    nLevels, &
    swX, &
    UseTiling

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: CreateFineMask
  PUBLIC :: DestroyFineMask
  PUBLIC :: IsNotLeafElement

  ! https://amrex-codes.github.io/amrex/docs_html/Basics.html#fine-mask
  INTEGER, PARAMETER :: iLeaf    = 0
  INTEGER, PARAMETER :: iNotLeaf = 1

  ! https://amrex-codes.github.io/amrex/docs_html/Basics.html#point-mask

  ! Ghost elements covered by valid elements
  INTEGER, PARAMETER :: iCovered = 2
  ! Ghost elements NOT covered by valid elements
  INTEGER, PARAMETER :: iNotCovered = 3
  ! Outside physical domain
  INTEGER, PARAMETER :: iPhysicalBoundary = 4
  ! Interior elements (i.e., valid elements)
  INTEGER, PARAMETER :: iInterior = 5

  INTERFACE

    SUBROUTINE amrex_fi_createfinemask_thornado &
      ( iMF_FineMask, CrseBA, CrseDM, FineBA, &
        iCoarse, iFine, swXX, geom ) BIND(c)

        IMPORT
        IMPLICIT NONE

        TYPE(c_ptr)           :: iMF_FineMask
        TYPE(c_ptr)   , VALUE :: CrseBA
        TYPE(c_ptr)   , VALUE :: CrseDM
        TYPE(c_ptr)   , VALUE :: FineBA
        INTEGER(c_int), VALUE :: iCoarse
        INTEGER(c_int), VALUE :: iFine
        INTEGER(c_int), VALUE :: swXX
        TYPE(c_ptr)   , VALUE :: geom

    END SUBROUTINE amrex_fi_createfinemask_thornado


    SUBROUTINE amrex_fi_destroyfinemask_thornado( iMF_FineMask ) BIND(c)

        IMPORT
        IMPLICIT NONE

        TYPE(c_ptr) :: iMF_FineMask

    END SUBROUTINE amrex_fi_destroyfinemask_thornado


    SUBROUTINE amrex_fi_createpointmask_thornado &
      ( iMF_PointMask, Geom, &
        iCovered, iNotCovered, iPhysicalBoundary, iInterior ) BIND(c)

        IMPORT
        IMPLICIT NONE

        TYPE(c_ptr)           :: iMF_PointMask
        TYPE(c_ptr)   , VALUE :: Geom
        INTEGER(c_int), VALUE :: iCovered
        INTEGER(c_int), VALUE :: iNotCovered
        INTEGER(c_int), VALUE :: iPhysicalBoundary
        INTEGER(c_int), VALUE :: iInterior

    END SUBROUTINE amrex_fi_createpointmask_thornado


  END INTERFACE

CONTAINS


  SUBROUTINE CreateFineMask &
    ( iLevel, iMF_FineMask, BA, DM )

    INTEGER              , INTENT(in)    :: iLevel ! Coarse level
    TYPE(amrex_imultifab), INTENT(inout) :: iMF_FineMask
    TYPE(amrex_boxarray) , INTENT(in)    :: BA(0:)
    TYPE(amrex_distromap), INTENT(in)    :: DM(0:)

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    INTEGER, CONTIGUOUS, POINTER :: FineMask(:,:,:,:)

    INTEGER :: iX1, iX2, iX3
    INTEGER :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)

    IF( nLevels .GT. 1 .AND. iLevel .LT. nLevels-1 )THEN

      ! This uses swX(1) because amrex uses IntVect(int), which creates
      ! a vector of dimension `nDimsX` and fills each value with `int`
      iMF_FineMask % owner = .TRUE.
      iMF_FineMask % nc    = 1
      CALL amrex_fi_createfinemask_thornado &
             ( iMF_FineMask % p, BA(iLevel) % p, &
               DM(iLevel) % p, BA(iLevel+1) % p, &
               iLeaf, iNotLeaf, swX(1), amrex_geom(iLevel) % p )

#if defined( THORNADO_OMP )
      !$OMP PARALLEL &
      !$OMP PRIVATE( MFI, BX, FineMask, iX_B0, iX_E0, iX_B1, iX_E1 )
#endif

      CALL amrex_mfiter_build( MFI, iMF_FineMask, tiling = UseTiling )

      DO WHILE( MFI % next() )

        FineMask => iMF_FineMask % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = iX_B0 - swX
        iX_E1 = iX_E0 + swX

        ! --- Fix physical boundary ghost cells ---

        ! --- X1 ---

        IF( .NOT. amrex_is_periodic(1) )THEN

          IF( iX_B0(1) .EQ. amrex_geom(iLevel) % domain % lo( 1 ) )THEN

            DO iX3 = iX_B0(3), iX_E0(3)
            DO iX2 = iX_B0(2), iX_E0(2)

              IF( IsNotLeafElement( FineMask(iX_B0(1),iX2,iX3,1) ) ) &
                FineMask(iX_B1(1),iX2,iX3,1) = iNotLeaf

            END DO
            END DO

          END IF

          IF( iX_E0(1) .EQ. amrex_geom(iLevel) % domain % hi( 1 ) )THEN

            DO iX3 = iX_B0(3), iX_E0(3)
            DO iX2 = iX_B0(2), iX_E0(2)

              IF( IsNotLeafElement( FineMask(iX_E0(1),iX2,iX3,1) ) ) &
                FineMask(iX_E1(1),iX2,iX3,1) = iNotLeaf

            END DO
            END DO

          END IF

        END IF ! .NOT. amrex_is_periodic(1)

        ! --- X2 ---

        IF( iX_B0(2) .NE. iX_E0(2) )THEN

          IF( .NOT. amrex_is_periodic(2) )THEN

            IF( iX_B0(2) .EQ. amrex_geom(iLevel) % domain % lo( 2 ) )THEN

              DO iX3 = iX_B0(3), iX_E0(3)
              DO iX1 = iX_B0(1), iX_E0(1)

                IF( IsNotLeafElement( FineMask(iX1,iX_B0(2),iX3,1) ) ) &
                  FineMask(iX1,iX_B1(2),iX3,1) = iNotLeaf

              END DO
              END DO

            END IF

            IF( iX_E0(2) .EQ. amrex_geom(iLevel) % domain % hi( 2 ) )THEN

              DO iX3 = iX_B0(3), iX_E0(3)
              DO iX1 = iX_B0(1), iX_E0(1)

                IF( IsNotLeafElement( FineMask(iX1,iX_E0(2),iX3,1) ) ) &
                  FineMask(iX1,iX_E1(2),iX3,1) = iNotLeaf

              END DO
              END DO

            END IF

          END IF ! .NOT. amrex_is_periodic(2)

        END IF ! iX_B0(2) .NE. iX_E0(2)

        ! --- X3 ---

        IF( iX_B0(3) .NE. iX_E0(3) )THEN

          IF( .NOT. amrex_is_periodic(3) )THEN

            IF( iX_B0(3) .EQ. amrex_geom(iLevel) % domain % lo( 3 ) )THEN

              DO iX2 = iX_B0(2), iX_E0(2)
              DO iX1 = iX_B0(1), iX_E0(1)

                IF( IsNotLeafElement( FineMask(iX1,iX2,iX_B0(3),1) ) ) &
                  FineMask(iX1,iX2,iX_B1(3),1) = iNotLeaf

              END DO
              END DO

            END IF

            IF( iX_E0(3) .EQ. amrex_geom(iLevel) % domain % hi( 3 ) )THEN

              DO iX2 = iX_B0(2), iX_E0(2)
              DO iX1 = iX_B0(1), iX_E0(1)

                IF( IsNotLeafElement( FineMask(iX1,iX2,iX_E0(3),1) ) ) &
                  FineMask(iX1,iX2,iX_E1(3),1) = iNotLeaf

              END DO
              END DO

            END IF

          END IF ! .NOT. amrex_is_periodic(3)

        END IF ! iX_B0(3) .NE. iX_E0(3)

      END DO ! WHILE( MFI % next() )

      CALL amrex_mfiter_destroy( MFI )

#if defined( THORNADO_OMP )
      !$OMP END PARALLEL
#endif

    ELSE

      CALL amrex_imultifab_build( iMF_FineMask, BA(iLevel), DM(iLevel), 1, swX )
      CALL iMF_FineMask % SetVal( iLeaf )

    END IF

  END SUBROUTINE CreateFineMask


  SUBROUTINE DestroyFineMask( iMF_FineMask )

    TYPE(amrex_imultifab), INTENT(inout) :: iMF_FineMask

    CALL amrex_imultifab_destroy( iMF_FineMask )

  END SUBROUTINE DestroyFineMask


  SUBROUTINE CreatePointMask( iLevel, iMF_PointMask, BA, DM )

    INTEGER              , INTENT(in)    :: iLevel
    TYPE(amrex_imultifab), INTENT(inout) :: iMF_PointMask
    TYPE(amrex_boxarray) , INTENT(in)    :: BA(0:)
    TYPE(amrex_distromap), INTENT(in)    :: DM(0:)

    CALL amrex_imultifab_build &
           ( iMF_PointMask, BA(iLevel), DM(iLevel), 1, swX )

    CALL amrex_fi_createpointmask_thornado( iMF_PointMask % p, &
                                            amrex_geom(iLevel) % p, &
                                            iCovered, &
                                            iNotCovered, &
                                            iPhysicalBoundary, &
                                            iInterior )

  END SUBROUTINE CreatePointMask


  SUBROUTINE DestroyPointMask( iMF_PointMask )

    TYPE(amrex_imultifab), INTENT(inout) :: iMF_PointMask

    CALL amrex_imultifab_destroy( iMF_PointMask )

  END SUBROUTINE DestroyPointMask


  LOGICAL FUNCTION IsNotLeafElement( Element )

    INTEGER, INTENT(in) :: Element

    IF( Element .EQ. iNotLeaf )THEN
      IsNotLeafElement = .TRUE.
    ELSE
      IsNotLeafElement = .FALSE.
    END IF

    RETURN
  END FUNCTION IsNotLeafElement


  LOGICAL FUNCTION IsCoveredElement( Element )

    INTEGER, INTENT(in) :: Element

    IF( Element .EQ. iCovered )THEN
      IsCoveredElement = .TRUE.
    ELSE
      IsCoveredElement = .FALSE.
    END IF

    RETURN
  END FUNCTION IsCoveredElement


END MODULE MaskModule
