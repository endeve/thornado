#include <AMReX_Config.H>

MODULE thornado_amrex_fillpatch_module

  USE amrex_base_module
  USE amrex_fillpatch_module

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: thornado_amrex_fillpatch, thornado_amrex_fillcoarsepatch

  INTERFACE thornado_amrex_fillpatch
    MODULE PROCEDURE amrex_fillpatch_single
    MODULE PROCEDURE amrex_fillpatch_dgconservative_two
    MODULE PROCEDURE amrex_fillpatch_dgpointwise_two
  END INTERFACE thornado_amrex_fillpatch

  INTERFACE thornado_amrex_fillcoarsepatch
    MODULE PROCEDURE amrex_fillcoarsepatch_dgconservative
    MODULE PROCEDURE amrex_fillcoarsepatch_dgpointwise
  END INTERFACE thornado_amrex_fillcoarsepatch

  INTERFACE

    ! --- This is just copied from the amrex source code ---
    SUBROUTINE amrex_fi_fillpatch_single &
      ( pMF, Time, sMF, sTime, ns, sComp, dComp, nComp, &
        pGeom, fill ) BIND(c)
      IMPORT
      IMPLICIT NONE
      TYPE(C_PTR)     , VALUE      :: pMF, pGeom
      TYPE(C_PTR)     , INTENT(in) :: sMF(*)
      REAL(amrex_real), VALUE      :: Time
      REAL(amrex_real), INTENT(in) :: sTime(*)
      INTEGER(c_int)  , VALUE      :: sComp, dComp, nComp, ns
      TYPE(C_FUNPTR)  , VALUE      :: fill
    END SUBROUTINE amrex_fi_fillpatch_single


     SUBROUTINE amrex_fi_fillpatch_dgconservative_two &
       ( pMF, pMF_G, Time, &
         pCrseMF, pCrseMF_G, CrseTime, nCrse, &
         pFineMF, pFineMF_G, FineTime, nFine, &
         sComp, dComp, nComp, &
         pCrseGeom, pFineGeom, fpCrseFillPhysBC, fpFineFillPhysBC, &
         RefRatio, interp, pLoBC, pHiBC, nFineV, nDOFX, &
         vpCoarseToFineProjectionMatrix, &
         pre_interp, post_interp ) BIND(c)
       IMPORT
       IMPLICIT NONE
       TYPE(C_PTR)     , VALUE      :: pMF, pMF_G, pCrseGeom, pFineGeom, &
                                       vpCoarseToFineProjectionMatrix
       TYPE(C_PTR)     , INTENT(in) :: pCrseMF(*), pCrseMF_G(*), &
                                       pFineMF(*), pFineMF_G(*), &
                                       pLoBC(*), pHiBC(*)
       TYPE(C_FUNPTR)  , VALUE      :: fpCrseFillPhysBC, fpFineFillPhysBC, &
                                       pre_interp, post_interp
       REAL(amrex_real), VALUE      :: Time
       REAL(amrex_real), INTENT(in) :: CrseTime(*), FineTime(*)
       INTEGER         , VALUE      :: nCrse, nFine, sComp, dComp, nComp, &
                                       RefRatio, interp, nFineV, nDOFX
     END SUBROUTINE amrex_fi_fillpatch_dgconservative_two


     SUBROUTINE amrex_fi_fillpatch_dgpointwise_two &
       ( pMF, Time, &
         pCrseMF, CrseTime, nCrse, &
         pFineMF, FineTime, nFine, &
         sComp, dComp, nComp, &
         pCrseGeom, pFineGeom, fpCrseFillPhysBC, fpFineFillPhysBC, &
         RefRatio, interp, pLoBC, pHiBC, nFineV, nDOFX, &
         vpCoarseToFineProjectionMatrix, &
         pre_interp, post_interp ) BIND(c)
       IMPORT
       IMPLICIT NONE
       TYPE(C_PTR)     , VALUE      :: pMF, pCrseGeom, pFineGeom, &
                                       vpCoarseToFineProjectionMatrix
       TYPE(C_PTR)     , INTENT(in) :: pCrseMF(*), &
                                       pFineMF(*), &
                                       pLoBC(*), pHiBC(*)
       TYPE(C_FUNPTR)  , VALUE      :: fpCrseFillPhysBC, fpFineFillPhysBC, &
                                       pre_interp, post_interp
       REAL(amrex_real), VALUE      :: Time
       REAL(amrex_real), INTENT(in) :: CrseTime(*), FineTime(*)
       INTEGER         , VALUE      :: nCrse, nFine, sComp, dComp, nComp, &
                                       RefRatio, interp, nFineV, nDOFX
     END SUBROUTINE amrex_fi_fillpatch_dgpointwise_two


     SUBROUTINE amrex_fi_fillcoarsepatch_dgconservative &
       ( pMF, pMF_G, Time, &
         pCrseMF, pCrseMF_G, sComp, dComp, nComp, pCrseGeom, pFineGeom, &
         fpCrseFillPhysBC, fpFineFillPhysBC, RefRatio, interp, pLoBC, pHiBC, &
         nFineV, nDOFX, vpCoarseToFineProjectionMatrix, &
         pre_interp, post_interp ) BIND(c)
       IMPORT
       IMPLICIT NONE
       TYPE(C_PTR)     , VALUE :: pMF, pMF_G, pCrseMF, pCrseMF_G, &
                                  pCrseGeom, pFineGeom, &
                                  vpCoarseToFineProjectionMatrix
       TYPE(C_PTR), INTENT(in) :: pLoBC(*), pHiBC(*)
       TYPE(C_FUNPTR)  , VALUE :: fpCrseFillPhysBC, fpFineFillPhysBC, &
                                  pre_interp, post_interp
       REAL(amrex_real), VALUE :: Time
       INTEGER         , VALUE :: sComp, dComp, nComp, RefRatio, interp, &
                                  nFineV, nDOFX
     END SUBROUTINE amrex_fi_fillcoarsepatch_dgconservative


     SUBROUTINE amrex_fi_fillcoarsepatch_dgpointwise &
       ( pMF, Time, &
         pCrseMF, sComp, dComp, nComp, pCrseGeom, pFineGeom, &
         fpCrseFillPhysBC, fpFineFillPhysBC, RefRatio, interp, pLoBC, pHiBC, &
         nFineV, nDOFX, vpCoarseToFineProjectionMatrix, &
         pre_interp, post_interp ) BIND(c)
       IMPORT
       IMPLICIT NONE
       TYPE(C_PTR)     , VALUE :: pMF, pCrseMF, &
                                  pCrseGeom, pFineGeom, &
                                  vpCoarseToFineProjectionMatrix
       TYPE(C_PTR), INTENT(in) :: pLoBC(*), pHiBC(*)
       TYPE(C_FUNPTR)  , VALUE :: fpCrseFillPhysBC, fpFineFillPhysBC, &
                                  pre_interp, post_interp
       REAL(amrex_real), VALUE :: Time
       INTEGER         , VALUE :: sComp, dComp, nComp, RefRatio, interp, &
                                  nFineV, nDOFX
     END SUBROUTINE amrex_fi_fillcoarsepatch_dgpointwise

  END INTERFACE

CONTAINS


  ! --- This is just copied from the amrex source code ---
  SUBROUTINE amrex_fillpatch_single &
    ( MF, OldTime, OldMF, NewTime, NewMF, Geom, fill_physbc, &
      Time, sComp, dComp, nComp)

    TYPE(amrex_multifab), INTENT(inout) :: MF
    TYPE(amrex_multifab), INTENT(in)    :: OldMF, NewMF
    INTEGER             , INTENT(in)    :: sComp, dComp, nComp
    REAL(amrex_real)    , INTENT(in)    :: OldTime, NewTime, Time
    TYPE(amrex_geometry), INTENT(in)    :: Geom
    PROCEDURE(amrex_physbc_proc)        :: fill_physbc

    REAL(amrex_real) :: teps
    REAL(amrex_real) :: sTime(2)
    TYPE(C_PTR)      :: smf(2)
    INTEGER          :: ns

    teps = 1.0e-4_amrex_real * ABS( NewTime - OldTime )
    IF( ABS( Time - NewTime ) .LE. teps )THEN
       ns       = 1
       sMF  (1) = NewMF % p
       sTime(1) = NewTime
    ELSE IF( ABS( Time - OldTime ) .LE. teps )THEN
       ns       = 1
       sMF  (1) = OldMF % p
       sTime(1) = OldTime
    ELSE
       ns       = 2
       sMF  (1) = OldMF % p
       sMF  (2) = NewMF % p
       sTime(1) = OldTime
       sTime(2) = NewTime
    END IF

    ! sComp-1 and dComp-1 because of Fortran index starts with 1
    CALL amrex_fi_fillpatch_single &
           ( MF % p, Time, smF, sTime, ns, sComp-1, dComp-1, nComp, &
             Geom % p, C_FUNLOC( fill_physbc ) )

  END SUBROUTINE amrex_fillpatch_single


  SUBROUTINE amrex_fillpatch_dgconservative_two &
    ( MF, MF_G, &
      OldTimeCrse, OldMFCrse, OldMFCrse_G, NewTimeCrse, NewMFCrse, NewMFCrse_G, &
      GeomCrse, FillPhysBCCrse, &
      OldTimeFine, OldMFFine, OldMFFine_G, NewTimeFine, NewMFFine, NewMFFine_G, &
      GeomFine, FillPhysBCFine, &
      Time, sComp, dComp, nComp, RefRatio, interp, LoBC, HiBC, &
      nFineV, nDOFX, vpCoarseToFineProjectionMatrix, &
      pre_interp, post_interp )

    TYPE(amrex_multifab), INTENT(inout) :: MF
    TYPE(amrex_multifab), INTENT(in)    :: &
      OldMFCrse, OldMFCrse_G, NewMFCrse, NewMFCrse_G, &
      OldMFFine, OldMFFine_G, NewMFFine, NewMFFine_G, MF_G
    INTEGER             , INTENT(in)    :: &
      sComp, dComp, nComp, RefRatio, interp, nFineV, nDOFX
    INTEGER, TARGET     , INTENT(in)    :: &
      LoBC(amrex_spacedim,sComp+nComp-1), HiBC(amrex_spacedim,sComp+nComp-1)
    REAL(amrex_real)    , INTENT(in)    :: &
      OldTimeCrse, NewTimeCrse, OldTimeFine, NewTimeFine, Time
    TYPE(amrex_geometry), INTENT(in)    :: &
      GeomCrse, GeomFine
    TYPE(c_ptr)         , INTENT(in)    :: &
      vpCoarseToFineProjectionMatrix
    PROCEDURE(amrex_physbc_proc)        :: &
      FillPhysBCCrse, FillPhysBCFine
    PROCEDURE(amrex_interp_hook_proc), OPTIONAL :: &
      pre_interp, post_interp

    REAL(amrex_real) :: teps
    REAL(amrex_real) :: CrseTime (2), FineTime (2)
    TYPE(c_ptr)      :: pCrseMF  (2), pFineMF  (2)
    TYPE(c_ptr)      :: pCrseMF_G(2), pFineMF_G(2)
    TYPE(c_ptr)      :: pLoBC(sComp+nComp-1), pHiBC(sComp+nComp-1)
    TYPE(c_funptr)   :: pre_interp_ptr, post_interp_ptr
    INTEGER          :: nCrse, nFine, iComp

    ! Coarse level
    teps = 1.0e-4_amrex_real * ABS( NewTimeCrse - OldTimeCrse )
    IF( ABS( Time - NewTimeCrse ) .LE. teps )THEN

       nCrse        = 1
       pCrseMF  (1) = NewMFCrse   % p
       pCrseMF_G(1) = NewMFCrse_G % p
       CrseTime (1) = NewTimeCrse

    ELSE IF( ABS( Time - OldTimeCrse ) .LE. teps )THEN

       nCrse        = 1
       pCrseMF  (1) = OldMFCrse   % p
       pCrseMF_G(1) = OldMFCrse_G % p
       CrseTime (1) = OldTimeCrse

    ELSE

       nCrse        = 2
       pCrseMF  (1) = OldMFCrse   % p
       pCrseMF  (2) = NewMFCrse   % p
       pCrseMF_G(1) = OldMFCrse_G % p
       pCrseMF_G(2) = NewMFCrse_G % p
       CrseTime (1) = OldTimeCrse
       CrseTime (2) = NewTimeCrse

    END IF

    ! Fine level
    teps = 1.0e-4_amrex_real * ABS( NewTimeFine - OldTimeFine )
    IF( ABS( Time - NewTimeFine ) .LE. teps )THEN

       nFine        = 1
       pFineMF  (1) = NewMFFine   % p
       pFineMF_G(1) = NewMFFine_G % p
       FineTime (1) = NewTimeFine

    ELSE IF( ABS( Time - OldTimeFine ) .LE. teps )THEN

       nFine        = 1
       pFineMF  (1) = OldMFFine   % p
       pFineMF_G(1) = OldMFFine_G % p
       FineTime (1) = OldTimeFine

    ELSE

       nFine        = 2
       pFineMF  (1) = OldMFFine   % p
       pFineMF  (2) = NewMFFine   % p
       pFineMF_G(1) = OldMFFine_G % p
       pFineMF_G(2) = NewMFFine_G % p
       FineTime (1) = OldTimeFine
       FineTime (2) = NewTimeFine

    END IF

    DO iComp = 1, sComp-1

       pLoBC(iComp) = c_null_ptr
       pHiBC(iComp) = c_null_ptr

    END DO

    DO iComp = sComp, sComp+nComp-1

       pLoBC(iComp) = C_LOC( LoBC(1,iComp) )
       pHiBC(iComp) = C_LOC( HiBC(1,iComp) )

    END DO

    pre_interp_ptr = c_null_funptr
    IF( PRESENT( pre_interp ) ) pre_interp_ptr = C_FUNLOC( pre_interp )

    post_interp_ptr = c_null_funptr
    IF( PRESENT( post_interp ) ) post_interp_ptr = C_FUNLOC( post_interp )

    ! sComp-1 and dComp-1 because of Fortran index starts with 1
    CALL amrex_fi_fillpatch_dgconservative_two &
           ( MF % p, MF_G % p, Time, &
             pCrseMF, pCrseMF_G, CrseTime, nCrse, &
             pFineMF, pFineMF_G, FineTime, nFine, &
             sComp-1, dComp-1, nComp, GeomCrse % p, GeomFine % p, &
             C_FUNLOC( FillPhysBCCrse ), C_FUNLOC( FillPhysBCFine ), &
             RefRatio, interp, pLoBC, pHiBC, &
             nFineV, nDOFX, vpCoarseToFineProjectionMatrix, &
             pre_interp_ptr, post_interp_ptr )

  END SUBROUTINE amrex_fillpatch_dgconservative_two


  SUBROUTINE amrex_fillpatch_dgpointwise_two &
    ( MF, &
      OldTimeCrse, OldMFCrse, NewTimeCrse, NewMFCrse, &
      GeomCrse, FillPhysBCCrse, &
      OldTimeFine, OldMFFine, NewTimeFine, NewMFFine, &
      GeomFine, FillPhysBCFine, &
      Time, sComp, dComp, nComp, RefRatio, interp, LoBC, HiBC, &
      nFineV, nDOFX, vpCoarseToFineProjectionMatrix, &
      pre_interp, post_interp )

    TYPE(amrex_multifab), INTENT(inout) :: MF
    TYPE(amrex_multifab), INTENT(in)    :: &
      OldMFCrse, NewMFCrse, &
      OldMFFine, NewMFFine
    INTEGER             , INTENT(in)    :: &
      sComp, dComp, nComp, RefRatio, interp, nFineV, nDOFX
    INTEGER, TARGET     , INTENT(in)    :: &
      LoBC(amrex_spacedim,sComp+nComp-1), HiBC(amrex_spacedim,sComp+nComp-1)
    REAL(amrex_real)    , INTENT(in)    :: &
      OldTimeCrse, NewTimeCrse, OldTimeFine, NewTimeFine, Time
    TYPE(amrex_geometry), INTENT(in)    :: &
      GeomCrse, GeomFine
    TYPE(c_ptr)         , INTENT(in)    :: &
      vpCoarseToFineProjectionMatrix
    PROCEDURE(amrex_physbc_proc)        :: &
      FillPhysBCCrse, FillPhysBCFine
    PROCEDURE(amrex_interp_hook_proc), OPTIONAL :: &
      pre_interp, post_interp

    REAL(amrex_real) :: teps
    REAL(amrex_real) :: CrseTime (2), FineTime (2)
    TYPE(c_ptr)      :: pCrseMF  (2), pFineMF  (2)
    TYPE(c_ptr)      :: pLoBC(sComp+nComp-1), pHiBC(sComp+nComp-1)
    TYPE(c_funptr)   :: pre_interp_ptr, post_interp_ptr
    INTEGER          :: nCrse, nFine, iComp

    ! Coarse level
    teps = 1.0e-4_amrex_real * ABS( NewTimeCrse - OldTimeCrse )
    IF( ABS( Time - NewTimeCrse ) .LE. teps )THEN

       nCrse       = 1
       pCrseMF (1) = NewMFCrse % p
       CrseTime(1) = NewTimeCrse

    ELSE IF( ABS( Time - OldTimeCrse ) .LE. teps )THEN

       nCrse       = 1
       pCrseMF (1) = OldMFCrse % p
       CrseTime(1) = OldTimeCrse

    ELSE

       nCrse       = 2
       pCrseMF (1) = OldMFCrse  % p
       pCrseMF (2) = NewMFCrse  % p
       CrseTime(1) = OldTimeCrse
       CrseTime(2) = NewTimeCrse

    END IF

    ! Fine level
    teps = 1.0e-4_amrex_real * ABS( NewTimeFine - OldTimeFine )
    IF( ABS( Time - NewTimeFine ) .LE. teps )THEN

       nFine       = 1
       pFineMF (1) = NewMFFine % p
       FineTime(1) = NewTimeFine

    ELSE IF( ABS( Time - OldTimeFine ) .LE. teps )THEN

       nFine       = 1
       pFineMF (1) = OldMFFine % p
       FineTime(1) = OldTimeFine

    ELSE

       nFine       = 2
       pFineMF (1) = OldMFFine % p
       pFineMF (2) = NewMFFine % p
       FineTime(1) = OldTimeFine
       FineTime(2) = NewTimeFine

    END IF

    DO iComp = 1, sComp-1

       pLoBC(iComp) = c_null_ptr
       pHiBC(iComp) = c_null_ptr

    END DO

    DO iComp = sComp, sComp+nComp-1

       pLoBC(iComp) = C_LOC( LoBC(1,iComp) )
       pHiBC(iComp) = C_LOC( HiBC(1,iComp) )

    END DO

    pre_interp_ptr = c_null_funptr
    IF( PRESENT( pre_interp ) ) pre_interp_ptr = C_FUNLOC( pre_interp )

    post_interp_ptr = c_null_funptr
    IF( PRESENT( post_interp ) ) post_interp_ptr = C_FUNLOC( post_interp )

    ! sComp-1 and dComp-1 because of Fortran index starts with 1
    CALL amrex_fi_fillpatch_dgpointwise_two &
           ( MF % p, Time, &
             pCrseMF, CrseTime, nCrse, &
             pFineMF, FineTime, nFine, &
             sComp-1, dComp-1, nComp, GeomCrse % p, GeomFine % p, &
             C_FUNLOC( FillPhysBCCrse ), C_FUNLOC( FillPhysBCFine ), &
             RefRatio, interp, pLoBC, pHiBC, &
             nFineV, nDOFX, vpCoarseToFineProjectionMatrix, &
             pre_interp_ptr, post_interp_ptr )

  END SUBROUTINE amrex_fillpatch_dgpointwise_two


  SUBROUTINE amrex_fillcoarsepatch_dgconservative &
    ( MF, MF_G, &
      OldTimeCrse, OldMFCrse, OldMFCrse_G, &
      NewTimeCrse, NewMFCrse, NewMFCrse_G, &
      CrseGeom, FillPhysBCCrse, FineGeom, FillPhysBCFine, &
      Time, nComp, RefRatio, interp, LoBC, HiBC, &
      nFineV, nDOFX, vpCoarseToFineProjectionMatrix, &
      pre_interp, post_interp )

    TYPE(amrex_multifab), INTENT(inout) :: MF, MF_G
    TYPE(amrex_multifab), INTENT(in)    :: OldMFCrse, OldMFCrse_G, &
                                           NewMFCrse, NewMFCrse_G
    INTEGER             , INTENT(in)    :: nComp, RefRatio, &
                                           interp, nFineV, nDOFX
    INTEGER, TARGET     , INTENT(in)    :: &
      LoBC(amrex_spacedim,nComp), HiBC(amrex_spacedim,nComp)
    REAL(amrex_real)    , INTENT(in)    :: OldTimeCrse, NewTimeCrse, Time
    TYPE(amrex_geometry), INTENT(in)    :: CrseGeom, FineGeom
    PROCEDURE(amrex_physbc_proc)        :: FillPhysBCCrse, FillPhysBCFine
    PROCEDURE(amrex_interp_hook_proc), OPTIONAL :: pre_interp
    PROCEDURE(amrex_interp_hook_proc), OPTIONAL :: post_interp
    TYPE(c_ptr)         , INTENT(in)    :: vpCoarseToFineProjectionMatrix

    REAL(amrex_real) :: teps
    TYPE(c_ptr)      :: pCrseMF, pCrseMF_G
    TYPE(c_ptr)      :: pLoBC(nComp), pHiBC(nComp)
    TYPE(c_funptr)   :: pre_interp_ptr, post_interp_ptr
    INTEGER          :: iComp

    INTEGER :: sComp, dComp

    sComp = 1
    dComp = 1

    ! Coarse level
    teps = 1.0e-4_amrex_real * ABS( NewTimeCrse - OldTimeCrse )
    IF( ABS( Time - NewTimeCrse ) .LE. teps )THEN

       pCrseMF   = NewMFCrse   % p
       pCrseMF_G = NewMFCrse_G % p

    ELSE IF( ABS( Time - OldTimeCrse ) .LE. teps )THEN

       pCrseMF   = OldMFCrse   % p
       pCrseMF_G = OldMFCrse_G % p

    ELSE

       pCrseMF = NewMFCrse % p

       CALL amrex_abort( "amrex_fillcoarsepatch: how did this happen?" )

    END IF

    DO iComp = 1, sComp-1

       pLoBC(iComp) = c_null_ptr
       pHiBC(iComp) = c_null_ptr

    END DO

    DO iComp = sComp, sComp+nComp-1

       pLoBC(iComp) = C_LOC( LoBC(1,iComp) )
       pHiBC(iComp) = C_LOC( HiBC(1,iComp) )

    END DO

    pre_interp_ptr = c_null_funptr
    IF( PRESENT( pre_interp ) ) pre_interp_ptr = C_FUNLOC( pre_interp )

    post_interp_ptr = c_null_funptr
    IF( PRESENT( post_interp ) ) post_interp_ptr = C_FUNLOC( post_interp )

    ! sComp-1 and dComp-1 because of Fortran index starts with 1
    CALL amrex_fi_fillcoarsepatch_dgconservative &
           ( MF % p, MF_G % p, Time, pCrseMF, pCrseMF_G, &
             sComp-1, dComp-1, nComp, &
             CrseGeom % p, FineGeom % p, &
             C_FUNLOC( FillPhysBCCrse ), &
             C_FUNLOC( FillPhysBCFine ), &
             RefRatio, interp, pLoBC, pHiBC, &
             nFineV, nDOFX, vpCoarseToFineProjectionMatrix, &
             pre_interp_ptr, post_interp_ptr)

  END SUBROUTINE amrex_fillcoarsepatch_dgconservative


  SUBROUTINE amrex_fillcoarsepatch_dgpointwise &
    ( MF, &
      OldTimeCrse, OldMFCrse, &
      NewTimeCrse, NewMFCrse, &
      CrseGeom, FillPhysBCCrse, FineGeom, FillPhysBCFine, &
      Time, nComp, RefRatio, interp, LoBC, HiBC, &
      nFineV, nDOFX, vpCoarseToFineProjectionMatrix, &
      pre_interp, post_interp )

    TYPE(amrex_multifab), INTENT(inout) :: MF
    TYPE(amrex_multifab), INTENT(in)    :: OldMFCrse, &
                                           NewMFCrse
    INTEGER             , INTENT(in)    :: nComp, RefRatio, &
                                           interp, nFineV, nDOFX
    INTEGER, TARGET     , INTENT(in)    :: &
      LoBC(amrex_spacedim,nComp), HiBC(amrex_spacedim,nComp)
    REAL(amrex_real)    , INTENT(in)    :: OldTimeCrse, NewTimeCrse, Time
    TYPE(amrex_geometry), INTENT(in)    :: CrseGeom, FineGeom
    PROCEDURE(amrex_physbc_proc)        :: FillPhysBCCrse, FillPhysBCFine
    PROCEDURE(amrex_interp_hook_proc), OPTIONAL :: pre_interp
    PROCEDURE(amrex_interp_hook_proc), OPTIONAL :: post_interp
    TYPE(c_ptr)         , INTENT(in)    :: vpCoarseToFineProjectionMatrix

    REAL(amrex_real) :: teps
    TYPE(c_ptr)      :: pCrseMF
    TYPE(c_ptr)      :: pLoBC(nComp), pHiBC(nComp)
    TYPE(c_funptr)   :: pre_interp_ptr, post_interp_ptr
    INTEGER          :: iComp

    INTEGER :: sComp, dComp

    sComp = 1
    dComp = 1

    ! Coarse level
    teps = 1.0e-4_amrex_real * ABS( NewTimeCrse - OldTimeCrse )
    IF( ABS( Time - NewTimeCrse ) .LE. teps )THEN

       pCrseMF = NewMFCrse % p

    ELSE IF( ABS( Time - OldTimeCrse ) .LE. teps )THEN

       pCrseMF = OldMFCrse % p

    ELSE

       pCrseMF = NewMFCrse % p

       CALL amrex_abort( "amrex_fillcoarsepatch: how did this happen?" )

    END IF

    DO iComp = 1, sComp-1

       pLoBC(iComp) = c_null_ptr
       pHiBC(iComp) = c_null_ptr

    END DO

    DO iComp = sComp, sComp+nComp-1

       pLoBC(iComp) = C_LOC( LoBC(1,iComp) )
       pHiBC(iComp) = C_LOC( HiBC(1,iComp) )

    END DO

    pre_interp_ptr = c_null_funptr
    IF( PRESENT( pre_interp ) ) pre_interp_ptr = C_FUNLOC( pre_interp )

    post_interp_ptr = c_null_funptr
    IF( PRESENT( post_interp ) ) post_interp_ptr = C_FUNLOC( post_interp )

    ! sComp-1 and dComp-1 because of Fortran index starts with 1
    CALL amrex_fi_fillcoarsepatch_dgpointwise &
           ( MF % p, Time, pCrseMF, &
             sComp-1, dComp-1, nComp, &
             CrseGeom % p, FineGeom % p, &
             C_FUNLOC( FillPhysBCCrse ), &
             C_FUNLOC( FillPhysBCFine ), &
             RefRatio, interp, pLoBC, pHiBC, &
             nFineV, nDOFX, vpCoarseToFineProjectionMatrix, &
             pre_interp_ptr, post_interp_ptr)

  END SUBROUTINE amrex_fillcoarsepatch_dgpointwise

END MODULE thornado_amrex_fillpatch_module
