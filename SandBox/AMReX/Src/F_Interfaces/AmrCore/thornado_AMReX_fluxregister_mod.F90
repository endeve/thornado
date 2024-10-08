module thornado_amrex_fluxregister_module

  use iso_c_binding
  use amrex_base_module

  implicit none

  private

  public :: amrex_fluxregister_destroy ! List first to avoid XL compiler bug
  public :: amrex_fluxregister_build

  type, public :: amrex_fluxregister
     logical     :: owner = .false.
     type(c_ptr) :: p     = c_null_ptr
     integer     :: flev  = -1
   contains
     procedure :: fineadd_dg    => amrex_fluxregister_fineadd_dg
     procedure :: crseinit_dg   => amrex_fluxregister_crseinit_dg
     procedure :: reflux_dg     => amrex_fluxregister_reflux_dg
  end type amrex_fluxregister

  interface
     subroutine amrex_fi_new_fluxregister (fr, ba, dm, rr, flev, nc) bind(c)
       import
       implicit none
       type(c_ptr) :: fr
       type(c_ptr), value :: ba, dm
       integer, value :: rr, flev, nc
     end subroutine amrex_fi_new_fluxregister

     subroutine amrex_fi_delete_fluxregister (fr) bind(c)
       import
       implicit none
       type(c_ptr), value :: fr
     end subroutine amrex_fi_delete_fluxregister

     subroutine amrex_fi_fluxregister_fineadd (fr, flxs, scale) bind(c)
       import
       implicit none
       type(c_ptr), value :: fr
       type(c_ptr), intent(in) :: flxs(*)
       real(amrex_real), value :: scale
     end subroutine amrex_fi_fluxregister_fineadd

     subroutine amrex_fi_fluxregister_fineadd_1fab_1dir (fr, fabdata, flo,fhi, dir, boxno, zeroFirst, nfluxes, scale) bind(c)
       import
       implicit none
       type(c_ptr), value :: fr
       type(c_ptr), value,intent(in) :: fabdata
       integer(c_int),intent(in) :: flo(*), fhi(*)
       integer(c_int),intent(IN),value :: dir, boxno, zeroFirst, nfluxes
       real(amrex_real), value :: scale
     end subroutine amrex_fi_fluxregister_fineadd_1fab_1dir

     SUBROUTINE amrex_fi_fluxregister_fineadd_dg &
       ( FluxRegister, SurfaceFluxes, nFields, FaceRatio, &
         nDOFX_X1, nDOFX_X2, nDOFX_X3, &
         nFineX_X1, nFineX_X2, nFineX_X3, &
         WeightsX_X1, WeightsX_X2, WeightsX_X3, &
         vpLX_X1_Refined, vpLX_X2_Refined, vpLX_X3_Refined ) BIND(c)
       IMPORT
       IMPLICIT NONE
       TYPE(C_PTR)     , VALUE      :: &
         FluxRegister, vpLX_X1_Refined, vpLX_X2_Refined, vpLX_X3_Refined
       TYPE(C_PTR)     , INTENT(in) :: SurfaceFluxes(*)
       INTEGER         , VALUE      :: &
         nFields, &
         nDOFX_X1, nDOFX_X2, nDOFX_X3, &
         nFineX_X1, nFineX_X2, nFineX_X3
       REAL(amrex_real), VALUE      :: FaceRatio
       REAL(amrex_real), INTENT(in) :: &
         WeightsX_X1(*), WeightsX_X2(*), WeightsX_X3(*)
     END SUBROUTINE amrex_fi_fluxregister_fineadd_dg

     subroutine amrex_fi_fluxregister_crseinit (fr, flxs, scale) bind(c)
       import
       implicit none
       type(c_ptr), value :: fr
       type(c_ptr), intent(in) :: flxs(*)
       real(amrex_real), value :: scale
     end subroutine amrex_fi_fluxregister_crseinit

     SUBROUTINE amrex_fi_fluxregister_crseinit_dg &
       ( FluxRegister, SurfaceFluxes, nFields, &
         nDOFX_X1, nDOFX_X2, nDOFX_X3, &
         WeightsX_X1, WeightsX_X2, WeightsX_X3 ) BIND(c)
       IMPORT
       IMPLICIT NONE
       TYPE(C_PTR), VALUE             :: &
         FluxRegister
       TYPE(C_PTR),        INTENT(in) :: &
         SurfaceFluxes(*)
       INTEGER    , VALUE, INTENT(in) :: &
         nFields, nDOFX_X1, nDOFX_X2, nDOFX_X3
      REAL(amrex_real)   , INTENT(in) :: &
         WeightsX_X1(*), WeightsX_X2(*), WeightsX_X3(*)
     END SUBROUTINE amrex_fi_fluxregister_crseinit_dg

     subroutine amrex_fi_fluxregister_crseadd (fr, flxs, scale, geom) bind(c)
       import
       implicit none
       type(c_ptr), value :: fr, geom
       type(c_ptr), intent(in) :: flxs(*)
       real(amrex_real), value :: scale
     end subroutine amrex_fi_fluxregister_crseadd

     subroutine amrex_fi_fluxregister_setval (fr, val) bind(c)
       import
       implicit none
       type(c_ptr), value :: fr
       real(amrex_real), value :: val
     end subroutine amrex_fi_fluxregister_setval

     subroutine amrex_fi_fluxregister_reflux (fr, mf, scale, geom) bind(c)
       import
       implicit none
       type(c_ptr), value :: fr, mf, geom
       real(amrex_real), value :: scale
     end subroutine amrex_fi_fluxregister_reflux

     subroutine amrex_fi_fluxregister_reflux_dg &
       ( FluxRegister, MF_G, MF_dU, geom, &
         nDOFX, nDOFX_X1, nDOFX_X2, nDOFX_X3, &
         nFields, iGF_SqrtGm, &
         NodeNumberTableX_X1, NodeNumberTableX_X2, NodeNumberTableX_X3, &
         WeightsX_q, &
         LX_X1_Up, LX_X1_Dn, &
         LX_X2_Up, LX_X2_Dn, &
         LX_X3_Up, LX_X3_Dn, &
         dX1, dX2, dX3 ) bind(c)
       import
       implicit none
       type(c_ptr)     , value :: FluxRegister, MF_G, MF_dU, geom, &
                                  NodeNumberTableX_X1, NodeNumberTableX_X2, &
                                  NodeNumberTableX_X3, &
                                  WeightsX_q, &
                                  LX_X1_Up, LX_X1_Dn, &
                                  LX_X2_Up, LX_X2_Dn, &
                                  LX_X3_Up, LX_X3_Dn
       integer         , value :: nDOFX, nDOFX_X1, nDOFX_X2, nDOFX_X3, &
                                  nFields, iGF_SqrtGm
       real(amrex_real), value :: dX1, dX2, dX3
     end subroutine amrex_fi_fluxregister_reflux_dg

  end interface

contains

  subroutine amrex_fluxregister_build (fr, ba, dm, ref_ratio, fine_lev, ncomp)
    type(amrex_fluxregister) :: fr
    type(amrex_boxarray), intent(in) :: ba
    type(amrex_distromap), intent(in) :: dm
    integer, intent(in) :: ref_ratio, fine_lev, ncomp
    fr%owner = .true.
    fr%flev  = fine_lev
    call amrex_fi_new_fluxregister(fr%p, ba%p, dm%p, ref_ratio, fine_lev, ncomp)
  end subroutine amrex_fluxregister_build

  impure elemental subroutine amrex_fluxregister_destroy (this)
    type(amrex_fluxregister), intent(inout) :: this
    if (this%owner) then
       if (c_associated(this%p)) then
          call amrex_fi_delete_fluxregister(this%p)
       end if
    end if
    this%owner = .false.
    this%p = c_null_ptr
  end subroutine amrex_fluxregister_destroy


  SUBROUTINE amrex_fluxregister_fineadd_dg &
    ( this, SurfaceFluxes, nFields, FaceRatio, &
      nDOFX_X1, nDOFX_X2, nDOFX_X3, &
      nFineX_X1, nFineX_X2, nFineX_X3, &
      WeightsX_X1, WeightsX_X2, WeightsX_X3, &
      vpLX_X1_Refined, vpLX_X2_Refined, vpLX_X3_Refined )

    CLASS(amrex_fluxregister), INTENT(inout) :: &
      this
    TYPE(amrex_multifab)     , INTENT(in)    :: &
      SurfaceFluxes(amrex_spacedim)
    INTEGER                  , INTENT(in)    :: &
      nFields, &
      nDOFX_X1, nDOFX_X2, nDOFX_X3, &
      nFineX_X1, nFineX_X2, nFineX_X3
    REAL(amrex_real)         , INTENT(in)    :: &
      FaceRatio
    REAL(amrex_real)         , INTENT(in)    :: &
      WeightsX_X1(*), WeightsX_X2(*), WeightsX_X3(*)
    TYPE(c_ptr)              , INTENT(in)     :: &
      vpLX_X1_Refined, vpLX_X2_Refined, vpLX_X3_Refined

    INTEGER     :: iDimX
    TYPE(C_PTR) :: MF(amrex_spacedim)

    DO iDimX = 1, amrex_spacedim

       MF(iDimX) = SurfaceFluxes(iDimX) % p

    END DO

    CALL amrex_fi_fluxregister_fineadd_dg &
           ( this%p, mf, nFields, FaceRatio, &
             nDOFX_X1, nDOFX_X2, nDOFX_X3, &
             nFineX_X1, nFineX_X2, nFineX_X3, &
             WeightsX_X1, WeightsX_X2, WeightsX_X3, &
             vpLX_X1_Refined, vpLX_X2_Refined, vpLX_X3_Refined )

  END SUBROUTINE amrex_fluxregister_fineadd_dg


  SUBROUTINE amrex_fluxregister_crseinit_dg &
    ( this, SurfaceFluxes, nFields, &
      nDOFX_X1, nDOFX_X2, nDOFX_X3, &
      WeightsX_X1, WeightsX_X2, WeightsX_X3 )

    CLASS(amrex_fluxregister), INTENT(inout) :: &
      this
    TYPE(amrex_multifab)     , INTENT(in)    :: &
      SurfaceFluxes(amrex_spacedim)
    INTEGER                  , INTENT(in)    :: &
      nFields, nDOFX_X1, nDOFX_X2, nDOFX_X3
    REAL(amrex_real)         , INTENT(in)    :: &
      WeightsX_X1(*), WeightsX_X2(*), WeightsX_X3(*)

    INTEGER     :: iDimX
    TYPE(C_PTR) :: MF(amrex_spacedim)

    DO iDimX = 1, amrex_spacedim

       MF(iDimX) = SurfaceFluxes(iDimX) % p

    END DO

    CALL amrex_fi_fluxregister_crseinit_dg &
           ( this % p, MF, nFields, nDOFX_X1, nDOFX_X2, nDOFX_X3, &
             WeightsX_X1, WeightsX_X2, WeightsX_X3 )

  END SUBROUTINE amrex_fluxregister_crseinit_dg


  subroutine amrex_fluxregister_reflux_dg &
    ( this, MF_G, MF_dU, &
      nDOFX, nDOFX_X1, nDOFX_X2, nDOFX_X3, &
      nFields, iGF_SqrtGm, &
      pNodeNumberTableX_X1, pNodeNumberTableX_X2, pNodeNumberTableX_X3, &
      pWeightsX_q, &
      pLX_X1_Up, pLX_X1_Dn, pLX_X2_Up, pLX_X2_Dn, pLX_X3_Up, pLX_X3_Dn, &
      dX1, dX2, dX3 )
    use amrex_amrcore_module, only : amrex_geom
    class(amrex_fluxregister), intent(inout) :: this
    type(amrex_multifab)     , intent(in)    :: MF_G
    type(amrex_multifab)     , intent(in)    :: MF_dU
    integer                  , intent(in)    :: &
      nDOFX, nDOFX_X1, nDOFX_X2, nDOFX_X3, nFields, iGF_SqrtGm
    type(c_ptr)              , intent(in)    :: &
      pNodeNumberTableX_X1, pNodeNumberTableX_X2, pNodeNumberTableX_X3, &
      pWeightsX_q, &
      pLX_X1_Up, pLX_X1_Dn, pLX_X2_Up, pLX_X2_Dn, pLX_X3_Up, pLX_X3_Dn
    real(amrex_real)         , intent(in)    :: dX1, dX2, dX3

    call amrex_fi_fluxregister_reflux_dg &
           ( this%p, MF_G%p, MF_dU%p, amrex_geom(this%flev-1)%p, &
             nDOFX, nDOFX_X1, nDOFX_X2, nDOFX_X3, nFields, iGF_SqrtGm, &
             pNodeNumberTableX_X1, pNodeNumberTableX_X2, pNodeNumberTableX_X3, &
             pWeightsX_q, &
             pLX_X1_Up, pLX_X1_Dn, pLX_X2_Up, pLX_X2_Dn, pLX_X3_Up, pLX_X3_Dn, &
             dX1, dX2, dX3 )
  end subroutine amrex_fluxregister_reflux_dg

end module thornado_amrex_fluxregister_module
