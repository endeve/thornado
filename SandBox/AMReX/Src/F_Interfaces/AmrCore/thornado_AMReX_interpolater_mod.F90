module thornado_amrex_interpolater_module
  use amrex_interpolater_module
  implicit none
  ! THIS MUST BE CONSISTENT WITH amrex_fi_fillpatch_two in AMReX_fillpatch_fi.cpp!!!
  integer, parameter :: amrex_interp_dg = 0
  integer, parameter :: amrex_interp_cg = 1
end module thornado_amrex_interpolater_module
