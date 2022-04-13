MODULE AMReX_BoundaryConditionsModule

  USE amrex_fort_module, ONLY: &
    amrex_spacedim
  USE amrex_bc_types_module, ONLY: &
    amrex_bc_int_dir

  IMPLICIT NONE
  PRIVATE

  ! Periodic BCs. Second dimension is number of components
  INTEGER, PUBLIC, SAVE :: lo_bc(amrex_spacedim,1) = amrex_bc_int_dir
  INTEGER, PUBLIC, SAVE :: hi_bc(amrex_spacedim,1) = amrex_bc_int_dir

END MODULE AMReX_BoundaryConditionsModule

