MODULE MyRestartModule

  USE iso_c_binding
  USE amrex_fort_module
  USE amrex_base_module 

  ! --- thornado modules ---
  USE ProgramHeaderModule,  ONLY: &
    nDOFX
  USE GeometryFieldsModule, ONLY: &
    nGF
  USE FluidFieldsModule,    ONLY: &
    nCF, nPF, nAF

  IMPLICIT NONE

  INTERFACE

    SUBROUTINE WriteCheckpointFile &
                 ( StepNo, FinestLevel, dt, time, pBA, &
                   pMF_uGF, pMF_uCF, pMF_uPF, pMF_uAF ) BIND(c)
       IMPORT
       IMPLICIT NONE
       INTEGER(c_int),   INTENT(in) :: StepNo(*)
       INTEGER(c_int),   VALUE      :: FinestLevel
       REAL(amrex_real), INTENT(in) :: dt(*), time(*)
       TYPE(c_ptr),      INTENT(in) :: pBA(*)
       TYPE(c_ptr),      INTENT(in) :: pMF_uGF(*)
       TYPE(c_ptr),      INTENT(in) :: pMF_uCF(*)
       TYPE(c_ptr),      INTENT(in) :: pMF_uPF(*)
       TYPE(c_ptr),      INTENT(in) :: pMF_uAF(*)
    END SUBROUTINE WriteCheckpointFile

    SUBROUTINE ReadHeaderAndBoxArrayData &
                 ( FinestLevel, StepNo, dt, time, pBA, pDM ) BIND(c)
      IMPORT
      IMPLICIT NONE
      INTEGER(c_int),   INTENT(out) :: FinestLevel(*)
      INTEGER(c_int),   INTENT(out) :: StepNo(*)
      REAL(amrex_real), INTENT(out) :: dt(*),time(*)
      TYPE(c_ptr),      INTENT(out) :: pBA(*), pDM(*)
    END SUBROUTINE ReadHeaderAndBoxArrayData

    SUBROUTINE ReadMultiFabData( FinestLevel, pMF, iMF ) BIND(c)
      IMPORT
      IMPLICIT NONE
      INTEGER(c_int), VALUE       :: FinestLevel
      TYPE(c_ptr),    INTENT(out) :: pMF(*)
      INTEGER(c_int), VALUE       :: iMF
    END SUBROUTINE ReadMultiFabData

    SUBROUTINE amrex_fi_set_boxarray( iLevel, pBA, amrcore ) BIND(c)
      IMPORT
      IMPLICIT NONE
      TYPE(c_ptr),    VALUE :: pBA
      INTEGER(c_int), VALUE :: iLevel
      TYPE(c_ptr),    VALUE :: amrcore
    END SUBROUTINE amrex_fi_set_boxarray

    SUBROUTINE amrex_fi_set_distromap  (lev, pdm, amrcore) BIND(c)
      IMPORT
      IMPLICIT NONE
      TYPE(c_ptr),    VALUE :: pdm
      INTEGER(c_int), VALUE :: lev
      TYPE(c_ptr),    VALUE :: amrcore
    END SUBROUTINE amrex_fi_set_distromap

    SUBROUTINE amrex_fi_clone_boxarray (bao, bai) BIND(c)
      IMPORT
      IMPLICIT NONE
      TYPE(c_ptr)        :: bao
      TYPE(c_ptr), VALUE :: bai
    END SUBROUTINE amrex_fi_clone_boxarray

    SUBROUTINE amrex_fi_set_finest_level (lev, amrcore) BIND(c)
      IMPORT
      IMPLICIT NONE
      INTEGER(c_int), VALUE :: lev
      TYPE(c_ptr),    VALUE :: amrcore
    END SUBROUTINE amrex_fi_set_finest_level

  END INTERFACE

CONTAINS

  SUBROUTINE ReadCheckpointFile()

    USE amrex_amr_module
    USE MyAmrDataModule, ONLY: &
      t_new, MF_uGF, MF_uCF, MF_uPF, MF_uAF, &
      flux_reg, stepno_vec, do_reflux, dt_vec

    IMPLICIT NONE
    INTEGER               :: iLevel, FinestLevel(1), nComps
    TYPE(c_ptr)           :: pBA(0:amrex_max_level)
    TYPE(c_ptr)           :: pDM(0:amrex_max_level)
    TYPE(c_ptr)           :: MF(0:amrex_max_level)
    TYPE(c_ptr)           :: GF(0:amrex_max_level)
    TYPE(c_ptr)           :: CF(0:amrex_max_level)
    TYPE(c_ptr)           :: PF(0:amrex_max_level)
    TYPE(c_ptr)           :: AF(0:amrex_max_level)
    TYPE(amrex_boxarray)  :: BA(0:amrex_max_level)
    TYPE(amrex_distromap) :: DM(0:amrex_max_level)
    TYPE(c_ptr)           :: amrcore
    TYPE(amrex_box)       :: DOMAIN, DOM
    TYPE(amrex_geometry)  :: GEOM(0:amrex_max_level)

    amrcore = amrex_get_amrcore()

    ! Dummy variables 
    DOMAIN = amrex_box( [0,0,0], [1,1,1] )

    DO iLevel = 0, amrex_max_level
      CALL amrex_geometry_build( GEOM(iLevel), DOMAIN )
    END DO

    DOM = GEOM(0) % DOMAIN

    DO iLevel = 0, amrex_max_level
      CALL amrex_boxarray_build ( BA(iLevel), DOM )
      CALL amrex_distromap_build( DM(iLevel), BA(iLevel) )
    END DO

    pBA(0:amrex_max_level) = BA(0:amrex_max_level) % P
    pDM(0:amrex_max_level) = DM(0:amrex_max_level) % P

    CALL ReadHeaderAndBoxArrayData &
           ( FinestLevel, stepno_vec, dt_vec, t_new, &
             pBA(0:amrex_max_level), pDM(0:amrex_max_level) )

    DO iLevel=0, FinestLevel(1)
      BA(iLevel) = pBA(iLevel)
      DM(iLevel) = pDM(iLevel)
    END DO

    DO iLevel = 0, FinestLevel(1)
      CALL amrex_fi_set_boxarray ( iLevel, BA(iLevel) % P, amrcore )
      CALL amrex_fi_set_distromap( iLevel, DM(iLevel) % P, amrcore )
    END DO

    nComps = 1

    DO iLevel = 0, FinestLevel(1)
      CALL amrex_multifab_build &
             ( MF_uGF(iLevel), BA(iLevel), DM(iLevel), nDOFX * nGF, 0 )
      CALL amrex_multifab_build &
             ( MF_uCF(iLevel), BA(iLevel), DM(iLevel), nDOFX * nCF, 0 )
      CALL amrex_multifab_build &
             ( MF_uPF(iLevel), BA(iLevel), DM(iLevel), nDOFX * nPF, 0 )
      CALL amrex_multifab_build &
             ( MF_uAF(iLevel), BA(iLevel), DM(iLevel), nDOFX * nAF, 0 )

      IF( iLevel .GT. 0 .AND. do_reflux )THEN
        CALL amrex_fluxregister_build &
               ( flux_reg(iLevel), BA(iLevel), DM(iLevel), &
                 amrex_ref_ratio(iLevel-1), iLevel, nComps )
      END IF
    END DO

    GF(0:amrex_max_level) = MF_uGF(0:amrex_max_level) % P
    CF(0:amrex_max_level) = MF_uGF(0:amrex_max_level) % P
    PF(0:amrex_max_level) = MF_uGF(0:amrex_max_level) % P
    AF(0:amrex_max_level) = MF_uGF(0:amrex_max_level) % P
    CALL readmultifabdata( FinestLevel(1), GF(0:amrex_max_level), 0 )
    CALL readmultifabdata( FinestLevel(1), CF(0:amrex_max_level), 1 )
    CALL readmultifabdata( FinestLevel(1), PF(0:amrex_max_level), 2 )
    CALL readmultifabdata( FinestLevel(1), AF(0:amrex_max_level), 3 )

    CALL amrex_fi_set_finest_level( FinestLevel(1), amrcore )
	
  END SUBROUTINE ReadCheckpointFile 

END MODULE MyRestartModule
