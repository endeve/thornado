#ifndef THORNADO_AMREX_INTERPOLATER_H_
#define THORNADO_AMREX_INTERPOLATER_H_
#include <AMReX_Config.H>

#include <AMReX_Extension.H>
#include <AMReX_GpuControl.H>
#include <AMReX_InterpBase.H>

namespace amrex {

class Geometry;
class FArrayBox;
class IArrayBox;

/**
* \brief Virtual base class for interpolaters.
*
* Specifies interpolater interface for coarse-to-fine interpolation in space.
*/
class thornadoInterpolater
    : public InterpBase
{
public:

    /**
    * \brief Coarse to fine interpolation in space.
    *
    * \param CrseFab
    * \param CrseFab_G
    * \param FineFab
    * \param FineFab_G
    * \param nComp
    * \param fine_region
    * \param RefRatio
    * \param nDOFX
    * \param CoarseToFineProjectionMatrix
    * \param runon
    */
    virtual void interpConservative
      ( const FArrayBox &  /* CrseFab                      */ ,
        const FArrayBox &  /* CrseFab_G                    */ ,
        FArrayBox       &  /* FineFab                      */ ,
        const FArrayBox &  /* FineFab_G                    */ ,
        int                /* nComp                        */ ,
        const Box       &  /* fine_region                  */ ,
        const IntVect   &  /* RefRatio                     */ ,
        int                /* nDOFX                        */ ,
        Array4<Real const> /* CoarseToFineProjectionMatrix */ ,
        RunOn              /* runon                        */  ){}

    /**
    * \brief Coarse to fine interpolation in space.
    *
    * \param CrseFab
    * \param CrseFab_G
    * \param FineFab
    * \param FineFab_G
    * \param nComp
    * \param fine_region
    * \param RefRatio
    * \param nDOFX
    * \param CoarseToFineProjectionMatrix
    * \param runon
    */
    virtual void interpPointWise
      ( const FArrayBox &  /* CrseFab                      */ ,
        FArrayBox       &  /* FineFab                      */ ,
        int                /* nComp                        */ ,
        const Box       &  /* fine_region                  */ ,
        const IntVect   &  /* RefRatio                     */ ,
        int                /* nDOFX                        */ ,
        Array4<Real const> /* CoarseToFineProjectionMatrix */ ,
        RunOn              /* runon                        */  ){}

    /**
    * \brief Coarse to fine interpolation in space for face-based data.
    *
    * \param crse
    * \param crse_comp
    * \param fine
    * \param fine_comp
    * \param ncomp
    * \param fine_region
    * \param fine_values
    * \param fine_known
    * \param ratio
    * \param crse_geom
    * \param fine_geom
    * \param bcr
    * \param bccomp
    * \param runon
    */
    virtual void interp_face (const FArrayBox& /*crse*/,
                              const int        /*crse_comp*/,
                              FArrayBox&       /*fine*/,
                              const int        /*fine_comp*/,
                              const int        /*ncomp*/,
                              const Box&       /*fine_region*/,
                              const IntVect&   /*ratio*/,
                              const IArrayBox& /*solve_mask*/,
                              const Geometry&  /*crse_geom*/,
                              const Geometry&  /*fine_geom*/,
                              Vector<BCRec> const & /*bcr*/,
                              const int        /*bccomp*/,
                              RunOn            /*runon*/)
    { amrex::Abort("The version of this Interpolater for face-based data is not implemented or does not apply. Call 'interp' instead."); }

};


class DGInterp
    :
    public thornadoInterpolater
{
public:

    /**
    * \brief Returns coarsened box given fine box and refinement ratio.
    *
    * \param fine
    * \param ratio
    */
    Box CoarseBox (const Box& fine,
                   int        ratio) override;

    /**
    * \brief Returns coarsened box given fine box and refinement ratio.
    *
    * \param fine
    * \param ratio
    */
    Box CoarseBox (const Box&     fine,
                   const IntVect& ratio) override;

    /**
    * \brief Coarse to fine interpolation in space.
    *
    * \param CrseFab
    * \param CrseFab_G
    * \param FineFab
    * \param FineFab_G
    * \param nComp
    * \param fine_region
    * \param RefRatio
    * \param nDOFX
    * \param CoarseToFineProjectionMatrix
    * \param runon
    */
    void interpConservative
           ( const FArrayBox &  CrseFab                     ,
             const FArrayBox &  CrseFab_G                   ,
             FArrayBox       &  FineFab                     ,
             const FArrayBox &  FineFab_G                   ,
             int                nComp                       ,
             const Box       &  fine_region                 ,
             const IntVect   &  RefRatio                    ,
             int                nDOFX                       ,
             Array4<Real const> CoarseToFineProjectionMatrix,
             RunOn              runon                        ) override;

    /**
    * \brief Coarse to fine interpolation in space.
    *
    * \param CrseFab
    * \param FineFab
    * \param nComp
    * \param fine_region
    * \param RefRatio
    * \param nDOFX
    * \param CoarseToFineProjectionMatrix
    * \param runon
    */
    void interpPointWise
           ( const FArrayBox &  CrseFab                     ,
             FArrayBox       &  FineFab                     ,
             int                nComp                       ,
             const Box       &  fine_region                 ,
             const IntVect   &  RefRatio                    ,
             int                nDOFX                       ,
             Array4<Real const> CoarseToFineProjectionMatrix,
             RunOn              runon                        ) override;
};
/**
* \brief CG interpolation on nodal data with cell.
*/

class CGInterp
    :
    public thornadoInterpolater
{
public:

    /**
    * \brief Returns coarsened box given fine box and refinement ratio.
    *
    * \param fine
    * \param ratio
    */
    Box CoarseBox (const Box& fine,
                   int        ratio) override;

    /**
    * \brief Returns coarsened box given fine box and refinement ratio.
    *
    * \param fine
    * \param ratio
    */
    Box CoarseBox (const Box&     fine,
                   const IntVect& ratio) override;

    /**
    * \brief Coarse to fine interpolation in space.
    *
    * \param CrseFab
    * \param CrseFab_G
    * \param FineFab
    * \param FineFab_G
    * \param nComp
    * \param fine_region
    * \param RefRatio
    * \param nDOFX
    * \param CoarseToFineProjectionMatrix
    * \param runon
    */
    void interpConservative
           ( const FArrayBox &  CrseFab                     ,
             const FArrayBox &  CrseFab_G                   ,
             FArrayBox       &  FineFab                     ,
             const FArrayBox &  FineFab_G                   ,
             int                nComp                       ,
             const Box       &  fine_region                 ,
             const IntVect   &  RefRatio                    ,
             int                nDOFX                       ,
             Array4<Real const> CoarseToFineProjectionMatrix,
             RunOn              runon                        ) override;

    /**
    * \brief Coarse to fine interpolation in space.
    *
    * \param CrseFab
    * \param FineFab
    * \param nComp
    * \param fine_region
    * \param RefRatio
    * \param nDOFX
    * \param CoarseToFineProjectionMatrix
    * \param runon
    */
    void interpPointWise
           ( const FArrayBox &  CrseFab                     ,
             FArrayBox       &  FineFab                     ,
             int                nComp                       ,
             const Box       &  fine_region                 ,
             const IntVect   &  RefRatio                    ,
             int                nDOFX                       ,
             Array4<Real const> CoarseToFineProjectionMatrix,
             RunOn              runon                        ) override;
};

//! CONSTRUCT A GLOBAL OBJECT OF EACH VERSION.
extern AMREX_EXPORT DGInterp dg_interp;
extern AMREX_EXPORT CGInterp cg_interp;

}

#endif /*_INTERPOLATER_H_*/
