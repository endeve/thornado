#include <AMReX_DistributionMapping.H>
#include <AMReX_Print.H>

using namespace amrex;

extern "C" {
  bool amrex_fi_distromap_issame
    ( const DistributionMapping& dma, const DistributionMapping& dmb )
  {
    return dma == dmb;
  }
}
