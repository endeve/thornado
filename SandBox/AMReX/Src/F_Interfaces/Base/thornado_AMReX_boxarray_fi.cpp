#include <AMReX_BoxArray.H>
#include <AMReX_Print.H>
#include <AMReX_Utility.H>

using namespace amrex;

extern "C" {
  bool amrex_fi_boxarray_issame( const BoxArray& baa, const BoxArray& bab )
  {
    return baa == bab;
  }
}
