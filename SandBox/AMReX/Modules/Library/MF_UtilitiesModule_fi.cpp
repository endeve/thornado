#include <iostream>

#include <AMReX_BoxArray.H>

using namespace amrex;

extern "C"
{
    void print_boxarray
      ( const BoxArray& ba )
    {
        std::cout << ba;
    }

}
