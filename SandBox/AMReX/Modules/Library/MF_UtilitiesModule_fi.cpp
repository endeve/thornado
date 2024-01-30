#include <iostream>

#include <AMReX_BoxArray.H>

using namespace amrex;

extern "C"
{
    void print_boxarray
      ( const BoxArray& ba )
    {
        std::cout << ba << "\n";
        std::cout << "Total number of boxes    = " << ba.size() << "\n";
        std::cout << "Total number of elements = " << ba.numPts() << "\n";

    }

}
