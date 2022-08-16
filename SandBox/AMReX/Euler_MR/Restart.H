#include <string>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

extern "C"
{
  // Utility to skip to next line in Header
  static void GotoNextLine ( std::istream& is );
  std::string restart_chkfile = "";
  std::string chk_file = "";
}
