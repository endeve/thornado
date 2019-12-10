/*
Restart.cpp

Originally created by nataraj for one MultiFab

Modified by sjdunham for multiple MultiFabs (not working yet)

*/

#include <sstream>
#include <iomanip>

#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>

#ifdef BL_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

#include "Restart.H"

using namespace amrex;

extern "C"
{
  void writefieldsamrex_checkpoint
         ( int StepNo[], int nLevels,
           Real dt[], Real time[], Real t_wrt[],
           MultiFab** MF_uCR, MultiFab** MF_uPR )
  {

    // chk00010            Write a checkpoint file with this root directory
    // chk00010/Header     This contains information you need to save
    //                       (e.g., FinestLevel, time, etc.), and also
    //                       the BoxArrays at each level
    // chk00010/Level_0/
    // chk00010/Level_1/
    // etc.                These subdirectories will hold the MultiFab
    //                     data at each level of refinement

    // Checkpoint file name, e.g., chk00010


    const std::string& checkpointname
                         = amrex::Concatenate( chk_file, StepNo[0], 8 );

    if ( ParallelDescriptor::IOProcessor() )
      amrex::Print() << "\n    Writing checkpoint " << checkpointname << "\n";

    const int FinestLevel = nLevels-1;

    // ---- Prebuild a hierarchy of directories
    // ---- dirName is built first. If dirName exists, it is renamed. Then build
    // ----   dirName/subDirPrefix_0 .. dirName/subDirPrefix_nLevels-1
    // ---- If callBarrier is true, call ParallelDescriptor::Barrier()
    // ----   after all directories are built
    // ---- ParallelDescriptor::IOProcessor() creates the directories
    amrex::PreBuildDirectorHierarchy( checkpointname, "Level_", nLevels, true );

    // Write Header file
    if ( ParallelDescriptor::IOProcessor() )
    {

      std::string HeaderFileName( checkpointname + "/Header" );
      std::ofstream HeaderFile( HeaderFileName.c_str(), std::ofstream::out   |
				                        std::ofstream::trunc |
				                        std::ofstream::binary );
      if( ! HeaderFile.good() )
      {
        amrex::FileOpenFailed( HeaderFileName );
      }

      HeaderFile.precision(17);

      VisMF::IO_Buffer io_buffer( VisMF::IO_Buffer_Size );
      HeaderFile.rdbuf() -> pubsetbuf( io_buffer.dataPtr(), io_buffer.size() );

      // Write out title line
      HeaderFile << "Checkpoint file\n";

      // Write out FinestLevel
      HeaderFile << FinestLevel << "\n";

      // Write out array of StepNo
      for( int i = 0; i < nLevels; ++i )
      {
        HeaderFile << StepNo[i] << " ";
      }
      HeaderFile << "\n";

      // Write out array of dt
      for(int i = 0; i < nLevels; ++i)
      {
        HeaderFile << dt[i] << " ";
      }
      HeaderFile << "\n";

      // Write out array of time
      for( int i = 0; i < nLevels; ++i )
      {
        HeaderFile << time[i] << " ";
      }
      HeaderFile << "\n";

      // Write out t_wrt
      HeaderFile << t_wrt[0] << "\n";

      // Write the BoxArray at each level

    } // End of Header file writing

    // Write the MultiFab data to, e.g., chk00010/Level_0/
    for( int iLevel = 0; iLevel <= FinestLevel; ++iLevel )
    {
      MultiFab& MF_uCR1 = *MF_uCR[iLevel];
      MultiFab& MF_uPR1 = *MF_uPR[iLevel];
      VisMF::Write( MF_uCR1, amrex::MultiFabFileFullPrefix
                                      ( iLevel, checkpointname,
                                        "Level_", "Conserved" ) );
      VisMF::Write( MF_uPR1, amrex::MultiFabFileFullPrefix
                                      ( iLevel, checkpointname,
                                        "Level_", "Primitive" ) );
    }

  } // End of WriteCheckpointFile function
} // End of extern "C" block
