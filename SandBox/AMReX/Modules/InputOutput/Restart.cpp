/*
Restart.cpp

Originally created by nataraj for one MultiFab

Modified by sjdunham for multiple MultiFabs

*/

#include <cstddef> /* For NULL */
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
  void writecheckpointfile
         ( int StepNo[], int nLevels,
           Real dt[], Real time[],
           Real BaryonicMassArr   [],
           Real EulerMomentumX1Arr[],
           Real EulerMomentumX2Arr[],
           Real EulerMomentumX3Arr[],
           Real EulerEnergyArr    [],
           Real ElectronNumberArr [],
           Real ADMMassArr        [],
           BoxArray** pBA,
           int RwChkFields_uGF = 1,
           int RwChkFields_uCF = 1,
           int RwChkFields_uDF = 0,
           int RwChkFields_uCR = 0,
           MultiFab** pMF_uGF = NULL,
           MultiFab** pMF_uCF = NULL,
           MultiFab** pMF_uDF = NULL,
           MultiFab** pMF_uCR = NULL )
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

    bool WriteFields_uGF = true;
    if( RwChkFields_uGF == 0 ) WriteFields_uGF = false;

    bool WriteFields_uCF = true;
    if( RwChkFields_uCF == 0 ) WriteFields_uCF = false;

    bool WriteFields_uDF = false;
    if( RwChkFields_uDF == 1 ) WriteFields_uDF = true;

    bool WriteFields_uCR = false;
    if( RwChkFields_uCR == 1 ) WriteFields_uCR = true;

    ParmParse pp("thornado");
    chk_file = "chk";
    pp.query("CheckpointFileNameRoot",chk_file);

    const std::string& checkpointname
                         = amrex::Concatenate( chk_file, StepNo[0], 8 );

    bool Verbose = true;
    if ( Verbose && ParallelDescriptor::IOProcessor() )
    {
        amrex::Print() << "\n    Writing CheckpointFile "
                     << checkpointname << "\n\n";
    }

    const int FinestLevel = nLevels-1;

    bool callBarrier = true;

    // ---- Prebuild a hierarchy of directories
    // ---- dirName is built first. If dirName exists, it is renamed. Then build
    // ----   dirName/subDirPrefix_0 .. dirName/subDirPrefix_nLevels-1
    // ---- If callBarrier is true, call ParallelDescriptor::Barrier()
    // ----   after all directories are built
    // ---- ParallelDescriptor::IOProcessor() creates the directories
    amrex::PreBuildDirectorHierarchy
             ( checkpointname, "Level_", nLevels, callBarrier );

    // Write Header file
    if ( ParallelDescriptor::IOProcessor() )
    {

      VisMF::IO_Buffer io_buffer( VisMF::IO_Buffer_Size );
      std::ofstream HeaderFile;
      HeaderFile.rdbuf() -> pubsetbuf( io_buffer.dataPtr(), io_buffer.size() );
      std::string HeaderFileName( checkpointname + "/Header" );
      HeaderFile.open( HeaderFileName.c_str(), std::ofstream::out   |
                                               std::ofstream::trunc |
                                               std::ofstream::binary );

      if( ! HeaderFile.good() )
      {
        amrex::FileOpenFailed( HeaderFileName );
      }

      HeaderFile.precision(17);

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

      // Write out initial values for tally
      HeaderFile << BaryonicMassArr   [0] << "\n";
      HeaderFile << BaryonicMassArr   [1] << "\n";
      HeaderFile << EulerMomentumX1Arr[0] << "\n";
      HeaderFile << EulerMomentumX1Arr[1] << "\n";
      HeaderFile << EulerMomentumX2Arr[0] << "\n";
      HeaderFile << EulerMomentumX2Arr[1] << "\n";
      HeaderFile << EulerMomentumX3Arr[0] << "\n";
      HeaderFile << EulerMomentumX3Arr[1] << "\n";
      HeaderFile << EulerEnergyArr    [0] << "\n";
      HeaderFile << EulerEnergyArr    [1] << "\n";
      HeaderFile << ElectronNumberArr [0] << "\n";
      HeaderFile << ElectronNumberArr [1] << "\n";
      HeaderFile << ADMMassArr        [0] << "\n";
      HeaderFile << ADMMassArr        [1] << "\n";
      HeaderFile << ADMMassArr        [2] << "\n"; // ADMMass_Interior

      // Write the BoxArray at each level
      for( int iLevel = 0; iLevel <= FinestLevel; ++iLevel )
      {
        BoxArray& BA = *pBA[iLevel];
        BA.writeOn( HeaderFile );
        HeaderFile << "\n";
      }

    } // End of Header file writing

//    VisMF::SetNOutFiles( 1 ); // write MultiFab data in serial

    // Write the MultiFab data to, e.g., chk00010/Level_0/
    for( int iLevel = 0; iLevel <= FinestLevel; ++iLevel )
    {

      if( WriteFields_uGF )
      {
        MultiFab& MF_uGF = *pMF_uGF[iLevel];
        VisMF::Write( MF_uGF,
                      amrex::MultiFabFileFullPrefix
                        ( iLevel, checkpointname,
                          "Level_", "Geometry" ) );
      }

      if( WriteFields_uCF )
      {
        MultiFab& MF_uCF = *pMF_uCF[iLevel];
        VisMF::Write( MF_uCF,
                      amrex::MultiFabFileFullPrefix
                        ( iLevel, checkpointname,
                          "Level_", "Conserved_Euler" ) );
      }

      if( WriteFields_uDF )
      {
        MultiFab& MF_uDF = *pMF_uDF[iLevel];
        VisMF::Write( MF_uDF,
                      amrex::MultiFabFileFullPrefix
                        ( iLevel, checkpointname,
                          "Level_", "Diagnostic_Euler" ) );
      }

      if( WriteFields_uCR )
      {
        MultiFab& MF_uCR = *pMF_uCR[iLevel];
        VisMF::Write( MF_uCR,
                      amrex::MultiFabFileFullPrefix
                               ( iLevel, checkpointname,
                                 "Level_", "Conserved_TwoMoment" ) );
      }
    }

  } // End of writecheckpointfile

  void readheaderandboxarraydata
         ( int FinestLevelArr[], int StepNo[],
	   Real dt[], Real Time[],
           Real BaryonicMassArr   [],
           Real EulerMomentumX1Arr[],
           Real EulerMomentumX2Arr[],
           Real EulerMomentumX3Arr[],
           Real EulerEnergyArr    [],
           Real ElectronNumberArr [],
           Real ADMMassArr        [],
           BoxArray** pba, DistributionMapping** pdm, int iChkFile )
  {

    int FinestLevel;
    Real BaryonicMass_Initial;
    Real BaryonicMass_OffGrid;
    Real EulerMomentumX1_Initial;
    Real EulerMomentumX1_OffGrid;
    Real EulerMomentumX2_Initial;
    Real EulerMomentumX2_OffGrid;
    Real EulerMomentumX3_Initial;
    Real EulerMomentumX3_OffGrid;
    Real EulerEnergy_Initial;
    Real EulerEnergy_OffGrid;
    Real ElectronNumber_Initial;
    Real ElectronNumber_OffGrid;
    Real ADMMass_Initial;
    Real ADMMass_OffGrid;
    Real ADMMass_Interior;

    ParmParse pp("thornado");
    chk_file = "chk";
    pp.query("CheckpointFileNameRoot",chk_file);

    std::stringstream sChkFile;
    sChkFile << chk_file << std::setw(8) << std::setfill('0') << iChkFile;
    restart_chkfile = sChkFile.str();

    amrex::Print() << "\n      Restart from checkpoint " << restart_chkfile
                   << "\n";

    // Header
    std::string File( restart_chkfile + "/Header" );

    VisMF::IO_Buffer io_buffer( VisMF::GetIOBufferSize() );

    Vector<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile( File, fileCharPtr );
    std::string fileCharPtrString( fileCharPtr.dataPtr() );
    std::istringstream is( fileCharPtrString, std::istringstream::in );

    std::string line, word;

    // Read in title line
    std::getline( is, line );

    // Read in FinestLevel
    is >> FinestLevel;
    GotoNextLine( is );
    FinestLevelArr[0] = FinestLevel;

    // Read in array of StepNo
    std::getline( is, line );
    {
        std::istringstream lis( line );
        int i = 0;
        while( lis >> word )
	{
          StepNo[i++] = std::stoi(word);
        }
    }

    // Read in array of dt
    std::getline( is, line );
    {
      std::istringstream lis( line );
      int i = 0;
      while( lis >> word )
      {
        dt[i++] = std::stod( word );
      }
    }

    // Read in array of time
    std::getline( is, line );
    {
      std::istringstream lis( line );
      int i = 0;
      while( lis >> word )
      {
        Time[i++] = std::stod( word );
      }
    }

    // Read in initial values for tally

    is >> BaryonicMass_Initial;
    GotoNextLine( is );
    BaryonicMassArr[0] = BaryonicMass_Initial;
    is >> BaryonicMass_OffGrid;
    GotoNextLine( is );
    BaryonicMassArr[1] = BaryonicMass_OffGrid;

    is >> EulerMomentumX1_Initial;
    GotoNextLine( is );
    EulerMomentumX1Arr[0] = EulerMomentumX1_Initial;
    is >> EulerMomentumX1_OffGrid;
    GotoNextLine( is );
    EulerMomentumX1Arr[1] = EulerMomentumX1_OffGrid;

    is >> EulerMomentumX2_Initial;
    GotoNextLine( is );
    EulerMomentumX2Arr[0] = EulerMomentumX2_Initial;
    is >> EulerMomentumX2_OffGrid;
    GotoNextLine( is );
    EulerMomentumX2Arr[1] = EulerMomentumX2_OffGrid;

    is >> EulerMomentumX3_Initial;
    GotoNextLine( is );
    EulerMomentumX3Arr[0] = EulerMomentumX3_Initial;
    is >> EulerMomentumX3_OffGrid;
    GotoNextLine( is );
    EulerMomentumX3Arr[1] = EulerMomentumX3_OffGrid;

    is >> EulerEnergy_Initial;
    GotoNextLine( is );
    EulerEnergyArr[0] = EulerEnergy_Initial;
    is >> EulerEnergy_OffGrid;
    GotoNextLine( is );
    EulerEnergyArr[1] = EulerEnergy_OffGrid;

    is >> ElectronNumber_Initial;
    GotoNextLine( is );
    ElectronNumberArr[0] = ElectronNumber_Initial;
    is >> ElectronNumber_OffGrid;
    GotoNextLine( is );
    ElectronNumberArr[1] = ElectronNumber_OffGrid;

    is >> ADMMass_Initial;
    GotoNextLine( is );
    ADMMassArr[0] = ADMMass_Initial;
    is >> ADMMass_OffGrid;
    GotoNextLine( is );
    ADMMassArr[1] = ADMMass_OffGrid;
    is >> ADMMass_Interior;
    GotoNextLine( is );
    ADMMassArr[2] = ADMMass_Interior;

    // Read in level 'iLevel' BoxArray from Header
    for( int iLevel = 0; iLevel <= FinestLevel; ++iLevel )
    {

      BoxArray& ba = *pba[iLevel];
      ba = BoxArray();
      ba.readFrom( is );
      pba[iLevel] = &ba;
      GotoNextLine( is );

      // Create a distribution mapping
      DistributionMapping dm{ ba, ParallelDescriptor::NProcs() };
      *pdm[iLevel] = dm;

    }

  } // End of readheaderandboxarraydata function

  void readmultifabdata
    ( int FinestLevel, MultiFab** MF, int iMF, int iChkFile )
  {

    ParmParse pp("thornado");
    chk_file = "chk";
    pp.query("CheckpointFileNameRoot",chk_file);

    std::stringstream sChkFile;

    sChkFile << chk_file << std::setw(8) << std::setfill('0') << iChkFile;
    restart_chkfile = sChkFile.str();

    // Header
    std::string File( restart_chkfile + "/Header" );

    VisMF::IO_Buffer io_buffer( VisMF::GetIOBufferSize() );

    Vector<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile( File, fileCharPtr );
    std::string fileCharPtrString( fileCharPtr.dataPtr() );
    std::istringstream is( fileCharPtrString, std::istringstream::in );

    // Read in the MultiFab data
    std::string MF_Name;
    switch( iMF )
    {
      case 0:
        MF_Name = "Geometry";
	break;
      case 1:
	MF_Name = "Conserved_Euler";
	break;
      case 2:
	MF_Name = "Diagnostic_Euler";
	break;
      case 3:
	MF_Name = "Conserved_TwoMoment";
	break;
      default:
        std::cout << "Invalid." << std::endl;
    }

    for( int iLevel = 0; iLevel <= FinestLevel; ++iLevel )
    {
      VisMF::Read( *MF[iLevel],
        amrex::MultiFabFileFullPrefix
                 ( iLevel, restart_chkfile, "Level_", MF_Name ) );
    }

  } // End of readmultifabdata function

  void GotoNextLine ( std::istream& is )
  {
    constexpr std::streamsize bl_ignore_max { 100000 };
    is.ignore( bl_ignore_max, '\n' );
  }

} // End of extern "C" block
