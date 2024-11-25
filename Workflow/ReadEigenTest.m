

function...
    [D, T, Y, P, E, cs_T, alpha, beta, dPdDe, Gm_A, Gm_T, R_1_n, invR_1_n, dFdU_1_n, R_1_a, invR_1_a, dFdU_1_a] = ReadEigenTest( Directory )

if( exist( 'File', 'var' ) )
    DirName = File;
else
     DirName = './Output';
end

if( exist( 'Directory', 'var' ) )
    DirName = Directory;
else
    DirName = '../SandBox/dgExperiments_Euler_NonRelativistic_TABLE/Executables/';
end

 D = ReadVector([DirName 'D.dat']);
 T = ReadVector([DirName 'T.dat']);
 Y = ReadVector([DirName 'Y.dat']);
 P = ReadVector([DirName 'P.dat']);
 E = ReadVector([DirName 'E.dat']);
 cs_T = ReadVector([DirName 'cs_T.dat']);
 alpha = ReadVector([DirName 'alpha.dat']);
 beta = ReadVector([DirName 'beta.dat']);
 dPdDe = ReadVector([DirName 'dPdDe.dat']);
 Gm_A = ReadVector([DirName 'Gm_A.dat']);
 Gm_T = ReadVector([DirName 'Gm_T.dat']);
 R_1_n = ReadRank3Tensor([DirName 'R_1_n.dat']);
 invR_1_n = ReadRank3Tensor([DirName 'invR_1_n.dat']);
 dFdU_1_n = ReadRank3Tensor([DirName 'dFdU_1_n.dat']);
 R_1_a = ReadRank3Tensor([DirName 'R_1_a.dat']);
 invR_1_a = ReadRank3Tensor([DirName 'invR_1_a.dat']);
 dFdU_1_a = ReadRank3Tensor([DirName 'dFdU_1_a.dat']);


end
