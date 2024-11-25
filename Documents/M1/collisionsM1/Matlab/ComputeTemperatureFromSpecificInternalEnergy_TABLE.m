function T = ComputeTemperatureFromSpecificInternalEnergy_TABLE(D, E, Y)

    global g_D1D g_T1D g_Y1D g_E g_EOS_OS; 
    
    offset = g_EOS_OS(g_E.i);
    
    T = zeros(size(D));
    
    for i = 1 : length(D)
        error = 0;
        iD = find( g_D1D < D(i), 1, 'last' );
        iY = find( g_Y1D < Y(i), 1, 'last' );
%         if ( size(iD,1)*size(iY,1) ~= 1 )
%             disp( 'Error in interpolateEos: Outside table boundary.' )
        if ( size( iD, 1 ) == 0 )
            disp( 'Error in interpolateEos: Density out of left boundary.' )
            disp('Error in interpolateEos: Replaced by the boundary.')
            iD = 1;
        end
        if ( size( iY, 1 ) == 0 )
            disp( 'Error in interpolateEos: Electron fraction out of left boundary.' )
            disp('Error in interpolateEos: Replaced by the boundary.')
            iY = 1;
        end
        if ( iD == length(g_D1D) )
            disp( 'Error in interpolateEos: Density out of right boundary.' )
            disp('Error in interpolateEos: Replaced by the boundary.')
            iD = iD - 1;
        end
        if ( iY == length(g_Y1D) )
            disp( 'Error in interpolateEos: Electron fraction out of right boundary.' )
            disp('Error in interpolateEos: Replaced by the boundary.')
            iY = iY - 1;
        end        
%         break
%         end
        %logE = log10(E(i) + offset);
        tmpE = E(i);
        
        iE1 = find( g_E.v3D(iD,  :,iY  ) < tmpE, 1, 'last' );
        iE2 = find( g_E.v3D(iD+1,:,iY  ) < tmpE, 1, 'last' );        
        iE3 = find( g_E.v3D(iD  ,:,iY+1) < tmpE, 1, 'last' );
        iE4 = find( g_E.v3D(iD+1,:,iY+1) < tmpE, 1, 'last' );        
        
        iEa = max( min( [iE1,iE2,iE3,iE4] ) - 1, 1 ); % left
        iEb = min( max( [iE1,iE2,iE3,iE4] ) + 2, length(g_T1D) ); % right
        nPtsE = length( g_T1D(iEa:iEb) );
        ptsD = D(i)*ones(nPtsE,1);
        ptsT = g_T1D(iEa:iEb);
        ptsY = Y(i)*ones(nPtsE,1);
%         ptsE = zeros(nPtsE,1);
        
        % now interpolate ptsE
        
        ptsE = interpolateEos( ptsD, ptsT, ptsY, g_D1D(iD:iD+1), g_T1D(iEa:iEb),...
            g_Y1D(iY:iY+1), g_E.v3D(iD:iD+1, iEa:iEb, iY:iY+1), offset);
        
        
        
        % if E is not in the interpolant
        % Question: is E monotone w.r.t. T?
        if ((ptsE(1)-tmpE)*(ptsE(nPtsE)-tmpE) > 0.0)
            if(E(i) < ptsE(1))
                tmpE = 0.5 * ( ptsE(1) + ptsE(2) );
            else
                tmpE = 0.5 * ( ptsE(nPtsE-1) + ptsE(nPtsE) );
            end
            error = 1;
        end
        
        iE = find( ptsE<tmpE, 1, 'last' );

        T(i) = 10^(log10(ptsT(iE))...
                + log10( ptsT(iE+1)/ptsT(iE) ) ...
                    * log10( (tmpE+offset)/(ptsE(iE)+offset) ) ...
                        / log10( (ptsE(iE+1)+offset)/(ptsE(iE)+offset) ) );
        if (error == 1)
            warning('E(i) not in the interpolant');
        end
        
    end
%     
%         DO i = 1, SIZE( D )
% 
%       Error = 0
% 
%       ilD = Index1D( D(i), D_table, SIZE( D_table ) )
%       ilY = Index1D( Y(i), Y_table, SIZE( Y_table ) )
% 
%       logE = LOG10( E(i) + Offset )
%       ilE1 = Index1D( logE, E_table(ilD,  :,ilY  ), SIZE( T_Table ) )
%       ilE2 = Index1D( logE, E_table(ilD+1,:,ilY  ), SIZE( T_Table ) )
%       ilE3 = Index1D( logE, E_table(ilD,  :,ilY+1), SIZE( T_Table ) )
%       ilE4 = Index1D( logE, E_table(ilD+1,:,ilY+1), SIZE( T_Table ) )
% 
%       ilEa = MAX( MINVAL( [ilE1,ilE2,ilE3,ilE4] ) - 1, 1 )
%       ilEb = MIN( MAXVAL( [ilE1,ilE2,ilE3,ilE4] ) + 2, SIZE( T_Table ) )
% 
%       nPtsE = SIZE( T_Table(ilEa:ilEb) )
%       ALLOCATE( ptsD(nPtsE), ptsT(nPtsE), ptsY(nPtsE), ptsE(nPtsE) )
%       ptsD = D(i)
%       ptsT = T_Table(ilEa:ilEb)
%       ptsY = Y(i)
% 
%       CALL LogInterpolateSingleVariable &
%              ( ptsD, ptsT, ptsY, D_Table(ilD:ilD+1), T_Table(ilEa:ilEb), &
%                Y_Table(ilY:ilY+1), LogInterp, Offset, &
%                E_Table(ilD:ilD+1,ilEa:ilEb,ilY:ilY+1), ptsE )
% 
%       tmpE = E(i)
% 
%       IF( (ptsE(1)-E(i))*(ptsE(nPtsE)-E(i)) > 0.0_DP )THEN
% 
%         IF( Debug )THEN
%           WRITE(*,*)
%           WRITE(*,'(A4,A)') &
%             '', 'Warning: ComputeTempFromIntEnergy_Lookup'
%           WRITE(*,'(A6,A20,ES10.4E2,A9,ES10.4E2)') &
%             '', 'No Root Between T = ', ptsT(1), ' and T = ', ptsT(nPtsE)
%           WRITE(*,*)
%           WRITE(*,*) '  nPtsE = ', nPtsE
%           WRITE(*,*) '  ptsE  = ', ptsE
%           WRITE(*,*) '    ia  = ', ilEa
%           WRITE(*,*) '    ib  = ', ilEb
%           WRITE(*,*) '    Ta  = ', T_Table(ilEa)
%           WRITE(*,*) '    Tb  = ', T_Table(ilEb)
%           WRITE(*,*) '    Ea  = ', ptsE(1)
%           WRITE(*,*) '    Eb  = ', ptsE(nPtsE)
%           WRITE(*,*) '     i  = ', i
%           WRITE(*,*) '     E  = ', E(i)
%           WRITE(*,*) '     D  = ', D(i)
%           WRITE(*,*) '     Y  = ', Y(i)
%           WRITE(*,*)
%           STOP
%         END IF
% 
%         ! --- Reset Energy Density ---
% 
%         IF( E(i) < ptsE(1) .AND. E(i) < ptsE(nPtsE) )THEN
% 
%           tmpE = 0.5 * ( ptsE(1) + ptsE(2) )
% 
%         ELSE
% 
%           tmpE = 0.5 * ( ptsE(nPtsE-1) + ptsE(nPtsE) )
% 
%         END IF
% 
%         Error = 1
% 
%       END IF
% 
%       ilE = Index1D( tmpE, ptsE, nPtsE )
% 
%       T(i) &
%         = 10.d0**( LOG10( ptsT(ilE) ) &
%                    + LOG10( ptsT(ilE+1)/ptsT(ilE) ) &
%                      * LOG10( (tmpE+Offset)/(ptsE(ilE)+Offset) ) &
%                        / LOG10( (ptsE(ilE+1)+Offset)/(ptsE(ilE)+Offset) ) )
% 
%       IF( Error == 1 .AND. Debug )THEN
%         WRITE(*,*)
%         WRITE(*,*) '  T = ', T(i)
%         WRITE(*,*)
%       END IF
% 
%       DEALLOCATE( ptsD, ptsT, ptsY, ptsE )
% 
%     END DO
% 
%   END SUBROUTINE ComputeTempFromIntEnergy_Lookup

end