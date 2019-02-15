close all
clear all

MeVperK = 8.62d-11;

[D,T,Y,nD,nT,nY,P,S,E,Me,Mp,Mn,Xp,Xn,Xa,Xh,Zh,Ah,Eh,Eth,Gm,OS]...
  = readEosTable( 'wl-EOS-SFHo-15-25-50.h5' );

interpolateEos( 1.d12, 1.604d10, 0.4, D, T, Y, Xp, 0.0 )
interpolateEos( 1.d12, 1.604d10, 0.4, D, T, Y, Xn, 0.0 )
interpolateEos( 1.d12, 1.604d10, 0.4, D, T, Y, Xa, 0.0 )
interpolateEos( 1.d12, 1.604d10, 0.4, D, T, Y, Xh, 0.0 )

% minEta_e = 1.0d100;
% maxEta_e = 0.0d000;
% eta = zeros( nD*nT*nY, 1 ); i = 1;
% for iY = 1 : nY
%   for iT = 1 : nT
%     for iD = 1 : nD
%       minEta_e = min( [ Me(iD,iT,iY)/(T(iT)*MeVperK) minEta_e ] );
%       maxEta_e = max( [ Me(iD,iT,iY)/(T(iT)*MeVperK) maxEta_e ] );
%       eta(i) = Me(iD,iT,iY)/(T(iT)*MeVperK); i = i + 1;
%     end
%   end
% end
% eta = sort( eta );