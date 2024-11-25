function [ eGrid ] = CreateEnergyGrid( E_lo, E_hi, nE, nGhostE, dE_0, nUni )

  if ~exist( 'dE_0', 'var' )
    dE_0 = (E_hi-E_lo)/nE;
  end
  
  if ~exist( 'nUni', 'var' )
    nUni = 1;
  end
  
  eGrid.nE     = nE;
  eGrid.E_lo   = E_lo;
  eGrid.E_hi   = E_hi;
  eGrid.nGhost = nGhostE;
  
  eGrid.Center  = zeros( nE+2*nGhostE, 1 );
  eGrid.Width   = zeros( nE+2*nGhostE, 1 );
  eGrid.Edge_lo = zeros( nE+2*nGhostE, 1 );
  eGrid.Edge_hi = zeros( nE+2*nGhostE, 1 );
  eGrid.V       = zeros( nE+2*nGhostE, 1 );
  
  if( dE_0 >= (E_hi-E_lo)/nE )
    
    % --- Return Unigrid ---
    
    eGrid.Width(:) = ( E_hi - E_lo ) / nE;
  
    for i = 1 : nE + 2 * nGhostE

      j = i - nGhostE;
    
      eGrid.Edge_lo(i) = (j-1) * eGrid.Width(i);
      eGrid.Edge_hi(i) = (j-0) * eGrid.Width(i);
    
      eGrid.Center(i) = 0.5 * ( eGrid.Edge_lo(i) + eGrid.Edge_hi(i) );
      
    end
    
  else
      
    % --- Return Geometric Grid ---
    
    nUni = max( 1, nUni );
    
    FUN = @( x ) (x.^(nE-(nUni-1))-1.0).*dE_0-(x-1.0).*(E_hi-(E_lo+(nUni-1)*dE_0));
    
    delta = 1.0d-8;
    alpha = fzero( FUN, [(1.0+delta) 2.0] ); % --- Zoom Factor
    
    iOS = nGhostE;
    
    eGrid.Width  (iOS+1:iOS+nUni) = dE_0;
    eGrid.Edge_lo(iOS+1:iOS+nUni) = E_lo + ((1:nUni)-1) .* dE_0;
    eGrid.Edge_hi(iOS+1:iOS+nUni) = E_lo + ((1:nUni)-0) .* dE_0;
    eGrid.Center (iOS+1:iOS+nUni)...
      = 0.5 * ( eGrid.Edge_hi(iOS+1:iOS+nUni) + eGrid.Edge_lo(iOS+1:iOS+nUni) );
    
    for i = iOS + nUni + 1 : iOS + nE
      
      eGrid.Width  (i) = alpha * eGrid.Width(i-1);
      eGrid.Edge_lo(i) = eGrid.Edge_hi(i-1);
      eGrid.Edge_hi(i) = eGrid.Edge_lo(i) + eGrid.Width(i);
      eGrid.Center (i) = 0.5 * ( eGrid.Edge_hi(i) + eGrid.Edge_lo(i) );
      
    end
    
  end
  
  for i = 1 + nGhostE : nE + nGhostE
    eGrid.V(i) = (eGrid.Edge_hi(i)^3-eGrid.Edge_lo(i)^3)/3.0;
  end
  
end