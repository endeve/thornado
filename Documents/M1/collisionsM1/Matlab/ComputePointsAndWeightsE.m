function [E_N, W2, W3] = ComputePointsAndWeightsE(NumofEBins, NorderE)
% generate E grid on [0 E_max]

%     global PlanckConstant SpeedOfLight AtomicMassUnit;
%     hc = PlanckConstant * SpeedOfLight;

    
    
    E_min = 0; E_max = 300; 
    
    E_N = zeros(NumofEBins*NorderE, 1);
    W2 = zeros(NumofEBins*NorderE, 1);
    W3 = zeros(NumofEBins*NorderE, 1);
    eGrid = CreateEnergyGrid( E_min, E_max, NumofEBins, 0, (E_max-E_min)/NumofEBins/10);
    
%     E_edge = linspace(0, E_max, NumofEBins+1); %equidistance
%     distE = E_edge(2:end) - E_edge(1:end-1);
    
    
    for iE = 1 : NumofEBins
        
        [E_q, W_q]=lgwt(NorderE, eGrid.Edge_lo(iE), eGrid.Edge_hi(iE));
        
%         [E_q, W_q]=lgwt(NorderE, E_edge(iE), E_edge(iE+1));
        
        [E_q, idx] = sort(E_q,'ascend');
        W_q = W_q(idx);
        iStart = (iE-1)*length(E_q) + 1;
        iEnd = iStart + length(E_q) - 1;
        E_N(iStart:iEnd) = E_q;
        
%         W2(iStart:iEnd) = W_q .* eGrid.Width(iE) .* ( E_q ).^2;
        W2(iStart:iEnd) = W_q .* ( E_q ).^2;
        W3(iStart:iEnd) = W2(iStart:iEnd) .* ( E_q );
        
        % scaled version
%         W2(iStart:iEnd) = W_q .* ( distE(iE) / hc ) .* ( E_q / hc ).^2;
%         W3(iStart:iEnd) = W2(iStart:iEnd) .* ( E_q / AtomicMassUnit );
    end

end
