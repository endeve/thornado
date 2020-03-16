function...
    [Lambda_n, I_n, Lambda_a, I_a, R_1_cond] = AnalyzeEigensystem(D, T, Y, P, E, cs_T, alpha, beta, dPdDe, Gm_A, Gm_T, R_1_n, invR_1_n, dFdU_1_n, R_1_a, invR_1_a, dFdU_1_a)

    nPoints = size(D, 1)

    for iPoint = 1:nPoints

      Lambda_n(:,:,iPoint) = invR_1_n(:,:,iPoint) * dFdU_1_n(:,:,iPoint) * R_1_n(:,:,iPoint);
      I_n(:,:,iPoint) = invR_1_n(:,:,iPoint) * R_1_n(:,:,iPoint);

      Lambda_a(:,:,iPoint) = invR_1_a(:,:,iPoint) * dFdU_1_a(:,:,iPoint) * R_1_a(:,:,iPoint);
      I_a(:,:,iPoint) = invR_1_a(:,:,iPoint) * R_1_a(:,:,iPoint);

      R_1_cond(iPoint) = cond(R_1_n(:,:,iPoint));

   end

end
