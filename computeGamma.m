%% computeGamma
% FastNLS simplified implementation for MATLAB coder - by Søren Vøgg Lyster

function [Gamma1, Gamma2]= computeGamma(L, F, pitchBounds, ...
    crossCorrelationVectors, nPitches, validFftIndices)
    a1 = crossCorrelationVectors(2:L, :)+...
        [zeros(1,nPitches);...
        crossCorrelationVectors((1:L-2)+2, :)];
    a2 = crossCorrelationVectors(2:L, :)-...
        [zeros(1,nPitches);...
        crossCorrelationVectors((1:L-2)+2, :)];
    R1 = nan(2,nPitches);
    R2 = nan(2,nPitches);
    phi1 = nan(L,nPitches);
    psi1 = nan(L,nPitches);
    phi2 = nan(L,nPitches);
    psi2 = nan(L,nPitches);
    alpha1 = nan(1,nPitches);
    alpha2 = nan(1,nPitches);
    gammaNew1 = nan(L,nPitches);
    gammaNew2 = nan(L,nPitches);
    gammaOld1 = nan(L,nPitches);
    gammaOld2 = nan(L,nPitches);
    Gamma1 = cell(L,1);
    Gamma2 = cell(L,1);
    for i = 1:L
        Gamma1{i} = nan;
        Gamma2{i} = nan;
    end
    coder.varsize('Gamma1','Gamma2', [30, 5000]);
    coder.varsize('a1', 'a2', 'R1', 'alpha1',...
        'gammaOld1', 'gammaNew1', 'phi1', 'psi1', 'R2',...
        'alpha2', 'gammaOld2', 'gammaNew2', 'phi2', 'psi2');
    for l = 1:L
        maxFftIndex = floor(min(F*pitchBounds(2),F/(2*l)-1));
        nPitches = nPitches-(validFftIndices(end)-maxFftIndex);
        validPitchIndices = 1:nPitches;
        validFftIndices = validFftIndices(validPitchIndices);
        if l == 1
           [psiout, phiout, gammaOld1] = computeGammaSingleSinus(...
               crossCorrelationVectors(1:3,validPitchIndices),...
               a1(1,validPitchIndices), true);
           psi1(1,:) = psiout;
           phi1(1,:) = phiout;
           Gamma1{l} = gammaOld1(1,:);
           [psiout, phiout, gammaOld2] = computeGammaSingleSinus(...
               crossCorrelationVectors(1:3,validPitchIndices),...
               a2(1,validPitchIndices), false);
           psi2(1,:) = psiout;
           phi2(1,:) = phiout;
           Gamma2{l} = gammaOld2(1,:);
        elseif l == 2
            [R1, alpha1, gammaNewout] = computeGammaTwoSinus(...
                crossCorrelationVectors(1:5,validPitchIndices),...
                psi1(1,validPitchIndices),...
                gammaOld1(validPitchIndices), true);
            gammaNew1(1:2,:) = gammaNewout;
            Gamma1{l} = gammaNew1(1:2,:);
            [R2, alpha2, gammaNewout] = computeGammaTwoSinus(...
                crossCorrelationVectors(1:5,validPitchIndices),...
                psi2(1,validPitchIndices),...
                gammaOld2(validPitchIndices),false);
            gammaNew2(1:2,:) = gammaNewout;
            Gamma2{l} = gammaNew2(1:2,:);
        else 
            ll = l-1;
            [R1, phiout, psiout, alpha1, gammaOld1, gammaNew1] = ...
            computeGammaMultipleSinus(...
                R1(1:end-1, validPitchIndices), ll+1, ...
                crossCorrelationVectors(1:2*l+1, validPitchIndices), ...
                a1(ll, validPitchIndices), ...
                phi1(1:l-2, validPitchIndices), ...
                psi1(1:l-2, validPitchIndices), ...
                gammaOld1(1:l-2, validPitchIndices), ...
                gammaNew1(1:l-1, validPitchIndices), ...
                alpha1(validPitchIndices), ...
                true);
            psi1(1:l-1,1:size(psiout,2)) = psiout;
            phi1(1:l-1,1:size(phiout,2)) = phiout;
            [R2, phiout, psiout, alpha2, gammaOld2, gammaNew2] = ...
            computeGammaMultipleSinus(...
                R2(1:end-1,validPitchIndices), l, ...
                crossCorrelationVectors(1:2*l+1, validPitchIndices), ...
                a2(l-1, validPitchIndices), ...
                phi2(1:l-2, validPitchIndices), ...
                psi2(1:l-2, validPitchIndices), ...
                gammaOld2(1:l-2, validPitchIndices), ...
                gammaNew2(1:l-1, validPitchIndices), ...
                alpha2(validPitchIndices), ...
                false);
            psi2(1:l-1,1:size(psiout,2)) = psiout;
            phi2(1:l-1,1:size(psiout,2)) = phiout;
            Gamma1{l} = gammaNew1(1:l,:);
            Gamma2{l} = gammaNew2(1:l,:);
        end
    end
    
end