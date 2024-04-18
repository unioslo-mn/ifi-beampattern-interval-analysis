function value = getArrayCoeff(obj)  

% Calculate array coefficient from tapering, steering and coupling
%
% This function calculates the total effect of tapering, steeering
% and coupling and provides a single complex value or interval for 
% each element depending on whether the coupling coefficient is zero
% or non-zero.
%______________________________________________________________________
% USAGE        
%   
% _________________________________________________________________________
% NECESSARY ARGUMENT
%   obj       : object of biat.SensorArray type
% _________________________________________________________________________
% OPTIONS
% _________________________________________________________________________
% EXAMPLES
% _________________________________________________________________________
%
% Copyright (C) 2023 H. Arnestad and G. Gereb, BSD-3
% If you use this software, please cite it as in CITATION.cff
% Project: Beampattern Interval Analysis 
% Website: doi.org/10.5281/zenodo.6856232
% Contact: haavaarn@uio.no, gaborge@uio.no
% (More information in README.md and LICENSE.md.)
% _________________________________________________________________________

    arguments
        obj
    end

    cplR = obj.CouplingAmp;
    cplPh = obj.CouplingPhase;
    w = obj.TaperWeights;
    k = obj.SteeringVector;
    posX = obj.ElPosX;
    posY = obj.ElPosY;
    if obj.CouplingCoeff == 0
        % value  = zeros(obj.ElCount, 1); % Ac as complex number
        % for iCpl = 1:obj.ElCount % each row is one element + coupling to others
        %     for iEl = 1:obj.ElCount % each collumn is a coupling element
        %         value(iCpl) = value(iCpl) + cplR (iEl, iCpl) * w(iEl) * ...
        %                       exp(1j*(cplPh(iEl, iCpl)-...
        %                       (k(1) * posX(iEl) + k(2) * posY(iEl))));
        % 
        %     end 
        %     (cplR(:,iCpl).*w).' * exp(1j*(cplPh(:,iCpl)-[posX posY]*k));
        % end
        value = (cplR .* exp(1j*cplPh)).' * (w .* exp(-1j*[posX posY]*k));
    else
        % value  = repmat(ciat.CircularInterval( 0, 0), obj.ElCount,1);  
        % for iCpl = 1:obj.ElCount % each row is one element + idxing to others
        %     center = 0;
        %     radius = 0;
        %     for iEl = 1:obj.ElCount % each collumn is a idxing element
        %         if iCpl == iEl
        %             center = w(iEl) * exp(-1j*(k(1) * posX(iEl) + ...
        %                                        k(2) * posY(iEl)));
        %         else
        %             radius = radius + cplR(iEl, iCpl) * w(iEl);
        %         end
        %     end 
        % value(iCpl) = ciat.CircularInterval(center, radius);
        % end
        center = w .* exp(-1j * [posX posY] * k);
        radius = (cplR - diag(cplR).*eye(obj.ElCount)).' * w;
        value = ciat.CircularInterval(center, radius);
    end
end