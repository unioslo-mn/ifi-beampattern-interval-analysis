function [E_extr,C_extr] = backtrackComplexIntervals(obj,options)

% Backtrack element interval points to errors and coupling
%
% The element intervals are calculated as the product of 
% the complex tolerance error interval (annular sector) 
% and the coupling interval (circle). This function
% backtracks the selected element interval points 
% (the ones contributing to minimum or maximum array power)
% to the contributing points in the tolerance error interval
% and the array interval.
%______________________________________________________________________
% USAGE        
%   [E_extr,C_extr] = backtrackComplexIntervals(Name, Value) 
% _________________________________________________________________________
% NECESSARY ARGUMENT
%   obj       : object of biat.BeamPattern type
% _________________________________________________________________________
% OPTIONS
%   getMax     : "1" to backtrack maximum power (default),
%                 "0" to backtrack minimum power point.
% _________________________________________________________________________
% EXAMPLES
%   [E_extr,C_extr] = backtrackComplexIntervals;
%   [E_extr,C_extr] = backtrackComplexIntervals('getMax',0);
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
       options.getMax       (1,1)   {mustBeNumeric} = 1
    end

    % Initialize output variables
    M = obj.Arrays.ElCount;
    N = obj.ArrayCount;
    E_extr = zeros(M,N);
    C_extr = zeros(M,M,N);

    if strcmp(obj.Type,'polygonal') 
        % Get extreme point to backtrack
        if options.getMax
            extreme_points = obj.MaxPowerPoints;
        else
            extreme_points = obj.MinPowerPoints;
        end
        
        % Get complex intervals contributing to the extreme point
        for n = 1:N
           [E_extr(:,n),C_extr(:,:,n)] = getCInt(obj.Arrays(n),...
                                                obj.IncidenceAngle,...
                                                extreme_points);
        end
    else
        warning('No backtracking available for this type')
    end
end

%% 

function [E_extr,C_extr] = getCInt(array,incAngle,extreme_points)

    % Construct EcI
    k = 2*pi / array.WaveLength * ([sin(incAngle); cos(incAngle)] );
    dPsi = incAngle - [array.OrientInterval]; 
    EcI_ph = k(1) * [array.PosXInterval] + ...
             k(2) * [array.PosYInterval]  + ...
             [array.PhaseInterval] ;
    EcI_r = [array.GainInterval] .* ...
            array.getElDirectivity(dPsi);

    % Get nominal coupling
    AcI = array.ArrayCoeff;
    if isnumeric(AcI)
       AcI_center = AcI; 
    else
        AcI_center = [AcI.Center];
    end

    % Get other parameters
    M = array.ElCount;
    w_apod = array.TaperWeights;
    k_s = array.SteeringVector;
    posX = array.ElPosX;
    posY = array.ElPosY;
    cplR = array.CouplingAmp;
    
    % Initialize output arrays
    E_extr = zeros(M,1);
    C_extr = zeros(M);

    for iInt = 1:M
        z_marked = extreme_points(iInt);
        E_inv = ciat.PolarInterval( ...
                            ciat.RealInterval( ...
                                abs(z_marked)/(EcI_r(iInt).Supremum+eps), ...
                                abs(z_marked)/(EcI_r(iInt).Infimum+eps)), ...
                            ciat.RealInterval(...
                                angle(z_marked)-EcI_ph(iInt).Supremum,...
                                angle(z_marked)-EcI_ph(iInt).Infimum) ); 

        % Backtrack element intervals
        % this block could be made analytical... (CONSIDER REPLACE)
        circSmpCnt = 100; % Replace this with polygon tolerance
        E_inv_sample = E_inv.sample(circSmpCnt);
        E_inv_sample = E_inv_sample{:};
        [~, min_idx] = min( abs( E_inv_sample - AcI_center(iInt) ) );
        A_marked = E_inv_sample(min_idx);

        E_marked = z_marked / A_marked;
        E_extr(iInt) = E_marked;

        % Backtrack coupling intervals
        angle_c = angle( E_inv_sample(min_idx) - w_apod(iInt) * exp(-1j * ...
                        (k_s(1)*posX(iInt) + k_s(2)*posY(iInt))) ) +...
                        (k_s(1)*posX + k_s(2)*posY);
        C_extr(:,iInt) = cplR(:,iInt) .* exp(1j*angle_c);
        C_extr(iInt,iInt) = 1;  % one on the main diagonal
    end



end

