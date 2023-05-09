function value = backtrackBeamPattern(obj,options)

% Get beampattern for backtracked error pattern
%
% Calculate the beampattern from the backtracked errors at 
% a reference angle. The beampattern is computed for all 
% obj.BeamAngles specified in the beampattern object.
% _________________________________________________________________________
% USAGE        
%   p = obj.backtrackBeamPattern(Name, Value) 
% _________________________________________________________________________
% NECESSARY ARGUMENT
%   obj       : object of biat.BeamPattern type
% _________________________________________________________________________
% OPTIONS
%   getMax     : "1" to backtrack maximum power (default),
%                 "0" to backtrack minimum power point.
% _________________________________________________________________________
% EXAMPLES
%   points = obj.backtrackBeamPattern;
%   points = obj.backtrackBeamPattern('getMax',0);
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

    % Preallocate memory
    value = zeros(length(obj.BeamAngles), 1);

    % Get other parameters
    M = obj.Arrays.ElCount;
    k_s = obj.Arrays.SteeringVector;
    posX = obj.Arrays.ElPosX;
    posY = obj.Arrays.ElPosY;
    w_apod = obj.Arrays.TaperWeights;
    [E_extr,C_extr] = backtrackComplexIntervals(obj,'getMax',options.getMax);
    
    %% Calculate complete beampattern for the backtracked errors
    thetas = obj.BeamAngles;
    theta_ref = obj.IncidenceAngle;
    k_abs = 2*pi / obj.Arrays(1).WaveLength;
    k_ref = k_abs * ([sin(theta_ref), cos(theta_ref)] );

    for iBeam = 1 : length(thetas) % Loop over angles
        k = k_abs * ([sin(thetas(iBeam)), cos(thetas(iBeam))] );
        dk = (k - k_ref);
        phase_corr = sum(dk.*[posX, posY],2);

        BP = 0;

        for c = 1:M
            sum_m = 0;
            for m = 1:M
                sum_m = sum_m + w_apod(m) * C_extr(m,c) * exp( -1j * ...
                                        ( k_s(1)*posX(m) + k_s(2)*posY(m) ));
            end
            BP = BP + E_extr(c)*exp(1j*phase_corr(c)) * sum_m;
        end
        value(iBeam) = abs(BP).^2;
    end
end

