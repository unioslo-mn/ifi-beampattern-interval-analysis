function value = getTotalInterval(obj)

% Calculate total interval as the sum of array intervals
%
% This function calculates the sum of array intervals (or
% nominal values) to get the total complex interval (value)
%______________________________________________________________________
% USAGE        
%   value = obj.getTotalInterval();
% _________________________________________________________________________
% NECESSARY ARGUMENT
%   obj       : object of biat.BeamPattern type
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

    % Extract array intervals
    arrInt = obj.ArrayInterval;
    tol = obj.PolygonTolerance;
    
    % Nominal
    if strcmp(obj.Type,"nominal") 
        if ~isnumeric(arrInt)
            arrInt = arrInt.Midpoint;
        end
        value = sum(arrInt);
    else 
    % Interval
        value = arrInt(1);
        for n = 2:obj.ArrayCount
            value = value + arrInt(n);
        end
    end
end

