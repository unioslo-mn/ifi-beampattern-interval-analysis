function value = getArrayInterval(obj)

% Calculate array interval(s) (SubarrayInterval × BlockInterval)
%
% This function calculates the array intervals for each given 
% array by multiplying the sub-array interval with the block
% interval. If there is no block array given, the array 
% intervals equal the sub-array intervals
%______________________________________________________________________
% USAGE        
%   intervals = obj.getArrayinterval
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


    % Get the sub-array interval as input for multiplication
    value = obj.SubArrayInterval;
    
    % Return NaN if input is NaN
    if isnumeric(value(1)) && isnan(value(1))
        return
    end
    
    % Calculate the product of each sub-array interval with the
    % corresponding block element interval
    if ~isempty(obj.Block)
        blInt = obj.BlockInterval;
        N = obj.ArrayCount;
        
        % Calculate product
        for n = 1:N
            value(n) = value(n) .* blInt(n);
        end
    end
end

