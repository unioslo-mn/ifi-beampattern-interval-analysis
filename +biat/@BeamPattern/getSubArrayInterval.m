function value = getSubArrayInterval(obj)

% Calculate sub-array interval(s) from element intervals
%
% This function calculates the sub-array intervals for each
% in the Arrays property. If no block orientation error is
% present, the sub-array interval is the sum of the element
% intervals. If there is a block orientation error then the
% sum of element intervals are calculated at each beam angle
% within the orientation-error-wide vicinity of the given
% incidence angle, and the union of these sum is taken to
% get the sub-array interval.
%______________________________________________________________________
% USAGE        
%   value = obj.getSubArrayInterval();
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


    elInt = obj.ElementIntervals;
    nArr = obj.ArrayCount;
    
    % If input is NaN return NaN
    if isnumeric(elInt(1)) && isnan(elInt(1))
        return
    end
    
    % Calculate sum of element intervals
    if strcmp(obj.Type,"nominal") || isempty(obj.Block) || ...
                                     obj.Block.OrientError == 0
        % If there is no block orientation error calculate the sum
        % of element intervals
        value = getElIntSum(elInt);
    else
    % If there is block orientation error, calculate element interval sum
    % for each angle in the interval and then get their union
       
        % Find beam index interval
        cntBeam = obj.BeamCount;
        if obj.BeamIndex ~= 0
            idxBeam = obj.BeamIndex;
        else
            [~,idxBeam] = min(abs(obj.BeamAngles - obj.BeamAngle));
        end
        difBeam = ceil(obj.Block.OrientError / obj.BeamResolution);
        intBeam = max(1,idxBeam-difBeam) : min(idxBeam+difBeam,cntBeam);    

        % Calculate sub-array intervals for the index interval
        subInt(1:length(intBeam),1:nArr) = elInt(1);
        for idx = 1:length(intBeam)
            obj.BeamIndex = intBeam(idx);
            subInt(idx,:) = getElIntSum(obj.ElementIntervals);
        end

        % Calculate the union of array intervals along the beam indeces
        value(1:nArr) = ciat.PolygonalInterval(0);
        for iArr = 1:nArr
            value(iArr) = union( subInt(:,iArr) );
        end
    end
end

%% 

function value = getElIntSum(elInt)
    [M,N] = size(elInt);
    value(1:N) = elInt(1,1);
    for n = 1:N
        value(n) = elInt(1,n);
        for m = 2:M
            value(n) = value(n) + elInt(m,n);
        end
    end
    
end

