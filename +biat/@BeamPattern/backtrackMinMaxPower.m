function value = backtrackMinMaxPower(obj,options)

% Backtrack array interval point to element interval points
%
% Array intervals are calculated as the sum of element 
% intervals, this function the boundary point on each 
% element interval that contributed to the minimum 
% or maximum power point.
% _________________________________________________________________________
% USAGE        
%   p = obj.backtrackMinMaxPower(Name, Value) 
% _________________________________________________________________________
% NECESSARY ARGUMENT
%   obj       : object of biat.BeamPattern type
% _________________________________________________________________________
% OPTIONS
%   getMax     : "1" to backtrack maximum power (default),
%                 "0" to backtrack minimum power point.
% _________________________________________________________________________
% EXAMPLES
%   points = obj.backtrackMinMaxPower;
%   points = obj.backtrackMinMaxPower('getMax',0);
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
    
    if strcmp(obj.Type,'polygonal')
        M = obj.Arrays(1).ElCount;
        N = obj.ArrayCount;
        value = zeros(M, N);
        arrayInt = obj.ArrayInterval;
        elemInt = obj.ElementIntervals;
        for n = 1:N
            points = arrayInt(n).Points;
            if options.getMax
                [~, extreme_idx] = max(abs(points));
                extreme_angle = wrapTo2Pi( angle( points (extreme_idx)));
            else
                if inpolygon(0,0,real(points),...
                                 imag(points))
                    value = [];
                    warning('Minimum power cannot be backtracked.')
                    return
                end
                [~, extreme_idx] = min(abs(points));
                extreme_angle = wrapTo2Pi( angle( ...
                                    points(extreme_idx)) - pi );
            end

            cpxIntv = elemInt(:,n);
            for m = 1:M % loop over all obj.ElementIntervals
                value(m,n) = findPoint(cpxIntv(m),extreme_angle);
            end 
        end
    else
        value = nan(obj.Arrays(1).ElCount,obj.ArrayCount);
        warning('No backtracking available for this type')
    end
end

%%
function point = findPoint(cpxIntv,extreme_angle)
    eps10 = 10*eps;
    v = cpxIntv.Points;
    prev_angle = wrapTo2Pi( angle( v(1) - v(end) )- pi/2); % init prev angle
    point = v(end); % just assume it is the final point to avoid extra condition later
    for j = 1 : length(v)-1
        current_angle = wrapTo2Pi( angle( v(j+1) - v(j) )- pi/2);
        if (current_angle - prev_angle) < 0 % we crossed 0
            cond1 = (extreme_angle >= prev_angle - eps10) && ...
                (extreme_angle <= 2*pi+current_angle+eps10);
            cond2 = (extreme_angle >= prev_angle -2*pi-eps10) && ...
                    (extreme_angle <= current_angle+eps10);
            cond = cond1 || cond2;
        else
            cond = (extreme_angle >= prev_angle-eps10) && ...
                   (extreme_angle <= current_angle+eps10);
        end

        if cond == true
            point = v(j);
            break;
        end
        prev_angle = current_angle;
    end
end