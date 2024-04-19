function value = getElDirectivity(obj,angle)

% Calculate element directivity for a given incident angle
%
% This function calculates the element directivity value 
% or interval for a given angle or angle interval.
%______________________________________________________________________
% USAGE        
%   directivity = obj.getElDirectivity(angle)
% _________________________________________________________________________
% NECESSARY ARGUMENT
%   obj       : object of biat.SensorArray type
%   angle     : incidence angle value or real interval [rad]
% _________________________________________________________________________
% OPTIONS
% _________________________________________________________________________
% EXAMPLES
%   dir = array.getElDirectivity(pi/4);
%   dir = array.getElDirectivity(bp_pol.IncidenceAngle - ...
%                                bp_pol.Arrays.OrientInterval);
% _________________________________________________________________________
%
% Copyright (C) 2023 H. Arnestad and G. Gereb, BSD-3
% If you use this software, please cite it as in CITATION.cff
% Project: Beampattern Interval Analysis 
% Website: doi.org/10.5281/zenodo.6856232
% Contact: haavaarn@uio.no, gaborge@uio.no
% (More information in README.md and LICENSE.md.)
% _________________________________________________________________________
    
    % If input angles are intervals, split to lower and upper bounds
    % and initialize output accordingly
    N = length(obj);
    M = size(angle,1);
    if isnumeric(angle)
        cDim = 1;
        value = zeros(M,N);
    else
        cDim = 2;
        angle = [angle.Infimum , angle.Supremum];
        value(M,N) = ciat.RealInterval(0,0);
    end
    
    % Calculate element directivities for each array in the input object
    for iArr = 1:length(obj)
        value(:,iArr) = getElDir(obj(iArr),angle,cDim,M);
    end
end


%% Local function to get element directivity of an array

function elDir = getElDir(array,angle,cDim,M)

    elRad = array.ElDiameter / 2;
        
    if array.ElDiameter > 0
        % Calculate element directivity for all elements
        elDir = zeros(M,cDim);
        for iDim = 1:cDim
            % Calculate element directivity from -90 to 90 degrees
            kxy = 2 * pi * sin(angle(:,iDim)) / array.WaveLength;
            elDir(:,iDim) = real(2 * pi * elRad * ...
                            besselj(1, kxy * elRad) ./ kxy);
            elDir((kxy == 0),iDim) = pi * elRad^2;

            % Normalize element directivity to one
            elDir(:,iDim) = elDir(:,iDim) / (pi * elRad^2); 

            % Add tapering: set to zero beyond taper range, 
            % get taperint indices
            elDir(abs(angle(:,iDim)) > array.TaperAngles(2),iDim)=0;     
            taperInds = find(...
                    abs(angle(:,iDim)) >= array.TaperAngles(1) ...
                    & abs(angle(:,iDim)) <= array.TaperAngles(2));      
            elDir(taperInds,iDim) = elDir(taperInds,iDim) .* ...
                            (1 + cos(pi * (abs(angle(taperInds,iDim)) - ...
                            array.TaperAngles(1)) / ...
                            diff(array.TaperAngles)))/2;
        end

        % Reshape output and form output interval if necessary
        if cDim > 1
            elDirInf = min( elDir,[],2 );
            elDirSup = max( elDir, [], 2 );
            elDirSup( angle(:,1)<=0 & angle(:,2)>=0 )  = 1;
            elDir = ciat.RealInterval(elDirInf,elDirSup);
        end
    else
        % Element directivity is one for point-like elements
        elDir = ones(M,1);
    end
end