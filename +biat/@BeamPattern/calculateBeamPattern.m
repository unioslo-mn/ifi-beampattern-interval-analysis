function value = calculateBeamPattern(obj)

% Calculate beampattern of the Arrays at the BeamAngles using Type
%
% This function samples the beampattern by cycling through the 
% BeamAngles and storing the power interval calculated by the 
% object. It also calculates an approximate beampatter if the
% type is nominal. Also if there is a block orientation error,
% the function calculates the element intervals for all angles
% before it would form the unions for the subarray intervals to
% avoid redundant calculation.
%______________________________________________________________________
% USAGE        
%   power = obj.calculateBeamPattern
% _________________________________________________________________________
% NECESSARY ARGUMENT
%   obj       : object of biat.BeamPattern type
% _________________________________________________________________________
% OPTIONS
% _________________________________________________________________________
% EXAMPLES
%   P_nom = bp_nom.calculateBeamPattern;
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

    cntBeam = obj.BeamCount;
    nArr = obj.ArrayCount;
    tol = obj.PolygonTolerance;
    value = zeros(cntBeam, 1+(obj.Type ~= 'nominal') ); 
    
    % Block calculation for non-implemented cases
    if (strcmp(obj.Type,'rectangular') || strcmp(obj.Type,'circular')) && ...
        obj.Arrays.CouplingCoeff ~= 0
        warning('Circular type is not implemented for coupling.')
        return
    end
    
    % Calculate beampattern according to type
    if strcmp(obj.Type,'nominal')
        % Nominal beampattern
        for idxBeam = 1 : cntBeam  
            obj.BeamIndex = idxBeam;
            value(idxBeam) = obj.PowerInterval;
        end
    else
        % Beampattern interval
        if isempty(obj.Block) || obj.Block.OrientError == 0
            % Without block orientation error
            for idxBeam = 1 : cntBeam  
                obj.BeamIndex = idxBeam;
                value(idxBeam,:) = obj.PowerInterval.Bounds';
            end
        else
            % With block orientation error
            if strcmp(obj.Type,'polygonal')
                % Calculate sub-array intervals for each angle
                % (switch of block orientation time during this phase)
                blockPsiErr = obj.Block.OrientError;
                obj.Block.OrientError = 0;
                subInt(1:cntBeam,1:nArr) = obj.ElementIntervals(1);
                for idxBeam = 1 : cntBeam  
                    obj.BeamIndex = idxBeam;
                    subInt(idxBeam,:) = obj.SubArrayInterval;
                end
                obj.Block.OrientError = blockPsiErr;
                % Calculate the sub-array interval unions within the 
                % block orientation error
                difBeam = ceil(obj.Block.OrientError / obj.BeamResolution);
                for idxBeam = 1 : cntBeam
                    % Get block intervals for the selected beam index
                    obj.BeamIndex = idxBeam;
                    blInt = obj.BlockInterval;
                    
                    % Find beam index interval for block orientation error
                    intBeam = max(1,idxBeam-difBeam) : ...
                              min(idxBeam+difBeam,cntBeam);
                          
                    % Calculate the union, multiply by the block interval
                    % and sum the intervals to get the total
                    totInt = ciat.PolygonalInterval(0);
                    for iArr = 1:nArr
                        tiltInt = union( subInt(intBeam,iArr) );
                        arrInt = tiltInt .* blInt(iArr);
                        totInt = totInt + arrInt;
                    end
                    ampInt = abs(totInt);
                    powerInt = ampInt * ampInt;
                    value(idxBeam,:) = powerInt.Bounds';
                end
            else
               error('Block orientation error is not implemented for this type') 
            end
        end
    end
    
    % Calculate approximate beampattern for nominal type
    if strcmp(obj.Type,'nominal')
        circle_bound = sqrt( mean([obj.Arrays.GainError],'all')^2 + ...
                             mean([obj.Arrays.PhaseError],'all')^2);
        r_tot_bound = circle_bound + 2 * mean([obj.Arrays.CouplingCoeff]);
        P_apx_max = (sqrt(value) + r_tot_bound).^2;
        A_apx_min = (sqrt(value) - r_tot_bound);
        A_apx_min(A_apx_min <0 ) = 0;
        P_apx_min = A_apx_min.^2;
            
        % Sort nominal and approximate beampatterns into a single array 
        value = [value , P_apx_min , P_apx_max];
    end

end

