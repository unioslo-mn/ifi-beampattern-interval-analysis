function value = getElementIntervals(obj,options)

% Calculate element intervals from tolerance errors and coupling
%
% Calculate element interval by forming the tolerance error interval
% from the given gain, phase, position and orientation errors
% (an annular sector represented by the selected interval type). 
% If the array has a coupling, then the tolerance error interval is
% multiplied by the coupling interval (only in polygonal mode).
%______________________________________________________________________
% USAGE        
%   obj.getElementIntervals(Name,Value)
% _________________________________________________________________________
% NECESSARY ARGUMENT
%   obj       : object of biat.BeamPattern type
% _________________________________________________________________________
% OPTIONS
%   getBlock  : if "0" the intervals are calculated for the elements
%               of the array(s) (default), if "1" it is calculated
%               for the elements of the block to get block intervals
% _________________________________________________________________________
% EXAMPLES
%   value = obj.getArrayInterval;
%   value = obj.getElementIntervals('getBlock',1);
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
           options.getBlock    (1,1)   {mustBeNumeric} = 0
        end

        % Extract array parameters according to given option
        if ~options.getBlock
            arrays = obj.Arrays;
            M = arrays(1).ElCount;
            N = length(arrays);
        else
            arrays = obj.Block;         
            M = obj.Block.ElCount;
            N = 1;
        end
        AcI = [arrays.ArrayCoeff];
        w_apod = [arrays.TaperWeights];
        
        % Extract parameters
        gI = [arrays.GainInterval];
        theta = obj.IncidenceAngle;
        psiI = [arrays.OrientInterval];
        
        % If array element intervals are calculated and there is a block
        % orientation error, add that to the element orientation error
        if ~options.getBlock && ~isempty(obj.Block)
            theta = theta - obj.Block.ElOrient';
            psiI = psiI + obj.Block.OrientError;
        end
        k = 2*pi / obj.Arrays(1).WaveLength * ([sin(theta); cos(theta)] );
        
        % Calculate amplitude and phase intervals
        EcI_r(1:M,1:N) = ciat.RealInterval();
        EcI_ph(1:M,1:N) = ciat.RealInterval();        
        posX = [arrays.PosXInterval];
        posY = [arrays.PosYInterval];
        phI = [arrays.PhaseInterval];
        for n = 1:N
            dTheta = theta(n) - psiI(:,n);
            EcI_r(:,n) = gI(:,n) .* arrays(n).getElDirectivity(dTheta);
            EcI_ph(:,n) = k(1,n)*posX(:,n) + k(2,n)*posY(:,n) + phI(:,n); 
        end
        
        % Combine with coupling
        if isnumeric(AcI)
            r_I = EcI_r .* abs(AcI);
            ph_I = EcI_ph + angle(AcI);
        end
        
        %         
        switch obj.Type
            case 'nominal'
                dk = k - [arrays.SteeringVector];
                phase = dk(1,:) .* [arrays.ElPosX] + ...
                        dk(2,:) .* [arrays.ElPosY];
                dTheta = obj.IncidenceAngle - [arrays.ElOrient]; % check this
                dir = arrays.getElDirectivity(dTheta);
                value = (dir .* w_apod) .* exp(1j*phase);
            case 'rectangular'
                if isnumeric(AcI)
                    value = repmat(ciat.RectangularInterval( ...
                        ciat.RealInterval(0), ciat.RealInterval(0)), M,N);
                    for n = 1:N
                        for m = 1:M
                            value(m,n) = ciat.RectangularInterval(...
                                            ciat.PolarInterval( ...
                                                r_I(m,n), ph_I(m,n) ) );
                        end
                    end
                else
                   value = nan();
                   warning('Rectangular type is not implemented for coupling.') 
                end
            case 'circular'
                if isnumeric(AcI)
                    value = repmat(ciat.CircularInterval( 0, 0), M,N);  
                    for n = 1:N
                        for m = 1:M
                            value(m,n) = ciat.CircularInterval( ...
                                            ciat.PolarInterval( ...
                                                r_I(m,n), ph_I(m,n) ) );
                        end
                    end
                else
                   value = nan();
                   warning('Circular type is not implemented for coupling.') 
                end
            case 'polygonal'
                inclusive = 1;
                tolerance = obj.PolygonTolerance;
                if isnumeric(AcI)
                    value = repmat(ciat.PolygonalInterval(0), M,N);
                    for n = 1:N
                        for m = 1:M
                            value(m,n) = ciat.PolygonalInterval( ...
                                            ciat.PolarInterval( ...
                                                r_I(m,n), ph_I(m,n) ) ,...
                                            'tolerance',tolerance ) ;
                        end
                    end
                else
                    value = repmat(ciat.PolygonalInterval(0), M,N);
                    for n = 1:N
                        for m = 1:M
                            EcI = ciat.PolarInterval( ...
                                                EcI_r(m,n), EcI_ph(m,n) );
                            value(m,n) = ciat.PolygonalInterval(EcI,...
                                                            AcI(m,n), ...
                                            'tolerance', tolerance);
                        end
                    end
                end
                case 'polyarcular'
                if isnumeric(AcI)
                    value(M,N) = ciat.PolyarcularInterval;  
                    for n = 1:N
                        for m = 1:M
                            value(m,n) = ciat.PolyarcularInterval( ...
                                            ciat.PolarInterval( ...
                                                r_I(m,n), ph_I(m,n) ) );
                        end
                    end
                else
                   value = nan();
                   warning('Polyarcular type is not implemented for coupling.') 
                end
            otherwise
                value = nan();
                warning('No element interval is implemented for the give type.')
        end
end

