function plotComplexIntervals(obj,options)

% Plot complex intervals (element, sub-array, block, array, total)
%
% This function blocks the intervals calculated for the selected
% beam angle. The plotted intervals include the element, 
% sub-array, block, array and total intervals. It is possible
% to set the color and style of the lines and markers and to 
% turn backtracking on or off.
%______________________________________________________________________
% USAGE        
%   obj.plotComplexIntervals(Name,Value)
% _________________________________________________________________________
% NECESSARY ARGUMENT
%   obj       : object of biat.BeamPattern type
% _________________________________________________________________________
% OPTIONS
%   Color      : same as the Color property of the Matlab 
%                plot function
%   LineStyle  : same as the LineStyle property of the Matlab 
%                plot function
%   LineWidth  : same as the LineWidth property of the Matlab 
%                plot function
%   Marker     : same as the Marker property of the Matlab 
%                plot function
%   MarkerSize : same as the MarkerSize property of the Matlab 
%                plot function
%   Backtracking: set to "1" to plot backtracked beampattern
%                (default), and "0" to skip it   
% _________________________________________________________________________
% EXAMPLES
%   obj.plotComplexIntervals('Color',plotCol(iPlt),...
%                            'Backtracking',true);
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
       options.Color        (1,:)   {mustBeText}    = 'b'
       options.LineStyle    (1,:)   {mustBeText}    = '-'
       options.LineWidth    (1,1)   {mustBeNumeric} = 2
       options.Marker       (1,:)   {mustBeText}    = '.'
       options.MarkerSize   (1,1)   {mustBeNumeric} = 10
       options.Backtracking (1,1)   {mustBeNumericOrLogical} = 1
    end
    
    % Get dimensions
    [M,N] = size(obj.ElementIntervals);

    % Initialize figure
    hold on;
    title('Complex intervals')

    %% Plot element and block intervals
    set(gca,'LineWidth',2)
    
    elInt = obj.ElementIntervals;
    blInt = obj.BlockInterval;
    if isnumeric(elInt)
       plot(real(elInt),imag(elInt),'x','Color',options.Color);
    else
        for n = 1:N
            for m = 1:M
                elInt(m,n).plot('Color',options.Color,...
                                'LineWidth',options.LineWidth,...
                                'LineStyle',options.LineStyle,...
                                'Marker',options.Marker,...
                                'MarkerSize',options.MarkerSize);
            end
            if ~isempty(blInt)
                blInt(n).plot('Color',options.Color,...
                                'LineWidth',options.LineWidth+1,...
                                'LineStyle',options.LineStyle,...
                                'Marker',options.Marker,...
                                'MarkerSize',options.MarkerSize); 
            end
        end
    end

    %% Plot array interval
    arrInt = obj.ArrayInterval;
    if isnumeric(arrInt)
        quiver(0,0,real(arrInt),imag(arrInt),options.Color,...
                            'LineWidth',options.LineWidth,...
                            'AutoScale','off','DisplayName','nom.');
    else
        for n = 1:N
            arrInt(n).plot( 'Color',options.Color,...
                            'LineWidth',options.LineWidth+2,...
                            'LineStyle',options.LineStyle,...
                            'Marker',options.Marker,...
                            'MarkerSize',options.MarkerSize);
        end
    end

    %% Plot backtracked points
    if strcmp(obj.Type,'polygonal')
        if options.Backtracking == 1
            if obj.Arrays(1).OrientError == 0
                errPat = obj.MaxPowerPoints;
                maxPnt = sum(errPat);
                plot(real(errPat),imag(errPat),'o','Color',options.Color);
                plot(real(maxPnt),imag(maxPnt),'o','Color',options.Color);
            else
                warning(['Backtracking is not implemented for arrays ',...
                            'with non-zero orientation error.']) 
            end
        end
    end

    % Finalized figure
    grid on
    axis equal
    hold off
    
end

