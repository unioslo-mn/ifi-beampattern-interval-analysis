function plotErrorPattern()  

%PLOT_ERROR_PATTERN   Plot backtracked error pattern
%
%  USAGE:   plot_error pattern(obj)  
% 
% _________________________________________________________________________
%  NECESSARY ARGUMENT
%     obj          = physical array parameter object 
%                    (see SimulationParameters.m)
% _________________________________________________________________________
%
% Copyright (C) 2023 Haavard Arnestad and Gabor Gereb, BSD-3 (LICENSE.md)
% If you use this software, please cite it as in CITATION.cff
% Project: Beampattern Interval Analysis (doi.org/10.5281/zenodo.6856232)
% Contact: haavaarn@uio.no, gaborge@uio.no (more in README.md)
% ________________________________________________________________________

    % Get parameter struct 
    guicp = evalin('base','gui');    

    if guicp.array.OrientError== 0

         % Reset plot
        ax = subplot(2,2,2);
        delete(ax)
        guicp.pltPA = subplot(2,2,2);
        assignin('base','gui',guicp);

        guicp.beampattern(4).plotErrorPattern();
    end 
end