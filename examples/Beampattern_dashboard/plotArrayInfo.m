function plotArrayInfo()
%PLOT_ARRAY_INFO   Plot physical array layout
%
%  USAGE:   plot_array_info(obj)  
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

    % Reset plot
    ax = subplot(2,2,2);
    delete(ax)
    guicp.pltPA = subplot(2,2,2);
    assignin('base','gui',guicp);
        
    % Initialize plot
    guicp.array.plot();
    
    % Set the background of the beampattern and complex interval plots to
    % grey
    subplot(2,2,3);
    set(gca,'color',[.9 .9 .9]);
    subplot(2,2,4);
    set(gca,'color',[.9 .9 .9]);
end

