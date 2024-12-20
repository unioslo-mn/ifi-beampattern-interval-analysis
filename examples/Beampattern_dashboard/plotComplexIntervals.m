function plotComplexIntervals()

%PLOT_CINTERVALS   Plot complex intervals
%
%  USAGE:   plot_cIntervals(obj)
% 
% _________________________________________________________________________
%  NECESSARY ARGUMENT
%     obj          = physical array objameter object 
%                    (see SimulationParameters.m)
% _________________________________________________________________________
%
% Copyright (C) 2023 Haavard Arnestad and Gabor Gereb, BSD-3 (LICENSE.md)
% If you use this software, please cite it as in CITATION.cff
% Project: Beampattern Interval Analysis (doi.org/10.5281/zenodo.6856232)
% Contact: haavaarn@uio.no, gaborge@uio.no (more in README.md)
% ________________________________________________________________________

% Set colors
plotCol = 'kcrbm';

% Get parameter struct
guicp = evalin('base','gui'); 

% Update complex interval plot
subplot(2,2,4);cla;hold on  
set(gca,'color','w');

% Plot complex intervals
for iPlt = 1:5
    if guicp.tb(iPlt).Value % guicp.plt.sel(idx)
        btr = (guicp.array.OrientError == 0);
        guicp.beampattern(iPlt).plotComplexIntervals('Color',plotCol(iPlt),...
                                                     'Backtracking',btr);
    end
end

% Update theta indicator in the beampattern plot
th_slider = guicp.beampattern(1).IncidenceAngle;
subplot(2,2,3);
ax = gca;
ax.Children(1).XData = [sin(th_slider),sin(th_slider)];

end