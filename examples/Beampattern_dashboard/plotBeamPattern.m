function plotBeamPattern()

%PLOT_BEAMPATTERN   Plot beampattern
%
%  USAGE:   plot_beampattern_interval(obj)
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

% Get parameter struct
guicp = evalin('base','gui'); 

% Initialize plot
subplot(2,2,3);cla;
set(gca,'color','w');

% Set colors
plotCol = 'kcrb';
plotName = {'nom','rect','circ','poly'};

% Plot beampatterns
plt = [];
for iPlt = 1:4
    if guicp.tb(iPlt).Value % guicp.plt.sel(iPlt)
        btr = (guicp.array.OrientError == 0);
        plt = [plt; guicp.beampattern(iPlt).plotBeamPattern(...
                                        'Color',plotCol(iPlt),...
                                        'Backtracking',btr)]; 
    else
        for iLine = 1:3
           plt = [plt ; plot(0,0,'Color',plotCol(iPlt),...
                                 'DisplayName',plotName{iPlt})];
        end
    end
end

% Add bar indicating chose incident angle
th_slider = guicp.beampattern(1).IncidenceAngle;
max_value = max([max(plt(3).YData),...
                 max(plt(5).YData),...
                 max(plt(8).YData),...
                 max(plt(11).YData)]);
plot([sin(th_slider),sin(th_slider)],...
     [guicp.plt.dynamicRange,max_value],'k:','DisplayName','angle')

ylim([guicp.plt.dynamicRange,max_value])

if guicp.plt.legend
    hL = legend(plt([3,4,7,10]),'AutoUpdate','off');
    set(hL,'Units', 'normalized', 'Position', [0.505,0.365,0.005,0.0663]);
    guicp.plt.legend = 0;
end

end

