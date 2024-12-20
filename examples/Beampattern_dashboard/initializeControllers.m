function initializeControllers()

    % Get parameter struct
    guicp = evalin('base','gui');

    % Initialize array
    guicp.array = biat.SensorArray('TaperType','kaiser',...
                                   'GainError',0.02,...
                                   'PhaseError',deg2rad(5));

    % Initialize beampatterns
    guicp.beampattern(1) = biat.BeamPattern(guicp.array,'nominal');
    guicp.beampattern(2) = biat.BeamPattern(guicp.array,'rectangular');
    guicp.beampattern(3) = biat.BeamPattern(guicp.array,'circular');
    guicp.beampattern(4) = biat.BeamPattern(guicp.array,'polygonal');
    guicp.beampattern(5) = biat.BeamPattern(guicp.array,'polyarx');

    % Initialize plot selectors
    guicp.plt.sel = ones(1,5);
    guicp.plt.legend = 1;
    guicp.plt.dynamicRange = -60;
    
    % Set parameter struct
    assignin('base','gui',guicp);

end

