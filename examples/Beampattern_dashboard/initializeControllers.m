function initializeControllers()

    % Get parameter struct
    guicp = evalin('base','gui');

    % Initialize array
    guicp.array = biat.SensorArray('TaperType','kaiser');

    % Initialize beampatterns
    guicp.beampattern(1) = biat.BeamPattern(guicp.array,'nominal');
    guicp.beampattern(2) = biat.BeamPattern(guicp.array,'rectangular');
    guicp.beampattern(3) = biat.BeamPattern(guicp.array,'circular');
    guicp.beampattern(4) = biat.BeamPattern(guicp.array,'polygonal');

    % Initialize plot selectors
    guicp.plt.sel = ones(1,4);
    guicp.plt.legend = 1;
    guicp.plt.dynamicRange = -60;
    
    % Set parameter struct
    assignin('base','gui',guicp);

end

