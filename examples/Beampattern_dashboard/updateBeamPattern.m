function updateBeamPattern(h)

    % Get event details
    name = func2str(h.Callback);
    value = h.Value;

    % Get parameter struct
    guicp = evalin('base','gui');
    
    % Get axis of beampattern plot
    subplot(2,2,3);
    ax = gca();
    
    % Update selected parameter
    switch name
        % Tickbox panel
        case 'cbTbN'
            guicp.plt.sel(1) = value;
        case 'cbTbR'
            guicp.plt.sel(2) = value;
        case 'cbTbC'
            guicp.plt.sel(3) = value;
        case 'cbTbM'
            guicp.plt.sel(4) = value;
        case 'cbThi'
            theta = round((asin(value)/pi*2+1)/...
                           2*(guicp.beampattern(1).BeamCount-1)+1);
            [guicp.beampattern.BeamIndex] = deal(theta);
            guicp.val.Thi.String = num2str(rad2deg(...
                                    guicp.beampattern(1).IncidenceAngle));
        case 'cbDRn'
            [guicp.plt.dynamicRange] = deal(value);
            subplot(2,2,3);
            ax = gca;
            ax.YLim(1) = value;
    end
    
    % Update beampattern plot visibility
    for plt = 1:4
       for lin = 0:2 
            ax.Children(end-(plt-1)*3-lin).Visible = guicp.plt.sel(plt);
       end
    end
    
    % Set parameter struct
    assignin('base','gui',guicp);
    
    % Update complex interval plot
    subplot(2,2,4);cla;hold on
    plotComplexIntervals();
    plotErrorPattern();
end