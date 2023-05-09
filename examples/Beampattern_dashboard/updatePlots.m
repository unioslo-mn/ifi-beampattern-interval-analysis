function a = updatePlots(h)
    % Check button value
    if get(h,'Value') == 1
        % Switch button string to busy
        set(h,'Value',1);

        % Update beampattern plot
        plotBeamPattern()
        
        % Update complex interval plot
        plotComplexIntervals()
        
        % Update error pattern plot
        plotErrorPattern();

        % Restore button string to update
        set(h,'Value',0);
    end
end