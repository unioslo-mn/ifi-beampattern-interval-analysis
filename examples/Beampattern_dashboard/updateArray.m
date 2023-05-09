function updateArray(h)
    % Get event details
    name = func2str(h.Callback);
    value = h.Value;

    % Get parameter struct
    guicp = evalin('base','gui');
    
    % Update selected parameter
    switch name
        % Control panel
        case 'cbNel'
            guicp.array.ElCount = round(value);
            guicp.val.Nel.String = formatValueString(name,value);
        case 'cbRcv'
            if value == 0
                guicp.array.Curvature = 0;
                guicp.val.Rcv.String = formatValueString(name,value);
            else
                guicp.array.Curvature = guicp.array.ElCount * ...
                                        guicp.array.ElPitch/(2*pi)/value;
                guicp.val.Rcv.String = formatValueString(name,value);
            end
            
        case 'cbEld'
            guicp.array.ElPitchRatio = value;
            guicp.val.Eld.String = formatValueString(name,value);
        case 'cbEls'
            guicp.array.ElDiameterRatio = value;
            guicp.val.Els.String = formatValueString(name,value);
        case 'cbSte'
            guicp.array.SteeringAngle = value * pi/2;        
            guicp.val.Ste.String = formatValueString(name,value);
        case 'cbWgh'
            guicp.array.TaperParam = value;
            guicp.val.Wgh.String = formatValueString(name,value);
        case 'cbEGn'
            guicp.array.GainError = value;
            guicp.val.EGn.String = formatValueString(name,value);
        case 'cbEPh'
            guicp.array.PhaseError = value;
            guicp.val.EPh.String = formatValueString(name,value);
        case 'cbEPs'
            guicp.array.OrientError = value;
            guicp.val.EPs.String = formatValueString(name,value);
        case 'cbEPx'
            guicp.array.PosXError = value * guicp.array.ElPitch;
            guicp.val.EPx.String = formatValueString(name,value);
        case 'cbEPy'
            guicp.array.PosYError = value * guicp.array.ElPitch;
            guicp.val.EPy.String = formatValueString(name,value);
        case 'cbCpl'
            guicp.array.CouplingCoeff = value;
            guicp.val.Cpl.String = formatValueString(name,value);
            
    end
    
    % Set plot selectors according to non-implement constraints
    if guicp.array.CouplingCoeff ~= 0
       guicp.tb(2).Value = 0;
       guicp.tb(3).Value = 0;
    else
       guicp.tb(2).Value = guicp.plt.sel(2); 
       guicp.tb(3).Value = guicp.plt.sel(3); 
    end
    
    % Set parameter struct
    assignin('base','gui',guicp);
    
    % Update physical plot
    plotArrayInfo()
end