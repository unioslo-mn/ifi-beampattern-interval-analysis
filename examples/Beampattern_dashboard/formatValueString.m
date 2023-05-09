function valstr = formatValueString(ctrl,val)

    guicp = evalin('base','gui');

    switch ctrl
        case 'cbNel'
            valstr = [num2str(round(val))];
        case 'cbRcv'
            if val == 0
                valstr = 'Linear';
            elseif val == 1
                valstr = 'Circular';
            else
                valstr = [num2str(guicp.array.ElCount * ...
                                  guicp.array.ElPitch/(2*pi)/val,2), 'm'];
            end
        case 'cbEld'
            valstr = [num2str(val,2), 'm'];
        case 'cbEls'
            valstr = [num2str(round(val*100)), '%'];
        case 'cbThi'
            valstr = [num2str(rad2deg(val),2), '°'];
        case 'cbSte'
            valstr = [num2str(rad2deg(val* pi/2),2), '°'];
        case 'cbWgh'
            valstr = [num2str(round(val))];
        case 'cbEGn'
            valstr = [num2str(round(val*100)), '%'];
        case 'cbEPh'
            valstr = [num2str(round(rad2deg(val))), '°'];
        case 'cbEPs'
            valstr = [num2str(round(rad2deg(val))), '°'];
        case 'cbEPx'
            valstr = [num2str(val*guicp.array.ElPitch*1e3,2), 'mm'];
        case 'cbEPy'
            valstr = [num2str(val*guicp.array.ElPitch*1e3,2), 'mm'];
        case 'cbCpl'
            valstr = [num2str(round(val*100)), '%'];
    end
end