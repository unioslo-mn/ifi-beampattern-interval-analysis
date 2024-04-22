% Reset controllers to default
function resetControllers(h)
    % Get parameter struct
    guicp = evalin('base','gui');

    % Get controller handlers
    ctrl = h.Parent.Children;
    
    % Update controllers
    for idx = 1:length(ctrl)
        if ~isempty(ctrl(idx).Callback)
            name = func2str(ctrl(idx).Callback);
            switch name
                case 'cbNel'
                    ctrl(idx).Value = guicp.array.ElCount;
                    guicp.val.Nel.String = formatValueString(name,...
                                                        ctrl(idx).Value);
                case 'cbRcv'
                    ctrl(idx).Value = guicp.array.Curvature;
                    guicp.val.Rcv.String = formatValueString(name,...
                                                        ctrl(idx).Value);
                case 'cbEld'
                    ctrl(idx).Value = guicp.array.ElPitchRatio;
                    guicp.val.Eld.String = formatValueString(name,...
                                                        ctrl(idx).Value);
                case 'cbEls'
                    ctrl(idx).Value = guicp.array.ElDiameterRatio;
                    guicp.val.Els.String = formatValueString(name,...
                                                        ctrl(idx).Value);
                case 'cbThi'
                    ctrl(idx).Value = sin(((...
                                        guicp.beampattern(1).BeamIndex-1) / ...
                                       (guicp.beampattern(1).BeamCount-1)*...
                                                              2-1)*pi/2);
                    guicp.val.Thi.String = formatValueString(name,...
                                                        ctrl(idx).Value);
                case 'cbSte'
                    ctrl(idx).Value = guicp.array.SteeringAngle/pi*2;
                    guicp.val.Ste.String = formatValueString(name,...
                                                        ctrl(idx).Value);
                case 'cbWgh'
                    ctrl(idx).Value = guicp.array.TaperParam;
                    guicp.val.Wgh.String = formatValueString(name,...
                                                        ctrl(idx).Value);
                case 'cbEGn'
                    ctrl(idx).Value = guicp.array.GainError;
                    
                    guicp.val.EGn.String = formatValueString(name,...
                                                        ctrl(idx).Value);
                case 'cbEPh'
                    ctrl(idx).Value = guicp.array.PhaseError;
                    guicp.val.EPh.String = formatValueString(name,...
                                                        ctrl(idx).Value);
                case 'cbEPs'
                    ctrl(idx).Value = guicp.array.OrientError;
                    guicp.val.EPs.String = formatValueString(name,...
                                                        ctrl(idx).Value);
                case 'cbEPx'
                    ctrl(idx).Value = guicp.array.PosXError;
                    guicp.val.EPx.String = formatValueString(name,...
                                                        ctrl(idx).Value);
                case 'cbEPy'
                    ctrl(idx).Value = guicp.array.PosYError;
                    guicp.val.EPy.String = formatValueString(name,...
                                                        ctrl(idx).Value);
                case 'cbCpl'
                    ctrl(idx).Value = guicp.array.CouplingCoeff;
                    guicp.val.Cpl.String = formatValueString(name,...
                                                        ctrl(idx).Value);
            end
        end
    end
    
    % Get tickbox handlers
    tbx = h.Parent.Parent.Children(1).Children;
    
    % Update tickboxes
    for idx = 1:length(tbx)
        if ~isempty(tbx(idx).Callback)
            name = func2str(tbx(idx).Callback);
            switch name
                case 'cbTbN'
                    tbx(idx).Value = guicp.plt.sel(1);
                case 'cbTbR'
                    tbx(idx).Value = guicp.plt.sel(2);
                case 'cbTbC'
                    tbx(idx).Value = guicp.plt.sel(3);
                case 'cbTbM'
                    tbx(idx).Value = guicp.plt.sel(4);
                case 'cbTbA'
                    tbx(idx).Value = guicp.plt.sel(5);
                case 'cbTbB'
                    tbx(idx).Value = guicp.plt.sel(6);
            end
        end
    end
end

