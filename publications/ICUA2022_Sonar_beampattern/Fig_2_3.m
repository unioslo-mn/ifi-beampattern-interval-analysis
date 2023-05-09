% Start dashboard
BeamPatternDashboard

% Set array parameters
gui.ctrl.Nel.Value = 11;            updateArray(gui.ctrl.Nel);  % Number of elements
gui.ctrl.Rcv.Value = 0.13;          updateArray(gui.ctrl.Rcv);  % Curvature
gui.ctrl.Eld.Value = 0.5;           updateArray(gui.ctrl.Eld);  % Element distance ratio
gui.ctrl.Els.Value = 0.8;           updateArray(gui.ctrl.Els);  % Element diameter ratio
gui.ctrl.Ste.Value = 15/90;         updateArray(gui.ctrl.Ste);  % Steering angle (pi/2 ratio)
gui.ctrl.Wgh.Value = 1;             updateArray(gui.ctrl.Wgh);  % Tapering Kaiser beta
gui.ctrl.EGn.Value = 0.03;          updateArray(gui.ctrl.EGn);  % Gain error
gui.ctrl.EPh.Value = deg2rad(3);    updateArray(gui.ctrl.EPh);  % Phase error
gui.ctrl.EPs.Value = deg2rad(3);    updateArray(gui.ctrl.EPs);  % Orientation error
gui.ctrl.EPx.Value = 1/100;         updateArray(gui.ctrl.EPx);  % Pos-X error (pitch-ratio)
gui.ctrl.EPy.Value = 1/100;         updateArray(gui.ctrl.EPy);  % Pos-Y error (pitch-ratio)
gui.ctrl.Cpl.Value = 0;             updateArray(gui.ctrl.Cpl);  % Coupling
gui.ctrl.Thi.Value = deg2rad(0);    updateArray(gui.ctrl.Thi);  % Incidence angle index
gui.ctrl.DRn.Value = -60;           updateArray(gui.ctrl.DRn);  % Dynamic range floor

% Update plots
gui.btn.Update.Value = 1;
updatePlots(gui.btn.Update);
updateBeamPattern(gui.ctrl.Thi);