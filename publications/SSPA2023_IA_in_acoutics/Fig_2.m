%% Beampattern tolerance analysis

% This script invokes the Beampattern Dashboard and sets its control
% values to produce the beampattern bounds of a 3 element array.
% The top right plot shows the acoustic array in question.
% The bottom left plot shows the power pattern bounds using 
% rectangular (cyan) circular (red) and polygonal (blue) complex interval 
% representations. The bottom right plot shows the corresponding 
% element intervals and their sum on the complex plane.

% Requirement: examples/Beampattern_dashboard folder in the path

% Start dashboard
BeamPatternDashboard

% Set array parameters
gui.ctrl.Nel.Value = 3;            updateArray(gui.ctrl.Nel);  % Number of elements
gui.ctrl.Rcv.Value = 0.16;          updateArray(gui.ctrl.Rcv);  % Curvature
gui.ctrl.Eld.Value = 0.37;          updateArray(gui.ctrl.Eld);  % Element distance ratio
gui.ctrl.Els.Value = 0;             updateArray(gui.ctrl.Els);  % Element diameter ratio
gui.ctrl.Ste.Value = 0;             updateArray(gui.ctrl.Ste);  % Steering angle
gui.ctrl.Wgh.Value = 0;             updateArray(gui.ctrl.Wgh);  % Tapering Kaiser beta
gui.ctrl.EGn.Value = 0.11;          updateArray(gui.ctrl.EGn);  % Gain error
gui.ctrl.EPh.Value = deg2rad(12);   updateArray(gui.ctrl.EPh);  % Phase error
gui.ctrl.EPs.Value = deg2rad(2);    updateArray(gui.ctrl.EPs);  % Orientation error
gui.ctrl.EPx.Value = 0;             updateArray(gui.ctrl.EPx);  % Pos-X error
gui.ctrl.EPy.Value = 0;             updateArray(gui.ctrl.EPy);  % Pos-Y error
gui.ctrl.Cpl.Value = 0;             updateArray(gui.ctrl.Cpl);  % Coupling
gui.ctrl.Thi.Value = deg2rad(17);   updateArray(gui.ctrl.Thi);  % Incidence angle index
gui.ctrl.DRn.Value = -60;           updateArray(gui.ctrl.DRn);  % Dynamic range floor

% Update plots
gui.btn.Update.Value = 1;
updatePlots(gui.btn.Update);
updateBeamPattern(gui.ctrl.Thi);