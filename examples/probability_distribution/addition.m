clear
close all
%% 

% Array parameters
M = 10;

% Error parameters
errGain = 0.1;
errPha = deg2rad(5);
errPsi = deg2rad(1);

% Beamforming parameters
incAngle = deg2rad(-11);
AngResDeg = 1;
PolTol = 5e-4;
PolIncl = true;

% Generate array and beampattern objects
% Initialize array
array = biat.SensorArray(       'ElCount',          M,...
                                'ElPitchRatio',     0.5,...
                                'ElDiameterRatio',  0,...
                                'Curvature',        0,...
                                'TaperType',        'none',...
                                'GainError',        errGain,...
                                'PhaseError',       errPha,...
                                'PosXError',        0,...   
                                'PosYError',        0,...
                                'OrientError',      errPsi,...
                                'CouplingCoeff',    0,...
                                'SteeringAngle',    0);

% Calculate nominal beampattern with tapering
bp_rec = biat.BeamPattern(array,'rectangular');
bp_cir = biat.BeamPattern(array,'circular');
bp_pgn = biat.BeamPattern(array,'polygonal','BeamResolutionDeg',AngResDeg,...
                                        'PolygonTolerance',PolTol);


%% Addition of rectangular intervals in various ways

bp_rec.BeamIndex = 72;

A = bp_rec.ElementIntervals(1);
B = bp_rec.ElementIntervals(2);
C = bp_rec.ElementIntervals(3);
D = bp_rec.ElementIntervals(4);

A.ProbaGrid = ciat.ProbaGrid(A, "normal", 'nx', 100, 'ny', 100);
B.ProbaGrid = ciat.ProbaGrid(B, "normal", 'nx', 100, 'ny', 100);
C.ProbaGrid = ciat.ProbaGrid(C, "normal", 'nx', 100, 'ny', 100);
D.ProbaGrid = ciat.ProbaGrid(D, "normal", 'nx', 100, 'ny', 100);

AC = [A;C];
BD = [B;D];

ACpBD = AC + BD;
ABpCDsum = sum(ACpBD);

figure;hold on
AC(1).plot;
AC(2).plot;
BD(1).plot;
BD(2).plot;
ACpBD(1).plot;
ACpBD(2).plot;
ABpCDsum.plot

%% Addition of circular intervals

bp_cir.BeamIndex = 72;

A = bp_cir.ElementIntervals(1);
B = bp_cir.ElementIntervals(2);
C = bp_cir.ElementIntervals(3);
D = bp_cir.ElementIntervals(4);

A.ProbaGrid = ciat.ProbaGrid(A, "normal", 'nx', 100, 'ny', 100);
B.ProbaGrid = ciat.ProbaGrid(B, "normal", 'nx', 100, 'ny', 100);
C.ProbaGrid = ciat.ProbaGrid(C, "normal", 'nx', 100, 'ny', 100);
D.ProbaGrid = ciat.ProbaGrid(D, "normal", 'nx', 100, 'ny', 100);

AC = [A;C];
BD = [B;D];

ACpBD = AC + BD;
ABpCDsum = sum(ACpBD);

figure;hold on
AC(1).plot;
AC(2).plot;
BD(1).plot;
BD(2).plot;
ACpBD(1).plot;
ACpBD(2).plot;
ABpCDsum.plot

%% Addition of polygonal intervals

bp_pgn.BeamIndex = 72;

A = bp_pgn.ElementIntervals(1);
B = bp_pgn.ElementIntervals(2);
C = bp_pgn.ElementIntervals(3);
D = bp_pgn.ElementIntervals(4);

A.ProbaGrid = ciat.ProbaGrid(A, "polarnormal", 'nx', 100, 'ny', 100);
B.ProbaGrid = ciat.ProbaGrid(B, "polarnormal", 'nx', 100, 'ny', 100);
C.ProbaGrid = ciat.ProbaGrid(C, "polarnormal", 'nx', 100, 'ny', 100);
D.ProbaGrid = ciat.ProbaGrid(D, "polarnormal", 'nx', 100, 'ny', 100);

AC = [A;C];
BD = [B;D];

ACpBD = AC + BD;
ABpCDsum = sum(ACpBD);

figure;hold on
AC(1).plot;
AC(2).plot;
BD(1).plot;
BD(2).plot;
ACpBD(1).plot;
ACpBD(2).plot;
ABpCDsum.plot

