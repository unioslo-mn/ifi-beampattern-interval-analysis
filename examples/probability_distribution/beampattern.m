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
bp_pgn = biat.BeamPattern(array,'polygonal','BeamResolutionDeg',AngResDeg,...
                                        'PolygonTolerance',PolTol);

%% Show the complex probabilistic interval sum for a chosen steering angle

bp_pgn.BeamIndex = 72;
elInt = bp_pgn.ElementIntervals;

% Initialize figure
figure;clf;hold on

for m = 1:array.ElCount
    elInt(m).ProbaGrid = ciat.ProbaGrid(elInt(m), ...
                                "polarnormal", 'nx', 100, 'ny', 100);
    elInt(m).plot;
end
arInt = sum(elInt);
arInt.plot


%% Plot probabilistic beampattern

nx = 100;
ny = 100;

N = bp_pgn.BeamCount;
M = array.ElCount;

bp_probPdf = zeros(nx,N); 
bp_probX = zeros(nx,N); 

for n = 1:N
    bp_pgn.BeamIndex = n;
    elInt = bp_pgn.ElementIntervals;
    for m = 1:M
        elInt(m).ProbaGrid = ciat.ProbaGrid(elInt(m), ...
                                    "polarnormal", 'nx', nx, 'ny', ny);
    end
    arInt = sum(elInt);
    bp_prob = abs2(arInt.ProbaGrid);
    bp_probPdf(:,n) = bp_prob.Pdf; 
    bp_probX(:,n) = bp_prob.x; 
end

% Initialize figure
figure;clf;hold on

% Plot
bp_pgn.plotBeamPattern

%% 

bp_pgnInt = zeros(1e3,N);
xInt = linspace(-60,0,1e3);

for n = 1:N
    bp_pgnInt(:,n) = interp1( db(bp_probX(:,n)) , bp_probPdf(:,n) , xInt );
    bp_pgnInt(:,n) = bp_pgnInt(:,n) / sum(bp_pgnInt(~isnan(bp_pgnInt(:,n)),n) );
end


