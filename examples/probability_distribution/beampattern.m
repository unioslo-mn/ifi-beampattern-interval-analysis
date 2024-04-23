clear
close all
%% 

% Array parameters
M = 10;

% Error parameters
errGain = 0.1;
errPha = deg2rad(1);
errPsi = deg2rad(0);

% Beamforming parameters
incAngle = deg2rad(-11);
AngResDeg = 1/2;
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
bp = biat.BeamPattern(array,'polygonal','BeamResolutionDeg',AngResDeg,...
                                        'PolygonTolerance',PolTol);
bp.plotBeamPattern
%% Show the complex probabilistic interval sum for a chosen steering angle

bp.BeamIndex = 72;
elInt = bp.ElementIntervals;
elInt = elInt.setProbaGrid("polarnormal", 'nx', 100, 'ny', 100);
arInt = sum(elInt);
ppInt = abs2(arInt.ProbaGrid);

% Initialize figure
figure;clf;
subplot(2,1,1);hold on
elInt.plot;
arInt.plot;
subplot(2,1,2);hold on
plot(db(ppInt.x)',ppInt.Pdf)

%% Plot probabilistic beampattern

nx = 1e2;
ny = 1e2;
pCnt = 1e4;
dinRange = [-60,0];


beamCnt = bp.BeamCount;
pAxis = linspace(2*dinRange(1),2*dinRange(2),pCnt);
ppPdf = zeros(pCnt,beamCnt);

for beamIdx = 1:beamCnt
    % Set steering
    bp.BeamIndex = beamIdx;

    % Calculate element interval sum
    elInt = bp.ElementIntervals;
    elInt = elInt.setProbaGrid("polarnormal", 'nx', nx, 'ny', ny);
    arInt = sum(elInt);

    % Extract power pattern PDF
    powerInt = abs2(arInt.ProbaGrid);
    powerPdf = interp1(db(powerInt.x),powerInt.Pdf,pAxis);

    % Replace NaN values with zeros
    powerPdf(isnan(powerPdf)) = 0;

    % Assign to image
    ppPdf(:,beamIdx) = flip(powerPdf);
    fprintf('%0.0f/%0.0f\n',beamIdx,beamCnt)
end
%% Plot results

% Calculate beampattern
bpVal = bp.calculateBeamPattern;
bpVal(bpVal==0) = eps;

% Initialize figure
close all
figure;clf;hold on

% Plot
imagesc(bp.BeamAngles,pAxis/2,flipud(ppPdf))
plot(bp.BeamAngles,db(bpVal,'power'))
ylim(dinRange)
xlim([-pi/2,pi/2])
colormap(turbo)
set(gca,'ColorScale','log')
colorbar