clear 
% close all

%%
% Define array
array = biat.SensorArray(   'ElCount',5,...
                            'ElDiameterRatio',0,...
                            'Curvature',0.2,...
                            'TaperType','chebwin',...
                            'TaperParam',20,...
                            'GainError',5/100,...
                            'PhaseError',deg2rad(4),...
                            'SteeringAngle',deg2rad(5));

bp_cir = biat.BeamPattern(array,'circular','BeamResolutionDeg',0.35/2);
bp_rec = biat.BeamPattern(array,'rectangular','BeamResolutionDeg',0.35/2);
bp_gon = biat.BeamPattern(array,'polygonal','BeamResolutionDeg',0.35/2, ...
                                            'PolygonTolerance',1e-6);
bp_arc = biat.BeamPattern(array,'polyarcular','BeamResolutionDeg',0.35/2);
bp_arx = biat.BeamPattern(array,'polyarx','BeamResolutionDeg',0.35/2);

%% Calculate beampattern and measure time

% Calculate rectangular beampattern
tic
cirBeamPattern = bp_cir.calculateBeamPattern;
cirTime = toc;

% Calculate rectangular beampattern
tic
recBeamPattern = bp_rec.calculateBeamPattern;
recTime = toc;

% Calculate polygonal beampattern
tic
gonBeamPattern = bp_gon.calculateBeamPattern;
gonTime = toc;

% Calculate polyarx beampattern
tic
arxBeamPattern = bp_arx.calculateBeamPattern;
arxTime = toc;

% Placeholder for polyarcular beampattern calculation
arcTime = 0;
arcBeamPattern = zeros(bp_arc.BeamCount,2);

%% Optional: This takes about two minutes

% Calculate polyarx beampattern 
tic
arcBeamPattern = bp_arc.calculateBeamPattern;
arcTime = toc;

%% Show an example for the element intervals

% Set beam angle
beamAngle = asin(-0.6);
[~,beamIndex] = min(abs(bp_rec.BeamAngles - beamAngle));
bp_cir.BeamIndex = beamIndex;
bp_rec.BeamIndex = beamIndex;
bp_gon.BeamIndex = beamIndex;
bp_arx.BeamIndex = beamIndex;
bp_arc.BeamIndex = beamIndex;

% Calculate array interval
cI = bp_cir.ArrayInterval;
rI = bp_rec.ArrayInterval;
gI = bp_gon.ArrayInterval;
xI = bp_arx.ArrayInterval;
aI = bp_arc.ArrayInterval;

% Cast intervals to polar
cpI = ciat.PolarInterval(cI);
rpI = ciat.PolarInterval(rI);
gpI = ciat.PolarInterval(gI);
xpI = ciat.PolarInterval(xI);
apI = ciat.PolarInterval(aI);

%% Plot

%figure;clf
subplot(2,2,1);cla;hold on
plot(sin(bp_rec.BeamAngles),db(cirBeamPattern,'power'),'r-','LineWidth',2','DisplayName','cir');
plot(sin(bp_rec.BeamAngles),db(recBeamPattern,'power'),'g-','LineWidth',2','DisplayName','rec');
plot(sin(bp_rec.BeamAngles),db(gonBeamPattern,'power'),'k-','LineWidth',2','DisplayName','gon');
plot(sin(bp_rec.BeamAngles),db(arxBeamPattern,'power'),'c--','LineWidth',2','DisplayName','arx');
plot(sin(bp_rec.BeamAngles),db(arcBeamPattern,'power'),'b:','LineWidth',2','DisplayName','arc');
yLim = ylim;
plot(sin(beamAngle)*[1;1] , ylim.' ,'k:','DisplayName','rec')
legend(legendUnq)

subplot(2,2,3);cla;hold on;title('Beampattern supremum difference')
plot(sin(bp_rec.BeamAngles),db(gonBeamPattern(:,2),'power')-db(arcBeamPattern(:,2),'power'));
plot(sin(bp_rec.BeamAngles),db(arxBeamPattern(:,2),'power')-db(arcBeamPattern(:,2),'power'));
plot(sin(beamAngle)*[1;1] , ylim.' ,'k:')
yLim = ylim;
xlabel('sin(angle)')
ylabel('Amplitude [dB]')
legend('g-a','x-a')

% Plot example of element intervals
subplot(2,2,[2,4]);cla;hold on;axis equal;title('Element intervals')
cI.plot('r','LineWidth',2');
rI.plot('g','LineWidth',2');
gI.plot('k','LineWidth',2');
xI.plot('c--','LineWidth',2');
aI.plot('b--','LineWidth',2');
cpI.plot('r:','LineWidth',1');
rpI.plot('g:','LineWidth',1');
gpI.plot('k:','LineWidth',1');
xpI.plot('c:','LineWidth',1');
apI.plot('b:','LineWidth',1');

%% Report


sprintf(['Circular time: %0.1fs\n'...
         'Rectangular time: %0.1fs\n'...
         'Polygonal time: %0.1fs\n'...
         'Polyarx time: %0.1fs\n',...
         'Polyarc time: %0.1fs'], ...
         cirTime,...
         recTime,...
         gonTime,...
         arxTime,...
         arcTime)
