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
                            'SteeringAngle',rand*90);%deg2rad(5));

bp_rec = biat.BeamPattern(array,'rectangular','BeamResolutionDeg',0.35/2);
bp_gon = biat.BeamPattern(array,'polygonal','BeamResolutionDeg',0.35/2, ...
                                            'PolygonTolerance',1e-6);
bp_arc = biat.BeamPattern(array,'polyarcular','BeamResolutionDeg',0.35/2);
bp_arx = biat.BeamPattern(array,'polyarx','BeamResolutionDeg',0.35/2);

%% Calculate beampattern and measure time

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

%% Optional: Still under debugging

arcTime = 0;
arcBeamPattern = zeros(bp_arc.BeamCount,2);

% % Calculate polyarx beampattern 
% tic
% for idxBeam = 363 : bp_arc.BeamCount  
%     bp_arc.BeamIndex = idxBeam;
%     arcBeamPattern(idxBeam,:) = bp_arc.ArrayInterval.Abs.Bounds.';
% end
% arcTime = toc;

%% Plot

%figure;clf
cla;hold on
bp_rec.plotBeamPattern('BeamPattern',recBeamPattern,'Color','g','LineWidth',2);
bp_rec.plotBeamPattern('BeamPattern',gonBeamPattern,'Color','r','LineWidth',2);
bp_rec.plotBeamPattern('BeamPattern',arxBeamPattern,'Color','c','LineWidth',2,'LineStyle','--');
bp_rec.plotBeamPattern('BeamPattern',arcBeamPattern,'Color','b','LineWidth',2,'LineStyle','--');

%% Report


sprintf(['Rectangular time: %0.1fs\n'...
         'Polygonal time: %0.1fs\n'...
         'Polyarx time: %0.1fs\n',...
         'Polyarc time: %0.1fs'], ...
         recTime,...
         gonTime,...
         arxTime,...
         arcTime)
