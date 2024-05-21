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

%% Plot

%figure;clf
cla;hold on
bp_rec.plotBeamPattern('BeamPattern',recBeamPattern,'Color','g','LineWidth',2);
bp_rec.plotBeamPattern('BeamPattern',gonBeamPattern,'Color','r','LineWidth',2);
bp_rec.plotBeamPattern('BeamPattern',arxBeamPattern,'Color','c','LineWidth',2,'LineStyle','--');

%% Report


sprintf(['Rectangular time: %0.1fms\n'...
         'Polygonal time: %0.1fms\n'...
         'Polyarx time: %0.1fms\n'], ...
         recTime*1e3,...
         gonTime*1e3,...
         arxTime*1e3);
