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

%% Calculate sum and measure time

% Add rectangular intervals
recElementInt = bp_rec.ElementIntervals;
tic
recArrayInt = sum(recElementInt);
recTime = toc;

% Add polygonal intervals
gonElementInt = bp_gon.ElementIntervals;
tic
gonArrayInt = sum(gonElementInt);
gonTime = toc;

% Add concanve polyarcular intervals
arcElementInt = bp_arc.ElementIntervals;
tic
arcArrayInt = sum(arcElementInt);
arcTime = toc;

% Add convex polyarcular intervals
arxElementInt = bp_arx.ElementIntervals;
tic
arxArrayInt = sum(arxElementInt);
arxTime = toc;



%% Plot

% figure;clf
cla;hold on; axis equal

recElementInt.plot('g-');
recArrayInt.plot('g-','linewidth',2);

gonElementInt.plot('r-');
gonArrayInt.plot('r-','linewidth',2);

arcElementInt.plot('b-');
arcArrayInt.plot('b-','linewidth',2);

arxElementInt.plot('c-');
arxArrayInt.plot('c-','linewidth',2);

%% Report

recArea = recArrayInt.Area;
gonArea = gonArrayInt.Area;
arxArea = arxArrayInt.Area;
arcArea = arcArrayInt.Area;

sprintf(['Rectangular area: %0.4f (tightness: %0.1f%%), Time: %0.1fms\n'...
         'Polygonal area: %0.4f (tightness: %0.1f%%), Time: %0.1fms\n'...
         'Polyarcular (convex) area: %0.4f (tightness: %0.1f%%), Time: %0.1fms\n',...
         'Polyarcular (concave) area: %0.4f (tightness: %0.1f%%), Time: %0.1fms'], ...
         recArea, arcArea / recArea * 100, recTime*1e3,...
         gonArea, arcArea / gonArea * 100, gonTime*1e3,...
         arxArea, arcArea / arxArea * 100, arxTime*1e3,...
         arcArea, arcArea / arcArea * 100, arcTime*1e3)


