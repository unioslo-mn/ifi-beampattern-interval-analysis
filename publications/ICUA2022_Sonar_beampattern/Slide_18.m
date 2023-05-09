%% Presentation slide 18
close all; clear all;

%% Set simulation parameters (include update function)

array = biat.SensorArray(   'ElCount',          33,...
                            'ElPitchRatio',     0.5*1.015,...
                            'ElDiameterRatio',  0.95,...
                            'Curvature',        0,...
                            'SteeringAngle',    deg2rad(20),...
                            'TaperType',        'none',...
                            'GainError',        0,...
                            'PhaseError',       0,...
                            'PosXError',        5e-4,...
                            'PosYError',        5e-4);  
                        
bp_nom = biat.BeamPattern(array,'nominal','BeamResolutionDeg',0.5);
bp_pol = biat.BeamPattern(array,'polygonal','BeamResolutionDeg',0.5,...
                                            'PolygonTolerance',1e-6);
P_nom = bp_nom.calculateBeamPattern;
P_pol = bp_pol.calculateBeamPattern;

degs = rad2deg(bp_nom.BeamAngles);

%% Plot array and errors
fig1 = figure(1);
array.plot;

%% Plot BPs

figure('Position',[200 200 350 240]); 
set(gca,'DefaultLineLineWidth',2)

hold on; xlim([-90,90]); ylim([-40,2]); grid on; set(gcf, 'color', 'white');
xlabel('Angle (deg)')
ylabel('Power (dB)')
set(gca,'FontSize',9)
yticks([-40,-30,-20,-10, 0])
xticks([-90,-45,0,45,90])

plot(degs, db(P_pol(:,2))/2, 'r','DisplayName', 'Upper bound');
plot(degs, db(P_nom(:,1))/2, 'color',[0.1,0.75,0.1], 'DisplayName', 'Nominal');
plot(degs, db(P_pol(:,1) + eps)/2,'b', 'DisplayName', 'Lower bound');

legend()
