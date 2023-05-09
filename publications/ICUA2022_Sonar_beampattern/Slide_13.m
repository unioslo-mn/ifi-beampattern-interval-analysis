%% Presentation slide 13
close all; clear

%% Set simulation parameters and calculate beampattern element errors

array = biat.SensorArray(   'ElCount',          30,...
                            'ElPitchRatio',     0.5,...
                            'ElDiameterRatio',  0,...
                            'Curvature',        0,...
                            'SteeringAngle',    deg2rad(30),...
                            'TaperType',        'none',...
                            'GainError',        1/100,...
                            'PhaseError',       deg2rad(2),...
                            'PosXError',        5e-4,...
                            'PosYError',        5e-4);  
                        
bp_nom = biat.BeamPattern(array,'nominal','BeamResolutionDeg',0.35/2);
bp_pol = biat.BeamPattern(array,'polygonal','BeamResolutionDeg',0.35/2,...
                                            'PolygonTolerance',1e-6);
P_nom = bp_nom.calculateBeamPattern;
P_pol = bp_pol.calculateBeamPattern;


%% Calculate beampattern with block errors

array.GainError = 0;
array.PhaseError = 0;
array.PosXError = 0;
array.PosYError = 0;

block = biat.SensorArray(   'ElCount',          3,...
                            'ElPitchRatio',     0.5*10,...
                            'ElDiameterRatio',  0,...
                            'Curvature',        0,...
                            'SteeringAngle',    deg2rad(30),...
                            'PosXError',        0e-4,...
                            'PosYError',        0e-4);


bp_blk = biat.BeamPattern(array,'polygonal','Block',block,...
                                            'BlockType','split',...
                                            'BeamResolutionDeg',0.35/2,...
                                            'PolygonTolerance',1e-6);

P_blk = bp_blk.calculateBeamPattern;

%% Plot beampatterns

degs = rad2deg(bp_nom.BeamAngles);

figure('Position',[200 200 500 240]); 
set(gca,'DefaultLineLineWidth',2)
hold on; xlim([-60,30]); ylim([-40,0]); grid on; set(gcf, 'color', 'white');
xlabel('Angle (deg)')
ylabel('Power (dB)')
set(gca,'FontSize',9)
yticks([-60:10: 0])
xticks([-90:15:90])


plot(degs, db(P_pol(:,2))/2, 'm:','DisplayName', 'P_U (element)');
plot(degs, db(P_blk(:,2))/2, 'r-','linewidth', 3,'DisplayName', 'P_U (blocks)');
plot(degs, db(P_nom(:,1))/2, 'color',[0.1,0.75,0.1],'DisplayName', 'nominal');
plot(degs, db(P_blk(:,1))/2, 'b-','linewidth', 3,'DisplayName', 'P_L (blocks)');
plot(degs, db(P_pol(:,1))/2, 'c:','DisplayName', 'P_L (element)');

legend('location','northwest')

set(gca, 'fontsize', 14)
