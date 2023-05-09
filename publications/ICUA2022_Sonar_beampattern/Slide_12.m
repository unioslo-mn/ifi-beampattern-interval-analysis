%% Presentation slide 12
close all; clear all;

%% Set simulation parameters (include update function)

array = biat.SensorArray(   'ElCount',          30,...
                            'ElPitchRatio',     0.5,...
                            'ElDiameterRatio',  0,...
                            'Curvature',        0,...
                            'SteeringAngle',    deg2rad(20),...
                            'TaperType',        'chebwin',...
                            'TaperParam',       40,...
                            'GainError',        3/100,...
                            'PhaseError',       deg2rad(3));  

bp_nom = biat.BeamPattern(array,'nominal','BeamResolutionDeg',0.5);
bp_pol = biat.BeamPattern(array,'polygonal','BeamResolutionDeg',0.5,...
                                            'PolygonTolerance',0.001);
P_nom = bp_nom.calculateBeamPattern;
P_pol = bp_pol.calculateBeamPattern;


%% Plot array and errors
fig1 = figure(1);clf
array.plot

%% Nominal
bp_pol.BeamIndex = 122;
P_bt = bp_pol.backtrackBeamPattern;
W = bp_pol.MaxErrorPattern;

%% Plot BPs

degs = rad2deg(bp_nom.BeamAngles);

fig2 = figure(2);clf
fig2.Position = [200 200 600 240]; 
set(gca,'DefaultLineLineWidth',2)

hold on; xlim([-90,90]); ylim([-50,2]); grid on; set(gcf, 'color', 'white');
xlabel('Angle (deg)')
ylabel('Power (dB)')
set(gca,'FontSize',9)
yticks([-40,-30,-20,-10, 0])
xticks([-90,-45,0,45,90])

plot(degs, db(P_pol(:,2))/2, 'r','DisplayName', 'Upper bound');
plot(degs, db(P_nom(:,1))/2, 'color',[0.1,0.75,0.1], 'DisplayName', 'Nominal');
plot(degs, db(P_pol(:,1) + eps)/2,'b', 'DisplayName', 'Lower bound');

plot(degs, db(P_bt)/2,'k-.','linewidth',2.5, 'DisplayName', 'Worst realization');

legend()

exportgraphics(gcf,'Backtracking_example-1.pdf')%,'Resolution',300)

%%

M = array.ElCount;
gI = array.GainInterval;
phI = array.PhaseInterval;


fig3 = figure(3);clf
fig3.Position = [200 200 600 240]; 

set(gcf, 'color', 'white');

subplot(1,2,1); hold on; ylim([0.95,1.05])

xlabel('Element')
ylabel('Amplitude, 1+g_c ( - )')
set(gca,'FontSize',12)

bar(abs(W),'FaceColor',[.3 .3 .3],'EdgeColor',[0 0 0],'displayname', 'Realization')
plot(1:M,[gI.Supremum],'r:','linewidth',5, 'displayname', 'Upper bound')
plot(1:M,[gI.Infimum],'b:','linewidth',5, 'displayname', 'Lower bound')
legend('Location','northwest')

subplot(1,2,2); hold on;
xlabel('Element')
ylabel('Phase, \phi_c (deg)')
set(gca,'FontSize',12)
bar(rad2deg(angle(W)),'FaceColor',[.3 .3 .3],'EdgeColor',[0 0 0],'displayname', 'Realization')
plot(1:M,rad2deg([phI.Supremum]),'r:','linewidth',5, 'displayname', 'upper bound')
plot(1:M,rad2deg([phI.Infimum]),'b:','linewidth',5, 'displayname', 'lower bound')

exportgraphics(gcf,'Backtracking_example-2.pdf')%,'Resolution',300)