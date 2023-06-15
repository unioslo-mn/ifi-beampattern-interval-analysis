%% Initilize workspace
clear 
close all

%% Set simulation parameters (include update function)
% Note: there is a certain order in which things must be set as of now
% Initialize simularion parameters
array = biat.SensorArray(   'ElCount',5,...
                            'ElDiameterRatio',0,...
                            'Curvature',0.2,...
                            'TaperType','chebwin',...
                            'TaperParam',20,...
                            'GainError',5/100,...
                            'PhaseError',deg2rad(4),...
                            'SteeringAngle',deg2rad(5));

bp_nom = biat.BeamPattern(array,'nominal','BeamResolutionDeg',0.35/2);
bp_pol = biat.BeamPattern(array,'polygonal','BeamResolutionDeg',0.35/2,...
                                            'PolygonTolerance',1e-6);

% Set beam index
beamAngle = deg2rad(50.5);
[~,bp_nom.BeamIndex] = min(abs(bp_nom.BeamAngles - beamAngle));
bp_pol.BeamIndex = bp_nom.BeamIndex;

% Make phase intervals non-uniform
M = array.ElCount;
phErrW = (1 + 0.5 * ( ( (1:M) - (M+1)/2) ./ ((M-1) / 2)) .^2)';
array.PhaseError = array.PhaseError .* phErrW;

% declare some parameters
thetas = bp_pol.BeamAngles;
degs = rad2deg(thetas);
theta_s = array.SteeringAngle;
k_s = array.SteeringVector;

% backtracking angle
theta_probe = deg2rad(50.5); % probe 50 degrees (almost on sidelobe in plot)
th0 = rad2deg(theta_probe);
[~, theta_idx] = min(abs(thetas - theta_probe));

% compute expected beampattern
c_0 = array.SoundSpeed;
f_0 = array.CenterFrequency;
c_t = c_0;
w_apod = array.TaperWeights;
T_se = sum(w_apod.^2);
var_g = 1/12 * (array.GainInterval(1).Width)^2;
var_phase = 1/12 * ( mean([array.PhaseInterval.Width]) )^2;
var_pos_lambda = (2*pi / (c_0/f_0))^2 * 1/12 *...
                 (array.PosXInterval(1).Width)^2;

B1 = exp(-(var_phase + var_pos_lambda));
B2 = T_se * ( 1 + var_g - exp(-(var_phase + var_pos_lambda)));

% some parameters
lambda = c_0 / f_0;
elRad = array.ElDiameter/2;
taperAngs = array.TaperAngles;

%% Backtrack

P_nom = bp_nom.calculateBeamPattern;
shapes = bp_pol.ElementIntervals;
final_potato = bp_pol.ArrayInterval;
extreme_points_min = bp_pol.MinPowerPoints;
extreme_points_max = bp_pol.MaxPowerPoints;
element_error = bp_pol.MaxErrorPattern;
element_error_min = bp_pol.MinErrorPattern;
P_GENPOLY = bp_pol.calculateBeamPattern;
P_backtrack = bp_pol.backtrackBeamPattern('getMax',1);
P_backtrack_min = bp_pol.backtrackBeamPattern('getMax',0);

%% Backtracked errors vs. error bounds - Fig. 8

figure('Position',[200 200 600 240]); clf
set(gcf, 'color', 'white');
subplot(1,2,1); hold on; 

xlabel('Element {\itm}')
ylabel('Amplitude error (%)')
set(gca,'FontSize',14)

bar(100*(abs(element_error)-1),'FaceColor',[.8 .3 .3],'EdgeColor',[.8 .3 .3],'displayname', 'Realization, BP upper')
bar(100*(abs(element_error_min)-1),0.33,'FaceColor',[.1 .1 .7],'EdgeColor',[.1 .1 .7],'displayname', 'Realization, BP lower')
plot((1:M),100*([array.GainInterval.Supremum]-1),'r:o','linewidth',5, 'markersize', 10, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'r', 'displayname', 'Upper error bound')
plot((1:M),100*([array.GainInterval.Infimum]-1),'b:o','linewidth',5, 'markersize', 10, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'b', 'displayname', 'Lower error bound')
legend('Location','southwest'); ylim(100*([0.85,1.07]-1)); xlim([0.5, M+0.5]);
legend boxoff

subplot(1,2,2); hold on;
xlabel('Element {\itm}')
ylabel('Phase error (deg)')
set(gca,'FontSize',14)
bar(rad2deg(angle(element_error)),'FaceColor',[.8 .3 .3],'EdgeColor',[.8 .3 .3],'displayname', 'Realization')
bar(rad2deg(angle(element_error_min)),0.33,'FaceColor',[.1 .1 .7],'EdgeColor',[.1 .1 .7],'displayname', 'Realization, lower')
plot((1:M),rad2deg([array.PhaseInterval.Supremum]),'r:o','linewidth',5, 'markersize', 10, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'r', 'displayname', 'Upper error bound')
plot((1:M),rad2deg([array.PhaseInterval.Infimum]),'b:o','linewidth',5, 'markersize', 10, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'b', 'displayname', 'Lower error bound')
xlim([0.5, M+0.5]);  ylim([-7,7]); 

ax = gcf;
%exportgraphics(ax,'Fig8_BP_1a_backtracked_errors.pdf')%,'Resolution',300) 



%% Plot BPs - Fig. 6b
k_abs = array.WaveNumber;
theta_i = bp_pol.IncidenceAngle;
k = k_abs * ([sin(theta_i); cos(theta_i)] );
dk = (k - k_s);
phase = dk(1,:) .* [array.ElPosX] + dk(2,:) .* [array.ElPosY];
dPsi = bp_pol.IncidenceAngle - bp_pol.Arrays.OrientInterval;
dir = array.getElDirectivity(dPsi);
vecs =([dir.Midpoint]' .* w_apod) .* exp(1j*phase);


% Comment: in the article the origin of the array was defined differently,
% therefore the resulted complex intervals are rotated, it has no influence
% on the beampattern interval or its interpretation

figure(2); clf
set(gcf,'Position',[100 200 400 400])
set(gca, 'fontsize', 17)
hold on; axis equal; 
xlim([-0.2272 0.2846])
ylim([-0.2497 0.2622])
grid on; set(gcf, 'color', 'white');
xlabel('Real axis')
ylabel('Imaginary axis')
yticks([-0.2,-0.1,0,0.1,0.2])
xticks([-0.2,-0.1,0,0.1,0.2])

legend('AutoUpdate','off'); 
pgon = polyshape(real(final_potato.Points), imag(final_potato.Points));
plot(pgon,'FaceColor',[0.4,0.4,0.4],'FaceAlpha',0.3)

for idx=1:length(vecs)
    pgon = polyshape(real(shapes(idx).Points), imag(shapes(idx).Points));
    plot(pgon,'FaceColor',[0.4,0.4,0.4],'FaceAlpha',0.1, 'Displayname', 'Complex intervals')

    if idx == 4; legend('Location','northwest','AutoUpdate','on'); end
end

legend('Location','northwest','AutoUpdate','off');
plot(real(sum(vecs)), imag(sum(vecs)), 'ko','markersize',13,'MarkerFaceColor',[0.1 0.8 0.1]);

for idx=1:length(vecs)
    [x,y] = array.getArrowLine([0,0],[(real(vecs(idx))), (imag(vecs(idx)))], 0.02);
    plot(x, y,  'linewidth',3, 'color','k', 'displayname', 'Element response');
    if idx == 4; legend('Location','northwest','AutoUpdate','on'); end
end

legend('Location','northwest','AutoUpdate','off');
[xSum,ySum] = array.getArrowLine([0,0],[sum(real(vecs)), sum(imag(vecs))], 0.02);
plot(xSum, ySum,  'linewidth',4, 'color','k', 'HandleVisibility','off');
legend('Location','northwest','AutoUpdate','on');
plot(xSum, ySum,  'linewidth',3, 'color',[0.1,0.8,0.1], 'displayname', 'Array response');

plot((extreme_points_min), 'kv','markersize',10,'MarkerFaceColor',[0.3 0.3 1], 'linewidth', 1,'Displayname', 'z_c lower');
plot((extreme_points_max), 'k^','markersize',10,'MarkerFaceColor',[1   0   0], 'linewidth', 1,'Displayname', 'z_c upper');

plot(sum(extreme_points_min), 'kv','markersize',13,'MarkerFaceColor',[0.3 0.3 1]*0.5, 'linewidth', 1, 'HandleVisibility','off');
plot(sum(extreme_points_max), 'k^','markersize',13,'MarkerFaceColor',[1   0   0]*0.5, 'linewidth', 1, 'HandleVisibility','off');

legend('Location','northeast')
legend boxoff
ax = gcf;
%exportgraphics(ax,'Fig6b_BP_1a_backtracked.pdf')



%% Fig 6a
figure(3);clf

set(gcf,'Position',[600 200 800 400])
set(gca, 'fontsize', 18)
hold on; xlim([-90-0.1,90]); ylim([-37,2]); grid on; set(gcf, 'color', 'white'); set(gca, 'color', [.85 .85 .85]);
xlabel('\theta (deg)')
ylabel('Power (dB)')
yticks([-40 -30,-20,-10, 0])
xticks([-90,-45,0,45,90])

pgon = polyshape([degs; flip(degs)], [db(P_GENPOLY(:,1)+ eps)/2 ; flip(db(P_GENPOLY(:,2))/2 + eps)]);
plot(pgon, 'linewidth', 1, 'edgealpha', 0, 'edgecolor',[0 0 0], 'FaceColor',[1 1 1],'FaceAlpha',1, 'handlevisibility', 'off')

degs = rad2deg(thetas);
p = plot(degs, db(P_backtrack,'power'), 'linewidth',4,'linestyle','-.', 'color',[1,0,0], 'Displayname','Backtrack upper');
p = plot(degs, db(P_nom(:,1),'power'), 'linewidth',4, 'color',[0.1,0.75,0.1], 'Displayname','Nominal');
p = plot(degs, db(P_backtrack_min,'power'), 'linewidth',4,'linestyle',':', 'color',[0 ,0 ,1], 'Displayname','Backtrack lower');

plot([th0,th0],[-50,2], 'k--','linewidth',3, 'Displayname','Reference angle');
plot(rad2deg([theta_s,theta_s]),[-50,2], 'k:','linewidth',3, 'Displayname','Steering angle');
lgnd = legend('Location','northwest','AutoUpdate','off');
set(lgnd,'color','w');

plot(th0,db(P_nom(theta_idx))/2,           'ko','markersize',16,'MarkerFaceColor',[0.1 0.8 0.1], 'Displayname','P');
plot(th0,db(abs(sum(extreme_points_min))), 'kv','markersize',16,'MarkerFaceColor',[0.3 0.3 1]*0.5, 'Displayname','P');
plot(th0,db(abs(sum(extreme_points_max))), 'k^','markersize',16,'MarkerFaceColor',[1   0   0]*0.5, 'Displayname','P');

ax = gca;
%exportgraphics(ax,'Fig6a_Backtracked_beampattern.pdf')%,'Resolution',300) 


%% Bounds and Monte Carlo of above - Fig. 1

figure(4);clf
set(gcf,'Position',[600 200 800 400])
set(gca, 'fontsize', 18)
hold on; xlim([-90-0.1,90]); ylim([-37,2]); grid on; set(gcf, 'color', 'white'); set(gca, 'color', [.85 .85 .85]);
xlabel('\theta (deg)')
ylabel('Power (dB)')
yticks([-40 -30,-20,-10, 0])
xticks([-90,-45,0,45,90])

pgon = polyshape([degs; flip(degs)], [db(P_GENPOLY(:,1)+ eps)/2 ; flip(db(P_GENPOLY(:,2))/2 + eps)]);
plot(pgon, 'linewidth', 1, 'edgealpha', 0, 'edgecolor',[0 0 0], 'FaceColor',[1 1 1],'FaceAlpha',1, 'handlevisibility', 'off')

degs = rad2deg(thetas);
p = plot(degs, db(P_GENPOLY(:,2))/2, '-','linewidth',5, 'color',[1,0,0], 'Displayname','Upper bound');

N_MC = 1000;
for n = 1:N_MC
    P_MC = MonteCarlo_P(bp_pol);
    p = plot(degs, db(P_MC)/2, 'linewidth',1, 'color',[0 0 0, 0.125], 'handlevisibility', 'off');
end

p = plot(degs, db(P_nom(:,1))/2, 'linewidth',5, 'color',[1 1 1], 'handlevisibility', 'off');
p = plot(degs, db(P_nom(:,1))/2, 'linewidth',4, 'color',[0.1,0.75,0.1], 'Displayname','Nominal');
p = plot(degs, db(P_GENPOLY(:,1))/2, 'linewidth',3, 'color',[0,0,1], 'Displayname','Lower bound');
p = plot([0,0], [0,0], 'linewidth',1, 'color',[0 0 0], 'displayname', strcat(num2str(N_MC), ' realizations'));
p = plot(degs, db(P_nom(:,1) .* B1 + B2)/2, 'linewidth',4, 'linestyle', ':', 'color',[255, 172, 28]/255, 'Displayname','Expectation');

plot(rad2deg([theta_s,theta_s]),[-50,2], 'k:','linewidth',2, 'Displayname','Steering angle');
lgnd = legend('Location','northwest','AutoUpdate','off');
set(lgnd,'color','w');

ax = gca;
%exportgraphics(ax,'MC_beampattern.jpg', 'Resolution',300)

%% REPRESENTAIONS - Fig. 2
nominal = 0.7*exp(1j*deg2rad(45+15));

radial_interval = ciat.RealInterval(0.50,0.90);
angular_interval = ciat.RealInterval( deg2rad(30+15-5),deg2rad(60+15+5));

cPOL = ciat.PolarInterval(radial_interval, angular_interval);
cREC = ciat.RectangularInterval(cPOL);
cCIRC = ciat.CircularInterval(cPOL);
cPOLY = ciat.PolygonalInterval(cPOL,'tolerance',3e-2);

% nominal + rect
fig5 = figure(5);clf
fig5.Position = [200 200 400 400]; 
hold on; axis equal; xlim([0,1]); ylim([0,1]); 
grid on; set(gcf, 'color', 'white');
xlabel('Real axis')
ylabel('Imaginary axis')
set(gca,'FontSize',15)
yticks([0:5])
xticks([1:5])

quiver(-1,-1,-1-real(nominal), -1-imag(nominal), 'linewidth',4, 'color','k','Maxheadsize',0.3, 'autoscale','off', 'Displayname', 'Nominal');
plot([-1],[-1],'color', [0 0 0], 'linewidth',4, 'Displayname', 'Interval bounds')
legend('Location','southeast','AutoUpdate','off')
plot(cPOL, 'color',[0 0 0], 'linewidth',4);
q1 = quiver(0,0,real(nominal), imag(nominal), 'linewidth',5, 'color','w','Maxheadsize',0.3, 'autoscale','off');
q2 = quiver(0,0,real(nominal), imag(nominal), 'linewidth',4, 'color','k','Maxheadsize',0.3, 'autoscale','off');
legend('Location','southeast','AutoUpdate','on')

h1 = plot(cREC,'color',[0 0.4470 0.7410], 'linewidth',4);
h1.set('DisplayName', 'Rectangular (rIA)');
h1.set('LineStyle', ':');

% circ
legend('Location','southeast','AutoUpdate','on')
legend('Location','southeast','AutoUpdate','on')
h2 = plot(cCIRC, 'color',[0.8500 0.3250 0.0980], 'linewidth',4, 'Displayname', 'Cirular (cIA)');

% poly
legend('Location','southeast','AutoUpdate','on')
h3 = plot(cPOLY, 'color',[0.9290 0.6940 0.1250], 'linewidth',4, 'linestyle', '--', 'marker', 'o','markersize', 12,'MarkerEdgeColor','none', 'MarkerFaceColor', [0.9290 0.6940 0.1250]);
h3.set('Displayname', 'Polygonal (pIA)');
legend('Location','southeast','AutoUpdate','off');
legend boxoff

ax = gca;
%% 
%exportgraphics(ax,'IA_representation.pdf')

