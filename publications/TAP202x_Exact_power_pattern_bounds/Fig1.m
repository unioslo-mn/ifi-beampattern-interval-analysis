%% Initilize workspace
clear 
close all

%% Set simulation parameters (include update function)
% Note: there is a certain order in which things must be set as of now
% Initialize simularion parameters
array = biat.SensorArray(   'ElCount',8,...
                            'ElDiameterRatio',0,...
                            'Curvature',0.0,...
                            'TaperType','chebwin',...
                            'TaperParam',25,...
                            'GainError',2/100,...
                            'PhaseError',deg2rad(2),...
                            'SteeringAngle',deg2rad(0));

bp_nom = biat.BeamPattern(array,'nominal','BeamResolutionDeg',0.35/2);
bp_pol = biat.BeamPattern(array,'polygonal','BeamResolutionDeg',0.35/2,...
                                            'PolygonTolerance',1e-6);

% Set beam index
beamAngle = deg2rad(39.5);
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
theta_probe = deg2rad(39.5); % probe 50 degrees (almost on sidelobe in plot)
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

%% Bounds and Monte Carlo of above - Fig. 1

fig = figure(4);clf
set(gcf,'Position',[600 200 800 400])
set(gca, 'fontsize', 18)
hold on; xlim([-90-0.1,90]); ylim([-47,2]); grid on; set(gcf, 'color', 'white'); set(gca, 'color', [.85 .85 .85]);
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
%p = plot(degs, db(P_nom(:,1) .* B1 + B2)/2, 'linewidth',4, 'linestyle', ':', 'color',[255, 172, 28]/255, 'Displayname','Expectation');

p = plot(degs, db(P_backtrack)/2, 'linewidth',3, 'color',[1 1 1], 'handlevisibility', 'off');
p = plot(degs, db(P_backtrack)/2, '-.','linewidth',2, 'color',[0.6, 0, 0], 'Displayname','Upper bound realization');
p = plot(degs, db(P_backtrack_min)/2, 'linewidth',3, 'color',[1 1 1], 'handlevisibility', 'off');
p = plot(degs, db(P_backtrack_min)/2, '-.','linewidth',2, 'color',[0.1, 0.8, 1], 'Displayname','Lower bound realization');

p = plot(degs(theta_idx), db(P_backtrack_min(theta_idx))/2, 'o','handlevisibility','off','MarkerFaceColor',[0.1, 0.8, 1],'markeredgecolor','k','MarkerSize',10);
p = plot(degs(theta_idx), db(P_backtrack(theta_idx))/2, 'o','handlevisibility','off','MarkerFaceColor',[0.6, 0, 0],'markeredgecolor','k','MarkerSize',10);


%plot(rad2deg([theta_s,theta_s]),[-50,2], 'k:','linewidth',2, 'Displayname','Steering angle');
lgnd = legend('Location','northwest','AutoUpdate','off');
set(lgnd,'color','w');

ax = gca;


%exportgraphics(ax,'/Users/havarn/Documents/DocGit/PolyArc/ifi-beampattern-interval-analysis/publications/TAP_Tight_bounds_202x/Arnestad_plots/MC_beampattern.jpg','Resolution',300)
