
%% REQUIREMENTS
% k-wave in path to have access to the function "acousticFieldPropagator": http://www.k-wave.org/documentation.php
clearvars;

%% SIMULATION
% There are some blocks to run before the plotting is done

% define the computational grid
Nx = 192*3+1;           % number of grid points in the x (row) direction
Ny = 192*3+1;           % number of grid points in the y (column) direction
dx = 0.5;    	        % grid point spacing [m]

% define source frequency and sound speed
c0 = 1500;          % sound speed [m/s]
lambda = 1;
f0 = c0/lambda;     % source frequency [Hz]

% define source aperture as a line array at the top of the grid
aperture_width = 31; % elements

% calculate position of each grid point in the array
el_pos = dx * (0:aperture_width - 1) - (aperture_width-1)/2*dx;

% focus and steering
steer = deg2rad(-10);

% some settings
apod_window = chebwin(aperture_width, 30);
theta_probe = deg2rad(13.6);

% errors
gamma = 0.05;
error_g = 5/100;
error_phi = deg2rad(5);


%% Interval Arithmetic

array = biat.SensorArray(   'ElCount',aperture_width,...
                            'ElDiameterRatio',0,...
                            'SoundSpeed',c0,...
                            'CenterFrequency',f0,...
                            'ElPitchRatio',dx/lambda,...
                            'Curvature',0,...
                            'SteeringAngle',steer,...
                            'GainError',error_g,...
                            'PhaseError',error_phi,...
                            'CouplingCoeff',gamma,...
                            'TaperType','chebwin',...
                            'TaperParam',30);
                        
bp_nom = biat.BeamPattern(array,'nominal','BeamResolutionDeg',0.35/2);
bp_pol = biat.BeamPattern(array,'polygonal','BeamResolutionDeg',0.35/2,...
                                            'PolygonTolerance',1e-4);
                            
% two lists I unfortunately have made myself dependent on (plotting mostly, might fix)
thetas = bp_nom.BeamAngles;
degs = rad2deg(thetas);

% Calculate beampattern
P_nom = bp_nom.calculateBeamPattern;
P_GENPOLY = bp_pol.calculateBeamPattern;

% Backtrack beampattern
bp_pol.PolygonTolerance = 1e-6;
[~,bp_nom.BeamIndex] = min(abs(bp_nom.BeamAngles - theta_probe));
bp_pol.BeamIndex = bp_nom.BeamIndex;
th0 = rad2deg(bp_nom.IncidenceAngle);
extreme_points_max = bp_pol.MaxPowerPoints;
shapes = bp_pol.ElementIntervals;
element_error = bp_pol.MaxErrorPattern;
[~,C_extr] = bp_pol.backtrackComplexIntervals('getMax',1);
P_backtrack = bp_pol.backtrackBeamPattern('getMax',1);


%% COMPUTE HISTOGRAM

N_MC = 100000;
P_MC = zeros(N_MC, 1);
for n = 1:N_MC
    % NOTE: theta in argument makes us only compute for that angle
    P_MC(n) = MonteCarlo_P(bp_pol, 'theta', bp_pol.IncidenceAngle); 
end


%% COMPUTE 2D COLOR BPs
% open new figure window

% allocate empty amplitute and phase matrices
amp_in_NOM   = zeros(Nx, Ny);
phase_in_NOM = zeros(Nx, Ny);
amp_in       = zeros(Nx, Ny);
phase_in     = zeros(Nx, Ny);

y1 = Ny/2 - aperture_width/2 + 1;
y2 = Ny/2 + aperture_width/2;
x1 = 1;

rads_m90_top90 = deg2rad(-90:1:90);

axX =([0: Nx-1])*dx;
axY =([-Ny/2: Ny/2-1]+0.5)*dx;

steer_phase = -2 * pi * f0 * el_pos * sin(steer) / c0;

% sjekk om C_extr er riktig vei, fight c og m.
Aworst = (C_extr * (apod_window .* exp(1j*steer_phase.')));
err = Aworst .* element_error;

% assign constant amplitude across the line array
amp_in_NOM(x1, y1:y2) =  apod_window.';
amp_in(x1, y1:y2) = abs(err); %element_error) .* apod_window.';

% calculate phase offset for each grid point in the line array based on
% element position and steering angle, and assign to the line array
phase_in_NOM(x1, y1:y2) = steer_phase ; % 2 * pi * f0 * el_pos * sind(steering_angle) / c0 +
phase_in(x1, y1:y2) = angle(err); %steer_phase + angle(element_error); % 2 * pi * f0 * el_pos * sind(steering_angle) / c0 +

% compute beam pattern
[amp_out_NOM, phase_out_NOM] = acousticFieldPropagator(amp_in_NOM, phase_in_NOM, dx, f0, c0);
[amp_out, phase_out] = acousticFieldPropagator(amp_in, phase_in, dx, f0, c0);

%% COMPUTE Statistics / approximation

% Get parameters
g_I = array.GainInterval;
phi_I = array.PhaseInterval;
rx_I = array.PosXInterval;
theta_idx = bp_pol.BeamIndex;

T_se = sum(apod_window.^2);
var_g = 1/12 * (g_I(1).Width)^2;
var_phase = 1/12 * (phi_I(1).Width)^2;
var_pos_lambda = (2*pi / (c0/f0))^2  * 1/12 * (rx_I(1).Width)^2;

B1 = exp(-(var_phase + var_pos_lambda));
B2 = T_se * ( 1 + var_g - exp(-(var_phase + var_pos_lambda)));

%% Fig9 - BOUND PLOT BP

figure(2);clf
set(gcf,'Position',[600 200 800 500])
set(gca, 'fontsize', 23)
hold on; 
xlim([-30 25] ); 
ylim([-50,2]); grid on; set(gcf, 'color', 'white'); set(gca, 'color', [.85 .85 .85]);
xlabel('\theta (deg)')
ylabel('Power (dB)')
yticks([-50:10:0])
xticks(-50:10:50);

% Plot interval bounds as white background
pgon = polyshape([degs; flip(degs)], [db(P_GENPOLY(:,1)+ eps)/2 ; flip(db(P_GENPOLY(:,2))/2 + eps)]);
plot(pgon, 'linewidth', 1, 'edgealpha', 0, 'edgecolor',[0 0 0], 'FaceColor',[1 1 1],'FaceAlpha',1, 'handlevisibility', 'off')

% plot nominal, bounds, steering, and add legend
p = plot(degs, db(P_GENPOLY(:,2))/2, '-','linewidth',4, 'color',[1,0,0], 'Displayname','Upper bound');
p = plot(degs, db( P_nom(:,3),'power'), 'linewidth',4, 'linestyle', ':', 'color',[0.4,0.2,0.7], 'Displayname','Approximate bound');
p = plot(degs, db(P_nom(:,1))/2, 'linewidth',5.8, 'color',[1 1 1], 'handlevisibility', 'off');
p = plot(degs, db(P_nom(:,1))/2, 'linewidth',5, 'color',[0.1,0.75,0.1], 'Displayname','Nominal');
p = plot(degs, db(P_backtrack)/2, 'linewidth',4,'linestyle','-.', 'color',[1,0,0], 'Displayname','Worst case');

plot([th0,th0],[-50,2], 'k--','linewidth',3, 'Displayname','Reference angle');

lgnd = legend('Location','northwest','AutoUpdate','off');
lgnd.BoxFace.ColorType='truecoloralpha';
lgnd.BoxFace.ColorData=uint8(255*[1 1 1 0.85]');

%exportgraphics(gcf,'Fig9_bounds_example_B.pdf')


%% Fig 10 - Define plot

set(groot,'defaultAxesFontSize',15)

figure(5);
set(gcf,'Position',[300 200 1000 460])
ha = tight_subplot(2,4,[.1 .05],[.1 .1],[.06 .1]);

elAx = [1:aperture_width];
elAx2 = [0 : 10 : aperture_width];
elAx2(1) = 1;

% PANEL 1
axes(ha(1)); hold on;
annotation('textbox', [ha(1).Position(1:2),0,0] + [0.01 0.28 0.03 0.05] , 'String', "(a)",'BackgroundColor','w', 'FitBoxToText','off', 'FontSize', 14, 'Fontweight', 'bold', 'Facealpha', 0.85);
annotation('textbox', [ha(2).Position(1:2),0,0] + [0.01 0.28 0.03 0.05] , 'String', "(b)",'BackgroundColor','w', 'FitBoxToText','off', 'FontSize', 14, 'Fontweight', 'bold', 'Facealpha', 0.85);
annotation('textbox', [ha(3).Position(1:2),0,0] + [0.01 0.28 0.09 0.05] , 'String', "(c) - arg({\itC})",'BackgroundColor','w', 'FitBoxToText','off', 'FontSize', 14, 'Fontweight', 'bold', 'Facealpha', 0.85);
annotation('textbox', [ha(4).Position(1:2),0,0] + [0.01 0.28 0.09 0.05] , 'String', "(d) - abs({\itC})",'BackgroundColor','w', 'FitBoxToText','off', 'FontSize', 14, 'Fontweight', 'bold', 'Facealpha', 0.85);
annotation('textbox', [ha(5).Position(1:2),0,0] + [0.01 0.28 0.03 0.05] , 'String', "(e)",'BackgroundColor','w', 'FitBoxToText','off', 'FontSize', 14, 'Fontweight', 'bold', 'Facealpha', 0.85);
annotation('textbox', [ha(6).Position(1:2),0,0] + [0.01 0.28 0.08 0.05] , 'String', "(f) - Nom.",'BackgroundColor','w', 'FitBoxToText','off', 'FontSize', 14, 'Fontweight', 'bold', 'Facealpha', 0.85);
annotation('textbox', [ha(7).Position(1:2),0,0] + [0.01 0.28 0.09 0.05] , 'String', "(g) - Worst",'BackgroundColor','w', 'FitBoxToText','off', 'FontSize', 14, 'Fontweight', 'bold', 'Facealpha', 0.85);
annotation('textbox', [ha(8).Position(1:2),0,0] + [0.01 0.28 0.07 0.05] , 'String', "(h) - Diff.",'BackgroundColor','w', 'FitBoxToText','off', 'FontSize', 14, 'Fontweight', 'bold', 'Facealpha', 0.85);

set(gcf, 'color', 'white');
xlabel('Element {\itm}')
ylabel('Amplitude error (%)')
bar(100*(abs(element_error)-1),'FaceColor',[.8 .3 .3],'EdgeColor',[.8 .3 .3],'handlevisibility', 'off');%'displayname', 'Worst realization')
plot((1:aperture_width),100*([g_I.Supremum]-1),'r:','linewidth',5, 'displayname', 'Upper bound')
plot((1:aperture_width),100*([g_I.Infimum]-1),'b:','linewidth',5, 'displayname', 'Lower bound')
ylim(100*[(g_I(1).Infimum-1)*1.4, (g_I(1).Supremum-1)*1.4]); xlim([0.5, aperture_width+0.5]);
xticks(elAx2)

% PANEL 2
axes(ha(2)); hold on;

xlabel('Element {\itm}')
ylabel('Phase error (deg)')
bar(rad2deg(angle(element_error)),'FaceColor',[.8 .3 .3],'EdgeColor',[.8 .3 .3],'displayname', 'Realization')
plot((1:aperture_width),rad2deg([phi_I.Supremum]),'r:','linewidth',5, 'displayname', 'upper bound')
plot((1:aperture_width),rad2deg([phi_I.Infimum]),'b:','linewidth',5, 'displayname', 'lower bound')
xlim([0.5, aperture_width+0.5]); 
ylim([rad2deg(phi_I(1).Infimum)-2, rad2deg(phi_I(1).Supremum) + 2]); 
xticks(elAx2)

% PANEL 3
axes(ha(3)); hold on;

imagesc(elAx,elAx, angle(C_extr) )
xlabel('Element {\itc}')
ylabel('Element {\itm}')
yticks(elAx2)
xticks(elAx2)
axis image;
set(gca,'YDir','reverse')

h = rectangle('Position', [ 21 21 10.1 10.1], ...
                'Curvature', 0.1, ...
                'FaceColor', [1, 1, 1, 0.8], ...
                'EdgeColor', [1, 1, 1, 0.8]);

axs = phasebar('location','se','rad');
axs.Position = axs.Position + [-0.004 0.01 0 0];

% PANEL 4
axes(ha(4)); hold on;
imagesc(elAx,elAx, 20*log10(abs(C_extr)) );
caxis([-50, 0])
set(gca,'YDir','reverse')
xlabel('Element {\itc}')
ylabel('Element {\itm}')
yticks(elAx2)
xticks(elAx2)
axis image;

% PANEL 5
axes(ha(5)); 

plot_type = 'db';

if strcmp(plot_type, 'power')
    list_MC = (P_MC);
    sup = ([P_GENPOLY(theta_idx,2),P_GENPOLY(theta_idx,2)]);
    inf = ([P_GENPOLY(theta_idx,1),P_GENPOLY(theta_idx,1)]);
    nom = ([P_nom(theta_idx,1),P_nom(theta_idx,1)]);        
    expc = ([P_nom(theta_idx,1),P_nom(theta_idx,1)]) .* B1 + B2 ;  
    xlabel_text = 'Power';
elseif strcmp(plot_type, 'amplitude')
    list_MC = sqrt(P_MC);
    sup = sqrt([P_GENPOLY(theta_idx,2),P_GENPOLY(theta_idx,2)]);
    inf = sqrt([P_GENPOLY(theta_idx,1),P_GENPOLY(theta_idx,1)]);
    nom = sqrt([P_nom(theta_idx,1),P_nom(theta_idx,1)]);
    expc = sqrt( [P_nom(theta_idx,1),P_nom(theta_idx,1)] .* B1 + B2 );
    xlabel_text = 'Amplitude'; 
elseif strcmp(plot_type, 'db')
    list_MC = 0.5 * db((P_MC));
    list_MC_fit =sqrt((P_MC));
    sup = 0.5 * db([P_GENPOLY(theta_idx,2),P_GENPOLY(theta_idx,2)]);
    inf = 0.5 * db([P_GENPOLY(theta_idx,1),P_GENPOLY(theta_idx,1)]  + eps);
    nom = 0.5 * db([P_nom(theta_idx,1),P_nom(theta_idx,1)]);   
    expc = 0.5 * db(([P_nom(theta_idx,1),P_nom(theta_idx,1)]) .* B1 + B2 );
    xlabel_text = 'dB'; 
end


h = histogram(list_MC, 70,'Normalization', 'pdf', 'FaceColor',[.25 .25 .25], 'EdgeColor',[0 0 0], 'displayname', 'Hist.');% 'handlevisibility', 'off',
hold on;

pgon = polyshape([0, 100, 100, 0] + sup(1), [0 0 1 1]);
plot(pgon, 'linewidth', 1, 'edgealpha', 0, 'edgecolor',[0 0 0], 'FaceColor',[.85 .85 .85],'FaceAlpha',1, 'handlevisibility', 'off')

plot(sup, [-1,max(h.Values)]*1.1, 'r-.','linewidth', 3', 'displayname', 'Worst')
plot(nom, [-1,max(h.Values)]*1.1, 'color',[.1 .8 .1],'linestyle', '-', 'linewidth', 3, 'displayname', 'Nom.')

xlim([-50, -5])
xlabel('Power (dB)')
ylabel('Pdf')
ylim([0, 1.05*max(h.Values)])
legend

% PANEL 6
axes(ha(6)); % COLOR NOMINAL
imagesc(axY, axX, db(amp_out_NOM / max(amp_out_NOM, [], 'all')));
caxis([-50, 0])
axis image;
xlabel('{\it x}/\lambda (-)');
ylabel('{\it y}/\lambda (-)');

% PANEL 7
axes(ha(7)); % COLOR MAX
% MAX
imagesc(axY, axX,  db(amp_out / max(amp_out_NOM, [], 'all')));
caxis([-50, 0])
axis image;
%title(['Worst case']);
xlabel('{\it x}/\lambda (-)');
ylabel('{\it y}/\lambda (-)');

% PANEL 8
axes(ha(8)); % COLOR DIFF

diff_img     = abs(amp_out_NOM .* exp(1j* phase_out_NOM) - amp_out .* exp(1j* phase_out));
imagesc(axY, axX, db(diff_img / max(amp_out_NOM, [], 'all')));
caxis([-50, 0])
axis image;
xlabel('{\it x}/\lambda (-)');
ylabel('{\it y}/\lambda (-)');
diff_img_ph     = (amp_out_NOM .* exp(1j* phase_out_NOM) - amp_out .* exp(1j* phase_out));
colormap('hot');    
h = axes(gcf,'visible','off'); 
c = colorbar(h,'Position',[0.92 0.2 0.022 0.6]);  % attach colorbar to h
colormap(c,'hot')
caxis(h,[-50,0]);
c.Title.String = "(dB)";

colormap( ha(3), twilight)
colormap( ha(4), hot )
colormap( axs, twilight)

%exportgraphics(gcf,'Fig10_eight_plot.jpg','Resolution',300) 
