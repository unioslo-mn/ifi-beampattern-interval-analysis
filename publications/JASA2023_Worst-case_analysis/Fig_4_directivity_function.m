clear
close all;

% Set parameters
theta = deg2rad((-120:0.1:120)');
dth = deg2rad(10);
lambda = 1;
elRad = 0.5;
taperAngs = deg2rad([80,100]);

% Generate sensor array
array = biat.SensorArray('TaperAngles',deg2rad([80,100]));
array.CenterFrequency = array.SoundSpeed / lambda;
array.ElDiameterRatio = 2*elRad / array.ElPitch;

% Calculate directivity
DIR = array.getElDirectivity(theta);
theta_I = ciat.RealInterval(theta -dth, theta+dth);
temp = array.getElDirectivity(theta_I);
dir_I = [temp.Infimum, temp.Supremum];

% Plot results
fig = figure(1);
fig.Position=[200 200 500 300]; 
set(gca, 'fontsize', 16)
hold on; 
set(gcf, 'color', 'white'); xlim([-120,120]); ylim([0,1.05]);
xlabel('\Delta\theta (deg)')
ylabel('{\it d}(\Delta\theta) (-)')
theta = rad2deg(theta);
plot(theta, dir_I(:,2), 'Displayname', 'Upper (\pm10\circ)', 'linewidth',4, 'color',[1,0,0])
plot(theta, DIR, 'Displayname', 'Nominal', 'linewidth',3, 'color',[0.1,0.75,0.1])
plot(theta, dir_I(:,1), 'Displayname', 'Lower (\pm10\circ)', 'linewidth',2, 'color',[0,0,1])
legend boxoff

ax = gca;
%% 
%exportgraphics(ax,'directivity_function.pdf')%,'Resolution',300) 
