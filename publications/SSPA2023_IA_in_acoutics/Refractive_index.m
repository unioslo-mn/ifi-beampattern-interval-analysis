%% Measurement of the refractive index

% Let us assume a measurement setup where the incident angle and reflected
% angle are both measured with a give precision and the aim is to estimate
% the refraction index of an interface between two media

clear
close all

% Define refractive index function
nFunc = @(th1,th2) sin(th1)./sin(th2);

% Define angle measurement error
thErr = deg2rad(5);

% Define measured values
th1Nom = deg2rad(15);
th2Nom = deg2rad(18);

% Define intervals
th1Int = ciat.RealInterval(th1Nom+[-1 1]*thErr);
th2Int = ciat.RealInterval(th2Nom+[-1 1]*thErr);

% Calculate and plot result (sampled and interval arithmetic)
[nSmp nInt] = plotIntervalFunction(nFunc, th1Int, th2Int,...
                                            [-pi/2,pi/2,-pi/2,pi/2])

% Set custom parameters
view(-70,60)
xlabel('Incident angle')
ylabel('Reflected angle')
zlabel('Refractive index')
grid on
legend('Location','best')
