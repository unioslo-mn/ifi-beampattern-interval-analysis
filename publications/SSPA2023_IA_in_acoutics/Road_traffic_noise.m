%% Calculate road traffic noise level from hourly traffic volume only

% Calculate the result of a non-linear function over an interval using both 
% sampled approach and interval arithmetic

clear
close all

% Set traffic volume
qInt = ciat.RealInterval(2,4);

% Define function and calculate value
L10Func = @(q) 42.2 + 10*log10(q);
L10Val = L10Func(qInt);

% Plot 
figure(1);clf;hold on
fplot(L10Func,'DisplayName','function')
h = fplot(@(q) L10Func(q),[qInt.Infimum,qInt.Supremum],'r-',...
                         'LineWidth',4,...
                         'HandleVisibility','off');
L10Int = [min(h.YData),max(h.YData)];
plot([qInt.Infimum,qInt.Supremum],[30 30],'r-','LineWidth',4,...
                        'DisplayName','Sampled')
plot([0 0],L10Int,'r-','LineWidth',4,...
                         'HandleVisibility','off')
plot([0 0],[L10Val.Infimum,L10Val.Supremum],'r:','LineWidth',4,...
                        'DisplayName','Interval')
xlabel('q : Traffic volume [vehicle/hour]')
ylabel('f(q) : Basic noise level')
title('Noise level from traffic volume,',...
      'no dependency effect')
grid on
ylim([30 50])
legend()

%% Calculate noise level correction from speed and heavy vehicle ratio

% Calculate the result of a non-linear function over an interval using both 
% sampled approach and interval arithmetic

% Set speed and heavy vehicle ratio
VInt = ciat.RealInterval(60,80);
pInt = ciat.RealInterval(40,60);

% Define function and calculate value
corrFunc = @(p,V) 33*log10(V + 40 + 500./V) + 10.*log10(1+5*p./V) - 68.8;

% Plot result 
figure();clf
[corrSmp, corrInt] = plotIntervalFunction(corrFunc, pInt, VInt,[0 100 20 100])

% Set custom parameters
xlabel('p : Heavy vehicle ratio [%])')
ylabel('V : Traffic speed [km/h]')
zlabel('f(p,V)')
view([-35 20])
title('Noise correction from speed and heavy vehicle ratio,',...
      'dependency effect is apparent')
grid on
legend('Location','best')
