%% Non-monotonic independent bivariate function

clear
close all

% Define intervals
a_int = ciat.RealInterval(pi/3,3/4*pi);
b_int = ciat.RealInterval(3,6);

% Define function
f = @(a,r) r.*sin(a);

% Plot results
figure(1);hold on
[f_smp,f_int] = plotIntervalFunction(f,a_int,b_int,[-pi pi 0 10])


% Add custom labels
xlabel('a : Measured angle  [rad])')
ylabel('r : Measured distance [m]')
zlabel('f(a,r)=r*sin(a): Elevation')
view([25.6 27.5])
grid on
title('Non-monotonic independent bivariate function',...
      'no dependency effect')
legend('Location','best')

%% Non-monotonic independent bivariate monotonic function
clear

% Define intervals
a_int = ciat.RealInterval(-1,3);
b_int = ciat.RealInterval(-2, 2);

% Define function
f = @(a,b) a.*b;

% Plot results
figure(2);hold on
[f_smp,f_int] = plotIntervalFunction(f,a_int,b_int)

% Add custom labels
xlabel('a ')
ylabel('b ')
zlabel('f(a,b)=a*b')
view([-39.4 22.2])
grid on
title('Non-monotonic independent bivariate monotonic function',...
      'no dependency effect')
legend('Location','best')

%% Dependent bivariate non-monotonic function
clear

% Define intervals
a_int = ciat.RealInterval(1, 3);
b_int = ciat.RealInterval(2, 4);

%Define function
f = @(a,b) a./(a+b);

% Plot results
figure(3);hold on
[f_smp,f_int] = plotIntervalFunction(f,a_int,b_int,[1 10 1 10])


% Add custom labels
xlim([0 10])
ylim([0 10])                                            
xlabel('a ')
ylabel('b ')
zlabel('f(a,b)=a/(a+b)')
view([-111.733 27.533])
grid on
title('Dependent bivariate non-monotonic function',...
      'clear dependency effect')
legend('Location','best')
  
%% Dependent bivariate monotonic function with aligned partial derivatives
clear

% Define intervals
a_int = ciat.RealInterval(3, 7);
b_int = ciat.RealInterval(2, 4);

%Define function
f = @(a,b) log(a) .* (a+b);

% Plot results
figure(4);hold on
[f_smp,f_int] = plotIntervalFunction(f,a_int,b_int,[1 10 1 10])

% Add custom labels
xlim([1 10])
ylim([0 10])
xlabel('a ')
ylabel('b ')
zlabel('f(a,b)=log(a)*(a+b)')
view([-51.733 27.533])
grid on
title('Dependent bivariate monotonic function with aligned partial derivatives',...
      'no dependency effect')
legend('Location','best')

%% Dependent bivariate monotonic function with non-aligned partial derivatives
clear


% Define intervals
a_int = ciat.RealInterval(30, 60);
b_int = ciat.RealInterval(50, 80);

%Define function
f = @(a,b) 3*log10(b + 500./b) + log10(1+a./b);

% Plot results
figure(5);hold on
[f_smp,f_int] = plotIntervalFunction(f,a_int,b_int,[10 100 10 100])

% Add custom labels
xlim([10 100])
ylim([10 100])
xlabel('a ')
ylabel('b ')
zlabel('f(a,b)=3*log10(b+500/b)+log10(1+a/b)')
view([-70 20])
grid on
title('Dependent bivariate monotonic function with non-aligned partial derivatives',...
      'clear dependency effect')
legend('Location','best')
 
%% Simplest example of dependent multivariate function with dependency effect

clear


% Define intervals
a_int = ciat.RealInterval(1, 3);

%Define function
f = @(a1,a2) a1./a2;

% Plot results
figure(6);hold on
[f_smp,f_int] = plotIntervalFunction(f,a_int,a_int)
ax = gca();
ax.Children(2).ZData = [1 1];

% Plot constraned function
fimplicit3(@(a1,a2,f) a1-a2 ,'FaceColor','none',...
                             'DisplayName','Constraint');
line([-5 5],[-5 5],[1,1],'LineWidth',5,'Color','k',...
                            'DisplayName','Constrained');

% Add custom labels
% xlim([10 100])
% ylim([10 100])
xlabel('a1 ')
ylabel('a2 ')
zlabel('f(a1,a2)=a1/a2')
view([-80 30])
grid on
title('Simplest example of dependent multivariate function with dependency effect')
legend('Location','best')
