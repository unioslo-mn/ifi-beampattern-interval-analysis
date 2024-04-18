%% Dependent bivariate non-monotonic function

% Demonstrative example of the dependence problem in interval arithmetic. 
% On the left side a bivariate function is calculated with variable 
% $a^I=[1,3]$ appearing both in the numerator and denominator in a fully 
% correlated way and $b^I=[1,2]$ appearing only in the denominator. 
% On the left side the denominator is represented by variable 
% $c^I=a^I+b^I$ and the result is calculated by considering it independent 
% from the numerator. Both results provide reliable bounds, but the latter 
% provide relaxed bounds by containing the independent solutions.

clear
close all

% Define intervals
a_int = ciat.RealInterval(1, 3);

%Define function
f = @(a1,a2) a1./a2;

% Plot results
figure();
subplot(2,2,1);hold on
[f_smp,f_int] = plotIntervalFunction(f,a_int,a_int)

% Plot constraned function
line([-5 5],[-5 5],[1,1],'LineWidth',5,'Color','k',...
                            'DisplayName','Constraint');

                        
% Zoom in and reposition lines
ax = gca();
xlim([-1 5])
ylim([-1 5])
zlim([-1 10])
ax.Children(2).XData = [-1 -1];
ax.Children(3).ZData = [1 1];
ax.Children(3).Marker = 'o'
ax.Children(3).Color = 'k'
ax.Children(3).XData = [-1 -1];
ax.CLim = [-1 10]
plot3(a_int.Bounds,[-1 -1],[-1 -1],'r-','LineWidth',3,'HandleVisibility','off')
plot3([-1 -1],a_int.Bounds,[-1 -1],'r-','LineWidth',3,'HandleVisibility','off')
                        
% Add custom labels
xlabel('a_1 ')
ylabel('a_2 ')
zlabel('f(a1,a2)=a1/a2')
view([-80 30])
grid on
title('Dependence constraint on f(a)=a/a')
legend('Location','best')
set(gca,'FontSize',16)

%%

% Define functions
f = @(a,b) a./(a+b);

% Set intervals
a_int = ciat.RealInterval(1, 3);
b_int = ciat.RealInterval(1, 2);

% Set plot axis limits
a_lim = [0 5];
b_lim = [0 5];
f_lim = [0 2];

% Plot results
subplot(2,2,2);hold on
[f_smp,f_int] = plotIntervalFunction(f,a_int,a_int,[0.1 5 0.1 5])
     
% Add custom labels
ax = gca();
ax.Children(2).Color = 'k'
xlabel('a')
ylabel('b')
zlabel('f(a,b)=a/(a+b)')
title('Dependency effect on f(a,b)=a/(a+b)')
view([-80 30])
grid on
set(gca,'FontSize',16)
legend('Location','best')
