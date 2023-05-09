%% Presentation slide 5
clear 
close all

%% Set simulation parameters and calculate beampattern

array = biat.SensorArray(   'ElCount',          5,...
                            'ElPitchRatio',     0.5,...
                            'ElDiameterRatio',  0,...
                            'Curvature',        0.2,...
                            'SteeringAngle',    deg2rad(10));                            

bp_nom = biat.BeamPattern(array,'nominal','BeamResolutionDeg',0.35);
P_nom = bp_nom.calculateBeamPattern;

%% Plot sensor array

fig1 = figure(1);clf;
array.plot;

%% Plot BPs

fig2 = figure(2);clf
fig2.Position = [200 200 600 240]; 
subplot(1,2,1);
hold on; axis equal; xlim([-0.2,0.2]); ylim([-0.2,0.2]); grid on; set(gcf, 'color', 'white');
xlabel('Real axis')
ylabel('Imaginary axis')
set(gca,'FontSize',9)
yticks([-0.2,-0.1,0,0.1,0.2])
xticks([-0.2,-0.1,0,0.1,0.2])

subplot(1,2,2);
hold on; xlim([-90,90]); ylim([-30,2]); grid on; set(gcf, 'color', 'white');
xlabel('Angle (deg)')
ylabel('Power (dB)')
set(gca,'FontSize',9)
yticks([-30,-20,-10, 0])
xticks([-90,-45,0,45,90])

thetas = rad2deg(bp_nom.BeamAngles);
L1 = 0.2;
zs = zeros(5,1);
vecs=zeros(5,1);
vecsum=0;

subplot(1,2,1);
q1 = quiver(zs,zs,real(vecs), imag(vecs), 'linewidth',3, 'color','k','Maxheadsize',0.05/L1, 'autoscale','off', 'displayname', 'Element response');
q2 = quiver(0,0,real(vecsum/5), imag(vecsum/5), 'linewidth',3, 'color',[0.1,0.75,0.1],'Maxheadsize',0.05/L1, 'autoscale','off','displayname', 'Array response (/5)' );
legend('Location','northwest','AutoUpdate','off')
subplot(1,2,2);
p = plot([1],[1], 'linewidth',3, 'color',[0.1,0.75,0.1], 'Displayname','Power pattern');
plot([10,10],[-50,2], 'k:','linewidth',2, 'Displayname','Steering angle');
legend('Location','northwest','AutoUpdate','off')  


for t = 1 : length(thetas) % Loop over angles
    bp_nom.BeamIndex = t;
    vecs = bp_nom.ElementIntervals;
    vecsum = bp_nom.ArrayInterval;
    
    subplot(1,2,1);
    set(q1,'udata',real(vecs),'vdata',imag(vecs))
    set(q2,'udata',real(vecsum/5),'vdata',imag(vecsum/5))
    subplot(1,2,2);
    set(p, 'xdata', thetas(1:t), 'ydata', 10*log10(abs(P_nom(1:t))));
    drawnow;
end


%%
subplot(1,2,2)
xlabel('angle (degs)')
xlim([-90,90])
ylim([-60,5])

ylabel('Amplitude [dB]')
set(gca,'LineWidth',2)
set(gca,'DefaultLineLineWidth',2)
set(gca,'defaultAxesFontSize',20)

plot(thetas, db(P_nom(:,1))/2, 'k', 'DisplayName', 'Nominal');


legend();
