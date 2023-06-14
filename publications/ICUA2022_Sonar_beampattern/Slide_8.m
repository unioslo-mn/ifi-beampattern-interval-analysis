%% Presentation slide 8
clear 
close all

%% Set simulation parameters (include update function)

array = biat.SensorArray(   'ElCount',          5,...
                            'ElPitchRatio',     0.5,...
                            'ElDiameterRatio',  0,...
                            'Curvature',        0.2,...
                            'SteeringAngle',    deg2rad(10),...
                            'GainError',        0/100,...
                            'PhaseError',       deg2rad(3));  

bp_nom = biat.BeamPattern(array,'nominal','BeamResolutionDeg',0.35/2);
bp_pol = biat.BeamPattern(array,'polygonal','BeamResolutionDeg',0.35/2,...
                                            'PolygonTolerance',0.001);
P_nom = bp_nom.calculateBeamPattern;
P_pol = bp_pol.calculateBeamPattern;

thetas = rad2deg(bp_nom.BeamAngles);

%% Plot array

fig1 = figure(1);clf
array.plot;


%% Plot BPs
fig2 = figure(2);clf
fig2.Position = [200 200 600*2 240*2]; 
subplot(1,2,1);
hold on; axis equal; xlim([-0.25,0.25]); ylim([-0.25,0.25]); grid on; 
set(gcf, 'color', 'white');
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

L1 = 0.2;
zs = zeros(5,1);
vecs=zeros(5,1);


subplot(1,2,1);
q1 = quiver(zs,zs,real(vecs), imag(vecs), 'linewidth',1, 'color','k',...
                            'Maxheadsize',0.05/L1, 'autoscale','off' ,...
                            'displayname', 'Element response');
q2 = quiver(0,0,real(sum(vecs)), imag(sum(vecs)), 'linewidth',3,...
                             'color',[0.1,0.75,0.1],'Maxheadsize',0.05/L1,...
                            'autoscale','off','displayname', 'Array response' );
p_LM = plot(1,1, 'o', 'Color','k','markersize', 8, ...
                            'MarkerFaceColor', 'b',...
                            'Displayname','P^L');
p_UM = plot(1,1, 'o', 'Color','k','markersize', 8, ...
                            'MarkerFaceColor', 'r',...
                            'Displayname','P^U');

%legend('Location','northwest','AutoUpdate','off')
p1 = plot(1,1, 'linewidth',2, 'color','k', 'Displayname','Element Int.');
p2 = plot(1,1, 'linewidth',3, 'color','y', 'Displayname','Array Int.');
p3 = plot(1,1, 'linewidth',2, 'color','k', 'Displayname','Array Int.');

subplot(1,2,2);
p = plot(1,1, 'linewidth',2, 'color',[0.1,0.75,0.1], ...
                            'Displayname','P_{nom.}');
p_L = plot(1,1, 'linewidth',2, 'color','b','linestyle','-', ...
                            'Displayname','P^L');
p_U = plot(1,1, 'linewidth',2, 'color','r','linestyle','-', ...
                            'Displayname','P^U');
legend('Location','northwest','AutoUpdate','off')  
plot([10,10],[-50,2], 'k:','linewidth',2, 'Displayname','Steering angle');


gif('anim.gif')
for t = 1 : length(thetas) % Loop over angles
    bp_nom.BeamIndex = t;
    bp_pol.BeamIndex = t;
    vecs = bp_nom.ElementIntervals;
    vecsum = bp_nom.ArrayInterval;
    elInt = bp_pol.ElementIntervals;
    arInt = bp_pol.ArrayInterval;
    
    subplot(1,2,1);
    set(q1,'udata',real(vecs),'vdata',imag(vecs))
    set(q2,'udata',real(sum(vecs)),'vdata',imag(sum(vecs)))
    
    
    p1 = elInt.plot('k','LineWidth',2);
    p2 = arInt.plot('linewidth',3.7, 'color','y');
    p3 = arInt.plot('linewidth',3, 'color','k');
    
    [vU,iU] = max(abs(arInt.Points));
    [vL,iL] = min(abs(arInt.Points));
    if vL == 0
        p_LM = plot(0,0, 'o', 'Color','k','markersize', 10, ...
                            'MarkerFaceColor', 'b','Displayname','P^L');
    else
        p_LM = plot(real(arInt.Points(iL)),imag(arInt.Points(iL)), 'o', 'Color','k',...
                            'markersize', 10, 'MarkerFaceColor', 'b',...
                            'Displayname','P^L');
    end
    p_UM = plot(real(arInt.Points(iU)),imag(arInt.Points(iU)), 'o', 'Color','k',...
                            'markersize', 10, 'MarkerFaceColor', 'r',...
                            'Displayname','P^U');
    
    subplot(1,2,2);
    set(p, 'xdata', thetas(1:t), ...
           'ydata', 10*log10(abs(P_nom(1:t))));
    set(p_L, 'xdata', thetas(1:t),  ...
             'ydata', 10*log10(abs(P_pol(1:t,1))));
    set(p_U, 'xdata', thetas(1:t),  ...
             'ydata', 10*log10(abs(P_pol(1:t,2))));
    
    drawnow;
    gif
    if t < length(thetas)
        delete(p1);
        delete(p2);
        delete(p3);
        delete(p_LM);
        delete(p_UM);
    end
end


%% Plot nominal beampattern

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
