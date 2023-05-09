%% Paper figure 1, Presentation slide 6

close all; clear;
%%

nominal = 0.75*exp(1j*deg2rad(45+15));

radial_interval = ciat.RealInterval(0.6,0.9);
angular_interval = ciat.RealInterval( deg2rad(30+15),deg2rad(60+15));

cPOL = ciat.PolarInterval(radial_interval, angular_interval);
cREC = ciat.RectangularInterval(cPOL);
cCIRC = ciat.CircularInterval(cPOL);
cPOLY = ciat.PolygonalInterval(cPOL,'tolerance',1e-2);

%% nominal + rect
fig1 = figure('Position',[200 200 400 400]); hold on; axis equal; xlim([0,1]); ylim([0,1]); grid on; set(gcf, 'color', 'white');
xlabel('Real axis')
ylabel('Imaginary axis')
set(gca,'FontSize',14)
yticks([0:5]/5)
xticks([1:5]/5)
quiver(-1,-1,-1-real(nominal), -1-imag(nominal), 'linewidth',3, 'color','k','Maxheadsize',0.3, 'autoscale','off', 'Displayname', 'Nominal');
plot([-1],[-1], 'b', 'linewidth',3, 'Displayname', 'Interval bounds')
legend('Location','southeast','AutoUpdate','off')
plot(cPOL, 'color','b', 'linewidth',3);
q1 = quiver(0,0,real(nominal), imag(nominal), 'linewidth',4, 'color','w','Maxheadsize',0.3, 'autoscale','off');
q2 = quiver(0,0,real(nominal), imag(nominal), 'linewidth',3, 'color','k','Maxheadsize',0.3, 'autoscale','off');
legend('Location','southeast','AutoUpdate','on')

%saveas(gcf,'Elem_bounds_nominal.png')

h1 = plot(cREC, 'color','r', 'linewidth',3);
h1.set('Displayname', 'Rectangular repr.');

%saveas(gcf,'Elem_bounds_rect.png')

%% circ
fig2 = figure('Position',[200 200 400 400]); hold on; axis equal; xlim([0,1]); ylim([0,1]); grid on; set(gcf, 'color', 'white');
xlabel('Real axis')
ylabel('Imaginary axis')
set(gca,'FontSize',14)
yticks([0:5]/5)
xticks([1:5]/5)
quiver(-1,-1,-1-real(nominal), -1-imag(nominal), 'linewidth',3, 'color','k','Maxheadsize',0.3, 'autoscale','off', 'Displayname', 'Nominal');
plot([-1],[-1], 'b', 'linewidth',3, 'Displayname', 'Interval bounds')
plot([-1],[-1], 'r', 'linewidth',3, 'Displayname', 'Circular repr.') % dirty fix to below
legend('Location','southeast','AutoUpdate','off')
plot(cPOL, 'color','b', 'linewidth',3);
q1 = quiver(0,0,real(nominal), imag(nominal), 'linewidth',4, 'color','w','Maxheadsize',0.3, 'autoscale','off');
q2 = quiver(0,0,real(nominal), imag(nominal), 'linewidth',3, 'color','k','Maxheadsize',0.3, 'autoscale','off');
legend('Location','southeast','AutoUpdate','on')
h2 = plot(cCIRC, 'color','r', 'linewidth',3);
%h2.set('Displayname', 'Ciruclar repr.'); % not working.

%saveas(gcf,'Elem_bounds_circ.png')

%% polygon
fig = figure('Position',[200 200 400 400]); hold on; axis equal; xlim([0,1]); ylim([0,1]); grid on; set(gcf, 'color', 'white');
xlabel('Real axis')
ylabel('Imaginary axis')
set(gca,'FontSize',14)
yticks([0:5]/5)
xticks([1:5]/5)
quiver(-1,-1,-1-real(nominal), -1-imag(nominal), 'linewidth',3, 'color','k','Maxheadsize',0.3, 'autoscale','off', 'Displayname', 'Nominal');
plot([-1],[-1], 'b', 'linewidth',3, 'Displayname', 'Interval bounds')

legend('Location','southeast','AutoUpdate','off')
plot(cPOL, 'color','b', 'linewidth',3);
q1 = quiver(0,0,real(nominal), imag(nominal), 'linewidth',4, 'color','w','Maxheadsize',0.3, 'autoscale','off');
q2 = quiver(0,0,real(nominal), imag(nominal), 'linewidth',3, 'color','k','Maxheadsize',0.3, 'autoscale','off');
legend('Location','southeast','AutoUpdate','on')

h3 = plot(cPOLY, 'color','r', 'linewidth',3);
h3.set('Displayname', 'Polygon repr.');
legend('Location','southeast','AutoUpdate','off')
h4 = plot(cPOLY.Points,'o', 'Color','k','markersize', 10, 'MarkerFaceColor', 'r');

%saveas(gcf,'Elem_bounds_poly.png')
