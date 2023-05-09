%% Presentation slide 7
close all; clear; clc;
%%

radial_interval = ciat.RealInterval(0.3,0.4);
dl = 10;
du = 25;
d = 25;
angular_interval1 = ciat.RealInterval( deg2rad(dl),deg2rad(du));
angular_interval2 = ciat.RealInterval( deg2rad(dl+d),deg2rad(du+d));
angular_interval3 = ciat.RealInterval( deg2rad(dl+2*d),deg2rad(du+2*d));

cPOL1 = ciat.PolarInterval(radial_interval, angular_interval1);
cPOL2 = ciat.PolarInterval(radial_interval, angular_interval2);
cPOL3 = ciat.PolarInterval(radial_interval, angular_interval3);
cPOLY1 = ciat.PolygonalInterval(cPOL1,'tolerance',3e-4);
cPOLY2 = ciat.PolygonalInterval(cPOL2,'tolerance',3e-4);
cPOLY3 = ciat.PolygonalInterval(cPOL3,'tolerance',3e-4);

cPOLYs = cPOLY1 + cPOLY2;
cPOLYss = cPOLYs + cPOLY3;

% polygon
fig1 = figure(1);clf
fig1.Position = [200 200 400 400]; 
hold on; axis equal; xlim([0,1]); ylim([0,1]); grid on; set(gcf, 'color', 'white');
xlabel('Real axis')
ylabel('Imaginary axis')
set(gca,'FontSize',14)
yticks([0:5]/5)
xticks([1:5]/5)

h1 = plot(cPOLY1, 'color','c', 'linewidth',3);
h2 = plot(cPOLY2, 'color','m', 'linewidth',3);
h3 = plot(cPOLY3, 'color','y', 'linewidth',3);
h4 = plot(cPOLY1.Points,'o', 'Color','k','markersize', 10, ...
                             'MarkerFaceColor', 'c', 'Displayname','A');
h5 = plot(cPOLY2.Points,'o', 'Color','k','markersize', 10,  ...
                             'MarkerFaceColor', 'm', 'Displayname','B');
h6 = plot(cPOLY3.Points,'o', 'Color','k','markersize', 10,  ...
                             'MarkerFaceColor', 'y', 'Displayname','C');
h7 = plot(cPOLYs.Points,'-o', 'Color','k','markersize', 10,  ...
                             'MarkerFaceColor', 'b',  ...
                             'Displayname','A\oplusB');
h8 = plot(cPOLYss.Points,'-o', 'Color','k','markersize', 10,  ...
                             'MarkerFaceColor', 'k',  ...
                             'Displayname','A\oplusB\oplusC');

legend([h4,h5,h6,h7,h8],'Location','northwest','AutoUpdate','on')

%
points1 = 1j+[1+1j, 2+1.1j, 1.4+1.8j].';
poly1 = ciat.PolygonalInterval(points1);

points2 = 1+[1+1j, 1+2j,1.5+1j, 1.5+2j].';
poly2 = ciat.PolygonalInterval(points2);

poly3 = poly1+poly2;

%% Animation

idxmax = 60;

xr = 2*rand(1000,1) + 1;
yr = 2*rand(1000,1) + 1;

in1 = inpolygon(xr,yr,real(poly1.Points),imag(poly1.Points));
in2 = inpolygon(xr,yr,real(poly2.Points),imag(poly2.Points));
xp1 = (xr(in1));
yp1 = (yr(in1));
xp2 = (xr(in2));
yp2 = (yr(in2));

xp3 = xp1(1:idxmax) + xp2(1:idxmax);
yp3 = yp1(1:idxmax) + yp2(1:idxmax);
for idx = 1:idxmax
    clf;
    clear axis; hold on; axis equal; xlim([0,5]); ylim([0,5]); grid on;
    xlabel('Real axis')
    ylabel('Imaginary axis')
    set(gca,'FontSize',14)
    yticks([1:5])
    xticks([1:5])


    fill(real(poly1.Points),imag(poly1.Points),[0 0.4470 0.7410],  ...
                             'Displayname', 'Polygon A')
    fill(real(poly2.Points),imag(poly2.Points),[0.8500 0.3250 0.0980], ...
                             'Displayname', 'Polygon B')
    fill(real(poly3.Points),imag(poly3.Points),[0.4660 0.6740 0.1880], ...
                             'Displayname', 'A\oplusB')
    legend('Location','southeast','AutoUpdate','off')
    plot(poly1);
    plot(poly2);
    plot(poly3);

    scatter(xp1(1:idx), yp1(1:idx), 'filled', 'ko');
    scatter(xp2(1:idx), yp2(1:idx), 'filled','ko');
    scatter(xp3(1:idx), yp3(1:idx), 'filled' ,'ko');

    x1 = xp1(idx); y1 = yp1(idx);
    x2 = xp2(idx); y2 = yp2(idx);

    L1 = sqrt(x1^2 + y1^2);
    L2 = sqrt(x2^2 + y2^2);

    q = quiver(0,0,x1,y1, 'linewidth',4, 'color','w','Maxheadsize',1/L1,  ...
                             'autoscale','off');
    q = quiver(0,0,x1,y1, 'linewidth',3, 'color','k','Maxheadsize',1/L1,  ...
                             'autoscale','off');

    q = quiver(0, 0, x2, y2, 'linewidth',4, 'color','w','Maxheadsize',1/L2, ...
                              'autoscale','off');
    q = quiver(0, 0, x2, y2, 'linewidth',3, 'color','k','Maxheadsize',1/L2, ...
                              'autoscale','off');

    q = quiver(x1+x2, y1+y2, -x1, -y1, 'linewidth',4, 'color','w', ...
                             'Maxheadsize',0.5, 'autoscale','off',  ...
                             'showarrowhead','off','markersize',35);
    q.Marker = '.';
    q = quiver(x1+x2, y1+y2, -x2, -y2, 'linewidth',4, 'color','w', ...
                             'Maxheadsize',0.5, 'autoscale','off',  ...
                             'showarrowhead','off','markersize',35);
    q.Marker = '.';
    q = quiver(x1+x2, y1+y2, -x2, -y2, 'linewidth',3, 'color',[0.3,0.3,0.3], ...
                             'Maxheadsize',0.5, 'autoscale','off',  ...
                             'showarrowhead','off','markersize',30);
    q.Marker = '.';
    q = quiver(x1+x2, y1+y2, -x1, -y1, 'linewidth',3, 'color',[0.3,0.3,0.3], ...
                             'Maxheadsize',0.5, 'autoscale','off',  ...
                             'showarrowhead','off','markersize',30);
    q.Marker = '.';
    
    drawnow;
end

