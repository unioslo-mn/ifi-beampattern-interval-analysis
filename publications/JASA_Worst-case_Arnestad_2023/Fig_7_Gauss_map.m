clear

% Make points
points1 = [1, 1.3, 1.5] + 1j*[1, 1.4, 1] + 0.5 + 0.1j;
points2 = [1, 1.5, 1, 1.5] + 1j*[0, 0, 0.7,1]-1j*0.2;
points2 = points2*exp(0.1*2*pi*1j); 

poly1 = ciat.PolygonalInterval( points1.');
poly2 = ciat.PolygonalInterval( points2.');
poly3 = poly1+poly2;

% Find maximum point on the sum polygon
[~, P_U_idx] = max(abs(poly3.Points));
P_U = poly3.Points(P_U_idx);
P_U_ang = wrapTo2Pi( angle( poly3.Points(P_U_idx)));

scale = (abs(P_U)+0.5)/abs(P_U);
[x_dir,y_dir] = biat.SensorArray.getArrowLine( ...
                            [0,0], scale*[real(P_U), imag(P_U)], 0.1);
circ = ciat.CircularInterval(0,abs(P_U));

% Backtrack maximum point to the other two intervals
[select_z] = backtrack([poly1,poly2,poly3], P_U_ang);

new_scale = (0.5)/abs(P_U);
start = [real(select_z(1)), imag(select_z(1))];
[x1,y1] = biat.SensorArray.getArrowLine( ...
                    start, start + [real(P_U), imag(P_U)]*new_scale, 0.1);

start = [real(select_z(2)), imag(select_z(2))];
[x2,y2] = biat.SensorArray.getArrowLine( ...
                    start, start + [real(P_U), imag(P_U)]*new_scale, 0.1);


fig=figure(2);
fig.Position=[100,100,600, 400];
set(gca,'LineWidth',2)
set(gca,'DefaultLineLineWidth',2)
set(gca,'defaultAxesFontSize',20)
set(gca, 'fontsize', 15)
hold on;


% Good code below, but bloated copy/paste to plot outer angles for Gauss map
v = poly1.Points;
prev_angle = wrapTo2Pi( angle( v(1) - v(end) ) - pi/2); % init prev angle
v = [v; v(end);v(1)];
for j = 1 : length(v)-1
    current_angle = wrapTo2Pi( angle( v(j+1) - v(j) )- pi/2);
    
    start = v(j);
    stop = (0.2*(v(j+1) - v(j))/abs(v(j+1) - v(j))  *exp(-1j*pi/2)) + v(j);
    start_vec = [real(start), imag(start)];
    stop_vec = [real(stop), imag(stop)];
    [x,y] = biat.SensorArray.getArrowLine( start_vec, stop_vec, 0.1);
    plot(x,y,'k','HandleVisibility','off');
    
    start = v(j+1);
    stop = (0.2*(v(j+1) - v(j))/abs(v(j+1) - v(j))  *exp(-1j*pi/2)) + v(j+1);
    start_vec = [real(start), imag(start)];
    stop_vec = [real(stop), imag(stop)];
    [x,y] = biat.SensorArray.getArrowLine( start_vec, stop_vec, 0.1);
    plot(x,y, 'k','HandleVisibility','off');    
    
    prev_angle = current_angle;
end

v = poly2.Points;
prev_angle = wrapTo2Pi( angle( v(1) - v(end) ) - pi/2); % init prev angle
v = [v; v(end);v(1)];
for j = 1 : length(v)-1
    current_angle = wrapTo2Pi( angle( v(j+1) - v(j) )- pi/2);
    
    start = v(j);
    stop = (0.2*(v(j+1) - v(j))/abs(v(j+1) - v(j))  *exp(-1j*pi/2)) + v(j);
    start_vec = [real(start), imag(start)];
    stop_vec = [real(stop), imag(stop)];
    [x,y] = biat.SensorArray.getArrowLine( start_vec, stop_vec, 0.1);
    plot(x,y,'k','HandleVisibility','off');
    
    start = v(j+1);
    stop = (0.2*(v(j+1) - v(j))/abs(v(j+1) - v(j))  *exp(-1j*pi/2)) + v(j+1);
    start_vec = [real(start), imag(start)];
    stop_vec = [real(stop), imag(stop)];
    [x,y] = biat.SensorArray.getArrowLine( start_vec, stop_vec, 0.1);
    plot(x,y, 'k','HandleVisibility','off');    
    
    prev_angle = current_angle;
end

v = poly3.Points;
prev_angle = wrapTo2Pi( angle( v(1) - v(end) ) - pi/2); % init prev angle
v = [v; v(end);v(1)];
for j = 1 : length(v)-1
    current_angle = wrapTo2Pi( angle( v(j+1) - v(j) )- pi/2);
    
    start = v(j);
    stop = (0.2*(v(j+1) - v(j))/abs(v(j+1) - v(j))  *exp(-1j*pi/2)) + v(j);
    start_vec = [real(start), imag(start)];
    stop_vec = [real(stop), imag(stop)];
    [x,y] = biat.SensorArray.getArrowLine( start_vec, stop_vec, 0.1);
    plot(x,y,'k','HandleVisibility','off');
    
    start = v(j+1);
    stop = (0.2*(v(j+1) - v(j))/abs(v(j+1) - v(j))  *exp(-1j*pi/2)) + v(j+1);
    start_vec = [real(start), imag(start)];
    stop_vec = [real(stop), imag(stop)];
    [x,y] = biat.SensorArray.getArrowLine( start_vec, stop_vec, 0.1);
    plot(x,y, 'k','HandleVisibility','off');    
    
    prev_angle = current_angle;
end

% plot remaining
pgon = polyshape(real(poly1.Points), imag(poly1.Points));
plot(pgon, 'linewidth', 2, 'edgecolor',[0 0.4470 0.7410], 'FaceColor',[0 0.4470 0.7410],'FaceAlpha',0.1, 'Displayname', 'Poly_1')
plot(poly1, '-o', 'color', [0 0.4470 0.7410]*0.6, 'MarkerFaceColor', [0 0.4470 0.7410],'HandleVisibility','off');

pgon = polyshape(real(poly2.Points), imag(poly2.Points));
plot(pgon, 'linewidth', 2, 'edgecolor',[0.4660 0.6740 0.1880], 'FaceColor',[0.4660 0.6740 0.1880]*0.6,'FaceAlpha',0.3, 'Displayname', 'Poly_2')
plot(poly2,'-o', 'color', [0.4660 0.6740 0.1880], 'MarkerFaceColor', [0.4660 0.6740 0.1880],'HandleVisibility','off');

pgon = polyshape(real(poly3.Points), imag(poly3.Points));
plot(pgon, 'linewidth', 2, 'edgecolor',[0.9290 0.6940 0.1250], 'FaceColor',[0.9290 0.6940 0.1250],'FaceAlpha',0.1, 'Displayname', 'Poly_1+Poly_2')
plot(poly3,'o', 'color', [0.9290 0.6940 0.1250], 'MarkerFaceColor', [0.9290 0.6940 0.1250],'HandleVisibility','off');

plot(x1,y1,'r', 'linewidth', 3,'HandleVisibility','off');
plot(x2,y2,'r', 'linewidth', 3,'HandleVisibility','off');
plot(x_dir,y_dir,'r', 'linewidth', 3,'HandleVisibility','off');

plot(circ, 'k:','HandleVisibility','off');
plot(select_z, 'k^','markersize',12,'MarkerFaceColor',[1 0 0],'Displayname', 'Points \it z_c');
plot(P_U, 'k^','markersize',14,'MarkerFaceColor',[1 0 0]*0.5, 'Displayname', 'Sum \it z_c') 

axis equal; xlim([0,3.5]); ylim([0,3.5]);
xticks([0: 1 : 3.5])
yticks([0: 1 : 3.5])
xlabel('Real axis'); ylabel('Imaginary axis');

quiver(0,0,0,0, 'displayname', 'Extr. direction', 'color','r', 'linewidth', 3)

h = legend('Location','northwest','AutoUpdate','off'); 
legend boxoff;
ax = gcf;
%exportgraphics(ax,'Fig7_Backtrack_sum.pdf')