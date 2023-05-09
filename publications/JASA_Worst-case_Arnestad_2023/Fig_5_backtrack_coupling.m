
close all;

array = biat.SensorArray;

rL = 2.6; rU = 3.1; phL = -0.2*pi; phU = -0.1*pi;

tol = 1e-3;

polar = ciat.PolarInterval(ciat.RealInterval(rL, rU), ciat.RealInterval(phL,phU));
circle = ciat.CircularInterval(2+1j, 0.6);

polP = ciat.PolygonalInterval(polar,'tolerance',tol);
cirP = ciat.PolygonalInterval(circle,'tolerance',tol);
prodP = (cirP .* polP);

fig=figure(1);
fig.Position=[100,100,550, 350];
hold on
axis equal

pgon = polyshape(real(polP.Points), imag(polP.Points));
plot(pgon,'linestyle', '-', 'linewidth', 1.5,'FaceColor',[0 0.4470 0.7410],'FaceAlpha',0.3, 'displayname', '\it E^I')

pgon = polyshape(real(cirP.Points), imag(cirP.Points));
plot(pgon,'linestyle', '--', 'linewidth', 1.5,'FaceColor',[0.9290 0.6940 0.1250],'FaceAlpha',0.3, 'displayname', '\it A^I')

% Plot product
pgon = polyshape(real(prodP.Points), imag(prodP.Points));
plot(pgon,'linestyle', '-.', 'linewidth', 1.5,'FaceColor',[0.8500 0.3250 0.0980],'FaceAlpha',0.3, 'displayname', '\it E^I \cdot A^I')


%% Backtrack point in product

% Select point to backtrack
z = 8.7680 + 0.6876i; %points(randi(length(points)));

% Plot E_inv
polar2 = ciat.PolarInterval( ciat.RealInterval(abs(z)/rU,abs(z)/rL), ciat.RealInterval(angle(z)-phU,angle(z)-phL) );
Pol2 = polar2.sample(10);
pgon = polyshape(real(Pol2), imag(Pol2));
plot(pgon,'linestyle', ':', 'linewidth', 2, 'FaceColor',[0.3, 0.3, 0.3],'FaceAlpha',0.6, 'displayname', '{\it E}_{inv.}^{\it I}')
plot(z, 'k^', 'markerfacecolor','r', 'displayname', "\it z", 'markersize', 13, 'linewidth',2)

% Find intersection of E_inv and A^I
polar2arc = polar2.Abs.Infimum * exp(1j*linspace(polar2.Angle.Infimum,...
                                                 polar2.Angle.Supremum,10 ));
temp = abs((cirP.Points - Pol2.'));
[m, argmin] = min(min(temp,[],1), [], 2);
Az = Pol2(argmin);

plot(Az,'k', 'LineStyle', 'none' ,'marker', 'o', 'markerfacecolor', [0.9290 0.6940 0.1250], 'displayname', "\it A", 'markersize', 13, 'linewidth',2)

plot(z/Az,'k', 'LineStyle', 'none' ,'marker', 'square','markerfacecolor', [0 0.4470 0.7410], 'displayname', "\it E = z / A", 'markersize', 13, 'linewidth',2)

polP = ciat.PolygonalInterval(polar,'tolerance',tol);
cirP = ciat.PolygonalInterval(circle,'tolerance',tol);
prodP = (cirP .* polP);
plot([prodP.Points(1:end-2); prodP.Points(1)],'-','HandleVisibility','off', 'color', [1,0.3,0.3,0.5]*0.7, 'linewidth', 1)
plot([prodP.Points(end-1:end); prodP.Points(end)],'-','HandleVisibility','off', 'color', [1,0.3,0.3, 0.5]*0.9, 'linewidth', 1)

%circle = ciat.CircularInterval(2+1j, 0.6);
polP = ciat.PolygonalInterval(polar,'tolerance',tol);
polP_inside = polP.Points .* (2+1j);
plot([polP_inside; polP_inside(1)],'-','HandleVisibility','off', 'color', [0,0,0, 0.2]*0.9, 'linewidth', 1)

legend('location', 'southwest');
xlabel('Real axis')
ylabel('Imaginary axis')
%yticks([0 : 0.2 : 1])
%xticks([0 : 0.2 : 1])
set(gca,'FontSize',18)
xlim([-1, 9.5])
legend boxoff
%exportgraphics(gca,'Fig5_backtrack_coupling.pdf')%,'Resolution',300) 