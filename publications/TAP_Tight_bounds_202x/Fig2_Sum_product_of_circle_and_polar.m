clear
% close all

%% Calculate polar circular product

% Set intervals
A = ciat.CircularInterval(1+0.5i,0.10);
E = ciat.PolarInterval(1.6 , 2.0 , 0.05*pi, 0.10*pi);

% Wrap intervals
Ar = ciat.RectangularInterval(A);
Er = ciat.RectangularInterval(E);
Ec = ciat.CircularInterval(E);

% Calculate product intervals
EAg = ciat.PolygonalInterval(E,A,'tolerance',0.1);
EAa = ciat.PolyarcularInterval(E,A);
EAx = ciat.PolyarxInterval(E,A);
EAr = Ar * Er;
EAc = A * Ec;

% Calculate corner circles
Ap(4,1) = ciat.CircularInterval;
Ap(1) = A * E.Abs.inf*exp(1j*E.Angle.inf);
Ap(2) = A * E.Abs.sup*exp(1j*E.Angle.inf);
Ap(3) = A * E.Abs.inf*exp(1j*E.Angle.sup);
Ap(4) = A * E.Abs.sup*exp(1j*E.Angle.sup);

% Calculate rotated polar
Ep = E * A.Center;

% Sample operands and calculate Minkowski product
Esmp = E.sample(30);
Asmp = A.Center + A.Radius * exp(1j*linspace(-pi,pi,30));
EAgsmp = Esmp * Asmp;

%% Plot polar circular product

% Set parameters
cols = [200 165 0]/255;
fontSizeDef = 20;
fontSizeS = 17;
fontSizeL = 25;

% Intialize figure
figure(1);clf;hold on; axis equal

% Polar intervals
p1 = E.plot('m--','LineWidth',2,'DisplayName','Polar');
Ecorn = E.abs.Bounds  * exp(1j*E.angle.Bounds');
plot(real(Ecorn),imag(Ecorn),'mo')
plot(real(Ecorn*A.Center),imag(Ecorn*A.Center),'mo')
Ep.plot('m--','LineWidth',2);

% Circular intervals
p2 = A.plot('--','LineWidth',2,'color',cols(1,:), ...
                'DisplayName','Circular');
plot(real(A.Center),imag(A.Center),'+','color',cols(1,:))
plot(real(A.Center)+[0,sqrt(1/2)*A.Radius] ,...
     imag(A.Center)+[0,sqrt(1/2)*A.Radius], ...
                    '-','color',cols(1,:),'LineWidth',1)
Ap.plot('--','LineWidth',2,'color',cols(1,:));

% Wrapped intervals
Ar.plot('c--','LineWidth',2);
Er.plot('c--','LineWidth',2);
Ec.plot('--','LineWidth',2,'Color',cols(1,:));

% Product intervals
p3 = EAg.plot('b-','LineWidth',2,'DisplayName','Convex polygon');
p4 = EAa.plot('g-','LineWidth',2,'DisplayName','Concave polyarc');
p5 = EAx.plot('r-','LineWidth',2,'DisplayName','Convex polyarc');
p6 = EAr.plot('c--','LineWidth',2,'DisplayName','Rectangular');
EAc.plot('--','LineWidth',2,'color',cols(1,:));

% Circular interval labels
text(real(A.Center),imag(A.Center)-0.05, '$A^I_m$',...
    'fontsize',fontSizeL,'HorizontalAlignment','center','Interpreter','latex')
text(real(A.Center)-0.05,imag(A.Center), '$p_m$','color',cols(1,:),...
    'fontsize',fontSizeS,'HorizontalAlignment','center','Interpreter','latex')
text(real(A.Center),imag(A.Center)+0.05, '$R_m$','color',cols(1,:),...
    'fontsize',fontSizeS,'HorizontalAlignment','center','Interpreter','latex')
text(real(Ap(1).Center),imag(Ap(1).Center)-0.05, ...
    '$A_m^I \underline{a}_m e^{j\underline{\varphi}_m}$','color',cols(1,:),...
    'fontsize',fontSizeS,'HorizontalAlignment','center','Interpreter','latex')
text(real(Ap(2).Center)+0.07,imag(Ap(2).Center)-0.06, ...
    '$A_m^I \overline{a}_m e^{j\underline{\varphi}_m}$','color',cols(1,:),...
    'fontsize',fontSizeS,'HorizontalAlignment','center','Interpreter','latex')
text(real(Ap(3).Center)-0.03,imag(Ap(3).Center)+0.07, ...
    '$A_m^I \underline{a}_m e^{j\overline{\varphi}_m}$','color',cols(1,:),...
    'fontsize',fontSizeS,'HorizontalAlignment','center','Interpreter','latex')
text(real(Ap(4).Center),imag(Ap(4).Center)+0.05, ...
    '$A_m^I \overline{a}_m e^{j\overline{\varphi}_m}$','color',cols(1,:),...
    'fontsize',fontSizeS,'HorizontalAlignment','center','Interpreter','latex')


% Polar interval labels
pnt = E.Abs.mid * exp(1j*E.Angle.mid);
text(real(pnt),imag(pnt), '$E^I_m$', ...
    'fontsize',fontSizeL,'HorizontalAlignment','center','Interpreter','latex')
text(real(pnt)-0.02,imag(pnt)+0.1, '$a^I_m$', 'color','m', ...
    'fontsize',fontSizeS,'HorizontalAlignment','center','Interpreter','latex')
text(real(pnt)+0.15,imag(pnt)+0.02, '$\varphi^I_m$', 'color','m', ...
    'fontsize',fontSizeS,'HorizontalAlignment','center','Interpreter','latex')
text(real(Ecorn(1))-0.05,imag(Ecorn(1))-0.05, ...
    '$\underline{a}_m e^{j\underline{\varphi}_m}$', 'color','m',...
    'fontsize',fontSizeS,'HorizontalAlignment','center','Interpreter','latex')
text(real(Ecorn(2))+0.10,imag(Ecorn(2))-0.03, ...
    '$\overline{a}_m e^{j\underline{\varphi}_m}$', 'color','m',...
    'fontsize',fontSizeS,'HorizontalAlignment','center','Interpreter','latex')
text(real(Ecorn(3))-0.08,imag(Ecorn(3)), ...
    '$\underline{a}_m e^{j\overline{\varphi}_m}$', 'color','m',...
    'fontsize',fontSizeS,'HorizontalAlignment','center','Interpreter','latex')
text(real(Ecorn(4))+0.07,imag(Ecorn(4))+0.04, ...
    '$\overline{a}_m e^{j\overline{\varphi}_m}$', 'color','m',...
    'fontsize',fontSizeS,'HorizontalAlignment','center','Interpreter','latex')
pnt = Ep.Abs.mid * exp(1j*Ep.Angle.mid);
text(real(pnt),imag(pnt), '$E^I_m p_m$', 'color','m',...
    'fontsize',fontSizeS,'HorizontalAlignment','center','Interpreter','latex')

% Product interval labels
text(1.95,0.90, '$E^I_m A^I_m$', ...
    'fontsize',fontSizeL,'HorizontalAlignment','center','Interpreter','latex')

% Tightness values
annotText = {['$\tau_{EA}^{(rect)}=' num2str(100*EAa.Area ./ EAr.Area,3),'\%$'],...
             ['$\tau_{EA}^{(circ)}=' num2str(100*EAa.Area ./ EAc.Area,3),'\%$'],...
             ['$\tau_{EA}^{(pgon)}=' num2str(100*EAa.Area ./ EAg.Area,3),'\%$'],...
             ['$\tau_{EA}^{(parc)}=' num2str(100*EAa.Area ./ EAx.Area,3),'\%$']};
annotation('textbox',[0.16 0.50 0.20 0.17],'String',annotText, ...
           'FontSize',fontSizeDef,...
           'BackgroundColor','w',...
           'HorizontalAlignment','right', ...
           'Interpreter','latex');

% Misc
xlim([0.55,2.15])
legend([p1(1),p2(1),p6(1),p3(1),p5(1),p4(1)],'location','northwest', ...
        'FontSize',fontSizeDef)

%% Calculate polar circular sum-product

% Set intervals
% A = ciat.CircularInterval(1+0.5i,0.10);
% E = ciat.PolarInterval(1.6 , 2.0 , 0.05*pi, 0.10*pi);
A(2) = ciat.CircularInterval(1-0.5i,0.10);
E(2) = ciat.PolarInterval(1.6 , 2.0 , -0.05*pi, -0.10*pi);

% Wrap intervals in various types
Er = ciat.RectangularInterval(E); 
Ec = ciat.CircularInterval(E);    
Eg = ciat.PolygonalInterval(E,'Tolerance',0.1); 
Ex = ciat.PolyarxInterval(E);
Ea = ciat.PolyarcularInterval(E);
Ar = ciat.RectangularInterval(A); 
Ac = A;    

% Multiply E and A
EAr = Er .* Ar; 
EAc = Ec .* Ac; 
EAg = ciat.PolygonalInterval(E,A,'Tolerance',0.1); 
EAx = ciat.PolyarxInterval(E,A); 
EAa = ciat.PolyarcularInterval(E,A); 

% Sum the EA intervals
Br = sum(EAr);
Bc = sum(EAc);
Bg = sum(EAg);
Bx = sum(EAx);
Ba = sum(EAa);

% Set parameters
Nsmp = 100;

% Sample E and A intervals
Es = E.sample(Nsmp);
As = A.sample(Nsmp);

% Multiply E and A intervals
EAs{1} = Es{1} * As{1}.';
EAs{2} = Es{2} * As{2}.';

% Sum EA intervals
Bs = EAs{1}(:) + EAs{2}(:).';

% Get convex hulls
idx = convhull(real(EAs{1}),imag(EAs{1}));
EAs{1} = EAs{1}(idx);
idx = convhull(real(EAs{2}),imag(EAs{2}));
EAs{2} = EAs{2}(idx);
% idx = convhull(real(Bs),imag(Bs));
% Bs = Bs(idx);


%% Plot polar circular sum-product
figure(2);clf;hold on;axis equal
set(0,'DefaultLineLineWidth',2)

% Plot wraps of E intervals
p1 = E.plot('m--','DisplayName','Polar');
p2 = Er.plot('c--','DisplayName','Rectangular');
p3 = Ec.plot('--','color',cols(1,:),'DisplayName','Circular');
% p3 = Eg.plot('b','DisplayName','Convex polygon');
% p4 = Ex.plot('r','lineWidth',2,'DisplayName','Convex polyarc');
% p5 = Ea.plot('g--','lineWidth',2,'DisplayName','Concave polyarc');

% Plot wraps of A intervals
Ar.plot('c--');
Ac.plot('--','color',cols(1,:));

% Plot origin
plot(0,0,'k+')

% Plot EA intervals
EAr.plot('c--');
EAc.plot('--','color',cols(1,:));
p4 = EAg.plot('b-','DisplayName','Convex polygon');
p5 = EAx.plot('r-','lineWidth',2,'DisplayName','Convex polyarc');
p6 = EAa.plot('g--','lineWidth',2,'DisplayName','Concave polyarc');
% plot(real(EAs{1}),imag(EAs{1}),'k--')
% plot(real(EAs{2}),imag(EAs{2}),'k--')

% Plot B interval
Br.plot('c--');
Bc.plot('--','color',cols(1,:));
Bg.plot('b');
Bx.plot('r','lineWidth',2);
Ba.plot('g--','lineWidth',2);

% Interval labels
text(real(Ec(1).Center),imag(Ec(1).Center),'$E^I_1$', ...
    'HorizontalAlignment', 'center','Interpreter','latex')
text(real(Ec(2).Center),imag(Ec(2).Center),'$E^I_2$', ...
    'HorizontalAlignment', 'center','Interpreter','latex')
text(real(EAc(1).Center),imag(EAc(1).Center),'$E^I_1 A^I_1$', ...
    'HorizontalAlignment', 'center','Interpreter','latex')
text(real(EAc(2).Center),imag(EAc(2).Center),'$E^I_2 A^I_2$', ...
    'HorizontalAlignment', 'center','Interpreter','latex')
text(real(Ac(1).Center),imag(Ac(1).Center),'$A^I_1$', ...
    'HorizontalAlignment', 'center','Interpreter','latex')
text(real(Ac(2).Center),imag(Ac(2).Center),'$A^I_2$', ...
    'HorizontalAlignment', 'center','Interpreter','latex')
text(real(Bc.Center),imag(Bc.Center)+0.1,'$B^I=E^I_1 A^I_1+E^I_2 A^I_2$', ...
    'HorizontalAlignment', 'center','Interpreter','latex')

% Tightness values
annotText = {['$\tau_B^{(rect)}=' num2str(100*Ba.Area ./ Br.Area,4),'\%$'],...
             ['$\tau_B^{(circ)}=' num2str(100*Ba.Area ./ Bc.Area,4),'\%$'],...
             ['$\tau_B^{(pgon)}=' num2str(100*Ba.Area ./ Bg.Area,4),'\%$'],...
             ['$\tau_B^{(parc)}=' num2str(100*Ba.Area ./ Bx.Area,4),'\%$']};
annotation('textbox',[0.65 0.75 0.21 0.17],'String',annotText, ...
           'HorizontalAlignment','right', ...
           'Interpreter','latex');


% Plot axes
ylim([-2.5,2.5])
xL = xlim;
yL = ylim;
plot(xL,[0,0],'color',0.7*ones(1,3),'linewidth',0.1)
plot([0,0],yL,'color',0.7*ones(1,3),'linewidth',0.1)

% Misc
fontsize(20,'points')
legend([p1(1),p2(1),p3(1),p4(1),p5(1),p6(1)],'location','southeast')
