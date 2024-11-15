clear
% close all

%% Anselmi 2013 and Hu 2017 (amplitude errors)

% Set parameters
conf = getFig5conf(1);
M = conf.M;
w = conf.w;
theta = conf.theta;
ampErr = conf.ampErr;

% Calculate nominal element phase 
phi = ((0:M-1)-(M-1)/2)' * pi*sin(theta);

% Calculate amplitude interval
A = w * (1 + ciat.RealInterval(-ampErr/2,ampErr/2));

% Calculate nominal beampattern and power pattern
E_nom = w .* exp(1j*phi);
B_nom = sum(E_nom);
P_nom = abs(B_nom)^2;

% Define and cast intervals
E_p = ciat.PolarInterval(A,ciat.RealInterval(phi));
tic;E_g = ciat.PolygonalInterval(E_p,'tolerance',conf.tol);T_g(1) = toc;
tic;E_x = ciat.PolyarxInterval(E_p);T_x(1) = toc;
tic; E_r = ciat.RectangularInterval(E_g); T_r(1) = toc;
E_a = ciat.PolyarcularInterval(E_p);

% Sum intervals
tic;B_r = sum(E_r);T_r(2)=toc;
tic;B_g = sum(E_g);T_g(2)=toc;
tic;B_x = sum(E_x);T_x(2)=toc;
B_a = sum(E_a);

% Get power intervals
tic;P_g = abs(B_g)^2;T_g(3)=toc;
tic;P_x = abs(B_x)^2;T_x(3)=toc;
tic;P_r = abs(B_r)^2;T_r(3)=toc;
P_a = abs(B_a)^2;

% Hu's Taylor based approximation
tic;
A_R = sum(A .* cos(phi));
A_I = sum(A .* sin(phi));
P_Hu_d = 2*cos(phi) .* A_R + 2*sin(phi) .* A_I;
P_Hu_0 = abs(sum(w .* exp(1j*phi)))^2;
P_Hu_U = P_Hu_0 + sum( abs(P_Hu_d.sup) .* A.width/2 );
P_Hu_I = P_Hu_0 - sum( abs(P_Hu_d.inf) .* A.width/2 );
P_Hu_I (P_Hu_I<0) = 0;
P_Hu = ciat.RealInterval(P_Hu_I,P_Hu_U);
T_Hu(3) = toc;

% He's matrix method
tic;
a_mid = A.mid';
a_rad = A.width'/2;
Theta = cos(phi - phi');
P_He_mid = abs(sum(w .* exp(1j*phi)))^2;
P_He_sup = P_He_mid + 2 * abs(a_mid * Theta) * a_rad' + a_rad * abs(Theta)*a_rad';
P_He_inf = P_He_mid - 2 * abs(a_mid * Theta) * a_rad';
P_He_inf(P_He_inf<0) = 0;
P_He = ciat.RealInterval(P_He_inf,P_He_sup);
T_He(3) = toc;

% Calculate tightness
tau_x = P_a.Width ./ P_x.Width;
tau_g = P_a.Width ./ P_g.Width;
tau_r = P_a.Width ./ P_r.Width;
tau_Hu = P_a.Width ./ (P_a.Width + abs(P_Hu.sup-P_x.sup) + ...
                                   abs(P_Hu.inf-P_x.inf));
tau_He = P_a.Width ./ P_He.Width;

%% Plot

% Set parameters
lineWidthL = 4;
lineWidthM = 3;
lineWidthS = 2;
cList = getColorList(conf.cID);
fBox = [0 2 B_x.Imag.inf B_x.Imag.sup];
xL = conf.xL;
yL = conf.yL;

% Get beampattern interval absolute values
B_abs = sqrt([P_r,P_g,P_x,P_Hu,P_He]);

figure(1);clf;
hold on;axis equal
set(gca,'DefaultLineLineWidth',3)
plot(0,0,'k+')
xlabel('Real')
ylabel('Imag')

% Plot operand intervals
E_g.plot('b','linewidth',lineWidthS);
E_x.plot('r','linewidth',lineWidthS);
E_r.plot('color',cList(1,:),'linewidth',lineWidthS);

% Plot sum intervals
lC = B_x.plot('r','linewidth',lineWidthS,'DisplayName','Polyarcular');
lB = B_g.plot('b','linewidth',lineWidthL,'DisplayName','Polygonal');
lA = B_r.plot('color',cList(1,:),'linewidth',lineWidthL,'DisplayName','Rectangular');

% Interval label
for m = 1:M
    text(real(E_nom(m))+0.015,imag(E_nom(m)),['$E_{' num2str(m) '}^I$'], ...
                'HorizontalAlignment','center', 'Interpreter','latex')
end
text(B_r.real.mid,B_r.imag.mid,'$B^I$',...
                    'HorizontalAlignment','center', 'Interpreter','latex')
text(B_r.real.inf-0.01,B_r.imag.mid,'$\underline{|B^I|}$',...
                    'HorizontalAlignment','right', 'Interpreter','latex')
text(B_r.real.sup+0.01,B_r.imag.mid,'$\overline{|B^I|}$',...
                    'HorizontalAlignment','left', 'Interpreter','latex')

% Anselmi
l3 = fimplicit(@(x,y) x.^2+y.^2-P_r.inf,fBox,'-.',...
            'color',cList(1,:),'linewidth',lineWidthM,'DisplayName','Anselmi');
fimplicit(@(x,y) x.^2+y.^2-P_r.sup,fBox,'-.',...
            'color',cList(1,:),'linewidth',lineWidthM)
% Hu
l2 = fimplicit(@(x,y) x.^2 + y.^2 - P_Hu.inf,fBox,':', ...
            'color',cList(3,:),'linewidth',lineWidthM,'DisplayName','Hu');
fimplicit(@(x,y) x.^2 + y.^2 - P_Hu.sup,fBox,':', ...
            'color',cList(3,:), 'linewidth',lineWidthM)
%He
l1 = fimplicit(@(x,y) x.^2 + y.^2 - P_He.inf,fBox,':',...
            'color',cList(4,:),'linewidth',lineWidthM,'DisplayName','He');
fimplicit(@(x,y) x.^2 + y.^2 - P_He.sup,fBox,':',...
            'color',cList(4,:),'linewidth',lineWidthM)
% Polygonal
l4 = fimplicit(@(x,y) x.^2+y.^2-P_g.sup,fBox,'b--','linewidth',lineWidthM,...
                                            'DisplayName','Tenuti');
fimplicit(@(x,y) x.^2+y.^2-P_g.inf,fBox,'b--','linewidth',lineWidthM)

% Set figure limits
xlim(xL); ylim(yL)


% % Add two legend windows
% leg1 = legend([lA(1),lB(1),lC(1)],'Location','NorthWest');
% title(leg1,'$E^I \mathrm{ wrappers}$','Interpreter','latex');
% ax2 = axes('position',get(gca,'position'),'visible','off');
% leg2 = legend(ax2, [l1,l2,l3,l4,l5], 'Location','northeast');
% title(leg2,'$|B^I| \mathrm{ methods}$','Interpreter','latex');

% Add zoom window for the operand interval
axes('position',[0.09,0.18,0.3,0.3]); hold on; box on
E_p.plot('b','linewidth',lineWidthL);
E_r.plot('color',cList(1,:),'linewidth',lineWidthL);
E_x.plot('r','linewidth',lineWidthS);
axis equal
set(gca, 'XAxisLocation', 'top')
m = 6;
maxWidth = max([E_r(m).Real.Width , E_r(m).Imag.Width]);
xlim(E_r(m).Real.Midpoint + [-1 1]*maxWidth/2)
ylim(E_r(m).Imag.Midpoint + [-1 1]*maxWidth/2)
xticks([]);
yticks([]);
set(gca, 'YAxisLocation', 'right')
text(real(E_nom(m)),imag(E_nom(m))+0.001,['$E_{' num2str(m) '}^I$'], ...
                'HorizontalAlignment','right', 'Interpreter','latex')

% Add zoom window for the infimum
axes('position',[0.73,0.178,0.14,0.3]); hold on; box on
fimplicit(@(x,y) x.^2+y.^2-P_r.inf,'-.', ...
                        'color',cList(1,:),'linewidth',lineWidthL)
fimplicit(@(x,y) x.^2 + y.^2 - P_Hu.inf,':', ...
                        'color',cList(3,:),'linewidth',lineWidthL)
fimplicit(@(x,y) x.^2 + y.^2 - P_He.inf,':', ...
                        'color',cList(4,:),'linewidth',lineWidthL)
fimplicit(@(x,y) x.^2+y.^2-(inf(abs(B_g)))^2,'b--','linewidth',lineWidthL)
fimplicit(@(x,y) x.^2+y.^2-P_x.inf,fBox,'r:','linewidth',lineWidthM);
axis equal
xlim([min(B_abs.inf)-1e-3 , max(B_abs.inf)+1e-3])
xticks([])
xtickangle(90)
set(gca, 'XAxisLocation', 'top')
yticks([])
text(mean(B_abs.inf),0,'$\underline{|B^I|}$',...
            'HorizontalAlignment','center', 'Interpreter','latex')


% Add zoom window for the supremum
axes('position',[0.730,0.62,0.14,0.30]); hold on; box on
fimplicit(@(x,y) x.^2+y.^2-P_r.sup,'-.', ...
                    'color',cList(1,:),'linewidth',lineWidthL)
fimplicit(@(x,y) x.^2 + y.^2 - P_Hu.sup,':',...
                    'color',cList(3,:),'linewidth',lineWidthL)
fimplicit(@(x,y) x.^2 + y.^2 - P_He.sup,':',...
                    'color',cList(4,:),'linewidth',lineWidthL)
fimplicit(@(x,y) x.^2+y.^2-P_g.sup,'b--','linewidth',lineWidthL)
fimplicit(@(x,y) x.^2+y.^2-P_x.sup,fBox,'r:','linewidth',lineWidthM);
axis equal
xlim([min(B_abs.sup)-1e-3 , max(B_abs.sup)+1e-3])
xticks([])
set(gca, 'XAxisLocation', 'top')
xtickangle(90)
yticks([])
text(mean(B_abs.sup),0,'$\overline{|B^I|}$',...
        'HorizontalAlignment','center', 'Interpreter','latex')

fontsize(40,'point')


% Tightness values
annotText = {join(['\color{red}\tau^{(a)}=' num2str(100*tau_x,3) '%' ...
                    ', T^{(a)}=' join(compose("%0.0f",1e3*T_x),'+') 'ms'],''),...
            join(['\color{blue}\tau^{(g)}=' num2str(100*tau_g,3) '%' ...
                    ', T^{(g)}=' join(compose("%0.0f",1e3*T_g),'+') 'ms'],''),...
            join(['\color[rgb]{' num2str(cList(1,:)) '}\tau^{(r)}=' ...
                        num2str(100*tau_r,3) '%' ...
                    ', T^{(r)}=' join(compose("%0.0f",1e3*T_r),'+') 'ms'],''),...
            join(['\color[rgb]{' num2str(cList(3,:)) '}\tau^{(T)}=' ...
                        num2str(100*tau_He,3) '%' ...
                    ', T^{(T)}=' join(compose("%0.0f",1e3*T_Hu),'+') 'ms'],''),...
            join(['\color[rgb]{' num2str(cList(4,:)) '}\tau^{(M)}=' ...
                        num2str(100*tau_Hu,3) '%' ...
                    ', T^{(M)}=' join(compose("%0.0f",1e3*T_He),'+') 'ms'],'')};
annotation('textbox',[0.164,0.553,0.280 0.367],'String',annotText, ...
           'BackgroundColor','w','VerticalAlignment','top','fontSize',30,...
           'HorizontalAlignment','center','FitBoxToText','on');
