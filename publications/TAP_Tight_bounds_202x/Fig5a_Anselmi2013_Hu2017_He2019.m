clear
% close all

%% Anselmi 2013 and Hu 2017 (amplitude errors)

% Set parameters
conf = getFig5conf(1);
M = conf.M;
w = conf.w;
theta = conf.theta;
ampErr = conf.ampErr;
xL = conf.xL;
yL = conf.yL;

% Calculate nominal element phase 
phi = ((0:M-1)-(M-1)/2)' * pi*sin(theta);

% Calculate amplitude interval
A = w * (1 + ciat.RealInterval(-ampErr/2,ampErr/2));

% Calculate nominal beampattern and power pattern
AF_nom = w .* exp(1j*phi);
B_nom = sum(AF_nom);
P_nom = abs(B_nom)^2;

% Define and cast intervals
AF_p = ciat.PolarInterval(A,ciat.RealInterval(phi));
AF_r = ciat.RectangularInterval(AF_p);
AF_g = ciat.PolygonalInterval(AF_p);
AF_a = ciat.PolyarcularInterval(AF_p);
AF_x = ciat.PolyarxInterval(AF_p);

% Sum intervals
B_r = sum(AF_r);
B_g = sum(AF_g);
B_a = sum(AF_a);
B_x = sum(AF_x);

% Get power intervals
P_g = abs(B_g)^2;
P_a = abs(B_a)^2;
P_x = abs(B_x)^2;
P_r = abs(B_r)^2;

% Hu's Taylor based approximation
A_R = sum(A .* cos(phi));
A_I = sum(A .* sin(phi));
P_Hu_d = 2*cos(phi) .* A_R + 2*sin(phi) .* A_I;
P_Hu_0 = abs(sum(w .* exp(1j*phi)))^2;
P_Hu_U = P_Hu_0 + sum( abs(P_Hu_d.sup) .* A.width/2 );
P_Hu_I = P_Hu_0 - sum( abs(P_Hu_d.inf) .* A.width/2 );
P_Hu_I (P_Hu_I<0) = 0;
P_Hu = ciat.RealInterval(P_Hu_I,P_Hu_U);

% He's matrix method
a_mid = A.mid';
a_rad = A.width'/2;
Theta = cos(phi - phi');
P_He_mid = abs(sum(w .* exp(1j*phi)))^2;
P_He_sup = P_He_mid + 2 * abs(a_mid * Theta) * a_rad' + a_rad * abs(Theta)*a_rad';
P_He_inf = P_He_mid - 2 * abs(a_mid * Theta) * a_rad';
P_He_inf(P_He_inf<0) = 0;
P_He = ciat.RealInterval(P_He_inf,P_He_sup);

% Power values
B_abs = sqrt([P_r,P_g,P_x,P_Hu,P_He]);

%% Plot

% Set parameters
lineWidthL = 4;
lineWidthM = 3;
lineWidthS = 2;
cList = getColorList(1);
fBox = [0 2 B_x.Imag.inf B_x.Imag.sup];

figure(1);clf;
hold on;axis equal
set(gca,'DefaultLineLineWidth',3)
plot(0,0,'k+')
xlabel('Real')
ylabel('Imag')

% Plot operand intervals
AF_r.plot('c','linewidth',lineWidthL);
AF_g.plot('b','linewidth',lineWidthL);
AF_x.plot('r','linewidth',lineWidthS);

% Plot sum intervals
lC = B_x.plot('r','linewidth',lineWidthS,'DisplayName','Polyarcular');
lB = B_g.plot('b','linewidth',lineWidthL,'DisplayName','Polygonal');
lA = B_r.plot('color',cList(1,:),'linewidth',lineWidthL,'DisplayName','Rectangular');

% Interval label
for n = 1:M
    text(real(AF_nom(n))+0.01,imag(AF_nom(n)),['$E_{' num2str(n) '}^I$'], ...
                'HorizontalAlignment','center', 'Interpreter','latex')
end
text(B_r.real.mid,B_r.imag.mid,'$B^I$',...
                    'HorizontalAlignment','center', 'Interpreter','latex')
text(B_r.real.mid-0.05,B_r.imag.mid,'$\underline{|B^I|}$',...
                    'HorizontalAlignment','center', 'Interpreter','latex')
text(B_r.real.mid+0.05,B_r.imag.mid,'$\overline{|B^I|}$',...
                    'HorizontalAlignment','center', 'Interpreter','latex')

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

% Tightness values
annotText = {['\color{red}\tau^{(a)}=' num2str(100*P_a.Width ./ P_x.Width,3) '%'],...
            ['\color{blue}\tau^{(g)}=' num2str(100*P_a.Width ./ P_g.Width,3) '%'],...
            ['\color[rgb]{' num2str(cList(1,:)) '}\tau^{(r)}=' ...
                        num2str(100*P_a.Width ./ P_r.Width,3) '%'],...
            ['\color[rgb]{' num2str(cList(3,:)) '}\tau^{(T)}=' ...
                        num2str(100*P_a.Width ./ P_Hu.Width,3) '%'],...
            ['\color[rgb]{' num2str(cList(4,:)) '}\tau^{(M)}=' ...
                        num2str(100*P_a.Width ./ P_He.Width,3) '%']};
annotation('textbox',[0.14 0.66 .14 0.25],'String',annotText, ...
           'BackgroundColor','w','VerticalAlignment','top',...
           'HorizontalAlignment','center','FitBoxToText','on');


% Add zoom window for the operand interval
axes('position',[.14 .10 .09 .3]); hold on; box on
    % Plot operand intervals
AF_p.plot('b','linewidth',lineWidthL);
AF_r.plot('c','linewidth',lineWidthL);
AF_x.plot('r','linewidth',lineWidthS);
axis equal
xlim([0.064 0.073])
xticks(-0.5)
xtickangle(90)
set(gca, 'XAxisLocation', 'top')
ylim([0.115 0.128])
yticks(0.65)
set(gca, 'YAxisLocation', 'right')
n=6;
text(real(AF_nom(n)),imag(AF_nom(n))+0.001,['$E_{' num2str(n) '}^I$'], ...
                'HorizontalAlignment','right', 'Interpreter','latex')

% Add zoom window for the infimum
axes('position',[.63 .13 .07 .15]); hold on; box on
fimplicit(@(x,y) x.^2+y.^2-P_r.inf,'c-.','linewidth',lineWidthL)
fimplicit(@(x,y) x.^2 + y.^2 - P_Hu.inf,':','color',cList(3,:),...
                            'linewidth',lineWidthL)
fimplicit(@(x,y) x.^2 + y.^2 - P_He.inf,':','color',cList(4,:),...
                            'linewidth',lineWidthL)
fimplicit(@(x,y) x.^2+y.^2-(inf(abs(B_g)))^2,'b--','linewidth',lineWidthL)
fimplicit(@(x,y) x.^2+y.^2-P_x.inf,fBox,'r:','linewidth',lineWidthM);
axis equal
xlim([min(B_abs.inf)-1e-3 , max(B_abs.inf)+1e-3])
xticks(1.2)
xtickangle(90)
set(gca, 'XAxisLocation', 'top')
yticks(0)
text(mean(B_abs.inf),0,'$\underline{|B^I|}$',...
            'HorizontalAlignment','center', 'Interpreter','latex')


% Add zoom window for the supremum
axes('position',[.76 .13 .07 .15]); hold on; box on
fimplicit(@(x,y) x.^2+y.^2-P_r.sup,'c-.','linewidth',lineWidthL)
fimplicit(@(x,y) x.^2 + y.^2 - P_Hu.sup,':','color',cList(3,:),...
                            'linewidth',lineWidthL)
fimplicit(@(x,y) x.^2 + y.^2 - P_He.sup,':','color',cList(4,:),...
                            'linewidth',lineWidthL)
fimplicit(@(x,y) x.^2+y.^2-P_g.sup,'b--','linewidth',lineWidthL)
fimplicit(@(x,y) x.^2+y.^2-P_x.sup,fBox,'r:','linewidth',lineWidthM);
axis equal
xlim([min(B_abs.sup)-1e-3 , max(B_abs.sup)+1e-3])
xticks(1.9)
set(gca, 'XAxisLocation', 'top')
xtickangle(90)
yticks(0)
text(mean(B_abs.sup),0,'$\overline{|B^I|}$',...
        'HorizontalAlignment','center', 'Interpreter','latex')


fontsize(30,'point')


