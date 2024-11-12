clear
% close all

%% Zhang 2017 (amplitude and phase errors)
M = 8;
w = taylorwin(M,3,-20);w = w / sum(w);
theta = deg2rad(40);

% Set errors
ampErr = 0.1;
phaErr = deg2rad(5);

% Calculate nominal element phase 
phi = ((0:M-1)-(M-1)/2)' * pi*sin(theta);

% Calculate amplitude interval
A = w * (1 + ciat.RealInterval(-ampErr/2,ampErr/2));

% Define and cast intervals
AF_nom = w .* exp(1i*phi);
AF_p = ciat.PolarInterval(w * (1 + ciat.RealInterval(-ampErr/2,ampErr/2)),...
                          ciat.RealInterval(phi + [-1 1]*phaErr/2));
AF_r = ciat.RectangularInterval(AF_p);
AF_g = ciat.PolygonalInterval(AF_p,'tolerance',1e-3);
AF_x = ciat.PolyarxInterval(AF_p);
AF_a = ciat.PolyarcularInterval(AF_p);

% Sum intervals
B_r = sum(AF_r);
B_g = sum(AF_g);
B_x = sum(AF_x);
B_a = sum(AF_a);

% Power intervals
P_r = abs(B_r).^2;
P_g = abs(B_g).^2;
P_x = abs(B_x).^2;
P_a = abs(B_a).^2;

%% Plot

% Set parameters
lineWidthL = 4;
lineWidthM = 3;
lineWidthS = 2; 
fBox = [0 3 B_r.Imag.inf B_r.Imag.sup];

% Plot
figure(1);clf;hold on;axis equal;
set(gca,'DefaultLineLineWidth',lineWidthL)
plot(0,0,'k+')

% Plot operand intervals
AF_r.plot('c','linewidth',lineWidthL);
AF_g.plot('b','linewidth',lineWidthL);
AF_x.plot('r','linewidth',lineWidthM);

% Plot sum intervals
lA = B_r.plot('c','linewidth',lineWidthL,'DisplayName','Rectangular');
lB = B_g.plot('b','linewidth',lineWidthL,'DisplayName','Polygonal');
lC = B_x.plot('r','linewidth',lineWidthM,'DisplayName','Polyarcular');

% Axis labels
xlabel('Real')
ylabel('Imag')

% Get figure limits
xL = xlim(); yL = ylim();

% Plot supremum
l1 = fimplicit(@(x,y) x.^2+y.^2-P_r.sup,fBox,'c-.','linewidth',lineWidthM,...
                                    'DisplayName','Zhang');
l2 = fimplicit(@(x,y) x.^2+y.^2-P_g.sup,fBox,'b:','linewidth',lineWidthM,...
                                    'DisplayName','Tenuti');
l3 = fimplicit(@(x,y) x.^2+y.^2-P_x.sup,fBox,'r:','linewidth',lineWidthM,...
                                    'DisplayName','Geréb');

% Plot Infimum
fimplicit(@(x,y) x.^2+y.^2-P_r.inf,fBox,'c-.','linewidth',lineWidthM);
fimplicit(@(x,y) x.^2+y.^2-P_g.inf,fBox,'b--','linewidth',lineWidthM);
fimplicit(@(x,y) x.^2+y.^2-P_x.inf,fBox,'r:','linewidth',lineWidthM);


% Set figure limits
xlim(xL-0.03); ylim(yL)


% Interval label
for m = 1:M
    text(real(AF_nom(m))+0.02,imag(AF_nom(m)),['$E_{' num2str(m) '}^I$'], ...
                'HorizontalAlignment','center', 'Interpreter','latex')
end
text(B_r.real.mid,B_r.imag.mid,'$B^I\!=\!\sum_n E_n^I$',...
                    'HorizontalAlignment','center', 'Interpreter','latex')
text(sqrt(P_g.inf)-0.01,0,'$\underline{|B^I|}$',...
                    'HorizontalAlignment','center', 'Interpreter','latex')
text(sqrt(P_g.sup)+0.01,0,'$\overline{|B^I|}$',...
                    'HorizontalAlignment','center', 'Interpreter','latex')

% Add two legend windows
legend([lA(1),lB(1),lC(1)],'Location','NorthWest')
ax2 = axes('position',get(gca,'position'),'visible','off');
legend(ax2, [l1,l2,l3], 'Location','northeast');

% Tightness values
annotText = {['$\tau_{P^I}^{\mathrm{(Zh)}}=' ...
                num2str(100*P_a.Width ./ P_r.Width,3) '\%$'],...
               ['$\tau_{P^I}^{\mathrm{(Te)}}=' ...
               num2str(100*P_a.Width ./ P_g.Width,3) '\%$'],...
               ['$\tau_{P^I}^{\mathrm{(Ge)}}=' ...
               num2str(100*P_a.Width ./ P_x.Width,3) '\%$']};
annotation('textbox',[0.30 0.45 .15 0.20],'String',annotText, ...
           'BackgroundColor','w',...
           'HorizontalAlignment','right', ...
           'Interpreter','latex');

% Add zoom window for the operand interval
axes('position',[.14 .13 .11 .25]); hold on; box on
AF_r.plot('c','linewidth',lineWidthL);
AF_g.plot('b','linewidth',lineWidthL);
AF_x.plot('r','linewidth',lineWidthS);
axis equal
xlim([0.0585 0.0775])
xticks(-0.5)
xtickangle(90)
set(gca, 'XAxisLocation', 'top')
ylim([0.112 0.133])
yticks(0.63)
set(gca, 'YAxisLocation', 'right')
m=6;
text(real(AF_nom(m)),imag(AF_nom(m)),['$E_{' num2str(m) '}^I$'], ...
                'HorizontalAlignment','center', 'Interpreter','latex')

% Add zoom window for the infimum
axes('position',[.56 .13 .07 .15]); hold on; box on
fimplicit(@(x,y) x.^2+y.^2-P_r.inf,fBox,'c-.','linewidth',lineWidthL);
fimplicit(@(x,y) x.^2+y.^2-P_g.inf,fBox,'b--','linewidth',lineWidthL);
fimplicit(@(x,y) x.^2+y.^2-P_x.inf,fBox,'r:','linewidth',lineWidthL);
axis equal
xlim([0.195 0.196])
xticks([])
xtickangle(90)
set(gca, 'XAxisLocation', 'top')
yticks(0)
yticklabels('0')
text(sqrt(P_g.inf+0.0002),0,'$\underline{|B^I|}$',...
                    'HorizontalAlignment','center', 'Interpreter','latex')


% Add zoom window for the supremum
axes('position',[.78 .13 .07 .15]); hold on; box on
fimplicit(@(x,y) x.^2+y.^2-P_r.sup,fBox,'c-.','linewidth',lineWidthL);
fimplicit(@(x,y) x.^2+y.^2-P_g.sup,fBox,'b:','linewidth',lineWidthL);
fimplicit(@(x,y) x.^2+y.^2-P_x.sup,fBox,'r:','linewidth',lineWidthL);
axis equal
xlim([0.315 0.325])
xticks(1.9)
set(gca, 'XAxisLocation', 'top')
xtickangle(90)
yticks(0)
text(sqrt(P_g.sup)+0.003,0,'$\overline{|B^I|}$',...
                    'HorizontalAlignment','center', 'Interpreter','latex')

% Set font size
fontsize(30,'point')

