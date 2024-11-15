clear
% close all

%% Poli 2015 (phase errors)

% Set parameters
conf = getFig5conf(1);
M = conf.M;
w = conf.w;
theta = conf.theta;
phaErr = conf.phaErr;

% Calculate nominal element phase 
phi = ((0:M-1)-(M-1)/2)' * pi*sin(theta);

% Define and cast intervals
E_nom = w .* exp(1i*phi);
E_p = ciat.PolarInterval(w , ciat.RealInterval(phi + [-1 1]*phaErr/2));
tic;E_g = ciat.PolygonalInterval(E_p,'tolerance',conf.tol);T_g(1)=toc;
tic;E_x = ciat.PolyarxInterval(E_p);T_x(1)=toc;
tic;E_r = ciat.RectangularInterval(E_g);T_r(1)=toc;
E_a = ciat.PolyarcularInterval(E_p);

% Sum intervals
tic;B_r = sum(E_r);T_r(2)=toc;
tic;B_g = sum(E_g);T_g(2)=toc;
tic;B_x = sum(E_x);T_x(2)=toc;
B_a = sum(E_a);

% Calculate power intervals
tic;P_r = abs(B_r).^2;T_r(3)=toc;
tic;P_g = abs(B_g).^2;T_g(3)=toc;
tic;P_x = abs(B_x).^2;T_x(3)=toc;
P_a = abs(B_a).^2;

% Calculate tightness
tau_x = P_a.Width ./ P_x.Width;
tau_g = P_a.Width ./ P_g.Width;
tau_r = P_a.Width ./ P_r.Width;

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
B_abs = sqrt([P_r,P_g,P_x]);

% Plot
figure(2);clf;hold on;axis equal;
set(gca,'DefaultLineLineWidth',lineWidthL)
plot(0,0,'k+')

% Plot operand intervals
E_g.plot('b','linewidth',lineWidthS);
E_x.plot('r','linewidth',lineWidthS);
E_r.plot('color',cList(1,:),'linewidth',lineWidthS);

% Plot sum intervals
lB = B_g.plot('b','linewidth',lineWidthL,'DisplayName','Polygonal');
lC = B_x.plot('r','linewidth',lineWidthS,'DisplayName','Polyarcular');
lA = B_r.plot('color',cList(1,:),'linewidth',lineWidthL,...
                                        'DisplayName','Rectangular');

% Axis labels
xlabel('Real')
ylabel('Imag')

% Plot supremum
l1 = fimplicit(@(x,y) x.^2+y.^2-P_r.sup,fBox,'-.','color',cList(1,:),...
                            'linewidth',lineWidthM,'DisplayName','Poli');
l2 = fimplicit(@(x,y) x.^2+y.^2-P_g.sup,fBox,'b--','linewidth',lineWidthM,...
                                    'DisplayName','Tenuti');
l3 = fimplicit(@(x,y) x.^2+y.^2-P_x.sup,fBox,'r:','linewidth',lineWidthM,...
                                    'DisplayName','Geréb');

% Plot Infimum
fimplicit(@(x,y) x.^2+y.^2-P_r.inf,fBox,'-.','linewidth',lineWidthM,...
                                            'color',cList(1,:));
fimplicit(@(x,y) x.^2+y.^2-P_g.inf,fBox,'b--','linewidth',lineWidthM);
fimplicit(@(x,y) x.^2+y.^2-P_x.inf,fBox,'r:','linewidth',lineWidthM);


% Set figure limits
xlim(xL); ylim(yL)

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


% % Add two legend windows
% legend([lA(1),lB(1),lC(1)],'Location','NorthWest')
% ax2 = axes('position',get(gca,'position'),'visible','off');
% legend(ax2, [l1,l2,l3], 'Location','northeast');

% Add zoom window for the operand interval
axes('position',[0.09,0.18,0.3,0.3]); hold on; box on
E_p.plot('b','linewidth',lineWidthL);
E_r.plot('color',cList(1,:),'linewidth',lineWidthL);
E_x.plot('r','linewidth',lineWidthS);
axis equal
set(gca, 'XAxisLocation', 'top')
m=6;
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
fimplicit(@(x,y) x.^2+y.^2-P_r.inf,fBox,'-.', ...
                          'color',cList(1,:),'linewidth',lineWidthL);
fimplicit(@(x,y) x.^2+y.^2-P_g.inf,fBox,'b--','linewidth',lineWidthL);
fimplicit(@(x,y) x.^2+y.^2-P_x.inf,fBox,'r:','linewidth',lineWidthL);
axis equal
xlim([min(B_abs.inf)-1e-3 , max(B_abs.inf)+1e-3])
xticks([])
xtickangle(90)
set(gca, 'XAxisLocation', 'top')
yticks([])
yticklabels('0')
text(mean(B_abs.inf),0,'$\underline{|B^I|}$',...
        'HorizontalAlignment','center', 'Interpreter','latex')


% Add zoom window for the supremum
axes('position',[0.730,0.62,0.14,0.30]); hold on; box on
fimplicit(@(x,y) x.^2+y.^2-P_r.sup,fBox,'-.', ...
                        'color',cList(1,:),'linewidth',lineWidthL);
fimplicit(@(x,y) x.^2+y.^2-P_g.sup,fBox,'b--','linewidth',lineWidthL);
fimplicit(@(x,y) x.^2+y.^2-P_x.sup,fBox,'r:','linewidth',lineWidthL);
axis equal
xlim([min(B_abs.sup)-1e-3 , max(B_abs.sup)+1e-3])
xticks([])
set(gca, 'XAxisLocation', 'top')
xtickangle(90)
yticks([])
text(mean(B_abs.sup),0,'$\overline{|B^I|}$',...
        'HorizontalAlignment','center', 'Interpreter','latex')

% Set font size
fontsize(40,'point')

% Tightness values
annotText = {join(['\color{red}\tau^{(a)}=' num2str(100*tau_x,3) '%' ...
                    ', T^{(a)}=' join(compose("%0.0f",1e3*T_x),'+') 'ms'],''),...
            join(['\color{blue}\tau^{(g)}=' num2str(100*tau_g,3) '%' ...
                    ', T^{(g)}=' join(compose("%0.0f",1e3*T_g),'+') 'ms'],''),...
            join(['\color[rgb]{' num2str(cList(1,:)) '}\tau^{(r)}=' ...
                        num2str(100*tau_r,3) '%' ...
                    ', T^{(r)}=' join(compose("%0.0f",1e3*T_r),'+') 'ms'],'')};
annotation('textbox',[0.164,0.553,0.280 0.367],'String',annotText, ...
           'BackgroundColor','w','VerticalAlignment','top','fontSize',30,...
           'HorizontalAlignment','center','FitBoxToText','on');