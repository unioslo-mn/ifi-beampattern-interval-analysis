clear
% close all

%% Zhang 2017 (amplitude and phase errors)

% Set parameters
conf = getFig5conf(1);
M = conf.M;
w = conf.w;
theta = conf.theta;
ampErr = conf.ampErr;
phaErr = conf.phaErr;

% Calculate nominal element phase 
phi = ((0:M-1)-(M-1)/2)' * pi*sin(theta);

% Calculate amplitude interval
A = w * (1 + ciat.RealInterval(-ampErr/2,ampErr/2));

% Define and cast intervals
EA_nom = w .* exp(1i*phi);
EA_p = ciat.PolarInterval(w * (1 + ciat.RealInterval(-ampErr/2,ampErr/2)),...
                          ciat.RealInterval(phi + [-1 1]*phaErr/2));
tic;EA_g = ciat.PolygonalInterval(EA_p,'tolerance',conf.tol);T_g(1)=toc;
tic;EA_x = ciat.PolyarxInterval(EA_p);T_x(1)=toc;
tic;EA_r = ciat.RectangularInterval(EA_g);T_r(1)=toc;
EA_a = ciat.PolyarcularInterval(EA_p);

% Sum intervals
tic;B_r = sum(EA_r);T_r(2)=toc;
tic;B_g = sum(EA_g);T_g(2)=toc;
tic;B_x = sum(EA_x);T_x(2)=toc;
B_a = sum(EA_a);

% Power intervals
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
fBox = [0 3 B_r.Imag.inf B_r.Imag.sup];
xL = conf.xL;
yL = conf.yL;

% Get beampattern interval absolute values
B_abs = sqrt([P_r,P_g,P_x]);

% Plot
figure(3);clf;hold on;axis equal;
set(gca,'DefaultLineLineWidth',lineWidthL)
plot(0,0,'k+')

% Plot operand intervals
EA_r.plot('c','linewidth',lineWidthS);
EA_g.plot('b','linewidth',lineWidthS);
EA_x.plot('r','linewidth',lineWidthS);

% Plot sum intervals
lA = B_r.plot('c','linewidth',lineWidthL,'DisplayName','Rectangular');
lB = B_g.plot('b','linewidth',lineWidthL,'DisplayName','Polygonal');
lC = B_x.plot('r','linewidth',lineWidthM,'DisplayName','Polyarcular');

% Axis labels
xlabel('Real')
ylabel('Imag')

% Plot supremum
l1 = fimplicit(@(x,y) x.^2+y.^2-P_r.sup,fBox,'-.','linewidth',lineWidthM,...
                                'color',cList(1,:),'DisplayName','Zhang');
l2 = fimplicit(@(x,y) x.^2+y.^2-P_g.sup,fBox,'b:','linewidth',lineWidthM,...
                                    'DisplayName','Tenuti');
l3 = fimplicit(@(x,y) x.^2+y.^2-P_x.sup,fBox,'r:','linewidth',lineWidthM,...
                                    'DisplayName','Geréb');

% Plot Infimum
fimplicit(@(x,y) x.^2+y.^2-P_r.inf,fBox,'-.', ...
                            'color',cList(1,:),'linewidth',lineWidthM);
fimplicit(@(x,y) x.^2+y.^2-P_g.inf,fBox,'b--','linewidth',lineWidthM);
fimplicit(@(x,y) x.^2+y.^2-P_x.inf,fBox,'r:','linewidth',lineWidthM);


% Set figure limits
xlim(xL); ylim(yL)


% Interval label
for m = 1:M
    text(real(EA_nom(m))+0.01,imag(EA_nom(m)), ...
                ['$E_{' num2str(m) '}^I A_{' num2str(m) '}^I$'], ...
                'HorizontalAlignment','left', 'Interpreter','latex')
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

% % Tightness values
% annotText = {['\color{red}\tau^{(a)}=' num2str(100*tau_x,3) '%'],...
%             ['\color{blue}\tau^{(g)}=' num2str(100*tau_g,3) '%'],...
%             ['\color[rgb]{' num2str(cList(1,:)) '}\tau^{(r)}=' ...
%                         num2str(100*tau_r,3) '%']};
% annotation('textbox',[0.14 0.66 .14 0.25],'String',annotText, ...
%            'BackgroundColor','w','VerticalAlignment','top',...
%            'HorizontalAlignment','center','FitBoxToText','on');

% Add zoom window for the operand interval
axes('position',[0.09,0.18,0.3,0.3]); hold on; box on
EA_p.plot('b','linewidth',lineWidthL);
EA_r.plot('color',cList(1,:),'linewidth',lineWidthL);
EA_x.plot('r','linewidth',lineWidthS);
axis equal
set(gca, 'XAxisLocation', 'top')
n=6;
maxWidth = max([EA_r(n).Real.Width , EA_r(n).Imag.Width]);
xlim(EA_r(n).Real.Midpoint + [-1 1]*maxWidth/2)
ylim(EA_r(n).Imag.Midpoint + [-1 1]*maxWidth/2)
xticks([]);
yticks([]);
set(gca, 'YAxisLocation', 'right')
text(real(EA_nom(n)),imag(EA_nom(n)), ...
     ['$E_{' num2str(n) '}^I A_{' num2str(n) '}^I$'], ...
                'HorizontalAlignment','center', 'Interpreter','latex')

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