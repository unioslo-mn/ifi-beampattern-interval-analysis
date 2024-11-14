clear
% close all

%%

% Set parameters
conf = getFig5conf(1);
M = conf.M;
w = conf.w;
theta = conf.theta;
ampErr = conf.ampErr;
phaErr = conf.phaErr;
mtlCpl = conf.mtlCpl;


% Calculate nominal element phase 
phi = ((0:M-1)-(M-1)/2)' * pi*sin(theta);
v = exp(1j*phi);

% Calculate nominal beampattern interval
B = w' * v;

% Generate error vectors
alpha = ones(M,1) * sqrt(ampErr^2 + sin(phaErr)^2)/2;
beta = ones(M-1,1)*mtlCpl;

% Generate polar and circular intervals
Eint = ciat.PolarInterval(w * (1 + ciat.RealInterval(-ampErr/2,ampErr/2)),...
                          ciat.RealInterval(phi + [-1 1]*phaErr/2));
Aint = ciat.CircularInterval(zeros(M,1) , [0;beta]+[beta;0]);

% Element intervals (coupled)
AF_nom = w .* v;
AF_a = ciat.PolyarcularInterval( Eint , ones(M,1)+Aint );
tic;AF_g = ciat.PolygonalInterval( Eint , ones(M,1)+Aint ,'tolerance',conf.tol);T_g(1)=toc;
tic;AF_x = ciat.PolyarxInterval( Eint , ones(M,1)+Aint );T_x(1)=toc;
tic;AF_r = ciat.RectangularInterval(AF_x);T_r(1)=toc;
tic;AF_c = ciat.CircularInterval(Eint) .* (ones(M,1)+Aint);T_c(1)=toc;


% Beampattern interval
tic;B_r = sum(AF_r);T_r(2)=toc;
tic;B_g = sum(AF_g);T_g(2)=toc;
tic;B_x = sum(AF_x);T_x(2)=toc;
tic;B_c = sum(AF_c);T_c(2)=toc;
B_a = sum(AF_a);

% Schmid method
tic;
wL2 = sqrt(sum(abs(w)*.2));
maxDeltaB = sqrt(M) * wL2 * (w'*alpha + w(1:end-1)'*beta + w(2:end)'*beta);
B_Sch = ciat.CircularInterval(B,maxDeltaB);
T_Sc(2) = toc;

% Anselmi method
tic;
diagC0 = diag(alpha);
diagC1 = diag(beta,1);
diagCm1 = diag(beta,-1);
I = eye(M);
Ca = ciat.CircularInterval(zeros(M) , diagC0);
Cb = ciat.CircularInterval(zeros(M) , diagC1 + diagCm1);
B_Ans = w' * (Ca + Cb + I) * v;
T_An(2) = toc;
AF_An = (w.*v)' * (Ca + Cb + I);


% Power intervals
tic;P_Sch = abs(B_Sch).^2;T_Sc(3)=toc;
tic;P_Ans = abs(B_Ans).^2;T_An(3)=toc;
tic;P_g = abs(B_g).^2;T_g(3)=toc;
tic;P_x = abs(B_x).^2;T_x(3)=toc;
tic;P_r = abs(B_r).^2;T_r(3)=toc;
tic;P_c = abs(B_c).^2;T_c(3)=toc;
P_a = abs(B_a).^2;

% Calculate tightness
tau_x = P_a.Width ./ P_x.Width;
tau_g = P_a.Width ./ P_g.Width;
tau_r = P_a.Width ./ P_r.Width;
tau_c = P_a.Width ./ P_c.Width;
tau_Sc = P_a.Width ./ P_Sch.Width;
tau_An = P_a.Width ./ P_Ans.Width;


%% Plot

% Set parameters
lineWidthM = 4;
lineWidthS = 3;
lineWidthXS = 1; 
cList = getColorList(conf.cID);
fBox = [0 3 B_r.Imag.inf B_r.Imag.sup];
xL = conf.xL;
yL = conf.yL;

% Get beampattern interval absolute values
B_abs = sqrt([P_Sch,P_Ans,P_g,P_x]);

% Plot
figure(4);clf;hold on;axis equal;
set(gca,'DefaultLineLineWidth',lineWidthM)
plot(0,0,'k+')

% Plot operand intervals
AF_g.plot('b','linewidth',lineWidthS);
AF_x.plot('r','linewidth',lineWidthS);
% AF_c.plot('color',cList(2,:),'linewidth',lineWidthS);
AF_An.plot('color',cList(4,:),'linewidth',lineWidthS);

% Plot sum intervals
lA = B_Sch.plot('color',cList(3,:),'linewidth',lineWidthM,'DisplayName','Schmid');
lB = B_Ans.plot('color',cList(4,:),'linewidth',lineWidthM,'DisplayName','Anselmi');
lC = B_g.plot('b','linewidth',lineWidthM,'DisplayName','Polygonal');
lD = B_x.plot('r','linewidth',lineWidthS,'DisplayName','Polyarcular');
lE = B_g.plot('color',cList(2,:),'linewidth',lineWidthM,'DisplayName','Circular');

% Axis labels
xlabel('Real')
ylabel('Imag')

% Plot supremum
l1 = fimplicit(@(x,y) x.^2+y.^2-P_Sch.sup,fBox,'-.','linewidth',lineWidthS,...
                             'color',cList(3,:),'DisplayName','Schmid');
l2 = fimplicit(@(x,y) x.^2+y.^2-P_Ans.sup,fBox,'-.','linewidth',lineWidthS,...
                             'color',cList(4,:),'DisplayName','Anselmi');
% l3 = fimplicit(@(x,y) x.^2+y.^2-P_c.sup,fBox,':','linewidth',lineWidthS,...
%                                 'color',cList(2,:),'DisplayName','Circular');
l4 = fimplicit(@(x,y) x.^2+y.^2-P_g.sup,fBox,'b:','linewidth',lineWidthS,...
                                    'DisplayName','Tenuti');
l5 = fimplicit(@(x,y) x.^2+y.^2-P_x.sup,fBox,'r:','linewidth',lineWidthS,...
                                    'DisplayName','Geréb');


% Plot Infimum
fimplicit(@(x,y) x.^2+y.^2-P_Sch.inf,fBox,'-.', ...
                    'color',cList(3,:),'linewidth',lineWidthS);
fimplicit(@(x,y) x.^2+y.^2-P_Ans.inf,fBox,'-.', ...
                    'color',cList(4,:),'linewidth',lineWidthS);
% fimplicit(@(x,y) x.^2+y.^2-P_c.inf,fBox,':', ...
%                     'color',cList(2,:),'linewidth',lineWidthS);
fimplicit(@(x,y) x.^2+y.^2-P_g.inf,fBox,'b--','linewidth',lineWidthS);
fimplicit(@(x,y) x.^2+y.^2-P_x.inf,fBox,'r:','linewidth',lineWidthS);

% Set figure limits
xlim(xL); ylim(yL)

% Interval label
for n = 1:M
    text(real(AF_nom(n))+0.02,imag(AF_nom(n)),['$E_{' num2str(n) '}^I$'], ...
            'HorizontalAlignment','center', 'Interpreter','latex')
end
text(B_r.real.mid,B_r.imag.mid,'$B^I$',...
                    'HorizontalAlignment','center', 'Interpreter','latex')
text(B_r.real.inf-0.03,B_r.imag.mid,'$\underline{|B^I|}$',...
                    'HorizontalAlignment','right', 'Interpreter','latex')
text(B_r.real.sup+0.03,B_r.imag.mid,'$\overline{|B^I|}$',...
                    'HorizontalAlignment','left', 'Interpreter','latex')

% % Add two legend windows
% legend([lA(1),lB(1),lC(1),lD(1)],'Location','NorthWest')
% ax2 = axes('position',get(gca,'position'),'visible','off');
% legend(ax2, [l1,l2,l3,l4], 'Location','northeast');

% % Tightness values
% annotText = {['\color{red}\tau^{(a)}=' num2str(100*tau_x,3) '%'],...
%             ['\color{blue}\tau^{(g)}=' num2str(100*tau_g,3) '%'],...
%             ['\color[rgb]{' num2str(cList(4,:)) '}\tau^{(A)}=' ...
%                         num2str(100*tau_Ans,3) '%'],...
%             ['\color[rgb]{' num2str(cList(3,:)) '}\tau^{(S)}=' ...
%                         num2str(100*tau_Sch,3) '%']};
% annotation('textbox',[0.14 0.66 .14 0.25],'String',annotText, ...
%            'BackgroundColor','w','VerticalAlignment','top',...
%            'HorizontalAlignment','center','FitBoxToText','on');

% Add zoom window for the operand interval
axes('position',[0.09,0.18,0.3,0.3]); hold on; box on
AF_g.plot('b','linewidth',lineWidthM);
AF_x.plot('r','linewidth',lineWidthS);
% AF_c.plot('color',cList(2,:),'linewidth',lineWidthS);
AF_An.plot('color',cList(4,:),'linewidth',lineWidthS);
axis equal
set(gca, 'XAxisLocation', 'top')
n=6;
maxWidth = max([AF_c(n).Real.Width , AF_c(n).Imag.Width]);
xlim(AF_c(n).Real.Midpoint + [-1 1]*maxWidth/2)
ylim(AF_c(n).Imag.Midpoint + [-1 1]*maxWidth/2)
xticks([]);
yticks([]);
set(gca, 'YAxisLocation', 'right')
text(real(AF_nom(n)),imag(AF_nom(n)),['$E_{' num2str(n) '}^I$'], ...
                'HorizontalAlignment','right', 'Interpreter','latex')


% Add zoom window for the infimum
axes('position',[0.73,0.178,0.14,0.3]); hold on; box on
fimplicit(@(x,y) x.^2+y.^2-P_Sch.inf,fBox,'-.', ...
                    'color',cList(3,:),'linewidth',lineWidthS);
fimplicit(@(x,y) x.^2+y.^2-P_Ans.inf,fBox,'-.', ...
                    'color',cList(4,:),'linewidth',lineWidthS);
% fimplicit(@(x,y) x.^2+y.^2-P_c.inf,fBox,':', ...
%                     'color',cList(2,:),'linewidth',lineWidthS);
fimplicit(@(x,y) x.^2+y.^2-P_g.inf,fBox,'b--','linewidth',lineWidthS);
fimplicit(@(x,y) x.^2+y.^2-P_x.inf,fBox,'r:','linewidth',lineWidthS);
axis equal
xlim([min(B_abs.inf)-1e-3 , max(B_abs.inf)+1e-3])
xticks([])
set(gca, 'XAxisLocation', 'top')
xtickangle(90)
yticks([])
text(mean(B_abs.inf),0,'$\underline{|B^I|}$',...
        'HorizontalAlignment','center', 'Interpreter','latex')

% Add zoom window for the supremum
axes('position',[0.730,0.62,0.14,0.30]); hold on; box on
fimplicit(@(x,y) x.^2+y.^2-P_Sch.sup,fBox,'-.', ...
                        'color',cList(3,:),'linewidth',lineWidthS);
fimplicit(@(x,y) x.^2+y.^2-P_Ans.sup,fBox,'-.', ...
                        'color',cList(4,:),'linewidth',lineWidthS);
% fimplicit(@(x,y) x.^2+y.^2-P_c.sup,fBox,':', ...
%                         'color',cList(2,:),'linewidth',lineWidthS);
fimplicit(@(x,y) x.^2+y.^2-P_g.sup,fBox,'b--','linewidth',lineWidthS);
fimplicit(@(x,y) x.^2+y.^2-P_x.sup,fBox,'r:','linewidth',lineWidthS);
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
            join(['\color[rgb]{' num2str(cList(4,:)) '}\tau^{(A)}=' ...
                        num2str(100*tau_An,3) '%' ...
                    ', T^{(A)}=' join(compose("%0.0f",1e3*T_An),'+') 'ms'],''),...
            join(['\color[rgb]{' num2str(cList(3,:)) '}\tau^{(S)}=' ...
                        num2str(100*tau_Sc,3) '%' ...
                    ', T^{(S)}=' join(compose("%0.0f",1e3*T_Sc),'+') 'ms'],'')};
annotation('textbox',[0.164,0.553,0.280 0.367],'String',annotText, ...
           'BackgroundColor','w','VerticalAlignment','top','fontSize',30,...
           'HorizontalAlignment','center','FitBoxToText','on');
