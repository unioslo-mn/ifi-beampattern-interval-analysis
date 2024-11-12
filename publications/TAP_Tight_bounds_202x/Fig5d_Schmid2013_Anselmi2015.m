clear
% close all

%%

% Set parameters
M = 8;
SNR = 20;
% alpha = [5 5 10 10 10 10 5 5]'*1e-2;
% beta = [1 2 3 4 3 2 1]' * 1e-2;
dTheta = 0.85;
theta = linspace(-pi,pi,M)' * dTheta;
w = chebwin(M,SNR) ; w = w / sum(w);
v = exp(1j*theta);
B = w' * v;

% Set error
ampErr = 0.1;
phaErr = deg2rad(10);
gamma = 0.05;

% Generate polar and circular intervals
% Eint = ciat.PolarInterval(w.*(1+[-1 1].*alpha*0.75) , angle(v')' + [-1 1].*phaErr*pi);
Eint = ciat.PolarInterval(w * (1 + ciat.RealInterval(-ampErr/2,ampErr/2)),...
                          ciat.RealInterval(theta + [-1 1]*phaErr/2));
% alpha = ciat.CircularInterval(Eint).Radius ./ w;
alpha = ones(M,1) * sqrt(ampErr^2 + sin(phaErr)^2)/2;
beta = ones(M-1,1)*gamma;
Aint = ciat.CircularInterval(zeros(M,1) , [0;beta]+[beta;0]);

% Element intervals (coupled)
AF_nom = w .* v;
AF_g = ciat.PolygonalInterval( Eint , ones(M,1)+Aint ,'tolerance',1e-3);
AF_x = ciat.PolyarxInterval( Eint , ones(M,1)+Aint );
AF_a = ciat.PolyarcularInterval( Eint , ones(M,1)+Aint );

% Beampattern interval
B_g = sum(AF_g);
B_x = sum(AF_x);
B_a = sum(AF_a);

% Schmid method
wL2 = sqrt(sum(abs(w)*.2));
maxDeltaB = sqrt(M) * wL2 * (w'*alpha + w(1:end-1)'*beta + w(2:end)'*beta);
B_Sch = ciat.CircularInterval(B,maxDeltaB);

% Anselmi method
diagC0 = diag(alpha);
diagC1 = diag(beta,1);
diagCm1 = diag(beta,-1);
I = eye(M);
% CaRad = radius(ciat.CircularInterval(Eint))./w;
Ca = ciat.CircularInterval(zeros(M) , diagC0);
Cb = ciat.CircularInterval(zeros(M) , diagC1 + diagCm1);
AF_An = (w.*v)' * (Ca + Cb + I);
% AF_An = ciat.CircularInterval(Eint) .* (1+Aint);
% AF_An = ciat.CircularInterval(AF_a);
B_Ans = w' * (Ca + Cb + I) * v;
% B_Ans = sum(AF_An);

% Power intervals
P_Sch = abs(B_Sch).^2;
P_Ans = abs(B_Ans).^2;
P_g = abs(B_g).^2;
P_x = abs(B_x).^2;
P_a = abs(B_a).^2;


%% Plot

% Set parameters
lineWidthM = 4;
lineWidthS = 3;
lineWidthXS = 1; 
fBox = [0 0.5 -1.1 1.1];

% Plot
figure(1);clf;hold on;axis equal;
set(gca,'DefaultLineLineWidth',lineWidthM)
plot(0,0,'k+')

% Plot operand intervals
AF_g.plot('b','linewidth',lineWidthM);
AF_x.plot('r','linewidth',lineWidthS);
AF_An.plot('m','linewidth',lineWidthS);

% Plot sum intervals
lA = B_Sch.plot('c','linewidth',lineWidthM,'DisplayName','Circular');
lB = B_Ans.plot('m','linewidth',lineWidthM,'DisplayName','Circular');
lC = B_g.plot('b','linewidth',lineWidthM,'DisplayName','Polygonal');
lD = B_x.plot('r','linewidth',lineWidthS,'DisplayName','Polyarcular');

% Axis labels
xlabel('Real')
ylabel('Imag')

% Get figure limits
xL = xlim(); yL = ylim();

% Plot supremum
l1 = fimplicit(@(x,y) x.^2+y.^2-P_Sch.sup,fBox,'c-.','linewidth',lineWidthS,...
                                     'DisplayName','Schmid');
l2 = fimplicit(@(x,y) x.^2+y.^2-P_Ans.sup,fBox,'m-.','linewidth',lineWidthS,...
                                     'DisplayName','Anselmi');
l3 = fimplicit(@(x,y) x.^2+y.^2-P_g.sup,fBox,'b:','linewidth',lineWidthS,...
                                    'DisplayName','Tenuti');
l4 = fimplicit(@(x,y) x.^2+y.^2-P_x.sup,fBox,'r:','linewidth',lineWidthS,...
                                    'DisplayName','Geréb');

% Plot Infimum
fimplicit(@(x,y) x.^2+y.^2-P_Sch.inf,fBox,'c-.','linewidth',lineWidthS);
fimplicit(@(x,y) x.^2+y.^2-P_Ans.inf,fBox,'m-.','linewidth',lineWidthS);
fimplicit(@(x,y) x.^2+y.^2-P_g.inf,fBox,'b--','linewidth',lineWidthS);
fimplicit(@(x,y) x.^2+y.^2-P_x.inf,fBox,'r:','linewidth',lineWidthS);


% Set figure limits
xlim(xL-0.03); ylim(yL)


% Interval label
for n = 1:M
    if real(AF_nom(n)) > 0
        text(real(AF_nom(n)),imag(AF_nom(n)),['$E_{' num2str(n) '}^I$'], ...
                'HorizontalAlignment','center', 'Interpreter','latex')
    else
        text(real(AF_nom(n))+0.015,imag(AF_nom(n)),['$E_{' num2str(n) '}^I$'], ...
                'HorizontalAlignment','center', 'Interpreter','latex')
    end
end
text(B_g.real.mid,B_g.imag.mid,'$B^I\!=\!\sum_n E_n^I$',...
                    'HorizontalAlignment','center', 'Interpreter','latex')
text(sqrt(P_Ans.inf)-0.02,0,'$\underline{|B^I|}$',...
                    'HorizontalAlignment','center', 'Interpreter','latex')
text(sqrt(P_Ans.sup)+0.02,0,'$\overline{|B^I|}$',...
                    'HorizontalAlignment','center', 'Interpreter','latex')

% Add two legend windows
legend([lA(1),lB(1),lC(1),lD(1)],'Location','NorthWest')
ax2 = axes('position',get(gca,'position'),'visible','off');
legend(ax2, [l1,l2,l3,l4], 'Location','northeast');

% Tightness values
annotText = {   ['$\tau_{P^I}^{\mathrm{(Sc)}}=' ...
                num2str(100*P_a.Width ./ P_Sch.Width,3) '\%$'],...
                ['$\tau_{P^I}^{\mathrm{(An)}}=' ...
                num2str(100*P_a.Width ./ P_Ans.Width,3) '\%$'],...
               ['$\tau_{P^I}^{\mathrm{(Te)}}=' ...
               num2str(100*P_a.Width ./ P_g.Width,3) '\%$'],...
               ['$\tau_{P^I}^{\mathrm{(Ge)}}=' ...
               num2str(100*P_a.Width ./ P_x.Width,3) '\%$']};
annotation('textbox',[0.14 0.13 .14 0.25],'String',annotText, ...
           'BackgroundColor','w',...
           'HorizontalAlignment','right', ...
           'Interpreter','latex');

% Add zoom window for the operand interval
axes('position',[.14 .43 .25 .25]); hold on; box on
AF_g.plot('b','linewidth',lineWidthM);
AF_x.plot('r','linewidth',lineWidthS);
AF_An.plot('m','linewidth',lineWidthS);
axis equal
xlim([-0.075 -0.045])
xticks(-0.5)
xtickangle(90)
set(gca, 'XAxisLocation', 'top')
ylim([0.6 0.72])
yticks(0.63)
set(gca, 'YAxisLocation', 'right')
n=6;
text(real(AF_nom(n)),imag(AF_nom(n))+0.01,['$E_{' num2str(n) '}^I$'], ...
                'HorizontalAlignment','center', 'Interpreter','latex')

% Add zoom window for the supremum
axes('position',[.82 .13 .09 .3]); hold on; box on
B_Sch.plot('c','linewidth',lineWidthM);
B_Ans.plot('m','linewidth',lineWidthM);
B_g.plot('b','linewidth',lineWidthM);
B_x.plot('r','linewidth',lineWidthS);
fimplicit(@(x,y) x.^2+y.^2-P_Sch.sup,fBox,'c-.','linewidth',lineWidthS);
fimplicit(@(x,y) x.^2+y.^2-P_Ans.sup,fBox,'m-.','linewidth',lineWidthS);
fimplicit(@(x,y) x.^2+y.^2-P_g.sup,fBox,'b--','linewidth',lineWidthS);
fimplicit(@(x,y) x.^2+y.^2-P_x.sup,fBox,'r:','linewidth',lineWidthS);
axis equal
xlim([0.285 0.292])
xticks(0.29)
xtickangle(90)
set(gca, 'XAxisLocation', 'top')
yticks(0)
ylim([-0.01 0.01])

% Finalize figure
fontsize(30,'point')


