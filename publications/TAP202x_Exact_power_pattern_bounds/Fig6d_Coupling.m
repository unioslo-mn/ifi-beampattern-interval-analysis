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
cplCoeff = conf.cplCoeff;
nT = conf.nT;

% Calculate nominal element phase 
phi = ((0:M-1)-(M-1)/2)' * pi*sin(theta);
v = exp(1j*phi);

% Calculate nominal beampattern interval
B = w' * v;

% Generate error vectors
alpha = ones(M,1) * sqrt(ampErr^2 + sin(phaErr)^2)/2;
beta = ones(M-1,1)*cplCoeff;

% Generate the element intervals
Eint = ciat.PolarInterval((ones(M,1) + ciat.RealInterval(-ampErr/2,ampErr/2)),...
                          ciat.RealInterval(phi + [-1 1]*phaErr/2));

% Generate the coupling interval (A)
p_m = w;
R_m = [0;w(1:M-1)].*[0;beta] + [w(2:M);0].*[beta;0];
Aint = ciat.CircularInterval(p_m , R_m);

% Element intervals (without coupling)
E_a = ciat.PolyarxInterval(Eint);

% Measure running time on many occurences and average the result
T_x    = zeros(3,1);
T_g    = zeros(3,1);
T_r    = zeros(3,1);
T_c    = zeros(3,1);
T_Sc    = zeros(3,1);
T_An    = zeros(3,1);
T_Ar    = zeros(3,1);

for iT = 1:nT
    fprintf("%0i,",iT)
    % Element intervals (coupled)
    EA_nom = w .* v;
    B_nom = sum(EA_nom);
    EA_a = ciat.PolyarcularInterval( Eint , Aint );
    tic;EA_g = ciat.PolygonalInterval( Eint , Aint ,'tolerance',conf.tol);T_g(1)=T_g(1)+toc;
    tic;EA_x = ciat.PolyarxInterval( Eint , Aint );T_x(1)=T_x(1)+toc;
    tic;EA_r = ciat.RectangularInterval(EA_x);T_r(1)=T_r(1)+toc;
    tic;EA_c = ciat.CircularInterval(Eint) .* (Aint);T_c(1)=T_c(1)+toc;
    
    
    % Beampattern interval
    tic;B_r = sum(EA_r);T_r(2)=T_r(2)+toc;
    tic;B_g = sum(EA_g);T_g(2)=T_g(2)+toc;
    tic;B_x = sum(EA_x);T_x(2)=T_x(2)+toc;
    tic;B_c = sum(EA_c);T_c(2)=T_c(2)+toc;
    B_a = sum(EA_a);
    
    % Schmid method
    tic;
    wL2 = sqrt(sum(abs(w)*.2));
    maxDeltaB = sqrt(M) * wL2 * (w'*alpha + w(1:end-1)'*beta + w(2:end)'*beta);
    B_Sch = ciat.CircularInterval(B,maxDeltaB);
    T_Sc(2) = T_Sc(2) + toc;
    
    % Anselmi method
    tic;
    diagC0 = diag(alpha);
    diagC1 = diag(beta,1);
    diagCm1 = diag(beta,-1);
    I = eye(M);
    Ca = ciat.CircularInterval(zeros(M) , diagC0);
    Cb = ciat.CircularInterval(zeros(M) , diagC1 + diagCm1);
    B_Ans = w' * (Ca + Cb + I) * v;
    T_An(2) = T_An(2)+toc;
    EA_An = (w.*v)' * (Ca + Cb + I);
    
    
    % Power intervals
    tic;P_Sc = abs(B_Sch).^2;T_Sc(3)=T_Sc(3)+toc;
    tic;P_An = abs(B_Ans).^2;T_An(3)=T_An(3)+toc;
    tic;P_g = abs(B_g).^2;T_g(3)=T_g(3)+toc;
    tic;P_x = abs(B_x).^2;T_x(3)=T_x(3)+toc;
    tic;P_r = abs(B_r).^2;T_r(3)=T_r(3)+toc;
    tic;P_c = abs(B_c).^2;T_c(3)=T_c(3)+toc;
    P_a = abs(B_a).^2;

    % Arnestad's method
    tic;
    P_Ar_sup = ( abs(B_nom) + sqrt((ampErr/2)^2 + (phaErr/2)^2) + 2*cplCoeff )^2;
    P_Ar_inf = ( abs(B_nom) - sqrt((ampErr/2)^2 + (phaErr/2)^2) - 2*cplCoeff )^2;
    P_Ar_inf(P_Ar_inf<0) = 0;
    P_Ar = ciat.RealInterval(P_Ar_inf,P_Ar_sup);
    T_Ar(3) = T_Ar(3)+ toc;
end

% Calculate average time
T_g    = T_g/nT;
T_r    = T_r/nT;
T_x    = T_x/nT;
T_c    = T_c/nT;
T_Sc   = T_Sc/nT;
T_An   = T_An/nT;
T_Ar   = T_Ar/nT;

% Calculate tightness
tau_x = P_a.Width ./ P_x.Width;
tau_g = P_a.Width ./ P_g.Width;
tau_r = P_a.Width ./ P_r.Width;
tau_c = P_a.Width ./ P_c.Width;
tau_Sc = P_a.Width ./ P_Sc.Width;
tau_An = P_a.Width ./ P_An.Width;
tau_Ar = P_a.Width ./ P_Ar.Width;


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
B_abs = sqrt([P_Sc,P_An,P_g,P_x]);

% Plot
figure(4);clf;hold on;axis equal;
set(gca,'DefaultLineLineWidth',lineWidthM)

% Grey boxes for the zoom windows
rectangle('Position',[-0.1510  -0.0278 0.0243   0.0243],...
    'EdgeColor','none','FaceColor',ones(3,1)*0.95);
rectangle('Position',[0.1722  -0.0152  0.0290   0.0304],...
    'EdgeColor','none','FaceColor',ones(3,1)*0.95);
rectangle('Position',[-0.0010  -0.0096  0.0183 0.0192],...
    'EdgeColor','none','FaceColor',ones(3,1)*0.95);

plot(0,0,'k+')

% Plot operand intervals
% Eint.plot('k-','linewidth',lineWidthXS);
% Aint.plot('k-','linewidth',lineWidthXS);
EA_g.plot('b','linewidth',lineWidthS);
EA_x.plot('r','linewidth',lineWidthS);
% EA_c.plot('color',cList(2,:),'linewidth',lineWidthS);
EA_An.plot('color',cList(4,:),'linewidth',lineWidthS);

% Plot sum intervals
B_Sch.plot('color',cList(3,:),'linewidth',lineWidthM,'DisplayName','Schmid');
B_Ans.plot('color',cList(4,:),'linewidth',lineWidthM,'DisplayName','Anselmi');
B_g.plot('b','linewidth',lineWidthM,'DisplayName','Polygonal');
B_x.plot('r','linewidth',lineWidthS,'DisplayName','Polyarcular');
B_g.plot('color',cList(2,:),'linewidth',lineWidthM,'DisplayName','Circular');

% Axis labels
xlabel('Real')
ylabel('Imag')

% Plot supremum
fimplicit(@(x,y) x.^2+y.^2-P_Sc.sup,fBox,'-.','linewidth',lineWidthS,...
                        'color',cList(3,:),'DisplayName','Schmid');
fimplicit(@(x,y) x.^2+y.^2-P_An.sup,fBox,'-.','linewidth',lineWidthS,...
                             'color',cList(4,:),'DisplayName','Anselmi');
fimplicit(@(x,y) x.^2+y.^2-P_Ar.sup,fBox,'-.','linewidth',lineWidthS,...
                             'color',cList(2,:),'DisplayName','Arnestad');
fimplicit(@(x,y) x.^2+y.^2-P_g.sup,fBox,'b:','linewidth',lineWidthS,...
                               'DisplayName','Tenuti');
fimplicit(@(x,y) x.^2+y.^2-P_x.sup,fBox,'r:','linewidth',lineWidthS,...
                                    'DisplayName','Geréb');


% Plot Infimum
fimplicit(@(x,y) x.^2+y.^2-P_Sc.inf,fBox,'-.', ...
                    'color',cList(3,:),'linewidth',lineWidthS);
fimplicit(@(x,y) x.^2+y.^2-P_An.inf,fBox,'-.', ...
                    'color',cList(4,:),'linewidth',lineWidthS);
fimplicit(@(x,y) x.^2+y.^2-P_Ar.inf,fBox,':', ...
                    'color',cList(2,:),'linewidth',lineWidthS);
fimplicit(@(x,y) x.^2+y.^2-P_g.inf,fBox,'b--','linewidth',lineWidthS);
fimplicit(@(x,y) x.^2+y.^2-P_x.inf,fBox,'r:','linewidth',lineWidthS);

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
text(B_r.real.inf-0.03,B_r.imag.mid,'$\underline{|B^I|}$',...
                    'HorizontalAlignment','right', 'Interpreter','latex')
text(B_r.real.sup+0.03,B_r.imag.mid,'$\overline{|B^I|}$',...
                    'HorizontalAlignment','left', 'Interpreter','latex')

% Add zoom window for the operand interval
axes('position',[0.09,0.18,0.3,0.3]); hold on; box on
set(gca,'Color',ones(3,1)*0.95)
EA_g.plot('b','linewidth',lineWidthM);
EA_x.plot('r','linewidth',lineWidthS);
EA_An.plot('color',cList(4,:),'linewidth',lineWidthS);
axis equal
set(gca, 'XAxisLocation', 'top')
m = 3;
maxWidth = max([EA_c(m).Real.Width , EA_c(m).Imag.Width]);
xlim(EA_c(m).Real.Midpoint + [-1 1]*maxWidth/2)
ylim(EA_c(m).Imag.Midpoint + [-1 1]*maxWidth/2)
xticks([]);
yticks([]);
set(gca, 'YAxisLocation', 'right')
text(real(EA_nom(m)),imag(EA_nom(m)), ...
    ['$E_{' num2str(m) '}^I A_{' num2str(m) '}^I$'], ...
                'HorizontalAlignment','center', 'Interpreter','latex')


% Add zoom window for the infimum
axes('position',[0.73,0.178,0.14,0.3]); hold on; box on
set(gca,'Color',ones(3,1)*0.95)
fimplicit(@(x,y) x.^2+y.^2-P_Sc.inf,fBox,'-.', ...
                    'color',cList(3,:),'linewidth',lineWidthS);
fimplicit(@(x,y) x.^2+y.^2-P_An.inf,fBox,'-.', ...
                    'color',cList(4,:),'linewidth',lineWidthS);
fimplicit(@(x,y) x.^2+y.^2-P_Ar.inf,fBox,':', ...
                    'color',cList(2,:),'linewidth',lineWidthS);
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
set(gca,'Color',ones(3,1)*0.95)
fimplicit(@(x,y) x.^2+y.^2-P_Sc.sup,fBox,'-.', ...
                        'color',cList(3,:),'linewidth',lineWidthS);
fimplicit(@(x,y) x.^2+y.^2-P_An.sup,fBox,'-.', ...
                        'color',cList(4,:),'linewidth',lineWidthS);
fimplicit(@(x,y) x.^2+y.^2-P_Ar.sup,fBox,':', ...
                        'color',cList(2,:),'linewidth',lineWidthS);
fimplicit(@(x,y) x.^2+y.^2-P_Ar.sup,fBox,'-.', ...
                        'color','y','linewidth',lineWidthS);
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
                    ', T^{(a)}=' join(compose("%0.0f",1e3*T_x),'+')  , ...
                    '=' num2str(sum(T_x,1)*1e3,'%0.0f') 'ms'],'')...
            join(['\color{blue}\tau^{(g)}=' num2str(100*tau_g,3) '%' ...
                    ', T^{(g)}=' join(compose("%0.0f",1e3*T_g),'+')  , ...
                    '=' num2str(sum(T_g,1)*1e3,'%0.0f') 'ms'],'')...
            join(['\color[rgb]{' num2str(cList(4,:)) '}\tau^{(C)}=' ...
                        num2str(100*tau_An,3) '%' ...
                    ', T^{(C)}=' join(compose("%0.0f",1e3*T_An),'+')  , ...
                    '=' num2str(sum(T_An,1)*1e3,'%0.0f') 'ms'],'')...
            join(['\color[rgb]{' num2str(cList(2,:)) '}\tau^{(A)}=' ...
                        num2str(100*tau_Ar,3) '%' ...
                    ', T^{(A)}=' join(compose("%0.0f",1e3*T_Ar),'+')  , ...
                    '=' num2str(sum(T_Ar,1)*1e3,'%0.0f') 'ms'],'')...
            join(['\color[rgb]{' num2str(cList(3,:)) '}\tau^{(S)}=' ...
                        num2str(100*tau_Sc,3) '%' ...
                    ', T^{(S)}=' join(compose("%0.0f",1e3*T_Sc),'+')  , ...
                    '=' num2str(sum(T_Sc,1)*1e3,'%0.0f') 'ms'],'')};
annotation('textbox',[0.152,0.555,0.280 0.367],'String',annotText, ...
           'BackgroundColor','w','VerticalAlignment','top','fontSize',25,...
           'HorizontalAlignment','center','FitBoxToText','on');

annotation('textbox',[0.18 0.52 0.1 0.1],'String','d)',...
                        'FontSize',60,'EdgeColor','none')


