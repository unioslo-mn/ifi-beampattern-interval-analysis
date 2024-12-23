    clear
% close all

%% Amplitude errors

% Set parameters
conf = getFig5conf(1);
M = conf.M;
w = conf.w;
theta = conf.theta;
ampErr = conf.ampErr;
nT =1% conf.nT;

% Calculate nominal element phase 
phi = ((0:M-1)-(M-1)/2)' * pi*sin(theta);

% Calculate amplitude interval
A = w * (1 + ciat.RealInterval(-ampErr/2,ampErr/2));

% Calculate nominal beampattern and power pattern
EA_nom = w .* exp(1j*phi);
B_nom = sum(EA_nom);
P_nom = abs(B_nom)^2;

% Measure running time on many occurences and average the result
T_g    = zeros(3,1);
T_r    = zeros(3,1);
T_x    = zeros(3,1);
T_He   = zeros(3,1);
T_Hu   = zeros(3,1);
T_Ar   = zeros(3,1);


for iT = 1:nT
    fprintf("%0i,",iT)
    % Define and cast intervals
    EA_p = ciat.PolarInterval(A,ciat.RealInterval(phi));
    tic;EA_g = ciat.PolygonalInterval(EA_p,'tolerance',conf.tol);T_g(1) = T_g(1)+toc;
    tic;EA_x = ciat.PolyarxInterval(EA_p);T_x(1) = T_x(1)+toc;
    tic; EA_r = ciat.RectangularInterval(EA_g); T_r(1) = T_r(1)+toc;
    EA_a = ciat.PolyarcularInterval(EA_p);
    
    % Sum intervals
    tic;B_r = sum(EA_r);T_r(2)=T_r(2)+toc;
    tic;B_g = sum(EA_g);T_g(2)=T_g(2)+toc;
    tic;B_x = sum(EA_x);T_x(2)=T_x(2)+toc;
    B_a = sum(EA_a);
    
    % Get power intervals
    tic;P_g = abs(B_g)^2;T_g(3)=T_g(3)+toc;
    tic;P_x = abs(B_x)^2;T_x(3)=T_x(3)+toc;
    tic;P_r = abs(B_r)^2;T_r(3)=T_r(3)+toc;
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
    T_Hu(3) = T_Hu(3)+toc;
    
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
    T_He(3) = T_He(3)+ toc;

    % % Arnestad's method
    % tic;
    % P_Ar_sup = ( abs(B_nom) + sqrt((ampErr/2)^2) )^2;
    % P_Ar_inf = ( abs(B_nom) - sqrt((ampErr/2)^2) )^2;
    % P_Ar_inf(P_Ar_inf<0) = 0;
    % P_Ar = ciat.RealInterval(P_Ar_inf,P_Ar_sup);
    % T_Ar(3) = T_Ar(3)+ toc;
end

% Calculate average time
T_g    = T_g/nT;
T_He   = T_He/nT;
T_Hu   = T_Hu/nT;
T_r    = T_r/nT;
T_x    = T_x/nT;
T_Ar    = T_Ar/nT;

% Calculate tightness
tau_x = P_a.Width ./ P_x.Width;
tau_g = P_a.Width ./ P_g.Width;
tau_r = P_a.Width ./ P_r.Width;
tau_Hu = P_a.Width ./ (P_a.Width + abs(P_Hu.sup-P_x.sup) + ...
                                   abs(P_Hu.inf-P_x.inf));
tau_He = P_a.Width ./ P_He.Width;
% tau_Ar = P_a.Width ./ P_Ar.Width;

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

% Grey boxes for the zoom windows
rectangle('Position',[-0.1456  -0.0226  0.0139   0.0139],...
    'EdgeColor','none','FaceColor',ones(3,1)*0.95);
rectangle('Position',[0.1246  -0.0033  0.0063   0.0066],...
    'EdgeColor','none','FaceColor',ones(3,1)*0.95);
rectangle('Position',[0.0505  -0.0086 0.0164 0.0172],...
    'EdgeColor','none','FaceColor',ones(3,1)*0.95);


% Plot operand intervals
EA_g.plot('b','linewidth',lineWidthS);
EA_x.plot('r','linewidth',lineWidthS);
EA_r.plot('color',cList(1,:),'linewidth',lineWidthS);

% Plot sum intervals
lC = B_x.plot('r','linewidth',lineWidthS,'DisplayName','Polyarcular');
lB = B_g.plot('b','linewidth',lineWidthL,'DisplayName','Polygonal');
lA = B_r.plot('color',cList(1,:),'linewidth',lineWidthL,'DisplayName','Rectangular');

% Interval label
for m = 1:M
    text(real(EA_nom(m))+0.01,imag(EA_nom(m))-0.005, ...
        ['$E_{' num2str(m) '}^I A_{' num2str(m) '}^I$'], ...
                'HorizontalAlignment','left', 'Interpreter','latex')
end
text(B_r.real.mid,B_r.imag.mid,'$B^I$',...
                    'HorizontalAlignment','center', 'Interpreter','latex')
text(B_r.real.inf-0.01,B_r.imag.mid,'$\underline{|B^I|}$',...
                    'HorizontalAlignment','right', 'Interpreter','latex')
text(B_r.real.sup+0.01,B_r.imag.mid,'$\overline{|B^I|}$',...
                    'HorizontalAlignment','left', 'Interpreter','latex')

% Plot infimum
fimplicit(@(x,y) x.^2+y.^2-P_r.inf,fBox,'-.',...
            'color',cList(1,:),'linewidth',lineWidthM,'DisplayName','Anselmi');
% fimplicit(@(x,y) x.^2+y.^2-P_Ar.inf,fBox,'-.',...
%             'color',cList(2,:),'linewidth',lineWidthM,'DisplayName','Arnestad');
fimplicit(@(x,y) x.^2 + y.^2 - P_Hu.inf,fBox,':', ...
            'color',cList(3,:),'linewidth',lineWidthM,'DisplayName','Hu');
fimplicit(@(x,y) x.^2 + y.^2 - P_He.inf,fBox,':',...
            'color',cList(4,:),'linewidth',lineWidthM,'DisplayName','He');
fimplicit(@(x,y) x.^2+y.^2-P_g.inf,fBox,'b--','linewidth',lineWidthM,...
                                            'DisplayName','Tenuti');
l3 = fimplicit(@(x,y) x.^2+y.^2-P_x.inf,fBox,'r:', ...
                                'linewidth',lineWidthM,'DisplayName','Geréb');

% Plot supremum
fimplicit(@(x,y) x.^2+y.^2-P_r.sup,fBox,'-.',...
            'color',cList(1,:),'linewidth',lineWidthM)
% fimplicit(@(x,y) x.^2+y.^2-P_Ar.sup,fBox,'-.',...
%             'color',cList(2,:),'linewidth',lineWidthM)
fimplicit(@(x,y) x.^2 + y.^2 - P_Hu.sup,fBox,':', ...
            'color',cList(3,:), 'linewidth',lineWidthM)
fimplicit(@(x,y) x.^2 + y.^2 - P_He.sup,fBox,':',...
            'color',cList(4,:),'linewidth',lineWidthM)
fimplicit(@(x,y) x.^2+y.^2-P_g.sup,fBox,'b--',...
                                'linewidth',lineWidthM)
l3 = fimplicit(@(x,y) x.^2+y.^2-P_x.sup,fBox,'r:', ...
                                'linewidth',lineWidthM);

% Set figure limits
xlim(xL); ylim(yL);

% Add zoom window for the operand interval
axes('position',[0.09,0.18,0.3,0.3]); hold on; box on
set(gca,'Color',ones(3,1)*0.95)
EA_p.plot('b','linewidth',lineWidthL);
EA_r.plot('color',cList(1,:),'linewidth',lineWidthL);
EA_x.plot('r','linewidth',lineWidthS);
axis equal
set(gca, 'XAxisLocation', 'top')
m = 3;
maxWidth = max([EA_r(m).Real.Width , EA_r(m).Imag.Width]);
xlim(EA_r(m).Real.Midpoint + [-1 1]*maxWidth/2)
ylim(EA_r(m).Imag.Midpoint + [-1 1]*maxWidth/2)
xticks([]);
yticks([]);
set(gca, 'YAxisLocation', 'right')
text(real(EA_nom(m)),imag(EA_nom(m))+0.003, ...
            ['$E_{' num2str(m) '}^I A_{' num2str(m) '}^I$'], ...
                'HorizontalAlignment','right', 'Interpreter','latex')

% Add zoom window for the infimum
axes('position',[0.73,0.178,0.14,0.3]); hold on; box on
set(gca,'Color',ones(3,1)*0.95)
fimplicit(@(x,y) x.^2+y.^2-P_r.inf,'-.', ...
                        'color',cList(1,:),'linewidth',lineWidthL)
% fimplicit(@(x,y) x.^2+y.^2-P_Ar.inf,'-.', ...
%                         'color',cList(2,:),'linewidth',lineWidthL)
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
set(gca,'Color',ones(3,1)*0.95)
fimplicit(@(x,y) x.^2+y.^2-P_r.sup,'-.', ...
                    'color',cList(1,:),'linewidth',lineWidthL)
% fimplicit(@(x,y) x.^2+y.^2-P_Ar.sup,'-.', ...
%                     'color',cList(2,:),'linewidth',lineWidthL)
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
                    ', T^{(a)}=' join(compose("%0.0f",1e3*T_x),'+'), ...
                    '=' num2str(sum(T_x,1)*1e3,'%0.0f') 'ms'],'')...
            join(['\color{blue}\tau^{(g)}=' num2str(100*tau_g,3) '%' ...
                    ', T^{(g)}=' join(compose("%0.0f",1e3*T_g),'+'), ...
                    '=' num2str(sum(T_g,1)*1e3,'%0.0f') 'ms'],'')...
            join(['\color[rgb]{' num2str(cList(1,:)) '}\tau^{(r)}=' ...
                        num2str(100*tau_r,3) '%' ...
                    ', T^{(r)}=' join(compose("%0.0f",1e3*T_r),'+'), ...
                    '=' num2str(sum(T_r,1)*1e3,'%0.0f') 'ms'],'')...
            join(['\color[rgb]{' num2str(cList(4,:)) '}\tau^{(M)}=' ...
                        num2str(100*tau_He,3) '%' ...
                    ', T^{(M)}=' join(compose("%0.0f",1e3*T_He),'+'), ...
                    '=' num2str(sum(T_He,1)*1e3,'%0.0f') 'ms'],'')...
            join(['\color[rgb]{' num2str(cList(3,:)) '}\tau^{(T)}=' ...
                        num2str(100*tau_Hu,3) '%' ...
                    ', T^{(T)}=' join(compose("%0.0f",1e3*T_Hu),'+'), ...
                    '=' num2str(sum(T_Hu,1)*1e3,'%0.0f') 'ms'],'')...
            % join(['\color[rgb]{' num2str(cList(2,:)) '}\tau^{(A)}=' ...
            %             num2str(100*tau_Ar,3) '%' ...
            %         ', T^{(M)}=' join(compose("%0.0f",1e3*T_Ar),'+'), ...
            %         '=' num2str(sum(T_Ar,1)*1e3,'%0.0f') 'ms'],'')  ...
                    };
annotation('textbox',[0.152,0.555,0.280 0.367],'String',annotText, ...
           'BackgroundColor','w','VerticalAlignment','top','fontSize',25,...
           'HorizontalAlignment','center','FitBoxToText','on');


annotation('textbox',[0.18 0.52 0.1 0.1],'String','a)',...
                        'FontSize',60,'EdgeColor','none')



