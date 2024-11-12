clear
% close all

% Set parameters
M = 8;
K = 10;
SNR = 20;
w = chebwin(M,SNR) ; w = w / sum(w);
phi = 0.03;
theta = linspace(-pi/2,pi/2,K)';

% Set error
ampErr = 0.1;
phaErr = deg2rad(10);
gamma = 0.05;
alpha = ones(M,1) * sqrt(ampErr^2 + sin(phaErr)^2)/2;
beta = ones(M-1,1)*gamma;

% Initialize beampattern intervals
B_nom = zeros(K,1);
B_r(K,4) = ciat.RectangularInterval;
B_c(K,4) = ciat.CircularInterval;
B_g(K,4) = ciat.PolygonalInterval;
B_x(K,4) = ciat.PolyarxInterval;
B_a(K,4) = ciat.PolyarcularInterval;


% Calculate beampattern intervals for each angle
tRect = zeros(1,4);
tCirc = zeros(1,4);
tPgon = zeros(1,4);
tParc = zeros(1,4);

for k=1:K
    v = exp(1j*pi*sin(theta(k))*(1:M)');
    
    % Calculate error intervals
    E_amp = ciat.PolarInterval(w.*(1+[-1 1].*alpha*0.75) , angle(v')' * [1 1]);
    E_pha = ciat.PolarInterval(w.* [1,1] , angle(v')' + [-1 1].*phi*pi);
    E_cal = ciat.PolarInterval(w.*(1+[-1 1].*alpha*0.75) , angle(v')' + [-1 1].*phi*pi);
    A_mcp = ciat.CircularInterval(zeros(M,1) , [0;beta]+[beta;0]);

    % Calculate element intervals
    AF_nom = w .* v;
    AF_a = [ciat.PolyarcularInterval([E_amp, E_pha, E_cal]) , ...
            ciat.PolyarcularInterval(E_cal , ones(M,1) + A_mcp)];
    tic
    AF_r = ciat.RectangularInterval(AF_a);
    tRect(1) = tRect(1) + toc;
    tic;
    AF_c = ciat.CircularInterval(AF_a);
    tCirc(1) = tCirc(1) + toc;
    tic
    AF_g = [ciat.PolygonalInterval([E_amp, E_pha, E_cal]) , ...
            ciat.PolygonalInterval(E_cal , ones(M,1) + A_mcp)];
    tPgon(1) = tPgon(1) + toc;
    tic
    AF_x = [ciat.PolyarxInterval([E_amp, E_pha, E_cal]) , ...
            ciat.PolyarxInterval(E_cal , ones(M,1) + A_mcp)];
    tParc(1) = tParc(1) + toc;

    % Calculate beampattern intervals
    B_nom(k,:) = sum(AF_nom);
    tic
    B_r(k,:) = sum(AF_r,1);
    tRect(2) = tRect(2) + toc;
    tic
    B_c(k,:) = sum(AF_c,1);
    tCirc(2) = tCirc(2) + toc;
    tic;
    B_g(k,:) = sum(AF_g,1);
    tPgon(2) = tPgon(2) + toc;
    tic
    B_x(k,:) = sum(AF_x,1);
    tParc(2) = tParc(2) + toc;
end

% Calculate power pattern intervals
P_nom = abs(B_nom).^2;
tic
P_r = abs(B_r).^2; P_r.Infimum(P_r.inf==0)=eps;
tRect(3) = tRect(3) + toc;
tic
P_c = abs(B_c).^2; P_c.Infimum(P_c.inf==0)=eps;
tCirc(3) = tCirc(3) + toc;
tic
P_g = abs(B_g).^2; P_g.Infimum(P_g.inf==0)=eps;
tPgon(3) = tPgon(3) + toc;
tic
P_x = abs(B_x).^2; P_x.Infimum(P_x.inf==0)=eps;
tParc(3) = tParc(3) + toc;

% Calculate pattern tolerance
D_r = sum(P_r.width,1) / sum(P_nom);
D_c = sum(P_c.width,1) / sum(P_nom);
D_g = sum(P_g.width,1) / sum(P_nom);
D_x = sum(P_x.width,1) / sum(P_nom);




%% Plot


% Set parameters
kSegments = [0 0.22 0.5 0.77];
xRange = [-90 90];
yRange = [-60 5];

figure(1);clf;hold on
set(gca,'DefaultLineLineWidth',3)

% Set axis
thAxis = rad2deg(theta);

% Plot nominal response
p0 = plot(thAxis,db(P_nom),'k-','DisplayName','Nominal');

% Amplitude error
for seg = 1:4
    % Select k range
    if seg < 4
        kRange = (round(K*kSegments(seg))+1):(round(K*kSegments(seg+1))+ 1);
    else
        kRange = (round(K*kSegments(seg))+1):K;
    end
    % kRange = 1:K;

    % Plot vertical splitting line
    if seg>1
        plot(thAxis(kRange(1))*[1,1],yRange,'k:')
    end

    % Rectangular
    p1 = plot(thAxis(kRange),db(P_r(kRange,seg).sup),'c-','DisplayName','Rectangle (sup)');
    p2 = plot(thAxis(kRange),db(P_r(kRange,seg).inf),'c-.','DisplayName','Rectangle (inf)');
    
    % Circular
    p3 = plot(thAxis(kRange),db(P_c(kRange,seg).sup),'m-','DisplayName','Circle (sup)');
    p4 = plot(thAxis(kRange),db(P_c(kRange,seg).inf),'m-.','DisplayName','Circle (inf)');

    % Polygonal
    p5 = plot(thAxis(kRange),db(P_g(kRange,seg).sup),'b-','DisplayName','Polygon (sup)');
    p6 = plot(thAxis(kRange),db(P_g(kRange,seg).inf),'b--','DisplayName','Polygon (inf)');

    % Polyarcular (convex)
    p7 = plot(thAxis(kRange),db(P_x(kRange,seg).sup),'r--','DisplayName','Polyarc (sup)');
    p8 = plot(thAxis(kRange),db(P_x(kRange,seg).inf),'r:','DisplayName','Polyarc (inf)');
end

% Pattern tolerance index
    % Amplitude errors
annotText = {   'Amplitude',...
                ['$\Delta^{\mathrm{(r)}}=' num2str(D_r(1),'%4.4f') '$'],...
                ['$\Delta^{\mathrm{(c)}}=' num2str(D_c(1),'%4.4f') '$'],...
                ['$\Delta^{\mathrm{(g)}}=' num2str(D_g(1),'%4.4f') '$'],...
                ['$\Delta^{\mathrm{(a)}}=' num2str(D_x(1),'%4.4f') '$']};
annotation('textbox',[0.14 0.70 .11 0.22],'String',annotText, ...
           'BackgroundColor','w',...
           'HorizontalAlignment','center', ...
           'Interpreter','latex');
    % Amplitude errors
annotText = {   'Phase',...
                ['$\Delta^{\mathrm{(r)}}=' num2str(D_r(2),'%4.4f') '$'],...
                ['$\Delta^{\mathrm{(c)}}=' num2str(D_c(2),'%4.4f') '$'],...
                ['$\Delta^{\mathrm{(g)}}=' num2str(D_g(2),'%4.4f') '$'],...
                ['$\Delta^{\mathrm{(a)}}=' num2str(D_x(2),'%4.4f') '$']};
annotation('textbox',[0.31 0.70 .11 0.22],'String',annotText, ...
           'BackgroundColor','w',...
           'HorizontalAlignment','center', ...
           'Interpreter','latex');
    % Calibration errors
annotText = {   'Calibration',...
                ['$\Delta^{\mathrm{(r)}}=' num2str(D_r(3),'%4.4f') '$'],...
                ['$\Delta^{\mathrm{(c)}}=' num2str(D_c(3),'%4.4f') '$'],...
                ['$\Delta^{\mathrm{(g)}}=' num2str(D_g(3),'%4.4f') '$'],...
                ['$\Delta^{\mathrm{(a)}}=' num2str(D_x(3),'%4.4f') '$']};
annotation('textbox',[0.61 0.70 .11 0.22],'String',annotText, ...
           'BackgroundColor','w',...
           'HorizontalAlignment','center', ...
           'Interpreter','latex');
    % Mutual coupling 
annotText = {   'Coupling',...
                ['$\Delta^{\mathrm{(r)}}=' num2str(D_r(4),'%4.4f') '$'],...
                ['$\Delta^{\mathrm{(c)}}=' num2str(D_c(4),'%4.4f') '$'],...
                ['$\Delta^{\mathrm{(g)}}=' num2str(D_g(4),'%4.4f') '$'],...
                ['$\Delta^{\mathrm{(a)}}=' num2str(D_x(4),'%4.4f') '$']};
annotation('textbox',[0.76 0.70 .11 0.22],'String',annotText, ...
           'BackgroundColor','w',...
           'HorizontalAlignment','center', ...
           'Interpreter','latex');

% Format figure
xlim(xRange)
ylim(yRange)
xlabel('Beam angle [deg]')
ylabel('Power')
%legend([p0 p1 p2 p3 p4 p5 p6 p7 p8],'Position',[0.517 0.33 0 0])
fontsize(25,'point')