clear
% close all

% Set parameters
M = 8;
K = 101;
w = taylorwin(M,3,-20); w = w / sum(w);
theta = linspace(-pi/2,pi/2,K)';

% Set errors
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

% Initaliz power pattern intervals
P_Hu(K,1) = ciat.RealInterval;
P_He(K,1) = ciat.RealInterval;
P_Sch(K,1) = ciat.RealInterval;
P_Ans(K,1) = ciat.RealInterval;

% Calculate beampattern intervals for each angle
T_r = zeros(3,4);
T_c = zeros(3,4);
T_g = zeros(3,4);
T_x = zeros(3,4);
T_He = 0;
T_Hu = 0;
T_Sc = 0;
T_An = 0;

for k=1:K
    phi = ((0:M-1)-(M-1)/2)' * pi*sin(theta(k));
    v = exp(1j*phi);
    
    % Calculate error intervals
    E_amp = ciat.PolarInterval(w.*(1+[-1 1].*ampErr/2) , ...
                                phi * [1 1]);
    E_pha = ciat.PolarInterval(w.* [1,1] , ...
                                phi + [-1 1].*phaErr/2);
    E_cal = ciat.PolarInterval(w.*(1+[-1 1].*ampErr/2) , ...
                                phi + [-1 1].*phaErr/2);
    E_int = [E_amp, E_pha, E_cal];
    A_int = ciat.CircularInterval(zeros(M,1) , [0;beta]+[beta;0]);
    AF_a = [ciat.PolyarcularInterval([E_amp, E_pha, E_cal]) , ...
            ciat.PolyarcularInterval(E_cal , ones(M,1) + A_int)];

    % Calculate nominal beampattern
    B_nom(k,:) = w' * v;
    
    % Calculate beampattern intervals
    for s=1:4
        % Rectangular
        tic;
        if s < 4
            AF_r = ciat.RectangularInterval(E_int(:,s));
        else
            E_r = ciat.RectangularInterval(E_cal);
            A_r = ciat.RectangularInterval(ones(M,1)+A_int);
            AF_r = E_r .* A_r;
        end
        T_r(1,s) = T_r(1,s) + toc;
        tic;
        B_r(k,s) = sum(AF_r);
        T_r(2,s) = T_r(2,s) + toc;

        % Circular
        tic;
        if s < 4
            AF_c = ciat.CircularInterval(E_int(:,s));
        else
            E_c = ciat.CircularInterval(E_cal);
            A_c = ones(M,1) + A_int;
            AF_c = E_c .* A_c;
        end
        T_c(1,s) = T_c(1,s) + toc;
        tic;
        B_c(k,s) = sum(AF_c);
        T_c(2,s) = T_c(2,s) + toc;

        % Polygonal
        tic
        if s < 4
            AF_g = ciat.PolygonalInterval(E_int(:,s));
        else
            AF_g = ciat.PolygonalInterval(E_cal , ones(M,1) + A_int);
        end
        T_g(1,s) = T_g(1,s) + toc;
        tic;
        B_g(k,s) = sum(AF_g);
        T_g(2,s) = T_g(2,s) + toc;

        % Polyarcular
        tic
        if s < 4
            AF_x = ciat.PolyarxInterval(E_int(:,s));
        else
            AF_x = ciat.PolyarxInterval(E_cal , ones(M,1) + A_int);
        end
        T_x(1,s) = T_x(1,s) + toc;
        tic;
        B_x(k,s) = sum(AF_x);
        T_x(2,s) = T_x(2,s) + toc;
    end

    % Hu's Taylor based approximation
        tic;
        A = w * (1 + ciat.RealInterval(-ampErr/2,ampErr/2));
        A_R = sum(A .* cos(phi));
        A_I = sum(A .* sin(phi));
        P_Hu_d = 2*cos(phi) .* A_R + 2*sin(phi) .* A_I;
        P_Hu_0 = abs(sum(w .* exp(1j*phi)))^2;
        P_Hu_U = P_Hu_0 + sum( abs(P_Hu_d.sup) .* A.width/2 );
        P_Hu_I = P_Hu_0 - sum( abs(P_Hu_d.inf) .* A.width/2 );
        P_Hu_I (P_Hu_I<0) = 0;
        P_Hu(k) = ciat.RealInterval(P_Hu_I,P_Hu_U);
        T_Hu = T_Hu + toc;
        
        % He's matrix method
        tic;
        A = w * (1 + ciat.RealInterval(-ampErr/2,ampErr/2));
        a_mid = A.mid';
        a_rad = A.width'/2;
        Theta = cos(phi - phi');
        P_He_mid = abs(sum(w .* exp(1j*phi)))^2;
        P_He_sup = P_He_mid + 2 * abs(a_mid * Theta) * a_rad' + ...
                              a_rad * abs(Theta)*a_rad';
        P_He_inf = P_He_mid - 2 * abs(a_mid * Theta) * a_rad';
        P_He_inf(P_He_inf<0) = 0;
        P_He(k) = ciat.RealInterval(P_He_inf,P_He_sup);
        T_He = T_He + toc;

        % Schmid method
        tic
        wL2 = sqrt(sum(abs(w)*.2));
        maxDeltaB = sqrt(M) * wL2 * (w'*alpha + w(1:end-1)'*beta + ...
                                     w(2:end)'*beta);
        B_Sch = ciat.CircularInterval(B_nom(k),maxDeltaB);
        P_Sch(k) = abs(B_Sch).^2;
        T_Sc = T_Sc + toc;
        
        % Anselmi method
        tic;
        diagC0 = diag(alpha);
        diagC1 = diag(beta,1);
        diagCm1 = diag(beta,-1);
        Ca = ciat.CircularInterval(zeros(M) , diagC0);
        Cb = ciat.CircularInterval(zeros(M) , diagC1 + diagCm1);
        B_Ans = w' * (Ca + Cb + eye(M)) * v;
        P_Ans(k) = abs(B_Ans).^2;
        T_An = T_An + toc;
end

% Calculate power pattern intervals
P_nom = abs(B_nom).^2;
tic
P_r = abs(B_r).^2; P_r.Infimum(P_r.inf==0)=eps;
T_r(3) = T_r(3) + toc;
tic
P_c = abs(B_c).^2; P_c.Infimum(P_c.inf==0)=eps;
T_c(3) = T_c(3) + toc;
tic
P_g = abs(B_g).^2; P_g.Infimum(P_g.inf==0)=eps;
T_g(3) = T_g(3) + toc;
tic
P_x = abs(B_x).^2; P_x.Infimum(P_x.inf==0)=eps;
T_x(3) = T_x(3) + toc;

% Calculate pattern tolerance
D_r = sum(P_r.width,1) / sum(P_nom);
D_c = sum(P_c.width,1) / sum(P_nom);
D_g = sum(P_g.width,1) / sum(P_nom);
D_x = sum(P_x.width,1) / sum(P_nom);




%% Plot


% Set parameters
kSegments = [0 0.23 0.5 0.76];
xRange = [-90 90];
yRange = [-60 5];

figure(2);clf;hold on
set(gca,'DefaultLineLineWidth',3)

% Set axis
thAxis = rad2deg(theta);

% Plot nominal response
p0 = plot(thAxis,db(P_nom),'k-','DisplayName','Nominal');

% Amplitude error
for s = 1%:4
    % Select k range
    if s < 4
        kRange = (round(K*kSegments(s))+1):(round(K*kSegments(s+1))+ 1);
    else
        kRange = (round(K*kSegments(s))+1):K;
    end
    kRange = 1:K;

    % Plot vertical splitting line
    if s>1
        plot(thAxis(kRange(1))*[1,1],yRange,'k:')
    end

    % Rectangular
    pRs = plot(thAxis(kRange),db(P_r(kRange,s).sup),'c-', ...
                                    'DisplayName','Rectangle (sup)');
    pRi = plot(thAxis(kRange),db(P_r(kRange,s).inf),'c-.','DisplayName','Rectangle (inf)');
    
    % Circular
    pCs = plot(thAxis(kRange),db(P_c(kRange,s).sup),'m-', ...
                                    'DisplayName','Circle (sup)');
    pCi = plot(thAxis(kRange),db(P_c(kRange,s).inf),'m-.', ...
                                    'DisplayName','Circle (inf)');

    % Polygonal
    pGs = plot(thAxis(kRange),db(P_g(kRange,s).sup),'b-', ...
                                    'DisplayName','Polygon (sup)');
    pGi = plot(thAxis(kRange),db(P_g(kRange,s).inf),'b--', ...
                                    'DisplayName','Polygon (inf)');

    % Polyarcular (convex)
    pAs = plot(thAxis(kRange),db(P_x(kRange,s).sup),'r--', ...
                                    'DisplayName','Polyarc (sup)');
    pAi = plot(thAxis(kRange),db(P_x(kRange,s).inf),'r:', ...
                                    'DisplayName','Polyarc (inf)');

    if s == 1
        % Hu
        pHuS = plot(thAxis(kRange),db(P_Hu(kRange).sup),'--', ...
                            'color',0.7*ones(1,3),'DisplayName','Hu (sup)');
        pHuI = plot(thAxis(kRange),db(P_Hu(kRange).inf),'-', ...
                            'color',0.7*ones(1,3),'DisplayName','Hu (inf)');
    
        % He
        pHeS = plot(thAxis(kRange),db(P_He(kRange).sup),'--', ...
                            'color',0.5*ones(1,3),'DisplayName','He (sup)');
        pHeI = plot(thAxis(kRange),db(P_He(kRange).inf),'-', ...
                            'color',0.5*ones(1,3),'DisplayName','He (inf)');
    end

    if s == 4
        % Schmid
        pScS = plot(thAxis(kRange),db(P_Sch(kRange).sup),'--', ...
                            'color',0.7*ones(1,3),'DisplayName','Schmid (sup)');
        pScI = plot(thAxis(kRange),db(P_Sch(kRange).inf),'-', ...
                            'color',0.7*ones(1,3),'DisplayName','Schmid (inf)');
    
        % Anselmi
        pAnS = plot(thAxis(kRange),db(P_Ans(kRange).sup),'--', ...
                            'color',0.5*ones(1,3),'DisplayName','Anselmi (sup)');
        pAnI = plot(thAxis(kRange),db(P_Ans(kRange).inf),'-', ...
                            'color',0.5*ones(1,3),'DisplayName','Anselmi (inf)');
    end
end

% Pattern tolerance index
    % Amplitude errors
annotText = {   'Amp. error',...
                ['\color{cyan}\Delta^{(r)}=' num2str(D_r(1),'%4.4f')],...
                ['\color{magenta}\Delta^{(c)}=' num2str(D_c(1),'%4.4f')],...
                ['\color{blue}\Delta^{(g)}=' num2str(D_g(1),'%4.4f')],...
                ['\color{red}\Delta^{(a)}=' num2str(D_x(1),'%4.4f')]};
annotation('textbox',[0.14 0.70 .11 0.22],'String',annotText, ...
           'BackgroundColor','w',...
           'HorizontalAlignment','center');
    % Phase errors
annotText = {   'Phase error',...
                ['\color{cyan}\Delta^{(r)}=' num2str(D_r(2),'%4.4f')],...
                ['\color{magenta}\Delta^{(c)}=' num2str(D_c(2),'%4.4f')],...
                ['\color{blue}\Delta^{(g)}=' num2str(D_g(2),'%4.4f')],...
                ['\color{red}\Delta^{(a)}=' num2str(D_x(2),'%4.4f')]};
annotation('textbox',[0.33 0.70 .11 0.22],'String',annotText, ...
           'BackgroundColor','w',...
           'HorizontalAlignment','center');
    % Calibration errors
annotText = {   'A & P errors',...
                ['\color{cyan}\Delta^{(r)}=' num2str(D_r(3),'%4.4f') ''],...
                ['\color{magenta}\Delta^{(c)}=' num2str(D_c(3),'%4.4f') ''],...
                ['\color{blue}\Delta^{(g)}=' num2str(D_g(3),'%4.4f') ''],...
                ['\color{red}\Delta^{(a)}=' num2str(D_x(3),'%4.4f') '']};
annotation('textbox',[0.61 0.70 .11 0.22],'String',annotText, ...
           'BackgroundColor','w',...
           'HorizontalAlignment','center');
    % Mutual coupling 
annotText = {   'A&P & Coupl.',...
                ['\color{cyan}\Delta^{(r)}=' num2str(D_r(4),'%4.4f')],...
                ['\color{magenta}\Delta^{(c)}=' num2str(D_c(4),'%4.4f')],...
                ['\color{blue}\Delta^{(g)}=' num2str(D_g(4),'%4.4f')],...
                ['\color{red}\Delta^{(a)}=' num2str(D_x(4),'%4.4f')]};
annotation('textbox',[0.76 0.70 .12 0.22],'String',annotText, ...
           'BackgroundColor','w',...
           'HorizontalAlignment','center');

% Format figure
xlim(xRange)
ylim(yRange)
xlabel('Beam angle [deg]')
ylabel('Power')
%legend([p0 p1 p2 p3 p4 p5 p6 p7 p8],'Position',[0.517 0.33 0 0])
fontsize(25,'point')