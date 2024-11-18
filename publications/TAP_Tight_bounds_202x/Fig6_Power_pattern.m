clear
% close all

% Set parameters
conf = getFig5conf(1);
M = conf.M;
w = conf.w;
ampErr = conf.ampErr;
phaErr = conf.phaErr;
mtlCpl = conf.mtlCpl;
K = 1001;
theta = linspace(-pi/2,pi/2,K)';
kSeg = [0.23 0.5 0.77]; kSeg = [1 ceil(kSeg * K) K];

% Generate error vectors
alpha = ones(M,1) * sqrt(ampErr^2 + sin(phaErr)^2)/2;
beta = ones(M-1,1)*mtlCpl;

% Initialize beampattern intervals
B_nom = zeros(K,1);
B_r(K,4) = ciat.RectangularInterval;
B_c(K,4) = ciat.CircularInterval;
B_g(K,4) = ciat.PolygonalInterval;
B_x(K,4) = ciat.PolyarxInterval;
B_a(K,4) = ciat.PolyarcularInterval;

% Initaliz power pattern intervals
P_r(K,4) = ciat.RealInterval;
P_c(K,4) = ciat.RealInterval;
P_g(K,4) = ciat.RealInterval;
P_x(K,4) = ciat.RealInterval;
P_Hu(K,1) = ciat.RealInterval;
P_He(K,1) = ciat.RealInterval;
P_Sc(K,1) = ciat.RealInterval;
P_An(K,1) = ciat.RealInterval;

% Calculate beampattern intervals for each angle
T_r = zeros(3,4);
T_c = zeros(3,4);
T_g = zeros(3,4);
T_x = zeros(3,4);
T_He = 0;
T_Hu = 0;
T_Sc = 0;
T_An = 0;


%%  Calculate beampattern

% Load data or calculate
if 1
    load Fig6_data.mat
else

% Calculate the beampattern bounds for the left half of the angle range
for k=1:K
    phi = ((0:M-1)-(M-1)/2)' * pi*sin(theta(k));
    v = exp(1j*phi);
    
    % Calculate element intervals without coupling
    EA_amp = ciat.PolarInterval(w.*(1+[-1 1].*ampErr/2) , ...
                                phi * [1 1]);
    EA_pha = ciat.PolarInterval(w.* [1,1] , ...
                                phi + [-1 1].*phaErr/2);
    EA_cal = ciat.PolarInterval(w.*(1+[-1 1].*ampErr/2) , ...
                                phi + [-1 1].*phaErr/2);
    EA_all = [EA_amp, EA_pha, EA_cal];

    % Calculate element interval with coupling
    E_cal = ciat.PolarInterval((ones(M,1) + ciat.RealInterval(-ampErr/2,ampErr/2)),...
                          ciat.RealInterval(phi + [-1 1]*phaErr/2));
    p_m = w;
    R_m = [0;w(1:M-1)].*[0;beta] + [w(2:M);0].*[beta;0];
    A_cpl = ciat.CircularInterval(p_m , R_m);
    EA_a = [ciat.PolyarcularInterval([EA_amp, EA_pha, EA_cal]) , ...
            ciat.PolyarcularInterval(E_cal , A_cpl)];

    % Calculate nominal beampattern
    B_nom(k,:) = w' * v;
    
    % Calculate beampattern intervals
    for s=1:4
        % if k >= kSeg(s) && k <= kSeg(s+1)
            % Rectangular
            tic;
            if s < 4
                EA_r = ciat.RectangularInterval(EA_all(:,s));
            else
                E_r = ciat.RectangularInterval(E_cal);
                A_r = ciat.RectangularInterval(A_cpl);
                EA_r = E_r .* A_r;
            end
            T_r(1,s) = T_r(1,s) + toc;
            tic;
            B_r(k,s) = sum(EA_r);
            T_r(2,s) = T_r(2,s) + toc;
    
            % Circular
            tic;
            if s < 4
                EA_c = ciat.CircularInterval(EA_all(:,s));
            else
                E_c = ciat.CircularInterval(E_cal);
                EA_c = E_c .* A_cpl;
            end
            T_c(1,s) = T_c(1,s) + toc;
            tic;
            B_c(k,s) = sum(EA_c);
            T_c(2,s) = T_c(2,s) + toc;
    
            % Polygonal
            tic
            if s < 4
                EA_g = ciat.PolygonalInterval(EA_all(:,s),...
                                                    'tolerance',conf.tol);
            else
                EA_g = ciat.PolygonalInterval(E_cal , A_cpl,...
                                                    'tolerance',conf.tol);
            end
            T_g(1,s) = T_g(1,s) + toc;
            tic;
            B_g(k,s) = sum(EA_g);
            T_g(2,s) = T_g(2,s) + toc;
    
            % Polyarcular
            tic
            if s < 4
                EA_x = ciat.PolyarxInterval(EA_all(:,s));
            else
                EA_x = ciat.PolyarxInterval(E_cal , A_cpl);
            end
            T_x(1,s) = T_x(1,s) + toc;
            tic;
            B_x(k,s) = sum(EA_x);
            T_x(2,s) = T_x(2,s) + toc;
        % end
    end

    % s = 1;
    % if k >= kSeg(s) && k <= kSeg(s+1)
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
    % end

    % s = 4;
    % if k >= kSeg(s) && k <= kSeg(s+1)
        % Schmid method
        tic
        wL2 = sqrt(sum(abs(w)*.2));
        maxDeltaB = sqrt(M) * wL2 * (w'*alpha + w(1:end-1)'*beta + ...
                                     w(2:end)'*beta);
        B_Sc = ciat.CircularInterval(B_nom(k),maxDeltaB);
        P_Sc(k) = abs(B_Sc).^2;
        T_Sc = T_Sc + toc;
        
        % Anselmi method
        tic;
        diagC0 = diag(alpha);
        diagC1 = diag(beta,1);
        diagCm1 = diag(beta,-1);
        Ca = ciat.CircularInterval(zeros(M) , diagC0);
        Cb = ciat.CircularInterval(zeros(M) , diagC1 + diagCm1);
        B_An = w' * (Ca + Cb + eye(M)) * v;
        P_An(k) = abs(B_An).^2;
        T_An = T_An + toc;
    % end
end

% Calculate power pattern intervals
P_nom = abs(B_nom).^2;
for s=1:4
    tic
    P_r(:,s) = abs(B_r(:,s)).^2; 
    T_r(3,s) = T_r(3,s) + toc;
    tic
    P_c(:,s) = abs(B_c(:,s)).^2; 
    T_c(3,s) = T_c(3,s) + toc;
    tic
    P_g(:,s) = abs(B_g(:,s)).^2; 
    T_g(3,s) = T_g(3,s) + toc;
    tic
    P_x(:,s) = abs(B_x(:,s)).^2; 
    T_x(3,s) = T_x(3,s) + toc;
end

% Replace zeros with eps so logarithm computes
P_r.Infimum(P_r.inf==0)=eps;
P_c.Infimum(P_c.inf==0)=eps;
P_g.Infimum(P_g.inf==0)=eps;
P_x.Infimum(P_x.inf==0)=eps;
P_Hu.Infimum(P_Hu.inf==0)=eps;
P_He.Infimum(P_He.inf==0)=eps;
P_Sc.Infimum(P_Sc.inf==0)=eps;
P_An.Infimum(P_An.inf==0)=eps;

% Calculate pattern tolerance
D_r = sum(P_r.width,1) / sum(P_nom);
D_c = sum(P_c.width,1) / sum(P_nom);
D_g = sum(P_g.width,1) / sum(P_nom);
D_x = sum(P_x.width,1) / sum(P_nom);
D_Hu = sum((P_x(:,1).Width + abs(P_Hu.sup-P_x(:,1).sup) + ...
            abs(P_Hu.inf-P_x(:,1).inf)),1) / sum(P_nom);
D_He = sum(P_He.width,1) / sum(P_nom);
D_Sc = sum(P_Sc.width,1) / sum(P_nom);
D_An = sum(P_An.width,1) / sum(P_nom);

end

%% Plot

% Set parameters
xRange = [-90 90];
yRange = [-60 5];
conf = getFig5conf(1);

% Set colors
cList = getColorList(conf.cID);

figure(5);clf;hold on
set(gca,'DefaultLineLineWidth',3)

% Set axis
thAxis = rad2deg(theta);

% Plot nominal response
p0 = plot(thAxis,db(P_nom),'k-','DisplayName','Nominal');

% Plot vertical splitting lines
plot(thAxis(kSeg(2:4))'.*[1;1],yRange'.*ones(1,3),'k:')

for s = 1:4
    
    kRange = kSeg(s):kSeg(s+1);
     % kRange = 1:K;

    % Polygonal
    pGs = plot(thAxis(kRange),db(P_g(kRange,s).sup),'-','color','b', ...
                                    'DisplayName','Polygon (sup)');
    pGi = plot(thAxis(kRange),db(P_g(kRange,s).inf),'-','color','b', ...
                                    'DisplayName','Polygon (inf)');

    % Polyarcular (convex)
    pAs = plot(thAxis(kRange),db(P_x(kRange,s).sup),'--', 'color','r',...
                                    'DisplayName','Polyarc (sup)');
    pAi = plot(thAxis(kRange),db(P_x(kRange,s).inf),'--', 'color','r',...
                                    'DisplayName','Polyarc (inf)');

    % Rectangular
    if s<4
        pRs = plot(thAxis(kRange),db(P_r(kRange,s).sup),':','color',cList(1,:), ...
                                        'DisplayName','Rectangle (sup)');
        pRi = plot(thAxis(kRange),db(P_r(kRange,s).inf),':','color',cList(1,:),...
                                        'DisplayName','Rectangle (inf)');
    end
        
    % % Circular
    % pCs = plot(thAxis(kRange),db(P_c(kRange,s).sup),'m-','color',cList(2,:),...
    %                                 'DisplayName','Circle (sup)');
    % pCi = plot(thAxis(kRange),db(P_c(kRange,s).inf),'m-','color',cList(2,:), ...
    %                                 'DisplayName','Circle (inf)');

    
    if s == 1
        % Hu
        pHuS = plot(thAxis(kRange),db(P_Hu(kRange).sup),'-','color',cList(3,:), ...
                                    'DisplayName','Hu (sup)');
        pHuI = plot(thAxis(kRange),db(P_Hu(kRange).inf),'-', 'color',cList(3,:),...
                                    'DisplayName','Hu (inf)');
    
        % He
        pHeS = plot(thAxis(kRange),db(P_He(kRange).sup),'-','color',cList(4,:),...
                                    'DisplayName','He (sup)');
        pHeI = plot(thAxis(kRange),db(P_He(kRange).inf),'-','color',cList(4,:), ...
                                    'DisplayName','He (inf)');
    end

    if s == 4
        % Schmid
        pScS = plot(thAxis(kRange),db(P_Sc(kRange).sup),'-', 'color',cList(3,:), ...
                                    'DisplayName','Schmid (sup)');
        pScI = plot(thAxis(kRange),db(P_Sc(kRange).inf),'-', 'color',cList(3,:), ...
                                    'DisplayName','Schmid (inf)');
    
        % Anselmi
        pAnS = plot(thAxis(kRange),db(P_An(kRange).sup),'-','color',cList(4,:), ...
                                    'DisplayName','Anselmi (sup)');
        pAnI = plot(thAxis(kRange),db(P_An(kRange).inf),'-','color',cList(4,:), ...
                                    'DisplayName','Anselmi (inf)');
    end
end

% Pattern tolerance index
    % Amplitude errors
annotText = {   'Amp. error',...
                ['\color{red}' ...
                    '\Delta^{(a)}=' num2str(100*D_x(1)/D_x(1),3) '%' ...
                    ', T^{(a)}=' num2str(sum(T_x(:,1)),'%0.1f') 's'],...
                ['\color{blue}' ...
                    '\Delta^{(g)}=' num2str(100*D_x(1)/D_g(1),3) '%'...
                    ', T^{(g)}=' num2str(sum(T_g(:,1)),'%0.1f') 's'],...
                ['\color[rgb]{' num2str(cList(1,:)) '}' ...
                    '\Delta^{(r)}=' num2str(100*D_x(1)/D_r(1),3) '%'...
                    ', T^{(r)}=' num2str(sum(T_r(:,1)),'%0.1f') 's'],...
                ['\color[rgb]{' num2str(cList(4,:)) '}' ...
                    '\Delta^{(M)}=' num2str(100*D_x(1)/D_He,3) '%'...
                    ', T^{(M)}=' num2str(T_He,'%0.1f') 's'],...
                ['\color[rgb]{' num2str(cList(3,:)) '}' ...
                    '\Delta^{(T)}=' num2str(100*D_x(1)/D_Hu,3) '%'...
                    ', T^{(T)}=' num2str(T_Hu,'%0.1f') 's']};
annotation('textbox',[0.165 0.42 .11 0.5],'String',annotText, ...
           'BackgroundColor','w',...
           'HorizontalAlignment','center','FitBoxToText','on');
    % Phase errors
annotText = {   'Phase error',...
                ['\color{red}' ...
                    '\Delta^{(a)}=' num2str(100*D_x(2)/D_x(2),3) '%' ...
                    ', T^{(a)}=' num2str(sum(T_x(:,2)),'%0.1f') 's'],...
                ['\color{blue}' ...
                    '\Delta^{(g)}=' num2str(100*D_x(2)/D_g(2),3) '%' ...
                    ', T^{(g)}=' num2str(sum(T_g(:,2)),'%0.1f') 's'],...
                ['\color[rgb]{' num2str(cList(1,:)) '}' ...
                    '\Delta^{(r)}=' num2str(100*D_x(2)/D_r(2),3) '%' ...
                    ', T^{(r)}=' num2str(sum(T_r(:,2)),'%0.1f') 's'],...
                ['\color[rgb]{' num2str(cList(2,:)) '}']};
annotation('textbox',[0.312,0.731,0.154,0.188],'String',annotText, ...
           'BackgroundColor','w','HorizontalAlignment','center');

    % Calibration errors
annotText = {   'Both errors',...
                ['\color{red}' ...
                    '\Delta^{(a)}=' num2str(100*D_x(3)/D_x(3),3) '%' ...
                    ', T^{(a)}=' num2str(sum(T_x(:,3)),'%0.1f') 's'],...
                ['\color{blue}' ...
                    '\Delta^{(g)}=' num2str(100*D_x(3)/D_g(3),3) '%' ...
                    ', T^{(g)}=' num2str(sum(T_g(:,3)),'%0.1f') 's'],...
                ['\color[rgb]{' num2str(cList(1,:)) '}' ...
                    '\Delta^{(r)}=' num2str(100*D_x(3)/D_r(3),3) '%' ...
                    ', T^{(r)}=' num2str(sum(T_r(:,3)),'%0.1f') 's'],...
                ['\color[rgb]{' num2str(cList(2,:)) '}']};
annotation('textbox',[0.568,0.726,0.154,0.194],'String',annotText, ...
           'BackgroundColor','w','HorizontalAlignment','center');


    % Mutual coupling 
annotText = {   'Errors & coupling',...
                ['\color{red}' ...
                    '\Delta^{(a)}=' num2str(100*D_x(4)/D_x(4),3) '%' ...
                    ', T^{(a)}=' num2str(sum(T_x(:,4)),'%0.1f') 's'],...
                ['\color{blue}' ...
                    '\Delta^{(g)}=' num2str(100*D_x(4)/D_g(4),3) '%' ...
                    ', T^{(g)}=' num2str(sum(T_g(:,4)),'%0.1f') 's'],...
                ['\color[rgb]{' num2str(cList(4,:)) '}' ...
                    '\Delta^{(A)}=' num2str(100*D_x(4)/D_An,3) '%' ...
                    ', T^{(A)}=' num2str(T_An,'%0.1f') 's'],...
                ['\color[rgb]{' num2str(cList(3,:)) '}' ...
                    '\Delta^{(S)}=' num2str(100*D_x(4)/D_Sc,3) '%' ...
                    ', T^{(S)}=' num2str(T_Sc,'%0.1f') 's'],...
                ['\color[rgb]{' num2str(cList(1,:)) '}']};
annotation('textbox',[0.739,0.681,0.162,0.239],'String',annotText, ...
           'BackgroundColor','w','HorizontalAlignment','center');


% Format figure
xlim(xRange)
ylim(yRange)
xlabel('Beam angle [deg]')
ylabel('Power')
%legend([p0 p1 p2 p3 p4 p5 p6 p7 p8],'Position',[0.517 0.33 0 0])
fontsize(20,'point')