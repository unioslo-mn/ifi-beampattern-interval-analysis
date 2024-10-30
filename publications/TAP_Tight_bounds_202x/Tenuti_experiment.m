clear
% close all

%% Set parameters

% Array
M = 8;
N = M;

% Tapering
apod = taylorwin(M,3,-25);
apod = apod .* apod.' / sum(apod)^2;

% Errors
errAmp = 0.01;
errAng = deg2rad(3);

% Plot
dinRange = [-60 5];

%% Calculate beampattern

% Angular axes
K = 100;
L = K;
u = linspace(-1,1,K);
v = u;

% Nominal beampattern
nomBP = zeros(K,L);
for k = 1:K
    for l = 1:L
        delay = pi*((1:M)*u(k) + (1:N)'*v(l));
        W = apod .* exp(1j*delay);
        nomBP(k,l) = abs(sum(W,'all'));
    end
end

%% Calculate beampattern interval for a given intersection


recBP(K,1) = ciat.RealInterval;
gonBP(K,1) = ciat.RealInterval;
arxBP(K,1) = ciat.RealInterval;
gonTime = 0;
arxTime = 0;
for k = K/2%1:K
    for l = 1:L
        % Calculate phase delay
        delay = pi*((1:M)*u(k) + (1:N)'*v(l));

        % Generate element intervals in polar type
        pI = ciat.PolarInterval(apod*(1-errAmp) , apod*(1+errAmp),...
                                delay-errAng , delay + errAng);

        % Cast element intervals to rectangular and calculate sum
        rI = ciat.RectangularInterval(pI);
        rIsum = sum(rI,'all');
        recBP(l) = abs(rIsum);

        % Cast element intervals to polygon and calculate sum
        tic
        gI = ciat.PolygonalInterval(pI,'tolerance',1e-4);
        gIsum = sum(gI,'all');
        gonBP(l) = abs(gIsum);
        gonTime = gonTime + toc;

        % Cast element intervals to polygon and calculate sum
        tic
        xI = ciat.PolyarxInterval(pI);
        xIsum = sum(xI(:));
        arxBP(l) = abs(xIsum);
        arxTime = arxTime + toc;
    end
end

recBP.Infimum(recBP.Infimum==0) = eps;
gonBP.Infimum(gonBP.Infimum==0) = eps;
arxBP.Infimum(arxBP.Infimum==0) = eps;
%% Plot
figure(1);clf;
cla
imagesc(u,v,db(nomBP))
clim(dinRange)
colormap(turbo)
colorbar

figure(2);clf;hold on
plot(u,db(recBP.Infimum),'r-')
plot(u,db(recBP.Supremum),'r-')
plot(u,db(gonBP.Infimum),'b-')
plot(u,db(gonBP.Supremum),'b-')
plot(u,db(arxBP.Infimum),'g--')
plot(u,db(arxBP.Supremum),'g--')
ylim(dinRange)

sprintf('Polygonal time: %0.3fms\n Polyarx time: %0.3f',gonTime*1e3,arxTime*1e3)