clear
% close all

%% Generate a random annular sector

% Set parameters
N = 100;

% Generate random radius and angles
rad = sort(rand(2,1));
ang = (2*rand-1)*pi;
ang = [ang ; ang + rand*pi];

% Sample interval
smp = [rad(1) * exp(1j*linspace(ang(1),ang(2),N)') ; ...
       rad(2) * exp(1j*flip(linspace(ang(1),ang(2),N))')];
smp = [smp;smp(1)];

figure(1);clf;hold on
plot(0,0,'k+')
plot(real(smp),imag(smp),'g.-')
axis equal 

%% Wrap arc with polygon

% Set parameters
L = 3;

% Calculate parameters
angStep = diff(ang)/2/L;
radExt = rad(2) / cos(angStep);

% Set vertex locations
vert = zeros(L+4,1);
for l = 0:L
    vert(l+1) = radExt * exp(1j*(ang(1)+l*2*angStep));
end
vert(L+2) = rad(1)*exp(1j*ang(2)); 
vert(L+3) = rad(1)*exp(1j*ang(1));
vert(L+4) = vert(1);

% Plot polygon
plot(real(vert),imag(vert),'ko-')