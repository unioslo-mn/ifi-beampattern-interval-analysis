clear
% close all

%%

% Set parameters
M = 8;
SNR = 20;
alpha = [5 5 10 10 10 10 5 5]'*1e-2;
beta = [1 2 3 4 3 2 1]' * 1e-2;
dTheta = 0.85;
theta = linspace(-pi,pi,M)' * dTheta;
w = chebwin(M,SNR) ; w = w / sum(w);
v = exp(1j*theta);
B = w' * v;

% Generate polyarcular intervals
phi = 0.01;
Eint = ciat.PolarInterval(w.*(1+[-1 1].*alpha*0.75) , angle(v')' + [-1 1].*phi*pi);
Aint = ciat.CircularInterval(zeros(M,1) , [0;beta]+[beta;0]);
aI = ciat.PolyarxInterval( Eint , ones(M,1)+Aint );
aIsum = sum(aI);

% Generate circular intervals
diagC0 = diag(alpha);
diagC1 = diag(beta,1);
diagCm1 = diag(beta,-1);
I = eye(M);
% Ca = ciat.CircularInterval(zeros(M) , diagC0);
CaRad = radius(ciat.CircularInterval(Eint))./w
Ca = ciat.CircularInterval(zeros(M) , diag(CaRad));
Cb = ciat.CircularInterval(zeros(M) , diagC1 + diagCm1);
AFint = (w.*v)' * (Ca + Cb + I);
AFsum = sum(AFint);
Bint = w' * (Ca + Cb + I) * v;

% Without coupling
AFa = (w.*v)' * (Ca + I);
AFint2 = AFa.' .* (ones(M,1)+Aint);
AFsum2 = sum(AFint);

% L2 norm and upper bound on max delta B 
wL2 = sqrt(sum(abs(w)*.2));
maxDeltaB = sqrt(M) * wL2 * (w'*alpha + w(1:end-1)'*beta + w(2:end)'*beta);

%% Plot
figure(1);clf;hold on;axis equal
set(gca,'DefaultLineLineWidth',2)
plot(0,0,'k+')
plot(real(w.*v),imag(w.*v),'bx')
plot(real(B),imag(B),'x')
% AFint.plot('b-');
AFint2.plot('b-');
Bint.plot('b--'); %==AFsum.plot('b--') == AFsum2.plot('k--')
fimplicit(@(x,y) (x-real(B)).^2 + (y-imag(B)).^2 - maxDeltaB^2,'r--','linewidth',2)
Eint.plot('g-');
AFa.plot('b-');
aI.plot('g-');
aIsum.plot('g-');
xlabel('Real')
ylabel('Imag')
text(real(B),imag(B),'B')
fontsize(20,'point')
