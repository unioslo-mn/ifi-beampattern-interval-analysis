clear
% close all

%% Zhang 2017 (amplitude and phase errors)
N = 10;
ampErr = 0.1;
phaErr = 0.1;
w = taylorwin(N,3,-20);
% w = ones(N,1);
A = w * (1 + ciat.RealInterval(-ampErr/2,ampErr/2));
dTheta = 0.9;
theta = linspace(-pi,pi,N)' * dTheta;

% Define and cast intervals
AFp = ciat.PolarInterval(A,ciat.RealInterval(theta + [-1 1]*phaErr/2));
AFr = ciat.RectangularInterval(AFp);
AFg = ciat.PolygonalInterval(AFp);

% Sum intervals
AFr_sum = sum(AFr);
AFg_sum = sum(AFg);

% Plot
figure(1);clf;hold on;axis equal
set(gca,'DefaultLineLineWidth',2)
plot(0,0,'k+')
AFp.plot('g')
% AFg.plot('g')
AFr.plot('b')
AFr_sum.plot('b')
AFg_sum.plot('g')
xlabel('Real')
ylabel('Imag')
text(AFr_sum.real.mid,AFr_sum.imag.mid,'AF_\Sigma','fontsize',20)
xL = xlim(); yL = ylim();
fimplicit(@(x,y) x.^2+y.^2-(AFr_sum.real.sup^2+AFr_sum.imag.sup^2),'b--','linewidth',2)
fimplicit(@(x,y) x.^2+y.^2-(sup(abs(AFg_sum)))^2,'g--','linewidth',2)
xlim(xL); ylim(yL)
fontsize(20,'point')


%% Poli 2015 (phase errors)

% Define and cast intervals
AFp = ciat.PolarInterval(w , ciat.RealInterval(theta + [-1 1]*phaErr/2));
AFr = ciat.RectangularInterval(AFp);
AFg = ciat.PolygonalInterval(AFp);

% Sum intervals
AFr_sum = sum(AFr);
AFg_sum = sum(AFg);

%% Plot
% figure(1);clf;hold on;axis equal;
set(gca,'DefaultLineLineWidth',2)
plot(0,0,'k+')
AFp.plot('g')
% AFg.plot('g')
AFr.plot('b')
AFr_sum.plot('b')
AFg_sum.plot('g')
xlabel('Real')
ylabel('Imag')
text(AFr_sum.real.mid,AFr_sum.imag.mid,'AF_\Sigma','fontsize',20)
xL = xlim(); yL = ylim();
fimplicit(@(x,y) x.^2+y.^2-(AFr_sum.real.sup^2+AFr_sum.imag.sup^2),'b--','linewidth',2)
fimplicit(@(x,y) x.^2+y.^2-(sup(abs(AFg_sum)))^2,'g--','linewidth',2)
xlim(xL); ylim(yL)
fontsize(20,'point')


