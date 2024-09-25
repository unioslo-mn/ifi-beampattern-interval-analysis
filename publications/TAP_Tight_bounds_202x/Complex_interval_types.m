clear
% close all

%% Generate intervals

% Set intervals
E = ciat.PolarInterval(2.0 , 2.5 , 0.05*pi, 0.10*pi);
E(2) = ciat.PolarInterval(2.0 , 2.5 , -0.1*pi, -0.2*pi);
A = ciat.CircularInterval(1.8,0.15);
A(2) = ciat.CircularInterval(2.2,0.15);

% Wrap intervals in various types
Er = ciat.RectangularInterval(E); 
Ec = ciat.CircularInterval(E);    
Eg = ciat.PolygonalInterval(E,'Tolerance',0.1); 
Ex = ciat.PolyarxInterval(E);
Ar = ciat.RectangularInterval(A); 
Ac = A;    

% Multiply E and A
EAr = Er .* Ar; 
EAc = Ec .* Ac; 
EAg = ciat.PolygonalInterval(E,A,'Tolerance',0.1); 
EAx = ciat.PolyarxInterval(E,A); 


% Sum the EA intervals
Br = sum(EAr);
Bc = sum(EAc);
Bg = sum(EAg);
Bx = sum(EAx);


%% Sample intervals

% Set parameters
Nsmp = 100;

% Sample E and A intervals
Es = E.sample(Nsmp);
As = A.sample(Nsmp);

% Multiply E and A intervals
EAs{1} = Es{1} * As{1}.';
EAs{2} = Es{2} * As{2}.';

% Sum EA intervals
Bs = EAs{1}(:) + EAs{2}(:).';

% Get convex hulls
idx = convhull(real(EAs{1}),imag(EAs{1}));
EAs{1} = EAs{1}(idx);
idx = convhull(real(EAs{2}),imag(EAs{2}));
EAs{2} = EAs{2}(idx);
idx = convhull(real(Bs),imag(Bs));
Bs = Bs(idx);


%% Plot
figure(1);clf;hold on;axis equal
set(0,'DefaultLineLineWidth',2)

% Plot wraps of E intervals
Er.plot('b','DisplayName','Rectangular');
Ec.plot('c','DisplayName','Circular');
Eg.plot('r','DisplayName','Polygonal (convex)');
Ex.plot('g','DisplayName','Polyarcular (convex)');
Ar.plot('b');
Ac.plot('c');

% Plot true intervals
E.plot('k--','DisplayName','Exact interval');
A.plot('k--');

% Add legend without update
legend(legendUnq(gcf),'Location','southwest','AutoUpdate',false)

% Plot origin
plot(0,0,'k+')

% Plot EA intervals
EAr.plot('b');
EAc.plot('c');
EAg.plot('r');
EAx.plot('g');
plot(real(EAs{1}),imag(EAs{1}),'k--')
plot(real(EAs{2}),imag(EAs{2}),'k--')

% Plot B interval
Br.plot('b');
Bc.plot('c');
Bg.plot('r');
Bx.plot('g');
plot(real(Bs),imag(Bs),'k--')

% Plot samples
% scatter(real(EAs{1}),imag(EAs{1}),1,'k.')
% scatter(real(EAs{2}),imag(EAs{2}),1,'k.')
% scatter(real(Bs),imag(Bs),1,'k.')

% Add texts
text(real(Ec(1).Center),imag(Ec(1).Center),'E_1','HorizontalAlignment', 'center')
text(real(Ec(2).Center),imag(Ec(2).Center),'E_2','HorizontalAlignment', 'center')
text(real(EAc(1).Center),imag(EAc(1).Center),'E_1 A','HorizontalAlignment', 'center')
text(real(EAc(2).Center),imag(EAc(2).Center),'E_2 A','HorizontalAlignment', 'center')
text(real(Ac(1).Center),imag(Ac(1).Center),'A_1','HorizontalAlignment', 'center')
text(real(Ac(2).Center),imag(Ac(2).Center),'A_2','HorizontalAlignment', 'center')
text(real(Bc.Center),imag(Bc.Center),'B=E_1 A + E_2 A','HorizontalAlignment', 'center')
fontsize(20,'points')

% Plot axes
xL = xlim;
yL = ylim;
plot(xL,[0,0],'color',0.7*ones(1,3),'linewidth',0.1)
plot([0,0],yL,'color',0.7*ones(1,3),'linewidth',0.1)
