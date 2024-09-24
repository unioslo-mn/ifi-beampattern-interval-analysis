clear
% close all

%% Generate three polar intervals
% Set parameters
r = [3.0 , 2.5 , 0.6 ; ...
     4.0 , 3.0 , 0.8];
a = [1.2 , 0.5 , 0.8 ; ...
     1.3 , 0.7 , 1.0]*pi;

% Set intervals
pI = ciat.PolarInterval(r(1,:), r(2,:), a(1,:), a(2,:));

%% Wrap intervals in various types
prI = ciat.RectangularInterval(pI); 
pcI = ciat.CircularInterval(pI);    
pgI = ciat.PolygonalInterval(pI,'Tolerance',0.1); 
paI = ciat.PolyarcularInterval(pI); 
pxI = ciat.PolyarxInterval(pI);

%% Sum the first two intervals
prIs = prI(1) + prI(2); 
pcIs = pcI(1) + pcI(2); 
pgIs = pgI(1) + pgI(2); 
paIs = paI(1) + paI(2); 
pxIs = pxI(1) + pxI(2);

%% Multiply the sum with the third interval
prIp = prIs * prI(3); 
pcIp = pcIs * pcI(3); 
pgIp = pgIs * pgI(3); 
pxIp = pxIs * pI(3);

%% Sample true sum and multiply with polar samples
paIs_smp = paIs.sample(100);
pI3_smp = pI(3).sample(100);
paIp_smp = paIs_smp{:} * pI3_smp{:}.';
paIp_ind = boundary(real(paIp_smp(:)),imag(paIp_smp(:)),0.5);
paIp = paIp_smp(paIp_ind);

%% Plot
figure(1);clf;hold on;axis equal
set(0,'DefaultLineLineWidth',2)

% Plot intervals
pI.plot('k--','DisplayName','Interval');

% Plot wraps
prI.plot('b','DisplayName','Rectangular');
pcI.plot('c','DisplayName','Circular');
pgI.plot('r','DisplayName','Polygonal');
pxI.plot('g','DisplayName','Polyarcular');

% Add legend without update
legend(legendUnq(gcf),'AutoUpdate',0)

% Plot origin
plot(0,0,'k+')

% Plot sums
prIs.plot('b');
pcIs.plot('c');
pgIs.plot('r');
pxIs.plot('g');
paIs.plot('k--');

% Plot products
prIp.plot('b');
pcIp.plot('c');
pgIp.plot('r');
pxIp.plot('g');
plot(real(paIp),imag(paIp),'k--')

% Add texts
text(real(pcI(1).Center),imag(pcI(1).Center),'A')
text(real(pcI(2).Center),imag(pcI(2).Center),'B')
text(real(pcI(3).Center),imag(pcI(3).Center),'C')
text(real(pcIs.Center),imag(pcIs.Center),'A+B')
text(real(pcIp.Center),imag(pcIp.Center),'(A+B)xC')
fontsize(20,'points')

% Plot axes
xL = xlim;
yL = ylim;
plot(xL,[0,0],'color',0.7*ones(1,3),'linewidth',0.1)
plot([0,0],yL,'color',0.7*ones(1,3),'linewidth',0.1)
