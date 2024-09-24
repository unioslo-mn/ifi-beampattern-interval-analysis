clear

% Define operand intervals
pI = ciat.PolarInterval(1,2,3,4);
cI = ciat.CircularInterval(1,0.3);

% Calculate product intervals
gI = ciat.PolygonalInterval(pI,cI);
aI = ciat.PolyarcularInterval(pI,cI);
xI = ciat.PolyarxInterval(pI,cI);

% Calculate corner circles
cIp(4,1) = ciat.CircularInterval;
cIp(1) = cI * pI.Abs.inf*exp(1j*pI.Angle.inf);
cIp(2) = cI * pI.Abs.sup*exp(1j*pI.Angle.inf);
cIp(3) = cI * pI.Abs.inf*exp(1j*pI.Angle.sup);
cIp(4) = cI * pI.Abs.sup*exp(1j*pI.Angle.sup);

% Sample operands and calculate Minkowski product
pIsmp = pI.sample(30); pIsmp = pIsmp{:};
cIsmp = cI.Center + cI.Radius * exp(1j*linspace(-pi,pi,30));
gIsmp = pIsmp * cIsmp;

% Plot
figure(1);clf;hold on; axis equal
pI.plot('b-')
cI.plot('r-')
gI.plot('g-')
aI.plot('c-')
xI.plot('k--')
cIp.plot('r-')
plot(real(pIsmp),imag(pIsmp),'b.')
plot(real(cIsmp),imag(cIsmp),'r.')
plot(real(gIsmp),imag(gIsmp),'g.')
plot(0,0,'+')