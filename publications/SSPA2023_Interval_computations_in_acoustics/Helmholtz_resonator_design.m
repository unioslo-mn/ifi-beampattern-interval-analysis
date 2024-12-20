%% Tolerance analysis of a Helmholtz resonator design

% This is an example how acoustic design errors can add up in the 
% worst-case to a large performance error. Tolerance analysis enables
% to provide guaranteed bounds on the performance.

% Define resonance frequency function
f = @(c,A,L,V) c/(2*pi) * sqrt( A/(V*L) );

% Define nominal parameters
cNom = 340;     % Sound speed [m/s]
ANom = 0.2;     % Surface of the neck [m^2]
LNom = 0.15;    % Length of the neck [m]
VNom = 0.04;    % Volume of the cavity [m^3]

% Relative errors of the parameters
cErr = 0.01;
AErr = 0.01;
LErr = 0.01;
VErr = 0.01;

% Calculate parameter intervals
cInt = ciat.RealInterval( cNom*(1-cErr) , cNom*(1+cErr) );
AInt = ciat.RealInterval( ANom*(1-AErr) , ANom*(1+AErr) );
LInt = ciat.RealInterval( LNom*(1-LErr) , LNom*(1+LErr) );
VInt = ciat.RealInterval( VNom*(1-VErr) , VNom*(1+VErr) );

% Calculate nominal frequency and frequency interval 
fNom = f(cNom,ANom,LNom,VNom);
fInt = f(cInt,AInt,LInt,VInt);

% Calculate relative error of the resonance frequency
fErr = fInt.Width / fNom

