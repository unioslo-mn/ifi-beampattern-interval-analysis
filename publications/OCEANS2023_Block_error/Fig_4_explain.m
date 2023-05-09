clear

%% Calculate beampattern intervap using our packages

% Set parameters
M = 16;
N = 2;
w_apod = chebwin(M,40);
w_apod = w_apod / sum(w_apod);
errPsi = deg2rad(1);
AngRes = 0.35/2;
PolTol = 1e-6;
PolIncl = true;

% Animation parameters
doAnimation = 0;
tPause = 0.2;
animFrmCnt = 10;

% Calculated parameters
Mb = M/N;
la = 1/4;
d = la/2; 
D = d * Mb;
r = repmat( (-(Mb-1)/2:(Mb-1)/2)' * d, 1,N);
rb = (-(N-1)/2:(N-1)/2)*D;
c = max(abs(w_apod));

% Define array and calculate beampattern
array = biat.SensorArray(       'ElCount',          M,...
                                'ElPitchRatio',     0.5,...
                                'ElDiameterRatio',  0,...
                                'Curvature',        0,...
                                'TaperType',        'custom',...
                                'TaperParam',       w_apod,...
                                'GainError',        0,...
                                'PhaseError',       0,...
                                'PosXError',        0,...
                                'PosYError',        0,...
                                'OrientError',      errPsi,...
                                'CouplingCoeff',    0,...
                                'SteeringAngle',    0);

% Calculate nominal beampattern                            
bp_nom = biat.BeamPattern(array,'nominal','BeamResolutionDeg',AngRes);
P_nom = bp_nom.calculateBeamPattern;


% Calculate beampattern interval
subArray(1:N) = biat.SensorArray();

for nb = 1:N
    % Create subarrays
    subArray(nb) = copy(array);
    subArray(nb).ElCount = Mb;
    subArray(nb).TaperParam = w_apod((1:Mb)+(nb-1)*Mb);
    subArray(nb).OrientError = 0;
end

% Create block
block = copy(array);
block.ElCount = N;
block.ElPitchRatio= 0.5*Mb;
block.TaperParam = ones(N,1);
block.ElDiameterRatio = 0;
block.OrientError = errPsi;

% Calculate beampattern
bp_int = biat.BeamPattern(subArray(1:N),'polygonal',...
                                        'Block',block,...
                                        'BeamResolutionDeg',0.35/2,...
                                        'PolygonTolerance',1e-6);
P_int = bp_int.calculateBeamPattern;

%% Find worst-case (maximum) response for each incidence angle and plot it

thCount = bp_nom.BeamCount;
thetas = bp_nom.BeamAngles;
wcAngles = zeros(thCount,2);
wcPowers = zeros(thCount,1);
for idx = 1:thCount
    th = thetas(idx);
    kx = 2*pi/la * sin(th);
    ky = 2*pi/la * cos(th);
    [wcAngles(idx,:), wcPowers(idx)] = findMaxPower(r,rb,kx,ky,errPsi,w_apod,0);    
end

%% Initialize figures %% Animate response with orientation error

% Parameters
th = deg2rad(61);

% Calculate nominal response
kx = 2*pi/la * sin(th);
ky = 2*pi/la * cos(th);
rx = r + rb;
ry = rx * 0;
s_nom = w_apod .* exp(1j*( kx*rx(:) ));

% Initialize figure
figure(1);clf;hold on
subplot(2,2,1);cla;hold on
subplot(2,2,2);cla; hold on
subplot(2,2,3:4);cla;hold on

% Plot unit circle and nominal complex response
subplot(2,2,1)
fimplicit(@(x,y) x.^2+y.^2-1,'k:')
p = plot(rx,ry,'k-');
scatter(real(s_nom),imag(s_nom),100,'k.')
scatter(real(s_nom/c),imag(s_nom/c),100,'ko')
quiver(0,0,real(sum(s_nom)),imag(sum(s_nom)),'k','AutoScale','off')
axis equal

% Plot beampattern and selected incidence angle
subplot(2,2,3:4)
plot(rad2deg(bp_nom.BeamAngles),db(P_nom(:,1),'power'),'k-')
plot(rad2deg(bp_nom.BeamAngles),db(P_int,'power'),'r-')
plot(rad2deg(th)*ones(1,2),[-50 0],'k--')
scatter(rad2deg(th),db(abs(sum(s_nom)).^2,'power'),50,'ko')
plot(rad2deg(thetas),db(wcPowers,'power'),'b-')
ylim([-50 0])

% Initialize plot handlers
p = plot(0,0,'+');

%

% Find maximum power
subplot(2,2,2)
[tiltAngles,~] = findMaxPower(r,rb,kx,ky,errPsi,w_apod,1);
view(130,30)
xlabel('dPsi-1')
ylabel('dPsi-2')
zlabel('Power')

% Define orientation deviation
if doAnimation
    dPsi = linspace(-errPsi,errPsi,animFrmCnt)' .* ones(1,2);
else
    dPsi = tiltAngles;
end

% Plot response
for idx = 1:size(dPsi,1)

    % Clear plots
    delete(p)
    
    % Calculate complex response
    rx = [ r(:,1) .* cos(dPsi(idx,1)) , ...
           r(:,2) .* cos(dPsi(idx,2)) ] + rb ;
    ry = [ r(:,1) .* sin(dPsi(idx,1)) , ...
           r(:,2) .* sin(dPsi(idx,2))];
    s_int = w_apod .* exp(1j*( kx*rx(:) + ky*ry(:) ));

    % Plot array and complex response
    subplot(2,2,1);
    p = plot(rx,ry,'r+');
    p = [p; scatter(real(s_int),imag(s_int),'r.')];
    p = [p; scatter(real(s_int/c),imag(s_int/c),100,'ro','filled')];
    p = [p; quiver(0,0,real(sum(s_int)),imag(sum(s_int)),'r',...
                                                    'AutoScale','off')];
    axis equal
    
    % Plot power value
    subplot(2,2,2)
    p = [p; scatter3(dPsi(idx,1),dPsi(idx,2),abs(sum(s_int)).^2,100,...
                                                        'ro','filled')];
    
    % Plot power response
    subplot(2,2,3:4);
    p = [p;scatter(rad2deg(th),db(abs(sum(s_int,1).^2),'power'),...
                                    50,'ro','filled')];
    
    drawnow
    pause(tPause)    
end

%%
% clf;hold on;
% ylim([-50,0])
% rx = r + rb + la*1.5;
% 
% s_nom = zeros(length(thetas),1);
% for idx = 1:length(thetas)
%     th = thetas(idx);
%     kx = 2*pi/la * sin(th);
%     s_nom(idx) = sum(w_apod .* exp(1j*( kx*rx(:) )));
% end
% plot(rad2deg(thetas),db(abs((s_nom)).^2,'power'),'k.-')
% plot(rad2deg(thetas),db(abs(real(s_nom)).^2,'power'),'b.-')
% plot(rad2deg(thetas),db(abs(imag(s_nom)).^2,'power'),'c.-')

%% Function to find maximum power

function [aOptGlob,powerOptGlob] = findMaxPower(r,rb,kx,ky,errPsi,w_apod,plotSurf)

    frx1 = @(a1) r(:,1) .* cos(a1) + rb(1);
    fry1 = @(a1) r(:,1) .* sin(a1);
    frx2 = @(a2) r(:,2) .* cos(a2) + rb(2);
    fry2 = @(a2) r(:,2) .* sin(a2);

    fph1 = @(a1) kx.*frx1(a1) + ky.*fry1(a1);
    fph2 = @(a2) kx.*frx2(a2) + ky.*fry2(a2);

    fp = @(a1,a2) ( sum(w_apod.*sin([fph1(a1);fph2(a2)])).^2 + ...
                     sum(w_apod.*cos([fph1(a1);fph2(a2)])).^2 );
    fOpt = @(a) 1-fp(a(1),a(2));

    x0 = errPsi/2 * ones(1,2) .* [1 1 ; 1 -1 ; -1 1 ; -1 -1 ];
    powerOptGlob = 0;
    for idx = 1:4
        aOpt = fminsearch(fOpt,x0(idx,:));
        if ~(	aOpt(1) > -errPsi && aOpt(1) < errPsi && ...
                aOpt(2) > -errPsi && aOpt(2) < errPsi)
            aOpt = fmincon(fOpt, aOpt , [], [], [], [], [-errPsi,-errPsi] , ...
                                                        [errPsi,errPsi]);
        end
        powerOpt = fp(aOpt(1),aOpt(2));
        if powerOpt > powerOptGlob
            aOptGlob = aOpt;
            powerOptGlob = powerOpt;
        end
    end
    
    if plotSurf
        fsurf(fp,[-errPsi,errPsi,-errPsi,errPsi])
        hold on
        scatter3(aOpt(1),aOpt(2),powerOpt,200,'ko','filled')
    end
    
end