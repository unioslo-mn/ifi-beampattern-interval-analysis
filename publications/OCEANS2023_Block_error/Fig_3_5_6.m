%%
close all; clear;

%% Set simulation parameters (include update function)
% Note: there is a certain order in which things must be set as of now


% Set parameters
M = 16;
N = [1 2 4 8];
w_apod = chebwin(M,40);
w_apod = w_apod / sum(w_apod);
errPosX = 5e-4;
errPosY = 5e-4;
AngRes = 0.35/2;
PolTol = 1e-6;
PolIncl = true;

% Initialize array
array = biat.SensorArray(       'ElCount',          M,...
                                'ElPitchRatio',     0.5,...
                                'ElDiameterRatio',  0,...
                                'Curvature',        0,...
                                'TaperType',        'custom',...
                                'TaperParam',       w_apod,...
                                'GainError',        0,...
                                'PhaseError',       0,...
                                'PosXError',        errPosX,...
                                'PosYError',        errPosY,...
                                'CouplingCoeff',    0,...
                                'SteeringAngle',    0);
% Calculate nominal beampattern
bp_nom = biat.BeamPattern(array,'nominal','BeamResolutionDeg',0.35/2);
P_nom = bp_nom.calculateBeamPattern;

%%
% Generate subarrays and block for each case and calculate beampattern
P_pol = cell(length(N),1);
for iCase = 1:length(N)
    Nb = N(iCase);
    Mb = M/Nb;
    subArray(1:Nb) = biat.SensorArray();
    
    for nb = 1:Nb
        % Create subarrays
        subArray(nb) = copy(array);
        subArray(nb).ElCount = Mb;
        subArray(nb).TaperParam = w_apod((1:Mb)+(nb-1)*Mb);
        subArray(nb).PosXError = errPosX * (Nb==1);
        subArray(nb).PosYError = errPosY * (Nb==1);
    end

    % Create block
    block = copy(array);
    block.ElCount = Nb;
    block.ElPitchRatio= 0.5*Mb;
    block.TaperParam = ones(Nb);
    block.PosXError = errPosX * (Nb>1);
    block.PosYError = errPosY * (Nb>1);

    % Calculate beampattern
    bp_pol = biat.BeamPattern(subArray(1:Nb),'polygonal',...
                                        'Block',block,...
                                        'BeamResolutionDeg',0.35/2,...
                                        'PolygonTolerance',1e-6);
    P_pol{iCase} = bp_pol.calculateBeamPattern;

end



%% Figure 3.

degs = rad2deg(bp_nom.BeamAngles);

% figure('Position',[200 200 500 300]); 
fig1 = figure(1);clf
fig1.Position = [200 200 700 420]; 
set(gca,'DefaultLineLineWidth',2)
hold on; xlim([-90,20]); ylim([-50,2]); grid on; set(gcf, 'color', 'white');
xlabel('\theta (deg)')
ylabel('Power (dB)')
set(gca,'FontSize',9)
yticks(-60:10: 0)
xticks(-90:15:90)

plot( degs, db(P_pol{2}(:,2))/2, 'linestyle','-','color',  0.0*[1 1 1],...
            'linewidth', 5,'DisplayName', 'Upper bound (2 blocks of 8)');
plot(degs, db(P_pol{3}(:,2))/2, 'linestyle','-','color',  0.5*[1 1 1],...
            'linewidth', 4,'DisplayName', 'Upper bound (4 blocks of 4)');
plot(degs, db(P_pol{4}(:,2))/2, 'linestyle','-','color',  0.8*[1 1 1],...
            'linewidth', 3,'DisplayName', 'Upper bound (8 blocks of 2)');
plot(degs, db(P_pol{1}(:,2))/2, 'r:',...
            'linewidth', 3, 'DisplayName', 'Upper bound (16 elements)');
plot(degs, db(P_nom(:,1))/2, 'color',[0.1,0.75,0.1],...
            'DisplayName', 'Nominal beampattern');

legend('location','northwest')

set(gca, 'fontsize', 14)
legend boxoff

%%
exportgraphics(gca,'fig3.pdf')%,'Resolution',300)

%% Figure 5
fig2 = figure(2);clf
fig2.Position = [200 200 700 420]; 
set(gca,'DefaultLineLineWidth',2)
hold on; grid on; set(gcf, 'color', 'white');
xlabel('Angle (deg)')
ylabel('Power (dB)')
set(gca,'FontSize',9)


A_el =  (P_pol{1}(:,2));
A_2 =  (P_pol{2}(:,2));
A_4 =  (P_pol{3}(:,2));
A_8 =  (P_pol{4}(:,2));

dpl = (degs >= 0);

plot( (degs(dpl)), db(A_2(dpl))/2, 'linestyle','-','color',  0.0*[1 1 1],...
            'linewidth', 5,'DisplayName', 'Upper bound (2 blocks of 8)');
plot( (degs(dpl)), db(A_4(dpl))/2, 'linestyle','-','color',  0.5*[1 1 1],...
            'linewidth', 4,'DisplayName', 'Upper bound (4 blocks of 4)');
plot( (degs(dpl)), db(A_8(dpl))/2, 'linestyle','-','color',  0.8*[1 1 1],...
            'linewidth', 3,'DisplayName', 'Upper bound (8 blocks of 2)');
plot( (degs(dpl)), db(A_el(dpl))/2, 'r:','DisplayName', ...
            'Upper bound (16 elements)', 'linewidth', 3);
plot( (degs(dpl)), db(P_nom(dpl))/2, 'color',[0.1,0.75,0.1],...
            'DisplayName', 'Nominal beampattern');

legend('location','northeast')

set(gca, 'fontsize', 14)
legend boxoff

lambda = array.WaveLength;
pitch = array.ElPitch;
kabs = 2*pi/lambda;

for Ni = [1,2,4,8]
    for ni = [1,3,5,7,9,11,13,15]
        ang = ni*pi/(Ni*pitch*kabs);
        if ang > 1; continue;
        else
            if ang == 1; GL = 1; % Grating lobe
            else; GL = 0; end;
            
            switch Ni
                case 1
                    colr = 1.0*[1 1 1];
                case 2
                    colr = 0.8*[1 1 1];
                case 4
                    colr = 0.5*[1 1 1];
                case 8
                    colr = 0*[1 1 1];
            end
            
            plot( asind(ang) ,db( ((1+GL)*[4*(errPosX*ang + ...
                                              errPosY*sqrt(1-(ang)^2)) / ...
                                              (ni*lambda)]).^2 )/2, 'o', ...
                                              'linewidth', 2, ...
                                              'markerfacecolor', colr, ...
                                              'markersize', 10, ...
                                              'markeredgecolor', 'r', ...
                                              'handlevisibility', 'off')
        end
    end
    
end



plot( [-1] ,[-1], 'o', 'linewidth', 1, 'markerfacecolor', 'r', ...
                      'markersize', 10, 'markeredgecolor', 'w', ...
                      'displayname', '(|B_0| + \beta_n peaks)^2')
%us = [0:0.01:1];
us = sind(degs(dpl));
plot( degs(dpl), db( (sqrt(P_nom(dpl)) +4*(errPosX*us + ...
                                          errPosY*sqrt(1-(us).^2)) / ...
                                          (lambda)).^2)/2,'b-', ...
                                            'linewidth', 1,...
                                            'displayname','(|B_0| + E_1)^2')

ylim([-50,1])
xlim([0,90])
xticks([-90:15:90])
legend('location', 'northeast');
xlabel('\theta (deg)')
ylabel('Power (dB)')
set(gca,'FontSize',12)
legend boxoff

%exportgraphics(gca,'fig5.pdf')%,'Resolution',300)

%% Figure 6

fig3 = figure(3);clf
fig3.Position = [200 200 700 420]; 

A_el =  sqrt(P_pol{1}(:,2)) - sqrt(P_nom);
A_2 =  sqrt(P_pol{2}(:,2)) - sqrt(P_nom);
A_4 =  sqrt(P_pol{3}(:,2)) - sqrt(P_nom);
A_8 =  sqrt(P_pol{4}(:,2)) - sqrt(P_nom);

dpl = (degs >= 0);
area((degs(dpl)), A_el(dpl), 'facecolor', 1*[1 0 0],...
                'DisplayName', '16 elem.', 'linewidth', 1); hold on;
area((degs(dpl)), A_8(dpl),  'facecolor', 0.8*[1 1 1],...
                'DisplayName', '8 bl. of 2', 'linewidth', 1);
area((degs(dpl)), A_4(dpl),  'facecolor', 0.5*[1 1 1],...
                'DisplayName', '4 bl. of 4', 'linewidth', 1);
area((degs(dpl)), A_2(dpl),  'facecolor', 0.0*[1 1 1],...
                'DisplayName', '2 bl. of 8', 'linewidth', 1);


plot( (degs(dpl)), A_2(dpl), 'linestyle','-','color',  0.0*[1 1 1],...
                'linewidth', 6,'DisplayName', 'Upper bound (2 blocks of 8)',...
                'handlevisibility', 'off')
plot( (degs(dpl)), A_4(dpl), 'linestyle','-','color',  0.5*[1 1 1],...
                'linewidth', 4,'DisplayName', 'Upper bound (8 blocks of 2)',...
                'handlevisibility', 'off')
plot( (degs(dpl)), A_8(dpl), 'linestyle','-','color',  0.8*[1 1 1],...
                'linewidth', 2,'DisplayName', 'Upper bound (4 blocks of 4)',...
                'handlevisibility', 'off')
plot( (degs(dpl)), A_el(dpl), 'r:','DisplayName', 'Upper bound (16 elements)', ...
                                'linewidth', 3,'handlevisibility', 'off')

for Ni = [1,2,4,8]
    for ni = [1,3,5,7,9,11,13,15]
        ang = ni*pi/(Ni*pitch*kabs);
        if ang > 1; continue;
        else
            if ang == 1; GL = 1; % Grating lobe
            else; GL = 0; end;
            
            switch Ni
                case 1
                    colr = 1.0*[1 1 1];
                case 2
                    colr = 0.8*[1 1 1];
                case 4
                    colr = 0.5*[1 1 1];
                case 8
                    colr = 0*[1 1 1];
            end
            
            plot( asind(ang) ,( ((1+GL)*[4*(errPosX*ang + ...
                                            errPosY*sqrt(1-(ang)^2)) / ...
                                            (ni*lambda)]) ), 'o', ...
                'linewidth', 2, 'markerfacecolor', colr, 'markersize', 10, ...
                'markeredgecolor', 'r', 'handlevisibility', 'off')
        end
    end
    
end



plot( [-1] ,[-1], 'o', 'linewidth', 1, 'markerfacecolor', 'r', ...
                    'markersize', 10, 'markeredgecolor', 'w', ...
                    'displayname', '\beta_n peaks')
us = [0:0.01:1];
plot( asind(us), 4*(errPosX*us + errPosY*sqrt(1-(us).^2))/(lambda),...
                        'b-', 'linewidth', 1,'displayname','E_1')

xticks([-90:15:90])
ylim([0,0.058])
xlim([0,90])
legend('location', 'northwest');
xlabel('\theta (deg)')
ylabel('$\sqrt{\overline{P}} - \sqrt{{P_0}}$','interpreter','latex')
set(gca,'FontSize',12)
legend boxoff
%% 

%exportgraphics(gca,'fig6.pdf')%,'Resolution',300)