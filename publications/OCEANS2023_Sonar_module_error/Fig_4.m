clear
close all
%%

% Set parameters
M = 16;
N = [1 2 4 8 16];
w_apod = chebwin(M,40);
w_apod = w_apod / sum(w_apod);
errPsi = deg2rad(0.5);
AngRes = 0.35/2;
PolTol = 1e-6;

% Generate array and beampattern objects
% Initialize array
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

% Calculate nominal beampattern with tapering
bp_nom = biat.BeamPattern(array,'nominal','BeamResolutionDeg',AngRes);
P_nom = bp_nom.calculateBeamPattern;
bp_pol = biat.BeamPattern(array,'polygonal','BeamResolutionDeg',AngRes,...
                                            'PolygonTolerance',PolTol);
P_pol = bp_pol.calculateBeamPattern;

%% Calculate bounds for various block sizes

cntCase = length(N);
subArray = cell(1,cntCase);
block = cell(1,cntCase);
bp_blk = cell(1,cntCase);
P_blk = cell(1,cntCase);

for idxCase = 1:cntCase
    Nb = N(idxCase);
    Mb = M/Nb;
    subArray{idxCase}(1:Nb) = biat.SensorArray();

    for nb = 1:Nb
        % Create subarrays
        subArray{idxCase}(nb) = copy(array);
        subArray{idxCase}(nb).ElCount = Mb;
        subArray{idxCase}(nb).TaperParam = w_apod((1:Mb)+(nb-1)*Mb);
        subArray{idxCase}(nb).OrientError = 0;
    end

    % Create block
    block{idxCase} = copy(array);
    block{idxCase}.ElCount = Nb;
    block{idxCase}.ElPitchRatio= 0.5*Mb;
    block{idxCase}.TaperParam = ones(Nb,1);
    block{idxCase}.ElDiameterRatio = 0;
    block{idxCase}.OrientError = errPsi;

    % Calculate beampattern
    bp_blk{idxCase} = biat.BeamPattern(subArray{idxCase}(1:Nb),'polygonal',...
                                        'Block',block{idxCase},...
                                        'BeamResolutionDeg',AngRes,...
                                        'PolygonTolerance',PolTol);
    P_blk{idxCase} = bp_blk{idxCase}.calculateBeamPattern;
end

%% Plot arrays


figure(1);clf;

subplot(2,3,1);
array.plot;
title(sprintf('Array of %0.0f elements',M))

for idxCase = 1:cntCase
    subplot(2,3,idxCase+1);
    hold on
    subArray{idxCase}.plot('Block',block{idxCase})
    title(sprintf("Block of %0.0fx%0.0f elements",N(idxCase),M/N(idxCase)))
end

%% Plot complex intervals

incAngle = deg2rad(-15);

figure(2);clf;
sgtitle(sprintf("Complex intervals at %0.1f˚ angle",rad2deg(incAngle)))

subplot(2,3,1);
bp_pol.BeamAngle = incAngle;
bp_pol.BeamIndex = 0;
bp_pol.ElementIntervals.plot('bo-','LineWidth',2);
bp_pol.SubArrayInterval.plot('co-','LineWidth',2);
bp_pol.ArrayInterval.plot('go-','LineWidth',2);
bp_pol.TotalInterval.plot('ko-','LineWidth',2);
grid on
axis equal
title(sprintf('Array of %0.0f elements',M))

for idxCase = 1:cntCase
    subplot(2,3,idxCase+1);hold on
    bp_blk{idxCase}.BeamAngle = incAngle;
    bp_blk{idxCase}.BeamIndex = 0;
    bp_blk{idxCase}.ElementIntervals.plot('bo-','LineWidth',2);
    bp_blk{idxCase}.SubArrayInterval.plot('co-','LineWidth',2);
    if ~isempty(bp_blk{idxCase}.Block)
        bp_blk{idxCase}.BlockInterval.plot('rs-','LineWidth',2);
    end
    bp_blk{idxCase}.ArrayInterval.plot('go-','LineWidth',2);
    bp_blk{idxCase}.TotalInterval.plot('ko-','LineWidth',2);
    grid on
    axis equal
    title(sprintf("Block of %0.0fx%0.0f elements",N(idxCase),M/N(idxCase)))
end

%% Plot beampatterns
degs = rad2deg(bp_nom.BeamAngles);

fig = figure(3);clf
fig.Position = [200 200 700 420]; 
set(gca,'DefaultLineLineWidth',2)
hold on; xlim([-90,20]); ylim([-50,2]); grid on; set(gcf, 'color', 'white');
xlabel('\theta (deg)')
ylabel('Power (dB)')
set(gca,'FontSize',9)
yticks(-60:10: 0)
xticks(-90:15:90)

lineWidth = [1 5 4 3 1];
lineColor = repmat([0.0 0.0 0.5 0.8 0.8]',1,3);
plural = ["","s"];

for idxCase = [1 2 3 4]
    plot(degs, db(P_blk{idxCase}(:,2),'power'),'-',...
                            'Color',lineColor(idxCase,:),...
                            'linewidth', lineWidth(idxCase),...
                'DisplayName',sprintf("Upper bound (%0.0f block%s of %0.0f)",...
                         N(idxCase),plural((N(idxCase)>1)+1),M/N(idxCase)));
end

plot(degs, db(P_pol(:,2))/2, 'r:','linewidth', 5,...
                  'DisplayName', sprintf("Upper bound (%0.0f elements)",M));
plot(degs, db(P_nom(:,1))/2, 'color',[0.1,0.75,0.1],...
                'DisplayName', 'Nominal beampattern');
legend('location','northwest')

set(gca, 'fontsize', 14)
legend boxoff
%%
% exportgraphics(gca,'fig4.pdf')%,'Resolution',300)

%%
figure(4);clf;

subplot(2,1,1);
hold on; xlim([-90,20]); ylim([-50,2]); grid on; set(gcf, 'color', 'white');
xlabel('\theta (deg)')
ylabel('Power (dB)')
set(gca,'FontSize',9)
yticks(-60:10: 0)
xticks(-90:15:90)
for idxCase = [1 2 3 4]
    plot(degs, db(P_blk{idxCase}(:,2),'power'),'-',...
                            'Color',lineColor(idxCase,:),...
                            'linewidth', lineWidth(idxCase),...
                'DisplayName',sprintf("Upper bound (%0.0f block%s of %0.0f)",...
                         N(idxCase),plural((N(idxCase)>1)+1),M/N(idxCase)));
end

plot(degs, db(P_pol(:,2))/2, 'r:','linewidth', 5,...
                  'DisplayName', sprintf("Upper bound (%0.0f elements)",M));
plot(degs, db(P_nom(:,1))/2, 'color',[0.1,0.75,0.1],...
                'DisplayName', 'Nominal beampattern');

subplot(2,1,2);
hold on; xlim([-90,20]); grid on; set(gcf, 'color', 'white');
for idxCase = [1 2 3 4]
    plot(degs, sqrt(P_blk{idxCase}(:,2)) - sqrt(P_nom(:,1)),'-',...
                            'Color',lineColor(idxCase,:),...
                            'linewidth', lineWidth(idxCase),...
                'DisplayName',sprintf("Upper bound (%0.0f block%s of %0.0f)",...
                         N(idxCase),plural((N(idxCase)>1)+1),M/N(idxCase)));
end

