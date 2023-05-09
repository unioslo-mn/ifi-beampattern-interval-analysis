clearvars;

% Code works by having one loop for each case. 
% Loops through PSLL to find worst case.
% This code takes some time to run.


% Set nominal PSLL
PSLL_nom = 10 : 2 : 50;
                                        
% Get parameters
L = length(PSLL_nom);

% Set error parameter arrays
elDiameterRatio = [0 , 0 , 0 , 0.95];
gainError = [ 5/100 , 5/100 , 1/100 , 0 ]; 
phaseError = deg2rad([ 5 , 1 , 1 , 0 ]);
orientError = deg2rad([ 0 , 0 , 0 , 2 ]);
couplingCoeff = [ 0.05 , 0.00 , 0.00 , 0.01 ];

%% Calculate results

% Preallocate memory
PSLL_bounds = zeros(L, 4);
PSLL_approx = zeros(L, 4);
r_tot_bound = zeros(1,4);

for iCase = 1:4
    
    % Generate array object
    array = biat.SensorArray(   'ElCount',30,...
                                'ElDiameterRatio',elDiameterRatio(iCase),...
                                'SteeringAngle',0,...
                                'TaperType','chebwin',...
                                'TaperParam',PSLL_nom(1),...
                                'GainError', gainError(iCase),...
                                'PhaseError', phaseError(iCase),...
                                'OrientError',orientError(iCase),...
                                'CouplingCoeff', couplingCoeff(iCase));
                            
	% Generate beampattern objects
    bp_nom = biat.BeamPattern(array,'nominal','BeamResolutionDeg',0.35);
    bp_pol = biat.BeamPattern(array,'polygonal','BeamResolutionDeg',0.35,...
                                            'PolygonTolerance',1e-4);
       
	% Get angles
    thetas = bp_nom.BeamAngles;
    
    for sl = 1:L

        % focus, window and steering
        array.TaperParam = PSLL_nom(sl);

        % Calculate beampattern
        P_nom = bp_nom.calculateBeamPattern;
        P_pol = bp_pol.calculateBeamPattern;

        % Extract peak side-lobe level
        [PSLL_bounds(sl,iCase),~] = PSLL( P_pol(:,2), P_nom(:,1), thetas); 
        [PSLL_approx(sl,iCase),~] = PSLL( P_nom(:,3), P_nom(:,1), thetas); 
    end
    circle_bound = sqrt( (array.GainInterval(1).Supremum - 1)^2 + ...
                          array.PhaseInterval(1).Supremum^2);
    r_tot_bound(iCase) = circle_bound + 2 * couplingCoeff(iCase);
end

%% Plot
set(gcf,'Position',[600 200 500 300])
set(gca, 'fontsize', 16)
set(groot,'defaultlinelinewidth', 3)
hold on; 
grid on; set(gcf, 'color', 'white'); %set(gca, 'color', [.85 .85 .85]);
xlabel('Nominal {\it|}PSLL{\it|} (dB)')
ylabel('PSLL (dB)')
yticks([-50:10:0])
xticks(-50:10:50);%[-90,-45,0,45,90]/2)
xlim([10,50])
ylim([-40,-5])

% Set labels
caseName = {'Case 1','Case 2','Case 3','Case 4'};
caseCol = [[1 .2 .2] ; [1 .2 .2]*2/4 ; [1 .2 .2]*1/4 ; [.3 .6 1]];
lineVis = {'off', 'off', 'off';'off', 'off', 'on'};
lineCol = {'w','k'};
lineWid = {2,1};

figure(1); hold on;
plot( abs(PSLL_nom), -PSLL_nom, 'linestyle', '-', ...
                                'color', [.1 .8 .1],...
                                'displayname', 'Nom.');
for iCase = 1:4
    plot( abs(PSLL_nom), db(PSLL_bounds(:,iCase),'power'),...
                                'linewidth', 4.5, ...
                                'linestyle', '-.', ...
                                'color', caseCol(iCase,:), ...
                                'displayname', caseName{iCase});
    if iCase < 4
        for iLine = 1:2
            plot( abs(PSLL_nom), db(PSLL_approx(:,iCase),'power'),...
                                'linewidth', lineWid{iLine}, ...
                                'linestyle', '-', ...
                                'color', lineCol{iLine}, ...
                                'handlevisibility', lineVis{iLine,iCase},...
                                'displayname', 'Approximations');
        end
        plot( abs(PSLL_nom)+30, db(r_tot_bound(iCase))*ones(1, L),...
                                'linewidth', 2, ...
                                'linestyle', ':', ...
                                'color', [.4 .4 .4], ...
                                'handlevisibility', lineVis{2,iCase},...
                                'displayname', 'Asymptotic limits');
    end
end

legend('location','southwest')
legend boxoff

%exportgraphics(gcf,'Fig11_PSLL.jpg','Resolution',300) 
