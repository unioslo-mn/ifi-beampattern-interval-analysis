clear
close all

%% Calculate

% Parameters
M = 16;
w_apod = chebwin(M,40);
w_apod = w_apod / sum(w_apod);
angResDeg = 0.35/2;
N = [2,4,8];
polTol = 1e-6;
elPosXErr = [0, 5e-4, 0, 0, 5e-4, 0]; 
elPosYErr = [0, 5e-4, 0, 0, 5e-4, 0]; 
blPosXErr = [5e-4, 5e-4, 0, 5e-4, 5e-4, 0]; 
blPosYErr = [5e-4, 5e-4, 0, 5e-4, 5e-4, 0]; 
blOrientErr = deg2rad([0, 0, 0.5, 0, 0, 0.5]);
elDiameter = ones(1,6)*0.90;
curvature = [0 0 0 0.5 0.5 0.5];

% Indexes
cntCase = length(N);
cntScen = length(curvature);

% Generate nominal array 
array = biat.SensorArray(       'ElCount',          M,...
                                'ElDiameterRatio',  0,...
                                'Curvature',        0,...
                                'TaperType',        'custom',...
                                'TaperParam',       w_apod,...
                                'GainError',        0,...
                                'PhaseError',       0,...
                                'PosXError',        0,...   
                                'PosYError',        0,...
                                'OrientError',      0,...
                                'CouplingCoeff',    0,...
                                'SteeringAngle',    0);

% Generate array with errors and block array
block = copy(array);
block.ElCount = 1;
block.TaperType = 'custom';
block.TaperParam = 1;        
                            
% Generate nominal beampattern
bp_nom = biat.BeamPattern(array,'nominal',...
                              'BeamResolutionDeg',angResDeg);
                          
% Generate polygonal beampattern
bp_err = biat.BeamPattern(array,'polygonal',...
                                  'BeamResolutionDeg',angResDeg,...
                                  'Block',block,...
                                  'BlockType','split',...
                                  'PolygonTolerance',polTol);

%% Plot arrays and calculate beampatterns

% Initialize figure
fig1 = figure(1);clf
fig1.Position = [200 200 750 450]; 

% Loop throught the scenarios
P_nom = zeros(bp_nom.BeamCount,cntScen);
P_err = zeros(bp_nom.BeamCount,cntCase,cntScen);
for idxScen = 1:cntScen
    
    % Adjust curvature and element diameter in all arrays
    array.ElDiameterRatio = elDiameter(idxScen);
    array.Curvature = curvature(idxScen);
    array.PosXError = elPosXErr(idxScen);
    array.PosYError = elPosYErr(idxScen);
    block.Curvature = curvature(idxScen);
    block.PosXError = blPosXErr(idxScen);
    block.PosYError = blPosYErr(idxScen);
    block.OrientError = blOrientErr(idxScen);
                            
    % Calculate nominal beampattern
    P_temp = bp_nom.calculateBeamPattern;
    P_nom(:,idxScen) = P_temp(:,1);

    % Loop through the module number variants
    for idxCase = 1:cntCase
        
        % Adjust block array
        block.ElCount = N(idxCase);
        block.ElPitchRatio= 0.5*M/N(idxCase); 
        block.ElDiameterRatio = 0;
        block.TaperType = 'custom';
        block.TaperParam = ones(N(idxCase),1);

        % Plot array
        subplot(2,3,idxScen);
        set(gca,'DefaultLineLineWidth',2)
        xlabel('Position-X')
        ylabel('Position-Y')
        set(gca,'FontSize',9)
        set(gcf, 'color', 'white');
        grid on; 
        if idxCase == 1
            bp_nom.plotArrays('PlotArrow',0)
            bp_err.plotArrays('PlotArrow',0)
        end

        % Calculate beampattern with various module numbers
    	P_temp = bp_err.calculateBeamPattern;
    	P_err(:,idxCase,idxScen) = P_temp(:,2);
    end   
end

%% Plot beampatterns
degs = rad2deg(bp_nom.BeamAngles);
lineWidth = [5 4 3 1];
lineColor = repmat([0.0 0.5 0.8 0.8]',1,3);
plural = ["","s"];
titles = ["a","b","c","d","e","f",];

% Initialize figure
fig2 = figure(2);clf
fig2.Position = [200 200 [750 450]*2.2]; 
ha = tight_subplot(2,3,[.1 .05],[.1 .1],[.06 .1]);

% Loop through error scenarios
for idxScen = 1:cntScen
    axes(ha(idxScen)); hold on; %subplot(2,3,idxScen);
    hold on; xlim([-90,20]); ylim([-50,2]); grid on; set(gcf, 'color', 'white');
    set(gca,'DefaultLineLineWidth',2)
    xlabel('\theta (deg)')
    ylabel('Power (dB)')
    set(gca,'FontSize',9)
    yticks(-60:10: 0)
    xticks(-90:30:90)

    % Loop through cases of module counts
    for idxCase = 1:cntCase
        plot(degs, db(P_err(:,idxCase,idxScen),'power'),'-',...
                                'Color',lineColor(idxCase,:),...
                                'linewidth', lineWidth(idxCase),...
                    'DisplayName',sprintf("Upper bound (%0.0f block%s of %0.0f)",...
                         N(idxCase),plural((N(idxCase)>1)+1),M/N(idxCase)));
    end
    plot(degs, db(P_nom(:,idxScen))/2, 'color',[0.1,0.75,0.1],...
                'DisplayName', 'Nominal beampattern');
            
	set(gca, 'fontsize', 14)
end

legend('location','northwest')
legend boxoff

% Plot figure indexes
% annotation('textbox', [ha(1).Position(1:2),0,0] + [0.21 0.28 0.03 0.03] , 'String', "(a)",'BackgroundColor','w', 'FitBoxToText','off', 'FontSize', 14, 'Fontweight', 'bold', 'Facealpha', 0.85);
% annotation('textbox', [ha(2).Position(1:2),0,0] + [0.21 0.28 0.03 0.03] , 'String', "(b)",'BackgroundColor','w', 'FitBoxToText','off', 'FontSize', 14, 'Fontweight', 'bold', 'Facealpha', 0.85);
% annotation('textbox', [ha(3).Position(1:2),0,0] + [0.21 0.28 0.03 0.03] , 'String', "(c)",'BackgroundColor','w', 'FitBoxToText','off', 'FontSize', 14, 'Fontweight', 'bold', 'Facealpha', 0.85);
% annotation('textbox', [ha(4).Position(1:2),0,0] + [0.21 0.28 0.03 0.03] , 'String', "(d)",'BackgroundColor','w', 'FitBoxToText','off', 'FontSize', 14, 'Fontweight', 'bold', 'Facealpha', 0.85);
% annotation('textbox', [ha(5).Position(1:2),0,0] + [0.21 0.28 0.03 0.03] , 'String', "(e)",'BackgroundColor','w', 'FitBoxToText','off', 'FontSize', 14, 'Fontweight', 'bold', 'Facealpha', 0.85);
% annotation('textbox', [ha(6).Position(1:2),0,0] + [0.21 0.28 0.03 0.03] , 'String', "(f)",'BackgroundColor','w', 'FitBoxToText','off', 'FontSize', 14, 'Fontweight', 'bold', 'Facealpha', 0.85);

annotation('textbox', [ha(1).Position(1:2),0,0] + [0.21 0.28 0.03 0.03] , 'String', "(a)",'BackgroundColor','w', 'FitBoxToText','off', 'FontSize', 14, 'Facealpha', 0.85);
annotation('textbox', [ha(2).Position(1:2),0,0] + [0.21 0.28 0.03 0.03] , 'String', "(b)",'BackgroundColor','w', 'FitBoxToText','off', 'FontSize', 14, 'Facealpha', 0.85);
annotation('textbox', [ha(3).Position(1:2),0,0] + [0.21 0.28 0.03 0.03] , 'String', "(c)",'BackgroundColor','w', 'FitBoxToText','off', 'FontSize', 14, 'Facealpha', 0.85);
annotation('textbox', [ha(4).Position(1:2),0,0] + [0.21 0.28 0.03 0.03] , 'String', "(d)",'BackgroundColor','w', 'FitBoxToText','off', 'FontSize', 14, 'Facealpha', 0.85);
annotation('textbox', [ha(5).Position(1:2),0,0] + [0.21 0.28 0.03 0.03] , 'String', "(e)",'BackgroundColor','w', 'FitBoxToText','off', 'FontSize', 14, 'Facealpha', 0.85);
annotation('textbox', [ha(6).Position(1:2),0,0] + [0.21 0.28 0.03 0.03] , 'String', "(f)",'BackgroundColor','w', 'FitBoxToText','off', 'FontSize', 14, 'Facealpha', 0.85);

% Plot array insets
% Adjust curvature and element diameter in all arrays
boxPos = [[0.10,0.73,0.12,0.14];[0.10,0.29,0.12,0.14]];
idxScens = [1,4];
for idx = 1:2
    idxScen = idxScens(idx);
    array.ElDiameterRatio = elDiameter(idxScen);
    array.Curvature = curvature(idxScen);
    array.PosXError = elPosXErr(idxScen);
    block.Curvature = curvature(idxScen);
    block.PosXError = blPosXErr(idxScen);
    block.OrientError = blOrientErr(idxScen);

    axes('Position',boxPos(idx,:))
    box on
    set(gca,'DefaultLineLineWidth',2)
    xlabel('X [m]')
    ylabel('Y [m]')
    set(gca,'FontSize',11)
    set(gcf, 'color', 'white');
    grid on; 
    bp_nom.plotArrays('PlotArrow',0)
    ax=gca();
    for idxC = 1:length(ax.Children)
        if strcmp(class(ax.Children(idxC)),"matlab.graphics.chart.primitive.Line")
            ax.Children(idxC).Color='b';
        end
    end
    ylim([-0.2,0.2])
    ax.Title.Position = [0 0.14 0];
    ax.Title.FontWeight = 'normal';
    ax.XLabel.Position = [0 -0.14 -4];
    ax.YLabel.Position = [-0.22 .1 -4];
    
end
%%
exportgraphics(gcf,'fig7.pdf')%,'Resolution',300)      
                          