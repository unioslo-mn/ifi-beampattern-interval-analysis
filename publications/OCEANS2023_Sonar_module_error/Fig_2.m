clear

% Array parameters
M = 2;
N = 2;

% Error parameters
errGain = 0.1;
errPha = deg2rad(5);
errPsi = deg2rad(1);
errPosBl = 5e-4;
errPsiBl = deg2rad(1);

% Beamforming parameters
incAngle = deg2rad(-11);
AngResDeg = 1;
PolTol = 5e-4;
PolIncl = true;

% Generate array and beampattern objects
% Initialize array
array = biat.SensorArray(       'ElCount',          M,...
                                'ElPitchRatio',     0.5,...
                                'ElDiameterRatio',  0,...
                                'Curvature',        0,...
                                'TaperType',        'none',...
                                'GainError',        errGain,...
                                'PhaseError',       errPha,...
                                'PosXError',        0,...   
                                'PosYError',        0,...
                                'OrientError',      errPsi,...
                                'CouplingCoeff',    0,...
                                'SteeringAngle',    0);

% Generate block
block = copy(array);
block.ElCount = N;
block.ElPitchRatio= 0.5*M;
block.ElDiameterRatio = 0;
block.GainError = 0;
block.PhaseError = 0;
block.OrientError = errPsiBl;
block.PosXError = errPosBl;
block.PosYError = errPosBl;
                            
% Calculate nominal beampattern with tapering
bp_pol = biat.BeamPattern(array,'polygonal','Block',block,...
                                        'BlockType','repeating',...
                                        'BeamResolutionDeg',AngResDeg,...
                                        'PolygonTolerance',PolTol);

%% Get sub-array intervals without orientation error for the angles within

% Get beam indexes of incidence angles within the orientation error
cntBeam = bp_pol.BeamCount;
[~,idxBeam] = min(abs(bp_pol.BeamAngles - incAngle));
difBeam = ceil(bp_pol.Block.OrientError / bp_pol.BeamResolution);
intBeam = max(1,idxBeam-difBeam) : min(idxBeam+difBeam,cntBeam);  

% Calculate sub-array itervals (without block orientation error)
bp_pol.Block.OrientError = 0;
elInt(1:M,1:length(intBeam)) = ciat.PolygonalInterval(0);
sumInt(1:length(intBeam)) = ciat.PolygonalInterval(0);
for idx = 1:length(intBeam)
	bp_pol.BeamIndex = intBeam(idx);
    elInt(:,idx) = bp_pol.ElementIntervals(:,1);
    sumInt(idx) = bp_pol.SubArrayInterval(1);
end
bp_pol.Block.OrientError = errPsiBl;
bp_pol.BeamIndex = idxBeam;

%%
% Plot array
% figure(1);clf;hold on
% bp_pol.plotArrays

fontSize = 13;
lineColors = [  117,112,179;...
                27,158,119;...
                217,95,2;...
                0 0 0;...
                0 0 0]/255;
            
% Plot complex intervals
fig = figure(2);clf;hold on
fig.Position = [200 200 750 450]; 
bp_pol.ElementIntervals.plot('o-','Color',lineColors(1,:),...
                               'LineWidth',3,'MarkerSize',3,...
                              'DisplayName','Element (El)');
elInt(:,1).plot('-','Color',lineColors(1,:),...
                                'LineWidth',0.5,'DisplayName','Element (El)');
elInt(:,3).plot('--','Color',lineColors(1,:),...
                                'LineWidth',0.5,'DisplayName','Element (El)');
bp_pol.SubArrayInterval.plot('o-','Color',lineColors(2,:),...
                                'LineWidth',2,'MarkerSize',3,...
                              'DisplayName','Array (Ar)');
sumInt(1).plot('-','Color',lineColors(2,:),...
                                'LineWidth',0.5,'DisplayName','Array (Ar)');
sumInt(2).plot('--','Color',lineColors(2,:),...
                                'LineWidth',0.5,'DisplayName','Array (Ar)');
sumInt(3).plot('-.','Color',lineColors(2,:),...
                                'LineWidth',0.5,'DisplayName','Array (Ar)');
bp_pol.BlockInterval.plot('s-','Color',lineColors(3,:),...
                                'LineWidth',3,'MarkerSize',3,...
                              'DisplayName','Block (Bl)');
bp_pol.ArrayInterval.plot('o:','Color',lineColors(4,:),...
                                'LineWidth',3,'MarkerSize',3,...
                              'DisplayName','Module (Mo)');
bp_pol.TotalInterval.plot('o-','Color',lineColors(5,:),...
                                'LineWidth',3,'MarkerSize',3,...
                              'DisplayName','Beamsum (B)');
grid on
axis equal
legend(legendUnq(fig))

xlim([0.25,1.38])

text(0.55,0.16,'El^I_1(\alpha)',...
                    'FontSize',fontSize,'Color',lineColors(1,:))
text(0.55,-0.16,'El^I_2(\alpha)',...
                    'FontSize',fontSize,'Color',lineColors(1,:))
text(1.09,0.00,'Ar^I=\cup_\alpha(El^I_1 +El^I_2)',...
                    'FontSize',fontSize,'Color',lineColors(2,:))
text(0.36,0.26,'Bl^I_1',...
                    'FontSize',fontSize,'Color',lineColors(3,:))
text(0.36,-0.26,'Bl^I_2',...
                    'FontSize',fontSize,'Color',lineColors(3,:))
text(0.49,0.29,'Mo^I_1=Bl_1^I\timesAr^I',...
                    'FontSize',fontSize,'Color',lineColors(4,:))
text(0.49,-0.29,'Mo^I_2=Bl_2^I\timesAr^I',...
                    'FontSize',fontSize,'Color',lineColors(4,:))
text(0.65,0,'B^I(\theta)=Mo_1^I+Mo_2^I',...
                    'FontSize',fontSize,'Color',lineColors(5,:))

xlabel('Real')
ylabel('Imaginary')
set(gca, 'fontsize', 14)
legend boxoff
%%
% exportgraphics(gcf,'fig2.pdf')%,'Resolution',300)