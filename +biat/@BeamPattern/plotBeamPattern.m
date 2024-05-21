function plt = plotBeamPattern(obj,options)

% Calculates (if not given) and plots the beampattern
%
% This function calculates the beampattern and plots it 
% according to the default or given line style. If a 
% previously calculated beampattern is given as an 
% optional argument, the function skips the calculation
% only plots the given beampattern. It also returns an
% array of handlers to the plots.
%______________________________________________________________________
% USAGE        
%   plt = obj.plotBeamPattern(Name,Value);
% _________________________________________________________________________
% NECESSARY ARGUMENT
%   obj       : object of biat.BeamPattern type
% _________________________________________________________________________
% OPTIONS
%   BeamPattern: Numeric array of 3 columns with the 
%                number of rows equal to object's BeamCount
%   Color      : same as the Color property of the Matlab 
%                plot function
%   LineWidth  : same as the LineWidth property of the Matlab 
%                plot function
%   LineStyle  : same as the LineStyle property of the Matlab 
%                plot function
%   Backtracking: set to "1" to plot backtracked beampattern
%                (default), and "0" to skip it
% _________________________________________________________________________
% EXAMPLES
%   plt = obj.beampattern();
%   plt = obj.beampattern('BeamPattern',P_nom);
%   plt = obj.plotBeamPattern('Color',[0.5,1.0,0.0],...
%                             'Backtracking',1)]; 
% _________________________________________________________________________
%
% Copyright (C) 2023 H. Arnestad and G. Gereb, BSD-3
% If you use this software, please cite it as in CITATION.cff
% Project: Beampattern Interval Analysis 
% Website: doi.org/10.5281/zenodo.6856232
% Contact: haavaarn@uio.no, gaborge@uio.no
% (More information in README.md and LICENSE.md.)
% _________________________________________________________________________


arguments
   obj
   options.BeamPattern  (:,:) {mustBeNumeric} = []
   options.Color        (1,:) = 'b'
   options.LineWidth    (1,1) = 2
   options.LineStyle    (1,:) = '-'
   options.Backtracking (1,1) = 1
end

% Initialize plot
hold on;
xlabel('sin(angle) [-]')
xlim([-1,1])
ylim([-60,0])
ylabel('Amplitude [dB]')

% Calculatebeampattern
degs = sin(obj.BeamAngles);
if isempty(options.BeamPattern)
    P = obj.calculateBeamPattern;
else
   P = options.BeamPattern;
   assert(size(P,1) == obj.BeamCount)
   assert(size(P,2) == 2 || size(P,2) == 3)
end

P(P==0) = eps;

% Plot beampattern
switch obj.Type 
    case 'nominal'
        plt(3) = plot(degs, db(P(:,1),'power') ,'-',...
                                    'Color',options.Color,...
                                   'LineWidth',options.LineWidth,...
                                   'LineStyle',options.LineStyle,...
                                   'DisplayName','nom');
        plt(1:2) = plot(degs, db(P(:,[2,3]),'power') ,':',...
                                    'Color',options.Color,...
                                   'LineWidth',options.LineWidth,...
                                   'LineStyle',options.LineStyle,...
                                   'DisplayName','apx');
    case 'rectangular'
        plt(3) = plot(0,0,'.','MarkerSize',0.1,'DisplayName','init');
        plt(1:2) = plot(degs, db(P,'power') ,'-',...
                                    'Color',options.Color,...
                                   'LineWidth',options.LineWidth,...
                                   'LineStyle',options.LineStyle,...
                                   'DisplayName','rect');
    case 'circular'
        plt(3) = plot(0,0,'.','MarkerSize',0.1,'DisplayName','init');
        plt(1:2) = plot(degs, db(P,'power') ,'-',...
                                    'Color',options.Color,...
                                   'LineWidth',options.LineWidth,...
                                   'LineStyle',options.LineStyle,...
                                   'DisplayName','circ');
    case 'polygonal'
        if options.Backtracking == 1
            if obj.Arrays(1).OrientError == 0   && ...
               (isempty(obj.Block) || obj.Block.OrientError == 0)
                P_btr = obj.backtrackBeamPattern;
                plt(3) = plot(degs, db(P_btr,'power') ,'--',...
                                    'Color',options.Color,...
                                   'LineWidth',options.LineWidth,...
                                   'LineStyle',options.LineStyle,...
                                       'DisplayName','back');
            else
                plt(3) = plot(0,0,'.','MarkerSize',0.1,'DisplayName','init');
                warning(['Backtracking is not implemented for arrays ',...
                        'with non-zero orientation error.']) 
            end
        else
            plt(3) = plot(0,0,'.','MarkerSize',0.1,'DisplayName','init');
        end
        plt(1:2) = plot(degs, db(P,'power') ,'-',...
                                'Color',options.Color,...
                               'LineWidth',options.LineWidth,...
                               'LineStyle',options.LineStyle,...
                                   'DisplayName','poly');
    otherwise
        warning('no beampattern plot function implement for this Type')
end

plt = plt(:);

% Deinitialize plot
grid on
%     hold off

end