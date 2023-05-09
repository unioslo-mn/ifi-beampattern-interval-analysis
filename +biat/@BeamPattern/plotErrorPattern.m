function plotErrorPattern(obj,options)

% Plot backtracked error pattern (min or max power)
%
% This function plots the backtracked error pattern
% contributing to the minimium or maximum total
% response power as one bar plot for the amplitude
% of the backtracked complex points and one bar plot
% for their phase values.
%______________________________________________________________________
% USAGE        
%   obj.plotErrorPattern(Name,Value)
% _________________________________________________________________________
% NECESSARY ARGUMENT
%   obj       : object of biat.BeamPattern type
% _________________________________________________________________________
% OPTIONS
%   getMax     : "1" to backtrack maximum power (default),
%                 "0" to backtrack minimum power point.
% _________________________________________________________________________
% EXAMPLES
%   points = obj.plotErrorPattern;
%   points = obj.plotErrorPattern('getMax',0);
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
       options.getMax       (1,1)   {mustBeNumeric} = 1
    end

    if strcmp(obj.Type,'polygonal')
        if obj.Arrays.OrientError == 0
            M = obj.Arrays.ElCount;
            if options.getMax
                errPat = obj.MaxErrorPattern;
            else
                errPat = obj.MinErrorPattern;
            end
            I_g = obj.Arrays.GainInterval;
            I_phi = obj.Arrays.PhaseInterval;

            % Initialize figure
            hold on; 
            title('Backtracked element errors')
            xlabel('Element {\itm}')
            xlim([0.5, M+0.5]);xticks([0 : 10 : M])
            sep = 1.2; eps = 1e-3;
            yyaxis left
            ylabel('Amplitude error (%)')

            % Plot error amplitudes
            bar(100*(abs(errPat)-1),'FaceColor','b','EdgeColor','b')
            plot(1:M , 100*([I_g.Supremum]'-1),'b--')
            plot(1:M , 100*([I_g.Infimum]'-1),'b--')
            ylim(100*[(I_g(1).Infimum-1)*sep-eps, (I_g(1).Supremum-1)*sep+eps]); 
            yl=ylim();ylim([2*yl(1)-yl(2) yl(2)])

            % Plot error phases
            yyaxis right
            ylabel('Phase error (deg)')
            bar(rad2deg(angle(errPat)),'FaceColor','r','EdgeColor','r')
            plot([1:M],rad2deg([I_g.Supremum]'),'r--')
            plot([1:M],rad2deg([I_g.Infimum]'),'r--')
            ylim([rad2deg(I_phi(1).Infimum)*sep-eps, ...
                  rad2deg(I_phi(1).Supremum)*sep+eps]); 
            yl=ylim();ylim([yl(1) 2*yl(2)-yl(1)])

            % Release figure hold
            hold off
        else
            warning(['Backtracking is not implemented for arrays ',...
                        'with non-zero orientation error.']) 
        end
    else
        warning('No backtracking available for this type.')
    end
end

