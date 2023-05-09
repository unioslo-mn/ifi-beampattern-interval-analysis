function plotArrays(obj,optional)

% Plot (sub-)arrays and block array (if given)
%
% This function plots the given arrays offset by the block 
% element positions and the block array itself. 
% (See biat.SensorArray.plot for further details.)
%______________________________________________________________________
% USAGE        
%   obj.plotArrays();
% _________________________________________________________________________
% NECESSARY ARGUMENT
%   obj       : object of biat.BeamPattern type
%   PlotArrow
% _________________________________________________________________________
% OPTIONS
%   PlotArrow : Boolean to switch arrow plot on (default) and off
% _________________________________________________________________________
% EXAMPLES
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
           optional.PlotArrow   (1,1)  {mustBeNumericOrLogical} = 1
    end

    if isempty(obj.Block)
        obj.Arrays.plot('PlotArrow',optional.PlotArrow);
    else
        obj.Arrays.plot('Block',obj.Block,'PlotArrow',optional.PlotArrow);
    end
    
end

