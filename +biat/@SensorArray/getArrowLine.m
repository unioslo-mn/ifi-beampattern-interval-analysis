function [x,y] = getArrowLine(r1,r2,L)

% Create and arrow out of lines (similar to the quiver function)
%
% This function creates an arrow out of lines between two points
% (r1,r2) in a 2-D plane (x,y) with an adjustable head size (L).
% This is a static function that can be called from a SensorArray
% object or directly from the ciat.SensorArray class.
%______________________________________________________________________
% USAGE        
%   [x,y] = obj.getArrowLine(r1,r2,L);
% _________________________________________________________________________
% NECESSARY ARGUMENT
%   r1       : coordinates of the first point as a two-element array
%   r2       : coordinates of the second point as a two-element array
%   L        : absolute length of the arrow head
% _________________________________________________________________________
% OPTIONS
% _________________________________________________________________________
% EXAMPLES
%   [x,y] = biat.SensorArray.getArrowLine([0,0],[1,1],0.2)
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
       r1       (1,2)   {mustBeNumeric}
       r2       (1,2)   {mustBeNumeric}
       L        (1,1)   {mustBeNumeric}
    end

    n = (r2-r1)/norm(r2-r1);
    
    th0 = atan2(n(2),n(1));
    
    theta = 25; % to rotate 90 counterclockwise
    
    th1 = 180-theta;
    th2 = -180+theta;
    R1 = [cosd(th1) -sind(th1); sind(th1) cosd(th1)]';
    R2 = [cosd(th2) -sind(th2); sind(th2) cosd(th2)]';
    d1 = (L*R1*n')';
    d2 = (L*R2*n')';
    
    %atan2(y,x);
    r3 = r2+d1;
    r4 = r2;
    r5 = r2+d2;
    r6 = r2;
    
    
    x = [r1(1), r2(1), r3(1), r4(1),r5(1),r6(1)];
    y = [r1(2), r2(2), r3(2), r4(2),r5(2),r6(2)];
end
