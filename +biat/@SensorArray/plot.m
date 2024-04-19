function plot(obj,optional)  

% Plot (sub-)arrays and block array (if given)
%
% This function plots the given arrays offset by the block 
% element positions and the block array itself. 
%______________________________________________________________________
% USAGE        
%   obj.plotArrays();
% _________________________________________________________________________
% NECESSARY ARGUMENT
%   obj       : object of biat.BeamPattern type
% _________________________________________________________________________
% OPTIONS
%   Block     : Block array from class biat.SensorArray (default empty)
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
    optional.Block  (:,:)  biat.SensorArray  = biat.SensorArray.empty
    optional.PlotArrow (1,1)  {mustBeNumericOrLogical} = 1
end

    hold on
    set(gca,'LineWidth',2)
    set(gca,'DefaultLineLineWidth',2)
    set(gca,'defaultAxesFontSize',20)
    title('Physical array')

    % Plot each array
    if isempty(optional.Block)
        for iArr = 1:length(obj)
            plotArr(obj(iArr),0,0,0,0,optional.PlotArrow);
        end
    else
        for iArr = 1:length(obj)
            plotArr(obj(iArr),  optional.Block.ElPosX(iArr),...
                                optional.Block.ElPosY(iArr),...
                                optional.Block.ElOrient(iArr),0,...
                                optional.PlotArrow);
        end
        plotArr(optional.Block,0,0,0,1,optional.PlotArrow);
    end
end

%%

function plotArr(obj, posXOffs, posYOffs, psiOffs, isBlock, plotArrow)

    % Shift element positions and orientation according to block offset
    posX = obj.ElPosX * cos(-psiOffs) - obj.ElPosY * sin(-psiOffs) + posXOffs;
    posY = obj.ElPosX * sin(-psiOffs) + obj.ElPosY * cos(-psiOffs) + posYOffs;    
    psi  = obj.ElOrient + psiOffs;
    psiI = obj.OrientInterval + psiOffs;

    % Plot array element positions
    if obj.ElDiameter == 0
        for elem = 1:obj.ElCount
            x = posX(elem) - obj.PosXError + obj.PosXBias;
            y = posY(elem) - obj.PosYError + obj.PosYBias;
            w = 2*obj.PosXError;
            h = 2*obj.PosYError;
            rectangle('Position',[x,y,w,h],'FaceColor',[1 .5 .5])
            
            if ~isBlock
                plot(posX, posY,'ko');
            else
                plot(posX, posY,'kx');
            end
        end
    else
        for elem = 1:obj.ElCount
            x = posX(elem) - obj.PosXError - obj.PosXBias;
            y = posY(elem) - obj.PosYError - obj.PosYBias;
            w = 2*obj.PosXError;
            h = 2*obj.PosYError;
            rectangle('Position',[x,y,w,h],'FaceColor',[1 .5 .5])
            
            x1 = posX(elem)  - ...
                 obj.ElDiameter/2 * cos(psi(elem));
            x2 = posX(elem) +...
                 obj.ElDiameter/2 * cos(psi(elem));
            y1 = posY(elem) -  ...
                 obj.ElDiameter/2 *sin(-psi(elem));
            y2 = posY(elem) +  ...
                 obj.ElDiameter/2 *sin(-psi(elem));
            plot([x1,x2],[y1,y2], 'k');
        end
    end
    
    % Plot array element orientations
    if plotArrow
        quiver(posX, posY, obj.ElPitch .* sin(psiI.Infimum), ...
                           obj.ElPitch .* cos(psiI.Infimum), ...
                                    'off','c');
        quiver(posX, posY, obj.ElPitch .* sin(psiI.Supremum), ...
                           obj.ElPitch .* cos(psiI.Supremum), ...
                                    'off','c');
        quiver(posX, posY, obj.ElPitch .* sin(psi),  ...
                           obj.ElPitch .* cos(psi), ...
                                    'off','b');
    end
    
       
    % Plot steering direction
    if plotArrow
        [x,y] = obj.getArrowLine([posXOffs,posYOffs],...
                    [obj.WaveLength * sin(obj.SteeringAngle + psiOffs ) + posXOffs, ...
                     obj.WaveLength * cos(obj.SteeringAngle + psiOffs) + posYOffs],...
                                obj.WaveLength/10);
        plot(x,y,'k-');
    end
    
    axis('equal')
end