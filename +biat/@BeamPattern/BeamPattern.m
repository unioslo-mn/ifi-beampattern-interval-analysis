
classdef BeamPattern
    
% BeamPattern class for beampatter interval analysis
%
% This is a class of the Beampattern Interval Analysis Toolbox.
% It allows the calculation of the nominal beampattern and
% beampattern tolerance interval of the given array(s) 
% affected by the given errors. The errors can be represented 
% as one of the available complex interval types.
% _________________________________________________________________________
%
% Copyright (C) 2023 H. Arnestad and G. Gereb, BSD-3
% If you use this software, please cite it as in CITATION.cff
% Project: Beampattern Interval Analysis 
% Website: doi.org/10.5281/zenodo.6856232
% Contact: haavaarn@uio.no, gaborge@uio.no
% (More information in README.md and LICENSE.md.)
% _________________________________________________________________________

    properties
        ArrayHandle         % Handle to the biat.SensorArray object(s) given during the construction (used to generate Arrays)
        Block               % A biat.SensorArray object defining the sub-array positions and errors
        BlockType           % An optional selector that allows defining split or repeating blocks (see constructor function)
        Type                % Type of complex interval representation
        BeamResolutionDeg   % Beam angle resolution in degrees(!) used to determine the BeamAngles array for sampling the beampattern
        BeamIndex           % Index of the selected beam out of the BeamAngles array
        BeamAngle           % Custom beam angle given in radians, works when BeamIndex is "0"
        PolygonTolerance    % Precision parameter for polygonal interval calculations (see ciat.Polygonal class)
    end

    properties (Dependent)
        Arrays              % One or more biat.SensorArray objects defining the element positions and errors
        ArrayCount          % Number of arrays (should be the same as the block element count)
        BeamResolution      % Beam resolution in radians [rad]
        BeamAngles          % Beam angles to sample the beampattern at [rad]
        BeamCount           % Number of beam angles
        IncidenceAngle      % Incidence angle used for the array response calculation (it is the BeamIndex-th element of BeamAngles or equals to BeamAngle if BeamIndex is "0") [rad]
        IncidenceVector     % Incidence vector used for the array reponse calculation (calculated from the Incidenceangle and array WaveLength) [rad/m]
        PhaseIntervals      % Steered element phase error intervals [rad]
        AmplitudeIntervals  % Steered element amplitude error intervals
        ElementIntervals    % Complex sensitivity intervals of the elements
        SubArrayInterval    % Complex sensitivity interval of each sub-array (sum of element intervals)
        BlockInterval       % Complex sensitivity intervals of the blocks elements
        ArrayInterval       % Complex sensitivity interval of each array (including all block errors)
        TotalInterval       % Complex sensitivity interval of the system (sum of array intervals)
        PowerInterval       % Total output power interval (absolute value square of the array interval)
        MaxPowerPoints      % Element contributions to array maximum power point
        MinPowerPoints      % Element contributions to array minimum power point
        MaxErrorPattern     % Error pattern resulting maximum array power
        MinErrorPattern     % Error pattern resulting minimum array power
    end

    methods
        %% Define class constructor
        function obj = BeamPattern(ArrayHandle, Type, optional)
            % BEAMPATTERN Construct an instance of this class
            %
            % Construct a biat.BeamPattern type object to calculate the 
            % beampattern of the given SensorArray objects using one of the
            % available complex interval arithmetic types (rectangular,
            % circular, polygonal). If no block array is defined, the
            % object will calculate the array interval of each array in the
            % Arrays property and then sum them to get the total interval.
            % If a block array is given, the constructor will check if the
            % number of block elements agrees with the number of arrays in
            % the input. If the BlockType is 'custom', then these have to
            % agree, if is is 'repeating' or 'split', then only one array
            % has to be given, which is then repeated the number of block
            % element numbers, or split into block element parts.
            %__________________________________________________________________________
            % USAGE        
            %   bp = biat.BeamPattern(Arrays, Type, Name, Value) 
            % _________________________________________________________________________
            % NECESSARY ARGUMENT
            %   Arrays       : array of object(s) of biat.SensorArray type
            %   Type         : complex interval type 
            %                 (rectangular, circular, polygonal)
            % _________________________________________________________________________
            % OPTIONS
            %   Block       : a biat.SensorArray type object
            %   BlockType   : block handling method 
            %                 (custom, repeating, split)
            %   BeamResolutionDeg: Beam angle resolution in degrees(!) 
            %                      used to determine the BeamAngles array 
            %                      for sampling the beampattern
            %   BeamIndex   :Index of the selected beam out of the 
            %                BeamAngles array
            %   BeamAngle   :Custom beam angle given in radians, 
            %                works when BeamIndex is "0"
            %   PolygonTolerance: Precision parameter for polygonal 
            %                     interval calculations 
            %                     (see ciat.Polygonal class)
            % _________________________________________________________________________
            % EXAMPLES
            %   bp_nom = biat.BeamPattern(array,'nominal',...
            %                            'BeamResolutionDeg',0.35/2);
            %   bp_pol = biat.BeamPattern(array,'rectangular',...
            %                            'BeamResolutionDeg',0.35/2);
            %   bp_pol = biat.BeamPattern(array,'circular',...
            %                            'BeamResolutionDeg',0.35/2);
            %   bp_pol = biat.BeamPattern(array,'polygonal',...
            %                            'BeamResolutionDeg',0.35/2,...
            %                            'PolygonTolerance',1e-6);
            %   bp_pol = biat.BeamPattern(array,'polygonal',...
            %                            'Block',block,...
            %                            'BlockType','repeating',...
            %                            'BeamResolutionDeg',0.5,...
            %                            'PolygonTolerance',1e-6);
            % _________________________________________________________________________

            arguments
               ArrayHandle (1,:)  biat.SensorArray   = biat.SensorArray
               Type        (1,1)  string {mustBeMember(Type ,...
                                      {'nominal','rectangular',...
                                       'circular','polygonal'})} ...
                                                        = 'rectangular'
               optional.Block               (:,:)  biat.SensorArray  ...
                                                    = biat.SensorArray.empty
               optional.BlockType           (1,1)  string ...
                                    {mustBeMember(optional.BlockType ,...
                                    {'custom','repeating','split'})} ...
                                                            = 'custom'
               optional.BeamResolutionDeg   (1,1)  {mustBeNumeric} = 1
               optional.BeamIndex           (1,1)  {mustBeNumeric} = 1
               optional.BeamAngle           (1,1)  {mustBeNumeric} = 0
               optional.PolygonTolerance    (1,1)  {mustBeNumeric} = 1e-5
            end
            obj.ArrayHandle         = ArrayHandle;
            obj.Type                = Type;
            obj.Block               = optional.Block;
            obj.BlockType           = optional.BlockType;
            obj.BeamResolutionDeg   = optional.BeamResolutionDeg;
            obj.BeamIndex           = optional.BeamIndex;
            obj.BeamAngle           = optional.BeamAngle;
            obj.PolygonTolerance    = optional.PolygonTolerance;


            % Forcing block array structure
            assert(range([obj.ArrayHandle.WaveLength]) == 0)
            assert(range([obj.ArrayHandle.ElCount]) == 0);
            assert(range([obj.ArrayHandle.SteeringAngle]) == 0)
            assert(range([obj.ArrayHandle.Curvature]) == 0);
            if ~isempty(obj.Block)
                assert(obj.ArrayHandle(1).WaveLength == obj.Block.WaveLength)
                assert(obj.ArrayHandle(1).SteeringAngle == obj.Block.SteeringAngle)
                assert(obj.ArrayHandle(1).Curvature == obj.Block.Curvature)
                obj.Block.ElDiameterRatio = 0;
                
                % Check block size and number of arrays
                N = length(obj.ArrayHandle);
                Nb = obj.Block.ElCount;
                switch obj.BlockType
                    case 'custom'
                        if N ~= Nb
                            error(['Number of arrays incompatible with ',...
                                    'the element number of the block.'])
                        end
                    case 'repeating'
                        if N == 1
                            obj.Block.TaperType = 'custom';
                            obj.Block.TaperParam = ones(Nb,1)/Nb;
                        else
                           error(['Only one array can be given in ',... 
                                 'repeating block mode'])
                        end
                    case 'split'
                        if N == 1
                            obj.Block.TaperType = 'custom';
                            obj.Block.TaperParam = ones(Nb,1);
                        else
                            error(['Only one array can be given in ',... 
                                 'split block mode'])                   
                        end
                end
            end
        end

        %% Define get methods for dependent properties

        % Arrays
        function value = get.Arrays(obj)
            % Default case
            value = obj.ArrayHandle;

            % Exception
            if ~isempty(obj.Block) && ~strcmp(obj.BlockType,'custom')
                M = obj.ArrayHandle.ElCount;
                Nb = obj.Block.ElCount;
                switch obj.BlockType
                    case 'repeating'
                        value(1:Nb) = biat.SensorArray;
                        for n = 1:Nb
                            value(n) = copy(obj.ArrayHandle);
                            value(n).SteeringAngle = obj.Block.SteeringAngle - ...
                                                    obj.Block.ElOrient(n); 
                        end
                    case 'split'
                        Mb = M / Nb;
                        w_apod = obj.ArrayHandle.TaperWeights;
                        value(1:Nb) = biat.SensorArray;
                        for n = 1:Nb
                            value(n) = copy(obj.ArrayHandle);
                            value(n).ElCount = Mb;
                            value(n).TaperType = 'custom';
                            value(n).TaperParam = w_apod((1:Mb) + (n-1)*Mb);
                            value(n).SteeringAngle = obj.Block.SteeringAngle - ...
                                                    obj.Block.ElOrient(n);
                        end
                end
            end
        end
        
        % Number of arrays
        function value = get.ArrayCount(obj)
           value = length(obj.Arrays); 
        end
        
        % Beam resolution (in radians)
        function value = get.BeamResolution(obj)
            value = deg2rad(obj.BeamResolutionDeg)';
        end

        % List of beam angles
        function value = get.BeamAngles(obj)
            value = ( -pi/2 : obj.BeamResolution : pi/2 )';
        end
        

        % Count of beam angles
        function value = get.BeamCount(obj)
            value = length(obj.BeamAngles);
        end

        % Incidence angle of selected beam
        function value = get.IncidenceAngle(obj)
            if obj.BeamIndex > 0
                value = obj.BeamAngles(obj.BeamIndex);
            else
                value = obj.BeamAngle;
            end
        end

        % Incidence vector of selected beam
        function value = get.IncidenceVector(obj)
            value = 2*pi / obj.Arrays(1).WaveLength * ...
                        ([sin(obj.IncidenceAngle); cos(obj.IncidenceAngle)] );
        end

        % Array element intervals (annular sectors or polar intervals)
        function value = get.ElementIntervals(obj)
            value = obj.getElementIntervals();
        end
        
        % Sub-array intervals (sum of array element intervals)
        function value = get.SubArrayInterval(obj)
            value = obj.getSubArrayInterval();
        end
        
        % ArrayInterval
        function value = get.ArrayInterval(obj)
            value = obj.getArrayInterval;
        end
        
        % BlockInterval
        function value = get.BlockInterval(obj)
            if ~isempty(obj.Block)
                value = obj.getElementIntervals('getBlock',1);
            else
               value = [];
            end
        end
        
        % TotalInterval
        function value = get.TotalInterval(obj)
            value = obj.getTotalInterval;
        end

        % PowerInterval
        function value = get.PowerInterval(obj)
            if strcmp(obj.Type,"nominal") 
                value = abs(obj.TotalInterval)^2;
            else
                absInt = abs(obj.TotalInterval);
                value = absInt * absInt ;
            end 
        end

        % Element contributions backtracked from maximum array power point
        function value = get.MaxPowerPoints(obj)
            value = obj.backtrackMinMaxPower('getMax',1);
        end

        % Element contributions backtracked from minimum array power point
        function value = get.MinPowerPoints(obj)
            value = obj.backtrackMinMaxPower('getMax',0);
        end

        % Error pattern for maximum array power
        function value = get.MaxErrorPattern(obj)
            k = obj.IncidenceVector;
            posX = obj.Arrays.ElPosX;
            posY = obj.Arrays.ElPosY;
            [E_extr,~] = obj.backtrackComplexIntervals('getMax',1);
            value = E_extr .* exp(-1j*( k(1)*posX + k(2)*posY));
        end

        % Error pattern for minimum array power
        function value = get.MinErrorPattern(obj)
            if ~ isempty(obj.MinPowerPoints)
                k = obj.IncidenceVector;
                posX = obj.Arrays.ElPosX;
                posY = obj.Arrays.ElPosY;
                [E_extr,~] = obj.backtrackComplexIntervals('getMax',0);
                value = E_extr .* exp(-1j*( k(1)*posX + k(2)*posY));
            else
                value = [];
                warning('Minimum error pattern cannot be backtracked')
            end
        end


        %% Define method signatures
        value = getElementIntervals(obj,options)
        value = getSubArrayInterval(obj)
        value = getArrayInterval(obj)
        value = getTotalInterval(obj)
        value = calculateBeamPattern(obj)
        value = backtrackMinMaxPower(obj,options)
        [E_extr,C_extr] = backtrackComplexIntervals(obj,options)
        value = backtrackBeamPattern(obj,options)
        plotComplexIntervals(obj,options)
        plt = plotBeamPattern(obj,options)
        plotErrorPattern(obj,options)
        plotArrays(obj)

    end

end

