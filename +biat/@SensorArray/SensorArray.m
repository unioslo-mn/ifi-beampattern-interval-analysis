classdef SensorArray < matlab.mixin.Copyable % handle
% SensorArray class for beampatter interval analysis
%
% This is a class of the Beampattern Interval Analysis Toolbox.
% It allows the definition of medium parameters, element 
% poisitions, tapering, calibration errors and coupling, to 
% enable the tolerance analysis of the array's beampattern. 
% The object automatically calculates several properties 
% required for the calculation of beampattern bounds.
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
    SoundSpeed      % Assumed speed of sound (m/s)
    CenterFrequency % Center frequency (monchromatic assumption) [Hz]
    ElCount         % Number of elements in array
    ElPitchRatio    % Ratio of element pitch and wavelength
    ElDiameterRatio % Ratio of element diameter and pitch
    Curvature       % Curvature radius (zero for linear array) [m]
    TaperType       % Type of tapering (Hamming, Kaiser, Chebyshev,custom)
    TaperParam      % Tapering parameter (Kaiser beta, Chebyshev PSLL, custom weights)
    TaperAngles     % Element directivity taper angles (outside this angle the element has zero sensitivity) [rad]
    CouplingCoeff   % Coupling coefficient (between elements)
    SteeringAngle   % Array steering angle in radians [rad]
    GainError       % Element gain error (+/-) as percentage of the nominal gain
    PhaseError      % Element phase error (+/-) in radians [rad]
    OrientError     % Element orientation error (+/-) in radians [rad]
    PosXError       % Element along-array-axis position error (+/-) in meters [m]
    PosYError       % Element across-array-axis position error (+/-) in meters [m]
    GainBias        % Element gain bias as percentage of the nominal gain
    PhaseBias       % Element phase bias (+/-) in radians [rad]
    OrientBias      % Element orientation bias (+/-) in radians [rad]
    PosXBias        % Element along-array-axis position bias (+/-) in meters [m]
    PosYBias        % Element across-array-axis position bias (+/-) in meters [m]
end

properties (Dependent)
    WaveLength      % Assumed wavelength [m]
    WaveNumber      % Assumed wavenumber [rad/m]
    SteeringVector  % Steering vector calculated from the wavenumber and steering angle [rad/m]
    ArrayCoeff      % Array coefficient representing the total effect of tapering, steering and coupling
    ElPitch         % Element pitch calculated as a ratio of the wavelength [m]
    ElDiameter      % Element diameter calculated as a ratio of the element pitch [m]
    ElPosX          % Element along-array-axis position [m]
    ElPosY          % Element accross-array-axis position [m]
    ElOrient        % Element orientations relative to the array-axis [rad]
    TaperWeights    % Tapering weights
    CouplingAmp     % Element amplitude coupling
    CouplingPhase   % Element phase coupling [rad]
    GainInterval    % Element gain interval calculated as the nominal gain * (1 +/- gain error) + gain bias 
    PhaseInterval   % Element phase interval calculated as the nominal phase +/- phase error + phase bias [rad]
    OrientInterval  % Element orientation interval calculated as the nominal orientation +/- orientation error + orientation bias [rad]
    PosXInterval    % Element along-array-axis position interval calculated as the nominal position +/- position error + position bias [m]
    PosYInterval    % Element accross-array-axis position interval calculated as the nominal position +/- position error + position bias [m]
end

methods
    %% Define class constructor
    function obj = SensorArray(options)
        %SENSORARRAY Construct an instance of this class
        %
        % This function generates a sensor array object based on  
        % the optional input argument, each having a default value,
        % therefore if no argument is given, the generated
        % object is a default array with 11 elements half-lambda
        % pitch (assuming a sound speed of 1500 m/s and a central 
        % frequency of 20000 Hz, giving a wavelength of 75 mm),
        % no tapering and 1% gain error and 2 degrees 
        % phase error, and no other errors or coupling.
        % This class is a copyable handle class. This means that 
        % the constructor returns a handle to the object, which
        % can be passed on as a parameter providing multiple 
        % reference to the same object (for example multiple 
        % beampattern object using the same array object). If a
        % property is changed via one handle the property read 
        % via another handle will be also changed. It is however
        % possible to create an independent copy using the copy()
        % function, which will act as a seperate object.
        %__________________________________________________________________________
        % USAGE        
        %   biat.SensorArray(Name,Value)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : object of biat.BeamPattern type
        % _________________________________________________________________________
        % OPTIONS
        %   all optional arguments have the same name as the 
        %   corresponding object property, give the property
        %   name and selected value in the arguments to set 
        %   it different from the default value.
        % _________________________________________________________________________
        % EXAMPLES
        %   array = biat.SensorArray();
        %   array = biat.SensorArray(   'ElCount',5,...
        %                               'ElDiameterRatio',0,...
        %                               'Curvature',0.2,...
        %                               'TaperType','chebwin',...
        %                               'TaperParam',20,...
        %                               'GainError',5/100,...
        %                               'PhaseError',deg2rad(4),...
        %                               'SteeringAngle',deg2rad(5));
        %   sameArray = array;          % (handle copy)
        %   otherArray = copy(array);   % (object copy)
        % _________________________________________________________________________

        arguments
            options.SoundSpeed      (1,1)   {mustBeNumeric} = 1500
            options.CenterFrequency (1,1)   {mustBeNumeric} = 20000
            options.ElCount         (1,1)   {mustBeInteger} = 11
            options.ElPitchRatio    (1,1)   {mustBeNumeric} = 0.5
            options.ElDiameterRatio (1,1)   {mustBeNumeric} = 0
            options.Curvature       (1,1)   {mustBeNumeric} = 0
            options.TaperType       (1,1)   string {mustBeMember(...
                                                    options.TaperType ...
                                                    ,{'none',...
                                                      'hamming',...
                                                      'kaiser',...
                                                      'chebwin',...
                                                      'custom'})} ...
                                                = 'none'
            options.TaperParam      (:,1)   {mustBeNumeric} = 0   
            options.TaperAngles     (1,2)   {mustBeNumeric} = deg2rad([80,100])
            options.CouplingCoeff   (1,1)   {mustBeNumeric} = 0
            options.GainError       (:,1)   {mustBeNumeric} = 0.01
            options.PhaseError      (:,1)   {mustBeNumeric} = deg2rad(2)
            options.OrientError     (:,1)   {mustBeNumeric} = 0
            options.PosXError       (:,1)   {mustBeNumeric} = 0
            options.PosYError       (:,1)   {mustBeNumeric} = 0
            options.GainBias        (:,1)   {mustBeNumeric} = 0
            options.PhaseBias       (:,1)   {mustBeNumeric} = 0
            options.OrientBias      (:,1)   {mustBeNumeric} = 0
            options.PosXBias        (:,1)   {mustBeNumeric} = 0
            options.PosYBias        (:,1)   {mustBeNumeric} = 0
            options.SteeringAngle   (1,1)   {mustBeNumeric} = 0
        end
        obj.SoundSpeed      = options.SoundSpeed;
        obj.CenterFrequency = options.CenterFrequency;
        obj.ElCount         = options.ElCount;
        obj.ElPitchRatio    = options.ElPitchRatio;
        obj.ElDiameterRatio = options.ElDiameterRatio;
        obj.Curvature       = options.Curvature;
        obj.TaperType       = options.TaperType;
        obj.TaperParam      = options.TaperParam;
        obj.TaperAngles     = options.TaperAngles;
        obj.CouplingCoeff   = options.CouplingCoeff;
        obj.GainError       = options.GainError;    
        obj.PhaseError      = options.PhaseError;   
        obj.OrientError     = options.OrientError;  
        obj.PosXError       = options.PosXError;    
        obj.PosYError       = options.PosYError;
        obj.GainBias        = options.GainBias;    
        obj.PhaseBias       = options.PhaseBias;   
        obj.OrientBias      = options.OrientBias;  
        obj.PosXBias        = options.PosXBias;    
        obj.PosYBias        = options.PosYBias;
        obj.SteeringAngle   = options.SteeringAngle;
    end

    %% Define get methods for dependent variables
    
    % Wave-lenght
    function value = get.WaveLength(obj)
        value = obj.SoundSpeed / obj.CenterFrequency;
    end
    
    % Wave-number
    function value = get.WaveNumber(obj)
        value = 2*pi / obj.WaveLength;
    end
    
    %Elelement pitch
    function value = get.ElPitch(obj)   
        value = obj.WaveLength * obj.ElPitchRatio;
    end
    
    % Element diameter
    function value = get.ElDiameter(obj)
        value = obj.ElPitch * obj.ElDiameterRatio;
    end
    
    % Element position (X)
    function value = get.ElPosX(obj)
        if obj.Curvature == 0
            arrayLength = obj.ElPitch * (obj.ElCount - 1);
            value = linspace(-arrayLength/2, arrayLength/2, obj.ElCount)';
        else
            value = obj.Curvature*sin(obj.ElOrient);
            value = value - mean(value);
        end
    end
    
    % Element positition (Y)
    function value = get.ElPosY(obj) 
        if obj.Curvature == 0
            value = zeros(obj.ElCount,1);
        else
            value = obj.Curvature * (cos(obj.ElOrient) - 1);
        end
    end
    
    % Element orientation
    function value = get.ElOrient(obj)  
        if obj.Curvature == 0
            value = zeros(obj.ElCount,1);
        else
            sectorCoverage = 2*pi* (obj.ElPitch * ...
                               (obj.ElCount-1))...
                               /(2*pi*obj.Curvature);
            value = linspace(-0.5,0.5,obj.ElCount).' * sectorCoverage;
        end
    end
    
    % Tapering weights
    function value = get.TaperWeights(obj)
        switch obj.TaperType
            case 'none'
                value = ones(obj.ElCount,1)/obj.ElCount;
            case 'hamming'
                value = hamming(obj.ElCount);
                value = value / sum(value);
            case 'kaiser'
                if length(obj.TaperParam) == 1
                    value = kaiser(obj.ElCount,obj.TaperParam);
                    value = value / sum(value);
                else
                    error('Invalid parameter given for kaiser tapering')
                end
            case 'chebwin'
                if length(obj.TaperParam) == 1
                    value = chebwin(obj.ElCount,obj.TaperParam);
                    value = value / sum(value);
                else
                    error('Invalid parameter given for chebwin tapering')
                end
            case 'custom'
                if length(obj.TaperParam) == obj.ElCount
                    value = obj.TaperParam(:);
                else
                    error('Invalid parameter given for custom tapering')
                end
            otherwise
                error('Invalid taper type given')
        end                
    end
    
    % Coupling amplitude
    function value = get.CouplingAmp(obj)    
        if obj.CouplingCoeff == 0
            value = eye(obj.ElCount);
        else
            value = toeplitz(exp(((1:obj.ElCount)-1) * ...
                                 log(obj.CouplingCoeff)));
        end
    end
    
    % Coupling phase
    function value = get.CouplingPhase(obj)  
        value = zeros(obj.ElCount);
    end
    
    % Element gain error interval
    function value = get.GainInterval(obj)   
        value = ciat.RealInterval(  1-obj.GainError(:) , ...
                                    1+obj.GainError(:) ) .* ...
                boxcar(obj.ElCount) + obj.GainBias(:);
    end
    
    % Element phase error interval
    function value = get.PhaseInterval(obj)  
        value = ciat.RealInterval( -obj.PhaseError(:) , ...
                    obj.PhaseError(:) ) .* ...
                boxcar(obj.ElCount) + obj.PhaseBias(:);
    end
    
    % Element orientation error interval
    function value = get.OrientInterval(obj)
        value = ciat.RealInterval(obj.ElOrient - obj.OrientError(:), ...
                                  obj.ElOrient + obj.OrientError(:) ) .* ...
                                  boxcar(obj.ElCount) + obj.OrientBias(:);
    end
    
    % Element position (X) error interval
    function value = get.PosXInterval(obj)   
        value = ciat.RealInterval(  obj.ElPosX - obj.PosXError(:), ...
                                    obj.ElPosX + obj.PosXError(:) ) .* ...
                                    boxcar(obj.ElCount) + obj.PosXBias(:);
    end
    
    % Element position (Y) error interval
    function value = get.PosYInterval(obj)   
        value = ciat.RealInterval(  obj.ElPosY - obj.PosYError(:), ...
                                    obj.ElPosY + obj.PosYError(:) ) .* ...
                                    boxcar(obj.ElCount) + obj.PosYBias(:);
    end
    
    % Array steering vector
    function value = get.SteeringVector(obj)   
        value = 2*pi / obj.WaveLength * ...
                [sin(obj.SteeringAngle); cos(obj.SteeringAngle)];
    end
    
    % Array coefficient (includes coupling, tapering and steering 
    function value = get.ArrayCoeff(obj)  
        value = getArrayCoeff(obj);
    end
    
    %% DEFINE METHOD SIGNATURES FOR FUNCTIONS IN EXTERNAL FILES
    
    value = getElDirectivity(obj,angle)
    value = getArrayCoeff(obj)  
    plot(obj,optional)  
    end % methods

    
    %% DEFINE STATIC METHODS 
    methods (Static)
        %% Define method signatures
        [x,y] = getArrowLine(r1,r2,L)
    end % static methods    
    
end % class

