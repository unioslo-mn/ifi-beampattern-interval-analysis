% GUI Run graphical user interface for the Beampattern_IA library
%
%  USAGE:   GUI()
% 
% _________________________________________________________________________
%
% Copyright (C) 2023 Haavard Arnestad and Gabor Gereb, BSD-3 (LICENSE.md)
% If you use this software, please cite it as in CITATION.cff
% Project: Beampattern Interval Analysis (doi.org/10.5281/zenodo.6856232)
% Contact: haavaarn@uio.no, gaborge@uio.no (more in README.md)
% _________________________________________________________________________


%% Initialize figures

% Close figures
clear
close all

% GUI figure
gui.figGui = figure();clf
gui.figGui.Position=[0,0,1280, 720];

% Physical array subplot
gui.plt.PA = subplot(2,2,2);

% Beampattern subplot
gui.plt.BP = subplot(2,2,3);
title('Beampattern')
xlabel('Steering angle [rad]');
ylabel('Power [dB]');

% Array factor subplot            
gui.plt.CI = subplot(2,2,4);
title('Complex intervals')
xlabel('Real');
ylabel('Imag');

%% Initialize variables

% Initialize sensor array
initializeControllers();

%% Initialize control panels

% Panel-1: physical parameters
gui.uip1 = uipanel('Units','normalized',...
    'Position',[0.1,0.4875,0.4,0.4125]);
            
% Panel-2: plotting parameters
gui.uip2 = uipanel('Units','normalized',...
               'Position',[0.468,0.10875,0.03,0.34475]);
           
%% Add title and information text to Panel-1
gui.text.Par = uicontrol('Style','text',...
                'String',['Backtracking is not implemented for non-zero ',...
                          'orientation error, and ambiguous for position ',...
                          'and phase errors. ',...
                          'Beampattern with rectangular and circular ',...
                          'methods are not implemented for non-zero ',... 
                          'coupling coefficient.'],...
                'Parent',gui.uip1,...
                'Units','normalized',...
                'Position',[0.03,0.23,0.93,0.26]);


%% Create physical array controllers in Panel-1 (left column)

% Initialize position variables
dY = 0.07;
textPos = [0,0.99,0.16,0.06]; 
ctrlPos = [0.16,0.99,0.24,0.06];
valPos = [0.40,0.99,0.08,0.06];

%   Slider: Number of elements
textPos(2)=textPos(2)-dY; ctrlPos(2)=ctrlPos(2)-dY; valPos(2)=valPos(2)-dY;

gui.text.Nel = uicontrol('Style','text',...
                'String','N_element',...
                'Parent',gui.uip1,...
                'Units','normalized',...
                'Position',textPos);
gui.ctrl.Nel = uicontrol( 'Style', 'slider',...
                'Parent',gui.uip1,...
                'Units','normalized',...
                'Position',ctrlPos,...
                'Min',2,...
                'Max',22,...
                'SliderStep',[0.05,0.1],...
                'Value',gui.array.ElCount,...
                'ToolTip','Number of elements', ...
                'Callback',{@cbNel});
gui.val.Nel = uicontrol('Style','text',...
                'String',formatValueString('cbNel',gui.ctrl.Nel.Value),...
                'Parent',gui.uip1,...
                'Units','normalized',...
                'Position',valPos);
          
%   Slider: Curvature
textPos(2)=textPos(2)-dY; ctrlPos(2)=ctrlPos(2)-dY; valPos(2)=valPos(2)-dY;
gui.text.Rcv = uicontrol('Style','text',...
                'String','Curvature',...
                'Parent',gui.uip1,...
                'Units','normalized',...
                'Position',textPos);
gui.ctrl.Rcv = uicontrol( 'Style', 'slider',...
                'Parent',gui.uip1,...
                'Units','normalized',...
                'Position',ctrlPos,...
                'Min',0,...
                'Max',1,...
                'SliderStep',[0.01,0.1],...
                'Value',gui.array.Curvature,...
                'ToolTip','Number of elements',...
                 'Callback',{@cbRcv});
gui.val.Rcv = uicontrol('Style','text',...
                'String',formatValueString('cbRcv',gui.ctrl.Rcv.Value),...
                'Parent',gui.uip1,...
                'Units','normalized',...
                'Position',valPos);
            
%   Slider: Element distance
textPos(2)=textPos(2)-dY; ctrlPos(2)=ctrlPos(2)-dY; valPos(2)=valPos(2)-dY;
gui.text.Eld = uicontrol('Style','text',...
                'String','El. dist',...
                'Parent',gui.uip1,...
                'Units','normalized',...
                'Position',textPos);
gui.ctrl.Eld = uicontrol( 'Style', 'slider',...
                'Parent',gui.uip1,...
                'Units','normalized',...
                'Position',ctrlPos,...
                'Min',0,...
                'Max',1,...
                'SliderStep',[0.01,0.1],...
                'Value',gui.array.ElPitchRatio,...
                'ToolTip','Number of elements', ...
                'Callback',{@cbEld});
gui.val.Eld = uicontrol('Style','text',...
                'String',formatValueString('cbEld',gui.ctrl.Eld.Value),...
                'Parent',gui.uip1,...
                'Units','normalized',...
                'Position',valPos);
            
%   Slider: Element size
textPos(2)=textPos(2)-dY; ctrlPos(2)=ctrlPos(2)-dY; valPos(2)=valPos(2)-dY;
gui.text.Els = uicontrol('Style','text',...
                'String','El. size%',...
                'Parent',gui.uip1,...
                'Units','normalized',...
                'Position',textPos);
gui.ctrl.Els = uicontrol( 'Style', 'slider',...
                'Parent',gui.uip1,...
                'Units','normalized',...
                'Position',ctrlPos,...
                'Min',0,...
                'Max',1,...
                'SliderStep',[0.01,0.1],...
                'Value',gui.array.ElDiameterRatio,...
                'ToolTip','Number of elements', ...
                'Callback',{@cbEls});
gui.val.Els = uicontrol('Style','text',...
                'String',formatValueString('cbEls',gui.ctrl.Els.Value),...
                'Parent',gui.uip1,...
                'Units','normalized',...
                'Position',valPos);
            
%   Slider: Steering
textPos(2)=textPos(2)-dY; ctrlPos(2)=ctrlPos(2)-dY; valPos(2)=valPos(2)-dY;
gui.text.Ste = uicontrol('Style','text',...
                'String','Steering',...
                'Parent',gui.uip1,...
                'Units','normalized',...
                'Position',textPos);
gui.ctrl.Ste = uicontrol( 'Style', 'slider',...
                'Parent',gui.uip1,...
                'Units','normalized',...
                'Position',ctrlPos,...
                'Min',-1,...
                'Max',1,...
                'SliderStep',[0.01,0.1],...
                'Value',gui.array.SteeringAngle/pi*2,...
                'ToolTip','Number of elements', ...
                'Callback',{@cbSte});
gui.val.Ste = uicontrol('Style','text',...
                'String',formatValueString('cbSte',gui.ctrl.Ste.Value),...
                'Parent',gui.uip1,...
                'Units','normalized',...
                'Position',valPos);
            
%   Slider: Weighting (Kaiser beta selector)
textPos(2)=textPos(2)-dY; ctrlPos(2)=ctrlPos(2)-dY; valPos(2)=valPos(2)-dY;
gui.text.Wgh = uicontrol('Style','text',...
                'String','Weighting',...
                'Parent',gui.uip1,...
                'Units','normalized',...
                'Position',textPos);
gui.ctrl.Wgh = uicontrol( 'Style', 'slider',...
                'Parent',gui.uip1,...
                'Units','normalized',...
                'Position',ctrlPos,...
                'Min',0,...
                'Max',10,...
                'SliderStep',[0.1,0.1],...
                'Value',gui.array.TaperParam,...
                'ToolTip','Weighting (Kaiser beta)', ...
                'Callback',{@cbWgh}); 
gui.val.Wgh = uicontrol('Style','text',...
                'String',formatValueString('cbWgh',gui.ctrl.Wgh.Value),...
                'Parent',gui.uip1,...
                'Units','normalized',...
                'Position',valPos);
            
%% Create array error controllers in Panel-1 (right column)

% Initialize position variables
dY = 0.07;
textPos = [0.5,0.99,0.16,0.06]; 
ctrlPos = [0.66,0.99,0.24,0.06];
valPos = [0.90,0.99,0.08,0.06];

%   Slider: Gain error
textPos(2)=textPos(2)-dY; ctrlPos(2)=ctrlPos(2)-dY; valPos(2)=valPos(2)-dY;
gui.text.EGn = uicontrol('Style','text',...
                'String','Gain error',...
                'Parent',gui.uip1,...
                'Units','normalized',...
                'Position',textPos);
gui.ctrl.EGn = uicontrol( 'Style', 'slider',...
                'Parent',gui.uip1,...
                'Units','normalized',...
                'Position',ctrlPos,...
                'Min',0,...
                'Max',1,...
                'SliderStep',[0.01,0.1],...
                'Value',gui.array.GainError,...
                'ToolTip','Gain error',...
                 'Callback',{@cbEGn});
gui.val.EGn = uicontrol('Style','text',...
                'String',formatValueString('cbEGn',gui.ctrl.EGn.Value),...
                'Parent',gui.uip1,...
                'Units','normalized',...
                'Position',valPos);
            
%   Slider: Phase error
textPos(2)=textPos(2)-dY; ctrlPos(2)=ctrlPos(2)-dY; valPos(2)=valPos(2)-dY;
gui.text.EPh = uicontrol('Style','text',...
                'String','Phase error',...
                'Parent',gui.uip1,...
                'Units','normalized',...
                'Position',textPos);
gui.ctrl.EPh = uicontrol( 'Style', 'slider',...
                'Parent',gui.uip1,...
                'Units','normalized',...
                'Position',ctrlPos,...
                'Min',0,...
                'Max',pi,...
                'SliderStep',[1/180,1/18],...
                'Value',gui.array.PhaseError,...
                'ToolTip','Phase error', ...
                'Callback',{@cbEPh});      
gui.val.EPh = uicontrol('Style','text',...
                'String',formatValueString('cbEPh',gui.ctrl.EPh.Value),...
                'Parent',gui.uip1,...
                'Units','normalized',...
                'Position',valPos);
            
%   Slider: Orientation error
textPos(2)=textPos(2)-dY; ctrlPos(2)=ctrlPos(2)-dY; valPos(2)=valPos(2)-dY;
gui.text.EPs = uicontrol('Style','text',...
                'String','Orient error',...
                'Parent',gui.uip1,...
                'Units','normalized',...
                'Position',textPos);
gui.ctrl.EPs = uicontrol( 'Style', 'slider',...
                'Parent',gui.uip1,...
                'Units','normalized',...
                'Position',ctrlPos,...
                'Min',0,...
                'Max',pi,...
                'SliderStep',[1/180,1/18],...
                'Value',gui.array.OrientError,...
                'ToolTip','Orientation error', ...
                'Callback',{@cbEPs});    
gui.val.EPs = uicontrol('Style','text',...
                'String',formatValueString('cbEPs',gui.ctrl.EPs.Value),...
                'Parent',gui.uip1,...
                'Units','normalized',...
                'Position',valPos);
            
%   Slider: Position error (X)
textPos(2)=textPos(2)-dY; ctrlPos(2)=ctrlPos(2)-dY; valPos(2)=valPos(2)-dY;
gui.text.EPx = uicontrol('Style','text','String',...
                'Pos(x) error',...
                'Parent',gui.uip1,...
                'Units','normalized',...
                'Position',textPos);
gui.ctrl.EPx = uicontrol( 'Style', 'slider',...
                'Parent',gui.uip1,...
                'Units','normalized',...
                'Position',ctrlPos,...
                'Min',0,...
                'Max',1,...
                'SliderStep',[0.01,0.1],...
                'Value',gui.array.PosXError,...
                'ToolTip','Position error (X)', ...
                'Callback',{@cbEPx});  
gui.val.EPx = uicontrol('Style','text',...
                'String',formatValueString('cbEPx',gui.ctrl.EPx.Value),...
                'Parent',gui.uip1,...
                'Units','normalized',...
                'Position',valPos);
            
%   Slider: Position error (Y)
textPos(2)=textPos(2)-dY; ctrlPos(2)=ctrlPos(2)-dY; valPos(2)=valPos(2)-dY;
gui.text.EPy = uicontrol('Style','text',...
                'String','Pos(y) error',...
                'Parent',gui.uip1,...
                'Units','normalized',...
                'Position',textPos);
gui.ctrl.EPy = uicontrol( 'Style', 'slider',...
                'Parent',gui.uip1,...
                'Units','normalized',...
                'Position',ctrlPos,...
                'Min',0,...
                'Max',1,...
                'SliderStep',[0.01,0.1],...
                'Value',gui.array.PosYError,...
                'ToolTip','Position error (Y)', ...
                'Callback',{@cbEPy});  
gui.val.EPy = uicontrol('Style','text',...
                'String',formatValueString('cbEPy',gui.ctrl.EPy.Value),...
                'Parent',gui.uip1,...
                'Units','normalized',...
                'Position',valPos);
            
%   Slider: Coupling coefficient
textPos(2)=textPos(2)-dY; ctrlPos(2)=ctrlPos(2)-dY; valPos(2)=valPos(2)-dY;
gui.text.Cpl = uicontrol('Style','text',...
                'String','Coupling coeff',...
                'Parent',gui.uip1,...
                'Units','normalized',...
                'Position',textPos);
gui.ctrl.Cpl = uicontrol( 'Style', 'slider',...
                'Parent',gui.uip1,...
                'Units','normalized',...
                'Position',ctrlPos,...
                'Min',0,...
                'Max',1,...
                'SliderStep',[0.01,0.1],...
                'Value',gui.array.CouplingCoeff,...
                'ToolTip','Coupling coefficient', ...
                'Callback',{@cbCpl});  
gui.val.Cpl = uicontrol('Style','text',...
                'String',formatValueString('cbCpl',gui.ctrl.Cpl.Value),...
                'Parent',gui.uip1,...
                'Units','normalized',...
                'Position',valPos);
            
            
%% Create miscellaneous controllers in Panel-1 (bottom)

%   Button: Update beampattern
gui.btn.Update = uicontrol( 'Style', 'togglebutton',...
                'Units','normalized',...
                'Position',[0.78,0.082,0.15,0.14],...
                'Parent',gui.uip1,...
                'String','Update',...
                'ToolTip','Update beampatter', ...
                'Callback',{@cbUpdate});
            
%   Button: Set default values
gui.btn.Default = uicontrol( 'Style', 'pushbutton',...
                'Units','normalized',...
                'Position',[0.6,0.08,0.15,0.14],...
                'Parent',gui.uip1,...
                'String','Default',...
                'ToolTip','Set default values', ...
                'Callback',{@cbDefault});

% Theta index selector
gui.ctrl.Thi = uicontrol( 'Style', 'slider',...
                'Parent',gui.uip1,...
                'Units','normalized',...
                'Position',[0.03,0.01,0.9,0.06],...
                'Min',-1,'Max',1,...
                'Value',sin(((gui.beampattern(1).BeamIndex-1) / ...
                             (gui.beampattern(1).BeamCount-1)*2-1)*pi/2),...
                'SliderStep',[1/gui.beampattern(1).BeamCount,...
                             10/gui.beampattern(1).BeamCount],...
                'ToolTip','Number of elements', ...
                'Callback',{@cbThi});
gui.val.Thi = uicontrol('Style','text',...
                'String',num2str(gui.ctrl.Thi.Value),...
                'Parent',gui.uip1,...
                'Units','normalized',...
                'Position',[0.925,0.01,0.05,0.06]);
            
%% Create plot controllers in Panel-2

% Initialize position variables
dY = 0.055;
tbPos = [0.34,0.98,0.36,0.055];

% Add checkboxes
tbPos(2)=tbPos(2) - dY;
gui.tb.N = uicontrol('Style','checkbox',...
                'Units','normalized',...
                'Position',tbPos,...
                'Parent',gui.uip2,...
                'Value',gui.plt.sel(1),...
                'Callback',@cbTbN);
tbPos(2)=tbPos(2) - dY;
gui.tb.R = uicontrol('Style','checkbox',...
                'Units','normalized',...
                'Position',tbPos,...
                'Parent',gui.uip2,...
                'Value',gui.plt.sel(2),...
                'Callback',@cbTbR);
tbPos(2)=tbPos(2) - dY;
gui.tb.C = uicontrol('Style','checkbox',...
                'Units','normalized',...
                'Position',tbPos,...
                'Parent',gui.uip2,...
                'Value',gui.plt.sel(3),...
                'Callback',@cbTbC);
tbPos(2)=tbPos(2) - dY;
gui.tb.G = uicontrol('Style','checkbox',...
                'Units','normalized',...
                'Position',tbPos,...
                'Parent',gui.uip2,...
                'Value',gui.plt.sel(4),...
                'Callback',@cbTbG);
tbPos(2)=tbPos(2) - dY;
gui.tb.X = uicontrol('Style','checkbox',...
                'Units','normalized',...
                'Position',tbPos,...
                'Parent',gui.uip2,...
                'Value',gui.plt.sel(5),...
                'Callback',@cbTbX);
gui.tb = [gui.tb.N;gui.tb.R;gui.tb.C;gui.tb.G;gui.tb.X];
            
% Add power dynamic range (power floor) selector slider            
% Theta index selector
gui.ctrl.DRn = uicontrol( 'Style', 'slider',...
                'Parent',gui.uip2,...
                'Units','normalized',...
                'Position',[0.3,0.007,0.49,0.60],...
                'Min',-110,...
                'Max',-10,...
                'Value',-60,...
                'SliderStep',[0.01,0.1],...
                'ToolTip','Number of elements', ...
                'Callback',{@cbDRn});
            
%% Initialize plots

plotArrayInfo();
plotBeamPattern();
plotComplexIntervals();
            
            
%% Callback functions of physical array controllers (Panel-1 left)

%   Slider: Number of elements
function a = cbNel(h,~)
    updateArray(h);
end

%   Slider: Curvature
function a = cbRcv(h,~)
    updateArray(h);
end

%   Slider: Element distance
function a = cbEld(h,~)
    updateArray(h);
end

%   Slider: Element size
function a = cbEls(h,~)
    updateArray(h);
end

%   Slider: Element size
function a = cbSte(h,~)
    updateArray(h);
end

%   Slider: Weighting
function a = cbWgh(h,~)
    updateArray(h);
end

%% Callback functions of array error controllers (Panel-1 right)

%   Slider: Gain error
function a = cbEGn(h,~)
    updateArray(h);
end

%   Slider: Phase error
function a = cbEPh(h,~)
    updateArray(h);
end

%   Slider: Orientation error
function a = cbEPs(h,~)
    updateArray(h);
end

%   Slider: Position error (X)
function a = cbEPx(h,~)
    updateArray(h);
end

%   Slider: Position error (Y)
function a = cbEPy(h,~)
    updateArray(h);
end

%   Slider: Coupling coefficient
function a = cbCpl(h,~)
    updateArray(h);
end

%% Callback functions of miscellaneous controllers (Panel-1 bottom)

%   Button: Update beampattern
function a = cbUpdate(h,~)
    updatePlots(h);
end

%   Button: Set default parameters
function a = cbDefault(h,~)
    initializeControllers();
    resetControllers(h);
    plotArrayInfo();
end

%   Slider: Theta index
function a = cbThi(h,~)
    updateBeamPattern(h);
end

%   Slider: Dynamic range
function a = cbDRn(h,~)
    updateBeamPattern(h);
end

%% Callback functions of plot controllers (Panel-2)

function a = cbTbN(h,~)
    updateBeamPattern(h);
end
function a = cbTbR(h,~)
    updateBeamPattern(h);
end
function a = cbTbC(h,~)
    updateBeamPattern(h);
end
function a = cbTbG(h,~)
    updateBeamPattern(h);
end
function a = cbTbX(h,~)
    updateBeamPattern(h);
end
