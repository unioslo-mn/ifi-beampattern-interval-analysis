function [h,c] = dcolor(varargin)
% Plot a complex-valued matrix or function using domain coloring
%
% SYNTAX
%
% DCOLOR(W)
% DCOLOR(X,Y,W)
% DCOLOR(@FZ)
% DCOLOR(D,@FZ)
% DCOLOR(X,Y,@FZ)
% DCOLOR(...,'cfun',CFUN)
% DCOLOR(...,'cfun',@CFUN)
% DCOLOR('cfun',CFUN)
% DCOLOR('cfun',@CFUN)
% DCOLOR(...,'grid')
% DCOLOR(AH,...)
% H = DCOLOR(...)
% [H,C] = DCOLOR(...)
%
% DESCRIPTION 
%
% DCOLOR(W) 
%
% Produce a plot of the complex matrix W.  DCOLOR uses one of several
% coloring functions to display both the angle and modulus of the data on
% the same plot.  This technique is known as domain coloring, and is useful
% to visualize the behavior and characteristics of complex-valued
% functions.
%
% DCOLOR(X,Y,W)
%
% Make a plot using the grid defined by X and Y.  The inputs X and Y can be
% vectors or matrices.  Their sizes must correspond to the size of W.
%
% DCOLOR(@FZ)
% DCOLOR(D,@FZ)
% DCOLOR(X,Y,@FZ)
%
% Calling DCOLOR with a function handle input, @FZ, will evaluate the
% specified function and plot the results.  The function specified by @FZ
% must accept a single complex input matrix, and return a complex matrix as
% output. By default, the function will be evaluated over the interval
% [-1-i,1+i].  If the argument D is supplied in the form of a 1x2 complex
% vector, @FZ will be evaluated over the domain specified by D.  If grid
% points X and Y are specified, either as vectors (which will be expanded
% using MESHGRID(X,Y)) or as MxN arrays, @FZ will be evaluated at the
% points X+iY.
%
% DCOLOR(...,'cfun',CFUN)
% DCOLOR(...,'cfun',@CFUN)
%
% Specifying the property/value pair 'cfun',CFUN, where CFUN is a number
% corresponding to one of the built-in coloring functions, will use the
% specified built-in coloring function to generate the plot.  
%
% The built-in coloring functions are numbered 1-7 and provide alternative
% representations of the complex data.  If CFUN is not specified or out of
% range, coloring function 1 is used.  See the Remarks section for a
% description of the available coloring functions.
% 
% If a function handle @CFUN is supplied, the specified function is
% evaluated and the result is used to form the image.  The coloring
% function must accept a single MxN complex array and return an MxNx3 array
% of RGB triplets (normalized to [0,1]).  
%
% DCOLOR('cfun',CFUN)
% DCOLOR('cfun',@CFUN)
%
% Calling DCOLOR with the 'cfun',CFUN argument pair but no W data will
% produce a plot of the coloring function evaluated over [-1-i,1+i]. This
% can be useful for comparing and understanding the coloring functions.
%
% DCOLOR(...,'grid')
%
% Supplying the optional argument, 'grid', will superimpose a grid in the
% w-domain by creating contours of the real and imaginary parts of W.  The
% real contours are shown with solid lines and the imaginary contours with 
% dotted lines.
% 
% DCOLOR(AH,...)
%
% Make the plot in the axes referenced by the handle AH instead of gca().
%
% H = DCOLOR(...)
%
% Return a handle to the image or surface object created by DCOLOR.
%
% [H,C] = DCOLOR(...)
%
% In addition to the graphics handle H, return the image color data in C.
%
% REMARKS 
% 
% If X and Y are matrices or non-regularly spaced vectors, DCOLOR will plot
% the data as a surface object.  Otherwise, the data will be plotted as an
% image object.  
%
% DCOLOR attempts to calculate an optimal number of contour levels when
% called with the 'grid' option.  However, the results may not be
% satisfactory for some inputs.  If the resulting grid is too sparse or too
% dense, create contours separately as follows:
%   1) Call DCOLOR without the 'grid' option
%   2) Set the axes hold property to 'on'
%   3) Call CONTOUR(X,Y,REAL(W)) and CONTOUR(X,Y,IMAG(W)).  If desired,
%   supply a list of contour levels to the CONTOUR function.
%
% Description of built-in coloring functions:
%
% The built-in coloring functions take an MxN complex array as input, and
% output an MxNx3 array of RGB triples.  The algorithms for calculating the
% output start by computing the modulus and angle of the complex data.
% These values are then mapped into the HSV color space.  The final step is
% to transform the HSV color data to the RGB color space. The available
% coloring functions highlight different aspects of the underlying complex
% data.  A brief description of each algorithm follows:
%
% (1) Hue is determined by angle(W), with red corresponding to 0 angle. The
% hue continues through yellow, green, blue and violet before returning to
% red at an angle of 2*pi. The saturation and value (brightness) are
% determined by a logarithmic scaling of the modulus, with small values
% tending to black and large values fading to white. Thus, zeros will
% appear dark with the color changing in a counterclockwise direction,
% while poles will appear white with clockwise rotation of the color bands.
%
% (2) Similar to (1) except the log-scaled modulus is compressed outside of
% the 2nd and 98th percentiles.  This option may give better results than
% (1) for some functions with a large dynamic range.
%
% (3) Hue is determined by angle using the same algorithm as in (1).
% Saturation and value exhibit dark-to-light bands where the modulus
% crosses a power of 10.
%
% (4),(5) Hue is determined by a normalized logarithmic transformation of
% the modulus.  Pure red represents min(abs(w)), while violet represents
% max(abs(w)).  Saturation and value are determined by the angle, which is
% always limited to a 2*pi interval and wraps around at a particular limit.
% As angle increases towards the limit, the color fades to white, while
% approaching the limit with decreasing angle fades to black. Function (4)
% wraps the angle at -pi, while (5) wraps at 2*pi.
%
% (6) Hue is determined by angle using the same algorithm as in (1).  The
% saturation and value are fixed at 1, so this coloring function displays
% only the angle of the data.
%
% (7) Hue is determined by modulus using the algorithm of (4) and (5).
% Saturation and brightness are fixed at 1.  Thus, this coloring function
% displays only the modulus of the data.
%
% EXAMPLES
%
% % Example 1: A Rational Complex Function
% % Plot a rational function with three zeros and one third-order pole
% x = linspace(-30,30,1001);
% y = x;
% [X,Y] = meshgrid(x,y);
% Z = X+1i*Y;
% W = 10*(Z-(12-5i)).*(Z-(14-4.6i)).*(Z+(7+14i)).*(1./(Z-15i)).^3;
% figure
% dcolor(x,y,W,'CFun',1,'grid')
% axis square
% figure
% dcolor(x,y,W,'cfun',3)
% axis square
% figure
% dcolor(x,y,W,'cfun',4)
% axis square
%
% % Example 2: Evaluate a Complex-valued Function Over a Specified Domain
% % Plot w = f(z) = exp(1/z) over D: [(-1-i),(1+i)]/sqrt(2)
% D = 1/sqrt(2)*[(-1 - 1*i),(1 + 1*i)];
% fz = @(z)(exp(1./z));
% figure
% dcolor(D,fz)
% axis square
%
% % Example 3: A Complex Function in Polar Coordinates
% [th,r] = meshgrid((0:.5:360)*pi/180,0:.005:2);
% [X,Y] = pol2cart(th,r);
% Z = X+1i*Y;
% W =.5*(Z.^8-1).^(1/4).*(1./Z.^4).*(Z-1+1i).*(Z+1-1i).*(Z-1-1i).*(Z+1+1i);
% figure
% dcolor(X,Y,W)
% axis square
%
% REFERENCES
%
% See the following for more information about visualizing complex-valued
% functions using domain coloring, including many example images:
%
% http://www.maa.org/pubs/amm_complements/complex.html
% http://www.math.liu.se/~halun/complex/domain_coloring-unicode.html
% http://www.mathworks.com/company/newsletters/news_notes/clevescorner/summer98.cleve.html 
% http://commons.wikimedia.org/wiki/User:Jan_Homann/Mathematics 
% http://www1.american.edu/academic.depts/cas/mathstat/People/lcrone/ComplexPlot.html
% http://en.wikipedia.org/wiki/Domain_coloring
%

% dcolor.m 
% Copyright (c) 2010-2011 by John Barber
% Distribution and use of this software is governed by the terms in the 
% accompanying license file.

% Release History:
% v 1.0 : 2010-Oct-14
%       - Initial release
% v 2.0 : 2010-Nov-30
%       - Added ability to supply function handle to calculate W
%       - Added new coloring functions
%       - Added ability to supply user-defined coloring function
%       - Added grid option
%       - Changed default plot type from surface to image
%       - Reduced memory usage
%       - Renamed from cplxpcolor.m to dcolor.m
% v 2.1 : 2011-Feb-15
%       - Added image color data output

%% Parse inputs

% Default values
cFun = 1;
gridFlag = 0;
minX = -1;
maxX = 1;
minY = -1;
maxY = 1;
nPts = 512;
surfaceFlag = 0;

% If the first input is an axes handle, use it
if (isscalar(varargin{1}) && ...
        ishghandle(varargin{1}) && ...
        strcmp(get(varargin{1},'Type'),'axes'))
    hAx = newplot(varargin{1});
    varargin(1) = [];
else
    hAx = newplot;
end

% Look for the strings 'CFun' and 'grid' in the inputs 
charIdx = cellfun(@ischar,varargin);
for k = find(charIdx)
    switch lower(varargin{k})
        case 'cfun'
            cFun = varargin{k+1};
        case 'grid'
            gridFlag = 1;
    end
end
varargin(find(charIdx,1):end) = [];

% Get the x, y and w data
if isempty(varargin)
    % Plot a color function over (-1 -1i) to (1 + 1i)
    x = linspace(-1,1,nPts);
    y = x;
    [X Y] = meshgrid(x);
    w = X + 1i*Y;
elseif length(varargin) == 1 
    if isnumeric(varargin{1}) && all(size(varargin{1}) == [1 2])
        % Plot a color function over the specified domain
        minX = real(varargin{1}(1));
        maxX = real(varargin{1}(2));
        minY = imag(varargin{1}(1));
        maxY = imag(varargin{1}(2));
        x = linspace(minX,maxX,nPts);
        y = linspace(minY,maxY,nPts);
        [X,Y] = meshgrid(x,y);
        w = X + 1i*Y;
    elseif isnumeric(varargin{1})
        % Plot w with integer-spaced domain of size(w)
        w = varargin{1};
        minX = 1;
        maxX = size(w,2);
        minY = 1;
        maxY = size(w,1);
        x = minX:maxX;
        y = minY:maxY;
    else
        % Plot w = @fz(z) over (-1 -1i) to (1 + 1i)
        x = linspace(-1,1,nPts);
        y = x;
        [X,Y] = meshgrid(x);
        w = feval(varargin{1},(X+1i*Y));
    end
elseif length(varargin) == 2
    % Plot w = @fz(z) over the specified domain
    minX = real(varargin{1}(1));
    maxX = real(varargin{1}(2));
    minY = imag(varargin{1}(1));
    maxY = imag(varargin{1}(2));
    x = linspace(minX,maxX,nPts);
    y = linspace(minY,maxY,nPts);
    [X,Y] = meshgrid(x,y);
    w = feval(varargin{2},(X+1i*Y));
elseif length(varargin) == 3
    % Domain is [min(x)+i*min(y), max(x)+i*max(y)], with grid points at x
    % and y.  x and y can be vector or matrix (meshgrid-style).
    if min(size(varargin{1})) == 1;
        x = varargin{1};
        y = varargin{2};
        minX = x(1);
        maxX = x(end);
        minY = y(1);
        maxY = y(end);
        % Check for non-uniform grid
        if max(diff(x))-min(diff(x)) > ...
                4*max(eps(minX),eps(maxX)) || ...
           max(diff(y))-min(diff(y)) > ...
                4*max(eps(minY),eps(maxY))
            % Non-uniform grid, use surface instead of image
            surfaceFlag = 1;
        end
        [X,Y] = meshgrid(x,y);
    else
        % If a matrix of gridpoints was given, use surface instead of image
        X = varargin{1};
        Y = varargin{2};
        surfaceFlag = 1;
    end

    if isnumeric(varargin{3})
        % w data was supplied as input
        w = varargin{3};
    else
        % Calculate w=@fz(z) using @fz supplied as input
        w = feval(varargin{3},(X+1i*Y));
    end
else
    msgid = [mfilename ':InvalidInputs'];
    msgtext = ['Invalid input.  Type ''help ' mfilename ...
               ''' for correct syntax'];
    error(msgid,msgtext)
end

if surfaceFlag == 0
    clear X Y
end

%% Compute the color data and make the plot

% Create a grid if necessary
if gridFlag == 1
    % Compute contour levels
    levels = DColorGrid(w);
    % Find and plot contours for real(w)
    if surfaceFlag == 1
        C = contours(X,Y,real(w),levels{1});
    else
        C = contourc(x,y,real(w),levels{1});
    end
    hCReal = hggroup('Parent',hAx);
    pos = 1;
    while pos < size(C,2)
        linepts = C(2,pos);
        line('XData',C(1,(pos+1):(pos+linepts)),...
            'YData',C(2,(pos+1):(pos+linepts)),...
            'Parent',hCReal,...
            'LineStyle','-',...
            'Color',[0 0 0],...
            'LineWidth',0.5)
        pos = pos+linepts+1;
    end
    % Find and plot contours for imag(w)
    if surfaceFlag == 1
        C = contours(X,Y,imag(w),levels{2});
    else
        C = contourc(x,y,imag(w),levels{2});
    end
    hCImag = hggroup('Parent',hAx);
    pos = 1;
    while pos < size(C,2)
        linepts = C(2,pos);
        line('XData',C(1,(pos+1):(pos+linepts)),...
            'YData',C(2,(pos+1):(pos+linepts)),...
            'Parent',hCImag,...
            'LineStyle',':',...
            'Color',[0 0 0],...
            'LineWidth',0.5)
        pos = pos+linepts+1;
    end
    % Save handles so we can make sure grid is displayed on top
    hC = [hCImag; hCReal];
end

% Calculate color data using a coloring function
if isa(cFun,'function_handle')
    % Call the user-supplied function handle
    CData = feval(cFun,w);
else
    % Use one of the built-in coloring functions
    CData = DColorColorFunctions(w,cFun);
end

clear w

if surfaceFlag == 0
    YDir = get(hAx,'YDir');
    hPlot = image('XData',[minX maxX],...
                  'YData',[minY maxY],...
                  'CData',CData,...
                  'Parent',hAx);
    set(hAx,'YDir',YDir)
else
    hPlot = surface(X,Y,zeros(size(X)),...
                    'CData',CData,...
                    'EdgeColor','none',...
                    'FaceColor','flat',...
                    'Parent',hAx);
end

if gridFlag == 1
    % Make sure grid is displayed on top of image
    hh = get(hAx,'Children');
    idxC = (hh==hC(2));
    idxI = (hh==hPlot);
    hh(idxC) = hPlot;
    hh(idxI) = hC(2);
    set(hAx,'Children',hh)
end

if strcmp(get(hAx,'XLimMode'),'auto') && strcmp(get(hAx,'YLimMode'),'auto')
    axis(hAx,'tight')
end
    
% If asked for, return the handle to the image/surface
if nargout ~= 0
    h = hPlot;
end

% If asked for, return the image data
if nargout == 2
    c = CData;
end

end % End of function dcolor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = iplog10(x)
% In-place implementation of log10 for reduced memory use
x = log2(x);
x = x/6.64385618977472436 + x/6.64385618977472525;
end % End of function iplog10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  CData = DColorColorFunctions(w,cFunction)
%Domain coloring functions for plotting complex data
%
% Transforms the complex data in the MxN matrix W into an MxNx3 array of
% RGB values using the algorithm selected by cFunction.  Subfunction used
% by dcolor.m.
%

%% Calculate the color data matrix in HSV space

% Get page size for ugly indexing operations
ps = numel(w); 

if cFunction == 1
    CData = cat(3,angle(w),abs(w),abs(w));
    % Shift H so that a 0 degree angle has an HSV hue of 0 (red)
    CData(CData(:,:,1)<0) = CData(CData(:,:,1)<0)+2*pi;
    % Normalize H to [0,1]
    CData(:,:,1) = CData(:,:,1)/(2*pi);
    % Compute saturation and brightness from modulus
    CData(:,:,2) = 1./(1+0.33*iplog10(CData(:,:,2)./10+1));
    CData(:,:,3) = 1 - 1./(1.1+25*iplog10(2*CData(:,:,3)+1));
    
elseif cFunction == 2
    CData = cat(3,angle(w),iplog10(abs(w)),zeros(size(w)));
    % Shift H so that a 0 degree angle has an HSV hue of 0 (red)
    CData(CData(:,:,1)<0) = CData(CData(:,:,1)<0)+2*pi;
    % Normalize H to [0,1]
    CData(:,:,1) = CData(:,:,1)/(2*pi);
    % Compute saturation and brightness from normalized modulus. Start with
    % log10(modulus).  Compute 2 and 98 percentiles, then compress data
    % outside of this range and normalize to [0,1].  This is then fed into
    % the S and V calculations.
    CData(ps+find(isinf(CData(:,:,2)))) = NaN;
    CData(:,:,3) = reshape(sort(CData((ps+1):(2*ps))),size(CData(:,:,2)));
    lastdat = 2*ps+find(isnan(CData(:,:,3)),1)-1;
    lastdat(isempty(lastdat)) = 3*ps;
    q = [0 100*(0.5:(lastdat-2*ps+1-0.5))./(lastdat-2*ps+1) 100]';
    x = [CData(2*ps+1); CData((2*ps+1):lastdat)'; CData(lastdat)];
    % 2 calls to interp1q with scalar xi are faster than 1 call with vector
    y(2) = interp1q(q,x,98);
    y(1) = interp1q(q,x,2);
    idx = ps+find(CData(:,:,2)<y(1));
    CData(idx) = y(1)+2/pi*(atan(CData(idx)-max(CData(idx))-.5)-atan(-.5));
    idx = ps+find(CData(:,:,2)>y(2));
    CData(idx) = y(2)+2/pi*(atan(CData(idx)-min(CData(idx))+.5)-atan(.5));
    CData(:,:,2) = (CData(:,:,2)-min(min(CData(:,:,2))))/...
                   (max(max(CData(:,:,2)))-min(min(CData(:,:,2))));
    % Calculate S and V using atan(), then normalize to [0.1,1)
    CData(:,:,3) = CData(:,:,2);
    CData(:,:,2) = 1 - atan(CData(:,:,2).^3 - 1);
    CData(:,:,2) = 0.89999*(CData(:,:,2)-min(min(CData(:,:,2))))/...
                   (max(max(CData(:,:,2)))-min(min(CData(:,:,2)))) + 0.1;
    CData(:,:,3) = atan((CData(:,:,3)+1).^3);
    CData(:,:,3) = 0.79999*(CData(:,:,3)-min(min(CData(:,:,3))))/...
                   (max(max(CData(:,:,3)))-min(min(CData(:,:,3)))) + 0.2;
    
elseif cFunction == 3
    CData = cat(3,angle(w),abs(w),abs(w));
    % Shift H so that a 0 degree angle has an HSV hue of 0 (red)
    CData(CData(:,:,1)<0) = CData(CData(:,:,1)<0)+2*pi;
    % Normalize H to [0,1]
    CData(:,:,1) = CData(:,:,1)/(2*pi);
    % Make dark:normal:light bands for every decade of modulus
    CData(ps+find(isinf(CData(:,:,2)))) = NaN;
    CData(:,:,2) = iplog10(CData(:,:,2)) - floor(iplog10(CData(:,:,2)));
    CData(:,:,2) = 10.^(CData(:,:,2));
    CData(:,:,3) = CData(:,:,2);
    CData(:,:,2) = atan(CData(:,:,2).^5);
    CData(:,:,2) = 0.59999*(CData(:,:,2)-min(min(CData(:,:,2))))/...
                   (max(max(CData(:,:,2)))-min(min(CData(:,:,2)))) + 0.4;
    CData(:,:,3) = 1-atan(((CData(:,:,3)/10).^10-1));
    CData(:,:,3) = 0.59999*(CData(:,:,3)-min(min(CData(:,:,3))))/...
                   (max(max(CData(:,:,3)))-min(min(CData(:,:,3)))) + 0.4;
    
elseif cFunction == 4
    CData = cat(3,abs(w),angle(w),angle(w));
    % Compute hue from normalized modulus
    CData(isinf(CData(:,:,1))) = NaN;
    CData(:,:,1) = 0.84999*(CData(:,:,1)-min(min(CData(:,:,1))))/...
                   (max(max(CData(:,:,1)))-min(min(CData(:,:,1))))+0.00001;
    % Use angle: [-pi,pi)
    CData(:,:,2) = (CData(:,:,2) + pi)/(2*pi);
    CData(:,:,3) = CData(:,:,2);
    CData(:,:,2) = 1-1./(1.01+19*iplog10(-CData(:,:,2)+2));
    CData(:,:,3) = 1-1./(1.01+19*iplog10(CData(:,:,3)+1));
    
elseif cFunction == 5
    CData = cat(3,abs(w),angle(w),angle(w));
    % Compute hue from normalized modulus
    CData(isinf(CData(:,:,1))) = NaN;
    CData(:,:,1) = 0.84999*(CData(:,:,1)-min(min(CData(:,:,1))))/...
                   (max(max(CData(:,:,1)))-min(min(CData(:,:,1))))+0.00001;
    % Use angle: [0,2*pi) 
    CData(ps+find(CData(:,:,2)<0)) = CData(ps+find(CData(:,:,2)<0))+2*pi;
    CData(:,:,2) = CData(:,:,2)/(2*pi);
    CData(:,:,3) = CData(:,:,2);
    CData(:,:,2) = 1-1./(1.0+39*iplog10(-CData(:,:,2)+2));
    CData(:,:,3) = 1-1./(1.0+39*iplog10(CData(:,:,3)+1));
               
elseif cFunction == 6
    CData = cat(3,angle(w),ones(size(w)),ones(size(w)));
    % Saturation and brightness are fixed at 1
    % Shift H so that a 0 degree angle has an HSV hue of 0 (red)
    CData(CData(:,:,1)<0) = CData(CData(:,:,1)<0)+2*pi;
    % Normalize H to [0,1]
    CData(:,:,1) = CData(:,:,1)/(2*pi);
               
elseif cFunction == 7             
    CData = cat(3,abs(w),ones(size(w)),ones(size(w)));
    % Saturation and brightness are fixed at 1
    % Compute hue from normalized modulus
    CData(:,:,1) = iplog10(CData(:,:,1)+1);
    CData(isinf(CData(:,:,1))) = NaN;
    CData(:,:,1) = 0.84999*(CData(:,:,1)-min(min(CData(:,:,1))))/...
                   (max(max(CData(:,:,1)))-min(min(CData(:,:,1))))+0.00001;

else
    % For invalid cFunction, call ourselves recursively with cFunction = 1
    CData = DColorColorFunctions(w,1);
    return
end

%%
% Now convert the HSV values to RGB.  Do it this way (with ugly code)
% because it is much more memory efficient than calling hsv2rgb().
%
% The transform is:
% Inputs:
% h: [0,6)
% s and v: [0,1]
% Intermediate values:
% f =  h-floor(h),   floor(h) odd
%      1-h+floor(h), floor(h) even
% m = v*(1-s)
% n = v*(1-s*f)
% Outputs (all [0,1]):
% (r,g,b) = (v,n,m),   0<=h<1
%           (n,v,m),   1<=h<2
%           (m,v,n),   2<=h<3
%           (m,n,v),   3<=h<4
%           (n,m,v),   4<=h<5
%           (v,m,n),   5<=h<6

% Renormalize hue from [0,1] to [0,6)
CData(:,:,1) = CData(:,:,1)*6;
CData(CData(:,:,1)==6) = 0;

% Group hue values into 6 bins
hp = uint8(floor(CData(:,:,1))); 

% Compute intermediate values: f, n+v, m-v
CData(:,:,1) = CData(:,:,1)-double(hp); % f=h-floor(h)
CData(:,:,1) = CData(:,:,3).*(1 - CData(:,:,2).*(1 - double(mod(hp,2))+...
           CData(:,:,1).*(-1).^double(1+mod(hp,2)))) + CData(:,:,3); %n+v  
CData(:,:,2) = -CData(:,:,3).*CData(:,:,2); %m-v

% CData now has (n+v,m-v,v).  The rest of the code is to shift these values
% around to get permutations of (v,n,m), depending on the hue.

% 0 <= h < 1 : (n+v,m-v,v)->(v,n,m)
idx = find(hp==0);
CData(idx+2*ps) =             +CData(idx+ps) +CData(idx+2*ps); %(n+v,m-v,m)
CData(idx+ps)   = +CData(idx) +CData(idx+ps) -CData(idx+2*ps); %(n+v,n,m)
CData(idx)      = +CData(idx) -CData(idx+ps)                 ; %(v,n,m)

% 1 <= h < 2 :(n+v,m-v,v)->(n,v,m)
idx = find(hp==1);
CData(idx)      = +CData(idx)                -CData(idx+2*ps); %(n,m-v,v)
CData(idx+2*ps) =             +CData(idx+ps) +CData(idx+2*ps); %(n,m-v,m)
CData(idx+ps)   =             -CData(idx+ps) +CData(idx+2*ps); %(n,v,m)

% 2 <= h < 3 :(n+v,m-v,v))->(m,v,n)
idx = find(hp==2);
CData(idx+2*ps) = +CData(idx)                -CData(idx+2*ps); %(n+v,m-v,n)
CData(idx)      = +CData(idx) +CData(idx+ps) -CData(idx+2*ps); %(m,m-v,n)
CData(idx+ps)   = +CData(idx) -CData(idx+ps)                 ; %(m,v,n)

% 3 <= h < 4 :(n+v,m-v,v)->(m,n,v)
idx = find(hp==3);
CData(idx+ps)   = +CData(idx) +CData(idx+ps)                 ; %(n+v,m+n,v)
CData(idx)      = -CData(idx) +CData(idx+ps) +CData(idx+2*ps); %(m,m+n,v)
CData(idx+ps)   = -CData(idx) +CData(idx+ps)                 ; %(m,n,v)

% 4 <= h < 5 :(n+v,m-v,v)->(n,m,v)
idx = find(hp==4);
CData(idx)      = +CData(idx)                -CData(idx+2*ps); %(n,m-v,v)
CData(idx+ps)   =             +CData(idx+ps) +CData(idx+2*ps); %(n,m,v)

% 5 <= h < 6 :(n+v,m-v,v)->(v,m,n)
idx = find(hp==5);
CData(idx+2*ps) = +CData(idx)                -CData(idx+2*ps); %(n+v,m-v,n)
CData(idx)      = +CData(idx)                -CData(idx+2*ps); %(v,m-v,n)
CData(idx+ps)   = +CData(idx) +CData(idx+ps)                 ; %(v,m,n)

end % End of function DColorColorFunctions 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gridLevels = DColorGrid(w)
% Get contour levels for a grid overlay on a DCOLOR plot.
%
% The grid overlay on a DCOLOR plot is created by calling CONTOUR() twice
% with the real and imaginary parts of W.  Because DCOLOR plots often have
% a large dynamic range, DCOLORGRID tries to calculate log-scale contour
% levels (for both positive and negative data) which are then supplied to
% the CONTOUR routine in place of the default values.  It doesn't always
% work out perfectly, but usually looks better than the default output of
% CONTOUR.

gridLevels = cell(1,2);

% Loop twice, once for real(w), once for imag(w)
for m = 1:2

    % Zero-crossing flag
    zeroFlag = 0;
    
    if m == 1
        minW = min(real(w(isfinite(real(w(:))))));
        maxW = max(real(w(isfinite(real(w(:))))));
        [modeWP fP] = mode(floor(iplog10(real(w(real(w(:))>0)))));
        fZ = length(find(real(w(:))==0));
        [modeWN fN] = mode(floor(iplog10(-real(w(real(w(:))<0)))));
    else
        minW = min(imag(w(isfinite(imag(w(:))))));
        maxW = max(imag(w(isfinite(imag(w(:))))));
        [modeWP fP] = mode(floor(iplog10(imag(w(imag(w(:))>0)))));
        fZ = length(find(imag(w(:))==0));
        [modeWN fN] = mode(floor(iplog10(-imag(w(imag(w(:))<0)))));
    end

    modeWDec = [modeWP 0 modeWN];
    f = [fP fZ fN];
    modeSgn = [1 0 -1];
    if length(find(f==max(f))) == 1
        modeWDec = modeWDec(find(f==max(f),1));
        modeSgn = modeSgn(find(f==max(f),1));
    else
        modeWDec = 0;
        modeSgn = 1;
    end
   
    if maxW - minW < 100
        % If the range is small, try to get 'nice' contour levels
        normRange = (maxW-minW)/10^floor(iplog10(maxW-minW));
        if normRange  < 1.2
            stepVal = (maxW-minW)/(normRange*10);
        elseif normRange < 2.4
            stepVal = (maxW-minW)/(normRange*5);
        else
            stepVal = (maxW-minW)/(normRange*4);
        end
        
        if minW < 0 && maxW > 0
            levels = [-stepVal:-stepVal:minW 0:stepVal:maxW];
        elseif minW < 0
            startVal = minW - (stepVal - mod(-minW,stepVal));
            levels = (startVal+stepVal):stepVal:maxW;
        else
            startVal = minW + (stepVal - mod(minW,stepVal));
            levels = startVal:stepVal:maxW;
        end

    else
        % If range is large, do log-scale contours centered about the mode
        minWDec = floor(iplog10(abs(minW)));
        minSgn = sign(minW);
        maxWDec = floor(iplog10(abs(maxW)));
        maxSgn = sign(maxW);

        % Depending on sign and zero-crossing, find a span of 5 decades
        % about the mode to do contour levels on, plus outer contour levels
        % extending another two decades
        if minSgn == -1 && maxSgn == -1
            % No zero-crossing, all negative values
            zeroFlag = 0;
            startDec = max(min(minWDec,modeWDec+2),min(maxWDec+5,minWDec));
            endDec = max(maxWDec,startDec-5);
            decades = startDec:-1:endDec;
            signs = -1*ones(size(decades));
            levels = -1*10.^((startDec+2):-1:(endDec-2));
        elseif minSgn == -1
            % Data range crosses zero
            minWDec = -max(minWDec,0);
            maxWDec = max(maxWDec,0);
            modeWDec = max(modeWDec,0)*modeSgn;
            startDec = min(max(minWDec,modeWDec-2),max(maxWDec-5,minWDec));
            endDec = min(maxWDec,startDec+5);
            if sign(startDec+.1) == sign(endDec-.1)
                decades = abs(startDec:endDec);
                signs = sign(startDec+.1)*ones(size(decades));
                levelVecStart = min(startDec,endDec);
                levelVecEnd = max(startDec,endDec);
                levelVec = (levelVecStart-2):(levelVecEnd+2);
            else
                zeroFlag = 1;
                decades = [-startDec:-1:0 0:endDec];
                signs = [-1*ones(1,-startDec+1) ones(1,endDec+1)];
                levelVec = [(startDec-2):startDec endDec:(endDec+2)];
            end
            levels = sign(levelVec+.1).*10.^abs(levelVec);
        else
            % No zero-crossing, all positive values
            zeroFlag = 0;
            startDec = min(max(minWDec,modeWDec-2),max(maxWDec-5,minWDec));
            endDec = min(maxWDec,startDec+5);
            decades = startDec:endDec;
            signs = ones(size(decades));
            levels = 10.^((startDec-2):(endDec+2));
        end

        % Now add 'nice' values for levels over the span given by decades
        if length(decades) <= 4;
            logVec = iplog10([1 2 3 4 5 6 7 8 9 10]);
        else
            logVec = iplog10([1 2 5 10]);
        end
        
        for k = 1:length(decades)
            levels = [levels signs(k)*10.^(decades(k)+logVec)]; %#ok
        end
        
        % If the data crosses 0, add a contour at 0
        if zeroFlag == 1;
            levels = [levels 0]; %#ok
        end
        
        % Filter out redundant and out-of-range levels
        levels = unique(levels);
        levels(levels<minW) = [];
        levels(levels>maxW) = [];
        
    end
    
     gridLevels{m} = levels;
end   

end % end of function DColorGrid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
