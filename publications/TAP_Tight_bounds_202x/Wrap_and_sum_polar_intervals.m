clear
% close all

%% Generate a random annular sector

% Generate polar intervals
r = [6.0 , 5.0 ; ...
     7.0 , 6.0];
a = [0.7 , -0.1 ; ...
     1.1 , 0.2]*pi;
pI = ciat.PolarInterval(r(1,:), r(2,:), a(1,:), a(2,:));

% Wrap in polygon and polyarc
gI = ciat.PolygonalInterval(pI,'tolerance',0.3);
aI = ciat.PolyarcularInterval(pI);
xI = ciat.PolyarxInterval(pI);

% Sum the intervals
gI = [gI sum(gI)];
aI = [aI sum(aI)];
xI = [xI sum(xI)];

% Sample intervals
pIsmp = pI.sample(100);
pIsmp = [pIsmp {cell2mat(pIsmp(1)) + cell2mat(pIsmp(2)).'}];
pIsmp = cell2mat(pIsmp);

% Extract arcs and vertices
gaI = ciat.PolyarcularInterval(gI);
gaI_vert = vertcat(gaI.Vertices{:});
xI_arc = vertcat(xI.Arcs{:});
xI_vert = vertcat(xI.Vertices{:});

%% PLOT RESULTS

% Initialize figure
figure(1);clf;hold on
aI.plot('g','linewidth',1.5);
gI.plot('b.--','linewidth',1.5);
xI.plot('r--','linewidth',1.5);

% Plot vertices
gaI_vert.plot('bo','MarkerSize',40);
xI_arc.plot('r-','linewidth',3);
xI_vert.plot('r','MarkerSize',30);
% scatter(real(xI_vert.Center),imag(xI_vert.Center),50,'ro')

% Gauss maps
gaI_vert.plotGaussMap(0.3,'b');
xI_arc.plotGaussMap(0.3,'r');
xI_vert.plotGaussMap(0.3,'r');

% Samples
% scatter(real(pIsmp),imag(pIsmp),1,'g.')

% Labels
xlabel('Real','fontsize',20)
ylabel('Imag','fontsize',20)
    % Intervals
labInt = {'A','B','C'};
textLoc = pI.Abs.mid.*exp(1j*pI.Angle.mid);
text( real(textLoc(1)) , imag(textLoc(1)) , 'A'  ,'fontsize',20)
text( real(textLoc(2)) , imag(textLoc(2))-0.5 , 'B'  ,'fontsize',20)
text( aI(3).Real.mid-1 , aI(3).Imag.mid , 'A+B'  ,'fontsize',20)
    % Polygonal Vertices
vert = gaI.Vertices;
for idxI = 1:2
    pnt = vert{idxI}.Center + 0.5*exp(1j*vert{idxI}.GaussMap.mid) - 0.2;
    for idxP = 1:length(vert{idxI})
        text(real(pnt(idxP)) , imag(pnt(idxP)), ...
             [labInt{idxI} '_' num2str(idxP)],...
             'color','b','fontsize',15)
    end
end
idxI = 3;
pnt = vert{idxI}.Center + 0.5*exp(1j*vert{idxI}.GaussMap.mid) - 0.5;
capGauss31 = width(cap(vert{3}.GaussMap,vert{1}.GaussMap.'))>10*eps;
capGauss32 = width(cap(vert{3}.GaussMap,vert{2}.GaussMap.'))>10*eps;
for idxP = 1:length(vert{idxI})
    if any(capGauss31(idxP,:)) && any(capGauss32(idxP,:))
        text(real(pnt(idxP)) , imag(pnt(idxP)), ...
             ['A_' , num2str(find(capGauss31(idxP,:)),1) , ...
              '+B_' , num2str(find(capGauss32(idxP,:)),1)],...
             'color','b','fontsize',15)
    end
end
    % Polyarx arcs
arc = xI.Arcs;
for idxI = 1:2
    pnt = arc{idxI}(end).Midpoint - 0.3*exp(1j*arc{idxI}(end).GaussMap.mid) - 0.2;
    text(real(pnt) , imag(pnt), [labInt{idxI} '_0'],...
             'color','r','fontsize',15)
end
idxI = 3;
pnt = arc{idxI}.Midpoint - 0.5*exp(1j*arc{idxI}.GaussMap.mid) - 0.5;
text(real(pnt(1)) , imag(pnt(1)), 'A_0+B_4','color','r','fontsize',15)
text(real(pnt(2)) , imag(pnt(2)), 'A_0+B_1','color','r','fontsize',15)
text(real(pnt(3)) , imag(pnt(3)), 'A_3+B_0','color','r','fontsize',15)
% text(real(pnt(4)) , imag(pnt(4)), 'A_0+B_0','color','r','fontsize',15)

    % Polyarx vertices
vert = xI.Vertices;
vert{2} = circshift(vert{2}(1:4),-1);
for idxI = 1:2
    pnt = vert{idxI}.Center - 0.5*exp(1j*vert{idxI}.GaussMap.mid) - 0.2;
    for idxP = 1:4
        text(real(pnt(idxP)) , imag(pnt(idxP)), ...
             [labInt{idxI} '_' num2str(idxP)],...
             'color','r','fontsize',15)
    end
end
idxI = 3;
pnt = vert{idxI}.Center - 0.5*exp(1j*vert{idxI}.GaussMap.mid) - 0.5;
capGauss31 = width(cap(vert{3}.GaussMap,vert{1}.GaussMap.'))>10*eps;
capGauss32 = width(cap(vert{3}.GaussMap,vert{2}.GaussMap.'))>10*eps;
for idxP = 1:length(vert{idxI})
    if any(capGauss31(idxP,:)) && any(capGauss32(idxP,:))
        text(real(pnt(idxP)) , imag(pnt(idxP)), ...
             ['A_' , num2str(find(capGauss31(idxP,:)),1) , ...
              '+B_' , num2str(find(capGauss32(idxP,:)),1)],...
             'color','r','fontsize',15)
    end
end

axis equal 
% fontsize(20,'point')
