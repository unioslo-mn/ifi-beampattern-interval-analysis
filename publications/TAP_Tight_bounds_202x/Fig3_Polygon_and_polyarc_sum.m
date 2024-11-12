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

%% Plot polygonal sum

% Set plot parameters
defLineWidth = 1.5;
defFontSize = 25;
pntFontSize = 20;
vectorSize = 0.5;

% Initialize figure
figure(1);
clf;hold on;axis equal;fontsize(20,'point')

% Plot intervals, vertices and arcs
p1 = aI.plot('g','linewidth',defLineWidth,'DisplayName','Concave polyarc');
p2 = gI.plot('b.--','linewidth',defLineWidth,'DisplayName','Convex polygon');
p3 = gaI_vert.plot('bo','MarkerSize',40,'DisplayName','Polygonal vertex');

% Plot trimmed parts of the concave polyarcular sum
arc1 = [aI(1).Arcs ; aI(1).Vertices];
arc2 = [aI(2).Arcs ; aI(2).Vertices];
edge1 = aI(1).Edges;
edge2 = aI(2).Edges;    
p4 = plot(arc1 + arc2.','g--','DisplayName','Trimmed segments');
plot(arc1 + edge2.','g--');
plot(edge1 + arc2.','g--');

% Gauss maps
p5 = gaI_vert.plotGaussMap(vectorSize,'b','DisplayName','Gauss map');

% Axis Labels
xlabel('Real','fontsize',defFontSize,'HorizontalAlignment','center')
ylabel('Imag','fontsize',defFontSize,'HorizontalAlignment','center')

% Interval labels
textLoc = pI.Abs.mid.*exp(1j*pI.Angle.mid);
text( real(textLoc(1)) , imag(textLoc(1)) , '$A^I$'  , ...
                    'fontsize',defFontSize, ...
                    'HorizontalAlignment','center',...
                    'Interpreter','latex')
text( real(textLoc(2)) , imag(textLoc(2))-0.5 , '$B^I$'  , ...
                    'fontsize',defFontSize,'HorizontalAlignment','center',...
                    'Interpreter','latex')
text( aI(3).Real.mid , aI(3).Imag.mid , '$A^I+B^I$'  , ...
                    'fontsize',defFontSize,'HorizontalAlignment','center',...
                    'Interpreter','latex')

% Tightness values
annotText = {['$\tau_A=' num2str(100*aI(1).Area ./ gI(1).Area,2),'\%$'],...
             ['$\tau_B=' num2str(100*aI(2).Area ./ gI(2).Area,2),'\%$'],...
             ['$\tau_{A+B}=' num2str(100*aI(3).Area ./ gI(3).Area,2),'\%$']};
annotation('textbox',[0.73 0.81 0.16 0.10],'String',annotText, ...
           'fontSize',pntFontSize,'HorizontalAlignment','right', ...
           'Interpreter','latex');

% Vertex labels of the A interval
vert = gaI.Vertices;
idxI = 1;
pnt = vert{idxI}.Center - 0.5*exp(1j*vert{idxI}.GaussMap.mid);
for idxP = 1:length(vert{idxI})
    text(real(pnt(idxP)) , imag(pnt(idxP)), ...
         ['[' num2str(idxP) ']'],...
         'color','b','fontsize',pntFontSize...
         ,'HorizontalAlignment','center')
end

% Vertex labels of the A interval
vert = gaI.Vertices;
idxI = 2;
pnt = vert{idxI}.Center - 0.5*exp(1j*vert{idxI}.GaussMap.mid);
for idxP = 1:length(vert{idxI})
    text(real(pnt(idxP)) , imag(pnt(idxP)), ...
         ['(' num2str(idxP) ')'],...
         'color','b','fontsize',pntFontSize...
         ,'HorizontalAlignment','center')
end

% Vertex labels of the A+B interval
idxI = 3;
capGauss31 = width(cap(vert{3}.GaussMap,vert{1}.GaussMap.'));
capGauss32 = width(cap(vert{3}.GaussMap,vert{2}.GaussMap.'));
for idxP = 1:length(vert{idxI})
    if any(~isnan(capGauss31(idxP,:))) && any(~isnan(capGauss32(idxP,:)))
        if ciat.RealInterval(-3,7).isin( imag(vert{idxI}.Center(idxP)))
            pnt = vert{idxI}.Center(idxP) - 1*exp(1j*vert{idxI}.GaussMap(idxP).mid);
        else
            pnt = vert{idxI}.Center(idxP) + 1*exp(1j*vert{idxI}.GaussMap(idxP).mid);
        end
        [~,vIdx1] = max(capGauss31(idxP,:));
        [~,vIdx2] = max(capGauss32(idxP,:));
        text(real(pnt) , imag(pnt), ...
             ['[' , num2str(vIdx1) ']+(' , num2str(vIdx2) ')'],...
             'color','b','fontsize',pntFontSize...
         ,'HorizontalAlignment','center')
    end
end

xlim([-8 7])
legend([p1(1),p4(1),p2(1),p3(1),p5(1)],'Location','northwest')



%% Plot polyarcular sum

% Set plot parameters
defLineWidth = 1.5;
arcLineWidth = 3;
defFontSize = 25;
pntFontSize = 20;
vectorSize = 0.5;

% Initialize figure
figure(2);clf;hold on;axis equal;fontsize(20,'point')

% Plot intervals, vertices and arcs
p0 = aI.plot('g','linewidth',defLineWidth,'DisplayName','Concave polyarc');
p1 = gI.plot('b:','linewidth',defLineWidth,'DisplayName','Convex polygon');
p2 = xI.plot('r--','linewidth',defLineWidth,'DisplayName','Convex polyarc');
p3 = xI_arc.plot('r-','linewidth',arcLineWidth,'DisplayName','Polyarcular arc');
p4 = xI_vert.plot('r','MarkerSize',30,'DisplayName','Polyarcular vertex');

% Gauss maps
p5 = xI_arc.plotGaussMap(vectorSize,'r','DisplayName','Gauss map');
xI_vert.plotGaussMap(vectorSize,'r');

% Axis labels
xlabel('Real','fontsize',defFontSize,'HorizontalAlignment','center')
ylabel('Imag','fontsize',defFontSize,'HorizontalAlignment','center')

% Interval labels
labInt = {'A','B','C'};
textLoc = pI.Abs.mid.*exp(1j*pI.Angle.mid);
text( real(textLoc(1)) , imag(textLoc(1)) , '$A^I$'  , ...
                'fontsize',defFontSize, ...
                'HorizontalAlignment','center',...
                'Interpreter','latex')
text( real(textLoc(2)) , imag(textLoc(2))-0.5 , '$B^I$'  , ...
                'fontsize',defFontSize, ...
                'HorizontalAlignment','center',...
                'Interpreter','latex')
text( aI(3).Real.mid , aI(3).Imag.mid , '$A^I+B^I$'  , ...
                'fontsize',defFontSize, ...
                'HorizontalAlignment','center',...
                'Interpreter','latex')

% Tightness values
annotText = {['$\tau_A=' num2str(100*aI(1).Area ./ xI(1).Area,2),'\%$'],...
             ['$\tau_B=' num2str(100*aI(2).Area ./ xI(2).Area,2),'\%$'],...
             ['$\tau_{A+B}=' num2str(100*aI(3).Area ./ xI(3).Area,2),'\%$']};
annotation('textbox',[0.73 0.81 0.16 0.10],'String',annotText, ...
           'fontSize',pntFontSize,'HorizontalAlignment','right', ...
           'Interpreter','latex');
   
% Polyarcular arcs
arc = cell(3,1);
for idxI = 1:3
    arc{idxI} = [xI(idxI).Arcs{:} ; xI(idxI).Vertices{:}];
    [~,iSort] = sort(arc{idxI}.GaussMap.sup);
    arc{idxI} = arc{idxI}(iSort);
end

% Arc labels of interval A
idxI = 1;
pnt = arc{idxI}.Midpoint - 0.4*exp(1j*arc{idxI}.GaussMap.mid);
pntCorr = [0.56+4j,0,0.2j,0,0];
for idxP = 1:length(arc{idxI})-1
    text(real(pnt(idxP)+pntCorr(idxP)) , imag(pnt(idxP)+pntCorr(idxP)), ...
        ['[' num2str(idxP) ']'],...
         'color','r','fontsize',pntFontSize, ...
         'HorizontalAlignment','center')
end

% Arc labels of interval B
idxI = 2;
pnt = arc{idxI}.Midpoint - 0.4*exp(1j*arc{idxI}.GaussMap.mid);
pntCorr = [-0.2j,0,0,0,0];
for idxP = 1:length(arc{idxI})-1
    text(real(pnt(idxP)+pntCorr(idxP)) , imag(pnt(idxP)+pntCorr(idxP)), ...
        ['(' num2str(idxP) ')'],...
         'color','r','fontsize',pntFontSize, ...
         'HorizontalAlignment','center')
end

% Arc labels of interval A+B
capGauss31 = width(cap(arc{3}.GaussMap,arc{1}.GaussMap.'))>10*eps;
capGauss32 = width(cap(arc{3}.GaussMap,arc{2}.GaussMap.'))>10*eps;
cntListA = [1:length(arc{1})-1,1];
cntListB = [1:length(arc{2})-1,1];
for idxP = 2:length(arc{3})
    if any(capGauss31(idxP,:)) && any(capGauss32(idxP,:))
        % pnt = arc{3}.Midpoint - 1*exp(1j*arc{3}.GaussMap.mid);
        if ciat.RealInterval(-3,7).isin( imag(arc{3}.Midpoint(idxP) ))
            pnt = arc{3}.Midpoint(idxP) - 1*exp(1j*arc{3}.GaussMap(idxP).mid);
        else
            pnt = arc{3}.Midpoint(idxP) + 1*exp(1j*arc{3}.GaussMap(idxP).mid);
        end
        cntA = cntListA(find(capGauss31(idxP,:),1));
        cntB = cntListB(find(capGauss32(idxP,:),1));
        text(real(pnt) , imag(pnt), ...
             ['[' , num2str(cntA), ']+(' , num2str(cntB) ,')'],...
             'color','r','fontsize',pntFontSize,...
             'HorizontalAlignment','center')
    end
end

xlim([-8 7])
legend([p0(1),p1(1),p2(1),p3(1),p4(1),p5(1)],'Location','northwest')
