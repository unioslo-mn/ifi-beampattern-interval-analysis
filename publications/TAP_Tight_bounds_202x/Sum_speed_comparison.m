clear
% close all

%% Generate two random polar intervals

% Define random polar intervals
N = 400;
absMin = rand(N,1);
absMax = absMin + rand(N,1)/10;
angMin = 2*pi*rand(N,1);
angMax = angMin + 2*pi*rand(N,1)/20;
pI = ciat.PolarInterval(absMin,absMax,angMin,angMax);

% Convert them to various other formats
prI = ciat.RectangularInterval(pI);
pgI = ciat.PolygonalInterval(pI);
paI = ciat.PolyarcularInterval(pI);
pxI = ciat.PolyarxInterval(pI);

%% Sum rectangular intervals

recSeg = 4*ones(N,1);
tic
prIsum = sum(prI);
recTime = toc;
recTime = recTime * ones(N,1);


%% Sum polygonal intervals

% Measure time and count segments
gonSeg = zeros(N,1);
gonTime = zeros(N,1);
pgIsum = pgI(1);
gonSeg(1) = length(pgIsum.Points{:});
tic
for n = 2:N
    % tic
    pgIsum = gonPlus(pgIsum,pgI(n));
    gonSeg(n) = length(pgIsum.Points{:});
    gonTime(n) = toc;
end


%% Sum polyarx intervals

% Prepare objects
arx = pxI.Arx;

% Measure time
arxSeg = zeros(N,1);
arxTime = zeros(N,1);
pxIsum = arx{1};
arxSeg(1) = size(pxIsum,1);
tic
for n = 2:N
    % tic
    pxIsum = arxPlus(pxIsum,arx{n});
    arxSeg(n) = size(pxIsum,1);
    arxTime(n) = toc;
end
pxIsum = ciat.PolyarxInterval(pxIsum);


%% Sum polyarcular intervals

% Measure time
arcSeg = zeros(N,1);
arcTime = zeros(N,1);
paIsum = paI(1);
arcSeg(1) = length(paIsum.DefArcs{:});
tic
for n = 2:N
    % tic
    paIsum = arcPlus(paIsum,paI(n));
    arcSeg(n) = length(paIsum.DefArcs{:});
    arcTime(n) = toc;
end



%% Report

% Measure areas
recArea = prIsum.Area;
gonArea = pgIsum.Area;
arcArea = paIsum.Area;
arxArea = pxIsum.Area;

% Print report
sprintf(['Rectangular interval: Area: %0.4f (tightness: %0.1f%%), Time: %0.1fms\n'...
         'Polygonal interval area: %0.4f (tightness: %0.1f%%), Time: %0.1fms\n'...
         'Polyarcular (convex) area: %0.4f (tightness: %0.1f%%), Time: %0.1fms\n',...
         'Polyarcular (concave) area: %0.4f (tightness: %0.1f%%), Time: %0.1fms'], ...
         recArea, arcArea / recArea * 100, recTime(end)*1e3,...
         gonArea, arcArea / gonArea * 100, gonTime(end)*1e3,...
         arxArea, arcArea / arxArea * 100, arxTime(end)*1e3,...
         arcArea, arcArea / arcArea * 100, arcTime(end)*1e3)

%% Plot
% figure;clf
subplot(1,2,1);cla;hold on;axis equal;title('Complex plane')
paI.plot('g','DisplayName','pol');
paIsum.plot('b','linewidth',2,'DisplayName','rec');
pgIsum.plot('k','linewidth',2,'DisplayName','gon');
pxIsum.plot('c','linewidth',2,'DisplayName','arx');
paIsum.plot('r','linewidth',2,'DisplayName','arc');
legend(legendUnq)

subplot(1,2,2);cla;hold on;title('Number of interval boundary segments')
plot(recSeg,'b-')
plot(gonSeg,'k-')
plot(arxSeg,'c-')
plot(arcSeg,'r-')
% legend('recSeg','gonSeg','arxSeg','arcSeg')
xlabel('Number of interval additions')
ylabel('Segment count')
yyaxis right ; cla
plot(recTime,'b--')
plot(gonTime,'k--')
plot(arxTime,'c--')
plot(arcTime,'r--')
ylim([0,gonTime(end)])
ylabel('Time')
legend('recSeg','gonSeg','arxSeg','arcSeg','recTime','gonTime','arxTime','arcTime')


















%% Stripped down polygonal sum algorithm

function r = gonPlus(obj1,obj2)

    v = obj1.Points{:};
    w = obj2.Points{:};
    
    % Handle exception when one of the inputs is a degenerate interval
    if length(v)==1 || length(w)==1 
        points = reshape(v + w.' ,[],1);
        r = ciat.PolygonalInterval(points);
        return
    end
    
    i = 1; j = 1; 
    eps10 = 10*eps; % needed so avoid skipping vertices due to numerical precision
    
    I = length(obj1.Points{:}) + 1; 
    v = [v; v(1); v(2)];
    v_arg = ciat.wrapTo2Pi( angle( v(2:end) - v(1:end-1) ));
    v_arg(end) = v_arg(end) + 2*pi; % otherwise it wraps around
    
    J = length(obj2.Points{:}) + 1; 
    w = [w; w(1); w(2)];
    w_arg = ciat.wrapTo2Pi( angle( w(2:end) - w(1:end-1) ));
    w_arg(end) = w_arg(end) + 2*pi; % otherwise it wraps around
    
    p = zeros( I + J, 1);
    n = 0;
    
    % v_arg(i) = angle(v_i+1 - v_i), which is why we also
    % repeat when i=I and j=J (as opposed to in DeBerg2008). 
    % The loop breaks in the following turn.
    while (i <= I) && (j <= J) % continue finding more points
        n = n+1;
        p(n) = v(i) + w(j); % add what we found
    
        if     v_arg(i) < w_arg(j) + eps10 
            i = i + 1;
        elseif v_arg(i) > w_arg(j) + eps10
            j = j + 1;
        else
            i = i + 1;            j = j + 1;
        end            
    end
    
    r = ciat.PolygonalInterval(p(1:n-1));

end

%% Stripped down convex polyarcular plus algorithm

function r = arxPlus(arx1,arx2)

    % Extract parameters
    cen1 = complex(arx1(:,1),arx1(:,2));
    rad1 = arx1(:,3);
    ang1 = arx1(:,4);
    cen2 = complex(arx2(:,1),arx2(:,2));
    rad2 = arx2(:,3);
    ang2 = arx2(:,4);
    
    % Taken from the polygonal plus function
    N1 = size(ang1,1);
    N2 = size(ang2,1);
    N3 = N1 + N2;
    cen3 = zeros(N3,1);
    rad3 = zeros(N3,1);
    ang3 = zeros(N3,1);
    n1 = 1;
    n2 = 1;
    n3 = 0;
    eps10 = eps*10;
    while (n1 <= N1) && (n2 <= N2) % continue finding more points
        n3 = n3 + 1;
    
        % Sum arcs
        cen3(n3) = cen1(n1) + cen2(n2);
        rad3(n3) = rad1(n1) + rad2(n2);
        ang3(n3) = min(ang1(n1),ang2(n2));
        
        % Increment index
        if  ang1(n1) < ang2(n2) + eps10 
            n1 = n1 + 1;
        elseif ang1(n1) > ang2(n2) + eps10
            n2 = n2 + 1;
        else
            n1 = n1 + 1;
            n2 = n2 + 1;
        end            
    end
    cen3 = cen3(1:n3,:);
    rad3 = rad3(1:n3,:);
    ang3 = ang3(1:n3,:);
       
    % Generate polyarc
    r = [real(cen3),imag(cen3),rad3,ang3];
end

%% Stripped down concave polyarc plus algorithm

function r = arcPlus(obj1,obj2)
    
    % Extract curve segments by type
    %   - extract arcs including vertices
    arc1 = [obj1.Arcs{:} ; obj1.Vertices{:}];
    arc2 = [obj2.Arcs{:} ; obj2.Vertices{:}];
    %   - extract edges with non-zero length
    edge1 = obj1.Edges{:};
    edge2 = obj2.Edges{:};    
    
    % Add arcs and vertices
    arcPlusArc = arc1 + arc2.';
    arcPlusEdge = arc1 + edge2.';
    edgePlusArc = edge1 + arc2.';
    
    % Extract valid segments
    arc3 = arcPlusArc(~isnan(arcPlusArc));
    edge3 = [ arcPlusEdge(~isnan(arcPlusEdge)) ; ...
              edgePlusArc(~isnan(edgePlusArc)) ];

    % Extract non-vertex segments
    arc3 = arc3(abs(arc3.Length)>10*eps);
    edge3 = edge3(abs(edge3.Length)>10*eps);

    % Split segments
    [arc3,edge3] = ciat.PolyarcularInterval.splitSegments(arc3,edge3);
    arc3 = arc3(abs(arc3.Length)>10*eps); % This should be unnecessary
    edge3 = edge3(abs(edge3.Length)>10*eps);
    
    % Trim segments
    arc3 = ciat.PolyarcularInterval.trimSegments(arc3,edge3,1);
    
    % Generate polyarc
    r = ciat.PolyarcularInterval(arc3);
    r = joinSegments(r);

end