clear
% close all

%% Set parameters

% Array
M = 100;
N = 1;

% Tapering
%apodM = taylorwin(M,3,-25);
%apodN = taylorwin(N,3,-25);
apodM = ones(M,1);
apodN = ones(N,1);
apod = apodN .* apodM.' / (sum(apodM) * sum(apodN));


% Errors
errAmp = 0.1;
errAng = deg2rad(10);

gamma = 0.05;

% Plot
dinRange = [-60 5];

% Tolerance
polygon_tolerance = 1e-4;


polygon_tolerance1 = 0.5*1e-3;  % 1/(1*M);
polygon_tolerance2 = 1e-4;  % 1 /(100*M);
polygon_tolerance3 = 1e-5;  % 1 / (10000*M);

%% Calculate beampattern

% Angular axes
K = 1;
L = 1;

if K == 1
    u = 0;
else
    u = linspace(-1,1,K);
end

if L == 1
    v = 0;
else
    v = linspace(-1,1,L);
end

%u = 0.36;
u = 0.029;

% nomBP = zeros(K,L);
% for k = 1:K
%     for l = 1:L
%         delay = pi*((1:M)*u(k) + (1:N)'*v(l));
%         W = apod .* exp(1j*delay);
%         nomBP(k,l) = abs(sum(W,'all'));
%     end
% end


%plot(u, db(nomBP));


%% Calculate beampattern interval for a given intersection


arx_area = zeros(1,M);
arcular_area = zeros(1,M);
gon_area1 = zeros(1,M);
gon_area2 = zeros(1,M);
gon_area3 = zeros(1,M);

arx_range = zeros(1,M);
gon_range1 = zeros(1,M);
gon_range2 = zeros(1,M);
gon_range3 = zeros(1,M);

arx_count = zeros(1,M);
arcular_count = zeros(1,M);
gon_count1 = zeros(1,M); 
gon_count2 = zeros(1,M);
gon_count3 = zeros(1,M);

for m = 1:M
    
        % Calculate phase delay
        delay = pi*m*u;

        % Generate element intervals in polar type
        pI = ciat.PolarInterval(apod(m)*(1-errAmp) , apod(m)*(1+errAmp), delay-errAng , delay + errAng);

        % Cast element intervals to rectangular and calculate sum
        %rI = ciat.RectangularInterval(pI);
        %rIsum = sum(rI,'all');
        %recBP(k,l) = abs(rIsum);
        
        cI = ciat.CircularInterval(1,gamma);


        % Cast element intervals to polygon and calculate sum
        gI1 = ciat.PolygonalInterval(pI,cI,'tolerance', polygon_tolerance1);
        gI2 = ciat.PolygonalInterval(pI,cI,'tolerance', polygon_tolerance2);
        gI3 = ciat.PolygonalInterval(pI,cI,'tolerance', polygon_tolerance3);
        if m == 1
            gon_sum1 = ciat.PolygonalInterval(pI,cI,'tolerance', polygon_tolerance1);
            gon_sum2 = ciat.PolygonalInterval(pI,cI,'tolerance', polygon_tolerance2);
            gon_sum3 = ciat.PolygonalInterval(pI,cI,'tolerance', polygon_tolerance3);
        else 
            gon_sum1 = gon_sum1 + gI1;
            gon_sum2 = gon_sum2 + gI2;
            gon_sum3 = gon_sum3 + gI3;
        end

        gon_area1(m) = gon_sum1.Area;
        gon_count1(m) = gon_sum1.PointCount;
        gon_range1(m) = gon_sum1.Abs.Supremum;
        gon_area2(m) = gon_sum2.Area;
        gon_count2(m) = gon_sum2.PointCount;
        gon_range2(m) = gon_sum2.Abs.Supremum;
        gon_area3(m) = gon_sum3.Area;
        gon_count3(m) = gon_sum3.PointCount;     
        gon_range3(m) = gon_sum3.Abs.Supremum;
        
        % Cast element intervals to polygon and calculate sum
        xI = ciat.PolyarxInterval(pI);
        if m == 1
            arx_sum = ciat.PolyarxInterval(pI,cI);
        else
            arx_sum = arx_sum + ciat.PolyarxInterval(pI,cI);
        end
        arx_area(m) = arx_sum.Area;
        arx_count(m) = arx_sum.ArxCount;
        arx_range(m) = arx_sum.Abs.Supremum;
       

end



%%

fig=figure(1);clf
fig.Position=[100 179 495 290];%[100 179 495 321];
set(gcf,'color','w');
set(gca,'LineWidth',2)
set(gca,'DefaultLineLineWidth',2)
set(gca,'defaultAxesFontSize',20)
set(gca, 'fontsize', 15)
hold on;

tr = 0.0037;

plot(ciat.PolyarcularInterval(pI,cI), 'displayname', 'Annular sector','color','k','LineWidth',1,'handlevisibility','off'); hold on;
plot([0,0],[0,eps], 'displayname', 'Annular sector','color','k','LineWidth',1);
plot(ciat.PolyarxInterval(pI,cI),'k--', 'handlevisibility', 'off','LineWidth',2);
plot([0],[0], 'k--','displayname','Polyarc')
%axis equal

aI = ciat.PolyarcularInterval(pI,cI) + 1*tr; 
plot(aI, 'handlevisibility', 'off','color','k','LineWidth',1); hold on;
plot(gI3+1*tr,'-o', 'displayname', 'Polygon tol. 1e-5','Color',[0.9290 0.6940 0.1250]);

aI = ciat.PolyarcularInterval(pI,cI) + 2*tr; 
plot(aI, 'handlevisibility', 'off','color','k','LineWidth',1); hold on;
plot(gI2+2*tr,'-o', 'displayname', 'Polygon tol. 1e-4','Color',[0.8500 0.3250 0.0980]);

aI = ciat.PolyarcularInterval(pI,cI) + 3*tr; 
plot(aI, 'handlevisibility', 'off','color','k','LineWidth',1); hold on;
plot(gI1+ 3*tr,'-o','color', [1 0.2 1], 'displayname', 'Polygon tol. 5e-4','Color',[0 0.4470 0.7410]);

xlabel("Real")
ylabel("Imag")
ylim([4.8869e-04 0.0087])
xlim([-0.0121 0.0041])
legend('NumColumns',2)
legend boxoff;
%exportgraphics(fig,strcat('/Users/havarn/Documents/DocGit/PolyArc/ifi-beampattern-interval-analysis/publications/TAP_Tight_bounds_202x/Arnestad_plots/Arr_4_1c.pdf'))

%%
fig=figure(2);clf
set(gcf,'color','w');
fig.Position=[100 179 495 321];
set(gca,'LineWidth',2)
set(gca,'DefaultLineLineWidth',2)
set(gca,'defaultAxesFontSize',20)
set(gca, 'fontsize', 15)
hold on;

plot([0],[0], 'k--','displayname','Polyarc','linewidth',2)
plot(gon_sum3,'-o', 'displayname', 'Polygon tol. 1e-5','Color',[0.9290 0.6940 0.1250]);
plot(gon_sum2,'-o', 'displayname', 'Polygon tol. 1e-4','Color',[0.8500 0.3250 0.0980]);
plot(gon_sum1,'-o', 'displayname', 'Polygon tol. 5e-4','Color',[0 0.4470 0.7410]     ); hold on;

plot(arx_sum,'k--', 'handlevisibility', 'off','linewidth',2);

%xl [-0.1826 -0.1134]
%yl [0.3721 0.4024]
xlim([-0.1966 -0.1433])
ylim([0.3100 0.3594])
xlabel("Real")
ylabel("Imag")
legend('location','northwest')
legend boxoff;
%axis equal

axes('position',[.55 .2 .4 .4], 'NextPlot', 'add')
box on
plot(gon_sum1,'-', 'displayname', 'Polygon','Color',[0 0.4470 0.7410]     ,'linewidth',1.5); hold on;
plot(gon_sum2,'-', 'displayname', 'Polygon','Color',[0.8500 0.3250 0.0980],'linewidth',1.5);
plot(gon_sum3,'-', 'displayname', 'Polygon','Color',[0.9290 0.6940 0.1250],'linewidth',1.5);
plot(arx_sum,'k--', 'handlevisibility', 'off','linewidth',1.5);
axis equal
axis tight
plot([-0.1966 -0.1433 -0.1433 -0.1966 -0.1966], [0.3594  0.3594 0.3100  0.3100 0.3594],'k-', 'linewidth',2)
title("Full view")

%exportgraphics(fig,strcat('/Users/havarn/Documents/DocGit/PolyArc/ifi-beampattern-interval-analysis/publications/TAP_Tight_bounds_202x/Arnestad_plots/Arr_4_2c.pdf'))

%%
fig=figure(3);clf
set(gcf,'color','w');
fig.Position=[100 179 495 321];
set(gca,'LineWidth',2)
set(gca,'DefaultLineLineWidth',2)
set(gca,'defaultAxesFontSize',20)
set(gca, 'fontsize', 15)
hold on;

correction = 1:100;
plot(arx_count ,'k--', 'displayname', "Polyarc"); hold on;
plot(gon_count3     , 'Color',[0.9290 0.6940 0.1250],'displayname', "Polygon tol. 1e-5");
plot(gon_count2     , 'Color',[0.8500 0.3250 0.0980],'displayname', "Polygon tol. 1e-4");
plot(gon_count1     , 'Color',[0 0.4470 0.7410]     ,'displayname', "Polygon tol. 5e-4"); hold on;
plot(arx_count ,'k--', 'handlevisibility', 'off'); hold on;


%plot(arcular_count,'r--')

ylim([0,2000])
xlabel("Number of elements")
ylabel("Number of vertices or arcs")
legend('location','northwest')
legend boxoff;
%exportgraphics(fig,strcat('/Users/havarn/Documents/DocGit/PolyArc/ifi-beampattern-interval-analysis/publications/TAP_Tight_bounds_202x/Arnestad_plots/Arr_4_3c.pdf'))

%%
fig=figure(4);clf
set(gcf,'color','w');
fig.Position=[100 179 495 321];
set(gca,'LineWidth',2)
set(gca,'DefaultLineLineWidth',2)
set(gca,'defaultAxesFontSize',20)
set(gca, 'fontsize', 15)
hold on;


plot(1:M, (arx_area ./ arx_area   -1)*100, 'k--','displayname', "Polyarc")
plot(1:M, ((gon_area3 ./ arx_area) -1)*100, 'Color',[0.9290 0.6940 0.1250],'displayname', "Polygon tol. 1e-5"); hold on;
plot(1:M, ((gon_area2 ./ arx_area) -1)*100, 'Color',[0.8500 0.3250 0.0980],'displayname', "Polygon tol. 1e-4")
plot(1:M, ((gon_area1 ./ arx_area) -1)*100, 'Color',[0 0.4470 0.7410]     ,'displayname', "Polygon tol. 5e-4")


xlabel("Number of elements")
ylabel("Area relaxation (%)")
legend('location','east')
legend boxoff;
ylim([-0.01,inf])
%'markerfacecolor', [0.9290 0.6940 0.1250], 'displayname', "\it A", 'markersize', 13, 

plot([60 80 80 60 60], [0.5 0.5 -0.01 -0.01 0.5],'k-', 'linewidth',2,'HandleVisibility','off')
plot([65 50], [0.15 5],'k-', 'linewidth',2,'HandleVisibility','off')

axes('position',[.25 .4 .3 .3], 'NextPlot', 'add')
box on
plot(1:M, (arx_area ./ arx_area   -1)*100, 'k--','displayname', "Polyarc",'linewidth',2)
plot(1:M, ((gon_area3 ./ arx_area) -1)*100, 'Color',[0.9290 0.6940 0.1250],'displayname', "Polygon tol. 1e-5",'linewidth',2); hold on;


%axis equal
axis tight
title("Zoomed view")
ylim([-0.01, 0.5])
xlim([60 80])
set(gca,'LineWidth',2)

%exportgraphics(fig,strcat('/Users/havarn/Documents/DocGit/PolyArc/ifi-beampattern-interval-analysis/publications/TAP_Tight_bounds_202x/Arnestad_plots/Arr_4_4c.pdf'))

%%
% subplot(2,2,3)
% plot(u,db(recBP.Infimum),'r-','displayname','rect')
% plot(u,db(recBP.Supremum),'r-','displayname','rect')
% plot(u,db(gonBP.Infimum),'b-','displayname','poly')
% plot(u,db(gonBP.Supremum),'b-','displayname','poly')
% plot(u,db(arxBP.Infimum),'g--','displayname','arx')
% plot(u,db(arxBP.Supremum),'g--','displayname','arx')
% plot(u,db(nomBP),'k-','displayname','nom')
% ylim(dinRange)



%%
% tolerance_index(recBP.Supremum, nomBP, recBP.Infimum, u, v)