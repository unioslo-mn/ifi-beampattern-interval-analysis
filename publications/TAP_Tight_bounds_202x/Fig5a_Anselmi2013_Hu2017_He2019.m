clear
% close all

%% Anselmi 2013 and Hu 2017 (amplitude errors)
N = 10;
ampErr = 0.1;
w = taylorwin(N,3,-20);
% w = ones(N,1);
A = w * (1 + ciat.RealInterval(-ampErr/2,ampErr/2));
dTheta = 0.9;
theta = linspace(-pi,pi,N)' * dTheta;

% Define and cast intervals
AFp = ciat.PolarInterval(A,ciat.RealInterval(theta));
AFr = ciat.RectangularInterval(AFp);
AFg = ciat.PolygonalInterval(AFp);

% Sum intervals
AFr_sum = sum(AFr);
AFg_sum = sum(AFg);

% Anselmi power pattern bounds
P_An = abs(AFr_sum)^2;

% Hu's Taylor based approximation
A_R = sum(A .* cos(theta));
A_I = sum(A .* sin(theta));
P_d = 2*cos(theta) .* A_R + 2*sin(theta) .* A_I;
P_0 = abs(sum(w .* exp(1j*theta)))^2;
P_U = P_0 + sum( abs(P_d.sup) .* A.width/2 );
P_I = P_0 - sum( abs(P_d.inf) .* A.width/2 );
P_Hu = ciat.RealInterval(P_I,P_U);

% He's matrix method
a_mid = A.mid';
a_rad = A.width'/2;
Theta = cos(theta - theta');
P_mid = abs(sum(w .* exp(1j*theta)))^2;
P_inf = P_mid - 2 * abs(a_mid * Theta) * a_rad';
P_sup = P_mid + 2 * abs(a_mid * Theta) * a_rad' + a_rad * abs(Theta)*a_rad';
P_He = ciat.RealInterval(P_inf,P_sup);

% Plot
figure(1);clf;
% subplot(1,2,1);
hold on;axis equal
set(gca,'DefaultLineLineWidth',2)
plot(0,0,'k+')
AFp.plot('g');
AFr.plot('b');
AFr_sum.plot('b');
AFg_sum.plot('g');
xlabel('Real')
ylabel('Imag')
text(AFr_sum.real.mid,AFr_sum.imag.mid,'AF_\Sigma','fontsize',20)
xL = xlim(); yL = ylim();

    % Polygonal power bounds
l1 = fimplicit(@(x,y) x.^2+y.^2-(sup(abs(AFg_sum)))^2,'g--','linewidth',2,...
                                        'DisplayName','Polyarc');
fimplicit(@(x,y) x.^2+y.^2-(inf(abs(AFg_sum)))^2,'g--','linewidth',2)
    % Anselmi
l2 = fimplicit(@(x,y) x.^2+y.^2-P_An.inf,'b--','linewidth',2,...
                                        'DisplayName','Anselmi');
fimplicit(@(x,y) x.^2+y.^2-P_An.sup,'b--','linewidth',2)
    % Hu
l3 = fimplicit(@(x,y) x.^2 + y.^2 - P_Hu.inf,'r:','linewidth',2,...
                                        'DisplayName','Hu');
fimplicit(@(x,y) x.^2 + y.^2 - P_Hu.sup,'r:','linewidth',2)
    %He
l4 = fimplicit(@(x,y) x.^2 + y.^2 - P_He.inf,'r--','linewidth',2,...
                                        'DisplayName','He');
fimplicit(@(x,y) x.^2 + y.^2 - P_He.sup,'r--','linewidth',2)

legend([l1,l2,l3,l4])

xlim(xL); ylim(yL)
fontsize(20,'point')

%% Add zoom window

% create a new pair of axes inside current figure
axes('position',[.83 .18 .07 .5]); hold on; box on
    % Polygonal power bounds
fimplicit(@(x,y) x.^2+y.^2-(sup(abs(AFg_sum)))^2,'g--','linewidth',2)
fimplicit(@(x,y) x.^2+y.^2-(inf(abs(AFg_sum)))^2,'g--','linewidth',2)
    % Anselmi
fimplicit(@(x,y) x.^2+y.^2-P_An.inf,'b--','linewidth',2)
fimplicit(@(x,y) x.^2+y.^2-P_An.sup,'b--','linewidth',2)
    % Hu
fimplicit(@(x,y) x.^2 + y.^2 - P_Hu.inf,'r:','linewidth',2)
fimplicit(@(x,y) x.^2 + y.^2 - P_Hu.sup,'r:','linewidth',2)
    %He
fimplicit(@(x,y) x.^2 + y.^2 - P_He.inf,'r--','linewidth',2)
fimplicit(@(x,y) x.^2 + y.^2 - P_He.sup,'r--','linewidth',2)

axis equal
xlim([1.85 1.92])
fontsize(20,'point')


