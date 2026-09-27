% runAirfoil.m
% =====================================================================
%  Inviscid aerodynamic analysis of the Eppler E420 high-lift airfoil
%  using a 2-D vortex-panel method (Kuethe & Chow).
%
%  Produces:
%    (1) the surface pressure distribution  Cp(x/c)  at a chosen alpha
%    (2) the lift curve                     Cl(alpha)
%
%  Run in MATLAB or MATLAB Online:  >> runAirfoil
% =====================================================================

clear; clc; close all;

Nps       = 100;          % panels per surface (total 2*Nps)
alpha_cp  = 4.0;          % angle of attack for the Cp plot [deg]
alpha_sweep = -6:1:14;    % angles for the lift curve      [deg]

% --- load and panel the airfoil ---
[x, y] = loadAirfoil('eppler420.dat', Nps);
fprintf('Eppler 420 loaded: %d panels\n', numel(x)-1);

% --- (1) pressure distribution at alpha_cp ---
[Cp, Cl, xc, yc] = panelMethod(x, y, alpha_cp);
up = yc >= 0;             % split surfaces for plotting

figure('Color','w','Position',[100 100 640 420]);
plot(xc(up),  Cp(up),  '.-','Color',[0.75 0.22 0.17],'LineWidth',1,'MarkerSize',6); hold on;
plot(xc(~up), Cp(~up), '.-','Color',[0.17 0.43 0.73],'LineWidth',1,'MarkerSize',6);
set(gca,'YDir','reverse'); grid on; yline(0,'k-','Alpha',0.4);
xlabel('x/c'); ylabel('C_p');
legend('upper surface','lower surface','Location','northeast');
title(sprintf('Eppler 420 pressure distribution  (\\alpha = %.0f\\circ,  C_l = %.2f)', alpha_cp, Cl));

% --- (2) lift curve ---
Cls = zeros(size(alpha_sweep));
for k = 1:numel(alpha_sweep)
    [~, Cls(k)] = panelMethod(x, y, alpha_sweep(k));
end
p = polyfit(deg2rad(alpha_sweep), Cls, 1);   % slope (per rad) & intercept

figure('Color','w','Position',[100 100 600 420]);
plot(alpha_sweep, Cls, 'o-','Color',[0.75 0.22 0.17],'MarkerSize',5,'LineWidth',1.2); hold on;
aa = [alpha_sweep(1) alpha_sweep(end)];
plot(aa, polyval(p, deg2rad(aa)), '--','Color',[0.4 0.4 0.4]);
grid on; xlabel('angle of attack  \alpha (deg)'); ylabel('C_l');
legend('panel method', sprintf('linear fit  (%.2f /rad)', p(1)), 'Location','southeast');
title('Eppler 420 lift curve');

% --- console summary ---
Cl0 = interp1(alpha_sweep, Cls, 0);
fprintf('\nResults\n-------\n');
fprintf('  Cl at alpha = 0 deg : %+.3f\n', Cl0);
fprintf('  lift-curve slope    : %.2f /rad  (%.4f /deg)\n', p(1), p(1)*pi/180);
fprintf('  (thin-airfoil theory predicts 2*pi = %.2f /rad)\n', 2*pi);
