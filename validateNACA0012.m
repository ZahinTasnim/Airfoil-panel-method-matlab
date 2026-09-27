% validateNACA0012.m
% =====================================================================
%  Validation of panelMethod.m against a known case: the symmetric
%  NACA 0012 airfoil. Thin-airfoil theory predicts
%      Cl(alpha = 0) = 0     and     dCl/dalpha ~ 2*pi per radian
%  (panel methods give slightly more, ~6.9/rad, because of thickness).
% =====================================================================
clear; clc;
N = 160;                                   % total panels
b  = linspace(0, 2*pi, N+1).';
xc = (1 + cos(b))/2;                       % cosine spacing, TE->LE->TE
t  = 0.12;
yt = 5*t*(0.2969*sqrt(xc) - 0.1260*xc - 0.3516*xc.^2 + 0.2843*xc.^3 - 0.1015*xc.^4);
y  = yt;  y(b <= pi) = -yt(b <= pi);       % lower surface first -> clockwise

alphas = [0 4 8];  Cl = zeros(size(alphas));
for k = 1:numel(alphas)
    [~, Cl(k)] = panelMethod(xc, y, alphas(k));
end
slope = (Cl(3) - Cl(1)) / deg2rad(8);
fprintf('NACA 0012 validation\n');
fprintf('  Cl(0 deg) = %+.4f   (theory: 0)\n', Cl(1));
fprintf('  Cl(4 deg) = %+.4f\n', Cl(2));
fprintf('  Cl(8 deg) = %+.4f\n', Cl(3));
fprintf('  slope     = %.2f /rad (theory 2*pi = %.2f; thickness adds ~10%%)\n', slope, 2*pi);
