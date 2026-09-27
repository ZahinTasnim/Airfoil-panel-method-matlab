function [x, y] = loadAirfoil(datfile, Nps)
% LOADAIRFOIL  Read a Selig-format airfoil .dat file and re-panel it.
%   [x,y] = loadAirfoil(datfile, Nps) reads the airfoil coordinates,
%   splits them at the leading edge, and re-samples each surface with
%   Nps panels using cosine spacing (points clustered at the leading and
%   trailing edges, where the flow gradients are largest).
%
%   The returned (x,y) is a single CLOSED loop ordered CLOCKWISE
%   (trailing edge -> lower surface -> leading edge -> upper surface ->
%   trailing edge), with a single sharp trailing-edge node.  This is the
%   ordering expected by panelMethod.m.
%
%   Inputs :  datfile - path to a Selig .dat file (1 header line, then x y)
%             Nps     - panels per surface (total panels = 2*Nps)
%   Outputs:  x, y    - node coordinates of the closed airfoil loop

    if nargin < 2, Nps = 100; end

    % --- read coordinates (skip the 1-line header) ---
    M   = readmatrix(datfile, 'NumHeaderLines', 1, 'FileType', 'text');
    xin = M(:,1);  yin = M(:,2);

    % --- split at the leading edge (minimum x) ---
    [~, iLE] = min(xin);
    % Selig order starts at the TE on the upper side, runs to the LE, then
    % back along the lower side to the TE.
    xUp = xin(1:iLE);      yUp = yin(1:iLE);        % TE -> LE (upper)
    xLo = xin(iLE:end);    yLo = yin(iLE:end);      % LE -> TE (lower)

    % --- cosine-spaced re-sampling of each surface ---
    [xUp, yUp] = cosineResample(xUp, yUp, Nps);
    [xLo, yLo] = cosineResample(xLo, yLo, Nps);

    % --- assemble the loop, then order it CLOCKWISE (TE -> lower -> LE ->
    %     upper -> TE), which is the ordering panelMethod.m expects ---
    x = [xUp; xLo(2:end)];           % TE -> upper -> LE -> lower -> TE (CCW)
    y = [yUp; yLo(2:end)];
    x = flipud(x);  y = flipud(y);   % -> clockwise
    x(1)   = 1.0;  y(1)   = 0.0;     % sharp trailing edge
    x(end) = 1.0;  y(end) = 0.0;
end

function [xr, yr] = cosineResample(x, y, Nps)
% Parameterise the surface by normalised cumulative chord length, then
% evaluate a spline at cosine-spaced parameter values.
    t   = [0; cumsum(hypot(diff(x), diff(y)))];
    t   = t / t(end);
    beta = linspace(0, pi, Nps+1).';
    s    = (1 - cos(beta)) / 2;        % cosine spacing in [0,1]
    xr = interp1(t, x, s, 'spline');
    yr = interp1(t, y, s, 'spline');
end
