function [Cp, Cl, xc, yc, V] = panelMethod(x, y, alpha_deg)
% PANELMETHOD  2-D incompressible vortex-panel method (Kuethe & Chow).
%   [Cp,Cl,xc,yc,V] = panelMethod(x,y,alpha_deg) computes the inviscid
%   pressure distribution and lift coefficient of an airfoil using a
%   linearly-varying vortex-panel method with an explicit Kutta condition.
%
%   Method (Kuethe & Chow, "Foundations of Aerodynamics"):
%     * the airfoil surface is divided into flat panels;
%     * each node carries a vortex strength, varying linearly along a panel;
%     * flow tangency is enforced at every panel control point (midpoint);
%     * the Kutta condition gamma_1 + gamma_N = 0 closes the system.
%
%   Inputs :  x, y      - closed airfoil loop, ordered CLOCKWISE (see
%                         loadAirfoil.m). N panels => N+1 nodes.
%             alpha_deg - angle of attack in degrees.
%   Outputs:  Cp  - pressure coefficient at each panel control point
%             Cl  - lift coefficient
%             xc,yc - control-point coordinates
%             V   - tangential surface velocity / Vinf at each control point
%
%   Freestream speed and chord are normalised to 1.

    N  = numel(x) - 1;                 % number of panels
    X  = x(1:end-1);  Y = y(1:end-1);  % panel start nodes
    xc = 0.5*(x(1:end-1) + x(2:end));  % control points (midpoints)
    yc = 0.5*(y(1:end-1) + y(2:end));
    th = atan2(y(2:end) - y(1:end-1), x(2:end) - x(1:end-1));  % panel angles
    S  = hypot(x(2:end) - x(1:end-1), y(2:end) - y(1:end-1));  % panel lengths

    % --- influence coefficients ---
    CN1 = zeros(N); CN2 = zeros(N); CT1 = zeros(N); CT2 = zeros(N);
    for i = 1:N
        for j = 1:N
            if i == j
                CN1(i,j) = -1.0;   CN2(i,j) = 1.0;
                CT1(i,j) = pi/2;   CT2(i,j) = pi/2;
            else
                A = -(xc(i)-X(j))*cos(th(j)) - (yc(i)-Y(j))*sin(th(j));
                B =  (xc(i)-X(j))^2 + (yc(i)-Y(j))^2;
                C =  sin(th(i)-th(j));   D = cos(th(i)-th(j));
                E =  (xc(i)-X(j))*sin(th(j)) - (yc(i)-Y(j))*cos(th(j));
                F =  log(1 + (S(j)^2 + 2*A*S(j))/B);
                G =  atan2(E*S(j), B + A*S(j));
                P =  (xc(i)-X(j))*sin(th(i)-2*th(j)) + (yc(i)-Y(j))*cos(th(i)-2*th(j));
                Q =  (xc(i)-X(j))*cos(th(i)-2*th(j)) - (yc(i)-Y(j))*sin(th(i)-2*th(j));
                CN2(i,j) = D + 0.5*Q*F/S(j) - (A*C + D*E)*G/S(j);
                CN1(i,j) = 0.5*D*F + C*G - CN2(i,j);
                CT2(i,j) = C + 0.5*P*F/S(j) + (A*D - C*E)*G/S(j);
                CT1(i,j) = 0.5*C*F - D*G - CT2(i,j);
            end
        end
    end

    % --- assemble the linear system for the nodal vortex strengths ---
    al  = deg2rad(alpha_deg);
    AN  = zeros(N+1, N+1);   AT = zeros(N, N+1);   RHS = zeros(N+1, 1);
    for i = 1:N
        AN(i,1)   = CN1(i,1);   AN(i,N+1) = CN2(i,N);
        AT(i,1)   = CT1(i,1);   AT(i,N+1) = CT2(i,N);
        for j = 2:N
            AN(i,j) = CN1(i,j) + CN2(i,j-1);
            AT(i,j) = CT1(i,j) + CT2(i,j-1);
        end
        RHS(i) = sin(th(i) - al);       % freestream normal component
    end
    AN(N+1,1) = 1.0;  AN(N+1,N+1) = 1.0;  RHS(N+1) = 0.0;   % Kutta condition

    gamma = AN \ RHS;                   % nodal vortex strengths / (2*pi*Vinf)

    % --- surface velocity, pressure, and lift ---
    V  = cos(th - al) + AT*gamma;       % tangential velocity / Vinf
    Cp = 1 - V.^2;
    gpan = 0.5*(gamma(1:end-1) + gamma(2:end));
    Cl = 4*pi*sum(gpan .* S);           % 2*Gamma/(Vinf*c), Gamma = 2*pi*sum(g*S)
end
