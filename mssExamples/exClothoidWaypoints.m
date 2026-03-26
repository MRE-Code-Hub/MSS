% exClothoidWaypoints is compatible with MATLAB and GNU Octave (www.octave.org). 
% The script computes a curvature-continuous (G2) path composed of 
% clothoid–arc–clothoid segments. Only waypoint positions are required.
% Optional start and terminal path tangents can be specified.
%
% Outputs:
%   path.x, path.y     sampled path coordinates
%   path.s             arc length
%   path.kappa         curvature profile
%   path.pi_wpt        tangent angles used at waypoints

Lh = 10; % Length of arrow in plots

% Waypoints (North-East coordinates)
Pxy = [ ...
    0,   0
    15,  12
    36, -5
    40, 20
    20, 30
    ];

%% Start and terminal path tangents (optional)
% Use [] to allow automatic tangent selection based on chord direction
piPath0 = deg2rad(10);   
piPathN = [];            

% Trurning radius, maximum curvature, and maximum curvature rate
R_max = 10;
kappaMax = 1 / R_max;
sigmaMax = 0.01;

% Compute clothoid path through the waypoints
path = clothoidPathThroughWaypoints(Pxy,kappaMax,sigmaMax,piPath0,piPathN);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); clf; figure(gcf)

subplot(211)
plot(path.y, path.x, 'LineWidth', 2); hold on;
plot(path.Pxy(:,2), path.Pxy(:,1), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);

% Plot waypoint tangent directions
quiver(path.Pxy(:,2), path.Pxy(:,1), ...
    Lh*sin(path.pi_tangent), ...
    Lh*cos(path.pi_tangent), ...
    0, 'r', 'LineWidth', 1.2);
hold off;
axis equal; grid on;
xlabel('East (m)');
ylabel('North (m)');
title('G2 clothoid path through waypoints (North-East)');

subplot(212)
plot(path.s, path.kappa, 'LineWidth', 2); grid on;
xlabel('s'); ylabel('\kappa(s)');
title('Curvature');

