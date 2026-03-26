
function path = clothoidPathThroughWaypoints( ...
    Pxy, kappaMax, sigmaMax, piPath0, piPathN, opts)
% clothoidPathThroughWaypoints is compatible with MATLAB and GNU Octave 
% (www.octave.org). The function generates a planar G2-continuous reference path
% through a sequence of waypoints using composite clothoid-arc-clothoid segments.
% The input waypoints are specified by position only,
%
%     Pxy = [ x_1 y_1
%             x_2 y_2
%              ...
%             x_N y_N ]
%
% while the path tangent angle at the initial and final waypoint may be
% prescribed optionally as piPath0 and piPathN. Interior waypoint tangent
% angles are computed automatically from neighboring waypoints.
%
% The generated path satisfies:
%
%   - position continuity (G0),
%   - tangent-angle continuity (G1),
%   - curvature continuity (G2),
%
% while the curvature derivative d kappa / ds is piecewise constant inside
% each clothoid segment and may be discontinuous at segment junctions.
%
% The path is assembled in two passes:
%
%   1) Each waypoint interval is solved independently as a composite
%      clothoid-arc-clothoid segment.
%   2) The segment solutions are stitched into one global path with
%      duplicate nodes at the segment interfaces removed.
%
% INPUTS
%   Pxy        = nPts x 2 waypoint matrix, [x_i, y_i]
%   kappaMax   = maximum allowed absolute curvature
%   sigmaMax   = maximum allowed curvature slope, sigma = d kappa / ds
%   piPath0    = optional initial path tangent angle
%   piPathN    = optional final path tangent angle
%   opts       = optional struct passed to clothoid3arc
%
% OUTPUT
%   path       = struct with fields
%       .Pxy        waypoint coordinates (Nx2)
%       .pi_wpt     tangent angles specified by the user (NaN if free)
%       .pi_tangent tangent angles of the solved G2 path at the waypoints
%       .kappaMax   maximum curvature
%       .sigmaMax   maximum curvature slope
%       .sols       cell array of segment solutions
%       .x          global x-coordinates
%       .y          global y-coordinates
%       .s          global arc-length vector
%       .kappa      global curvature vector
%
% EXAMPLES
%   path = clothoidPathThroughWaypoints(Pxy, kappaMax, sigmaMax)
%   path = clothoidPathThroughWaypoints(Pxy, kappaMax, sigmaMax, piPath0)
%   path = clothoidPathThroughWaypoints(Pxy, kappaMax, sigmaMax, piPath0, piPathN)
%   path = clothoidPathThroughWaypoints(Pxy, kappaMax, sigmaMax, piPath0, piPathN, opts)

if nargin < 6
    opts = struct();
end

nPts = size(Pxy,1);

% Curvature values at waypoints: default zero
K = zeros(nPts,1);

% Build waypoint tangent-angle vector
pi_wpt = nan(nPts,1);

% Initial tangent
if ~isempty(piPath0)
    pi_wpt(1) = piPath0;
end

% Interior tangents from neighboring waypoints
if nPts >= 3
    for i = 2:nPts-1
        dx = Pxy(i+1,1) - Pxy(i-1,1);
        dy = Pxy(i+1,2) - Pxy(i-1,2);
        pi_wpt(i) = atan2(dy, dx);
    end
end

% Final tangent
if ~isempty(piPathN)
    pi_wpt(end) = piPathN;
end

% If start/end tangents are free, use chord directions
if isnan(pi_wpt(1)) && nPts >= 2
    dx = Pxy(2,1) - Pxy(1,1);
    dy = Pxy(2,2) - Pxy(1,2);
    pi_wpt(1) = atan2(dy, dx);
end

if isnan(pi_wpt(end)) && nPts >= 2
    dx = Pxy(end,1) - Pxy(end-1,1);
    dy = Pxy(end,2) - Pxy(end-1,2);
    pi_wpt(end) = atan2(dy, dx);
end

% ------------------------------------------------------------------------------
% Pass 1: solve all segments
% ------------------------------------------------------------------------------
sols = cell(nPts-1,1);

for i = 1:nPts-1
    p0 = [Pxy(i,1);   Pxy(i,2);   pi_wpt(i)];
    p1 = [Pxy(i+1,1); Pxy(i+1,2); pi_wpt(i+1)];

    if i == 1
        kStart = K(1);
    else
        kStart = sols{i-1}.kappa(end);
    end

    kEnd = K(i+1);
    piPathActive = ~isnan(pi_wpt(i+1));

    sols{i} = clothoid3arc(p0,p1,kStart,kEnd,kappaMax,sigmaMax,piPathActive,opts);
end

% ------------------------------------------------------------------------------
% Pass 2: preallocate and stitch the full path
% ------------------------------------------------------------------------------
Ntot = length(sols{1}.x);
for i = 2:nPts-1
    Ntot = Ntot + length(sols{i}.x) - 1;
end

xAll     = zeros(Ntot,1);
yAll     = zeros(Ntot,1);
sAll     = zeros(Ntot,1);
kappaAll = zeros(Ntot,1);

idx = 1;
sOffset = 0;

for i = 1:nPts-1
    sol = sols{i};

    if i == 1
        Ni = length(sol.x);

        xAll(idx:idx+Ni-1)     = sol.x;
        yAll(idx:idx+Ni-1)     = sol.y;
        sAll(idx:idx+Ni-1)     = sol.s;
        kappaAll(idx:idx+Ni-1) = sol.kappa;
    else
        Ni = length(sol.x) - 1;

        xAll(idx:idx+Ni-1)     = sol.x(2:end);
        yAll(idx:idx+Ni-1)     = sol.y(2:end);
        sAll(idx:idx+Ni-1)     = sOffset + sol.s(2:end);
        kappaAll(idx:idx+Ni-1) = sol.kappa(2:end);
    end

    idx = idx + Ni;
    sOffset = sAll(idx-1);
end

% Path tangents
path.pi_tangent = pi_wpt;  % start with constraints

% replace free ones by true solved tangents
if isempty(piPath0)
    path.pi_tangent(1) = sols{1}.piPath(1);
end

if isempty(piPathN)
    path.pi_tangent(end) = sols{end}.piPath(end);
end

path.Pxy = Pxy;
path.pi_wpt = pi_wpt;
path.kappaMax = kappaMax;
path.sigmaMax = sigmaMax;
path.sols = sols;
path.x = xAll;
path.y = yAll;
path.s = sAll;
path.kappa = kappaAll;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subfunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sol = clothoid3arc( ...
    p0, p1, kappa0, kappa1, kappaMax, sigmaMax, piPathActive, opts)
% G2 CLOTHOID CONNECTION USING A COMPOSITE 3-SEGMENT PATH
%
% This function computes a planar G2-continuous connection between two poses
%
%     p0 = [x0; y0; piPath0],     curvature kappa0
%     p1 = [x1; y1; piPath1],     curvature kappa1
%
% using a composite path consisting of:
%
%   1) Entry clothoid:
%        curvature varies linearly from kappa0 to kappa_m over length L1
%
%            kappa(s) = kappa0 + sigma1 s
%
%   2) Circular arc:
%        constant curvature kappa_m over length L2
%
%            kappa(s) = kappa_m
%
%   3) Exit clothoid:
%        curvature varies linearly from kappa_m to kappa1 over length L3
%
%            kappa(s) = kappa_m + sigma3 s
%
% The unknown parameter vector
%
%        z = [kappa_m, L1, L2, L3]
%
% is computed numerically using fminsearch by minimizing a cost function
% that penalizes:
%
%   - endpoint position error,
%   - optional endpoint tangent-angle error,
%   - total path length,
%   - curvature-limit violations,
%   - curvature-slope violations,
%   - negative segment lengths.
%
% The endpoint tangent-angle penalty is included only when piPathActive is
% true. This allows the terminal path tangent to be either prescribed or
% left free.
%
% INPUTS
%   p0           = [x0; y0; piPath0], start pose
%   p1           = [x1; y1; piPath1], end pose
%   kappa0       = start curvature
%   kappa1       = end curvature
%   kappaMax     = maximum allowed absolute curvature
%   sigmaMax     = maximum allowed curvature slope, sigma = d kappa / ds
%   piPathActive = logical flag, true if terminal tangent angle is active
%   opts         = optional struct with fields
%       .N1         number of samples in entry clothoid
%       .N2         number of samples in circular arc
%       .N3         number of samples in exit clothoid
%       .maxIter    fminsearch maximum iterations
%       .tolX       fminsearch TolX
%       .tolFun     fminsearch TolFun
%       .wPos       position error weight
%       .wPsi       path-tangent error weight
%       .wLen       total path length weight
%       .wKappa     curvature-limit penalty weight
%       .wSigma     curvature-slope penalty weight
%       .wNeg       negative-length penalty weight
%       .init       initial guess [kappaMid, L1, L2, L3]
%
% OUTPUT
%   sol          = struct with fields
%       .ok          true if final cost is finite
%       .cost        final objective value
%       .x           x-coordinates of the segment
%       .y           y-coordinates of the segment
%       .piPath      path tangent angle along the segment
%       .s           arc-length vector
%       .kappa       curvature along the segment
%       .kappaMid    intermediate curvature kappa_m
%       .L1          entry clothoid length
%       .L2          circular arc length
%       .L3          exit clothoid length
%       .sigma1      entry clothoid curvature slope
%       .sigma3      exit clothoid curvature slope
%       .pEnd        final pose [x_end; y_end; piPath_end]
%       .err         endpoint error [dx; dy; dPiPath]
%
% EXAMPLES
%   sol = clothoid3arc(p0, p1, kappa0, kappa1, kappaMax, sigmaMax, piPathActive)
%   sol = clothoid3arc(p0, p1, kappa0, kappa1, kappaMax, sigmaMax, piPathActive, opts)

if ~isfield(opts,'N1'),      opts.N1 = 80; end
if ~isfield(opts,'N2'),      opts.N2 = 80; end
if ~isfield(opts,'N3'),      opts.N3 = 80; end
if ~isfield(opts,'maxIter'), opts.maxIter = 4000; end
if ~isfield(opts,'tolX'),    opts.tolX = 1e-10; end
if ~isfield(opts,'tolFun'),  opts.tolFun = 1e-10; end
if ~isfield(opts,'wPos'),    opts.wPos = 1e8; end
if ~isfield(opts,'wPsi'),    opts.wPsi = 10; end
if ~isfield(opts,'wLen'),    opts.wLen = 1; end
if ~isfield(opts,'wKappa'),  opts.wKappa = 1e6; end
if ~isfield(opts,'wNeg'),    opts.wNeg = 1e8; end
if ~isfield(opts,'init'),    opts.init = []; end
if ~isfield(opts,'wSigma'),  opts.wSigma = 1e6; end

% Initial guess
D = hypot(p1(1)-p0(1), p1(2)-p0(2));
dpiPath = ssa(p1(3) - p0(3));

if isempty(opts.init)
    kappaMid0 = max(min(2*dpiPath/max(D,1e-3), kappaMax), -kappaMax);
    L10 = max(0.2*D, 0.5);
    L20 = max(0.6*D, 0.5);
    L30 = max(0.2*D, 0.5);
    z0 = [kappaMid0, L10, L20, L30];
else
    z0 = opts.init(:).';
    if numel(z0) ~= 4
        error('opts.init must be [kappaMid L1 L2 L3].');
    end
end

fobj = @(z) objective(z,p0,p1,kappa0,kappa1,kappaMax,sigmaMax,piPathActive,opts);

fmopts = optimset('Display','off', ...
    'MaxIter',opts.maxIter, ...
    'MaxFunEvals',10*opts.maxIter, ...
    'TolX',opts.tolX, ...
    'TolFun',opts.tolFun);

z = fminsearch(fobj, z0, fmopts);

% Build final solution
[J, out] = objective(z,p0,p1,kappa0,kappa1,kappaMax,sigmaMax,piPathActive,opts);

sol = out;
sol.ok   = isfinite(J);
sol.cost = J;

end

% ------------------------------------------------------------------------------
function [J, out] = objective(z,p0,p1,k0,k1,kMax,sigmaMax,piPathActive,opts)
% Cost function used by fminsearch.
%
% Propagates the 3-segment clothoid-arc-clothoid path defined by
% z = [kappaMid L1 L2 L3], computes endpoint error and penalties for
% curvature bounds, curvature-rate bounds and negative lengths.
%
% Returns scalar cost J and struct 'out' with sampled path solution.

kappaMid = z(1);
L1 = z(2);
L2 = z(3);
L3 = z(4);

% Penalize invalid lengths
negPenalty = sum(max(-[L1 L2 L3], 0).^2);

L1 = max(L1, 1e-6);
L2 = max(L2, 1e-6);
L3 = max(L3, 1e-6);

% Curvature slopes for clothoids
sigma1 = (kappaMid - k0) / L1;
sigma3 = (k1 - kappaMid) / L3;

% Sigma satuaration
sigmaViol1 = max(abs(sigma1) - sigmaMax, 0);
sigmaViol3 = max(abs(sigma3) - sigmaMax, 0);
pSigma = sigmaViol1^2 + sigmaViol3^2;

% Propagate 3 segments
[x1, y1, piPath1, s1, k1s] = propagateClothoid( ...
    p0(1), p0(2), p0(3), k0, sigma1, L1, opts.N1);

[x2, y2, piPath2, s2, k2s] = propagateArc( ...
    x1(end), y1(end), piPath1(end), kappaMid, L2, opts.N2);

[x3, y3, piPath3, s3, k3s] = propagateClothoid( ...
    x2(end), y2(end), piPath2(end), kappaMid, sigma3, L3, opts.N3);

% Full path
x = [x1; x2(2:end); x3(2:end)];
y = [y1; y2(2:end); y3(2:end)];
piPath = [piPath1; piPath2(2:end); piPath3(2:end)];
s = [s1; L1+s2(2:end); L1+L2+s3(2:end)];
kappa = [k1s; k2s(2:end); k3s(2:end)];

% Endpoint error
dx = x(end) - p1(1);
dy = y(end) - p1(2);
dpiPath = ssa(piPath(end) - p1(3));

% Curvature penalty
kappaViol = max(abs(kappa) - kMax, 0);
pKappa = sum(kappaViol.^2);

% Cost
if piPathActive
    piPathCost = opts.wPsi*(dpiPath^2);
else
    piPathCost = 0;
end

J = opts.wPos*(dx^2 + dy^2) + ...
    piPathCost + ...
    opts.wLen*(L1 + L2 + L3) + ...
    opts.wKappa*pKappa + ...
    opts.wSigma*pSigma + ...
    opts.wNeg*negPenalty;

out = struct();
out.x = x;
out.y = y;
out.piPath = piPath;
out.s = s;
out.kappa = kappa;
out.kappaMid = kappaMid;
out.L1 = L1;
out.L2 = L2;
out.L3 = L3;
out.sigma1 = sigma1;
out.sigma3 = sigma3;
out.pEnd = [x(end); y(end); piPath(end)];
out.err = [dx; dy; dpiPath];

end

% ------------------------------------------------------------------------------
function [x,y,piPath,s,kappa] = propagateClothoid(x0,y0,piPath0,kappa0,sigma,L,N)
% Numerically integrates a clothoid segment
%   kappa(s) = kappa0 + sigma*s ,  s ∈ [0,L]
% using trapezoidal integration of the heading angle.

s = linspace(0,L,N).';
ds = s(2)-s(1);

kappa = kappa0 + sigma*s;
piPath = piPath0 + kappa0*s + 0.5*sigma*s.^2;

x = zeros(N,1);
y = zeros(N,1);
x(1) = x0;
y(1) = y0;

for i = 2:N
    c1 = cos(piPath(i-1));
    c2 = cos(piPath(i));
    s1 = sin(piPath(i-1));
    s2 = sin(piPath(i));
    x(i) = x(i-1) + 0.5*ds*(c1 + c2);
    y(i) = y(i-1) + 0.5*ds*(s1 + s2);
end

end

% ------------------------------------------------------------------------------
function [x,y,piPath,s,kappa] = propagateArc(x0,y0,piPath0,kappa0,L,N)
% Numerically integrates a constant-curvature circular arc
% over arc length s ∈ [0,L].

s = linspace(0, L, N).';
ds = s(2)-s(1);

kappa = kappa0 * ones(N,1);
piPath = piPath0 + kappa0*s;

x = zeros(N,1);
y = zeros(N,1);
x(1) = x0;
y(1) = y0;

for i = 2:N
    c1 = cos(piPath(i-1));
    c2 = cos(piPath(i));
    s1 = sin(piPath(i-1));
    s2 = sin(piPath(i));
    x(i) = x(i-1) + 0.5*ds*(c1 + c2);
    y(i) = y(i-1) + 0.5*ds*(s1 + s2);
end

end

