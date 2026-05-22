
close all
addpath('utils')

% Mesh order
p = 8;

% Create sphere mesh
[x, y, z] = cubedSphere(p, 4);

% Sphere parameterization
X = @(u, v) [cos(u) .* sin(v) sin(u) .* sin(v) cos(v)];

% Number of inclusions
N_i = 5;

% Radius of inclusion in spherical coordinates
v0 = 0.05 * pi;

% Random set of points (spherical coordinates)
rng(1);
u = linspace(0, 2 * pi, N_i + 1)';
u = u(1 : end - 1);
v = pi / 8 + (pi - 2 * pi / 8) * rand(size(u));

% Inclusion coordinates
X_i = X(u, v);

% Generate mesh
for n_i = 1 : N_i

  clc, fprintf('forming inclusion %d of %d\n', n_i, N_i)

  % Get index of intersecting patches
  P = locatePatches(v0, X_i(n_i, :), x, y, z);

  % Preallocate
  xedge = [];
  yedge = []; 
  zedge = [];

  % Group edges of intersection
  for p = P
    [xe, ye, ze] = locateEdges(v0, X_i(n_i, :), x{p}, y{p}, z{p});
    xedge = [xedge; xe];
    yedge = [yedge; ye];
    zedge = [zedge; ze];
  end

  % Discretize inclusion
  [x_i, y_i, z_i] = inclusion(v0, X_i(n_i, :), xedge, yedge, zedge);

  % Delete original patch
  x(P) = []; 
  y(P) = []; 
  z(P) = [];

  % Current number of elements
  Np = length(x);

  % Append discretization
  for m_i = 1 : length(x_i)
    x{Np + m_i} = x_i{m_i};
    y{Np + m_i} = y_i{m_i};
    z{Np + m_i} = z_i{m_i};
  end

end

% Create surface mesh
S = surfacemesh(x, y, z);

% Solve a Poisson problem
pdo = struct('lap', -1);
L = surfaceop(S, pdo, 1);

% Visualize solution
u = L.solve(0);
plot(u)
colorbar
