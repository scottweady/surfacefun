
function [x, y, z, nx, ny, nz] = inclusion(v0, X, xN, yN, zN)
% Constructs an inclusion patch
%   v0 - angle of inclusion
%   X - coordinate of inclusion
%   xN, yN, zN - coordinates of connecting edge

  % Number of edges
  Ne = size(xN, 1);

  % Coordinates
  x = cell(Ne, 1);
  y = cell(Ne, 1);
  z = cell(Ne, 1);

  % Normal vector
  nx = cell(Ne, 1);
  ny = cell(Ne, 1);
  nz = cell(Ne, 1);

  % Move to pole
  [xN, yN, zN] = rotate(xN, yN, zN, X, [0 0 1]);

  % Loop over edges
  for ne = 1 : Ne

    % Discretize
    [tx, ty, tz, tnx, tny, tnz] = buildPatch(v0, xN(ne, :), yN(ne, :), zN(ne, :));

    % Rotate back
    [x{ne}, y{ne}, z{ne}] = rotate(tx, ty, tz, [0 0 1], X);
    [nx{ne}, ny{ne}, nz{ne}] = rotate(tnx, tny, tnz, [0 0 1], X);
    
  end

end

function [x, y, z, nx, ny, nz] = buildPatch(v0, xN, yN, zN)

  % Number of discretization points
  N = length(xN);

  % Get corners in spherical coordinates
  uL = atan2(yN(1), xN(1));
  vL = acos(zN(1));

  uR = atan2(yN(end), xN(end));
  vR = acos(zN(end));

  du = mod(uR - uL, 2 * pi);
  
  % Take shortest path
  if du < pi
    uR = uL + du;
  else
    uR = uL - (2 * pi - du);
  end

  ub = zeros(N);
  vb = zeros(N);

  % Interpolate
  s = chebpts(N, [0 1])';

  ub([1 end], :) = [uL; uR] + 0 * s;
  ub(:, 1) = chebpts(N, [uL uR]);

  vb([1 end], :) = (1 - s) * v0 + s .* [vL; vR];
  vb(:, 1) = v0 + 0 * s;

  % Convert to Cartesian coordinates
  xb = cos(ub) .* sin(vb);
  yb = sin(ub) .* sin(vb);
  zb =            cos(vb);

  % Apply edge boundary condition
  xb(:, end) = xN;
  yb(:, end) = yN;
  zb(:, end) = zN;
  
  % Construct mesh
  [x, y, z] = spherePatch(xb, yb, zb);

  % Normal vector to inclusion
  nx = cos(ub) .* cos(v0);
  ny = sin(ub) .* cos(v0);
  nz = -sin(v0);

  % Normal only defined on edge
  nx(:, 1 : end - 1) = nan;
  ny(:, 1 : end - 1) = nan;
  nz(:, 1 : end - 1) = nan;

end